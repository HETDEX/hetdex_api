import random
import numpy as np
import time
from scipy.spatial import cKDTree
from astropy.table import Table


# Author: Niv Drory
# NOTE, reqires scipy>1.4 due to a bug in KDTree.query_ball_point() in prior versions!
#


def mktree(x, y, z, dsky=3.0, dwave=5.0, euclidean=False):
    """
   Construct a kd-tree and calculate the linking length in normalized coordinates (where
   the distance is on a sphere)

   Parameters
   ----------
      x,y,z : numpy.array
         coordinates RA, DEC in degrees, WAVELENGTH in Angstroms
      dksy : float
         linking length on sky in arcsec, defaults to 3.0"
      dwave : float
         linking length in wavelength in Angstrom, defaults to 5A
      euclidean : bool
         coordinates are euclidean, not spherical, do not transform. 
         defaults to False

   Returns
   -------
      kd, r : scipy.spatial.cKDTree object, spherical linking length
   """
    # construct a kd-tree in normalized coordinates and calculate the linking length in that system
    # return cKDTree, r

    if euclidean == True:
        x_tree = x
        y_tree = y
        z_tree = z
    else:
        # assume x=ra/deg, y=dec/deg, z=wavelength/A
        x_tree = np.deg2rad(x)  # convert to standard spherical angle
        y_tree = np.deg2rad(90.0 - y)
        zf = (
            np.pi / 180 / (5500.0 - 3500.0) * dsky / dwave
        )  # scale wavelength to a "radians"-like scale
        z_tree = zf * (z - 3500.0)
        dsky = np.deg2rad(dsky / 3600.0)  # convert from arcsec to radians

    data = np.vstack((x_tree, y_tree, z_tree)).T
    kd = cKDTree(data)

    return kd, dsky


def frinds_of_friends(kdtree, r, Nmin=3):
    """
   Friends of Friends group finder

   For each particle, n, not yet assigned to any group:

       Create member list for the group
       Add particle n to group
       For each particle p in group
           Find neighbors of p within the linking length
              add those to the group if not already members of any group
       Retain the group if it has more than one member

   Example
   -------
      id, x, y, z, w = ...
      kdtree, r = mktree(x,y,z)
      group_lst = frinds_of_friends(kdtree, r, Nmin=3)

      # filter and sort groups (optional)
      group_lst = [l for l in group_lst if len(l) > 10]
      group_lst = sorted(group_lst, key=lambda a:len(a))

      # get group parameters:
      group_table = process_group_list(group_lst, x, y, z, w)  

      # save group table (load with load_groups())
      save_groups('groups', group_table)

   Parameters
   ----------
      kdtree : scipy.spatial.kdtree
         KD-Tree containing the coordinates of all the particles

      r : float
         linking length

      Nmin : int
         minimal number of members for a group

   Returns
   -------
      group_lst : list of lists
         List of the groups, where each group is represented by a list of the indices
         of all its members. The indices refer to the original coordinate arrays
         the KD Tree was constructed of, using the function mktree().
         The list is sorted by group size, with the largest group first.
   """
    # make a set with all particles, makes set ops below much faster
    id_lst = set(range(0, kdtree.n - 1))
    group_lst = []
    # iterate over all particles
    while len(id_lst) > 0:
        index = id_lst.pop()

        # start new group with current particle, keep it a list for iteration
        friends_lst = [index]

        i = 0
        # while there are more friends (of friends)
        while i < len(friends_lst):
            # find all friends of current particle
            friends = kdtree.query_ball_point(kdtree.data[friends_lst[i]], r)
            if len(friends) > 1:
                friends = set(friends)
                # restrict friends to those not already in groups (linking is commutative)
                friends = friends & id_lst
                # remove the friends from the particles list
                id_lst = id_lst - friends
                # add unique new friends to the group
                friends_lst.extend(friends)

            i = i + 1

        # finally, if we have found more than one particle in the group, add
        # to the group list
        if len(friends_lst) > 0:
            group_lst.append(friends_lst)

    # filter and sort groups
    group_lst = [l for l in group_lst if len(l) >= Nmin]
    group_lst = sorted(group_lst, key=lambda a: len(a))
    group_lst.reverse()  # largest group first

    return group_lst


def evaluate_group(x, y, z, f, euclidean=False):
    """
   Calculate group properties from a list of it's members
   x, y, z, and f (flux or weight) arrays.

   Return a tuple of group properties:
     (lum, icx, icy, icz, ixx, iyy, ixy, izz, a, b, ka, kb, pa)

     lum : total luminosity (sum of f)
     icx, icy, icz : f weighted centroid
     ixx, iyy, ixy : f weighted second-order moments in first two (sky) dimenstions
     izz : f weighted second order moment along third dimension
     a, b : f weighted semi-major and semi-minor axis in the x/y plane
     pa : position angle on sky
     a2, b2 : non-weighted semi-major and semi-minor axis in the x/y plane
     pa2 : position angle in degrees
   """
    # calculate some aggregate properties most useful for further processing

    # first and second order flux-weighed moments
    lum = np.sum(f)
    icy = np.sum(f * y) / lum
    c = 1.0 if euclidean else np.cos(np.deg2rad(icy))  # cos-dec
    icx = np.sum(f * x) / lum
    icz = np.sum(f * z) / lum
    dx = (x - icx) * c
    dy = y - icy
    ixx = np.sum(f * np.square(dx)) / lum
    iyy = np.sum(f * np.square(dy)) / lum
    ixy = np.sum(f * dx * dy) / lum
    izz = np.sum(f * np.square(z - icz)) / lum

    # semi-major and semi-minor axes
    a_ = 0.5 * (ixx + iyy)
    c_ = np.sqrt(0.25 * np.square(ixx - iyy) + np.square(ixy))
    a = np.sqrt(a_ + c_)
    b = np.sqrt(a_ - c_) if a_ - c_ > 0 else a_

    # position angle in degree
    pa = 0.5 * np.arctan2(2.0 * ixy, ixx - iyy) if ixx - iyy != 0.0 else np.pi / 4.0
    pa *= 180.0 / np.pi

    # kron-like aperture
    # t = ixx*iyy-ixy*ixy
    # cxx = ixx/t
    # cyy = iyy/t
    # cxy = -2.0*ixy/t
    # ir1 = np.sum(f*a*np.sqrt(cxx*dx*dx + cyy*dy*dy + cxy*dx*dy))/lum
    # fac = a/b if a/b < 4.0 else 4.0
    # ka = ir1*fac
    # kb = ir1/fac

    ixx2 = np.average(np.square(dx))
    iyy2 = np.average(np.square(dy))
    ixy2 = np.average(dx * dy)
    izz2 = np.average(np.square(z - icz))

    # semi-major and semi-minor axes
    a2_ = 0.5 * (ixx2 + iyy2)
    c2_ = np.sqrt(0.25 * np.square(ixx2 - iyy2) + np.square(ixy2))
    a2 = np.sqrt(a2_ + c2_)
    b2 = np.sqrt(a2_ - c2_) if a2_ - c2_ > 0 else a2_

    pa2 = (
        0.5 * np.arctan2(2.0 * ixy2, ixx2 - iyy2) if ixx2 - iyy2 != 0.0 else np.pi / 4.0
    )
    pa2 *= 180.0 / np.pi

    return (lum, icx, icy, icz, ixx, iyy, ixy, izz2, a, b, pa, a2, b2, pa2)


def process_group_list(group_lst, detectid, x, y, z, f):
    """
   Post-process the group memebership list by evaluating group properties using
   coordinates and fluxes (weights).

   Parameters
   ----------
      group_lst : list (of lists)
         group list returned by friends_of_friends(). This list contains
         indices into the original coordinate arrays.

      detectid : array like
         detect_ids corresponding to the original x,y,z inputs of friends_of_friends()

      x,y,z : array like
         coordinates corresponding to the detect-ids

      f : array like
         flux (or broadly the weight) of the points x,y,z

   Returns
   -------
      astropy.Table of processed groups with the fields
         ['id', 'size', 'lum', 'icx', 'icy', 'icz', 'ixx', 'iyy', 'ixy', 'izz',\                                          
          'a', 'b', 'pa', 'a2', 'b2', 'pa2', 'numpy array:members']
   """
    # create list of rows
    rows = []
    for i in np.arange(0, len(group_lst), 1):
        n = len(group_lst[i])
        # calculate some aggregate properties most useful for further processing
        m = group_lst[i]
        params = evaluate_group(x[m], y[m], z[m], f[m])
        # finally, an array with the list of all members' detectid
        members = np.array(detectid[m])
        # append to list. this is much quicker than iterating Table.append()
        rows.append((i, n, *params, members))
        # names=['id', 'size', 'lum', 'icx', 'icy', 'icz', 'ixx', 'iyy', 'ixy', 'izz',\
        #'a', 'b', 'pa', 'a2', 'b2', 'pa2', 'members'])
    return group_list_to_table(rows)


def print_groups(group_lst):
    """
   Print a summary of the groups found by friend_of_friends().

   Prints the number, size, average coordinates, and extent as
   rms of the coordinates for all groups in group_lst.

   Parameters
   ----------
      group_lst : list of lists
         group list returned by friend_of_friends().
   """
    for i in np.arange(0, len(group_lst), 1):
        m = group_lst[i]
        print(i, len(m), *m)


def group_list_to_table(group_lst):
    """
   Convert the list of lists representing the groups (output of process_group_list) to
   an astropy table, with one row per group.

   Parameters
   ----------
      group_lst : list (of lists)
         group list returned by process_group_list().

   Returns
   -------
      astropy.Table with one row per group.
   """
    return Table(
        rows=group_lst,
        names=[
            "id",
            "size",
            "lum",
            "icx",
            "icy",
            "icz",
            "ixx",
            "iyy",
            "ixy",
            "izz",
            "a",
            "b",
            "pa",
            "a2",
            "b2",
            "pa2",
            "members",
        ],
    )


def table_to_group_list(table):
    """
   Convert an astropy.Table of groups back to a list of lists

   Needed to be able to save, since astropy.Table can't deal with 
   saving an numpy array as a cell. It claims it can, but it's buggy.
   """
    return [[x for x in table[i]] for i in range(len(table))]


def save_groups(fname, gtable):
    """
   Save an astropy.Table with group data

   Parameters
   ----------
   fname : string
      file name

   gtable : astropy.Table
      table of groups
   """
    np.save(fname, table_to_group_list(gtable))


def load_groups(fname):
    """
   Load an astropy.Table with group data

   Parameters
   ----------
   fname : string
      file name

   Returns
   -------
   gtable : astropy.Table
      table of groups
   """
    return group_list_to_table(np.load(fname))


#
# T E S T
#


def test():
    x = np.random.uniform(0, 1000, 20000)
    y = np.random.uniform(0, 1000, 20000)
    z = x * 0
    print("building tree ...")
    kdtree, r = mktree(x, y, z, dsky=3.0, euclidean=True)
    t0 = time.time()
    print("starting fof ...")
    group_lst = frinds_of_friends(kdtree, r)

    group_lst = [l for l in group_lst if len(l) > 10]
    # sorting groups
    group_lst = sorted(group_lst, key=lambda a: len(a))
    group_lst.reverse()

    print(time.time() - t0, "seconds")

    print_groups(group_lst, x, y, z, z * 0 + 1.0)
    import matplotlib.pyplot as plt

    t0 = group_lst[0]
    plt.plot(x[t0], y[t0], "x")
    plt.show()


if __name__ == "__main__":
    test()
