import random
import numpy as np
import time
from astropy.io import ascii
from scipy.spatial import cKDTree
from astropy.table import Table, Column

# i=np.where(np.isin(clustering.D.detectid,m,assume_unique=True))[0]

#
# NOTE, reqires scipy>1.4 due to a bug in KDTree.query_ball_point() in prior versions!
#

def mktree(x, y, z, dsky=3.0, dwave=5.0, euclidean=False):
   '''
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
   '''
   # construct a kd-tree in normalized coordinates and calculate the linking length in that system
   # return cKDTree, r

   if euclidean==False:
      dsky = dsky/3600.0
      c = np.cos(y/180.0*np.pi)
      ra_tree = x*c
      dec_tree = y
      # need to be on a sphere
      zf = dsky/dwave          # scale wavelength so distance is in the same units as sky
      z_tree = zf * z
   else:
      ra_tree = x
      dec_tree = y
      z_tree = z

   data = np.vstack((ra_tree, dec_tree, z_tree)).T
   kd = cKDTree(data)

   r = dsky #np.sqrt(3*dsky*dsky)
   return kd, r


def frinds_of_friends(kdtree, r):
   '''
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
      group_lst = frinds_of_friends(kdtree, r)
      group_lst = [l for l in group_lst if len(l) > 10]
      group_lst = sorted(group_lst, key=lambda a:len(a))
      t = group_list_to_table(group_lst, id, x, y, z, w)
      t.write('groups.ecsv', overwrite=True)

   Parameters
   ----------
      kdtree : scipy.spatial.kdtree
         KD-Tree containing the coordinates of all the particles

      r : float
         linking length

   Returns
   -------
      group_lst : list of lists
         List of the groups, where each group is represented by a list of the indices
         of all its members. The indices refer to the original coordinate arrays
         the KD Tree was constructed of, using the function mktree().
   '''
   # make a set with all particles, makes set ops below much faster
   id_lst = set(range(0, kdtree.n-1))
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
      if len(friends_lst)>0:
         group_lst.append(friends_lst)

   return group_lst


def evaluate_group(x, y, z, f):
   '''
   Calculate group properties from a list of it's members
   x, y, z, and f (flux or weight) arrays.

   Return a tuple of group properties:
     (lum, icx, icy, icz, ixx, iyy, ixy, izz, a, b, ka, kb, pa)

     lum : total luminosity (sum of f)
     icx, icy, icz : f weighted centroid
     ixx, iyy, ixy : f weighted second-order moments in first two (sky) dimenstions
     izz : f weighted second order moment along third dimension
     a, b : semi-major and semi-minor axis in the x/y plane
     ka, kb : kron-like (first order moment) adaptive semi-major and minor axes in the x/y plane
     pa : position angle in degrees
   '''
   # calculate some aggregate properties most useful for further processing

   # first and second order flux-weighed moments
   lum = np.sum(f)
   icx = np.sum(f*x)/lum
   icy = np.sum(f*y)/lum
   icz = np.sum(f*z)/lum
   ixx = np.sum(f*np.square(x-icx))/lum
   iyy = np.sum(f*np.square(y-icy))/lum
   ixy = np.sum(f*(x-icx)*(y-icy))/lum
   izz = np.sum(f*np.square(z-icz))/lum

   # semi-major and semi-minor axes
   a_ = 0.5*(ixx+iyy);
   c_ = np.sqrt(0.25*(ixx-iyy)*(ixx-iyy)+ixy*ixy);
   a = np.sqrt (a_+c_);
   b = np.sqrt (a_-c_) if a_-c_>0 else a_

   # position angle in degree
   pa = 0.5*np.arctan2 (2.0*ixy, ixx-iyy) if ixx-iyy!=0.0 else np.pi/4.;
   pa *= 180./np.pi;

   # kron-like aperture
   t = ixx*iyy-ixy*ixy;
   cxx = ixx/t;
   cyy = iyy/t;
   cxy = -2.0*ixy/t;
   dx = x-icx
   dy = y-icy
   ir1 = np.sum(f*a*np.sqrt(cxx*dx*dx + cyy*dy*dy + cxy*dx*dy))/lum
   ka = ir1*a/b
   kb = ir1*b/a

   return (lum, icx, icy, icz, ixx, iyy, ixy, izz, a, b, ka, kb, pa)


def print_groups(group_lst, x, y, z, f):
   '''
   Print a summary of the groups found by friend_of_friends().

   Prints the number, size, average coordinates, and extent as
   rms of the coordinates for all groups in group_lst.

   Parameters
   ----------
      group_lst : list of lists
         group list returned by friend_of_friends().

      x, y, z: numpy.array
         arrays with coordinates used to construct the kd-tree
         using mktree().
   '''
   for i in np.arange(0, len(group_lst), 1):
      m = group_lst[i]
      print(i, len(m), *evaluate_group(x[m], y[m], z[m], f[m]))


def group_list_to_table(group_lst, detectid, x, y, z, f):
   '''
   Convert the list of lists representing the groups (output of friends_of_friends) to
   an astropy table, with one row per group.

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
      astropy.Table with one row per group.
   '''
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
      rows.append((i,n,*params,members))

   return Table(rows=rows, names=['id', 'size', 'lum', 'icx', 'icy', 'icz', 'ixx', 'iyy', 'ixy', 'izz',\
                                  'a', 'b', 'ka', 'kb', 'pa', 'members'])
#
# T E S T
#


def test():
   x = np.random.uniform(0,1000,20000)
   y = np.random.uniform(0,1000,20000)
   z = x*0
   print('building tree ...')
   kdtree, r = mktree(x,y,z,dsky=3.0,euclidean=True)
   t0 = time.time()
   print('starting fof ...')
   group_lst = frinds_of_friends(kdtree, r)

   group_lst = [l for l in group_lst if len(l) > 10]
   # sorting groups
   group_lst = sorted(group_lst, key=lambda a:len(a))
   group_lst.reverse()

   print(time.time() - t0, "seconds")

   print_groups(group_lst,x,y,z,z*0+1.0)
   import matplotlib.pyplot as plt
   t0 = group_lst[0]
   plt.plot(x[t0], y[t0], 'x')
   plt.show()

if __name__ == "__main__":
   test()
