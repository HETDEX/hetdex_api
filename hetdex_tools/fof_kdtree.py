import random
import numpy as np
import time
from astropy.io import ascii
from scipy.spatial import cKDTree


def read_catalog(filen):
   '''
   Read a hetdex ascii catalog of detections and return data selected on vis_class

   Parameters
   ----------
   filen : str
      name of catalog file
   
   Returns
   -------
   astropy.table.Table object representing the catalog
   '''
   data = ascii.read(filen)
   return data[data['vis_class'] >= 4]


def get_coordinates(table):
   '''
   Return the ra, dec, wavelength and ids of detections

   Parameters
   ----------
   table : astropy.table.Table
      Table object
   
   Returns
   -------
      x,y,z,id : numpy.array with RA, DEC, WAVELENGTH and DETECT_ID

   '''
   x = table['ra']
   y = table['dec']
   z = table['wave']
   id = table['detectid']
   return x, y, z, id


def mktree(x, y, z, dsky=3.0, dwave=20.0):
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
         linking length in wavelength in Angstrom, defaults to 20A
   Returns
   -------
      kd, r : scipy.spatial.cKDTree object, spherical linking length
   '''
   # construct a kd-tree in normalized coordinates and calculate the linking length in that system
   # return cKDTree, r

   dsky = dsky/3600.0

   c = np.cos(y/180.0*np.pi)
   ra_tree = x*c
   dec_tree = y
   # need to be on a sphere
   zf = dsky/dwave          # scale wavelength so distance is in the same units as sky
   z_tree = zf * z
   data = np.vstack((ra_tree, dec_tree, z_tree)).T
   kd = cKDTree(data)

   r = np.sqrt(3*dsky*dsky)
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

      kdtree, r = mktree(x,y,z)
      group_lst = frinds_of_friends(kdtree, r)
      group_lst = [l for l in group_lst if len(l) > 10]
      group_lst = sorted(group_lst, key=lambda a:len(a))

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


def print_groups(group_lst, x, y, z):
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
      n = len(group_lst[i])
      avg_ra = np.mean([x[j] for j in group_lst[i]])
      avg_dec = np.mean([y[j] for j in group_lst[i]])
      avg_z = np.mean([z[j] for j in group_lst[i]])
      rms_ra = np.std([x[j] for j in group_lst[i]])
      rms_dec = np.std([y[j] for j in group_lst[i]])
      rms_z = np.std([z[j] for j in group_lst[i]])

      print(i, len(
         group_lst[i]), avg_ra, avg_dec, avg_z, rms_ra, rms_dec, rms_z)

#
# T E S T
#


def test():
   print('reading catalog ...')
   x, y, z, detect_id = get_coordinates(read_catalog('hdr1_sngt6pt5_for.tab'))
   print('building tree ...')
   #x = x[1:5000]
   #y = y[1:5000]
   #z = z[1:5000]
   kdtree, r = mktree(x,y,z)
   t0 = time.time()
   print('starting fof ...')
   group_lst = frinds_of_friends(kdtree, r)

   group_lst = [l for l in group_lst if len(l) > 10]
   # sorting groups
   group_lst = sorted(group_lst, key=lambda a:len(a))
   group_lst.reverse()

   print(time.time() - t0, "seconds")

   print_groups(group_lst,x,y,z)


if __name__ == "__main__":
   test()
