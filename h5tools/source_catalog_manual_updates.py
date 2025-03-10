"""
definitions and helpers for the manual updates table for the source catalog

Note: permissions should work for ALL users on TACC, but will NOT work if the path is a remote ssh mount
(that is, you MUST be logged in to a TACC cluster or on the jupyter-hub)

!!! For Jupyter-hub users ...
you may need to run:   !pip install filelock
then restart your kernel. This will hold until you restart the server.

"""

import os
import numpy as np
from astropy.table import Table
try:
    from filelock import FileLock
except:
    print("For Jupyter-hub users ...")
    print("you may need to run:   !pip install filelock")
    print("then restart your kernel. This will hold until you restart the server.")

import copy
import traceback
import datetime
from datetime import timezone


#shared
SOURCECAT_UPDATES_ROOTPATHS_DICT = {"cluster":"/work/05350/ecooper/wrangler/team_classify/shared",
                                    "hub":"/home/jovyan/shared"}
SOURCECAT_UPDATES_SUBPATH = "hdrX/source_catalog_updates/"
SOURCECAT_UPDATES_FILE = "source_catalog_manual_updates.fits"
SOURCECAT_UPDATES_AUDIT = "source_catalog_manual_updates.audit"
SOURCECAT_UPDATES_LOCK = "source_catalog_manual_updates.lock"


#corral
# SOURCECAT_UPDATES_ROOTPATHS_DICT = {"cluster":"/corral-repl/utexas/", "hub":"/home/jovyan/"}
# SOURCECAT_UPDATES_SUBPATH = "Hobby-Eberly-Telesco/hdrX/source_catalog_updates/"
# SOURCECAT_UPDATES_FILE = "source_catalog_manual_updates.fits"
# SOURCECAT_UPDATES_AUDIT = "source_catalog_manual_updates.audit"
# SOURCECAT_UPDATES_LOCK = "source_catalog_manual_updates.lock"


#test
# SOURCECAT_UPDATES_ROOTPATHS_DICT = {"cluster":"/scratch/03261/polonius/temp/", "hub":"/home/jovyan/"}
# SOURCECAT_UPDATES_SUBPATH = "Hobby-Eberly-Telesco/hdrX/source_catalog_updates/"
# SOURCECAT_UPDATES_FILE = "source_catalog_manual_updates.fits"
# SOURCECAT_UPDATES_AUDIT = "source_catalog_manual_updates.audit"
# SOURCECAT_UPDATES_LOCK = "source_catalog_manual_updates.lock"




class SrcCatUpdateTable:
    def __init__(self): #,table_fqfn=None, audit_fqfn=None):

        # COULD allow user to specify paths/files, but for now at least, keep it simple and only use the defined paths

        self.table_fqfn = None
        self.audit_fqfn = None
        self.lock_fqfn = None
        self.rootpath = None

        self.status = 0
        self.status_msg = None

        self.the_table = None

        self.the_user = None #the user who instantiates

        self.updatable_cols = ['conf_status',
                               'z',
                               'parent_detectid',
                               'classification_labels',
                               'comments' ] #noteice: revision, revision_user, revision_date auto-updated


    def make_new_table(self):

        # the update table
        # key column explanation
        # revision = incrementing integer of how many revisions (updates) have been made to THIS row
        # revision_user = username of the most recent person to make an update
        # revision_date = date/time of the most recent update as YYYYMMDD HH:MM:SS in UTC
        # conf_status = confirmation status, e.g. confirmed_line, confirmed_absorber, false positive, etc
        # parent_detectid = the detectid (NOT the source_id) of a detection that represents the singular object to which
        #                  this detection belongs (and from which the redshift is determined)
        # classification_labels = string of 0 or more comma separated labels (e.g. AGN, LzG, LAB, etc)
        # comments = free form comment string from the last user (up to 256 characters)

        return Table(dtype=[('detectid', np.int64),('shotid', int), ('ra',float), ('dec',float),('obs_wave',float),
                 ('revision',int),('revision_user','S32'),('revision_date','S16'),
                 ('conf_status', int),('z',float),('parent_detectid', int),
                 ('classification_labels','S32'),('comments','S256')])

    def utc_now_str(self):
        return datetime.datetime.now(timezone.utc).strftime("%Y%m%d %H:%M:%S")


    def get_row_dict(self):
        dict = {}
        if self.the_table is not None:
            T = self.the_table
        else:
            T = self.make_new_table()

        for col in T.colnames:
            dict[col] = None

        return dict

    def reset_status(self):
        self.status = 0
        self.status_msg = ""

    def get_status(self):
        return self.status, self.status_msg

    def get_user(self):

        self.reset_status()
        try:
            userid = os.getuid()
        except:
            userid = None

        try:
            self.the_user = os.getenv('JUPYTERHUB_USER')
            if self.the_user is None and userid is not None:
                self.the_user = 'unknown'
                import pwd
                pwdinfo = pwd.getpwuid(userid)
                self.the_user = pwdinfo.pw_name
        except:
            try:
                if userid is not None and self.the_user != "unknown":
                    import pwd
                    pwdinfo = pwd.getpwuid(userid)
                    self.the_user = pwdinfo.pw_name
            except:
                self.the_user = None

        if self.the_user is None:
            self.status = -1
            self.status_msg = "Cannot identify the user"


    def set_paths(self):
        """

        Returns
        -------

        """

        try:

            if os.path.exists(os.path.join(SOURCECAT_UPDATES_ROOTPATHS_DICT['hub'],SOURCECAT_UPDATES_SUBPATH)):
                path = os.path.join(SOURCECAT_UPDATES_ROOTPATHS_DICT['hub'],SOURCECAT_UPDATES_SUBPATH)
                self.rootpath = SOURCECAT_UPDATES_ROOTPATHS_DICT['hub']
            elif os.path.exists(os.path.join(SOURCECAT_UPDATES_ROOTPATHS_DICT['cluster'],SOURCECAT_UPDATES_SUBPATH)):
                path = os.path.join(SOURCECAT_UPDATES_ROOTPATHS_DICT['cluster'],SOURCECAT_UPDATES_SUBPATH)
                self.rootpath = SOURCECAT_UPDATES_ROOTPATHS_DICT['cluster']
            else:
                path = None
                self.rootpath = None

            if path is not None:
                self.table_fqfn = os.path.join(path,SOURCECAT_UPDATES_FILE)
                self.audit_fqfn = os.path.join(path,SOURCECAT_UPDATES_AUDIT)
                self.lock_fqfn = os.path.join(path,SOURCECAT_UPDATES_LOCK)
                if not os.path.isfile(self.table_fqfn): #has priority
                    self.status = -1
                    self.status_msg = f"Cannot locate update table {self.table_fqfn}"
                elif not os.path.isfile(self.audit_fqfn): #has priority
                    self.status = -1
                    self.status_msg = f"Cannot locate audit {self.audit_fqfn}"
                elif not os.path.isfile(self.lock_fqfn):
                    self.status = -1
                    self.status_msg = f"Cannot locate lock {self.lock_fqfn}"

            else:
                self.status = -1
                self.status_msg = "Cannot locate path(s) to table, audit, lock"

        except:
            self.status = -1
            self.status_msg = "Cannot locate path(s) to table, audit, lock"
            self.status_msg += f"\n{traceback.format_exc()}"

    def load(self):
        """
        load the table

        """
        self.reset_status()
        try:
            if self.the_table:
                try:
                    del self.the_table
                except:
                    pass

            self.the_table = Table.read(self.table_fqfn,format="fits")

        except:
            self.status = -1
            self.status_msg = "Cannot load/refresh update table"
            self.status_msg += f"\n{traceback.format_exc()}"


    def get_row(self,detectid,refresh=False):
        """
        organized by detectid, though other columns are searchable

        todo: future allow get_row by other mechanisms (RA, Dec, etc?), but these might not be UNIQUE

        Parameters
        ----------
        detectid

        Returns row and its index
        -------

        """
        row = None
        index = -1
        self.reset_status()
        try:
            if not self.the_table or refresh is True:
                self.load()

            if not self.the_table:
                return row, index

            #this is a read only, no lock needed
            #NOTE: though, table COULD be stale ... user can force a refresh

            sel = np.array(self.the_table['detectid']==detectid)
            ct = np.count_nonzero(sel)
            if ct == 0:
                row = None
                index = -1
                self.status = 0 #this is not an error ... there just was no matching row
                self.status_msg = f"No matching record found for detectid {detectid}"
            elif ct > 1:
                row = None
                index = -1
                self.status = -1
                self.status_msg = f"Error! {ct} matching records found for detectid {detectid}"
            else:
                row = copy.copy(self.the_table[sel][0]) #simple copy is fine. no complex structures
                index = list(sel).index(True)

        except:
            row = None
            index = -1
            self.status = -1
            self.status_msg = f"Unable to search on detectid {detectid}"
            self.status_msg += f"\n{traceback.format_exc()}"

        return row, index



    def add_audit(self,row,action=0):
        """

        Parameters
        ----------
        row
        action: -1 = UNDO (e.g. undo the previous audit if an update failed)
                 0 = UPDATE existing row
                 1 = ADD new row

        Returns
        -------

        """
        try:
            with open(self.audit_fqfn,"a") as f:
                f.write(f"{action} ")
                for col in self.get_row_dict().keys():
                    f.write(f"{col}:{row[col]}\t")
                f.write("\n")

            return 0
        except:
            try:
                detectid = row['detectid']
            except:
                detectid = "???"

            #Audit does not overwrite the status, just appends
            #BUT really this should not happen
            self.status_msg += f"Unable to audit for {detectid}"
            self.status_msg += f"\n{traceback.format_exc()}"
            return -1

    def add_row(self,row):
        """
        add a new row

        Parameters
        ----------
        row

        Returns
        -------

        """
        try:
            row['revision'] = 1
            row['revision_date'] = self.utc_now_str()
            row['revision_user'] = self.the_user

            if self.add_audit(row,1) == 0:
                try:
                    self.the_table.add_row(row) #or do I need to break it down by column

                    try:
                        self.the_table.write(self.table_fqfn, format='fits', overwrite=True)
                    except:
                        self.status = -1
                        self.status_msg = f"Table file write failed."
                        self.status_msg += f"\n{traceback.format_exc()}"
                        self.add_audit(row, -1) #return the staus from the table write attempt
                except:
                    self.status = -1
                    self.status_msg = f"Add row failed."
                    self.status_msg += f"\n{traceback.format_exc()}"
                    self.add_audit(row, -1)  # return the staus from the table write attempt
            else:
                self.add_audit(row,-1)
                try:
                    detectid = row['detectid']
                except:
                    detectid = "???"
                self.status = -1
                #use the status message from add_audit() but in this case, set the status flag


        except:
            try:
                detectid = row['detectid']
            except:
                detectid = "???"
            self.status = -1
            self.status_msg = f"Unable to add row for {detectid}"
            self.status_msg += f"\n{traceback.format_exc()}"

    def update_row(self,row,idx):
        """
        update changed columns in an existing row
        update revision
        Parameters
        ----------
        row

        Returns
        -------

        """
        try:
            row['revision'] = self.the_table[idx]['revision'] + 1
            row['revision_date'] = self.utc_now_str()
            row['revision_user'] = self.the_user

            if self.add_audit(row) == 0:

                self.the_table[idx]['revision'] = row['revision']
                self.the_table[idx]['revision_date'] = row['revision_date']
                self.the_table[idx]['revision_user'] = row['revision_user']

                for col in self.updatable_cols:
                    self.the_table[idx][col] = row[col]

                try:
                    self.the_table.write(self.table_fqfn, format='fits', overwrite=True)
                except:
                    self.status = -1
                    self.status_msg = f"Table file write failed."
                    self.status_msg += f"\n{traceback.format_exc()}"

                    self.add_audit(row, -1)  # return the staus from the table write attempt
            else:
                self.add_audit(row,-1)
                try:
                    detectid = row['detectid']
                except:
                    detectid = "???"
                self.status = -1
                self.status_msg = f"Audit failed. Unable to add row for {detectid}"
                self.status_msg += f"\n{traceback.format_exc()}"

        except:
            try:
                detectid = row['detectid']
            except:
                detectid = "???"
            self.status = -1
            self.status_msg = f"Unable to update row for {detectid}"
            self.status_msg += f"\n{traceback.format_exc()}"


    def update_table(self, row):
        """
        Single row update
        may update an existing row or insert a new row
        row can be an actual table row or a dictionary with matching keys
        Returns
        -------

        """
        self.reset_status()
        try:

            if self.the_user is None:
                self.get_user()
                if self.the_user is None:
                    self.status = -1
                    self.status_msg = "Cannot identify the user"
                    return

            if row is None:
                self.status = -1
                self.status_msg = "Bad update"
                return

            #get a lock
            lock = FileLock(self.lock_fqfn)

            with lock:
                try:
                    #MUST reload the table after we get a lock
                    self.load()

                    if self.the_table is None:
                        self.status = -1
                        self.status_msg += f"\nCould not (re)load table to assert condition."
                        return  #should release the lock


                    current_row, current_index = self.get_row(row['detectid'],refresh=False)

                    if current_row is None:
                        if self.status != 0: #this is a problem
                            self.status = -1
                            self.status_msg += f"Error! Could not update table."
                            return
                        else: # this is just a new row to be added
                            self.reset_status()
                    else: #found a row, make sure we are current
                        if current_row['revision'] != row['revision']:
                            self.status = -1
                            self.status_msg += f"Error! Could not update table. Row record is stale. {current_row['revision']} != {row['revision']}"
                            return

                    #all good; either adding a row or updting an existing row
                    if current_row is None:
                        #this is a new insert
                        self.add_row(row)

                    else:
                        #this is an update
                        self.update_row(row,current_index)


                except:
                    self.status = -1
                    self.status_msg = "Could not modify the update table"
                    self.status_msg += f"\n{traceback.format_exc()}"

        except:
            self.status = -1
            self.status_msg = "Could not modify the update table"
            self.status_msg += f"\n{traceback.format_exc()}"




