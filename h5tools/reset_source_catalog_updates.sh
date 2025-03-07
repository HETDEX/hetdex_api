#!/usr/bin/env python

# archive existing source_catalog_manual_updates.fits and .audit
# create new empty table and audit and set permissions
# can only be successfully executed as the hetdex user

import os
import shutil
from h5tools import source_catalog_manual_updates
import traceback
from filelock import FileLock


fail  = False
sct = source_catalog_manual_updates.SrcCatUpdateTable()
T = sct.make_new_table()
sct.set_paths()


try:


  if os.path.isfile(sct.table_fqfn): #assum the audit also exists
    dirname = sct.utc_now_str()
    dirname = dirname.replace(":","").replace(" ","_")
    os.mkdir(dirname)
    #get fullpath to dirname
    dirname_fqfn = os.path.abspath(dirname)

    #move the tabel file and audit file to the new location
    #NOTE this WILL fail if anyone has an active lock?
    lock = FileLock(sct.lock_fqfn)
    with lock:
      try:
        shutil.move(sct.table_fqfn, os.path.join(dirname_fqfn,os.path.basename(sct.table_fqfn)))
        shutil.move(sct.audit_fqfn, os.path.join(dirname_fqfn,os.path.basename(sct.audit_fqfn)))
      except:
          print("move failed")
          print(traceback.format_exc())
          fail = True

  if not fail:
    #set permissions
    print("Creating new files")
    try:
      T.write(sct.table_fqfn)
      with open(sct.audit_fqfn,"w+") as f:
        f.write(f"new audit created: {sct.utc_now_str()}\n")

      os.chmod(sct.table_fqfn,0o666)
      os.chmod(sct.audit_fqfn,0o666)
    except:
      print("new file creation")
      print(traceback.format_exc())
      fail = True
  else:
    print("Failure. Will not create new files.")
except:
  print("directory creation failed")
  print(traceback.format_exc())


