#!/usr/bin/env python

# Import required python modules
import ftplib
import os
import tarfile

# Define the local directory name to put data in
ddir=r"C:\data\3_fire_data\burned_area\fire_cci"

# Change the local directory to where you want to put the data
os.chdir(ddir)

# login to FTP
f=ftplib.FTP("ftp.ceda.ac.uk", "nrietze", "FNrMg^e5Gn5M")

# loop through years
for year in range(2019,2020):
    
    # create a new directory for the year's data
    ydir = os.path.join(ddir, str(year))
    os.makedirs(ydir, exist_ok=True)
    
    # loop through months
    for month in range(4, 5):
        # change the remote directory
        f.cwd("/neodc/esacci/fire/data/burned_area/MODIS/pixel/v5.1/compressed/%.4d" % year )
        
        # get the remote files to the local directory
        filenames = f.nlst()
        for filename in filenames:
            if filename.endswith(".tar.gz") and ("AREA_1" in filename or "AREA_3" in filename or "AREA_4" in filename) and filename.startswith("%d%.2d" % (year, month)):
                f.retrbinary("RETR %s" % filename, open(filename, "wb").write)
                with tarfile.open(filename, "r:gz") as tar:
                    tar.extractall(path = ydir)
                os.remove(filename)

# Close FTP connection
f.close()