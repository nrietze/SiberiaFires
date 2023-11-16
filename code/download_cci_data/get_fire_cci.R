library(RCurl)

# Define the local directory name to put data in
ddir <- "C:/data/3_fire_data/burned_area/fire_cci"

HOSTNAME <- 'ftp.ceda.ac.uk'

USERNAME <- 'nrietze'
PASSWORD <- 'FNrMg^e5Gn5M'

# loop through years
for (year in 2001:2020) {
  
  # create a new directory for the year's data
  ydir <- file.path(ddir, as.character(year))
  
  if (dir.exists(ydir)){
    print('Folder exists. Interrupting download.')
    next
  }
  
  dir.create(ydir, showWarnings = FALSE)
  
  f.URL <- paste0("ftp://", USERNAME, ":", PASSWORD, "@", HOSTNAME,"/neodc/esacci/fire/data/burned_area/MODIS/pixel/v5.1/compressed/", year, "/")
  
  # get the remote files to the local directory
  filenames <- getURL(url = f.URL, 
                      # userpwd = paste(USERNAME, PASSWORD, sep = ":"),
                      ftp.use.epsv = FALSE, dirlistonly = TRUE)
  filenames <- strsplit(filenames, "\r*\n")[[1]]
  
  # loop through months
  for (month in 4:11) {
    print('Processing: ')
    
    cat("Downloading data from:", month, year, "\n")
    
    # define the remote file name pattern
    pattern <- sprintf("%d%02d01", year, month)
    filenames_filtered <- filenames[grepl(pattern, filenames) & grepl("AREA_[134]", filenames)]
    
    for (filename in filenames_filtered) {
      local_filename <- file.path(ydir, filename)
      download.file(paste0(f.URL, "/", filename), 
                    destfile = local_filename, 
                    mode = "wb",
                    quiet = TRUE)
      untar(local_filename, exdir = ydir)
      file.remove(local_filename)
    }
  }
  
  Sys.sleep(10)
}

