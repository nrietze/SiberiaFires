@echo off

set "outputFolder=./warped_tiles/"

:: Step 1: Warp and process individual TIF files
for /F %%i in ('dir /b /s "*.tif"') do (
   echo %%i
   gdalwarp -t_srs "+proj=laea +lat_0=90 +lon_0=-40 +x_0=0 +y_0=0" -tr 30 30 -r mode  %%i ./warped_tiles/%%~ni_warped.tif
)

:: Step 2: Merge warped TIF files using gdal_merge.py
cd /d ./warped_tiles/
echo Merging warped TIF files...

dir /b /s *.tif > tiff_list.txt
gdal_merge.py -o merged.tif -of GTiff --optfile tiff_list.txt

cd /d ..
