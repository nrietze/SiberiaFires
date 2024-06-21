import cdsapi
import tarfile
import os 

c = cdsapi.Client()


# %% Loop requests


for year in range(2019,2023):
    ydir = rf'C:\data\3_fire_data\burned_area\c3sba_11\{year}'
    
    os.makedirs(ydir, exist_ok=True)
    
    for month in range(4,12):
        
        if os.path.exists(os.path.join(ydir, f'{year}{month:02d}01-C3S-L3S_FIRE-BA-OLCI-AREA_1-fv1.1.nc')):
            print('This data already exists. Downloading next time frame...', end = '\n')
            
            continue
        
        filename = rf'C:\data\3_fire_data\burned_area\c3sba_11\c3sba_{year}_{month}.tar.gz'
        
        c.retrieve(
            'satellite-fire-burned-area',
            {
                'format': 'tgz',
                'origin': 'c3s',
                'sensor': 'olci',
                'variable': 'pixel_variables',
                'version': '1_1',
                'region': [
                    'asia', 'europe', 'north_america',
                ],
                'year': [
                    str(year),
                ],
                'month': [
                    '{:02d}'.format(month),
                ],
                'nominal_day': '01',
            },
            filename)
        
        with tarfile.open(filename, "r:gz") as tar:
            tar.extractall(path = ydir)
        os.remove(filename)