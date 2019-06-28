from setuptools import setup

req_packages = ['numpy',
                'scipy',
                'pyproj',
                'pandas',
                'GDAL',
                'geopandas',
                'fiona',
                'matplotlib>=2.0.0',
                'scikit-image',
                'Pillow',
                'joblib',
                'netCDF4',
                'shapely',
                'rasterio',
                'configobj',
                'pytest',
                'xarray',
                'progressbar2',
                'boto3',
                'requests',
                'salem']

setup(
	  # Project info
	  name='initialization',
      version='0.1',
      description='initialization of past glacier states',
	  # Github repostiory link
      url='http://github.com/OGGM/initialization',
	  # Authors details
      author='OGGM Developers',
      author_email='jeis@uni-bremen.de',
	  # License
      license='GPLv3+',
	  # What does the project relate to?
	  keywords=['geosciences', 'glaciers', 'initialization'],
      packages=['initialization'],
	  # Python 3 only 
	  python_requires='>=3.5',
	  # Install all dependencies --> same as for OGGM
	  install_requires=req_packages,
      zip_safe=False)

