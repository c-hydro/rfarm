"""
Library Features:

Name:          lib_io_dst_tiff
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20201102'
Version:       '3.0.0'
"""
#################################################################################
# Library
import logging
from osgeo import gdal, gdalconst

from rfarm.settings.lib_args import logger_name

# logging
log_stream = logging.getLogger(logger_name)
#################################################################################


# -------------------------------------------------------------------------------------
# Method to write geotiff file
def write_dset_tiff(file_name, map_data, map_geo_x, map_geo_y, map_n=1,
                    map_metadata=None, map_proj='EPSG:4326', map_format=gdalconst.GDT_Float32):

    # map metadata
    if map_metadata is None:
        map_metadata = {'description_field': 'data'}

    # map geotransform
    map_high = map_data.shape[0]
    map_wide = map_data.shape[1]
    map_x_min, map_y_min, map_x_max, map_y_max = [map_geo_x.min(), map_geo_y.min(), map_geo_x.max(), map_geo_y.max()]
    map_x_res = (map_x_max - map_x_min) / float(map_high)
    map_y_res = (map_y_max - map_y_min) / float(map_wide)
    map_geotransform = (map_x_min, map_x_res, 0, map_y_max, 0, -map_y_res)

    # map handle
    map_handle = gdal.GetDriverByName('GTiff').Create(
        file_name, map_wide, map_high, map_n, map_format, options=['COMPRESS=DEFLATE'])
    # map geographical info
    map_handle.SetGeoTransform(map_geotransform)
    map_handle.SetProjection(map_proj)
    # map data
    map_handle.GetRasterBand(1).WriteArray(map_data)
    map_handle.GetRasterBand(1).SetMetadata(map_metadata)

    del map_handle
# -------------------------------------------------------------------------------------
