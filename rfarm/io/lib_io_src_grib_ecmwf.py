
"""
Library Features:

Name:          lib_io_src_grib_ecmwf
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20210104'
Version:       '1.0.0'
"""

# -------------------------------------------------------------------------------------
# Libraries
import logging
import numpy as np
import pandas as pd
import xarray as xr

from rfarm.settings.lib_args import time_format, logger_name

# Logging
logging.getLogger("cfgrib").setLevel(logging.WARNING)
log_stream = logging.getLogger(logger_name)

# Debug
# import matplotlib.pylab as plt
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to scale time
def time_scaling(time_range_in, time_range_out):
    time_size_in = time_range_in.shape[0]
    time_size_out = time_range_out.shape[0]

    if time_size_out > time_size_in:
        time_scale = int(time_size_out / time_size_in)
    elif time_size_out == time_size_in:
        time_scale = 1
    else:
        log_stream.error(' ==> Time scale definition not implemented!')
        raise NotImplementedError('Time scale definition not implemented yet!')

    return time_scale

# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to read data for ecmwf 0100
def read_data_ecmwf_0100(file_name, data_filters=None, data_var='tp',
                         tag_var_time='valid_time', tag_var_geo_x='longitude', tag_var_geo_y='latitude',
                         tag_dim_time='step', tag_dim_geo_x='longitude', tag_dim_geo_y='latitude'):

    # Starting info
    log_stream.info(' --> Open file ' + file_name + ' ... ')

    # Define data filter for default configuration
    if data_filters is None:
        data_filters = {'filter_by_keys': {'dataType': 'fc'}}

    # Open datasets
    dst = xr.open_dataset(file_name, engine='cfgrib', backend_kwargs=data_filters)

    # Get variables ALL and DATA
    var_list_all = list(dst.variables)
    var_list_data = list(dst.data_vars)

    # Get time, geo x and geo y
    log_stream.info(' --->  Get time, geo_x and geo_y data ... ')
    if tag_var_time in var_list_all:
        da_time = dst[tag_var_time]
    else:
        log_stream.error(' ==> Time dimension name is not in the variables list of grib file')
        raise IOError(' ==> Check the time dimension!')
    if tag_var_geo_x in var_list_all:
        da_geo_x = dst[tag_var_geo_x]
    else:
        log_stream.error(' ==> GeoX dimension name is not in the variables list of grib file')
        raise IOError(' ==> Check the GeoX dimension!')
    if tag_var_geo_y in var_list_all:
        da_geo_y = dst[tag_var_geo_y]
    else:
        log_stream.error(' ==> GeoY dimension name is not in the variables list of grib file')
        raise IOError(' ==> Check the GeoY dimension!')
    log_stream.info(' --->  Get time, geo_x and geo_y data ... DONE')

    # Get data
    if data_var in var_list_data:
        da_var = dst[data_var]
    else:
        log_stream.error(' ==> Variable name ' + data_var + ' is not in the variables list of grib file')
        raise IOError(' ==> Check the variable datasets!')

    da_dims = list(da_var.dims)
    if tag_dim_time not in da_dims:
        log_stream.error(' ==> Dimension name ' + tag_dim_time + ' is not in the dimensions list of grib file')
        raise IOError(' ==> Check the dimension name!')
    if tag_dim_geo_x not in da_dims:
        log_stream.error(' ==> Dimension name ' + tag_dim_geo_x + ' is not in the dimensions list of grib file')
        raise IOError(' ==> Check the dimension name!')
    if tag_dim_geo_y not in da_dims:
        log_stream.error(' ==> Dimension name ' + tag_dim_geo_y + ' is not in the dimensions list of grib file')
        raise IOError(' ==> Check the dimension name!')

    # Ending info
    log_stream.info(' --> Open file ' + file_name + ' ... DONE')

    return da_var, da_time, da_geo_x, da_geo_y
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to adjust ecmwf 0100 datasets
def adjust_data_ecmwf_0100(da_var, da_time, da_geo_x, da_geo_y,
                           tag_var_time='valid_time', tag_var_geo_x='longitude', tag_var_geo_y='latitude',
                           tag_dim_time='step', tag_dim_geo_x='longitude', tag_dim_geo_y='latitude'
                           ):

    da_dims = list(da_var.dims)
    idx_dim_time = da_dims.index(tag_dim_time)

    if (list(da_geo_y.dims).__len__() == 1) and (list(da_geo_x.dims).__len__() == 1):
        tmp_point_x = da_geo_x.values
        tmp_point_y = da_geo_y.values
        tmp_geo_x, tmp_geo_y = np.meshgrid(tmp_point_x, tmp_point_y)
    elif (list(da_geo_y.dims).__len__() == 2) and (list(da_geo_x.dims).__len__() == 2):
        tmp_geo_x = da_geo_x.values
        tmp_geo_y = da_geo_y.values
    else:
        log_stream.error(' ==> Longitude and Latitude dimensions are not equal.')
        raise IOError(' ==> Check the dimensions of Longitude and Latitude')

    if idx_dim_time == 0:

        geo_y_upper = tmp_geo_y[0, 0]
        geo_y_lower = tmp_geo_y[-1, 0]
        if geo_y_lower > geo_y_upper:
            values_geo_y = np.flipud(tmp_geo_y)
            values_geo_x = tmp_geo_x
        else:
            values_geo_y = tmp_geo_y
            values_geo_x = tmp_geo_x

        values_raw = da_var.values
        values_var = np.zeros(shape=[values_raw.shape[1], values_raw.shape[2], values_raw.shape[0]])
        values_var[:, :, :] = np.nan
        for i in range(0, values_raw.shape[0]):
            if geo_y_lower > geo_y_upper:
                values_tmp = np.flipud(values_raw[i, :, :])
            else:
                values_tmp = values_raw[i, :, :]
            values_var[:, :, i] = values_tmp

    else:
        log_stream.error(' ==> Time dimension index is not in a allowed position.')
        raise IOError(' ==> Check the position of time dimension!')

    values_time = pd.DatetimeIndex(da_time.values)

    return values_var,  values_time, values_geo_x, values_geo_y
# -------------------------------------------------------------------------------------
