
"""
Library Features:

Name:          lib_io_src_grib
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20210104'
Version:       '1.0.0'
"""

# -------------------------------------------------------------------------------------
# Libraries
import logging
import numpy as np
# import numpy.ma as ma
import pandas as pd
import xarray as xr
# import pygrib

from scipy.interpolate import griddata

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


# -------------------------------------------------------------------------------------
# Method to adjust lami 2i datasets
def adjust_data_lami_2i(da_var, da_time, da_geo_x, da_geo_y,
                        tag_var_time='valid_time', tag_var_geo_x='longitude', tag_var_geo_y='latitude',
                        tag_dim_time='step', tag_dim_geo_x='x', tag_dim_geo_y='y'
                        ):

    da_dims = list(da_var.dims)
    idx_dim_time = da_dims.index(tag_dim_time)

    if idx_dim_time == 0:
        if isinstance(da_geo_y, xr.DataArray) and isinstance((da_geo_x, xr.DataArray)):
            geo_y_upper = da_geo_y.values[0, 0]
            geo_y_lower = da_geo_y.values[-1, 0]
            geo_y_values = da_geo_y.values
            geo_x_values = da_geo_x.values
        else:
            geo_y_upper = da_geo_y[0, 0]
            geo_y_lower = da_geo_y[-1, 0]
            geo_y_values = da_geo_y
            geo_x_values = da_geo_x

        if geo_y_lower > geo_y_upper:
            values_geo_y = np.flipud(geo_y_values)
            values_geo_x = geo_x_values

            if isinstance(da_var, xr.DataArray):
                values_raw = da_var.values
            else:
                values_raw = da_var

            values_var = np.zeros(shape=[values_raw.shape[1], values_raw.shape[2], values_raw.shape[0]])
            values_var[:, :, :] = np.nan
            for i in range(0, values_raw.shape[0]):
                values_tmp = np.flipud(values_raw[i, :, :])
                values_var[:, :, i] = values_tmp
        else:
            values_geo_y = da_geo_y.values
            values_geo_x = da_geo_x.values
            values_var = da_var.values

    else:
        log_stream.error(' ==> Time dimension index is not in a allowed position.')
        raise IOError(' ==> Check the position of time dimension!')

    values_time = pd.DatetimeIndex(da_time.values)

    return values_var,  values_time, values_geo_x, values_geo_y
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to read data for lami-2i
def read_data_lami_2i(file_name, data_filters=None, data_var='tp',
                          tag_var_time='valid_time', tag_var_geo_x='longitude', tag_var_geo_y='latitude',
                          tag_dim_time='step', tag_dim_geo_x='x', tag_dim_geo_y='y'):

    # Starting info
    log_stream.info(' --> Open file ' + file_name + ' ... ')

    # Define data filter for default configuration
    if data_filters is None:
        data_filters = {'filter_by_keys': {'typeOfLevel': 'surface'}}

    # Open datasets
    dst = xr.open_dataset(file_name, engine='cfgrib', backend_kwargs=data_filters)

    '''
    file_handle = pygrib.open(file_name)
    file_msg = file_handle.select(nameECMF="Total Precipitation")

    msg_n = file_msg.__len__()

    file_geo_x = None
    file_geo_y = None
    file_data = None
    for msg_id, msg_step in enumerate(file_msg):

        if (file_geo_x is None) and (file_geo_y is None):
            file_geo_y, file_geo_x = msg_step.latlons()

        data_step_tmp = msg_step.values
        if ma.is_masked(data_step_tmp):
            data_step_value = data_step_tmp.data
        else:
            data_step_value = deepcopy(data_step_tmp)

        if file_data is None:
            file_data = np.zeros([data_step_value.shape[0], data_step_value.shape[1], msg_n])
            file_data[:, :, :] = np.nan
        file_data[:, :, msg_id] = data_step_value
    '''

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

    interp_method = 'nearest'
    var_values = da_var.values
    var_geo_x = da_geo_x.values
    var_geo_y = da_geo_y.values
    var_time = da_time.values

    var_high = var_values.shape[2]
    var_wide = var_values.shape[1]
    var_x_min, var_y_min, var_x_max, var_y_max = [var_geo_x.min(), var_geo_y.min(), var_geo_x.max(), var_geo_y.max()]
    var_x_res = (var_x_max - var_x_min) / float(var_wide)
    var_y_res = (var_y_max - var_y_min) / float(var_high)
    var_res = np.mean([var_x_res, var_y_res])

    var_geo_x_tmp = np.arange(var_x_min, var_x_max, var_res, float)
    var_geo_y_tmp = np.arange(var_y_min, var_y_max, var_res, float)

    #var_geo_x_tmp = np.linspace(var_x_min, var_x_max, var_wide, endpoint=True)
    #var_geo_y_tmp = np.linspace(var_y_min, var_y_max, var_high, endpoint=True)

    var_geo_x_cmp, var_geo_y_cmp = np.meshgrid(var_geo_x_tmp, var_geo_y_tmp)
    #var_geo_y_cmp = np.flipud(var_geo_y_cmp)

    cmp_values = np.zeros([var_values.shape[0], var_geo_x_cmp.shape[0], var_geo_y_cmp.shape[1]])
    cmp_values[:, :, :] = np.nan
    for step in range(0, var_values.shape[0]):
        step_values = var_values[step, :, :]

        tmp_values = griddata((var_geo_x.ravel(), var_geo_y.ravel()),
                               step_values.ravel(),
                               (var_geo_x_cmp, var_geo_y_cmp), method='nearest', fill_value=-9999.0)

        cmp_values[step, :, :] = tmp_values

        '''
        # DEBUG
        plt.figure()
        plt.imshow(step_values)
        plt.colorbar()
        plt.figure()
        plt.imshow(tmp_values)
        plt.colorbar()
        plt.show()
        '''

    '''
    var_high_cmp = var_geo_x_cmp.shape[1]
    var_wide_cmp = var_geo_y_cmp.shape[0]

    var_x_min_cmp, var_y_min_cmp, var_x_max_cmp, var_y_max_cmp = [var_geo_x_cmp.min(), var_geo_y_cmp.min(), var_geo_x_cmp.max(), var_geo_y_cmp.max()]
    var_x_res_cmp = (var_x_max_cmp - var_x_min_cmp) / float(var_high_cmp)
    var_y_res_cmp = (var_y_max_cmp - var_y_min_cmp) / float(var_wide_cmp)
    var_res_cmp = np.mean([var_x_res_cmp, var_y_res_cmp])
    var_mask_cmp = np.ones([var_high_cmp, var_wide_cmp])

    var_geo_cmp = xr.DataArray(var_mask_cmp, name='mask', dims=[tag_var_geo_y, tag_var_geo_x],
                               coords={
                                   tag_var_geo_x: ([tag_var_geo_x], var_geo_x_cmp[:, 0]),
                                   tag_var_geo_y: ([tag_var_geo_y], var_geo_y_cmp[0, :])})
    interp_dict = {tag_var_geo_y: var_geo_cmp[tag_var_geo_y], tag_var_geo_x: var_geo_cmp[tag_var_geo_x]}

    da_var_interp = var_values_cmp.interp(interp_dict, method=interp_method)
    '''

    da_var_cmp = xr.DataArray(cmp_values, name=data_var, dims=[tag_dim_time, tag_dim_geo_y, tag_dim_geo_x],
                              coords={
                                    tag_dim_time: ([tag_dim_time], var_time),
                                    tag_dim_geo_x: ([tag_dim_geo_x], var_geo_x_cmp[0, :]),
                                    tag_dim_geo_y: ([tag_dim_geo_y], var_geo_y_cmp[:, 0])})

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

    return da_var_cmp, da_time, var_geo_x_cmp, var_geo_y_cmp
# -------------------------------------------------------------------------------------
