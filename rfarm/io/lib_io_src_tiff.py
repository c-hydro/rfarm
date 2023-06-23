
"""
Library Features:

Name:          lib_io_src_tiff
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20210104'
Version:       '1.0.0'
"""

# -------------------------------------------------------------------------------------
# Libraries
import copy
import logging
import os
import rasterio
from rasterio.crs import CRS

from copy import deepcopy

import numpy as np
import pandas as pd
import xarray as xr

from rfarm.dataset.lib_variables import compute_rain_moloc
from rfarm.utils.lib_utils_generic import reshape_var3d
from rfarm.settings.lib_args import logger_name

log_stream = logging.getLogger(logger_name)

# Debug
import matplotlib.pylab as plt
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to read data for moloc
def read_data_moloc(file_list, time_list,
                    var_name='rain', var_geox='longitude', var_geoy='latitude', var_time='time',
                    var_proj='EPSG:4326', var_type='float32',
                    var_no_data=-9999.0, var_limit_min=0, var_limit_max=None,
                    time_step_start=None, time_step_end=None):

    # Create list of filename
    if isinstance(file_list, str):
        file_list = [file_list]

    file_list_break = None
    file_list_n = file_list.__len__()

    if time_step_start is None:
        time_step_start = 0
    if time_step_end is None:
        time_step_end = file_list_n

    file_available = []
    for file_id, file_step in enumerate(file_list):
        if os.path.exists(file_step):
            file_available.append(file_step)
            file_list_break = False
        else:
            log_stream.warning(' ====> File not found [' + file_step + '] - Following files are not included!')
            file_list_break = True
            break

    file_available_n = file_available.__len__()
    if file_available_n != file_list_n:
        log_stream.warning(' ====> Available file(s) are ' + str(file_available_n) +
                           ' -- Expected file(s) are ' + str(file_list_n))
        log_stream.warning(' ====> Time steps are limited to ' + str(file_available_n))

        file_step_start = time_step_start
        file_step_end = file_available_n
    else:
        file_step_start = deepcopy(time_step_start)
        file_step_end = deepcopy(time_step_end)

    file_selected = file_list[file_step_start:file_step_end]
    time_selected = time_list[file_step_start:file_step_end]
    var_data_xyt, var_data_time = None, []
    for file_id, (file_step, time_step) in enumerate(zip(file_selected, time_selected)):
        file_handle = rasterio.open(file_step)

        file_bounds = file_handle.bounds
        file_res = file_handle.res
        file_transform = file_handle.transform
        file_data = file_handle.read()

        if file_handle.crs is None:
            file_crs = CRS.from_string(var_proj)
        else:
            file_crs = file_handle.crs

        if var_type == 'float32':
            var_data_xy = np.float32(file_data[0, :, :])
        else:
            log_stream.error(' ===> Data type is not allowed.')
            raise NotImplementedError('Case not implemented yet')

        if var_data_xyt is None:
            var_data_xyt = np.zeros(
                shape=[var_data_xy.shape[0], var_data_xy.shape[1], file_selected.__len__()])
            var_data_xyt[:, :, :] = np.nan

        var_data_xy[var_data_xy == var_no_data] = np.nan
        if var_limit_min is not None:
            var_limit_min = np.float32(var_limit_min)
            var_data_xy[var_data_xy < var_limit_min] = np.nan
        if var_limit_max is not None:
            var_limit_max = np.float32(var_limit_max)
            var_data_xy[var_data_xy > var_limit_max] = np.nan

        decimal_round_geo = 7

        file_center_right = file_bounds.right - (file_res[0] / 2)
        file_center_left = file_bounds.left + (file_res[0] / 2)
        file_center_top = file_bounds.top - (file_res[1] / 2)
        file_center_bottom = file_bounds.bottom + (file_res[1] / 2)

        if file_center_bottom > file_center_top:
            log_stream.warning(' ===> Coords "center_bottom": ' + str(file_center_bottom) +
                               ' is greater than "center_top": '
                               + str(file_center_top) + '. Try to inverse the bottom and top coords. ')
            file_center_tmp = file_center_top
            file_center_top = file_center_bottom
            file_center_bottom = file_center_tmp

        var_geox_1d = np.arange(file_center_left, file_center_right + np.abs(file_res[0] / 2),
                                np.abs(file_res[0]), float)
        var_geoy_1d = np.flip(np.arange(file_center_bottom, file_center_top + np.abs(file_res[1] / 2),
                                        np.abs(file_res[1]), float), axis=0)
        var_geox_2d, var_geoy_2d = np.meshgrid(var_geox_1d, var_geoy_1d)

        var_geoy_upper, var_geoy_lower = var_geoy_2d[0, 0], var_geoy_2d[-1, 0]
        if var_geoy_lower > var_geoy_upper:
            var_geox_2d = np.flipud(var_geox_2d)
            var_data_xy = np.flipud(var_data_xy)

        '''
        # debug
        plt.figure(); plt.imshow(var_data_xy); plt.colorbar()
        plt.figure(); plt.imshow(var_lats_2d); plt.colorbar()
        plt.show()
        '''

        var_lons_min_round = round(np.min(var_geox_2d), decimal_round_geo)
        var_lons_max_round = round(np.max(var_geox_2d), decimal_round_geo)
        var_lats_min_round = round(np.min(var_geoy_2d), decimal_round_geo)
        var_lats_max_round = round(np.max(var_geoy_2d), decimal_round_geo)

        file_center_right_round = round(file_center_right, decimal_round_geo)
        file_center_left_round = round(file_center_left, decimal_round_geo)
        file_center_bottom_round = round(file_center_bottom, decimal_round_geo)
        file_center_top_round = round(file_center_top, decimal_round_geo)

        assert var_lons_min_round == file_center_left_round
        assert var_lons_max_round == file_center_right_round
        assert var_lats_min_round == file_center_bottom_round
        assert var_lats_max_round == file_center_top_round

        var_data_xyt[:, :, file_id] = var_data_xy
        var_data_time.append(time_step)

    var_time_idx = pd.date_range(start=var_data_time[0], end=var_data_time[-1], periods=file_list_n)

    return var_data_xyt, var_time_idx, var_geox_2d, var_geoy_2d
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to convert data for moloc
def convert_data_moloc(var_data_raw, var_units, var_type_feat='accumulated'):

    var_data_cmp = compute_rain_moloc(var_data_raw, oVarUnits=[var_units],  oVarType=[var_type_feat])

    var_data_n = var_data_cmp.shape[2]

    var_data_def = np.zeros([var_data_cmp.shape[0], var_data_cmp.shape[1], var_data_cmp.shape[2]])
    for i in range(0, var_data_n):
        var_data_step = var_data_cmp[:, :, i]

        if np.any(var_data_step < 0):
            log_stream.warning(' ====> Some values are less then 0.0 at step ' + str(i + 1) + ' -- Set to 0.0')
            var_data_step[var_data_step <= 0.0] = 0.0

        if np.all(np.isnan(var_data_step)):
            log_stream.warning(' ====> All values are NaNs at step ' + str(i + 1) + ' -- Set to 0.0')
            var_data_step[:, :] = 0.0
        elif np.any(np.isnan(var_data_step)):
            log_stream.warning(' ====> Some values are NaNs at step ' + str(i + 1) + ' -- Set to 0.0')
            var_data_step[np.isnan(var_data_step)] = 0.0

        var_data_def[:, :, i] = var_data_step

    return var_data_def
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to convert time for wrf
def convert_time_wrf(var_time_cmp, var_time_exp=None, var_type_feat='accumulated'):

    if var_type_feat == 'accumulated':
        var_time_cmp = var_time_cmp[1:]
    elif var_type_feat == 'istantaneous':
        var_time_cmp = var_time_cmp[0:]
    else:
        log_stream.error(' ====> Data type are not defined correctly')
        raise TypeError

    if var_time_exp.equals(var_time_cmp):
        var_time_series = var_time_exp
    else:
        var_freq_inferred = var_time_cmp.inferred_freq

        if var_freq_inferred is not None:
            var_time_series = pd.date_range(start=var_time_cmp[0],
                                            periods=list(var_time_cmp.values).__len__(),
                                            freq=var_freq_inferred)
        else:
            log_stream.error(' ====> Data frequency is not defined correctly')
            raise TypeError(' Data frequency is equal to ' + str(var_freq_inferred))

    return var_time_series
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to read gfs025 forecast
def read_data_gfs_025(file_name, tag_time='time', tag_geo_x='lon', tag_geo_y='lat',
                      var_units='kg m**-2', var_step_type=None):

    # Parse args
    #var_name = list(var_name)[0]
    #var_units = var_units[0]
    #var_step_type = var_step_type[0]

    # Starting info
    log_stream.info(' --> Open file ' + file_name + ' ... ')

    # Open datasets
    dst_tmp = xr.open_dataset(file_name)

    if dst_tmp.dims.__len__() > 3:
        if 'height' in list(dst_tmp.dims):
            dst = dst_tmp.squeeze('height')
            dst = dst.drop('height')
        else:
            log_stream.warning(
                ' ==> Datasets has more then 3 dimensions. Add the name of extra-variable to remove it from datasets')
    else:
        dst = dst_tmp

    # Get variables ALL and DATA
    var_list_all = list(dst.variables)
    var_list_data = list(dst.data_vars)

    # Get time, geo x and geo y
    log_stream.info(' --->  Get time, geo_x and geo_y data ... ')
    if tag_time in var_list_all:
        da_time = dst[tag_time]
    else:
        log_stream.error(' ==> Time dimension name is not in the variables list of grib file')
        raise IOError(' ==> Check the time dimension!')
    if tag_geo_x in var_list_all:
        da_geo_x = dst[tag_geo_x]
    else:
        log_stream.error(' ==> GeoX dimension name is not in the variables list of grib file')
        raise IOError(' ==> Check the GeoX dimension!')
    if tag_geo_y in var_list_all:
        da_geo_y = dst[tag_geo_y]
    else:
        log_stream.error(' ==> GeoY dimension name is not in the variables list of grib file')
        raise IOError(' ==> Check the GeoY dimension!')
    log_stream.info(' --->  Get time, geo_x and geo_y data ... DONE')

    # Get data
    da_var = []
    for var_list_step in var_list_data:
        log_stream.info(' --->  Get ' + var_list_step + ' data ... ')
        da_step = dst[var_list_step]
        da_var.append(da_step)
        log_stream.info(' --->  Get ' + var_list_step + ' data ... DONE')

    # Ending info
    log_stream.info(' --> Open file ' + file_name + ' ... DONE')

    # Start Debug
    #mat = da_values[0].values
    #plt.figure()
    #plt.imshow(mat[0,:,:])
    #plt.colorbar()
    #plt.show()
    # End Debug

    if var_step_type is None:
        var_step_type = ['accum']
    if var_units is None:
        var_units = ['m']

    # Get values
    var_da_in = da_var[0]
    var_values_in = var_da_in.values
    var_dims_in = var_da_in.dims

    var_time = da_time
    var_geo_x = da_geo_x
    var_geo_y = da_geo_y

    if (var_units == 'kg m**-2') or (var_units == 'Kg m**-2'):
        var_units = 'mm'

    if var_units == 'm':
        var_scale_factor = 1000
    elif var_units == 'mm':
        var_scale_factor = 1
    else:
        log_stream.error(' ===> Rain components units are not allowed! Check your data!')
        raise IOError('Selected units are not allowed!')

    if (var_dims_in[0] == 'step') or (var_dims_in[0] == 'time'):
        var_values_in = reshape_var3d(var_values_in)
    var_shape_in = var_values_in.shape

    # Check attributes
    if not (var_units == 'mm') and not (var_units == 'm'):
        log_stream.error(' ===> Rain components units are not allowed! Check your data!')
        raise IOError('Data units is not allowed!')
    if (var_step_type == 'accum') or (var_step_type == 'accumulated'):

        var_values_step_start = None
        var_values_out = np.zeros([var_shape_in[0], var_shape_in[1], var_shape_in[2]])
        var_values_out[:, :, :] = np.nan

        for var_step in range(0, var_shape_in[2]):

            var_values_step_tmp = var_values_in[:, :, var_step]

            if var_values_step_start is None:
                var_values_step_end = var_values_step_tmp
                var_values_step = var_values_step_end
                var_values_step_start = var_values_step_end
            else:
                var_values_step_end = var_values_step_tmp
                var_values_step = var_values_step_end - var_values_step_start
                var_values_step_start = var_values_step_end

            var_values_step[var_values_step < 0.0] = 0.0
            var_values_out[:, :, var_step] = var_values_step / var_scale_factor

            [var_geox_2d, var_geoy_2d] = np.meshgrid(var_geo_x, var_geo_y, sparse=False, indexing='xy')

    elif (var_step_type == 'inst') or (var_step_type == 'instantaneous'):
        var_values_out = copy.deepcopy(var_values_in)
        var_values_in[var_values_in < 0.0] = 0.0
        [var_geox_2d, var_geoy_2d] = np.meshgrid(var_geo_x, var_geo_y, sparse=False, indexing='xy')

    # DEBUG START
    #import matplotlib.pylab as plt
    #plt.figure()
    #plt.imshow(var_geox_2d)
    #plt.colorbar();
    #plt.figure()
    #plt.imshow(var_geoy_2d)
    ##plt.colorbar()
    #plt.figure()
    #plt.imshow(var_values_out[:, :, 0])
    #plt.clim(0, 20)
    #plt.show(block=True)
    # DEBUG END

    return var_values_out, var_time, var_geox_2d, var_geoy_2d
# -------------------------------------------------------------------------------------