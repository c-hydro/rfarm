
"""
Library Features:

Name:          lib_io_src_tiff
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20210104'
Version:       '1.0.0'
"""

# -------------------------------------------------------------------------------------
# Libraries
import logging
import os
import rasterio
from rasterio.crs import CRS

from copy import deepcopy

import numpy as np
import pandas as pd

from rfarm.settings.lib_args import logger_name

log_stream = logging.getLogger(logger_name)

# Debug
# import matplotlib.pylab as plt
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to read data for moloc
def read_data_moloc(file_list, time_list,
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
        else:
            log_stream.warning(' ====> File not found [' + file_step + '] - Following files are not included!')
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
