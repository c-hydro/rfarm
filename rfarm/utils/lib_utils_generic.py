"""
Library Features:

Name:          lib_utils_generic
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20201202'
Version:       '1.0.0'
"""

#######################################################################################
# Library
import os
import xarray as xr

import numpy as np
from datetime import datetime
#######################################################################################


# -------------------------------------------------------------------------------------
# Method to reshape 3d variable
def reshape_var3d(var_values_in):
    var_shape = var_values_in.shape
    var_values_out = np.zeros([var_shape[1], var_shape[2], var_shape[0]])
    for id_step, var_step in enumerate(var_values_in):
        var_value_step = var_values_in[id_step, :, :]
        var_values_out[:, :, id_step] = var_value_step
    return var_values_out
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to create a data array
def create_darray_3d(data, time, geo_x, geo_y, geo_1d=True,
                     dim_key_x='west_east', dim_key_y='south_north', dim_key_time='time',
                     dim_name_x='west_east', dim_name_y='south_north', dim_name_time='time',
                     dims_order=None):

    if dims_order is None:
        dims_order = [dim_name_y, dim_name_x, dim_name_time]

    if geo_1d:
        if geo_x.shape.__len__() == 2:
            geo_x = geo_x[0, :]
        if geo_y.shape.__len__() == 2:
            geo_y = geo_y[:, 0]

        data_da = xr.DataArray(data,
                               dims=dims_order,
                               coords={dim_key_time: ([dim_name_time], time),
                                       dim_key_x: (dim_name_x, geo_x),
                                       dim_key_y: (dim_name_y, geo_y)})
    else:
        data_da = None
    return data_da
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to merge dictionaries (DA RIVEDERE)
def merge_dict(d1, d2):

    dd = {}
    for d in (d1, d2):  # you can list as many input dicts as you want here
        for key, value in iter(d.items()):
            dd[key] = value

    return dd
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to make folder
def make_folder(path_folder):
    if not os.path.exists(path_folder):
        os.makedirs(path_folder)
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to add time in a unfilled string (path or filename)
def fill_tags2string(string_raw, tags_format=None, tags_filling=None):

    apply_tags = False
    if string_raw is not None:
        for tag in list(tags_format.keys()):
            if tag in string_raw:
                apply_tags = True
                break

    if apply_tags:

        # string_filled_1 = string_raw.format(**tags_format)

        for tag_key, tag_value in tags_format.items():
            tag_key = '{' + tag_key + '}'
            if tag_value is not None:
                string_filled = string_raw.replace(tag_key, tag_value)
                string_raw = string_filled

        for tag_format_name, tag_format_value in list(tags_format.items()):

            if tag_format_name in list(tags_filling.keys()):
                tag_filling_value = tags_filling[tag_format_name]
                if tag_filling_value is not None:

                    if isinstance(tag_filling_value, datetime):
                        tag_filling_value = tag_filling_value.strftime(tag_format_value)

                    string_filled = string_filled.replace(tag_format_value, tag_filling_value)

        string_filled = string_filled.replace('//', '/')
        return string_filled
    else:
        return string_raw
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to get dictionary values using a key
def get_dict_values(d, key, value=[]):

    for k, v in iter(d.items()):

        if isinstance(v, dict):
            if k == key:

                for kk, vv in iter(v.items()):
                    temp = [kk, vv]
                    value.append(temp)

            else:
                vf = get_dict_values(v, key, value)

                if isinstance(vf, list):
                    if vf:
                        vf_end = vf[0]
                    else:
                        vf_end = None

                elif isinstance(vf, np.ndarray):
                    vf_end = vf.tolist()
                else:
                    vf_end = vf

                if vf_end not in value:
                    if vf_end:

                        if isinstance(value, list):
                            value.append(vf_end)
                        elif isinstance(value, str):
                            value = [value, vf_end]

                    else:
                        pass
                else:
                    pass

        else:
            if k == key:

                if isinstance(v, np.ndarray):
                    value = v
                else:
                    value = v
    return value
# -------------------------------------------------------------------------------------
