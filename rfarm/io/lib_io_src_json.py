
"""
Library Features:

Name:          lib_io_src_json
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20210104'
Version:       '1.0.0'
"""

# -------------------------------------------------------------------------------------
# Libraries
import logging
import os
# import numpy as np
import pandas as pd
# import xarray as xr

from copy import deepcopy

from rfarm.settings.lib_args import time_format, logger_name
log_stream = logging.getLogger(logger_name)

# Debug
# import matplotlib.pylab as plt
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to convert expert forecast time period
def convert_time_expert_forecast_OLD(var_time_idx, var_time_freq='H', var_time_period=12):

    var_time_idx = pd.DatetimeIndex(var_time_idx)

    var_timestamp_base_start = pd.date_range(start=var_time_idx[0], periods=2, freq=var_time_freq)[1]
    var_timestamp_base_end = var_time_idx[-1]
    var_timestamp_base_period = pd.date_range(start=var_timestamp_base_start, end=var_timestamp_base_end,
                                              freq=var_time_freq)

    var_time_ext_start = pd.date_range(var_timestamp_base_period[-1], periods=2, freq=var_time_freq)[1]
    var_timestamp_ext_period = pd.date_range(start=var_time_ext_start, freq=var_time_freq, periods=var_time_period)

    var_time_period = var_timestamp_base_period.union(var_timestamp_ext_period)

    return var_time_period
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# method to configure expert forecast data
def configure_data_expert_forecast(
        var_dataset, time_dataset,
        time_tag_base_start='time_start', time_tag_extended_start='start',
        data_tag_base_start='rain_start_run',
        info_tag_time_base_period='time_base_period',
        info_tag_time_base_start='time_base_start', info_tag_time_base_end='time_base_end',
        info_tag_time_base_length='time_base_len',
        info_tag_time_extended_period='time_extended_period',
        info_tag_time_extended_start='time_extended_start', info_tag_time_extended_end='time_extended_end',
        info_tag_time_extended_length='time_extended_len',
        info_tag_time_run_period='time_run_period',
        info_tag_time_run_start='time_run_start', info_tag_time_run_end='time_run_end',
        info_tag_time_run_length='time_run_len',
        info_tag_time_split='time_split', info_tag_data_base='rain_base_period',
        time_freq='H', time_period=12):

    # iterate over object(s)
    time_collections, data_collections, info_collections = {}, {}, {}
    for (time_key, time_obj), (data_key, data_obj) in zip(time_dataset.items(), var_dataset.items()):

        # time information
        time_base = time_obj[time_tag_base_start]
        time_forecast = time_obj[time_tag_extended_start]

        if time_tag_base_start == time_tag_extended_start:
            if time_base.__len__() > 1 and time_forecast.__len__() > 1:
                time_start_str = time_base[0]
                time_end_str = time_base[-1]
                time_split = False
            else:
                log_stream.error(' ===> Variables "time_base" and "time_extended" '
                                 'are the same and must be greater than 1')
                raise NotImplemented('Case not implemented yet')
        elif time_tag_base_start != time_tag_extended_start:
            if time_base.__len__() == 1 and time_forecast.__len__() == 1:
                time_start_str = time_base[0]
                time_end_str = time_forecast[0]
                time_split = True
            else:
                log_stream.error(' ===> Variables "time_base" and "time_forecast" '
                                 'are not the same and must be equal to 1')
                raise NotImplemented('Case not implemented yet')
        else:
            log_stream.error(' ===> Variables "time_base" and "time_forecast" are in unsupported format')
            raise NotImplemented('Case not implemented yet')

        # data information
        if time_tag_base_start != time_tag_extended_start:
            if data_tag_base_start in list(data_obj.keys()):
                rain_base = data_obj[data_tag_base_start][0]
            else:
                log_stream.error(' ===> Variables "rain_base" must be defined in the "data_obj"')
                raise NotImplemented('Case not implemented yet')
        else:
            rain_base = None

        # get time string start and end
        time_start_stamp = pd.Timestamp(time_start_str)
        time_end_stamp = pd.Timestamp(time_end_str)

        # compute time period base
        time_base_start = pd.date_range(start=time_start_stamp, periods=2, freq=time_freq)[1]
        time_base_end = deepcopy(time_end_stamp)
        time_base_period = pd.date_range(start=time_base_start, end=time_base_end, freq=time_freq)
        time_base_len = time_base_period.shape[0]
        # compute time period extended (add the last period of forecast)
        time_extended_start = pd.date_range(start=time_base_end, periods=2, freq=time_freq)[1]
        time_extended_period = pd.date_range(start=time_extended_start, freq=time_freq, periods=time_period)
        time_extended_end = time_extended_period[-1]
        time_extended_len = time_extended_period.shape[0]

        # merge time base and extended periods
        time_run_period = time_base_period.union(time_extended_period)
        time_run_start = time_run_period[0]
        time_run_end = time_run_period[-1]
        time_run_len = time_run_period.shape[0]

        # save time information in a common time workspace
        time_collections[time_key] = time_run_period

        # update data information
        if time_tag_base_start in list(data_obj.keys()):
            data_obj[time_tag_base_start] = [time_base_start]
        if time_tag_extended_start in list(data_obj.keys()):
            data_obj[time_tag_extended_start] = [time_extended_start]

        # save data information in a data workspace
        data_collections[time_key] = data_obj

        info_collections[time_key] = {}
        info_collections[time_key][info_tag_time_base_period] = time_base_period
        info_collections[time_key][info_tag_time_base_start] = time_base_start
        info_collections[time_key][info_tag_time_base_end] = time_base_end
        info_collections[time_key][info_tag_time_base_length] = time_base_len
        info_collections[time_key][info_tag_time_extended_period] = time_extended_period
        info_collections[time_key][info_tag_time_extended_start] = time_extended_start
        info_collections[time_key][info_tag_time_extended_end] = time_extended_end
        info_collections[time_key][info_tag_time_extended_length] = time_extended_len
        info_collections[time_key][info_tag_time_run_period] = time_run_period
        info_collections[time_key][info_tag_time_run_start] = time_run_start
        info_collections[time_key][info_tag_time_run_end] = time_run_end
        info_collections[time_key][info_tag_time_run_length] = time_run_len
        info_collections[time_key][info_tag_time_split] = time_split
        info_collections[time_key][info_tag_data_base] = rain_base

    return data_collections, time_collections, info_collections
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to read expert forecast datasets
def read_data_expert_forecast(
        file_list,
        time_tag_start_run_in=None, time_tag_start_forecast_in=None,
        var_tag_name_list_in=None, var_tag_reference_in=None,
        time_tag_start_run_out=None, time_tag_start_forecast_out=None,
        var_tag_name_list_out=None, var_tag_reference_out=None):

    if time_tag_start_run_in is None:
        time_tag_start_run_in = ['time_start']
    if time_tag_start_forecast_in is None:
        time_tag_start_forecast_in = ['time']

    if var_tag_reference_in is None:
        var_tag_reference_in = ['name']
    if var_tag_name_list_in is None:
        var_tag_name_list_in = ['time_start', 'time',
                                'rain_int_start', 'rain_average', 'rain_peak', 'slope_x', 'slope_y', 'slope_t', 'pth']

    if time_tag_start_run_out is None:
        time_tag_start_run_out = ['time_start_run']
    if time_tag_start_forecast_out is None:
        time_tag_start_forecast_out = ['time_start_forecast']

    if var_tag_reference_out is None:
        var_tag_reference_out = ['name']
    if var_tag_name_list_out is None:
        var_tag_name_list_out = ['time_start_run', 'time_start_forecast',
                                 'rain_start_run', 'rain_average', 'rain_peak', 'slope_x', 'slope_y', 'slope_t', 'pth']

    if not isinstance(time_tag_start_run_in, list):
        time_tag_start_run_in = [time_tag_start_run_in]
    if not isinstance(time_tag_start_forecast_in, list):
        time_tag_start_forecast_in = [time_tag_start_forecast_in]
    if not isinstance(var_tag_reference_in, list):
        var_tag_reference_in = [var_tag_reference_in]

    if not isinstance(time_tag_start_run_out, list):
        time_tag_start_run_out = [time_tag_start_run_out]
    if not isinstance(time_tag_start_forecast_out, list):
        time_tag_start_forecast_out = [time_tag_start_forecast_out]
    if not isinstance(var_tag_reference_out, list):
        var_tag_reference_out = [var_tag_reference_out]

    tag_var_fields_in = time_tag_start_run_in + time_tag_start_forecast_in + \
        var_tag_reference_in + var_tag_name_list_in
    tag_var_fields_out = time_tag_start_run_out + time_tag_start_forecast_out + \
        var_tag_reference_out + var_tag_name_list_out

    if not isinstance(var_tag_name_list_in, list):
        var_tag_name_list_in = [var_tag_name_list_in]
    if not isinstance(var_tag_name_list_out, list):
        var_tag_name_list_out = [var_tag_name_list_out]

    if not isinstance(file_list, list):
        file_list = [file_list]

    file_time, file_reference, file_datasets = {}, {}, {}
    for file_id, file_step in enumerate(file_list):

        file_time[file_id], file_reference[file_id], file_datasets[file_id] = {}, {}, {}
        if os.path.exists(file_step):
            file_data = pd.read_json(file_step)
            for var_name_in, var_name_out in zip(tag_var_fields_in, tag_var_fields_out):
                var_data = file_data[var_name_in].values
                if var_data.shape.__len__() == 1:
                    var_list_tmp = var_data.tolist()
                elif var_data.shape.__len__() > 1:
                    var_data = var_data.tolist()
                    if var_data.shape[1] == 1:
                        var_list_tmp = [var_value[0] for var_value in var_data]
                    else:
                        log_stream.error(' ===> Variable "var_data" is defined by unsupported dimensions')
                        raise NotImplementedError('Dataset dimensions are not allowed')
                else:
                    log_stream.error(' ===> Variable "var_data" is defined by unsupported format')
                    raise NotImplementedError('Variable format is not allowed')

                if var_name_in in var_tag_reference_in:
                    file_reference[file_id] = list(set(var_list_tmp))[0]

                if var_name_in in var_tag_name_list_in:
                    file_datasets[file_id][var_name_out] = var_list_tmp
                if var_name_in in time_tag_start_run_in:
                    file_time[file_id][var_name_out] = var_list_tmp
                if var_name_in in time_tag_start_forecast_in:
                    file_time[file_id][var_name_out] = var_list_tmp

    time_collections = {}
    for (file_ref_key, file_ref_name), file_time_fields in zip(file_reference.items(), file_time.values()):
        time_collections[file_ref_name] = file_time_fields

    datasets_collections = {}
    for (file_ref_key, file_ref_name), file_datasets_fields in zip(file_reference.items(), file_datasets.values()):
        datasets_collections[file_ref_name] = file_datasets_fields

    return datasets_collections, time_collections

# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to read expert forecast datasets
def read_data_expert_forecast_OLD(file_list, tag_var_time=None, tag_var_name_list=None, tag_var_reference=None):

    if tag_var_time is None:
        tag_var_time = ['time']

    if tag_var_reference is None:
        tag_var_reference = ['name']
    if tag_var_name_list is None:
        tag_var_name_list = ['rain_average', 'rain_peak', 'slope_x', 'slope_y', 'slope_t']

    tag_var_fields = tag_var_time + tag_var_reference + tag_var_name_list

    if not isinstance(tag_var_name_list, list):
        tag_var_name_list = [tag_var_name_list]
    if not isinstance(file_list, list):
        file_list = [file_list]

    file_reference = {}
    file_datasets = {}
    file_time = None
    for file_id, file_step in enumerate(file_list):
        file_reference[file_id] = {}
        file_datasets[file_id] = {}
        if os.path.exists(file_step):
            file_data = pd.read_json(file_step)
            for var_name in tag_var_fields:
                var_data = file_data[var_name].values
                if var_data.shape.__len__() == 1:
                    var_list_tmp = var_data.tolist()
                elif var_data.shape.__len__() > 1:
                    var_data = var_data.tolist()
                    if var_data.shape[1] == 1:
                        var_list_tmp = [var_value[0] for var_value in var_data]
                    else:
                        raise NotImplementedError('Dataset dimensions not allowed')
                else:
                    raise NotImplementedError('Variable dimensions not allowed')

                if var_name in tag_var_time:
                    if file_time is None:
                        file_time = var_list_tmp
                if var_name in tag_var_reference:
                    file_reference[file_id] = list(set(var_list_tmp))[0]
                if var_name in tag_var_name_list:
                    file_datasets[file_id][var_name] = var_list_tmp

    time_collections = file_time
    datasets_collections = {}
    for (file_ref_key, file_ref_name), file_datasets_fields in zip(file_reference.items(), file_datasets.values()):
        datasets_collections[file_ref_name] = file_datasets_fields

    return datasets_collections, time_collections

# -------------------------------------------------------------------------------------


