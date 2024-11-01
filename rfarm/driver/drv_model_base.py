"""
Library Features:

Name:          drv_model_base
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20190903'
Version:       '1.0.0'
"""
# -------------------------------------------------------------------------------------
# Library
import logging
import os
import re

import numpy as np
import pandas as pd
import xarray as xr

from copy import deepcopy

from rfarm.utils.lib_utils_generic import fill_tags2string
from rfarm.settings.lib_args import logger_name, time_units, time_format, time_calendar

from rfarm.driver.drv_model_io import RFarmData, RFarmResult
from rfarm.driver.drv_model_exec import RFarmModel

# Logging
log_stream = logging.getLogger(logger_name)

# Debug
# import matplotlib.pylab as plt
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Class data object
class DataObj(dict):
    pass
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Class to compute time model
class ModelTime:

    # -------------------------------------------------------------------------------------
    # Method to initialize class
    def __init__(self, **kwargs):
        self.time_step = kwargs['time_step']
        self.time_run = kwargs['time_run']
        self.time_settings = kwargs['time_settings']
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to compute model time
    def computeModelTime(self, time_ascending=True, time_closed='right'):

        log_stream.info(' ---> Compute model time ... ')

        time_step = self.time_step
        time_run = self.time_run
        time_settings = self.time_settings

        if 'time_observed_period' in time_settings and 'time_observed_frequency' in time_settings:
            time_observed_period = time_settings['time_observed_period']
            time_observed_frequency = time_settings['time_observed_frequency']
        else:
            time_observed_period = 0
            time_observed_frequency = 'H'

        if 'time_forecast_period' in time_settings and 'time_forecast_frequency' in time_settings:
            time_forecast_period = time_settings['time_forecast_period']
            time_forecast_frequency = time_settings['time_forecast_frequency']
        else:
            time_forecast_period = 0
            time_forecast_frequency = 'H'

        if time_observed_frequency == 'A-OFFSET':
            time_observed_frequency = pd.DateOffset(years=1)

        time_observed_range = pd.date_range(end=time_step, periods=time_observed_period + 1,
                                            freq=time_observed_frequency)
        time_forecast_range = pd.date_range(start=time_step, periods=time_forecast_period + 1,
                                            freq=time_forecast_frequency)

        if not time_observed_range.empty:
            time_start = time_observed_range[0]
        else:
            time_start = time_step
        if not time_forecast_range.empty:
            time_end = time_forecast_range[-1]
        else:
            time_end = time_step
        try:
            time_range = pd.date_range(start=time_start, end=time_end,
                                   #freq=time_observed_frequency, inclusive=time_closed)
                                   freq=time_observed_frequency, closed=time_closed)
        except:
            me_range = pd.date_range(start=time_start, end=time_end,
                                    freq=time_observed_frequency, inclusive=time_closed)

        time_range = time_range.sort_values(return_indexer=False, ascending=time_ascending)

        time_obj = DataObj
        time_obj.time_run = time_run
        time_obj.time_step = time_step
        time_obj.time_range = time_range
        time_obj.time_from = time_range[0]
        time_obj.time_to = time_range[-1]
        time_obj.time_format = time_format
        time_obj.time_units = time_units
        time_obj.time_calendar = time_calendar

        log_stream.info(' ---> Compute model time ... DONE!')

        return time_obj

# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Class to driver model
class ModelRunner:

    # -------------------------------------------------------------------------------------
    # Method to initialize class
    def __init__(self, time_step, time_range,
                 variable_in, variable_out, data_geo, template, parameters,
                 file_dict_ancillary_in=None, file_dict_ancillary_out=None,
                 file_dict_in=None, file_dict_out=None,
                 file_ancillary_in_updating=True, file_ancillary_out_updating=None,
                 file_out_updating=None, file_out_zipping=False, file_ext_zipping='.gz',
                 file_write_engine='netcdf4', file_domain_name='regional_domain',
                 file_reference_dim='time', file_reference_step=None,
                 tag_folder_name='folder', tag_file_name='filename',
                 tag_terrain_data='terrain_data', tag_alert_area_data='alert_area_data',
                 tag_model_algorithm='exec_nwp'):

        # Generic information
        self.time_step = time_step
        self.time_range = time_range
        self.var_info_in = variable_in
        self.var_info_out = variable_out
        self.data_geo = data_geo
        self.model_tags_template = template
        self.model_parameters = parameters

        self.file_dict_in = file_dict_in
        self.file_dict_ancillary_in = file_dict_ancillary_in
        self.file_dict_ancillary_out = file_dict_ancillary_out
        self.file_dict_out = file_dict_out

        self.tag_folder_name = tag_folder_name
        self.tag_file_name = tag_file_name
        self.tag_terrain_data = tag_terrain_data
        self.tag_alert_area_data = tag_alert_area_data

        self.tag_model_algorithm = tag_model_algorithm

        # add parameters not defined in the configuration file
        if self.tag_model_algorithm == 'exec_expert_forecast':
            if 'rain_max_thr' not in list(self.model_parameters.keys()):
                log_stream.warning(' ===> Parameter "rain_max_thr" is not defined in the configuration file')
                log_stream.warning(' ===> Parameter "rain_max_thr" set equal to 150 [mm]')
                self.model_parameters['rain_max_thr'] = float(150)
            else:
                self.model_parameters['rain_max_thr'] = float(self.model_parameters['rain_max_thr'])
        else:
            self.model_parameters['rain_max_thr'] = None

        folder_name_in = self.file_dict_in[self.tag_folder_name]
        file_name_in = self.file_dict_in[self.tag_file_name]
        file_in = os.path.join(folder_name_in, file_name_in)

        folder_ancillary_in = self.file_dict_ancillary_in[self.tag_folder_name]
        file_ancillary_in = self.file_dict_ancillary_in[self.tag_file_name]
        file_ancillary_in = os.path.join(folder_ancillary_in, file_ancillary_in)

        folder_ancillary_out = self.file_dict_ancillary_out[self.tag_folder_name]
        file_ancillary_out = self.file_dict_ancillary_out[self.tag_file_name]
        file_ancillary_out = os.path.join(folder_ancillary_out, file_ancillary_out)

        folder_name_out = self.file_dict_out[self.tag_folder_name]
        file_name_out = self.file_dict_out[self.tag_file_name]
        file_out = os.path.join(folder_name_out, file_name_out)

        self.data_geo_terrain = self.data_geo[self.tag_terrain_data]
        self.data_geo_alert_area = self.data_geo[self.tag_alert_area_data]

        if isinstance(file_domain_name, str):
            self.var_domain_name = [file_domain_name]
            self.var_domain_id = [1]
        elif isinstance(file_domain_name, dict):
            var_domain_name = []
            var_domain_id = []
            for domain_key, domain_fields in file_domain_name.items():
                domain_name = domain_fields['name']
                domain_id = domain_fields['id']
                var_domain_name.append(domain_name)
                var_domain_id.append(domain_id)
            self.var_domain_name = var_domain_name
            self.var_domain_id = var_domain_id
        else:
            log_stream.error(' ===> Domain name format not permitted')
            raise NotImplementedError('Case not implemented yet')

        self.file_reference_dim = file_reference_dim
        self.file_reference_step = file_reference_step

        # Data input file
        if file_in is not None:

            if self.file_reference_dim == 'time':
                model_tags_in_values = {'datetime_input': self.time_step, 'sub_path_time': self.time_step,
                                        'timestep_input': 1}

                file_in_raw = fill_tags2string(file_in, self.model_tags_template, model_tags_in_values)
                self.folder_in_raw, self.filename_in_raw = os.path.split(file_in_raw)
            elif self.file_reference_dim == 'domain':
                model_tags_in_values = {'datetime_input': self.time_step, 'sub_path_time': self.time_step,
                                        'domain': self.var_domain_name[0], 'timestep_input': 1}

                file_in_raw = fill_tags2string(file_in, self.model_tags_template, model_tags_in_values)
                self.folder_in_raw, self.filename_in_raw = os.path.split(file_in_raw)
            else:
                log_stream.error(' ===> File dimensions type is not allowed')
                raise NotImplementedError('Case not implemented yet')

            if self.file_reference_dim == 'time':
                self.folder_in_list, self.filename_in_list = [], []
                reference_id = None
                for time_id, time_step in enumerate(self.time_range):

                    if reference_id is None:
                        if time_step.hour == 1:
                            reference_id = 1
                        elif time_step.hour == 0:
                            reference_id = 0
                        else:
                            log_stream.error(' ===> Reference step "' + str(time_iter.hour) +
                                             '" is not supported')
                            raise NotImplementedError('Case not implemented yet')

                    if self.file_reference_step == 'first':
                        reference_time = self.time_range[0]
                    elif self.file_reference_step is None:
                        reference_time = deepcopy(time_step)
                    else:
                        log_stream.error(' ===> Reference time  "' + str(self.file_reference_step) +
                                         '" for "file_reference_dim == time" is not supported')
                        raise NotImplementedError('Case not implemented yet')

                    model_tags_list_values = {'datetime_input': reference_time, 'sub_path_time': self.time_step,
                                              'timestep_input': time_id + reference_id}
                    file_in_list = fill_tags2string(file_in, self.model_tags_template, model_tags_list_values)
                    folder_in_list, filename_in_list = os.path.split(file_in_list)
                    self.folder_in_list.append(folder_in_list)
                    self.filename_in_list.append(filename_in_list)

            elif self.file_reference_dim == 'domain':

                if self.file_reference_step is not None:
                    log_stream.error(' ===> Reference time  "' + str(self.file_reference_step) +
                                     '" for "file_reference_dim == domain" is not supported')
                    raise NotImplementedError('Case not implemented yet')

                self.folder_in_list, self.filename_in_list = [], []
                for domain_iter in self.var_domain_name:
                    model_tags_list_values = {'datetime_input': self.time_step, 'sub_path_time': self.time_step,
                                              'domain': domain_iter}
                    file_in_list = fill_tags2string(file_in, self.model_tags_template, model_tags_list_values)
                    folder_in_list, filename_in_list = os.path.split(file_in_list)
                    self.folder_in_list.append(folder_in_list)
                    self.filename_in_list.append(filename_in_list)

            else:
                log_stream.error(' ===> File dimensions type is not allowed')
                raise NotImplementedError('Case not implemented yet')
        else:
            log_stream.error(' ===> File input is defined by NoneType')
            raise TypeError('File input are not correctly defined')

        if not os.path.exists(self.folder_in_raw):
            log_stream.error(' ===> Input folder "' + self.folder_in_raw + '" does not exist')
            raise FileNotFoundError('Check your settings to search the right location.')

        # Ancillary input file
        if file_ancillary_in is not None:
            model_tags_values = {'datetime_input': self.time_step, 'sub_path_time': self.time_step}
            file_ancillary_in = fill_tags2string(file_ancillary_in, self.model_tags_template, model_tags_values)
            self.folder_ancillary_in_raw, self.filename_ancillary_in_raw = os.path.split(file_ancillary_in)
        else:
            self.folder_ancillary_in_raw = None
            self.filename_ancillary_in_raw = None

        if self.folder_ancillary_in_raw is not None:
            folder_ancillary_in_root = re.split('[{}]', self.folder_ancillary_in_raw)[0]
            if not os.path.exists(folder_ancillary_in_root):
                os.makedirs(folder_ancillary_in_root)

        # Ancillary outcome file
        if file_ancillary_out is not None:
            model_tags_template = deepcopy(self.model_tags_template)
            model_tags_template.pop('ensemble', None)
            model_tags_values = {'datetime_outcome': self.time_step, 'sub_path_time': self.time_step}
            file_ancillary_out = fill_tags2string(file_ancillary_out, model_tags_template, model_tags_values)
            self.folder_ancillary_out_raw, self.filename_ancillary_out_raw = os.path.split(file_ancillary_out)
        else:
            self.folder_ancillary_out_raw = None
            self.filename_ancillary_out_raw = None

        if self.folder_ancillary_out_raw is not None:
            folder_ancillary_out_root = re.split('[{}]', self.folder_ancillary_out_raw)[0]
            if not os.path.exists(folder_ancillary_out_root):
                os.makedirs(folder_ancillary_out_root)

        # Result outcome file
        if file_out is not None:
            model_tags_template = deepcopy(self.model_tags_template)
            model_tags_template.pop('ensemble', None)
            model_tags_values = {'datetime_outcome': self.time_step, 'sub_path_time': self.time_step}
            file_out = fill_tags2string(file_out, model_tags_template, model_tags_values)
            self.folder_out_raw, self.filename_out_raw = os.path.split(file_out)
        else:
            raise TypeError

        folder_out_root = re.split('[{}]', self.folder_out_raw)[0]
        if not os.path.exists(folder_out_root):
            os.makedirs(folder_out_root)

        # Model driver for collecting RFarm data
        self.driver_rf_data = RFarmData(
            self.folder_in_raw, self.filename_in_raw,
            folder_data_list=self.folder_in_list, filename_data_list=self.filename_in_list,
            var_name=self.var_info_in['id']['var_name'],
            var_dims=self.var_info_in['id']['var_type'][0],
            var_type=self.var_info_in['id']['var_type'][1],
            var_units=self.var_info_in['attributes']['units'],
            var_domain=self.var_domain_name,
            file_format=self.var_info_in['id']['var_format'],
            file_source=self.var_info_in['id']['var_source'],
            folder_tmp=self.folder_ancillary_in_raw, filename_tmp=self.filename_ancillary_in_raw,
        )

        # Model driver for running RFarm model
        self.driver_rf_model = RFarmModel(
            ensemble_n=self.model_parameters['ensemble'],
            ensemble_format=self.model_tags_template['ensemble'],
            ratio_s=self.model_parameters['ratio_s'],
            ratio_t=self.model_parameters['ratio_t'],
            slope_s=self.model_parameters['slope_s'],
            slope_t=self.model_parameters['slope_t'],
            cs_sf=self.model_parameters['cs_sf'],
            ct_sf=self.model_parameters['ct_sf'],
            multi_core=self.model_parameters['multi_core'],
            domain_extension=self.model_parameters['domain_extension'],
            folder_tmp=self.folder_ancillary_out_raw,
            filename_tmp=self.filename_ancillary_out_raw,
            model_var=self.var_info_out['id']['var_name'],
            model_algorithm=self.tag_model_algorithm,
            rain_max_thr=self.model_parameters['rain_max_thr'],
        )

        if 'var_frequency' in self.var_info_out['id']:
            var_frequency = self.var_info_out['id']['var_frequency']
        else:
            log_stream.warning(' ===> Result variable frequency is not set. Use default hourly frequency "H"')
            var_frequency = 'H'

        # Model driver of saving RFarm result(s)
        self.driver_rf_result = RFarmResult(
            self.driver_rf_model.ensemble_filename,
            ensemble_n=self.model_parameters['ensemble'],
            ensemble_format=self.model_tags_template['ensemble'],
            folder_out=self.folder_out_raw,
            filename_out=self.filename_out_raw,
            var_name=self.var_info_out['id']['var_name'],
            var_freq=var_frequency,
            var_dims=self.var_info_out['id']['var_type'][0],
            var_attrs=self.var_info_out['attributes'],
            file_format=self.var_info_out['id']['var_format'],
            ensemble_zip=file_out_zipping,
            ext_zip_type=file_ext_zipping,
            write_engine=file_write_engine
        )

        self.file_ancillary_in_updating = file_ancillary_in_updating
        self.file_ancillary_out_updating = file_ancillary_out_updating
        self.file_out_updating = file_out_updating
        self.file_out_zipping = file_out_zipping
        self.file_write_engine = file_write_engine

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to clean temporary data
    def clean(self, cleaning_dynamic_tmp=False):

        if cleaning_dynamic_tmp:

            # Clean ensemble group tmp input file(s)
            file_group_in_tmp = self.driver_rf_data.file_tmp
            if isinstance(file_group_in_tmp, str):
                file_group_in_tmp = [file_group_in_tmp]
            for file in file_group_in_tmp:
                if os.path.exists(file):
                    os.remove(file)

            # Clean ensemble group tmp outcome file(s)
            file_group_out_tmp = self.driver_rf_result.ensemble_group_in
            if isinstance(file_group_out_tmp, str):
                file_group_out_tmp = [file_group_out_tmp]
            for file in file_group_out_tmp:
                if os.path.exists(file):
                    os.remove(file)

        if self.file_out_zipping:
            file_group_out_unzip = self.driver_rf_result.ensemble_group_out
            if isinstance(file_group_out_unzip, str):
                file_group_out_unzip = [file_group_out_unzip]
            # Clean ensemble group outcome file(s)
            for file in file_group_out_unzip:
                if os.path.exists(file):
                    os.remove(file)
        else:
            file_group_out_zip = self.driver_rf_result.zip_group_out
            if isinstance(file_group_out_zip, str):
                file_group_out_zip = [file_group_out_zip]
            # Clean ensemble group outcome file(s)
            for file in file_group_out_zip:
                if os.path.exists(file):
                    os.remove(file)

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to collect data
    def collect(self, var_time_obj):

        # -------------------------------------------------------------------------------------
        # Starting info
        log_stream.info(' ---> Collect data ... ')

        # Get time steps
        var_time_range = var_time_obj.time_range

        # Check data availability
        if self.driver_rf_data.found_data_first:

            # Configure updating flag about ancillary in file
            if os.path.exists(self.driver_rf_data.file_tmp):
                if self.file_ancillary_in_updating:
                    os.remove(self.driver_rf_data.file_tmp)
            else:
                self.file_ancillary_in_updating = True

            if self.file_ancillary_in_updating:
                data_obj = self.driver_rf_data.organize_data(var_time_range)
                # Ending info
                log_stream.info(' ---> Collect data ... DONE')
            else:
                data_obj = self.driver_rf_data.callback_data()
                # Ending info
                log_stream.info(' ---> Collect data ... PREVIOUSLY DONE!')

        else:
            # Ending info
            log_stream.info(' ---> Collect data ... SKIPPED! DATA INPUT NOT FOUND!')
            data_obj = None

        return data_obj

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to configure model
    def configure(self, data_obj):

        # -------------------------------------------------------------------------------------
        # info start
        log_stream.info(' ---> Model configuration ... ')

        # Check data availability
        if self.driver_rf_data.found_data_first:

            # Method to set ancillary out flag
            self.file_ancillary_out_updating = self.__set_ancillary_out_flag(self.file_ancillary_out_updating)

            # Check ancillary outcome flag
            if self.file_ancillary_out_updating:

                if isinstance(data_obj, dict):

                    values_obj = data_obj['data']
                    lons_obj = self.data_geo_terrain['longitude']
                    lats_obj = self.data_geo_terrain['latitude']
                    time_obj = data_obj['time']
                    info_obj = data_obj['info']

                elif isinstance(data_obj, xr.DataArray):

                    values_obj = data_obj.values
                    lons_obj = data_obj['longitude'].values
                    lats_obj = data_obj['latitude'].values
                    time_obj = data_obj['time'].values
                    info_obj = None

                else:
                    log_stream.error(' ===> Data object not in supported format.')
                    raise NotImplementedError('Case not implemented yet')

                lons_geo = self.data_geo_terrain['longitude']
                lats_geo = self.data_geo_terrain['latitude']
                res_lon_geo = self.data_geo_terrain['res_lon']
                res_lat_geo = self.data_geo_terrain['res_lat']

                # configure model grid(s)
                self.driver_rf_model.configure_grid(lons_obj, lats_obj,
                                                    lons_geo, lats_geo, res_lon_geo, res_lat_geo,
                                                    self.driver_rf_model.domain_extension)
                # configure model time(s)
                self.driver_rf_model.configure_time(time_obj)

                # configure model info
                self.driver_rf_model.configure_info(info_obj)

                # configure model data
                self.driver_rf_model.configure_data(values_obj, info_obj)

                # info end
                log_stream.info(' ---> Model configuration ... DONE!')

            else:
                # info end
                log_stream.info(' ---> Model configuration ... SKIPPED! Results are saved during previously executions!')

        else:
            # info end
            log_stream.info(' ---> Model configuration ... SKIPPED! DATA INPUT NOT FOUND!')
        # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to execute model
    def exec(self):

        # -------------------------------------------------------------------------------------
        # Starting info
        log_stream.info(' ---> Model execution ... ')

        # Get domain information
        domain_id = self.var_domain_id
        domain_name = self.var_domain_name

        # Get alert area domain(s)
        if self.data_geo_alert_area is not None:
            domain_mask = self.data_geo_alert_area['values']
        else:
            terrain_data = self.data_geo_terrain['values']
            domain_mask = np.zeros(shape=[terrain_data.shape[0], terrain_data.shape[1]])
            domain_mask[:, :] = 1
            domain_mask[terrain_data < 0] = -9999

        # Check data availability
        if self.driver_rf_data.found_data_first:

            # Check ancillary outcome flag
            if self.file_ancillary_out_updating:
                # Method to execute model run(s)
                ensemble_obj = self.driver_rf_model.execute_run(domain_name=domain_name, domain_id=domain_id,
                                                                domain_mask=domain_mask)
                # Ending info
                log_stream.info(' ---> Model execution ... DONE!')
            else:
                # Method to collect model run(s)
                ensemble_obj = self.driver_rf_model.callback_run()
                # Ending info
                log_stream.info(' ---> Model execution ... SKIPPED! Results are collected using previously executions!')

        else:
            # Ending info
            log_stream.info(' ---> Model execution ... SKIPPED! DATA INPUT NOT FOUND!')
            ensemble_obj = None

        return ensemble_obj
        # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to save model
    def save(self, ensemble_obj):

        # -------------------------------------------------------------------------------------
        # Starting info
        log_stream.info(' ---> Save result(s) ... ')

        # Check data availability
        if self.driver_rf_data.found_data_first:

            # Geographical info
            values_geo = self.data_geo_terrain['values']
            lons_geo = self.data_geo_terrain['longitude']
            lats_geo = self.data_geo_terrain['latitude']

            # Method to set ancillary out flag
            self.file_out_updating = self.__set_ancillary_out_flag(self.file_out_updating)

            # Check ancillary outcome flag
            if self.file_out_updating:
                # Method to save model result(s)
                self.driver_rf_result.organize_result(ensemble_obj, values_geo, lons_geo, lats_geo)
                # Method to zip model result(s)
                self.driver_rf_result.zip_result()
                # Ending info
                log_stream.info(' ---> Save result(s) ... DONE!')
            else:
                # Ending info
                log_stream.info(' ---> Save result(s) ... SKIPPED! Results are collected using previously executions!')

        else:
            # Ending info
            log_stream.info(' ---> Save result(s) ... SKIPPED! DATA INPUT NOT FOUND!')

        # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to set ancillary flag
    def __set_ancillary_out_flag(self, file_out_updating):

        if self.file_out_zipping:
            for ensemble_filename in self.driver_rf_result.ensemble_group_out:
                if os.path.exists(ensemble_filename):
                    os.remove(ensemble_filename)
            for ensemble_filename_zip in self.driver_rf_result.zip_group_out:
                if os.path.exists(ensemble_filename_zip):
                    if file_out_updating:
                        os.remove(ensemble_filename_zip)
                else:
                    file_out_updating = True
                    break
        else:
            for ensemble_filename_zip in self.driver_rf_result.zip_group_out:
                if os.path.exists(ensemble_filename_zip):
                    os.remove(ensemble_filename_zip)
            for ensemble_filename in self.driver_rf_result.ensemble_group_out:
                if os.path.exists(ensemble_filename):
                    if file_out_updating:
                        os.remove(ensemble_filename)
                else:
                    file_out_updating = True
                    break
        return file_out_updating
    # -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
