"""
Library Features:

Name:          drv_model_io
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20241031'
Version:       '1.2.0'
"""
#################################################################################
# Library
import logging
import os
import re

from copy import deepcopy

import numpy as np
import pandas as pd

from rfarm.settings.lib_conventions import oVarConventions as var_def_conventions  # FX DA RIVEDERE
from rfarm.io.lib_io_generic_gzip import zip_filename

from rfarm.settings.lib_args import logger_name

from rfarm.dataset.lib_variables import compute_rain_lami_2i, compute_rain_icon_2i, compute_rain_ecmwf_0100

from rfarm.io.lib_io_generic_fx import create_darray_3d, write_obj, read_obj
from rfarm.io.lib_io_dst_netcdf import create_dset_nc, write_dset_nc
from rfarm.io.lib_io_dst_binary import write_dset_binary
from rfarm.io.lib_io_src_grib_lami import read_data_lami_2i, adjust_data_lami_2i
from rfarm.io.lib_io_src_grib_icon import read_data_icon_2i, adjust_data_icon_2i
from rfarm.io.lib_io_src_grib_ecmwf import read_data_ecmwf_0100, adjust_data_ecmwf_0100
from rfarm.io.lib_io_src_netcdf import read_data_wrf, read_data_gfs_025
from rfarm.io.lib_io_src_netcdf import convert_data_wrf, convert_time_wrf
from rfarm.io.lib_io_src_tiff import read_data_moloc

from rfarm.io.lib_io_src_json import read_data_expert_forecast, configure_data_expert_forecast

# from rfarm.core.lib_core_generic import plotResult
from rfarm.core.lib_core_generic import compute_ensemble

from rfarm.utils.lib_utils_generic import fill_tags2string
from rfarm.utils.lib_utils_generic import merge_dict    # FX DA RIVEDERE
from rfarm.utils.lib_utils_system import create_tmp     # FX DA RIVEDERE

# Logging
log_stream = logging.getLogger(logger_name)
# Debug
# import matplotlib.pylab as plt
#################################################################################


# -------------------------------------------------------------------------------------
# Class data object
class DataObj(dict):
    pass
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Class to organize result(s)
class RFarmResult:

    def __init__(self, ensemble_group_in,
                 ensemble_n=None, ensemble_format='{:03d}',
                 ensemble_zip=False, ext_zip_type='.gz',
                 folder_out=None, filename_out='rfarm_{ensemble}.nc',
                 var_name='Rain', var_dims='var3d', var_attrs=None, var_freq='H',
                 file_format='netcdf', file_dst=None,
                 dim_x_name='west_east', dim_y_name='south_north', write_engine='netcdf4'):

        if ensemble_n is None:
            ensemble_n = {'start': 1, 'end': 2}

        # define ensemble elements
        self.ensemble = compute_ensemble(ensemble_n['start'], ensemble_n['end'])

        if isinstance(ensemble_group_in, str):
            ensemble_group_in = [ensemble_group_in]

        self.ensemble_group_in = ensemble_group_in
        self.ensemble_n = ensemble_n
        self.ensemble_format = ensemble_format
        self.ensemble_zip = ensemble_zip

        self.var_name = var_name
        self.var_dims = var_dims
        self.var_freq = var_freq
        self.var_attrs_dict = {var_name: var_attrs}

        self.folder_out = folder_out
        self.filename_out = filename_out
        self.ext_zip_type = ext_zip_type

        self.dim_x_name = dim_x_name
        self.dim_y_name = dim_y_name

        self.write_engine = write_engine

        if self.folder_out is None:
            self.folder_out = create_tmp()

        if self.filename_out is None:
            raise TypeError

        self.ensemble_group_out = []
        self.zip_group_out = []
        for ensemble_id in self.ensemble:

            tags_tmpl = {'ensemble': ensemble_format}
            tags_values = {'ensemble':  ensemble_format.format(ensemble_id)}

            folder_out = fill_tags2string(self.folder_out, tags_tmpl, tags_values)
            filename_out = fill_tags2string(self.filename_out, tags_tmpl, tags_values)
            self.ensemble_group_out.append(os.path.join(folder_out, filename_out))

            if not os.path.exists(folder_out):
                os.makedirs(folder_out)

            filezip_out = filename_out + self.ext_zip_type
            self.zip_group_out.append(os.path.join(folder_out, filezip_out))

        self.file_dst = file_dst
        self.file_format = file_format

        if self.file_format == 'netcdf':
            if self.write_engine is None:
                log_stream.error(' ===> The "write_engine" must be defined with "file_format" == "netcdf')
                raise RuntimeError('Impossible to run the algorithm without "write_engine" information')

    def zip_result(self):

        # Zip result file according with zipping flag
        if self.ensemble_zip:
            for ensemble_name, filename_out, filezip_out in zip(self.ensemble, self.ensemble_group_out, self.zip_group_out):

                # Starting info
                log_stream.info(' ----> Zip ensemble ' + str(ensemble_name) + ' ... ')

                if os.path.exists(filezip_out):
                    os.remove(filezip_out)

                if os.path.exists(filename_out):

                    # zip filename to reduce disk usage (if selected)
                    zip_filename(filename_out, filezip_out)
                    # Starting info
                    log_stream.info(' ----> Zip ensemble ' + str(ensemble_name) + ' ... DONE')

                else:
                    # Starting info
                    log_stream.warning(' ----> Zip ensemble ' + str(ensemble_name) +
                                       ' ... FAILED! File unzipped not found!')

                if os.path.exists(filezip_out) and os.path.exists(filename_out):
                    os.remove(filename_out)

    def organize_result(self, ensemble_obj, terrain, lons, lats, terrain_name='Terrain',
                        compression_level=0):

        # Check ensemble expected with ensemble(s) saved
        if isinstance(ensemble_obj, str):
            ensemble_obj = [ensemble_obj]

        for ensemble_found in ensemble_obj:
            if ensemble_found not in self.ensemble_group_in:
                log_stream.warning(' -----> File not found: ' + ensemble_found)

        # Terrain attributes
        terrain_attr_dict = {}
        for var_def_key, var_def_value in var_def_conventions.items():
            if terrain_name.lower() == var_def_key.lower():
                terrain_attr_dict[terrain_name] = var_def_value

        # Var attributes
        var_attrs_dict = self.var_attrs_dict
        var_freq_expected = self.var_freq

        # Iterate over ensemble(s)
        for filename_id, (ensemble_name, filename_in, filename_out) in enumerate(
                zip(self.ensemble, self.ensemble_group_in, self.ensemble_group_out)):

            # Starting info
            log_stream.info(' ----> Dump ensemble ' + str(ensemble_name) + ' ... ')

            # Check file availability
            if os.path.exists(filename_in):

                # Get data
                log_stream.info(' -----> Get dataset object ... ')
                da_tmp = read_obj(filename_in)
                log_stream.info(' -----> Get dataset object ... DONE')

                # Adjust time frequency
                log_stream.info(' -----> Adjust dataset time frequency ... ')
                var_freq_tmp = pd.to_datetime(list(da_tmp.time.values)).inferred_freq
                if isinstance(var_freq_tmp, str):
                    if not re.findall('\d+', var_freq_expected):
                        var_freq_expected = str(1) + var_freq_expected
                    if not re.findall('\d+', var_freq_tmp):
                        var_freq_tmp = str(1) + var_freq_tmp
                else:
                    log_stream.info(' ----> Dump ensemble ' + str(ensemble_name) + ' ... FAILED.')
                    log_stream.error(' ===> Frequency format is not allowed')
                    raise NotImplementedError('Case not implemented yet')

                if pd.to_timedelta(var_freq_expected) > pd.to_timedelta(var_freq_tmp):
                    log_stream.info(' ------> Resample datasets from object frequency "' + var_freq_tmp +
                                    '" to the expected frequency "' + var_freq_expected + '" ... ')
                    da_out = da_tmp.resample(time=var_freq_expected, closed='right', label='right').sum()
                    log_stream.info(' ------> Resample datasets from object frequency "' + var_freq_tmp +
                                    '" to the expected frequency "' + var_freq_expected + '" ... DONE')
                elif pd.to_timedelta(var_freq_expected) == pd.to_timedelta(var_freq_tmp):
                    log_stream.info(' ------> Resample datasets is not activated. ' +
                                    'The object frequency is equal to the expected frequency.')
                    da_out = deepcopy(da_tmp)
                else:
                    log_stream.info(' ----> Dump ensemble ' + str(ensemble_name) + ' ... FAILED.')
                    log_stream.error(' ===> Resampling type is not allowed by the procedure to save result')
                    raise NotImplementedError('Case not implemented yet')

                log_stream.info(' -----> Adjust dataset time frequency ... DONE')

                # Extract values and time
                log_stream.info(' -----> Extract dataset values ... ')
                values_out_raw = da_out.values
                geo_x_out_raw = da_out['longitude'].values
                geo_y_out_raw = da_out['latitude'].values
                time_out = pd.to_datetime(list(da_out.time.values))
                log_stream.info(' -----> Extract dataset values ... DONE')

                # DEBUG START
                # import matplotlib.pylab as plt
                # plt.figure()
                # plt.imshow(values_out_raw[1, :, :])
                # plt.colorbar()
                # plt.clim(0, 10)
                # plt.show()
                # plotResult(values_out_raw[1, :, :], lons, lats)
                # plotResult(values_out_raw[1, :, :], geo_x_out_raw, np.flipud(geo_y_out_raw))
                # DEBUG END

                # Organizing values to save in a correct 3D format
                log_stream.info(' -----> Organize dataset values ... ')
                values_out_def = np.zeros([values_out_raw.shape[0], values_out_raw.shape[1], values_out_raw.shape[2]])
                values_out_def[:, :, :] = np.nan
                for id in range(0, values_out_raw.shape[0]):
                    values_out_def[id, :, :] = np.flipud(values_out_raw[id, :, :])

                    # DEBUG
                    # import matplotlib.pylab as plt
                    # plt.figure(1)
                    # plt.imshow(values_out_raw[id, :, :])
                    # plt.colorbar()
                    # plt.clim(0, 100)
                    # plt.figure(2)
                    # plt.imshow(values_out_def[id, :, :])
                    # plt.colorbar()
                    # plt.clim(0, 100)
                    # plt.show()
                    # DEBUG
                log_stream.info(' -----> Organize dataset values ... DONE')

                # Create and write dset
                log_stream.info(' -----> Save dataset values ... ')
                if self.file_format == 'netcdf':
                    dset_out = create_dset_nc(
                        time_out, values_out_def,
                        np.flipud(terrain), lons, np.flipud(lats),
                        var_name=self.var_name, terrain_name=terrain_name,
                        var_attrs=var_attrs_dict[self.var_name],
                        terrain_attrs=terrain_attr_dict[terrain_name],
                        dim_x_name=self.dim_x_name, dim_y_name=self.dim_y_name)

                    write_dset_nc(
                        filename_out, dset_out, attrs=merge_dict(var_attrs_dict, terrain_attr_dict),
                        compression=compression_level, engine=self.write_engine)

                elif self.file_format == 'binary':
                    write_dset_binary(filename_out, file_data=values_out_def)
                else:
                    log_stream.error(' ===> File format "' + file_format + '" is not supported')
                    raise NotImplemented('Case not implemented yet')

                log_stream.info(' -----> Save dataset values ... DONE')

                # Ending info
                log_stream.info(' ----> Dump ensemble ' + str(ensemble_name) + ' ... DONE')

            else:
                # Ending info
                log_stream.warning(' ----> Dump ensemble ' + str(ensemble_name) + ' ... FAILED! Data not available!')

# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Class to read data
class RFarmData:

    def __init__(self, folder_data_first, filename_data_first,
                 folder_data_list=None, filename_data_list=None,
                 var_dims='var3d', var_type='istantaneous', var_name='Rain', var_units='mm', var_domain=None,
                 file_format=None, file_source=None,
                 folder_tmp=None, filename_tmp='rfarm_data.nc'
                 ):

        self.folder_data_first = folder_data_first
        self.filename_data_first = filename_data_first

        self.var_dims_data = var_dims
        self.var_type_data = var_type
        self.var_name_data = var_name
        self.var_units_data = var_units
        if var_domain is None:
            var_domain = ['domain']
        self.var_domain = var_domain

        self.file_format_data = file_format
        self.file_source_data = file_source

        self.folder_data_list = folder_data_list
        self.filename_data_list = filename_data_list

        if self.file_format_data is None:
            raise TypeError
        if self.file_source_data is None:
            raise TypeError

        self.file_data_first = os.path.join(self.folder_data_first, self.filename_data_first)

        if os.path.exists(self.file_data_first):
            self.found_data_first = True
        else:
            self.found_data_first = False

        self.file_data_list = []
        if self.var_dims_data == 'var1d':
            for folder_data_step, filename_data_step in zip(self.folder_data_list, self.filename_data_list):
                self.file_data_list.append(os.path.join(folder_data_step, filename_data_step))
        elif self.var_dims_data == 'var2d':
            for folder_data_step, filename_data_step in zip(self.folder_data_list, self.filename_data_list):
                self.file_data_list.append(os.path.join(folder_data_step, filename_data_step))
        elif self.var_dims_data == 'var3d':
            self.file_data_list = None
        else:
            raise NotImplementedError('Variable dimensions not allowed.')

        self.folder_tmp = folder_tmp
        self.filename_tmp = filename_tmp

        if folder_tmp is None:
            self.folder_tmp = create_tmp()

        if self.filename_tmp is None:
            raise TypeError

        if not os.path.exists(self.folder_tmp):
            os.makedirs(self.folder_tmp)

        self.file_tmp = os.path.join(self.folder_tmp, self.filename_tmp)

    def callback_data(self):

        var_obj = read_obj(self.file_tmp)
        return var_obj

    def organize_data(self, file_time_data=None):

        # Starting info
        log_stream.info(' ----> Get data ... ')

        # Check file data availability
        var_obj = None
        if self.found_data_first:

            var_data_def, var_time_def, var_geox, var_geoy, var_info_def = None, None, None, None, None
            if self.file_format_data == 'grib':

                if self.file_source_data == 'lami_2i':

                    if self.var_dims_data == 'var2d':
                        log_stream.error(' ----> Get data ... FAILED! FILE SOURCE DIMS NOT ALLOWED!')
                        raise NotImplementedError('NWP lami 2D dimensions datasets not implemented yet')
                    elif self.var_dims_data == 'var3d':

                        [var_da, time_da, geox_obj, geoy_obj] = read_data_lami_2i(
                            self.file_data_first, data_var=self.var_name_data)

                        [var_data_raw,  var_time_def, var_geox, var_geoy] = adjust_data_lami_2i(
                            var_da, time_da, geox_obj, geoy_obj)

                        var_data_def = compute_rain_lami_2i(var_data_raw)

                    else:
                        log_stream.error(' ----> Get data ... FAILED! FILE SOURCE DIMS NOT ALLOWED!')
                        raise NotImplementedError('NWP lami case dimension datasets not implemented yet')

                elif self.file_source_data == 'icon_2i':

                    if self.var_dims_data == 'var2d':
                        log_stream.error(' ----> Get data ... FAILED! FILE SOURCE DIMS NOT ALLOWED!')
                        raise NotImplementedError('NWP icon 2D dimensions datasets not implemented yet')
                    elif self.var_dims_data == 'var3d':

                        [var_da, time_da, geox_obj, geoy_obj] = read_data_icon_2i(
                            self.file_data_first, data_var=self.var_name_data)

                        [var_data_raw, var_time_def, var_geox, var_geoy] = adjust_data_icon_2i(
                            var_da, time_da, geox_obj, geoy_obj)

                        var_data_def = compute_rain_icon_2i(var_data_raw)

                    else:
                        log_stream.error(' ----> Get data ... FAILED! FILE SOURCE DIMS NOT ALLOWED!')
                        raise NotImplementedError('NWP icon case dimension datasets not implemented yet')

                elif self.file_source_data == 'ecmwf0100':

                    if self.var_dims_data == 'var2d':
                        log_stream.error(' ----> Get data ... FAILED! FILE SOURCE DIMS NOT ALLOWED!')
                        raise NotImplementedError('NWP ecmwf0100 2D dimensions datasets not implemented yet')
                    elif self.var_dims_data == 'var3d':

                        [var_da, time_da, geox_da, geoy_da] = read_data_ecmwf_0100(
                            self.file_data_first, data_var=self.var_name_data)

                        [var_data_raw,  var_time_def, var_geox, var_geoy] = adjust_data_ecmwf_0100(
                            var_da, time_da, geox_da, geoy_da)

                        var_data_def = compute_rain_ecmwf_0100(var_data_raw)

                    else:
                        log_stream.error(' ----> Get data ... FAILED! FILE SOURCE DIMS NOT ALLOWED!')
                        raise NotImplementedError('NWP ecmwf0100 case dimension datasets not implemented yet')

                else:
                    log_stream.error(' ----> Get data ... FAILED! FILE SOURCE LIBRARY NOT ALLOWED!')
                    raise NotImplementedError('GRIB datasets type not implemented yet')

            elif self.file_format_data == 'tiff' or self.file_format_data == 'tif':

                if self.file_source_data == 'moloc_15':

                    if self.var_dims_data == 'var2d':

                        if self.var_type_data == 'accumulated_classic':
                            file_data_list = self.file_data_list
                        else:
                            log_stream.error(' ----> Get data ... FAILED! FILE SOURCE TYPE NOT ALLOWED!')
                            raise NotImplementedError('NWP WRF type datasets not implemented yet')

                        [var_data_def, var_time_def,
                         var_geox, var_geoy] = read_data_moloc(file_data_list, file_time_data)

                    elif self.var_dims_data == 'var3d':
                        log_stream.error(' ----> Get data ... FAILED! FILE SOURCE DIMS NOT ALLOWED!')
                        raise NotImplementedError('NWP MOLOC 3D dimensions datasets not implemented yet')
                    else:
                        log_stream.error(' ----> Get data ... FAILED! FILE SOURCE DIMS NOT ALLOWED!')
                        raise NotImplementedError('NWP MOLOC case dimension datasets not implemented yet')

            elif self.file_format_data == 'netcdf':

                if self.file_source_data == 'wrf':

                    if self.var_dims_data == 'var2d':

                        if self.var_type_data == 'accumulated_from_first_step':
                            file_data_list = [self.file_data_first] + self.file_data_list
                        elif self.var_type_data == 'accumulated_classic':
                            file_data_list = self.file_data_list
                        else:
                            log_stream.error(' ----> Get data ... FAILED! FILE SOURCE TYPE NOT ALLOWED!')
                            raise NotImplementedError('NWP WRF type datasets not implemented yet')

                        [var_data_raw, var_time_raw,
                         var_geox, var_geoy] = read_data_wrf(file_data_list)

                        var_time_def = convert_time_wrf(var_time_raw, file_time_data, self.var_type_data)
                        var_data_def = convert_data_wrf(var_data_raw, self.var_units_data, self.var_type_data)

                    elif self.var_dims_data == 'var3d':
                        log_stream.error(' ----> Get data ... FAILED! FILE SOURCE DIMS NOT ALLOWED!')
                        raise NotImplementedError('NWP WRF 3D dimensions datasets not implemented yet')
                    else:
                        log_stream.error(' ----> Get data ... FAILED! FILE SOURCE DIMS NOT ALLOWED!')
                        raise NotImplementedError('NWP WRF case dimension datasets not implemented yet')

                elif self.file_source_data == 'gfs025':

                    if self.var_dims_data == 'var2d':
                        log_stream.error(' ----> Get data ... FAILED! FILE SOURCE DIMS NOT IMPLEMENTED!')
                        raise NotImplementedError

                    elif self.var_dims_data == 'var3d':

                        [var_data_def, var_time_def, var_geox, var_geoy] = read_data_gfs_025(
                            self.file_data_first, var_step_type=self.var_type_data)
                else:
                    log_stream.error(' ----> Get data ... FAILED! FILE SOURCE LIBRARY NOT IMPLEMENTED!')
                    raise NotImplementedError('NetCDF datasets type not implemented yet')

            elif self.file_format_data == 'csv':

                if self.file_source_data == 'expert_forecast':
                    if self.var_dims_data == 'var1d':

                        if self.var_type_data == 'accumulated_over_domain':
                            file_data_list = self.file_data_list
                        else:
                            log_stream.error(' ----> Get data ... FAILED! FILE SOURCE TYPE NOT ALLOWED!')
                            raise NotImplementedError('Expert Forecast type datasets not implemented yet')

                        [var_data_raw, var_time_raw] = read_data_expert_forecast(
                            file_data_list,
                            time_tag_start_run_in='time', time_tag_start_forecast_in='time',
                            var_tag_name_list_in=
                                ['rain_average', 'rain_peak', 'slope_x', 'slope_y', 'slope_t'],
                            time_tag_start_run_out='time_start_forecast', time_tag_start_forecast_out='time_start_forecast',
                            var_tag_name_list_out=[
                                'rain_average', 'rain_peak', 'slope_x', 'slope_y', 'slope_t']
                        )

                        var_data_def, var_time_def, var_info_def = configure_data_expert_forecast(
                            var_data_raw, var_time_raw,
                            time_tag_base_start='time_start_forecast', time_tag_extended_start='time_start_forecast',
                            time_freq='H', time_period=12)

                    else:
                        log_stream.error(' ----> Get data ... FAILED! FILE SOURCE DIMS NOT ALLOWED!')
                        raise NotImplementedError('CSV case dimension datasets not implemented yet')
                else:
                    log_stream.error(' ----> Get data ... FAILED! FILE SOURCE LIBRARY NOT IMPLEMENTED!')
                    raise NotImplementedError('Expert Forecast  datasets type not implemented yet')

            elif self.file_format_data == 'json':

                if self.file_source_data == 'expert_forecast':
                    if self.var_dims_data == 'var1d':

                        if self.var_type_data == 'accumulated_over_domain':
                            file_data_list = self.file_data_list
                        else:
                            log_stream.error(' ----> Get data ... FAILED! FILE SOURCE TYPE NOT ALLOWED!')
                            raise NotImplementedError('Expert Forecast type datasets not implemented yet')

                        [var_data_raw, var_time_raw] = read_data_expert_forecast(
                            file_data_list,
                            time_tag_start_run_in='time_start', time_tag_start_forecast_in='time',
                            var_tag_name_list_in=[
                                'time_start', 'time', 'rain_int_start', 'rain_average', 'rain_peak',
                                'slope_x', 'slope_y', 'slope_t', 'pth'],
                            time_tag_start_run_out='time_start_run', time_tag_start_forecast_out='time_start_forecast',
                            var_tag_name_list_out=[
                                'time_start_run', 'time_start_forecast', 'rain_start_run', 'rain_average', 'rain_peak',
                                'slope_x', 'slope_y', 'slope_t', 'pth']
                        )

                        var_data_def, var_time_def, var_info_def = configure_data_expert_forecast(
                            var_data_raw, var_time_raw,
                            time_tag_base_start='time_start_run', time_tag_extended_start='time_start_forecast',
                            time_freq='H', time_period=12)

                    else:
                        log_stream.error(' ----> Get data ... FAILED! FILE SOURCE DIMS NOT ALLOWED!')
                        raise NotImplementedError('CSV case dimension datasets not implemented yet')
                else:
                    log_stream.error(' ----> Get data ... FAILED! FILE SOURCE LIBRARY NOT IMPLEMENTED!')
                    raise NotImplementedError('Expert Forecast  datasets type not implemented yet')

            else:
                log_stream.error(' ----> Get data ... FAILED! FILE TYPE LIBRARY NOT IMPLEMENTED!')
                raise NotImplementedError('NWP datasets format not implemented yet')

            # Create data array
            if var_data_def is not None:

                if self.var_dims_data == 'var2d' or self.var_dims_data == 'var3d':

                    # Create darray for 2d or 3d variable(s)
                    var_obj = create_darray_3d(var_data_def, var_time_def, var_geox, var_geoy)
                    # Dump tmp file
                    write_obj(self.file_tmp, var_obj)

                    '''
                    # DEBUG START (DUMP DATA IN TIFF NETCDF FORMAT)
                    from rfarm.lib_rfarm_utils_generic import writeGeoTiff
                    file_name = '/home/fabio/test/rfarm/lami_2i_model.tiff'
                    writeGeoTiff(file_name, var_data_cmp[:,:,1], var_geox, var_geoy)

                    file_test = self.file_tmp.split()[0] + '.nc'
                    dset_test = var_da.to_dataset(name='Rain')
                    dset_test.to_netcdf(path=file_test)
                    # DEBUG END
                    '''

                elif self.var_dims_data == 'var1d':

                    # Create dictionary for 1d variable(s)
                    var_obj = {'data': var_data_def, 'time': var_time_def, 'info': var_info_def}

                    # Dump tmp file
                    write_obj(self.file_tmp, var_obj)

                # Ending info
                log_stream.info(' ----> Get data ... DONE')

            else:
                # Ending info
                var_obj = None
                log_stream.warning(' ----> Get data ... FAILED! DATA NOT AVAILABLE!')

        else:
            # Ending info
            var_obj = None
            log_stream.warning(' ----> Get data ... FAILED! FIRST FILE OF DATA NOT FOUND')

        return var_obj

# -------------------------------------------------------------------------------------
