
"""
Library Features:

Name:          drv_model_exec
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20190903'
Version:       '1.0.0'
"""

# -------------------------------------------------------------------------------------
# Library
import logging
import os

import numpy as np
import pandas as pd

from copy import deepcopy

from rfarm.utils.lib_utils_system import create_tmp     # FX DA RIVEDERE

from rfarm.settings.lib_args import logger_name, time_format
from rfarm.utils.lib_utils_generic import fill_tags2string, convert_obj2dict

from rfarm.core.lib_core_generic import extend_grid, compute_grid, \
    compute_ensemble, compute_var, check_result, save_result_nwp, save_result_expert_forecast
from rfarm.core.lib_core_time import compute_time_steps

import rfarm.core.lib_core_app as lib_core_apps

# Debug
import matplotlib.pylab as plt

# Logging
log_stream = logging.getLogger(logger_name)
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Dictionary of model default parameter(s)
model_parameters_default = {
    'ensemble': {'start': 1, 'end': 2},         # ensemble n (min values EnsMin=1, EnsMax=1)
    'ratio_s': 4,          	                    # spatial disaggregated ratio (min value == 1)
    'ratio_t': 6,                               # time disaggregated ratio (min value == 1)
    'slope_s': None,                            # fft spatial slope (undefined == None)
    'slope_t': None,                            # fft temporal slope (undefined == None)
    'cs_sf': 4,                                 # reliable spatial scale (spatial/Csst)
    'ct_sf': 2,                                 # reliable time scale (time model aggregated Ctsf times)
    'multi_core': False,                        # multi core process (False or True)
    'domain_extension': 0,      	            # domain extended buffer (min value = 0) [km]
    'folder_tmp': None,                         # tmp folder to store data
    'filename_tmp': 'rf_{ensemble}.pkl',        # tmp filename to store data
    'rain_max_thr': 180                         # rain maximum threshold to limit unfair values due to wrong scales
}
# List of model available algorithm(s)
model_algorithm_type = ['exec_nwp', 'exec_expert_forecast']
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Class to configure and run Rainfarm model
class RFarmModel:

    # -------------------------------------------------------------------------------------
    # Initialize class variable(s)
    time_steps_in, time_steps_rel, time_steps_rf, time_steps_ref = None, None, None, None
    time_delta_in, time_delta_ref = None, None
    lons_ref, lats_ref = None, None

    lons_rf, lats_rf, index_rf = None, None, None

    ll_lon_rf, ll_lat_rf, i_min_rf, i_max_rf, j_min_rf, j_max_rf = None, None, None, None, None, None
    lon_min_rf, lon_max_rf, lat_min_rf, lat_max_rf = None, None, None, None

    data_rf, data_slopes = None, None

    nt, ndelta = None, None
    ns, nsl, nr, nas, ntl, nat = None, None, None, None, None, None

    info_common, info_domain = None, None

    ensemble_status = None
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to initialize class
    def __init__(self,
                 ensemble_n=model_parameters_default['ensemble'],
                 ensemble_format='{:03d}',
                 ratio_s=model_parameters_default['ratio_s'],
                 ratio_t=model_parameters_default['ratio_t'],
                 slope_s=model_parameters_default['slope_s'],
                 slope_t=model_parameters_default['slope_t'],
                 cs_sf=model_parameters_default['cs_sf'],
                 ct_sf=model_parameters_default['ct_sf'],
                 multi_core=model_parameters_default['multi_core'],
                 domain_extension=model_parameters_default['domain_extension'],
                 folder_tmp=model_parameters_default['folder_tmp'],
                 filename_tmp=model_parameters_default['filename_tmp'],
                 rain_max_thr=model_parameters_default['rain_max_thr'],
                 model_algorithm="exec_nwp",
                 model_var="Rain",
                 model_metagauss=None,
                 ):

        self.ratio_s = ratio_s
        self.ratio_t = ratio_t
        self.slope_s = slope_s
        self.slope_t = slope_t
        self.cs_sf = cs_sf
        self.ct_sf = ct_sf
        self.multi_core = multi_core
        self.domain_extension = domain_extension
        self.rain_max_thr = rain_max_thr

        if folder_tmp is None:
            self.folder_tmp = create_tmp(model_parameters_default['folder_tmp'])
        else:
            self.folder_tmp = folder_tmp
        if filename_tmp is None:
            self.filename_tmp = model_parameters_default['filename_tmp']
        else:
            self.filename_tmp = filename_tmp

        self.model_algorithm = model_algorithm
        self.model_var = model_var
        self.model_metagauss = model_metagauss

        # check the model algorithm mode
        if self.model_algorithm not in model_algorithm_type:
            log_stream.error(' ===> Algorithm of "' + self.model_algorithm + '" is not supported!')
            raise NotImplemented('The algorith mode is not implemented yet')

        # define ensemble elements
        self.ensemble = compute_ensemble(ensemble_n['start'], ensemble_n['end'])

        self.ensemble_filename = []
        for ensemble_id in self.ensemble:

            tags_tmpl = {'ensemble': ensemble_format}
            tags_values = {'ensemble':  ensemble_format.format(ensemble_id)}
            folder_tmp = fill_tags2string(self.folder_tmp, tags_tmpl, tags_values)
            filename_tmp = fill_tags2string(self.filename_tmp, tags_tmpl, tags_values)
            self.ensemble_filename.append(os.path.join(folder_tmp, filename_tmp))

            if not os.path.exists(folder_tmp):
                os.makedirs(folder_tmp)

        if hasattr(lib_core_apps, self.model_algorithm):
            self.lib_algorithm = getattr(lib_core_apps, self.model_algorithm)
        else:
            log_stream.error(' ===> Algorithm of "' + self.model_algorithm + '" is not available!')
            raise ModuleNotFoundError('Check if the method is currently implemented')

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to configure model grid(s)
    def configure_grid(self, lons_in, lats_in, lons_ref, lats_ref, res_lon_ref, res_lat_ref, domain_extension=0):

        # method to extend rainfarm model grid
        self.lons_ref, self.lats_ref = extend_grid(lons_in, lats_in,
                                                   lons_ref, lats_ref, res_lon_ref, res_lat_ref,
                                                   domain_extension)

        # method to compute rainfarm model grid
        [self.lons_rf, self.lats_rf, self.index_rf,
         self.ll_lon_rf, self.ll_lat_rf,
         self.i_min_rf, self.i_max_rf, self.j_min_rf, self.j_max_rf,
         self.lon_min_rf, self.lon_max_rf, self.lat_min_rf, self.lat_max_rf,
         ratio_s_upd, res_geo_rf,
         res_pixels_rf] = compute_grid(lons_in, lats_in,
                                       self.lons_ref, self.lats_ref, res_lon_ref, res_lat_ref,
                                       self.ratio_s)

        # define rainfarm parameter(s)
        self.ns = res_geo_rf
        self.nsl = res_pixels_rf
        self.nr = 1
        self.nas = res_pixels_rf

        self.ratio_s = ratio_s_upd

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to configure model info
    def configure_info(self, info_domain):

        if info_domain:

            info_time_run_nat_obj, info_time_run_dim_obj, info_time_period = {}, {}, None
            for info_key, info_fields in info_domain.items():

                time_ext_start_step = info_fields['time_extended_start']
                time_ext_end_step = info_fields['time_extended_end']
                time_run_period_step = info_fields['time_run_period']

                # compute time length and nat
                if info_fields['time_split']:
                    time_run_select_step = time_run_period_step[
                        (time_run_period_step >= time_ext_start_step) & (time_run_period_step <= time_ext_end_step)]
                else:
                    time_run_select_step = deepcopy(time_run_period_step)

                if info_time_period is None:
                    info_time_period = deepcopy(time_run_period_step)
                else:
                    info_time_period = info_time_period.union(time_run_period_step)

                info_time_run_dim_step = time_run_select_step.shape[0]
                info_time_run_nat_step = deepcopy(info_time_run_dim_step)

                info_domain[info_key]['dim_time'] = info_time_run_dim_step
                info_domain[info_key]['dim_nat'] = info_time_run_nat_step

            info_time_run_dim_list, info_time_run_nat_list = [], []
            for info_key, info_fields in info_domain.items():

                time_base_start_step = info_fields['time_base_start']
                time_base_end_step = info_fields['time_base_end']
                time_ext_start_step = info_fields['time_extended_start']
                time_ext_end_step = info_fields['time_extended_end']

                info_time_run_dim_list.append(info_domain[info_key]['dim_time'])
                info_time_run_nat_list.append(info_domain[info_key]['dim_nat'])

                idx_base_start_step = info_time_period.get_loc(time_base_start_step)
                idx_base_end_step = info_time_period.get_loc(time_base_end_step)
                idx_ext_start_step = info_time_period.get_loc(time_ext_start_step)
                idx_ext_end_step = info_time_period.get_loc(time_ext_end_step)

                info_domain[info_key]['idx_base_start'] = idx_base_start_step
                info_domain[info_key]['idx_base_end'] = idx_base_end_step
                info_domain[info_key]['idx_extended_start'] = idx_ext_start_step
                info_domain[info_key]['idx_extended_end'] = idx_ext_end_step

            # compute maximum values
            info_time_run_dim_max = max(info_time_run_dim_list)
            info_time_run_nat_max = max(info_time_run_nat_list)

            info_common = {'time_period': info_time_period, 'dim_time': info_time_run_dim_max,
                           'dim_nat': info_time_run_nat_max}
        else:
            info_common = None

        self.info_domain = info_domain
        self.info_common = info_common
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to configure model time(s)
    def configure_time(self, time_in, time_tag_default='generic'):

        # Configure time in
        if isinstance(time_in, pd.DatetimeIndex):   # case for json type 1 dataset (marche)
            time_in = convert_obj2dict(time_in, time_tag_default)
        elif isinstance(time_in, dict):             # case for json type 2 dataset (liguria)
            pass
        elif isinstance(time_in, np.ndarray):       # case for nwp datasets (lami-2i, ecmwf0100, wrf, gfs ... )
            time_in = pd.DatetimeIndex(time_in)
            time_in = convert_obj2dict(time_in, time_tag_default)
        else:
            log_stream.error(' ===> The "time_in" obj must be in datetime index or dict format.')
            raise NotImplemented('Case not implemented yet')

        # iterate on time key(s)
        time_steps_in, time_steps_rel, time_steps_rf, time_steps_ref = {}, {}, {}, {}
        nt, ndelta = {}, {}
        time_delta_collection_in, time_delta_collection_ref = {}, {}
        for time_key, time_value in time_in.items():

            # get time information
            time_start_in = pd.to_datetime(str(time_value[0]))
            time_end_in = pd.to_datetime(str(time_value[-1]))
            time_delta_in = ((time_end_in - time_start_in) / (time_value.__len__() - 1)).seconds

            # Compute time steps for input time length
            time_steps_in[time_key] = compute_time_steps(
                time_start_in.strftime(time_format), time_end_in.strftime(time_format),
                time_delta_in, time_delta_in, time_format)

            # Compute time steps for reliable time length
            time_steps_rel[time_key] = compute_time_steps(
                time_start_in.strftime(time_format), time_end_in.strftime(time_format),
                time_delta_in, time_delta_in * self.ct_sf, time_format)

            # Compute time steps for rainfarm disaggregation
            if self.model_algorithm == 'exec_nwp':
                time_steps_rf[time_key] = compute_time_steps(
                    time_start_in.strftime(time_format), time_end_in.strftime(time_format),
                    time_delta_in, time_delta_in / self.ratio_t, time_format)

            elif self.model_algorithm == 'exec_expert_forecast':
                time_ratio = self.ct_sf / self.ratio_t
                time_steps_rf[time_key] = compute_time_steps(
                    time_start_in.strftime(time_format), time_end_in.strftime(time_format),
                    time_delta_in, time_delta_in * time_ratio, time_format)

            # Compute time steps for output data
            if self.model_algorithm == 'exec_nwp':
                time_steps_ref[time_key] = compute_time_steps(
                    time_start_in.strftime(time_format), time_end_in.strftime(time_format),
                    time_delta_in, time_delta_in * self.ct_sf / (self.ratio_t * self.ct_sf), time_format)
            elif self.model_algorithm == 'exec_expert_forecast':
                time_ratio = self.ct_sf / self.ratio_t
                time_steps_ref[time_key] = compute_time_steps(
                    time_start_in.strftime(time_format), time_end_in.strftime(time_format),
                    time_delta_in, time_delta_in * time_ratio, time_format)

            time_start_ref = pd.to_datetime(time_steps_ref[time_key][0])
            time_end_ref = pd.to_datetime(time_steps_ref[time_key][-1])
            time_delta_ref = ((time_end_ref - time_start_ref) / (time_steps_ref[time_key] .__len__() - 1)).seconds

            if self.ct_sf > self.ratio_t:
                log_stream.warning(' ===> The parameter "ct_sf" [' + str(self.ct_sf) + '] is greater than "ratio_t" ['
                                   + str(self.ratio_t))
                log_stream.warning(
                    ' ===> The effect is a shift in the disaggregated values. To be investigate and correct.')
                # raise RuntimeError('Case not properly managed to under sample the dataset ') # mod francesco

            # Define RF variable(s)
            nt[time_key] = time_steps_ref[time_key].__len__()
            ndelta[time_key] = time_delta_in * self.ct_sf
            # define time delta variable(s)
            time_delta_collection_in[time_key] = time_delta_in
            time_delta_collection_ref[time_key] = time_delta_ref

        # store variable(s) as global dataset(s)
        self.time_steps_in = time_steps_in
        self.time_steps_rel = time_steps_rel
        self.time_steps_rf = time_steps_rf
        self.time_steps_ref = time_steps_ref

        self.nt = nt
        self.ndelta = ndelta

        self.time_delta_in = time_delta_collection_in
        self.time_delta_ref = time_delta_collection_ref

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to configure model data
    def configure_data(self, values_in, info_in):

        # iterate on time key(s)
        data_rf, data_slopes = {}, {}
        ntl, nat = {}, {}
        for (time_key_in, time_delta_in_step), (time_key_ref, time_delta_ref_step) in \
                zip(self.time_delta_in.items(), self.time_delta_ref.items()):

            # check keys
            assert time_key_in == time_key_ref, "the time key for the datasets must be the same"

            # time ratio to convert accumulated rain to instantaneous rain
            time_ratio_rf = int(time_delta_in_step/time_delta_ref_step)

            # compute input datases
            if self.model_algorithm == 'exec_nwp':
                data_rf[time_key_in] = compute_var(
                    values_in, time_ratio_rf, self.i_min_rf, self.i_max_rf, self.j_min_rf, self.j_max_rf)
                data_slopes[time_key_in] = None
            elif self.model_algorithm == 'exec_expert_forecast':

                if self.info_common is not None:
                    info_common_step = self.info_common
                else:
                    log_stream.error(' ===> Obj "info_common" is defined by NoneType')
                    raise RuntimeError('Obj "info_common" is needed by the procedure')
                if time_key_in in list(self.info_domain.keys()):
                    info_domain_step = self.info_domain[time_key_in]
                else:
                    log_stream.error(' ===> Obj "info_domain" is not available for key "' + time_key_in + '"')
                    raise RuntimeError('Obj "info_domain" is needed by the procedure')

                dim_geo_x, dim_geo_y = self.lons_ref.shape[0], self.lats_ref.shape[1]
                dim_time, dim_nat = info_domain_step['dim_time'], info_domain_step['dim_nat']

                data_rf[time_key_in] = np.zeros(shape=[dim_geo_x, dim_geo_y, dim_time])
                data_slopes[time_key_in] = values_in[time_key_in]

            # DEBUG
            # DEBUG START (DUMP DATA IN TIFF NETCDF FORMAT)
            # lon_rf_1d = np.linspace(self.lon_min_rf, self.lon_max_rf, self.nsl, endpoint=True)
            # lat_rf_1d = np.linspace(self.lat_min_rf, self.lat_max_rf, self.nsl, endpoint=True)
            # lon_rf_2d, lat_rf_2d = np.meshgrid(lon_rf_1d, lat_rf_1d)
            # lat_rf_2d = np.flipud(lat_rf_2d)
            # from src.hyde.model.rfarm.lib_rfarm_utils_generic import writeGeoTiff
            # file_name = '/home/fabio/test/rfarm/lami_2i_model_cut.tiff'
            # writeGeoTiff(file_name, self.data_rf[:, :, 1], lon_rf_2d, lat_rf_2d)

            # import matplotlib.pylab as plt
            # plt.figure(1)
            # plt.imshow(values_in[:, :, 0])
            # plt.colorbar()
            # plt.clim(0, 10)
            # plt.figure(2)
            # plt.imshow(self.data_rf[:, :, 0])
            # plt.colorbar()
            # plt.clim(0, 10)
            # plt.show()
            # DEBUG

            # define rainfarm variable(s)
            ntl[time_key_in] = data_rf[time_key_in].shape[2]
            nat[time_key_in] = data_rf[time_key_in].shape[2]

        # store variable(s) as global dataset(s)
        self.data_rf = data_rf
        self.data_slopes = data_slopes

        self.ntl = ntl
        self.nat = nat

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to recover run status
    def callback_run(self):
        ensemble_status = []
        for ensemble_id, (ensemble_filename, ensemble_n) in enumerate(zip(self.ensemble_filename, self.ensemble)):
            # Starting info
            log_stream.error(' ----> Callback ensemble ' + str(ensemble_n) + ' ... ')
            if os.path.exists(ensemble_filename):
                ensemble_status.append(ensemble_filename)
                # Ending info
                log_stream.info(' ----> Callback ensemble ' + str(ensemble_n) + ' ... DONE')
            else:
                # Ending info
                log_stream.info(' ----> Callback ensemble ' + str(ensemble_n) + ' ... FAILED! FILE NOT FOUND')
        return ensemble_status
    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to execute model
    def execute_run(self, domain_name='domain', domain_id=1, domain_mask=None):

        # Choose computing mode
        if self.multi_core:
            log_stream.error(' ===> Multicore application of RainFarm model is not available')
            raise NotImplemented('Case not implemented yet')
            # Define process(es) for RF model using multi core mode
            # ensemble_status = self.worker_multi_core()
        else:
            # Define process(es) for RF model using single core mode
            ensemble_status = self.worker_single_core(domain_name, domain_id, domain_mask)

        # Save data in global workspace
        return ensemble_status

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # method to run model in single core
    def worker_single_core(self, domain_name, domain_id, domain_mask):

        # get variable name
        ensemble_var = self.model_var

        # iterate over ensemble(s)
        ensemble_status = []
        info_common, info_domain = None, None
        for ensemble_id, (ensemble_filename, ensemble_n) in enumerate(zip(self.ensemble_filename, self.ensemble)):

            # info ensemble start
            log_stream.info(' ----> Compute ensemble "' + str(ensemble_n) + '" ... ')

            # select rainfarm mode
            if self.model_algorithm == 'exec_nwp':

                # info mode start
                log_stream.info(' -----> Algorithm mode "' + self.model_algorithm + '" ... ')

                # parse obj data to rainfarm model of nwp mode
                if isinstance(self.data_rf, dict):
                    data_rf = deepcopy(self.data_rf['generic'])
                elif isinstance(self.data_rf, np.ndarray):
                    data_rf = deepcopy(self.data_rf)
                else:
                    log_stream.error(' ===> Format of data obj is not supported')
                    raise NotImplemented('Case not implemented yet')

                # check source data (to evaluate data >=0)
                if np.any(data_rf):

                    # execute rainfarm model
                    [ensemble_rf_common, metagauss_rf] = self.lib_algorithm(
                        data_rf,
                        self.ratio_s, self.ratio_t,
                        cssf=self.cs_sf, ctsf=self.ct_sf,
                        f=self.model_metagauss,
                        celle3_rainfarm=None,
                        sx=self.slope_s, st=self.slope_t)

                    # check rainfarm results
                    check_result(data_rf, ensemble_rf_common, self.ratio_s, self.ratio_t)

                    '''
                    # debug
                    plt.figure(1); plt.imshow(self.data_rf[:, :, 1]); plt.colorbar(); plt.clim(0, 10)
                    plt.figure(2); plt.imshow(ensemble_rf_common[:, :, 1]); plt.colorbar(); plt.clim(0, 10)
                    plt.show()
                    '''

                    # info mode start
                    log_stream.info(' -----> Algorithm mode "' + self.model_algorithm + '" ... DONE')

                else:

                    # exit (all data are equal to zeros)
                    ensemble_rf_common = np.zeros(
                        [self.data_rf.shape[0]*self.ratio_s,
                         self.data_rf.shape[1]*self.ratio_s,
                         self.data_rf.shape[2]*int(int(self.nt)/int(self.nat))])

                    # info mode start
                    log_stream.info(' -----> Algorithm mode "' + self.model_algorithm +
                                    '" ... SKIPPED. All values are defined by zero')

            elif self.model_algorithm == 'exec_expert_forecast':

                # info mode start
                log_stream.info(' -----> Algorithm mode "' + self.model_algorithm + '" ... ')

                # get common and domain information
                info_common = self.info_common
                info_domain = self.info_domain
                # get rain maximum threshold
                rain_max_thr = self.rain_max_thr

                # get alert area index
                subdomain_mask = deepcopy(domain_mask)
                subdomain_shape = subdomain_mask.shape

                # get maximum simulation length
                subdomain_nat_common = info_common['dim_nat'] * self.ratio_t
                subdomain_time_period_common = info_common['time_period']

                # initialize common rf workspace
                ensemble_rf_common = np.zeros([subdomain_shape[0], subdomain_shape[1], subdomain_nat_common])
                # iterate over subdomain(s)
                ensemble_rf_collection, avg_rf_collection = {}, {}
                for area_tag, area_id, (slope_tag, slope_fields) in zip(domain_name, domain_id,
                                                                        self.data_slopes.items()):
                    # info domain start
                    log_stream.info(' ------> Domain "' + area_tag + '" ... ')
                    assert (area_tag == slope_tag)

                    # get simulation length
                    nat = self.nat[area_tag]
                    time_idx_max = int(nat / self.ratio_t)

                    # get domain idx
                    subdomain_idx = np.argwhere(subdomain_mask.ravel() == area_id)
                    if subdomain_idx.ndim == 2:
                        subdomain_idx = subdomain_idx[:, 0]
                    if subdomain_idx.size == 0:
                        logging.error(' ===> Domain "' + area_tag + '" with ID "' + str(area_id) +
                                      '" is not defined in the domain mask file. The domain indexes are empty.')
                        raise RuntimeError('Check the domains file for the domain and id availability')

                    subdomain_cell = np.argwhere(subdomain_mask == area_id)
                    subdomain_cell_x, subdomain_cell_y = subdomain_cell[:, 0], subdomain_cell[:, 1]

                    # iterate over time steps
                    time_idx_list, rain_avg_obj = [], {}
                    avg_rf_collection[area_tag] = {}
                    for time_idx_step in range(0, time_idx_max):

                        # get time start idx and time end idx
                        time_idx_start = time_idx_step * self.ratio_t
                        time_idx_end = time_idx_step * self.ratio_t + self.ratio_t
                        time_idx_list.append([time_idx_start, time_idx_end])

                        # initialize rf dataset
                        ensemble_rf_t = np.zeros([subdomain_shape[0] * subdomain_shape[1], nat])

                        # info idx start
                        log_stream.info(' -------> IdxMin: "' + str(time_idx_start) +
                                        '" -- IdxMax: "' + str(time_idx_end) + '" ... ')

                        # rainfarm parameters
                        rain_avg = slope_fields['rain_average'][time_idx_step]
                        slope_s = slope_fields['slope_x'][time_idx_step]
                        slope_t = slope_fields['slope_t'][time_idx_step]

                        nt = int(nat / self.ratio_t)
                        ratio_t = int(self.ratio_t/nt)
                        ct_sf = self.ratio_t / self.ct_sf

                        if 'pth' in slope_fields.keys():
                            pth = slope_fields['pth'][time_idx_step]
                        else:
                            log_stream.warning(' ===> Parameters "pth" defined by the default value (pth=0.25)')
                            pth = 0.25

                        # execute rainfarm model
                        """
                        Parameters example:
                        rain average: 0.0; slope s: 3.5; slope t: 0.5; ratio s: 1; ratio t/ nat: 3; 
                        cs_sf: 1; ct_sf: 1; nat: 4; ns: 644; ns: 644
                        """
                        [ensemble_rf_xyt, metagauss_rf] = self.lib_algorithm(
                            None,
                            self.ratio_s, ratio_t,
                            cssf=self.cs_sf, ctsf=ct_sf, pth=pth,
                            f=self.model_metagauss,
                            celle3_rainfarm=None,
                            sx=slope_s, st=slope_t,
                            nx=self.ns, ny=self.ns, nt=nt)

                        # post-process results
                        ensemble_rf_tmp = ensemble_rf_xyt[0:subdomain_shape[0], 0:subdomain_shape[1], :]
                        ensemble_rf_step = np.reshape(ensemble_rf_tmp, [subdomain_shape[0] * subdomain_shape[1], self.ratio_t])

                        # compute avg and weights
                        ensemble_rf_avg = np.nanmean(ensemble_rf_step[subdomain_idx, :])
                        if ensemble_rf_avg > 0.0:
                            ensemble_rf_weights = rain_avg / self.ratio_t / ensemble_rf_avg
                        else:
                            ensemble_rf_weights = 0.0

                        ensemble_rf_check = ensemble_rf_step[subdomain_idx, :] * ensemble_rf_weights

                        # debug
                        # rain_max_thr = 10.0

                        # check fields limits
                        log_stream.info(' --------> Check rainfarm fields limits ... ')
                        idx_check = np.where(ensemble_rf_check > rain_max_thr)
                        array_check = np.array(idx_check)
                        if np.size(array_check, 1):
                            ensemble_rf_check[ensemble_rf_check > rain_max_thr] = rain_max_thr
                            ensemble_rf_avg_n = np.nanmean(ensemble_rf_check)
                            ensemble_rf_weights_C = rain_avg / self.ratio_t / ensemble_rf_avg_n
                            ensemble_rf_check_N = ensemble_rf_check/ensemble_rf_weights_C
                            ensemble_rf_avg_nn = np.nanmean(ensemble_rf_check_N)
                            ensemble_rf_weights_C = rain_avg / self.ratio_t / ensemble_rf_avg_nn
                            ensemble_rf_check = ensemble_rf_check_N * ensemble_rf_weights_C
                            rain_rf_max = np.nanmax(np.nanmax(ensemble_rf_check))
                            rain_rf_mean = np.nanmean(np.nanmean(ensemble_rf_check))
                            ensemble_rf_check[ensemble_rf_check > 1.1 * rain_max_thr] = 1.1 * rain_max_thr

                            log_stream.info(' ::: Field Analysis ::: ')
                            log_stream.info(' ::: RainThr: "' + '{:.1f}'.format(rain_max_thr) +
                                            '" -- RainMax: "' + '{:.1f}'.format(rain_rf_max) +
                                            '" -- RainAvg: "' + '{:.1f}'.format(rain_rf_mean) + '" ::: ')
                            log_stream.info(' --------> Check rainfarm fields limits ...  DONE')

                        else:

                            rain_rf_max = np.nanmax(np.nanmax(ensemble_rf_check))
                            rain_rf_mean = np.nanmean(np.nanmean(ensemble_rf_check))
                            log_stream.info(' ::: Field Analysis ::: ')
                            log_stream.info(' ::: RainThr: "' + '{:.1f}'.format(rain_max_thr) +
                                            '" -- RainMax: "' + '{:.1f}'.format(rain_rf_max) +
                                            '" -- RainAvg: "' + '{:.1f}'.format(rain_rf_mean) + '" ::: ')
                            log_stream.info(' --------> Check rainfarm fields limits ... SKIPPED. NOTHING TO DO')

                        # normalize result(s) using weight(s)
                        ensemble_rf_t[subdomain_idx, time_idx_start:time_idx_end] = ensemble_rf_check

                        # reshape results in XYT format
                        ensemble_rf_domain = np.reshape(ensemble_rf_t, [subdomain_shape[0], subdomain_shape[1], nat])

                        # check the average value
                        tmp_rf_avg_act_cmp = np.nanmean(
                            ensemble_rf_domain[subdomain_cell_x, subdomain_cell_y, time_idx_start:time_idx_end])

                        # organize the collection datasets
                        if area_tag not in list(ensemble_rf_collection.keys()):
                            ensemble_rf_collection[area_tag] = deepcopy(ensemble_rf_domain)
                            time_idx_bnd = time_idx_list[time_idx_step]
                            tmp_rf_avg_act_save = np.nanmean(
                                ensemble_rf_domain[subdomain_cell_x, subdomain_cell_y, time_idx_bnd[0]:time_idx_bnd[1]])

                            assert tmp_rf_avg_act_cmp == tmp_rf_avg_act_save, "average values must be the same"

                        else:
                            # get data from the collections
                            ensemble_rf_tmp = deepcopy(ensemble_rf_collection[area_tag])

                            time_idx_bnd = time_idx_list[time_idx_step - 1]
                            tmp_rf_avg_check_prev = np.nanmean(
                                ensemble_rf_tmp[subdomain_cell_x, subdomain_cell_y, time_idx_bnd[0]:time_idx_bnd[1]])

                            # update values in the common workspace
                            ensemble_rf_tmp[
                                subdomain_cell_x, subdomain_cell_y, time_idx_start:time_idx_end] = deepcopy(
                                ensemble_rf_domain[subdomain_cell_x, subdomain_cell_y, time_idx_start:time_idx_end])

                            time_idx_bnd = time_idx_list[time_idx_step]
                            tmp_rf_avg_act_save = np.nanmean(
                                ensemble_rf_tmp[subdomain_cell_x, subdomain_cell_y, time_idx_bnd[0]:time_idx_bnd[1]])

                            assert tmp_rf_avg_act_cmp == tmp_rf_avg_act_save, "average values must be the same"

                            # put values in the collections
                            ensemble_rf_collection[area_tag] = ensemble_rf_tmp

                        # collect the average values
                        avg_rf_collection[area_tag][time_idx_step] = [tmp_rf_avg_act_cmp, tmp_rf_avg_act_save]

                        # info idx end
                        log_stream.info(' -------> IdxMin: "' + str(time_idx_start) +
                                        '" -- IdxMax: "' + str(time_idx_end) + '" ... DONE')
                        '''
                        # debug
                        plt.figure(1)
                        plt.imshow(ensemble_rf_domain[:, :, time_idx_start]); plt.colorbar()
                        plt.figure(2)
                        plt.imshow(ensemble_rf_domain[:, :, time_idx_start]); plt.colorbar()
                        plt.show()
                        '''

                    # info domain end
                    log_stream.info(' ------> Domain "' + area_tag + '" ... DONE')

                # initialize common rf workspace
                ensemble_rf_common = np.zeros([subdomain_shape[0], subdomain_shape[1], subdomain_time_period_common.__len__()])
                # iterate over subdomain(s)
                for area_tag, area_id, (info_tag, info_fields) in zip(domain_name, domain_id, info_domain.items()):

                    time_split = info_fields['time_split']
                    time_base_start, time_base_end = info_fields['time_base_start'], info_fields['time_base_end']
                    time_ext_start, time_ext_end = info_fields['time_extended_start'], info_fields['time_extended_end']

                    rain_base = None
                    if 'rain_base_period' in list(info_fields.keys()):
                        rain_base = info_fields['rain_base_period']

                    idx_base_start = pd.DatetimeIndex(subdomain_time_period_common).get_loc(time_base_start)
                    idx_base_end = pd.DatetimeIndex(subdomain_time_period_common).get_loc(time_base_end)
                    idx_ext_start = pd.DatetimeIndex(subdomain_time_period_common).get_loc(time_ext_start)
                    idx_ext_end = pd.DatetimeIndex(subdomain_time_period_common).get_loc(time_ext_end)

                    ensemble_rf_step = ensemble_rf_collection[area_tag]

                    # get domain idx
                    cell_step = np.argwhere(subdomain_mask == area_id)
                    x_step, y_step = cell_step[:, 0], cell_step[:, 1]

                    # check time splitting (base and forecast or only forecast)
                    if time_split:

                        # fill with rain fields (pass by user) for base period
                        if rain_base is not None:
                            ensemble_rf_common[x_step, y_step, idx_base_start:idx_base_end + 1] = rain_base
                        # fill with rainfarm fields for ext period
                        ensemble_rf_common[
                            x_step, y_step, idx_ext_start:idx_ext_end + 1] = deepcopy(ensemble_rf_step[x_step, y_step, :])

                    else:
                        # fill with rainfarm fields for ext period
                        ensemble_rf_common[
                            x_step, y_step, idx_base_start:idx_ext_end + 1] = deepcopy(ensemble_rf_step[x_step, y_step, :])

                    '''
                    # Debug
                    plt.figure()
                    plt.imshow(ensemble_rf_step[:, :, idx_base_start], interpolation=None); plt.colorbar()
                    plt.show()
                    plt.figure()
                    plt.imshow(ensemble_rf_common[:, :, idx_base_start], interpolation=None); plt.colorbar()
                    plt.show()
                    '''

                # info mode start
                log_stream.info(' -----> Algorithm mode "' + self.model_algorithm + '" ... DONE')

            else:

                log_stream.error(' ===> RainFarm application type is not correctly defined [' +
                                 self.model_algorithm + ']. Check your settings')
                raise NotImplemented('RainFarm application type not implemented yet')

            # store ensemble status
            ensemble_status.append(ensemble_filename)

            # print ensemble info
            ensemble_rf_avg = np.nanmean(ensemble_rf_common)
            ensemble_rf_min = np.nanmin(ensemble_rf_common)
            ensemble_rf_max = np.nanmax(ensemble_rf_common)
            log_stream.info(' -----> Evaluate'
                            ' -- Values Avg: ' + "{:.4f}".format(ensemble_rf_avg) +
                            ' -- Values Min: ' + "{:.4f}".format(ensemble_rf_min) +
                            ' -- Values Max: ' + "{:.4f}".format(ensemble_rf_max)
                            )

            # dump ensemble to free the system memory
            if self.model_algorithm == 'exec_expert_forecast':

                # save results for expert forecast mode
                save_result_expert_forecast(
                    ensemble_filename, var_name=ensemble_var, var_data_in=ensemble_rf_common,
                    var_time=info_common['time_period'],
                    var_geo_x=self.lons_ref, var_geo_y=self.lats_ref)

            elif self.model_algorithm == 'exec_nwp':

                # save results for nwp mode
                save_result_nwp(
                    ensemble_filename, ensemble_var, ensemble_rf_common,
                    self.lons_rf, self.lats_rf, self.time_steps_rf,
                    self.lons_ref, self.lats_ref, self.time_steps_ref,
                    geoindex_in=self.index_rf)

            # info ensemble end
            log_stream.info(' ----> Compute ensemble ' + str(ensemble_n) + ' ... DONE')

        # return ensemble status
        return ensemble_status
        # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------

    # -------------------------------------------------------------------------------------
    # Method to run model in multi cores mode
    def worker_multicore(self):

        # Iterate over ensemble(s)
        ensemble_status = []
        for ensemble_id, ensemble_filename in enumerate(self.ensemble_filename):
            ensemble_status.append(ensemble_filename)
        return ensemble_status
    # -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
