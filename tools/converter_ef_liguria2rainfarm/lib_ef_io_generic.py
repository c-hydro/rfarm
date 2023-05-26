"""
Library Features:

Name:          lib_ef_io_generic
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20210104'
Version:       '1.0.0'
"""
#######################################################################################
# Libraries
import logging
import tempfile
import os
import json
import pickle
import re
import xarray as xr
import pandas as pd
import numpy as np

import datetime as dtime
from datetime import datetime
from copy import deepcopy
from scipy.io import loadmat

import matplotlib.pylab as plt
#######################################################################################


# -------------------------------------------------------------------------------------
# Method to create a data array
def create_darray_2d(data, geo_x, geo_y, geo_1d=True, name='geo',
                     coord_name_x='west_east', coord_name_y='south_north',
                     dim_name_x='west_east', dim_name_y='south_north',
                     dims_order=None):

    if dims_order is None:
        dims_order = [dim_name_y, dim_name_x]

    if geo_1d:
        if geo_x.shape.__len__() == 2:
            geo_x = geo_x[0, :]
        if geo_y.shape.__len__() == 2:
            geo_y = geo_y[:, 0]

        data_da = xr.DataArray(data,
                               dims=dims_order,
                               coords={coord_name_x: (dim_name_x, geo_x),
                                       coord_name_y: (dim_name_y, geo_y)},
                               name=name)
        data_da.name = name
    else:
        logging.error(' ===> Longitude and Latitude must be 1d')
        raise IOError('Variable shape is not valid')

    return data_da
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to save data info in json format
def save_file_json(file_name, file_data_dict, file_attrs=None, file_indent=4):

    file_workspace = {}
    for file_key, file_value in file_data_dict.items():

        if isinstance(file_value, dict):
            file_time_list = []
            file_data_list = []
            for file_time_step, file_data_step in file_value.items():

                file_time_list.append(file_time_step.strftime('%Y-%m-%d %H:%M'))
                file_data_list.append(file_data_step)

                if 'time' not in list(file_workspace.keys()):
                    file_workspace['time'] = file_time_list

                file_workspace[file_key] = file_data_list

        else:
            logging.error(' ===> Error in getting datasets')
            raise RuntimeError('Datasets case not implemented yet')

    if file_attrs is not None:
        for attr_key, attr_data in file_attrs.items():
            file_workspace[attr_key] = attr_data

    file_data = json.dumps(file_workspace, indent=file_indent, ensure_ascii=False, sort_keys=False)
    #file_data = re.sub(r'",\s+', '", ', file_data)

    with open(file_name, "w", encoding='utf-8') as file_handle:
        file_handle.write(file_data)

    pass

# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Create default dataframe
def create_default_dframe(df_columns, df_shape, df_nodata=0.0):

    df_data = np.zeros(shape=df_shape)
    df_data[:, :] = df_nodata

    df_obj = pd.DataFrame(data=df_data, columns=df_columns)

    return df_obj
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to read csv file
def read_file_csv(file_name, file_sep=';', file_skiprows=0, tag_time_dim='time'):

    file_dframe = pd.read_table(file_name, sep=file_sep, skiprows=file_skiprows)

    file_dframe.columns = file_dframe.columns.str.strip()
    file_dframe = file_dframe.loc[:, ~file_dframe.columns.str.contains('^Unnamed')]
    file_dframe = file_dframe.replace(to_replace=',', value='', regex=True)

    file_dframe = file_dframe.rename(columns={'data previsione': tag_time_dim})

    return file_dframe

# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to read info file Liguria
def read_file_ForecastLiguria(file_name):

    file = open(file_name,'r')
    list_EF_liguria = file.readlines()
    file.close()

    return list_EF_liguria

# -------------------------------------------------------------------------------------


def OrganizeDataLiguria(list_EF_liguria,slope_y_data, slope_t_data,pth_data):
    a2oDATA = {}
    a2oDATA['AreaAllerta'] = {}
    sDateForecast=list_EF_liguria[1]
    sline=list_EF_liguria[7]
    a1dRain=np.array(sline.split(' '))
    sline = list_EF_liguria[9]
    a1dRain3H = np.array(sline.split(' '))
    sline=list_EF_liguria[11]
    a1dTwin=np.array(sline.split(' '))
    sline=list_EF_liguria[13]
    a1dRainStart=np.array(sline.split(' '))
    iNumAA=a1dRain.__len__()
    if(iNumAA==5):
        a1dList=['alert_area_a','alert_area_b','alert_area_c','alert_area_d','alert_area_e']
        a1dID=['1','3','5','2','4']
    else:
        a1dList = ['alert_area_a', 'alert_area_b', 'alert_area_cm', 'alert_area_d', 'alert_area_e', 'alert_area_m']
        a1dID = ['1', '3', '5', '2', '4','6']

    #
    sDatePP=str(list_EF_liguria[3])
    sDatePP=sDatePP[0:8]+'0000'
    #oDateIni = datetime.datetime(int(sDatePP[0:4]), int(sDatePP[4:6]), int(sDatePP[6:8]), \
    #                             int(sDatePP[8:10]), int(sDatePP[10:12]))
    oDateIni = dtime.datetime(int(sDatePP[0:4]), int(sDatePP[4:6]), int(sDatePP[6:8]), \
                                 int(sDatePP[8:10]), int(sDatePP[10:12]))

    #Initalize AA
    for i in range(0,iNumAA):

        a2oDATA['AreaAllerta'][str(a1dList[i])]={}


    for i in range(0,iNumAA):
        #ID
        a2oDATA['AreaAllerta'][str(a1dList[i])]['id']={}
        a2oDATA['AreaAllerta'][str(a1dList[i])]['id']=str(a1dID[i])
        dDD=float(a1dTwin[i])-12
        if(dDD<0):
            dDD=0
        dDS=dDD*3600
        #Date
        oDateWin = oDateIni + dtime.timedelta(seconds=dDS)
        sDateWin = str(oDateWin)
        sDeteIni = str(oDateIni)
        a2oDATA['AreaAllerta'][str(a1dList[i])]['time_start'] = {}
        a2oDATA['AreaAllerta'][str(a1dList[i])]['time_start'] = [sDeteIni]
        a2oDATA['AreaAllerta'][str(a1dList[i])]['time'] = {}
        a2oDATA['AreaAllerta'][str(a1dList[i])]['time'] = [sDateWin]
        dRainStart=float(a1dRainStart[i])/dDD
        #a2oDATA['AreaAllerta'][str(a1dList[i])]['time'][1] = sDateWin
        a2oDATA['AreaAllerta'][str(a1dList[i])]['rain_int_start'] = {}
        a2oDATA['AreaAllerta'][str(a1dList[i])]['rain_int_start'] = [float(dRainStart)]
        a2oDATA['AreaAllerta'][str(a1dList[i])]['rain_average'] = {}
        a2oDATA['AreaAllerta'][str(a1dList[i])]['rain_average'] = [int(a1dRain[i])]
        a2oDATA['AreaAllerta'][str(a1dList[i])]['rain_peak'] = {}
        a2oDATA['AreaAllerta'][str(a1dList[i])]['rain_peak'] =[int(a1dRain3H[i])]
        a2oDATA['AreaAllerta'][str(a1dList[i])]['slope_y'] = {}
        a2oDATA['AreaAllerta'][str(a1dList[i])]['slope_y'] =[float(slope_y_data[i])]
        a2oDATA['AreaAllerta'][str(a1dList[i])]['slope_x'] = {}
        a2oDATA['AreaAllerta'][str(a1dList[i])]['slope_x'] = [float(slope_y_data[i])]
        a2oDATA['AreaAllerta'][str(a1dList[i])]['slope_t'] = {}
        a2oDATA['AreaAllerta'][str(a1dList[i])]['slope_t'] = [float(slope_t_data[i])]
        a2oDATA['AreaAllerta'][str(a1dList[i])]['pth'] = {}
        a2oDATA['AreaAllerta'][str(a1dList[i])]['pth'] =[float(pth_data[i])]
        a2oDATA['AreaAllerta'][str(a1dList[i])]['name'] = {}
        a2oDATA['AreaAllerta'][str(a1dList[i])]['name'] = str(a1dList[i])
        a2oDATA['AreaAllerta'][str(a1dList[i])]['file_date'] = {}
        a2oDATA['AreaAllerta'][str(a1dList[i])]['file_date'] = sDateForecast[0:-1]
        a2oDATA['AreaAllerta'][str(a1dList[i])]['nodata_value'] = {}
        a2oDATA['AreaAllerta'][str(a1dList[i])]['nodata_value'] = 0

        #a2oDATA['AreaAllerta']['rain_average']
    return a2oDATA

# -------------------------------------------------------------------------------------
# Method to write Liguria ef
def write_ef_liguria(oInputData,a2oDATA):

    # Write on file (PORTARE STA PARTE IN UNA FUNZIONE
    #sPathOut = str(oInputData['data']['dynamic']['destination']['folder_name'])
    #dest_path_time = oInputData['algorithm']['template']['destination_sub_path_time']

    a2oTime = {}
    a2oTime = oInputData['algorithm']['template']
    sDateFor = a2oDATA['AreaAllerta']['alert_area_a']['file_date']
    stime = pd.Timestamp(sDateFor)
    sPathOut = str(oInputData['data']['dynamic']['destination']['folder_name'])
    sPathEF = define_file_pathout(a2oTime, stime, sPathOut)


    ditem = a2oDATA['AreaAllerta'].items()
    for key, value in ditem:
        print(key, '->', value)
        file_indent = 4
        file_data = json.dumps(value, indent=file_indent, ensure_ascii=False, sort_keys=False)

        #Use the date WHEN forecast is done
        sDateFor=a2oDATA['AreaAllerta'][key]['file_date']
        file_name = key + '.json'
        #sPathEF = os.path.join(sPathOut,sDateFor[0:4],sDateFor[4:6],sDateFor[6:8])
        #file_name = os.path.join(sPathOut,sDateFor[0:4],sDateFor[4:6],sDateFor[6:8], file_name)
        file_name = os.path.join(sPathEF, file_name)
        check = os.path.exists(sPathEF)
        if (check == False):
            os.makedirs(sPathEF)

        with open(file_name, "w", encoding='utf-8') as file_handle:
            file_handle.write(file_data)

    return


# -------------------------------------------------------------------------------------
# Method to define ancillary filename
def define_file_pathout(alg_template_tags, time, folder_name_raw):

    #alg_template_tags = self.alg_template_tags

    file_path_dict = {}
    #for domain_name in self.domain_name_list:

    alg_template_values = {'destination_sub_path_time': time}

    folder_name_def = fill_tags2string_2(folder_name_raw, alg_template_tags, alg_template_values)
    #file_name_def = fill_tags2string(file_name_raw, alg_template_tags, alg_template_values)

    file_path_def = os.path.join(folder_name_def)

    #file_path_dict[domain_name] = file_path_def

    return file_path_def


# -------------------------------------------------------------------------------------
# Method to add time in a unfilled string (path or filename)
def fill_tags2string_2(string_raw, tags_format=None, tags_filling=None):

    apply_tags = False
    if string_raw is not None:
        for tag in list(tags_format.keys()):
            if tag in string_raw:
                apply_tags = True
                break

    if apply_tags:

        tags_format_tmp = deepcopy(tags_format)
        for tag_key, tag_value in tags_format.items():
            tag_key_tmp = '{' + tag_key + '}'
            if tag_value is not None:
                if tag_key_tmp in string_raw:
                    string_filled = string_raw.replace(tag_key_tmp, tag_value)
                    string_raw = string_filled
                else:
                    tags_format_tmp.pop(tag_key, None)

        for tag_format_name, tag_format_value in list(tags_format_tmp.items()):

            if tag_format_name in list(tags_filling.keys()):
                tag_filling_value = tags_filling[tag_format_name]
                if tag_filling_value is not None:

                    if isinstance(tag_filling_value, datetime):
                        tag_filling_value = tag_filling_value.strftime(tag_format_value)

                    if isinstance(tag_filling_value, (float, int)):
                        tag_filling_value = tag_format_value.format(tag_filling_value)

                    string_filled = string_filled.replace(tag_format_value, tag_filling_value)

        string_filled = string_filled.replace('//', '/')
        return string_filled
    else:
        return string_raw
# -------------------------------------------
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# Method to read mat file
def read_file_mat(file_name, var_name='vm'):

    file_data = loadmat(file_name)
    if var_name in list(file_data.keys()):
        var_data = file_data[var_name]
    else:
        logging.warning(' ===> Variable not found in mat file. Return none value')
        var_data = None
    return var_data

# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to create a tmp name
def create_filename_tmp(prefix='tmp_', suffix='.tiff', folder=None):

    if folder is None:
        folder = '/tmp'

    with tempfile.NamedTemporaryFile(dir=folder, prefix=prefix, suffix=suffix, delete=False) as tmp:
        temp_file_name = tmp.name
    return temp_file_name
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to get file settings in json format
def read_file_settings(file_name_settings):
    if os.path.exists(file_name_settings):
        with open(file_name_settings) as file_handle:
            data_settings = json.load(file_handle)
    else:
        logging.error(' ===> Error in reading algorithm settings file')
        raise IOError('File not found')
    return data_settings
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to read data obj
def read_obj(filename):
    if os.path.exists(filename):
        data = pickle.load(open(filename, "rb"))
    else:
        data = None
    return data
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to write data obj
def write_obj(filename, data):
    if os.path.exists(filename):
        os.remove(filename)
    with open(filename, 'wb') as handle:
        pickle.dump(data, handle, protocol=pickle.HIGHEST_PROTOCOL)
# -------------------------------------------------------------------------------------
