"""
RFARM Processing Tool - Expert Forecast Rain Processing Tool

__date__ = '20210104'
__version__ = '2.0.0'
__author__ = 'Fabio Delogu (fabio.delogu@cimafoundation.org'
__library__ = 'rfarm'

General command line:
python RFARM_Converter_EF_Rain.py -settings_file configuration.json -time "YYYY-MM-DD HH:MM"

Version(s):
20210104 (2.0.0) --> Beta release for hyde package
20150923 (1.1.0) --> Previous final version
"""
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# Complete library
import logging
import json
import os
import sys
import pandas as pd
import datetime

from argparse import ArgumentParser
from time import time, strftime, gmtime

#from drv_configuration_algorithm_ef import DriverAlgorithm
#from drv_configuration_time_ef import DriverTime

#from drv_data_ef_geo import DriverGeo
#from drv_data_ef_io import DriverData

#Added for Liguria
import lib_ef_io_generic
import lib_ef_analysis

# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# Algorithm information
alg_project = 'RFarm'
alg_name = 'EXPERT FORECAST Rain Processing Tool'
alg_version = '2.0.0'
alg_release = '2021-01-04'
alg_type = 'Converter'
# Algorithm parameter(s)
time_format = '%Y-%m-%d %H:%M'
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Script Main
def main():


    # -------------------------------------------------------------------------------------
    # Get script argument(s)
    [file_script, file_settings, time_arg] = get_args()

    #************ case  Liguria *********************************************************

    oInputData = json.load(open(file_settings))

    """ Costants Variables """
    startime = datetime.datetime.now().strftime("%Y-%m-%d_%H_%M")

    """Definisco il logging"""
    # Definisco il "contenitore" dei logger
    logger = logging.getLogger(oInputData["log"]["docker_name"])
    logger.setLevel(logging.DEBUG)

    # Definisco l'handler che stampera' su terminale
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.DEBUG)

    # Definisco l'handler che scivera' su file
    file_handler = logging.FileHandler(oInputData["log"]["folder_name"] \
                                        + "Log_"+ startime + oInputData["log"]["file_name"])

    # Imposta un formato per i file di log
    formatter = logging.Formatter('%(asctime)s %(name)s %(levelname)s: %(message)s')

    # Imposto il forato dei file di log
    console_handler.setFormatter(formatter)
    file_handler.setFormatter(formatter)

    # Aggiungo gli handlers al logger
    logger.addHandler(console_handler)
    logger.addHandler(file_handler)

    logger.info("Logger Created")
    logger.info("Init")


    sDomain=str(oInputData['algorithm']['info']['domain'])

    # Read expert forecast
    sFileInput=str(oInputData['data']['dynamic']['source']['folder_name']) + \
                    str(oInputData['data']['dynamic']['source']['file_name'])

    list_EF_liguria=lib_ef_io_generic.read_file_ForecastLiguria(sFileInput)

    # Read slopes
    sFileSlopes=str(oInputData['data']['static']['source']['slopes']['folder_name']) + \
                    str(oInputData['data']['static']['source']['slopes']['file_name'])

    a2dTabSlope=lib_ef_io_generic.read_file_mat(sFileSlopes, var_name='Tabella_slopes')

    #Calculate slopes and pth
    [slope_y_data, slope_t_data, pth_data]=lib_ef_analysis.find_slopes_liguria(list_EF_liguria, a2dTabSlope)

    #Organize dictionary
    a2oDATA=lib_ef_io_generic.OrganizeDataLiguria(list_EF_liguria,slope_y_data, slope_t_data,pth_data)

    #Write on file
    lib_ef_io_generic.write_ef_liguria(oInputData,a2oDATA)


    #*******************************************************************************************

    logging.info('[' + alg_project + ' ' + alg_type + ' - ' + alg_name + ' (Version ' + alg_version + ')]')
    # -------------------------------------------------------------------------------------

    iExitStatus = -1;
    sys.exit(iExitStatus)  # exit the program
    #

# -------------------------------------------------------------------------------------
# Method to get script argument(s)
def get_args():

    parser_handle = ArgumentParser()
    parser_handle.add_argument('-settings_file', action="store", dest="alg_settings")
    parser_handle.add_argument('-time', action="store", dest="alg_time")
    parser_values = parser_handle.parse_args()

    alg_script = parser_handle.prog

    if parser_values.alg_settings:
        alg_settings = parser_values.alg_settings
    else:
        alg_settings = 'configuration.json'

    if parser_values.alg_time:
        alg_time = parser_values.alg_time
    else:
        alg_time = None

    return alg_script, alg_settings, alg_time

# -------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------
# Call script from external library
if __name__ == "__main__":
    main()
# ----------------------------------------------------------------------------
