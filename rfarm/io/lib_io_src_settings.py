"""
Library Features:

Name:          lib_io_src_settings
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20201102'
Version:       '3.0.0'
"""
#################################################################################
# Library
import logging
import os
import json

from rfarm.settings.lib_args import logger_name

log_stream = logging.getLogger(logger_name)
#################################################################################


# -------------------------------------------------------------------------------------
# Method to get file settings in json format
def read_file_settings(file_name_settings):
    if os.path.exists(file_name_settings):
        with open(file_name_settings) as file_handle:
            data_settings = json.load(file_handle)
    else:
        log_stream.error(' ===> Error in reading algorithm settings file "' + file_name_settings + '"')
        raise IOError('File not found')
    return data_settings
# -------------------------------------------------------------------------------------
