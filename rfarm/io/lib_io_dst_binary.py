"""
Library Features:

Name:          lib_io_binary
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20201102'
Version:       '3.0.0'
"""
#################################################################################
# Library
import logging
import numpy as np
from rfarm.settings.lib_args import logger_name

log_stream = logging.getLogger(logger_name)
#################################################################################


# --------------------------------------------------------------------------------
# method to write dset in binary format
def write_dset_binary(file_name, file_data, file_format='float32'):
    file_data_formatted = np.array(file_data, dtype=file_format)
    file_handle = open(file_name, "wb")  # open file in write mode
    file_data_formatted.tofile(file_handle)
    file_handle.close()
# --------------------------------------------------------------------------------
