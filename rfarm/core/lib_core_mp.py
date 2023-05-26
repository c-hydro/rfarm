"""
Library Features:

Name:           lib_core_regrid
Author(s):      Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:           '20170530'
Version:        '3.5.0'
"""

#######################################################################################
# Logging
import logging

import os
import multiprocessing as mp

from rfarm.settings.lib_args import logger_name

log_stream = logging.getLogger(logger_name)

# Debug
# import matplotlib.pylab as plt
#######################################################################################


# -------------------------------------------------------------------------------------
# method to get process id
def get_process_id(process_name):
    # Info
    log_stream.info(' ------> Info: ' + str(process_name) + ' ModuleName: ' + str(__name__))

    if hasattr(os, 'getppid'):  # only available on Unix
        log_stream.info(' -------> Parent process id: ' + str(os.getppid()))

    log_stream.info(' -------> Process id: ' + str(os.getppid()))

# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# method to get process signal start
def get_process_signal_start():
    # Info
    log_stream.info(' ------> Process: ' + str(mp.current_process().name) + ' ... START')
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# method to get process signal end
def get_process_signal_end(process_obj):
    # Info
    log_stream.info(' ------> Process: ' + str(process_obj.name) + ' ExitCode: ' +
                    str(process_obj.exitcode) + ' ... CLOSED')

# -------------------------------------------------------------------------------------
