"""
Library Features:

Name:           lib_core_time
Author(s):      Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:           '20230519'
Version:        '3.5.1'
"""

#######################################################################################
# Logging
import logging
import datetime
from copy import deepcopy

from rfarm.settings.lib_args import logger_name

log_stream = logging.getLogger(logger_name)

# Debug
# import matplotlib.pylab as plt
#######################################################################################


# -------------------------------------------------------------------------------------
# method to compute time steps
def compute_time_steps(sTimeFrom, sTimeTo, iTimeDelta_IN, iTimeDelta_OUT, sTimeFormat='%Y%m%d%H%M'):

    # Get time from and time to information
    oTimeFrom = datetime.datetime.strptime(sTimeFrom, sTimeFormat)
    oTimeTo = datetime.datetime.strptime(sTimeTo, sTimeFormat)

    oTimeDelta_IN = datetime.timedelta(seconds=iTimeDelta_IN)
    oTimeDelta_OUT = datetime.timedelta(seconds=iTimeDelta_OUT)

    # Compute initial step
    oTimeStep = oTimeFrom - oTimeDelta_IN
    oTimeStep = oTimeStep + oTimeDelta_OUT

    # Compute time steps OUT
    a1oTimeSteps = []
    while oTimeStep <= oTimeTo:
        a1oTimeSteps.append(oTimeStep.strftime(sTimeFormat))
        oTimeStep += oTimeDelta_OUT

    # Info
    sTimeFrom = a1oTimeSteps[0]
    sTimeTo = a1oTimeSteps[-1]
    log_stream.info(' -------> Time -- From: ' + sTimeFrom + ' To: ' + str(sTimeTo))

    # Return variable(s)
    return a1oTimeSteps
    # -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------

