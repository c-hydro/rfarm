
"""
Library Features:

Name:          lib_io_generic
Author(s):     Fabio Delogu (fabio.delogu@cimafoundation.org)
Date:          '20210104'
Version:       '1.0.0'
"""

# -------------------------------------------------------------------------------------
# Libraries
import logging
import os
import pickle

import xarray as xr

from rfarm.settings.lib_args import logger_name

# Logging
log_stream = logging.getLogger(logger_name)
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# Method to create a data array
def create_darray_2d(data, geox, geoy):
    darray = xr.DataArray(data,
                          dims=['south_north', 'west_east'],
                          coords={'longitude': (['south_north', 'west_east'], geox),
                                  'latitude': (['south_north', 'west_east'], geoy)})
    return darray
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# Method to create a data array
def create_darray_3d(data, time, geox, geoy):
    darray = xr.DataArray(data,
                          dims=['south_north', 'west_east', 'time'],
                          coords={'time': (['time'], time),
                                  'longitude': (['south_north', 'west_east'], geox),
                                  'latitude': (['south_north', 'west_east'], geoy)})
    return darray
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
