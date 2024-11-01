=========
Changelog
=========

Version 4.4.3 [2024-10-31]
**************************
APP: **RFARM_Model_RUN_Manager.py**
	- Add functions to manage and use different sources in "nwp" mode (icon-2i)
	- Fix bugs related to not initilized conditions


Version 4.4.2 [2023-11-27]
**************************
APP: **RFARM_Model_RUN_Manager.py**
	- Add functions to manage and use different sources in "nwp" mode (moloc)
	- Fix methods to organize disaggregation domain with a buffer zone (to extend the starting grid)
	- Fix methods to manage nans or zeros in the Fourier transform


Version 4.4.1 [2023-06-23]
**************************
APP: **RFARM_Model_RUN_Manager.py**
	- Add functions to manage different sources in "nwp" mode (lami-2i, ecmwf0100, gfs0.25)


Version 4.4.0 [2023-05-26]
**************************
APP: **RFARM_Model_RUN_Manager.py**
	- Add and refactor functions for "expert_forecast" and "nwp" mode
	- Add the "pth" parameter to the "expert_forecast" mode
	- Add the boundary limits to the rainfarm fieldss to the "expert_forecast" mode
	- Support for different time window in the "expert_forecast" mode
APP: **RFARM_Converter_EF_Rain.py**
	- Add the procedures, for Marche and Liguria region, to convert the "expert_forecast" original file
	to the file needed by RainFarm to execute the "expert_forecast" mode
FIX: **RFARM_Model_RUN_Manager.py**
	- Fix bugs in the grid definition
	- Fix bugs in the "expert_forecast" output fields
	- Fix bugs in the fx and methods 
	- Fix bugs related to the libraries versions

Version 4.3.0 [2022-12-19]
**************************
APP: **RFARM_Model_RUN_Manager.py**
	- RainFarm package (refactor from HyDe package previous versions)
FIX: **RFARM_Model_RUN_Manager.py**
	- Fix metagauss domain issue related to the expert-forecast mode

Version 4.2.1 [2021-08-01]
**************************
FIX: **HYDE_Model_RFarm.py**
	- Add resampling of outcome datasets according with the time delta settings
    - Fix 3d array vertical direction for the ecmwf0100 nwp
    - Fix time shift issue for the lami-2i nwp
    - Fix geographical position issue for the irregulared grid of lami-2i nwp

Version 4.2.0 [2021-05-03]
**************************
APP: **HYDE_Model_RFarm.py**
    - Add the support to GFS 0.25 products
    - Update WRF reader for supporting 3d files

Version 4.1.0 [2021-02-02]
**************************
APP: **HYDE_Model_RFarm.py**
	- Add expert forecast routine

FIX: **HYDE_Model_RFarm.py**
	- Adapt scripts and fix bugs

Version 4.0.3 [2021-01-25]
**************************
FIX: **HYDE_Model_RFarm.py**
	- Adapt scripts and fix bugs

Version 4.0.2 [2021-01-13]
**************************
FIX: **HYDE_Model_RFarm.py**
    - Fix and update the data reader of the lami-2i nwp and the ecmwf0100 nwp
    - Fix and update the data reader of the expert forecast datasets
    - Fix and update the model application for the expert forecast case
    
Version 4.0.1 [2020-05-22]
**************************
FIX: **HYDE_Model_RFarm.py**
	- Fix the incorrect slope in time estimation
    - Minor bugs in bash scripts and python scripts

Version 4.0.0 [2019-09-02]
**************************
APP: **HYDE_Model_RFarm.py**
    - Beta release for HyDE package and Python 3.x

DRV: **lib_rfarm_core.py**
	- Beta release for HyDE package for RainFarm model

Version 3.5.2 [2018-09-10]
**************************
APP: **FP_Model_RainFarm.py**
	- Beta release for FloodProofs library

Version 3.5.1 [2017-11-14]
**************************
FIX: **FP_Model_RainFarm.py**
	- Fix bugs (accumulated and istantaneous rain)

Version 3.5.0 [2017-05-30]
**************************
APP: **FP_Model_RainFarm.py**
	- Refactoring of the codes to update the libraries and the applications (pandas and xarray libraries)

Version 3.0.1 [2015-09-24]
**************************
APP: **FP_Model_RainFarm.py**
	- Final release for operational chain mode (i.e. Regional Operational Chain)

Version 3.0.0 [2015-08-23]
**************************
APP: **FP_Model_RainFarm.py**
	- Final release for experimental project (i.e. DRIHM)
	
Version 2.0.1 [2014-04-08]
**************************
APP: **RainFarm.py**
	- Final version for experimental mode (RainFarm Python Language based on RainFarm MatLab version

Version 2.0.0 [2014-01-22]
**************************
APP: **RainFarm.py**
	- Beta version in Python language based on RainFarm MatLab version

Version 1.0.0 [2013-05-02]
**************************
APP: **RainFarm.m**
	- RainFarm MatLab version

