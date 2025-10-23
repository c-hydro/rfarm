# Changelog

All notable changes to this project are documented here.

---

## Version 4.5.0 — *2025-10-23*
### APP: **RFARM_Model_RUN_Manager.py**
- Add functions to manage and use different sources in "nwp" mode (ecmwf-0100 1h/3h hours resolution)

---

## Version 4.4.3 — *2024-10-31*
### APP: **RFARM_Model_RUN_Manager.py**
- Add functions to manage and use different sources in "nwp" mode (icon-2i)
- Fix bugs related to not initialized conditions

---

## Version 4.4.2 — *2023-11-27*
### APP: **RFARM_Model_RUN_Manager.py**
- Add functions to manage and use different sources in "nwp" mode (moloc)
- Fix methods to organize disaggregation domain with a buffer zone (to extend the starting grid)
- Fix methods to manage nans or zeros in the Fourier transform

---

## Version 4.4.1 — *2023-06-23*
### APP: **RFARM_Model_RUN_Manager.py**
- Add functions to manage different sources in "nwp" mode (lami-2i, ecmwf0100, gfs0.25)

---

## Version 4.4.0 — *2023-05-26*
### APP: **RFARM_Model_RUN_Manager.py**
- Add and refactor functions for "expert_forecast" and "nwp" mode
- Add the "pth" parameter to the "expert_forecast" mode
- Add the boundary limits to the RainFARM fields in the "expert_forecast" mode
- Support for different time window in the "expert_forecast" mode

### APP: **RFARM_Converter_EF_Rain.py**
- Add procedures for Marche and Liguria regions to convert "expert_forecast" original files to the format required by RainFARM

### FIX: **RFARM_Model_RUN_Manager.py**
- Fix bugs in the grid definition
- Fix bugs in the "expert_forecast" output fields
- Fix function and method issues
- Fix library version compatibility problems

---

## Version 4.3.0 — *2022-12-19*
### APP: **RFARM_Model_RUN_Manager.py**
- RainFARM package refactored from HyDe package (previous versions)

### FIX: **RFARM_Model_RUN_Manager.py**
- Fix metagauss domain issue related to expert-forecast mode

---

## Version 4.2.1 — *2021-08-01*
### FIX: **HYDE_Model_RFarm.py**
- Add resampling of output datasets according to time delta settings
- Fix 3D array vertical direction for ecmwf0100 NWP
- Fix time shift issue for lami-2i NWP
- Fix geographical position issue for irregular grid of lami-2i NWP

---

## Version 4.2.0 — *2021-05-03*
### APP: **HYDE_Model_RFarm.py**
- Add support for GFS 0.25 products
- Update WRF reader for supporting 3D files

---

## Version 4.1.0 — *2021-02-02*
### APP: **HYDE_Model_RFarm.py**
- Add expert forecast routine

### FIX: **HYDE_Model_RFarm.py**
- Adapt scripts and fix bugs

---

## Version 4.0.3 — *2021-01-25*
### FIX: **HYDE_Model_RFarm.py**
- Adapt scripts and fix bugs

---

## Version 4.0.2 — *2021-01-13*
### FIX: **HYDE_Model_RFarm.py**
- Fix and update data readers for lami-2i and ecmwf0100 NWP
- Fix and update data readers for expert forecast datasets
- Fix and update model application for expert forecast case

---

## Version 4.0.1 — *2020-05-22*
### FIX: **HYDE_Model_RFarm.py**
- Fix incorrect slope in time estimation
- Minor bugs in bash and Python scripts

---

## Version 4.0.0 — *2019-09-02*
### APP: **HYDE_Model_RFarm.py**
- Beta release for HyDe package and Python 3.x

### DRV: **lib_rfarm_core.py**
- Beta release for HyDe package for RainFARM model

---

## Version 3.5.2 — *2018-09-10*
### APP: **FP_Model_RainFarm.py**
- Beta release for FloodProofs library

---

## Version 3.5.1 — *2017-11-14*
### FIX: **FP_Model_RainFarm.py**
- Fix bugs (accumulated and instantaneous rain)

---

## Version 3.5.0 — *2017-05-30*
### APP: **FP_Model_RainFarm.py**
- Refactoring of code to update libraries and applications (pandas and xarray)

---

## Version 3.0.1 — *2015-09-24*
### APP: **FP_Model_RainFarm.py**
- Final release for operational chain mode (Regional Operational Chain)

---

## Version 3.0.0 — *2015-08-23*
### APP: **FP_Model_RainFarm.py**
- Final release for experimental project (DRIHM)

---

## Version 2.0.1 — *2014-04-08*
### APP: **RainFarm.py**
- Final version for experimental mode (Python version based on RainFarm MATLAB implementation)

---

## Version 2.0.0 — *2014-01-22*
### APP: **RainFarm.py**
- Beta version in Python language based on RainFarm MATLAB version

---

## Version 1.0.0 — *2013-05-02*
### APP: **RainFarm.m**
- Initial MATLAB version of RainFARM
