# RainFARM (rfarm)

RainFARM is a stochastic downscaling method for precipitation fields.  
This repository provides a Python-based application and utilities to prepare environments and run workflows around the RainFARM algorithm.

> **License:** EUPL-1.2 (see [`LICENSE.md`](./LICENSE.md)).

---

## What is RainFARM?

RainFARM (Rainfall Filtered Autoregressive Model) downscales coarse precipitation fields by extrapolating their Fourier power spectrum to smaller scales and adding stochastic small-scale variability, followed by a nonlinear transformation to obtain realistic rainfall fields.  
The method is widely used in hydrological and climate modeling for generating high-resolution rainfall fields consistent with large-scale forcings.

---

## Repository Layout

- `rfarm/` — Python package source  
- `apps/` — Application drivers and workflow wrappers  
- `bin/` — Helper scripts and CLI utilities  
- `tools/` — Additional utilities or configuration tools  
- Documentation: `README.md`, `AUTHORS.md`, `CHANGELOG.md`, `CODEOWNERS`, `LICENSE.md`

---

## Requirements

- Linux (Debian/Ubuntu recommended)  
- Python 3.x  
- Standard scientific Python stack (NumPy, SciPy, etc.)  
- Optional: NetCDF utilities, QGIS, R, or other hydrometeorological tools for extended use

---

## Quick Start (Environment Setup)

This repository includes a setup script for automatic environment configuration:

```bash
# from the repository root
bash setup_rfarm_system_conda_python.sh
```

This script creates a Conda-based Python environment and installs the dependencies required by the `rfarm` application.

> After setup, activate the environment as instructed in the script output.

---

## Installation (Editable Mode)

Once your environment is active, install the package in editable mode:

```bash
pip install -e .
```

This allows you to modify and test the code without reinstalling after each change.

---

## Citing RainFARM

If you use RainFARM in scientific work, please cite the original publications describing the algorithm and its applications in hydrometeorological modeling and climate downscaling.

---

## Contributing

Contributions are welcome! Please:

1. Fork this repository and create a new branch  
2. Follow existing code style and documentation conventions  
3. Submit a Pull Request with a clear description of your change

See also:  
- [`AUTHORS.md`](./AUTHORS.md)  
- [`CHANGELOG.md`](./CHANGELOG.md)

---

## Related Links & References

- [RainFARM GitHub Repository](https://github.com/c-hydro/rfarm)  
- Related projects by [CIMA Research Foundation](https://github.com/c-hydro)  
- Method overview and community examples (e.g., *pysteps* RainFARM implementation and publications on stochastic downscaling)

---

## License

This project is licensed under the **European Union Public License 1.2 (EUPL-1.2)**.  
See [`LICENSE.md`](./LICENSE.md) for full terms.

---

### Notes

- This README intentionally avoids detailed code examples (as requested).  
- You can extract concise usage snippets directly from the included `.sh` scripts if you wish to document them later.
