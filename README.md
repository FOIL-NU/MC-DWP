# MC-DWP: A Monte-Carlo simulation package of Spectroscopic Single-Molecule Localization Microscopy
The MC-DWP package is a comprehensive suite of simulation tools for optimizing and understanding spectroscopic single-molecule localization microscopy (sSMLM) experiments.

## Main Function
The main function of the package is `ssmlm_sim`, which takes in parameters related to the experimental setup, the microscope, the camera, and the spectrometer, and outputs the fitted data of the localization and spectral information of the simulation.

## Included Functions
This package includes several helper functions:

- `checksavefile`: Checks if a savefile exists for the current simulation settings.
- `dwp_module`: Calculates the diffraction angle after passing through the DWP module in radians.
- `gaussian2dcirc`: Returns a symmetric 2D Gaussian function with offset.
- `grating_data`: Returns the dispersion and efficiency of a given grating.
- `load_efficiency`: Loads the efficiency data for the spectrometer.
- `load_fluor`: Loads fluorophore data from a file.
- `load_spectral`: Loads the spectral information from the fpbase database.
- `mle_fn`: Implements the MLE cost function for the 2D Gaussian.
- `photon_budgets`: Calculates the background photon budget for a given splitting ratio.
- `psfsample`: Implements the sampling of the PSF onto the imaging sensor for a given set of photon positions.
- `sellmeier`: Implements the Sellmeier equation for glass.

## Getting Started

To use the ssmlm-sim package, clone the repository:

```bash
git clone https://github.com/FOIL-NU/MC-DWP.git
```

Then in your MATLAB script, you can use the ssmlm-sim functions. Here's an example using the main `ssmlm_sim` function:

```MATLAB
addpath(genpath('path_to_ssmlm_sim_package'));
% define your input parameters here
outdata = ssmlm_sim(photons, sim_settings, specmtr, spectral, camera, microscope, n_iter, rng_init);
```

Please replace `'path_to_ssmlm_sim_package'` with the actual path to the cloned ssmlm-sim package.

Please refer to individual function documentation for usage instructions.

See `example.m` for an example script to start out.

## License

This project is licensed under the terms of the CC BY-NC license.

## Contact

If you have any questions, please reach out to Wei Hong Yeo from [Functional Optical Imaging Laboratory](foil.northwestern.edu), Northwestern University, 2023.