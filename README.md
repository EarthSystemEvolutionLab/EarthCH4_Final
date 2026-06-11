# EarthCH4_Final
A revised CH4 wetland emission model based on Beerling et al. (2009), with atmospheric chemistry scripts for Photochem (Wogan et al. 2023, 2024, 2025) to compute pCH4 over the last 375 Ma.

OVERVIEW
This model system has multiple components, packaged as subfolders: (C1, EarthCH4_Emission_Model) a CH4 emissions model, primarily focused on fluxes from "coal wetlands" since the end-Devonian, (C2, Photochemical_Scripts) photochemical scripts using the 1D Photochem model (Wogan et al. 2022, 2025) to compute steady-state pCH4 given the emissions calculated in C1, and (C3, CH4_ClimateCalcs_Plots) scripts to compute and plot radiative forcing and climate warming from CH4 over the last 375 Ma, as well as to generate other relevant plots from the manuscript. We also include our nominal model output data files in these folders. No non-standard hardware is required (i.e., laptops or desktop computers are fine).

SYSTEM REQUIREMENTS AND INSTALLATION
You should be able to run the Matlab scripts (we used R2023a herein) without special toolboxes. To use the Photochem scripts, you will need to install Photochem in a Python environment (see https://github.com/Nicholaswogan/photochem for details; conda-forge works well). We use Photochem version 0.8.4 in our work. In addition, make sure that the following packages are installed in your Python environment: numpy, scipy.io, matplotlib, and IPython.display. Total installation of these packages in an environment should be relatively brief (>30 minutes may signal a problem somewhere).

INSTRUCTIONS FOR DEMO AND USE
To run the Matlab scripts, simply click "Run" in the Editor tab of the Matlab software. We recommend keeping all input data files and model scripts together in the same folder, though they are packaged in separate folders for sorting purposes. All Matlab scripts (for the emissions model, climate calculations, and figure generation) should run within 10 seconds or less.

To run the Python scripts, activate your environment ("conda activate [env]", where [env] is your Python environment with Photochem installed) and enter the command "python [script_name].py > [script_name]_output.txt". This will run the relevant script and save key outputs as .mat files while printing secondary outputs (e.g., proof of successful model convergence) to a .txt file. Since the model runs over 375 million years at 1 myr resolution, it may take over 24 hours for a full time-series computation (>48 hours may be cause for concern). You can run these on a laptop (up to 2 at a time) but using several tmux seessions on a Linux cluster is helpful if running multiple different scenarios.

You can demo the Matlab figure generation scripts using the nominal data files supplied in the CH4_ClimateCalcs_Plots/pCH4_PhotochemOutputs folder, or run the emissions model to feed through the photochemical scripts to regenerate the nominal data files given. To demo the photochemical codes, you can test-run one of the CH4_Holmes...Vdep.py scripts, as these are relatively quick to complete (should be <6 hours) compared to a full 375 Ma model time-series.

Please feel free to contact the developers (John Herring: herring3@iastate.edu or Ben Johnson: bwj@iastate.edu) with any questions or to request alternative or updated versions of any files (e.g., Jupyter Notebook versions of the Python scripts).

Note: core code and input files for Photochem scripts herein were shared by Nicholas Wogan (https://github.com/Nicholaswogan/photochem). We also modify the geotimescale.m file created by Ben Mills (https://github.com/bjwmills/geotimescale) for use in our plots.
