# EarthCH4_Final
A revised CH4 wetland emission model based on Beerling et al. (2009), with atmospheric chemistry scripts for Photochem (Wogan et al. 2023, 2024, 2025) to compute pCH4 over the last 375 Ma.

This model system has multiple components: (C1) a CH4 emissions model, primarily focused on fluxes from "coal wetlands" since the end-Devonian, (C2) photochemical scripts using the 1D Photochem model (Wogan et al. 2022, 2025) to compute steady-state pCH4 given the emissions calculated in C1, and (C3) scripts to compute and plot radiative forcing and climate warming from CH4 over the last 375 Ma, as well as to generate other relevant plots from the manuscript. Each of these modules is packaged in a subfolder above.

You should be able to run the Matlab scripts (we used R2023a herein) without special toolboxes. To use the Photochem scripts, you will need to install Photochem in a Python environment (see https://github.com/Nicholaswogan/photochem for details; conda-forge works well). We use Photochem version 0.8.4 in our work.

To run the Matlab scripts, simply click "Run" in the Editor tab of the Matlab software. All Matlab scripts (for the emissions model, climate calculations, and figure generation) should run within 10 seconds or less. We recommend keeping all input data files and model scripts together in the same folder, though they are packaged in separate folders for sorting purposes.

To run the Python scripts, activate your environment (with Photochem installed) and enter the command "python [script_name].py > [script_name]_output.txt". This will run the relevant script and save key outputs as .mat files while printing secondary outputs (e.g., proof of successful model convergence) to a .txt file.

Please feel free to contact the developers (John Herring: herring3@iastate.edu or Ben Johnson: bwj@iastate.edu) with any questions or to request alternative or updated versions of any files (e.g., Jupyter Notebook versions of the Python scripts).
