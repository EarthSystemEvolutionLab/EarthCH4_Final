# EarthCH4_Final
A revised CH4 wetland emission model based on Beerling et al. (2009), with atmospheric chemistry scripts for Photochem (Wogan et al. 2023, 2024, 2025) to compute pCH4 over the last 375 Ma.

This model system has multiple components: (C1) a CH4 emissions model, primarily focused on fluxes from "coal wetlands" since the end-Devonian, (C2) photochemical scripts using the 1D Photochem model to compute steady-state pCH4 given the emissions calculated in C1, and (C3) scripts to compute and plot radiative forcing and climate warming from CH4 over the last 375 Ma, as well as to generate other relevant plots from the manuscript.

You should be able to run the Matlab scripts (R2023a or higher) without special toolboxes. To use the Photochem scripts, you will need to install Photochem in a Python environment (see https://github.com/Nicholaswogan/photochem for details; conda-forge works well). We use Photochem version 0.8.4 in our work.

Please feel free to contact the developers (John Herring: herring3@iastate.edu or Ben Johnson: bwj@iastate.edu) with any questions or to request alternative or updated versions of any files (e.g., Jupyter Notebook versions of the Python scripts).
