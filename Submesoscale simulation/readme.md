## Spectrum extraction
- nc_to_EDS.py create standardised csv files from netCDF files containing the EDS of each timestep.

- The standardised spectrum files are stored in the spectrums file 

- Resulting files can be opened to be compared in the spectrum_comparison.ipynb file.

## Vortex tracking : WIP
- nc_to_vortex.py create a standardised json file containing the path of the detected vortex in a netCDF file. The method used to detect the vortex require some tuning of a treshold depending on the flow characteristics and the algorithm to sort the vortex positions into tracks is very slow.

- Resulting files can be opened to be  compared in the vortex_comparison.ipynb file.

## Measures
- simulation_measures.ipynb computes the energy and enstrophy evolution of a given EDS file.





