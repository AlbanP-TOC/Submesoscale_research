from netCDF4 import *
import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd
import xarray as xr
import numpy as np
import sys

from matplotlib import rc
font = {"family": "serif", "weight": "normal", "size": 13}
rc("font", **font)
rc("text", usetex=True)


def EDS_fields(file, show = False, register = True, suffix = ""):
    """Converts a simulation output (series of velocity fields) into an energy density spectrum.
    The simulation can be a netCDF file containing all the fields in geo2D coordinates or a directory containing a netCDF file for each component in a cartesian grid.
    The result is registered into a csv file for further analysis.
    The input netCDF file must contain "uo" and "vo" variables. GLORYS files have a specific depth dimension that is handled by a flag. 
    A mask can also be applied to select an area of interest"""
    
    #Import the dataset
    path = Path(str(file))

    is_directory = path.is_dir()

    if is_directory:
        ux = xr.open_dataset(path.joinpath("ux.nc"), engine="netcdf4")
        uy = xr.open_dataset(path.joinpath("uy.nc"), engine="netcdf4")
        trange = range(len(ux["time"]))  

        #Definition of the size of the grid
        nx = len(ux["x"])
        ny = len(ux["y"])

        #Grid spacing in meters
        dx = np.abs(np.mean(np.diff(ux["x"])))
        dy = np.abs(np.mean(np.diff(ux["y"])))

    else :  
        dsr = xr.open_dataset(path, engine="netcdf4")
        trange = range(len(dsr["time"]))    

        #Definition of the size of the grid
        nx = len(dsr["longitude"])
        ny = len(dsr["latitude"])

        #Grid spacing in geo degrees
        dx_deg = np.abs(np.mean(np.diff(dsr["longitude"])))
        dy_deg = np.abs(np.mean(np.diff(dsr["latitude"])))

        # Earth radius in meters
        R_earth = 6371e3  

        #Computes the grid spacing in meters
        dy = float(np.deg2rad(dy_deg) * R_earth)
        dx = float(np.deg2rad(dx_deg) * R_earth * np.cos(np.deg2rad(dsr["latitude"].mean())))

    

        
    # Wavenumber axes
    kx = 2 * np.pi *np.fft.fftfreq(ny, d=dx)  
    ky = 2 * np.pi * np.fft.fftfreq(nx, d=dy) 

    # Radial wavenumber grid
    kx_grid, ky_grid = np.meshgrid(kx, ky, indexing='ij')
    k_mag = np.sqrt(kx_grid**2 + ky_grid**2)

    #Bins definition
    valid_k = k_mag > 0
    k_bins = np.logspace(np.log10(k_mag[valid_k].min())+1, np.log10(k_mag.max()), num=40)
    k_centers = (k_bins[:-1] * k_bins[1:])**0.5
   
    # Nyquist frequency mask
    nyq_mask = (k_centers <= 2*np.pi/(max(dx,dy)*2))

    dict_df = {"characteristic_length": 2*np.pi/k_centers[nyq_mask]}   

    print(list(dsr.keys()))

    #For each dataset
    for i in trange:
        # print(dsr.variables.names)
        # Adapt to the GLORYS format if needed
        if is_directory:
            u = ux["ux"][i]
            v = uy["uy"][i]
        elif "ugos" in list(dsr.keys()) :
            u = dsr["ugos"][i]
            v = dsr["vgos"][i]
        elif dsr["uo"].ndim == 4 : #GLORYS format with depth dimension
            u = dsr["uo"][i,0]
            v = dsr["vo"][i,0]
        else :
            u = dsr["uo"][i]
            v = dsr["vo"][i]

        # Compute Fourier transform of the velocity components
        E_uf = np.fft.fft2(u)
        E_vf = np.fft.fft2(v)

        # Compute power spectrum normalised to verify the Parseval's theorem
        E_f = 0.5 * (np.abs(E_uf)**2 + np.abs(E_vf)**2)/(nx*ny)**2       

        # Bin the energies   
        E_spectrum = np.zeros_like(k_centers)
        for k in range(len(k_bins)-1):
            mask_bin = (k_mag >= k_bins[k]) & (k_mag < k_bins[k+1])
            E_spectrum[k] = np.sum(E_f[mask_bin])

        #update the dictionnary for the data exportation
        dict_df.update({str(i) : E_spectrum[nyq_mask]})

        # Plot the spectrum if needed
        if show :
            ax = plt.gca()
            ax.plot(2*np.pi/(k_centers[nyq_mask]), E_spectrum[nyq_mask], '-',c = (0,0,0))
            
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_xlabel('Characteristic length (m) ')
            
            ax.set_ylabel('Energy Density Spectrum')
            ax.grid(visible= True)

    # Register the spectrum in a csv file if needed
    if register :
        
        df = pd.DataFrame(dict_df)

        #Extract the name of the file from the local path
        register_name = file.split("\\")[-1]

        #add an underscore if there is a suffix
        if suffix != "":
            suffix = "_"+suffix
        
        df.to_csv("spectrums\\EDS_"+str(register_name)+suffix+".csv",index = False)


if "__main__" == __name__ :
    if len(sys.argv) > 1 :  
     
        if len(sys.argv) > 2 :
            suf = str(sys.argv[2])
        else : 
            suf = ""

        EDS_fields(str(sys.argv[1]), show = False, register = True, suffix = suf)

    else : 
        print("nc_to_EDS.py <file_path> <suffix_for_output_file : '' for none>")