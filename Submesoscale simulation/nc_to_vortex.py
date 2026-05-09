from netCDF4 import *
from pathlib import Path
import xarray as xr
import numpy as np
import math
import cv2
import json
import sys

trange = range(0)
is_directory = None
file = ""
 

def vortex_detection() :
    # detect Okubo Weiss contours with OpenCV and filter large ones and noise
    
    path = Path(r"C:\\Users\\albou\\Desktop\\The_Ocean_Cleanup\\Work\\Github\\Submesoscale simulation\\"+file)

    


    global is_directory
    is_directory = not(path.is_file())

    global trange

    if is_directory:

        ux = xr.open_dataset(path.joinpath("ux.nc"), engine="netcdf4")
        uy = xr.open_dataset(path.joinpath("uy.nc"), engine="netcdf4")

        trange = range(len(ux["time"]))  

        #Grid spacing in meters
        dx = np.abs(np.mean(np.diff(ux["x"])))
        dy = np.abs(np.mean(np.diff(ux["y"])))

    else :  
        dsr = xr.open_dataset(path, engine="netcdf4")

        mask_loc =  maskloc = (dsr["longitude"] >= -145) & (dsr["longitude"] <= -135) & (dsr["latitude"] >= 25) & (dsr["latitude"] <= 35 )

        #apply the mask if needed
        if mask_loc is not None:
            dsr = dsr.where(mask_loc, drop=True)

        trange = range(len(dsr["time"]))    

        #Grid spacing in geo degrees
        dx_deg = np.abs(np.mean(np.diff(dsr["longitude"])))
        dy_deg = np.abs(np.mean(np.diff(dsr["latitude"])))

        # Earth radius in meters
        R_earth = 6371e3  

        #Computes the grid spacing in meters
        dy = float(np.deg2rad(dy_deg) * R_earth)
        dx = float(np.deg2rad(dx_deg) * R_earth * np.cos(np.deg2rad(dsr["latitude"].mean())))
    
    
    pos = []

    for t,time in enumerate(trange) :
        pos.append([])

        #if the dataset is from GLORYS, the depth dimension must be taken into account
        if is_directory :
            duy = (ux["ux"][time].shift(y=-1)-ux["ux"][time].shift(y=1))/(2*dx)
            dvx = (uy["uy"][time].shift(x=-1)-uy["uy"][time].shift(x=1))/(2*dy)
            dux = (ux["ux"][time].shift(x=-1)-ux["ux"][time].shift(x=1))/(2*dx)
            dvy = (uy["uy"][time].shift(y=-1)-uy["uy"][time].shift(y=1))/(2*dy)
        elif dsr["uo"].ndim == 4 :
            duy = (dsr["uo"][time,0].shift(latitude=-1)-dsr["uo"][time,0].shift(latitude=1))/(2*dx)
            dvx = (dsr["vo"][time,0].shift(longitude=-1)-dsr["vo"][time,0].shift(longitude=1))/(2*dy)
            dux = (dsr["uo"][time,0].shift(longitude=-1)-dsr["uo"][time,0].shift(longitude=1))/(2*dx)
            dvy = (dsr["vo"][time,0].shift(latitude=-1)-dsr["vo"][time,0].shift(latitude=1))/(2*dy)
        else : 
            duy = (dsr["uo"][time].shift(latitude=-1)-dsr["uo"][time].shift(latitude=1))/(2*dx)
            dvx = (dsr["vo"][time].shift(longitude=-1)-dsr["vo"][time].shift(longitude=1))/(2*dy)
            dux = (dsr["uo"][time].shift(longitude=-1)-dsr["uo"][time].shift(longitude=1))/(2*dx)
            dvy = (dsr["vo"][time].shift(latitude=-1)-dsr["vo"][time].shift(latitude=1))/(2*dy)
                    
        #compute the Okubo-Weiss parameter
        owp = (dux - dvy)**2 + (dvx + duy)**2 - (dvx - duy)**2

        #filter only positive Okubo Weiss parameter values
        owp[:,:] = np.where(owp[:,:] > 0, 0 , owp[:,:])

        #convert the Okubo Weiss values to shades of grey
        owp_fit_gray = np.array(owp[:,:]*255/np.min(owp[:,:]), dtype=np.uint8)
        
        owp_rgb = cv2.cvtColor(owp_fit_gray, cv2.COLOR_GRAY2RGB)

        #threshold the Okubo Weiss gray values
        _,owp_tresh = cv2.threshold(owp_fit_gray,5, 255, cv2.THRESH_BINARY)
        
        #detect the contours
        contours,_ = cv2.findContours(owp_tresh, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)

        contours_filtered = []
        if is_directory:
            height = ux["y"].shape[0]
            width = ux["x"].shape[0]
        else :
            height = dsr["latitude"].shape[0]
            width = dsr["longitude"].shape[0]

        #contour filtering
        for c,cnt in enumerate(contours) :
            #get the size of the rectangle bounding the contour
            rect = cv2.boundingRect(cnt)

            #filter noise (value is the maximum size of the noise contours)
            size_noise = 0.01
            
            #filter the noise and the large contours like the border of the domain
            if float(rect[2]/width) > size_noise and float(rect[3]/height) > size_noise and float(rect[2]/width) < 0.5 and float(rect[3]/height) < 0.5 :
                
                contours_filtered.append(cnt)

                #register the position of the detected vortex
                if is_directory:
                    pos[t].append([float(ux["y"][int(rect[1]+rect[3]/2)]),float(ux["x"][int(rect[0]+rect[2]/2)]),float(owp[int(rect[1]+rect[3]/2),int(rect[0]+rect[2]/2)])])
                else :    
                    pos[t].append([float(dsr["latitude"][int(rect[1]+rect[3]/2)]),float(dsr["longitude"][int(rect[0]+rect[2]/2)]),float(owp[int(rect[1]+rect[3]/2),int(rect[0]+rect[2]/2)])])

        cv2.drawContours(owp_rgb, contours, -1, (0, 255, 0), 1)
        cv2.imshow("Okubo Weiss contours", owp_rgb)
        cv2.waitKey(1)

    if is_directory:
        ux.close()
        uy.close()
    else :
        dsr.close()

    print("Detection done")
    return pos



def vortex_tracking(pos) :
    tracks = []

    #create tracks for each vortex detected at t = 0
    for vortex_num in range(len(pos[0])) :
        tracks.append([pos[0][vortex_num][0:2]])

    #loop remaining time
    for t in range(1,len(trange)) :
        print("Tracking in progress... "+str(int(t/len(trange)*100))+"%",end="\r")
        #range of the tracks currently registered
        former = range(len(tracks))
        
        #range of the vortexes to be assigned to tracks by next timestep
        current = range(len(pos[t]))

        #loop vortexes to be assigned
        for i in current :

            #find the closest last registered vortex among the existing tracks only if distance is less than 0.1 degrees (~11 km)
            min_dist = 100000000
            min_j = -1
            for j in former :
                if tracks[j][-1][0] is not None and tracks[j][-1][1] is not None :
                    dist = math.sqrt((pos[t][i][0]-tracks[j][t-1][0])**2+(pos[t][i][1]-tracks[j][t-1][1])**2)
                    if dist < min_dist and (dist < 11000 and is_directory) or (dist < 0.1 and not(is_directory)) :
                        min_dist = dist
                        min_j = j
            #if a minimum is found : extend the track and remove the track from the available tracks list      
            if min_j != -1 :
                tracks[min_j].append(pos[t][i][0:2])
                former = [j for j in former if j != min_j]
            #if no track exists to connect the vortex to : create a new one  
            else :
                tracks.append([[None,None] for k in range(t)] + [pos[t][i][0:2]])

        #if tracks are not extended they are discontinued by attributing a default position value        
        for i in former :
            tracks[i].append([None,None])
    print("Tracking done")
    return tracks



def filter_tracks(tracks,minimum_length_ratio = 0.02) :
    #filter tracks that last more than "minimum_length" timesteps
    tracks_filtered = []
    for track in tracks :
        if len(track)-track.count([None,None]) > minimum_length_ratio * len(trange) :
            tracks_filtered.append(track)
    print("Filtering done")
    return tracks_filtered



def register_tracks(tracks, suffix = "", domain_coords = None) :

    #Convert cartesian coordinates to geo degrees if a domain geocoordinates are provided
    # domain_coords should be an array of [longitude,latitude] from the center of the domain
    if domain_coords is not None :
        for i in range(len(tracks)) :
            for j in range(len(tracks[i])) :
                if tracks[i][j][0] is not None and tracks[i][j][1] is not None :
                    #convert the cartesian coordinates to geo degrees
                    lat = domain_coords[1] + (tracks[i][j][0]/110000)
                    lon = domain_coords[0] + (tracks[i][j][1]/(110000*np.cos(np.deg2rad(lat))))
                    tracks[i][j] = [lat,lon]
                else :
                    tracks[i][j] = [None,None]

    #Extract the name of the file from the local path
    register_name = file.split("\\")[-1]

    #add an underscore if there is a suffix
    if suffix != "":
        suffix = "_"+suffix

    with open("vortexes\\Tracks_"+str(register_name)+suffix+".json", "w") as f:
        json.dump(tracks, f, indent=4)
        print("Registering done")


"""
def circulation(ds,N,W,radius,time):
    #compute the circulation around a circle of given radius (in km) centered on (N,W) at given time index

    if W > 0 :
        W -= 360
    R = 6.371e6 #radius of the Earth
    pR = 2*np.pi*R

    #import velocity field
    u = np.array(ds["uo"])
    v = np.array(ds["vo"])
        
    #define number of points for circulation computation
    pts = 40

    #create coords for circle of the defined radius
    xcircle = np.cos(np.linspace(0,6.28,pts))*(radius*1000) 
    ycircle = np.sin(np.linspace(0,6.28,pts))*(radius*1000)

    #transform the circle to focus on the correct point
    ycircle /= pR/360
    ycircle += np.array(N)
    xcircle /= pR*np.cos(np.deg2rad(ycircle))/360
    xcircle += np.array(W) 

    #define interpolation coordinates
    Xc = np.linspace(np.min(ds["longitude"]),np.max(ds["longitude"]),u.shape[3])
    Yc = np.linspace(np.min(ds["latitude"]),np.max(ds["latitude"]),u.shape[2])

    #extract velocity field only with type float64
    ui = u[time,0,:,:].astype("float64")
    vi = v[time,0,:,:].astype("float64")

    #create interpolator
    fu = interp.RegularGridInterpolator((Xc,Yc), ui.T)
    fv = interp.RegularGridInterpolator((Xc,Yc), vi.T)

    circulation = 0

    #integrate u.dl over the circle
    for i in range(1,len(xcircle)) :
        
        if xcircle[i] < np.min(Xc) or xcircle[i] > np.max(Xc) or ycircle[i] < np.min(Yc) or ycircle[i] > np.max(Yc):
            return 0,0

        #circle length vector
        dlx = (xcircle[i] - xcircle[i-1])*pR*np.cos(np.deg2rad(ycircle[i]))/360
        dly = (ycircle[i] - ycircle[i-1])*pR/360

        #velocity vector
        ul = fu((xcircle[i],ycircle[i]))
        vl = fv((xcircle[i],ycircle[i]))

        #scalar product
        circulation += ul*dlx + vl*dly

    #tangent speed computation
    tg_spd = circulation/(2*3.1415*radius*1000)
    
    return circulation,tg_spd


def lamb_oseen_fit(spd_profile,circ_profile):
    #fit a lamb-oseen profile to the given speed and circulation profiles (1D arrays)
    spd_argmax = np.argmax(np.abs(spd_profile))
    rmax = 1+spd_argmax*2
    circ = circ_profile[spd_argmax]/0.71
    return rmax,circ


def flat(array) :
    array_ret = []
    for i in range(len(array)):
        array_ret.extend(array[i])
    return array_ret


def circulation_profile(dsr,tracks,trange) :
    
    #velocity profile
    spd_profile = []

    #circulation profile
    circulation_profile = []    
    coords = []
    rmax = []   
    circ = []

    for i in range(len(tracks)):
        for j in range(1,100,2) :
            if not(tracks[i][0][0] == 0 and tracks[i][0][1] == 0) :
                result = circulation(dsr,tracks[i][0][0],tracks[i][0][1],radius=j,time=trange[0])

                #increase the radius until the circulation starts decreasing or hit the edge of the domain
                if np.isnan(result[0])  or (len(spd_profile) > 0 and np.abs(spd_profile[-1]) >np.abs(result[1])) :

                    #extraction of the parameters
                    rmaxi,circi = lamb_oseen_fit(spd_profile,circulation_profile)
                    rmax.append(rmaxi)
                    circ.append(circi)
                    coords.append([tracks[i][0][0],tracks[i][0][1],j])
                    spd_profile = []
                    circulation_profile = []
                    break
              
            spd_profile.append(result[1])
            circulation_profile.append(result[0])
    print("Circulation profile done")
    return rmax,circ,coords


def register_IC(filename, circ , rmax , coords , long, lat):
    
    #check if path exists
    if not(os.path.exists("LO_IC_files/"+filename)):
        
        #create file
        with open("LO_IC_files/"+filename, "a") as f:
            for i in range(len(rmax)):

                #convert geo2D coords into domain coords
                coordsx = int((float(coords[i][1]) - long)*110000)
                coordsy = int((float(coords[i][0]) - lat)*110000)

                #write the coords, core radius and circulation of the Lamb-Oseen vortex to be used as initial conditions
                f.write(str(coordsx)+" "+str(coordsy)+" "+str(round(rmax[i]*1000/math.sqrt(1.25643),0))+" "+str(round(circ[i],0))+"\n")
        print("IC registered")
"""



if __name__ == "__main__" :
    if len(sys.argv) > 1 :
        file = str(sys.argv[1])
    
        if len(sys.argv) > 2 :
            suf = str(sys.argv[2])
        else :
            suf = ""

        if len(sys.argv) > 3 :
            #Domain geocoordinates for converting cartesian to degrees 
            domain = np.array(sys.argv[3].split(","))
            domain = [float(i) for i in domain]
        else :
            domain = None


        #vortex position detection algorithm
        pos = vortex_detection()

        #vortex tracking algorithm
        tracks = vortex_tracking(pos)

        #vortex longevity filtering
        tracks_filtered = filter_tracks(tracks)

        #resulting tracks registration in a comparable file
        register_tracks(tracks_filtered, suffix = suf, domain_coords = domain)

    else :
        print("nc_to_vortex.py <file_path> <suffix_for_output_file : '' for none> <domain_geocoordinates_for_cartesian_to_degrees_conversion : 'longitude,latitude'>")

##### Lamb-Oseen IC extraction #####

# rmax,circ,coords = circulation_profile(dsr,tracks_filtered,trange)

# register_IC("IC_GLORYS_aug2022", circ , rmax , coords , long = domain_coords[0], lat = domain_coords[1])

