import pandas as pd
from pathlib import Path
from math import *

def ReadParticles(prefix, folder):
    fname = f"{prefix}"
    args = {"header" : None, "engine" : "c"}

    ds = pd.read_table(folder.joinpath(f"{fname}_particle_velocity_x.dat"), **args)
    vx = pd.read_table(folder.joinpath(f"{fname}_particle_velocity_x.dat"), **args)
    vy = pd.read_table(folder.joinpath(f"{fname}_particle_velocity_y.dat"), **args)
    px = pd.read_table(folder.joinpath(f"{fname}_particle_position_x.dat"), **args)
    py = pd.read_table(folder.joinpath(f"{fname}_particle_position_y.dat"), **args)
    npart = ds.shape[1] - 1
    
    return [Particles(vx, vy, px, py, ip) for ip in range(npart)]
    
class Particles(object):
    def __init__(self, vx, vy, px, py, ip):
        self.vx = vx.values[:,ip + 1]
        self.vy = vy.values[:,ip + 1]
        self.px = px.values[:,ip + 1]
        self.py = py.values[:,ip + 1]
        self.time = py.values[:,0]