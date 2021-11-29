# generated using CFL simulator on Nov. 29, 2021, 
# Hualin Xiao
#~/FHNW/STIX/SolarFlareAnalysis/grids/pySimulator$ python apps/simulate_cfl.py 

import numpy as np
from pathlib import Path

filename = Path(__file__).resolve().parents[1] / 'data' / 'skylut.npz'
data= np.load(filename)
lut=data['lut']
X=data['x']
Y=data['y']
#the lut is a numpy array with dimensions ns x  14, ns is the number of random points
#2 of 14 elements are flare location in arsec and the other 12 elements are pixel area of unshaded area

