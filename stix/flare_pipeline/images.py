#!/usr/bin/env python3  
# -*- coding: utf-8 -*- 
#----------------------------------------------------------------------------
# Created By  : Hualin Xiao (hualin.xiao)
# Created Date: Sept 5, 2022
# version =1.0
"""
    various tools to manipulate stix images
"""

import os
import sys
import matplotlib
import numpy as np
from scipy import ndimage
from matplotlib import cm
from matplotlib import pyplot as plt
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
import sunpy
from sunpy.map import make_fitswcs_header
from sunpy.coordinates.frames import HeliocentricEarthEcliptic, HeliographicStonyhurst
