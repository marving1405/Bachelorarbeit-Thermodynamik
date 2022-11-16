"""
first link of refprop for python.
Documentation available: https://refprop-docs.readthedocs.io/en/latest/DLL/high_level.html#f/_/SETPATHdll
Tutorial: https://nbviewer.org/github/usnistgov/REFPROP-wrappers/blob/master/wrappers/python/notebooks/Tutorial.ipynb
(set-up extracted from the latter)

author: Alexandra Welp
"""


# Import the main class from the Python library
from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary
import os
# Imports from the standard library
import glob

# Imports from conda-installable packages
import pandas

# Import numpy
import numpy as np

# Import matplotlib for plotting
import matplotlib.pyplot as plt

# Now we instantiate the library, and use the environment variable to
# explicitly state which path we want to use. It was decided to make
# the path handling explicit (though more verbose), because explicit
# is almost always better than implicit
RP = REFPROPFunctionLibrary(os.environ['RPPREFIX'])

# Now we tell REFPROP what the root directory is that it should use.  This root directory should contain, at least:
# A) REFPRP64.DLL (or REFPROP.dll for 32-bit windows, or librefprop.so or librefprop.dylib, for linux or OSX respectively)
# B) FLUIDS folder (case sensitive)
# C) MIXTURES folder (case sensitive)
RP.SETPATHdll(os.environ['RPPREFIX'])

# Get the unit system we want to use (we will revisit this GETENUM function later)
MOLAR_BASE_SI = RP.GETENUMdll(0, "MOLAR BASE SI").iEnum

# first example
# input of variables used for later property call
p_Pa = 101.325
Q = 0

# Setup fluid/fluid mixture
RP.SETFLUIDSdll("Water")

# Example for 2 different functions that can be used in 2 phase region
l = RP.ABFLSHdll("PQ", p_Pa, Q, [1.0], 000) # faster, easier but limited properties available
print(l)

r = RP.REFPROPdll("","PQ","T;D;h",MOLAR_BASE_SI,0,0,p_Pa*10**3,Q,[1.0]).Output[0:3] # slower, complicated, but huge amount of properties
print(r)

