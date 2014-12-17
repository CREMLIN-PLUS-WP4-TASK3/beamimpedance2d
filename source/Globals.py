"""
FEniCS 2D beam coupling impedance simulation in frequency domain

by Uwe Niedermayer 2014
"""


import math
import os
import __builtin__



###Constants (SI-units)
__builtin__.mu0=4*math.pi*10**-7
__builtin__.nu0=1.0/mu0
__builtin__.c0=299792458
__builtin__.eps0=1/(mu0*c0**2)
__builtin__.Z0=1/(eps0*c0)
__builtin__.I=complex(0,1)


###Approximation Orders
__builtin__.div_long_order=1
__builtin__.curl_long_order=1

__builtin__.curl_order=1
__builtin__.div_order=1


###Result folder
__builtin__.resultfolder="./results/"
try:
    os.system('mkdir '+ resultfolder)
except:
    print ("results directory already exists")

###Options
__builtin__.logscale=True     #log or lin frequency stepping
__builtin__.dispersive=False  #Only for dispersive material update is required, set false if nondispersive-->faster
__builtin__.plot3Dflag=False  #Make all 3D plots (use only for one frequency point)
__builtin__.dipole=True       #Dipole transverse impedance
__builtin__.horizontal=True     # Only relevant for dipole=True
__builtin__.dataexport=False    # Write (Overwrite) file in the result folder
__builtin__.wallcurrent=False    #compute and export wall current
__builtin__.SIBC=True           # Use surface impedance boundary condition. Set conductivity in Curlsolver.py.
__builtin__.kappa_s=1e6         # Wall conductivity for surface impedance






