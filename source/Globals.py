"""
FEniCS 2D beam coupling impedance simulation in frequency domain

by Uwe Niedermayer 2014
"""


import math
import os
import builtins



###Constants (SI-units)
builtins.mu0=4*math.pi*10**-7
builtins.nu0=1.0/mu0
builtins.c0=299792458
builtins.eps0=1/(mu0*c0**2)
builtins.Z0=1/(eps0*c0)
builtins.I=complex(0,1)


###Approximation Orders
builtins.div_long_order=1
builtins.curl_long_order=1

builtins.curl_order=1
builtins.div_order=1


###Result folder
builtins.resultfolder="./SIS100-EmergencyKicker-RESULTS-roomtemp/"
try:
    os.system('mkdir '+ resultfolder)
except:
    print ("results directory already exists")

###Ferrite Data file
#builtins.FerriteDataFile="mischung43.dat"
builtins.FerriteDataFile="mu8c11.dat"

###Options
builtins.logscale=False     #log or lin frequency stepping
builtins.dispersive=False  #Only for dispersive material update is required, set false if nondispersive-->faster
builtins.plot3Dflag=False#Make all 3D plots (use only for one frequency point)
builtins.quadrupole=False       #Quadrupole transverse impedance (very preliminary!!!)
builtins.dipole=False       #Dipole transverse impedance
builtins.horizontal=True     # Only relevant for dipole=True
builtins.dataexport=True # Write (Overwrite) file in the result folder
builtins.wallcurrent=False    #compute and export wall current
builtins.SIBC=False           # Use surface impedance boundary condition. 
builtins.kappa_s=6e9 #cold copper         # Wall conductivity for single layer surface impedance
builtins.twolayer=False   





