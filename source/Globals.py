'''
Created on Dec 5, 2013

@author: niedermayer
'''
# constants: (SI units)
import math
import os
import __builtin__

#aaa=Constant(123.34)

###Constants
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
__builtin__.logscale=True
__builtin__.dispersive=True
__builtin__.plot3Dflag=False
__builtin__.dipole=False
__builtin__.horizontal=True
__builtin__.dataexport=True
__builtin__.wallcurrent=False    #compute and export wall current
__builtin__.SIBC=False







