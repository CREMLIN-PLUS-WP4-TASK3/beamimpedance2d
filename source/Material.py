"""
FEniCS 2D beam coupling impedance simulation in frequency domain

by Uwe Niedermayer 2014
"""

import numpy
import pylab


class MaterialProperties(object):
    def __init__(self, nur, nui, epsilon, kappa):
        self.Nur=nur
        self.Nui=nui
        self.Eps=epsilon
        self.Kappa=kappa

    def reluctivityUpdate(self,nur,nui):
        print("old reluctivity: ", self.Nur, " +i ",self.Nui)
        self.Nur=nur
        self.Nui=nui
        print("new reluctivity: ", self.Nur, " +i ",self.Nui)



###############################################
## Read dispersion data list
def PermeabilityRead():
    try:
        infilename = FerriteDataFile
        ifile = open( infilename, 'r')  # open file for reading
    except:
        raise FileNotFoundError("File not found or cannot be opened!")

    murArray=[]
    muiArray=[]
    fArray=[]

    re=0
    im=0
    fc=0

    for line in ifile:
        linetuple = line.split()
        try:
            fc=float(linetuple[0])
            re=float(linetuple[1])
            im=float(linetuple[2])
            fArray.append(fc)
            murArray.append(re)
            muiArray.append(im)
        except:
            print ("Line ommitted")

    return [fArray,murArray,muiArray]
###############################################################

##############################################################
# A linear interpolator to enable evaluation of material parameters at arbitrary frequency
def ReluctivityInterpolate(fd,fArray,murArray,muiArray):
    #fd=f[0]#desired f   NOOO!
    print ('Interpolating permemeability for f= ', fd, 'Hz')
    #print ('Length of fArray', len(fArray))
    #print (fArray)
    mur=0
    mui=0
    nur=0
    nui=0
    n=0
    notfound=True

    while notfound:
        fc=fArray[n]        #current f
        fn= fArray[n+1]     #next f
        #print (fc, ' ' , fd, ' ' ,fn)
        if fc <=fd and fd<= fn:
            mur=murArray[n]+(murArray[n+1]-murArray[n])/(fn-fc) *(fd-fc)
            mui=muiArray[n]+(muiArray[n+1]-muiArray[n])/(fn-fc) *(fd-fc)
            notfound=False
            #print ('great success!')

        n=n+1

        if n>len(fArray)-2:
            notfound=False
            raise ValueError("Problem with frequency range")

    print ('mu= ', mur, ' -i ',mui)

    if mur==0 and mui==0:
        print ("mu value not found!!!")
        mur=1
        mui=0

    nur=mur/(mur**2 + mui**2)
    nui=mui/(mur**2 + mui**2)

    return [nur,nui]
###################################################################################

def MaterialTest():
    [fArray,murArray,muiArray]=PermeabilityRead()

    nur=[]
    nui=[]
    f2Ar=[]

    print (len(fArray))

    for n in range(len(fArray)-10):
        print (n)
        f2Ar.append(fArray[n]*1.1)
        [nurttt,nuittt]=ReluctivityInterpolate(f2Ar[n],fArray,murArray,muiArray)
        nur.append(nurttt)
        nui.append(nuittt)

    fdysdf=pylab.figure(figsize=(20,12))
    pylab.loglog(fArray,murArray,fArray,muiArray)
    ffdfd=pylab.figure(figsize=(20,12))
    pylab.loglog(f2Ar,nur,f2Ar,nui)
    pylab.show()
    return
