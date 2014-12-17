"""
FEniCS 2D beam coupling impedance simulation in frequency domain

by Uwe Niedermayer 2014
"""

from dolfin import *
import numpy
import pylab
import os

#####################################################
##Some nice values for the plots
kappa_value=10**6
h=0.041
b=0.04
length=1.0
R_Ohm=length /(kappa_value*pi*(h**2-b**2))
f_skin=1.0 /(pi*mu0*kappa_value*(h-b)**2)
#####################################################

def Zlong(Elr,Eli,Jszr,Jszi,q):
    Jzrfct=Jszr
    Jzifct=Jszi
    Zl=(-1.0/q**2)*(assemble(inner(Elr,Jzrfct)*dx)+0*assemble(inner(Eli,Jzifct)*dx) \
                +I*(-0*assemble(inner(Elr,Jzifct)*dx)+assemble(inner(Eli,Jzrfct)*dx) ))
    return Zl


def Ztrans(Elr,Eli,Jszr,Jszi, omega, beta,xd,yd):
    Jzrfct=Jszr
    Jzifct=Jszi
    
    if horizontal:
        DipTest=Expression('x[0]-xd',xd=xd)
    else:
        DipTest=Expression('x[1]-yd',yd=yd)
    
    #x=Expression('x[0]')
    DM=assemble(DipTest*Jzrfct*dx)
    print "DM: ", DM
    Ztr=-beta*c0/(omega*DM**2) *(assemble(inner(Elr,Jzrfct)*dx)+0*assemble(inner(Eli,Jzifct)*dx) \
                +I*(-0*assemble(inner(Elr,Jzifct)*dx)+assemble(inner(Eli,Jzrfct)*dx) ))
    #check me! I'm just the power loss integral!
    return Ztr



def CplxImpExport(filename,f,Z):
    completeName=os.path.join(resultfolder,filename)
    Zfile=open(completeName,'w')
    for j in range(len(f)):
        Zfile.write(str(f[j])+"\t" +str(Z[j].real) +"\t " + str(Z[j].imag) + "\n")
    Zfile.close()
    return None

def PlotZlong(f,Zlong,Zlongloss,Zsc_ana):
    fig=pylab.figure(figsize=(20,12))
    pylab.suptitle('Longitudinal Impedance',fontsize=30)
    
    RElog=fig.add_subplot(221)
    try:
        pylab.loglog(f,numpy.fabs(numpy.asarray(Zlong).real),f,numpy.fabs(numpy.asarray(Zlongloss).real), linewidth=4 )
    except:
        print ("Zero loglogplot")
        #pylab.axhline(y=R_Ohm, color='k')
    pylab.axvline(x=f_skin, linewidth=4, color='k')
    pylab.axhline(y=R_Ohm, linewidth=4, color='k')
        ## legend
    pylab.legend((r'FEM real',r'Zloss real'), shadow = True, fontsize=25)

    ltext = pylab.gca().get_legend().get_texts()

    pylab.xlabel(r'$f \;$ [Hz]', fontsize=25)
    pylab.ylabel(r'Re $\{Z_\parallel \} [\Omega]$', fontsize=25)
    #pylab.grid()
    for label in RElog.get_xticklabels() + RElog.get_yticklabels():
        label.set_size(25)
    for line in RElog.get_xticklines() + RElog.get_yticklines():
        line.set_markersize(10)
        line.set_markeredgewidth(2)

    
    IMlog=fig.add_subplot(222)
    try:   
        pylab.loglog(f,abs(numpy.asarray(Zlong).imag),f,abs(numpy.asarray(Zsc_ana).imag), linewidth=4 )

        #pylab.axvline(x=f_skin, linewidth=4, color='k')
        ## legend
        pylab.legend(( r'FEM imag',r'Zsc_ana imag'), shadow = True, fontsize=25)
        ltext = pylab.gca().get_legend().get_texts()
        pylab.setp(ltext[0], fontsize = 25)
        pylab.xlabel(r'$f \;$ [Hz]', fontsize=25)
        pylab.ylabel(r'Im $\{Z_\perp \} [\Omega/\mathrm{m}]$', fontsize=25)
        #pylab.title(tit)
        #pylab.setp(ltext[1], fontsize = 20)
        for label in IMlog.get_xticklabels() + IMlog.get_yticklabels():
            label.set_size(20)   #Fontsize Labels
        for line in IMlog.get_xticklines() + IMlog.get_yticklines():
            line.set_markersize(10) #Lenth of ticks
            line.set_markeredgewidth(2) #Thickness of ticks
    except:
        print ("Zero loglogplot")
    
    RElin=fig.add_subplot(223)
    try:
        pylab.semilogx(f,numpy.asarray(Zlong).real,f,numpy.asarray(Zlongloss).real, linewidth=4 )
        #pylab.axhline(y=R_Ohm, color='k')
        #pylab.axvline(x=f_skin, linewidth=4, color='k')
        #pylab.axhline(y=R_Ohm, color='k')
        pylab.axvline(x=f_skin, linewidth=4, color='k')
        pylab.axhline(y=R_Ohm, linewidth=4, color='k')
        ## legend
        pylab.legend((r'Zlong real',r'Zloss real'), shadow = True, fontsize=25)

        ltext = pylab.gca().get_legend().get_texts()

        pylab.xlabel(r'$f \;$ [Hz]', fontsize=25)
        pylab.ylabel(r'Re $\{Z_\parallel \} [\Omega]$', fontsize=25)
        #pylab.grid()
        for label in RElin.get_xticklabels() + RElin.get_yticklabels():
            label.set_size(25)
        for line in RElin.get_xticklines() + RElin.get_yticklines():
            line.set_markersize(10)
            line.set_markeredgewidth(2)
    except:
        print ("hmmmm")
    
    IMlin=fig.add_subplot(224)
    try:   
        pylab.semilogx(f,numpy.asarray(Zlong).imag,f,numpy.asarray(Zlongloss).imag , linewidth=4 )
        #pylab.axhline(y=R_Ohm, color='k')
        #pylab.axvline(x=f_skin, linewidth=4, color='k')
        #pylab.axhline(y=R_Ohm, color='k')
        #pylab.axvline(x=f_skin, linewidth=4, color='k')
        ## legend
        pylab.legend(('FEM imag',r'Zsc_ana imag'), shadow = True, fontsize=25)
        ltext = pylab.gca().get_legend().get_texts()
        pylab.setp(ltext[0], fontsize = 25)
        pylab.xlabel(r'$f \;$ [Hz]', fontsize=25)
        pylab.ylabel(r'Im $\{Z_\parallel \} [\Omega]$', fontsize=25)
        #pylab.title(tit)
        #pylab.setp(ltext[1], fontsize = 20)
        for label in IMlin.get_xticklabels() + IMlin.get_yticklabels():
            label.set_size(20)   #Fontsize Labels
        for line in IMlin.get_xticklines() + IMlin.get_yticklines():
            line.set_markersize(10) #Lenth of ticks
            line.set_markeredgewidth(2) #Thickness of ticks
    except:
        print ("hmmmm")
    
    """
    pylab.subplot(221)    
    try:
        pylab.loglog(f,numpy.fabs(numpy.asarray(Zlong).real),f,numpy.fabs(numpy.asarray(Zlongloss).real))
        pylab.axhline(y=R_Ohm, color='k')
        pylab.axvline(x=f_skin, linewidth=4, color='k')
        ## legend
        pylab.legend((r'Zlong real',r'Zloss real'), shadow = True)
        ltext = pylab.gca().get_legend().get_texts()
        pylab.setp(ltext[0], fontsize = 20)
        pylab.setp(ltext[1], fontsize = 20)
    except:
        print ("Zero loglogplot")
        #raise
    pylab.subplot(222) 
    try:   
        pylab.loglog(f,abs(numpy.asarray(Zlong).imag),f,abs(numpy.asarray(Zsc_ana).imag))
        pylab.axhline(y=R_Ohm, color='k')
        pylab.axvline(x=f_skin, linewidth=4, color='k')
        ## legend
        pylab.legend(( r'Zlong imag',r'Zsc_ana imag'), shadow = True)
        ltext = pylab.gca().get_legend().get_texts()
        pylab.setp(ltext[0], fontsize = 20)
        pylab.setp(ltext[1], fontsize = 20)
    except:
        print ("Zero loglogplot")
    pylab.subplot(223) 
    pylab.semilogx(f,numpy.asarray(Zlong).real,f,numpy.asarray(Zlongloss).real)
    pylab.axhline(y=R_Ohm, color='k')
    pylab.axvline(x=f_skin, linewidth=4, color='k')
    ## legend
    pylab.legend((r'Zlong real',r'Zloss real'), shadow = True)
    ltext = pylab.gca().get_legend().get_texts()
    pylab.setp(ltext[0], fontsize = 20)
    pylab.setp(ltext[1], fontsize = 20)


    pylab.subplot(224) 
    pylab.semilogx(f,numpy.asarray(Zlong).imag,f,numpy.asarray(Zsc_ana).imag)
    pylab.axhline(y=R_Ohm, color='k')
    pylab.axvline(x=f_skin, linewidth=4, color='k')
    ## legend
    pylab.legend(('Zlong imag',r'Zsc_ana imag'), shadow = True)
    ltext = pylab.gca().get_legend().get_texts()
    pylab.setp(ltext[0], fontsize = 20)
    pylab.setp(ltext[1], fontsize = 20)
    """
    
    #pylab.show()
    return None


def PlotZtrans(f,Ztrans,ZscTr_ana,tit):
    fig=pylab.figure(figsize=(20,12))
    pylab.suptitle(tit, fontsize=30)
    RE=fig.add_subplot(121)
    #pylab.subplot(121)    
    try:
        RE=pylab.loglog(f,numpy.fabs(numpy.asarray(Ztrans).real),color='r' , linewidth=4 )
        #pylab.axhline(y=R_Ohm, color='k')
        pylab.axvline(x=f_skin, linewidth=4, color='k')
        ## legend
        pylab.legend((r'FEM real',r'$f_{skin}$'), shadow = True, fontsize=25)

        ltext = pylab.gca().get_legend().get_texts()

        pylab.xlabel(r'$f \;$ [Hz]', fontsize=25)
        pylab.ylabel(r'Re $\{Z_\perp \} [\Omega/\mathrm{m}]$', fontsize=25)
        #pylab.grid()
        for label in RE.get_xticklabels() + RE.get_yticklabels():
            label.set_size(25)
        for line in RE.get_xticklines() + RE.get_yticklines():
            line.set_markersize(10)
            line.set_markeredgewidth(2)
    except:
        print ("Zero loglogplot")
        #raise



    IM=fig.add_subplot(122)
    #pylab.subplot(122) 
    try:   
        pylab.loglog(f,numpy.fabs(numpy.asarray(Ztrans).imag),f,numpy.fabs(numpy.asarray(ZscTr_ana).imag), linewidth=4)
        #pylab.axhline(y=R_Ohm, color='k')
        #pylab.axvline(x=f_skin, linewidth=4, color='k')
        ## legend
        pylab.legend(( r'FEM imag','Zsc_ana imag' ), shadow = True)
        ltext = pylab.gca().get_legend().get_texts()
        pylab.setp(ltext[0], fontsize = 25)
        pylab.xlabel(r'$f \;$ [Hz]', fontsize=25)
        pylab.ylabel(r'Im $\{Z_\perp \} [\Omega/\mathrm{m}]$', fontsize=25)
        #pylab.title(tit)
        #pylab.setp(ltext[1], fontsize = 20)
        for label in IM.get_xticklabels() + IM.get_yticklabels():
            label.set_size(20)   #Fontsize Labels
        for line in IM.get_xticklines() + IM.get_yticklines():
            line.set_markersize(10) #Lenth of ticks
            line.set_markeredgewidth(2) #Thickness of ticks
    except:
        print ("Zero loglogplot")
    #######################
    #pylab.show()
    return None

def PlotZtranslinear(f,Ztrans,ZscTr_ana,tit):
    fig=pylab.figure(figsize=(20,12))
    pylab.suptitle(tit, fontsize=30)
    RE=fig.add_subplot(121)
    #pylab.subplot(121)    
    try:
        RE=pylab.semilogx(f,numpy.fabs(numpy.asarray(Ztrans).real),color='r' , linewidth=4 )
        #pylab.axhline(y=R_Ohm, color='k')
        pylab.axvline(x=f_skin, linewidth=4, color='k')
        ## legend
        pylab.legend((r'FEM real',r'$f_{skin}$'), shadow = True, fontsize=25)

        ltext = pylab.gca().get_legend().get_texts()

        pylab.xlabel(r'$f \;$ [Hz]', fontsize=25)
        pylab.ylabel(r'Re $\{Z_\perp \} [\Omega/\mathrm{m}]$', fontsize=25)
        #pylab.grid()
        for label in RE.get_xticklabels() + RE.get_yticklabels():
            label.set_size(25)
        for line in RE.get_xticklines() + RE.get_yticklines():
            line.set_markersize(10)
            line.set_markeredgewidth(2)
    except:
        print ("hmmmm")
        #raise



    IM=fig.add_subplot(122)
    #pylab.subplot(122) 
    try:   
        pylab.semilogx(f,numpy.fabs(numpy.asarray(Ztrans).imag),f,numpy.fabs(numpy.asarray(ZscTr_ana).imag), linewidth=4)
        #pylab.axhline(y=R_Ohm, color='k')
        #pylab.axvline(x=f_skin, linewidth=4, color='k')
        ## legend
        pylab.legend(( r'FEM imag','Zsc_ana imag' ), shadow = True)
        ltext = pylab.gca().get_legend().get_texts()
        pylab.setp(ltext[0], fontsize = 25)
        pylab.xlabel(r'$f \;$ [Hz]', fontsize=25)
        pylab.ylabel(r'Im $\{Z_\perp \} [\Omega/\mathrm{m}]$', fontsize=25)
        #pylab.title(tit)
        #pylab.setp(ltext[1], fontsize = 20)
        for label in IM.get_xticklabels() + IM.get_yticklabels():
            label.set_size(20)   #Fontsize Labels
        for line in IM.get_xticklines() + IM.get_yticklines():
            line.set_markersize(10) #Lenth of ticks
            line.set_markeredgewidth(2) #Thickness of ticks
    except:
        print ("hmmmm")
    #######################
    #pylab.show()
    return None




def PlotWallCurrent(f,TotalCurrent):
    pylab.figure(figsize=(20,12))
    pylab.title('Wall current')
    pylab.semilogx(f,numpy.asarray(TotalCurrent).real,f,numpy.asarray(TotalCurrent).imag)
    pylab.legend((r'Current1 real',r'Current1 imag'), shadow = True)
    ltext = pylab.gca().get_legend().get_texts()
    pylab.setp(ltext[0], fontsize = 20)
    pylab.setp(ltext[1], fontsize = 20)
    #pylab.show()
    return None


def RealToComplex(Re,Im):
    Cplx=range(len(Re))
    for idx in range(len(Re)):
        Cplx[idx]=Re[idx]+I*Im[idx]
    
    return Cplx

"""
#################################################################################################################
#alt



def pppPlotZtranslinear(f,Ztrans,ZscTr_ana,tit):
    fig0=pylab.figure(figsize=(20,12))
    pylab.suptitle(tit)
    pylab.subplot(121)    
    try:
        pylab.semilogx(f,numpy.asarray(Ztrans).real )
        #pylab.axhline(y=R_Ohm, color='k')
        pylab.axvline(x=f_skin, linewidth=4, color='k')
        ## legend
        pylab.legend((r'Ztrans real'), shadow = True)
        ltext = pylab.gca().get_legend().get_texts()
        pylab.setp(ltext[0], fontsize = 20)
        #pylab.setp(ltext[1], fontsize = 20)
    except:
        print ("Zero loglogplot")
        #raise


    pylab.subplot(122) 
    try:   
        p=pylab.semilogx(f,numpy.asarray(Ztrans).imag,f,numpy.asarray(ZscTr_ana).imag)
        #pylab.title(tit)
        #pylab.axhline(y=R_Ohm, color='k')
        #pylab.axvline(x=f_skin, linewidth=4, color='k')
        ## legend
        pylab.legend(( r'Ztrans imag','Zsc_ana imag' ), shadow = True)
        ltext = pylab.gca().get_legend().get_texts()
        pylab.setp(ltext[0], fontsize = 20)
        #pylab.setp(ltext[1], fontsize = 20)
    except:
        print ("Zero loglogplot")
    #######################
    #pylab.show()
    return None
def n_div(ftr_,fti_,flr_,fli_,omega,beta,mesh,tspace):
    Vdivt = tspace
    Vdivl = FunctionSpace(mesh, "CG", div_long_order) 
    ftr=project(ftr_,Vdivt)
    fti=project(fti_,Vdivt)
    flr=project(flr_,Vdivl)
    fli=project(fli_,Vdivl)
    
    nL2tr=norm(ftr,'L2')
    nL2lr=norm(flr,'L2')
    nL2ti=norm(fti,'L2')
    nL2li=norm(fli,'L2')
    
    ndivr=norm(ftr,'Hdiv0')
    ndivi=norm(fti,'Hdiv0')
    

    nHdivr=(ndivr**2+(omega/(beta*c0)*nL2li)**2)**0.5    #Note that im and re are exchanged in the long. comp.!
    nL2r=(nL2tr**2+nL2lr**2)**0.5 
    
    nHdivi=(ndivi**2+(omega/(beta*c0)*nL2lr)**2)**0.5
    nL2i=(nL2ti**2+nL2li**2)**0.5
    
    nr=nHdivr/nL2r
    ni=nHdivi/nL2i
    
    ncomplex=((nHdivr**2+nHdivi**2)/(nL2r**2+nL2i**2))
    
    print
    print "ncomplex: ", ncomplex
    print
    if (nr>0.01 or ni>0.01):
        print "||Real_perp||_Hdiv0 : ", ndivr
        print "||Imag_perp||_Hdiv0 : ", ndivi 
        print "||Real||_Hdiv0 : ", nHdivr
        print "||Imag||_Hdiv0 : ", nHdivi 
        print "||Real||_L2: ", nL2r
        print "||Imag||_L2: ", nL2i
        print "nr: ",nr
        print "ni: ",ni
        print
    else:
        print "Div Norm is nice and small!"
    
    return [nr,ni]
    
"""
    
