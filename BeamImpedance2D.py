#!/usr/bin/env python3

"""
FEniCS 2D beam coupling impedance simulation in frequency domain

by Uwe Niedermayer 2014

niedermayer@temf.tu-darmstadt.de
u.niedermayer@gsi.de

...developed on Python 2.7.6
"""

import os
from dolfin import *
#from dolfin_utils import meshconvert

import matplotlib
import pylab
import numpy
import scipy
import scipy.special

from source import Globals
from source import CurlSolver
from source import MeshGenerator
from source import PostProc
from source import PoissonSolver
from source import RHS
#from source import Boundary
from source import Material


# constants: (SI units)
zero=Constant(0.0)



#############################################################
#Import mesh
#MeshFileName='SIS100-TransferKickerVacGap'
#MeshFileName='SIS100-TransferKicker'
#MeshFileName='SIS100-EmergencyKickerVacGap'
# MeshFileName='SIS100-EmergencyKicker'
#MeshFileName='EllipticPipeDeltaZt'
#MeshFileName='simplepipeDeltaZt_S2'
#MeshFileName='collimator_shifted'
#####MeshFileName='simplepipeDeltaZt'
# MeshFileName='simplepipe'
MeshFileName='ThinShellPipeZt_eps'
#MeshFileName='fccCoating300'
#MeshFileName='fcc'
#MeshFileName='fccWithoutHole'
#MeshFileName='fccSIBC'
#MeshFileName='fccWithoutHoleSIBC'
#MeshFileName='ThinShellPipeGNDfine'
##MeshFileName='simplepipeZlspch'
#MeshFileName='extrudepipe'
#MeshFileName='extrudepipe3fineoutside'
#MeshFileName='FerriteRingZt'
[mesh,subdomains]= MeshGenerator.MeshImport(MeshFileName)
print ("Mesh imported and converted! File: " + MeshFileName )
#############################################################

output_file = File('mesh.pvd')
output_file << mesh
output_file << subdomains

######################################################################################
#Beam parameters (SI units)
px=0.02  #x-position
py=0.00  #y-position
a=0.01  #beam radius
s_mesh=1e-6   #Parameter for delta function representation
eps=s_mesh*0.1   #Apply Delta-function representation on little smaller epsilon (safety-factor)
#length=0.0254
length=1.0

gamma_in=50000.0
#gamma_in=100.0
#beta=(1.0-1.0/gamma_in**2)**0.5

#beta=0.947485301 #2GeV protons
beta=0.9999999 #roughly 7TeV
# beta=0.9
# beta=0.1

gamma=1/sqrt(1-beta**2)
bgsinv= 1.0/(beta**2 *gamma**2) # 1/(beta^2gamma^2)

#beta=1.0

print ("beta ",beta)
print ("gamma ",gamma)

# Charge should be normalized / does not matter
q=1.0  #Charge (As)
#######################################################################################


###########################################################################################################
# Variables: (SI units)
# Set these for analytic reference in the plots
b=0.04 #pipe inner radius
h=0.0403 #pipe outer radius
h2=0.1 #boundary radius
#length=1.0
g_ana=0.25+ln(b/a) #longitudinal space charge impedance geometry factor
#print('g_ana: ',g_ana)
###########################################################################################################


###############################################################
#Frequency stepping
maxfpoints=30 #Number of frequency points to compute
#logscale
startfexp=3.0
stopfexp=7.0
#linear scale
startf=5e9
stopf=5e9
###############################################################


#############################################################
#Import Material Data for specific Ferrite
if dispersive:
    #Material.MaterialTest()
    [fArray,murArray,muiArray]=Material.PermeabilityRead()
    #PostProc.PlotZtrans(fArray,PostProc.RealToComplex(murArray,muiArray),PostProc.RealToComplex(murArray,muiArray))
    #pylab.show()
#############################################################


##############################################################
#parameters["reorder_dofs_serial"] = False #Order dofs in the same way as vertices
##############################################################

######################################
if plot3Dflag:
    viz=plot(mesh, basename='mesh')
    #viz.write_png('mesh')
######################################

######################################################################################################################
# Function Spaces and Functions
#Vdivtr = VectorFunctionSpace(mesh, "CG", div_order)    #Space for transverse static fields
Hcurl= FunctionSpace(mesh, "Nedelec 1st kind H(curl)", curl_order)
H1= FunctionSpace(mesh, "Lagrange", curl_long_order)

Edivtr=Function(Hcurl)  #The solution of the Poisson solver (Edivt=-grad Phi)
Edivti=Function(Hcurl)

RHSvr=Function(Hcurl)   #The RHS for the curlcurl solver
RHSvi=Function(Hcurl)
RHSsr=Function(H1)
RHSsi=Function(H1)

Etr=Function(Hcurl)     #The final solution
Eti=Function(Hcurl)
Elr=Function(H1)
Eli=Function(H1)

divTesttr=Function(Hcurl)
divTestti=Function(Hcurl)

Jcondr=Function(H1)     #The wall current
Jcondi=Function(H1)
#####################################################################################################################



##############################################################
###Sources
if dipole:
    if quadrupole:
        Jszr=project(RHS.ExCurrentShiftedQuadrupole(mesh,subdomains,q,a,eps,px,py),H1)  ####!!!
    else:
        Jszr=project(RHS.ExCurrentShiftedDipole(mesh,subdomains,q,a,eps,px,py),H1)  ####!!!

    #Jszr=project(RHS.ExCurrentShifted(mesh,subdomains,q,a,0.0,d/2.0)-RHS.ExCurrentShifted(mesh,subdomains,q,a,0.0,-d/2.0),Vcurllr)
    Monointegral=assemble(Jszr*dx)
    if horizontal:
        DIP=Expression('x[0]-px',px=px, degree=2)
    else:
        DIP=Expression('x[1]-py',py=py, degree=2)

    Dipintegral=assemble((Jszr*DIP)*dx)

    Quad=Expression('(x[0]-px)*(x[0]-px)-(x[1]-py)*(x[1]-py)',px=px, py=py, degree=2)
    Quadintegral=assemble((Jszr*Quad)*dx)
    #mp=project(RHS.ExCurrentShifted(mesh,subdomains,q,a,0.0,0.0),Vcurllr)
    #Monointegral=assemble(mp*dx)
    print ("Quadintegral: ", Quadintegral)
    print ("Dipintegral: ", Dipintegral)
    print ("Monointegral: ", Monointegral)
    #print ("ratio:" , Dipintegral/Monointegral)
else:
    #monopole
    Jszr=project(RHS.ExCurrentShifted(mesh,subdomains,q,a,0.0,0.0),H1)             ####!!!
    Monointegral=assemble(Jszr*dx)
    print ("Monopole source integral: ", Monointegral)
    #Jszr=RHS.ExCurrent(mesh,subdomains,q,a)

Jszi=zero

output_file = File('Jszr.pvd')
output_file << Jszr

if plot3Dflag:
        viz=plot(Jszr, mesh=mesh,interactive=True,title='Jszr', basename='Jszr')
        #viz.write_png('Jszr.png')
        #interactive()

# tcur=assemble(Jszr*dx, mesh=mesh)
# print ('tcur= ', tcur)

rhor=Jszr*1.0/(beta*c0)
rhoi=Jszi*1.0/(beta*c0)
# plot(Jszr,mesh=mesh,interactive=True)
##############################################################


##############################################################################################
# Create frequency arrays
# Note that omega is just a scalar for the current frequency point
if logscale:
    f = numpy.logspace(startfexp, stopfexp, num=maxfpoints)
else:
    f = numpy.linspace(startf, stopf, num=maxfpoints)


if dipole:
    pass
    ZscTr_ana=numpy.zeros(maxfpoints, dtype="complex")
    ZscTr_ana_direct=numpy.zeros(maxfpoints, dtype="complex")
    ZscTr_ana_indirect=numpy.zeros(maxfpoints, dtype="complex")
    Ztrans_ind=numpy.zeros(maxfpoints, dtype="complex")
    Ztrans_full=numpy.zeros(maxfpoints, dtype="complex")
else:
    Zsc_ana=numpy.zeros(maxfpoints, dtype="complex")
    Zlong=numpy.zeros(maxfpoints, dtype="complex")
    # Zlongloss=numpy.zeros(maxfpoints, dtype="complex")
    # TotalCurrent=numpy.zeros(maxfpoints, dtype="complex")

################################################################################################
#initialize material
if dispersive:
    [nurFvalue,nuiFvalue]=Material.ReluctivityInterpolate(1.0e5,fArray,murArray,muiArray)
    epsFvalue=10
    Ferrite=Material.MaterialProperties(nurFvalue*nu0,nuiFvalue*nu0,epsFvalue*eps0,0.0)
else:
    Ferrite=Material.MaterialProperties(nu0,0.0,eps0,0.0)  #a dummy if there is no ferrite

BeamMaterial=Material.MaterialProperties(1.0*nu0,0.0,1.0*eps0,0.0)
VacuumMaterial=Material.MaterialProperties(1.0*nu0,0.0,1.0*eps0,0.0)
Steel=Material.MaterialProperties(1.0*nu0,0.0,1.0*eps0,1.4e8)
Copper=Material.MaterialProperties(1.0*nu0,0.0,1.0*eps0,5.8e7) #Equivalent kappa
Titanium=Material.MaterialProperties(1.0*nu0,0.0,1.0*eps0,1.8e8) #1e-6
Dielectric=Material.MaterialProperties(1.0*nu0,0.0,100.0*eps0,0.0)


[nur,nui,epsilon,kappa]= MeshGenerator.MaterialOnMesh(mesh,subdomains,BeamMaterial,VacuumMaterial,Steel,Copper,Ferrite,Titanium,Dielectric)

# plot(nui,mesh=mesh)
#interactive()
print ('Material initialized!')
#################################################################################################

output_file = File('material_props.pvd')
output_file << mesh
output_file << nur
output_file << nui
output_file << epsilon
output_file << kappa


#The big frequency loop
for i, fpoint in enumerate(f):

    omega=2*pi*fpoint
    print ("fpoint: ", fpoint,"  f= ",fpoint/1e6,"MHz" )

    ###############################################################################
    #Analytical references
    if dipole:
        if quadrupole:
            ZscTr_ana_direct[i]=-I*bgsinv*(beta*c0*mu0)*length/(pi*a**2)* \
                scipy.special.iv(2,(omega*a)/(beta*gamma*c0))*scipy.special.kn(2,(omega*a)/(beta*gamma*c0))
            ZscTr_ana_indirect[i]=I*bgsinv*(beta*c0*mu0)*length/(pi*a**2)* \
                (scipy.special.iv(2,(omega*a)/(beta*gamma*c0)))**2 \
                *scipy.special.kn(2,(omega*b)/(beta*gamma*c0))/scipy.special.iv(2,(omega*b)/(beta*gamma*c0))
            ZscTr_ana[i]=ZscTr_ana_direct[i]+ZscTr_ana_indirect[i]
        else:
            ZscTr_ana_direct[i]=-I*bgsinv*(beta*c0*mu0)*length/(pi*a**2)* \
                scipy.special.iv(1,(omega*a)/(beta*gamma*c0))*scipy.special.kn(1,(omega*a)/(beta*gamma*c0))  #check me
            ZscTr_ana_indirect[i]=I*bgsinv*(beta*c0*mu0)*length/(pi*a**2)* \
                (scipy.special.iv(1,(omega*a)/(beta*gamma*c0)))**2 \
                *scipy.special.kn(1,(omega*b)/(beta*gamma*c0))/scipy.special.iv(1,(omega*b)/(beta*gamma*c0))
            ZscTr_ana[i]=ZscTr_ana_direct[i]+ZscTr_ana_indirect[i]
    else:
        Zsc_ana[i]=-I*omega*mu0*bgsinv*g_ana/(2*pi)*length
        #ZscTr_ana[i]=-I*bgsinv*(beta*c0*mu0)*(1.0/a**2-1.0/b**2) *length/(2*pi)  #MQS
    ################################################################################

    ###################################################################################################
    #Determine frequency dependent material parameters, update, and then imprint again on the mesh
    if (dispersive):
        [nurFvalue,nuiFvalue]=Material.ReluctivityInterpolate(f[i],fArray,murArray,muiArray)
        epsFvalue=10
        print ("nur:", nurFvalue)
        print ("nui:", nuiFvalue)
        Ferrite.reluctivityUpdate(nurFvalue*nu0,nuiFvalue*nu0)
        [nur,nui,epsilon,kappa]= MeshGenerator.MaterialOnMesh(mesh,subdomains,BeamMaterial,VacuumMaterial,Steel,Copper,Ferrite,Titanium,Dielectric)
        #plot(nui,mesh=mesh,interactive=True)
        #print (nui.vector().array())
        print ('New material parameters imprinted on mesh')
    ####################################################################################################



    ####################################################################################################################
    ##Calculate Ediv
    [Phir,Phii]=PoissonSolver.CplxPoisson(mesh,omega, beta, epsilon, kappa, Jszr,Jszi)
    """
    $E_{div\;tr}^{\Re}=-\nabla\Phi^{\Re}$
    $E_{div\;tr}^{\Im}=-\nabla\Phi^{\Im}$
    $\left(\partial_z \rightarrow -i \omega/v\right)$
    $E_{div\;lon}^{\Re}=-\frac{\omega}{\beta c}\Phi^{\Im}$
    $E_{div\;lon}^{\Im}=\frac{\omega}{\beta c}\Phi^{\Re}$
    """
    Edivtr=project(-nabla_grad(Phir),Hcurl)
    Edivti=project(-nabla_grad(Phii),Hcurl)
    Edivlr= project(-omega/(beta *c0) * Phii,H1)
    Edivli= project(omega/(beta *c0) * Phir,H1)

    output_file = File('Phi.pvd')
    output_file << Phir
    output_file << Phii

    output_file = File('Ediv.pvd')
    output_file << Edivtr
    output_file << Edivti
    output_file << Edivlr
    output_file << Edivli


    # Jdivtr=project(-kappa*nabla_grad(Phir),Hcurl)
    # Jdivti=project(-kappa*nabla_grad(Phii),Hcurl)
    # Jdivlr= project(-omega/(beta *c0) *kappa*Phii,H1)
    # Jdivli= project(omega/(beta *c0) *kappa*Phir,H1)

    # output_file = File('Jdiv.pvd')
    # output_file << Jdivtr
    # output_file << Jdivti
    # output_file << Jdivlr
    # output_file << Jdivli

    ####################################################################################################################
    ####################################
    #Plot Ediv
    if(plot3Dflag):
        plot(Edivtr,title='Edivtr',basename='Edivtr')
        plot(Edivti,title='Edivti',basename='Edivti')
        plot(Edivlr,title='Edivlr',basename='Edivlr')
        plot(Edivli,title='Edivli',basename='Edivli')
        interactive()
        #plot(Jdivlr,'Jdivlr')
        #plot(Jdivli,'Jdivli')
        #plot(Jdivtr,'Jdivtr')
        #plot(Jdivti,'Jdivti')
        #interactive()
    ####################################

    ####################################################################################################################
    #calculate RHS
    """
    $R_v^\Re=\varepsilon\omega^2 E_{div\;tr}^{\Re}+\omega\kappa E_{div\;tr}^{\Im}$
    $R_v^\Im=\varepsilon\omega^2 E_{div\;tr}^{\Im}-\omega\kappa E_{div\;tr}^{\Re}$
    $R_s^\Re=\varepsilon\omega^2 E_{div\;lon}^{\Re}+\omega\kappa E_{div\;lon}^{\Im}$
    $R_s^\Im=\varepsilon\omega^2 E_{div\;lon}^{\Im}+\omega\kappa E_{div\;lon}^{\Re}-\omega J_{sz}^{\Re}$
    """
    RHSvr=omega*omega*epsilon*Edivtr+omega*kappa*Edivti
    RHSvi=omega*omega*epsilon*Edivti-omega*kappa*Edivtr
    RHSsr=omega*omega*epsilon*Edivlr+omega*kappa*Edivli
    RHSsi=omega*omega*epsilon*Edivli+omega*kappa*Edivlr-omega*Jszr
    ####################################################################################################################
    ####################################
    #Plot RHS
    if(plot3Dflag):
        plot(RHSvr,title='RHSvr')
        plot(RHSvi,title='RHSvi')
        plot(RHSsr,title='RHSsr')
        plot(RHSsi,title='RHSsi')
        interactive()

    # output_file = File('CurlCurlRHS.pvd')
    # output_file << RHSvr
    # output_file << RHSvi
    # output_file << RHSsr
    # output_file << RHSsi

    ####################################

    #######################################################
    #Check div norm of RHS
    #print
    #print ("Div-Norm of RHS")
    #PostProc.n_div(RHSvr,RHSvi,HSsr,RHSsi,omega,beta,mesh,Vcurltr)
    #########################################################

    ######
    #Run curlcurl solver
    # [Ecurltr,Ecurlti,Ecurllr,Ecurlli]=CurlSolver.CurlCurl(mesh,omega, beta, epsilon, kappa, nur, RHSsr/omega**2, RHSsi/omega**2, RHSvr/omega**2,RHSvi/omega**2)
    [Ecurltr,Ecurlti,Ecurllr,Ecurlli]=CurlSolver.CurlCurlCplxNu(mesh,omega, beta, epsilon, kappa, nur,nui, RHSsr/omega**2, RHSsi/omega**2, RHSvr/omega**2,RHSvi/omega**2)   #sign RHSvr
    #Rescaling with omega brings LF stabilization???
    ######

    output_file = File('CurlCurl.pvd')
    output_file << Ecurllr
    output_file << Ecurlli
    output_file << Ecurltr
    output_file << Ecurlti


    # #######################################################
    # #Check div norm of Ecurl
    # #print
    # #print ("Div-Norm of Ecurl")
    # #PostProc.n_div(Ecurltr,Ecurlti,Ecurllr,Ecurlli,omega,beta,mesh,Vcurltr)
    # #########################################################


    #These are both in H(curl)
    Etr=omega**2*Ecurltr+Edivtr
    Eti=omega**2*Ecurlti+Edivti
    Elr=omega**2*Ecurllr+Edivlr
    Eli=omega**2*Ecurlli+Edivli
    # Etr=Ecurltr+Edivtr
    # Eti=Ecurlti+Edivti
    # Elr=Ecurllr+Edivlr
    # Eli=Ecurlli+Edivli
    ######
    """
    filer = File("Edivtr.pvd")
    filer << Edivtr
    filei = File("Edivti.pvd")
    filei << Edivti
    print ("huhu")
    """

    # Etr = Function(Hcurl)
    # Etr.vector()[:] = Ecurltr.vector() + Edivtr.vector()
    # Eti = Function(Hcurl)
    # Eti.vector()[:] = Ecurlti.vector() + Edivti.vector()
    # Elr = Function(H1)
    # Elr.vector()[:] = Ecurllr.vector() + Edivlr.vector()
    # Eli = Function(H1)
    # Eli.vector()[:] = Ecurlli.vector() + Edivli.vector()

    # output_file = File('E.pvd')
    # output_file << Etr
    # output_file << Eti
    # output_file << Elr
    # output_file << Eli



    ##################################################################################################
    #Calculate Impedance
    if dipole:
        Ztrans_full[i]=length*PostProc.Ztrans(Elr, Eli, Jszr, Jszi, omega, beta, px,py)
        if quadrupole:
            Ztrans_ind[i]=a**2*Ztrans_full[i]-ZscTr_ana_direct[i] #test!
        else:
            Ztrans_ind[i]=Ztrans_full[i]-ZscTr_ana_direct[i] #test!
        print
        print ("Ztrans_full= " , Ztrans_full[i], " at ", omega/(2*pi),"Hz")
        print ("Ztrans_indirect= " , Ztrans_ind[i], " at ", omega/(2*pi),"Hz")
        print
    else:
        Zlong[i]=length*PostProc.Zlong(Elr, Eli, Jszr, Jszi, q)
        print
        print ("Zlong= " , Zlong[i], " at ", omega/(2*pi),"Hz")
        print
    ###################################################################################################

    ############################################################
    #calculate conduction current
    if wallcurrent:
        Jcondr=Elr*kappa
        Jcondi=Eli*kappa
        TotalCurrent[i]=assemble(Jcondr*dx)+I*assemble(Jcondi*dx)
        print ('Totalcurrent: ', TotalCurrent[i])
        Zlongloss[i]=(length/q**2)*(assemble(inner(Elr,Jcondr)*dx)+assemble(inner(Eli,Jcondi)*dx))
    ##############################################################


    if(plot3Dflag):
        plot(Etr,title='Etr')
        plot(Eti,title='Eti')
        plot(Elr,title='Elr')
        plot(Eli,title='Eli')
        interactive()

#End of f-loop
###########################################################################

##################################################################################################################

#############################################################
# Export impedance
if dataexport:
    ExportName=MeshFileName+'beta'+str(beta)
    if dipole:
        PostProc.CplxImpExport(ExportName+'_a'+str(a)+'_s'+str(s_mesh)+'horiz_'+str(horizontal)+'Ztr_full.dat',f,Ztrans_full)
        PostProc.CplxImpExport(ExportName+'_a'+str(a)+'_s'+str(s_mesh)+'horiz_'+str(horizontal)+'Ztr_ind.dat',f,Ztrans_ind)
        PostProc.CplxImpExport(ExportName+'_a'+str(a)+'_s'+str(s_mesh)+'Ztr_ana_full.dat',f,ZscTr_ana)
        PostProc.CplxImpExport(ExportName+'_a'+str(a)+'_s'+str(s_mesh)+'Ztr_ana_ind.dat',f,ZscTr_ana_indirect)
    else:
        PostProc.CplxImpExport(ExportName+'Zl.dat',f,Zlong)
    if wallcurrent:
        PostProc.CplxImpExport(ExportName+'ZlLoss.dat',f,Zlongloss)
        PostProc.CplxImpExport(ExportName+'Current.dat',f,TotalCurrent)

#############################################################



#############################################################
#Plot impedance
if dipole:
    if quadrupole:
        PostProc.PlotZtrans(f,numpy.multiply(a**2,Ztrans_full),ZscTr_ana,'Full Transverse Impedance')
        PostProc.PlotZtrans(f,Ztrans_full,numpy.multiply(1/(2*a**2),ZscTr_ana),'Full Transverse Impedance')
        PostProc.PlotZtrans(f,Ztrans_ind,ZscTr_ana_indirect,'Indirect Transverse Impedance')
    else:
        PostProc.PlotZtranslinear(f,Ztrans_full,ZscTr_ana,'Full Transverse Impedance')
        PostProc.PlotZtrans(f,Ztrans_full,ZscTr_ana,'Full Transverse Impedance')
        PostProc.PlotZtrans(f,Ztrans_ind,ZscTr_ana_indirect,'Indirect Transverse Impedance')
else:
    if wallcurrent:
        PostProc.PlotZlong(f,Zlong,Zlongloss,Zsc_ana)
    else:
        PostProc.PlotZlong(f,Zlong,Zlong,Zsc_ana)
##############################################################


# ################################################################################
# #Plot Conduction current
# if wallcurrent:
#     PostProc.PlotWallCurrent(f,TotalCurrent)
# ################################################################################
pylab.savefig('output.png')
