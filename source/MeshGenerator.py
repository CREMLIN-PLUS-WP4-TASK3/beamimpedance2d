"""
FEniCS 2D beam coupling impedance simulation in frequency domain

by Uwe Niedermayer 2014
"""

import os
from dolfin import *
from dolfin_utils import meshconvert

from source import Material


###################################################################################
# This runs dolfin-convert before every mesh import
def MeshImport(file):
    print("Importing mesh")
    #Importing mesh from GMSH
    os.system('dolfin-convert ' + file + '.msh '+ file + '.xml')
    mesh=Mesh(file+".xml")
    subdomains = MeshFunction("size_t", mesh, file+"_physical_region.xml")
    print("Mesh imported")
    return [mesh,subdomains]
###################################################################################


#################################################
# Materials are converted into constants on the mesh triangles, i.e. projected on the zeroth order Discontinuous Galerkin basis
def MaterialOnMesh(mesh,subdomains,BeamMaterial,VacuumMaterial,Steel,Copper,Ferrite,Titanium,Dielectric):
    print("Imprintintg Material properties on mesh, i.e. generating DG-0 functions")


    #print ('nu= ', nuivalue,' +i ', nurvalue)

    #Material parameters array over domains
    #epsilon_values=[eps0,eps0,epsvalue*eps0,eps0] #epsr=1+chi
    #nur_values=[nu0,nu0,nurvalue*nu0,nu0]
    #nui_values=[0.0,0.0,nuivalue*nu0,0.0]
    #kappa_values=[0.0,0.0,0.0,0.0] #[1e6,1e6,1e6,1e6]

    #epsilon_values=[eps0,eps0,eps0,eps0] #epsr=1+chi
    #nur_values=[nu0,nu0,nu0,nu0]
    #nui_values=[0.0,0.0,0.0,0.0]
    #kappa_values=[0.0,0.0,1e6,0.0] #[1e6,1e6,1e6,1e6]



    V0=FunctionSpace(mesh,'Discontinuous Lagrange',0)   # Space for material distribution functions
    epsilonF=Function(V0)
    kappaF=Function(V0)
    nurF=Function(V0)
    nuiF=Function(V0)

    epsilonF.vector()[subdomains.where_equal(1)]=BeamMaterial.Eps
    kappaF.vector()[subdomains.where_equal(1)]=BeamMaterial.Kappa
    nurF.vector()[subdomains.where_equal(1)]=BeamMaterial.Nur
    nuiF.vector()[subdomains.where_equal(1)]=BeamMaterial.Nui

    epsilonF.vector()[subdomains.where_equal(2)]=VacuumMaterial.Eps
    kappaF.vector()[subdomains.where_equal(2)]=VacuumMaterial.Kappa
    nurF.vector()[subdomains.where_equal(2)]=VacuumMaterial.Nur
    nuiF.vector()[subdomains.where_equal(2)]=VacuumMaterial.Nui

    epsilonF.vector()[subdomains.where_equal(3)]=Steel.Eps
    kappaF.vector()[subdomains.where_equal(3)]=Steel.Kappa
    nurF.vector()[subdomains.where_equal(3)]=Steel.Nur
    nuiF.vector()[subdomains.where_equal(3)]=Steel.Nui

    epsilonF.vector()[subdomains.where_equal(4)]=Copper.Eps
    kappaF.vector()[subdomains.where_equal(4)]=Copper.Kappa
    nurF.vector()[subdomains.where_equal(4)]=Copper.Nur
    nuiF.vector()[subdomains.where_equal(4)]=Copper.Nui

    epsilonF.vector()[subdomains.where_equal(5)]=Ferrite.Eps
    kappaF.vector()[subdomains.where_equal(5)]=Ferrite.Kappa
    nurF.vector()[subdomains.where_equal(5)]=Ferrite.Nur
    nuiF.vector()[subdomains.where_equal(5)]=Ferrite.Nui

    epsilonF.vector()[subdomains.where_equal(6)]=Titanium.Eps
    kappaF.vector()[subdomains.where_equal(6)]=Titanium.Kappa
    nurF.vector()[subdomains.where_equal(6)]=Titanium.Nur
    nuiF.vector()[subdomains.where_equal(6)]=Titanium.Nui

    epsilonF.vector()[subdomains.where_equal(7)]=Dielectric.Eps
    kappaF.vector()[subdomains.where_equal(7)]=Dielectric.Kappa
    nurF.vector()[subdomains.where_equal(7)]=Dielectric.Nur
    nuiF.vector()[subdomains.where_equal(7)]=Dielectric.Nui

    print("Material properties as DG-0 (Discontinuous Galerkin) functions on mesh available")
    return [nurF,nuiF,epsilonF,kappaF]
### The distribution of nur, nui, kappa and eps is returned
######################################################################
