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
    subdomains = MeshFunction("uint", mesh, file+"_physical_region.xml")
    print("Mesh imported")
    return [mesh,subdomains]
###################################################################################


#################################################
# Materials are converted into constants on the mesh triangles, i.e. projected on the zeroth order Discontinuous Galerkin basis
def MaterialOnMesh(mesh,subdomains,BeamMaterial,VacuumMaterial,Steel,Copper,Ferrite,Grounding,Dielectric):
    print("Imprintintg Material properties on mesh, i.e. generating DG-0 functions")
    

    #print 'nu= ', nuivalue,' +i ', nurvalue
    
    #Material parameters array over domains
    #epsilon_values=[eps0,eps0,epsvalue*eps0,eps0] #epsr=1+chi
    #nur_values=[nu0,nu0,nurvalue*nu0,nu0] 
    #nui_values=[0.0,0.0,nuivalue*nu0,0.0]
    #kappa_values=[0.0,0.0,0.0,0.0] #[1e6,1e6,1e6,1e6]
    
    #epsilon_values=[eps0,eps0,eps0,eps0] #epsr=1+chi
    #nur_values=[nu0,nu0,nu0,nu0] 
    #nui_values=[0.0,0.0,0.0,0.0]
    #kappa_values=[0.0,0.0,1e6,0.0] #[1e6,1e6,1e6,1e6]
    
    
    
    V0=FunctionSpace(mesh,'DG',0)   # Space for material distribution functions
    epsilonF=Function(V0)
    kappaF=Function(V0)
    nurF=Function(V0)
    nuiF=Function(V0)
        
        
    cell_markers0 = CellFunction("bool", mesh)    #beam
    cell_markers1 = CellFunction("bool", mesh)    #vacuum
    cell_markers2 = CellFunction("bool", mesh)    #Steel
    cell_markers3 = CellFunction("bool", mesh)    #Copper
    cell_markers4 = CellFunction("bool", mesh)    #Ferrite
    cell_markers5 = CellFunction("bool", mesh)    #Grounding
    cell_markers6 = CellFunction("bool", mesh)    #Dielectric
    cell_markers0.set_all(False)
    cell_markers1.set_all(False)
    cell_markers2.set_all(False)
    cell_markers3.set_all(False)
    cell_markers4.set_all(False)
    cell_markers5.set_all(False)
    cell_markers6.set_all(False)
    
    for cell in range(len(subdomains.array())):
        subdomain_no = subdomains.array()[cell]
        if subdomain_no==1:
            cell_markers0[cell]=True
        if subdomain_no==2:
            cell_markers1[cell]=True
        if subdomain_no==3:
            cell_markers2[cell]=True
        if subdomain_no==4:
            cell_markers3[cell]=True
        if subdomain_no==5:
            cell_markers4[cell]=True
        if subdomain_no==6:
            cell_markers5[cell]=True
        if subdomain_no==7:
            cell_markers6[cell]=True

            
                
        
    for cell_number in range(len(subdomains.array())):
        subdomain_no = subdomains.array()[cell_number]
        #print('subdomain_no cell_no',subdomain_no, cell_number)
        if subdomain_no==1:
            epsilonF.vector()[cell_number]=BeamMaterial.Eps
            kappaF.vector()[cell_number]=BeamMaterial.Kappa
            nurF.vector()[cell_number]=BeamMaterial.Nur
            nuiF.vector()[cell_number]=BeamMaterial.Nui
        if subdomain_no==2:
            epsilonF.vector()[cell_number]=VacuumMaterial.Eps
            kappaF.vector()[cell_number]=VacuumMaterial.Kappa
            nurF.vector()[cell_number]=VacuumMaterial.Nur
            nuiF.vector()[cell_number]=VacuumMaterial.Nui
        if subdomain_no==3:
            epsilonF.vector()[cell_number]=Steel.Eps
            kappaF.vector()[cell_number]=Steel.Kappa
            nurF.vector()[cell_number]=Steel.Nur
            nuiF.vector()[cell_number]=Steel.Nui
        if subdomain_no==4:
            epsilonF.vector()[cell_number]=Copper.Eps
            kappaF.vector()[cell_number]=Copper.Kappa
            nurF.vector()[cell_number]=Copper.Nur
            nuiF.vector()[cell_number]=Copper.Nui

        if subdomain_no==5:
            epsilonF.vector()[cell_number]=Ferrite.Eps
            kappaF.vector()[cell_number]=Ferrite.Kappa
            nurF.vector()[cell_number]=Ferrite.Nur
            nuiF.vector()[cell_number]=Ferrite.Nui
        if subdomain_no==6:
            epsilonF.vector()[cell_number]=Grounding.Eps
            kappaF.vector()[cell_number]=Grounding.Kappa
            nurF.vector()[cell_number]=Grounding.Nur
            nuiF.vector()[cell_number]=Grounding.Nui    
        if subdomain_no==7:
            epsilonF.vector()[cell_number]=Dielectric.Eps
            kappaF.vector()[cell_number]=Dielectric.Kappa
            nurF.vector()[cell_number]=Dielectric.Nur
            nuiF.vector()[cell_number]=Dielectric.Nui
    
    print("Material properties as DG-0 functions on mesh available")
    return [nurF,nuiF,epsilonF,kappaF]
### The distribution of nur, nui, kappa and eps is returned
######################################################################

