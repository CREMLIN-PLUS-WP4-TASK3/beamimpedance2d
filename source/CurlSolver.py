"""
FEniCS 2D beam coupling impedance simulation in frequency domain

by Uwe Niedermayer 2014
"""


from dolfin import *
import numpy as np
###################################################################
#operators

def my_cross(a,b):
    return a[0]*b[1]-a[1]*b[0]


def my_Tangential(n):          ##### ez cross
    return as_vector((-n[1],n[0]))



def Zmat(vec,omega,beta):
    rx,ry=variable(vec)
    rx=(omega/(beta*c0))*vec[1]
    ry=-(omega/(beta*c0))*vec[0]
    return as_vector((rx,ry))


def Amat(scal):
    rx=variable(scal)
    ry=variable(scal)
    rx=scal.dx(1)
    ry=-scal.dx(0)
    return as_vector((rx,ry))

def Bmat(vec):
    rz=variable(vec)
    rz=-vec[0].dx(1)+vec[1].dx(0)
    return rz


def ZAmat(scal,omega,beta):
    rx=variable(scal)
    ry=variable(scal)
    rx=-(omega/(beta*c0))*scal.dx(0)
    ry=-(omega/(beta*c0))*scal.dx(1)
    return as_vector((rx,ry))




###################################################################
def CurlCurlCplxNu(mesh,omega, beta, epsilon, kappa, nur,nui, RHSsr, RHSsi, RHSvr,RHSvi):
    domains = MeshFunction("size_t", mesh, mesh.topology().dim())
    domains.set_all(0)
    dx = Measure("dx")(subdomain_data=domains)

    ################################################################
    # SIBC
    class AdmittanceBoundary(SubDomain):
        def inside(self,x,on_boundary):
            return on_boundary

    AdmittanceBoundaryObject=AdmittanceBoundary()

    Abdy = MeshFunction("size_t", mesh, mesh.topology().dim()-1)
    Abdy.set_all(0)

    AdmittanceBoundaryObject.mark(Abdy,0)
    ds = Measure("ds")(subdomain_data=Abdy)
    n = -FacetNormal(mesh)  #Normal pointing inwards
    #plot(assemble(n,mesh=mesh),mesh=mesh,interactive=True)
    #plot(n,mesh=mesh)
    ################################################################

    Vtr = VectorElement("Nedelec 1st kind H(curl)", mesh.ufl_cell(), curl_order)
    Vlr = FiniteElement("Lagrange", mesh.ufl_cell(), curl_long_order)
    Vli=Vlr
    Vti=Vtr

    V= FunctionSpace(mesh, MixedElement([Vtr, Vti, Vlr, Vli])) # 3D complex vector space

    (uvr,uvi,usr,usi) = TrialFunctions(V)
    (vvr,vvi,vsr,vsi) = TestFunctions(V)

    Ecurl = Function(V)
    #Ecurltr=Function(Vtr)
    #Ecurlti=Function(Vti)
    #Ecurllr=Function(Vlr)
    #Ecurlli=Function(Vli)


    a_epsilon=-omega*omega*(inner(vvr,epsilon*uvr)*dx(0)  + inner(vvi,epsilon*uvi)*dx(0)  + inner(vsr,epsilon*usr)*dx(0)  + inner(vsi,epsilon*usi)*dx(0) )
    a_kappa=omega*(inner(vvi,kappa*uvr)*dx(0)  - inner(vvr,kappa*uvi)*dx(0)  + inner(vsi,kappa*usr)*dx(0)  - inner(vsr,kappa*usi)*dx(0) )


    #Only this part changes wrt the real permeabiliy case!
    a_curlcurl=  inner(Bmat(vvr),nur*Bmat(uvr))*dx(0)  -inner(vvr,-(omega/(beta*c0))**2 *nur*uvr)*dx(0) \
                -inner(vvr,nui*ZAmat(usr,omega,beta))*dx(0)  \
                -(inner(Bmat(vvr),nui*Bmat(uvi))*dx(0)  -inner(vvr,-(omega/(beta*c0))**2 *nui*uvi)*dx(0) ) \
                -inner(vvr,nur*ZAmat(usi,omega,beta))*dx(0)  \
                -inner(Amat(vsr),nui*Zmat(uvr,omega,beta))*dx(0)  \
                +inner(Amat(vsr),nur*Amat(usr))*dx(0)  \
                -inner(Amat(vsr),nur*Zmat(uvi,omega,beta))*dx(0)  \
                -inner(Amat(vsr),nui*Amat(usi))*dx(0)  \
                \
                +inner(Bmat(vvi),nur*Bmat(uvi))*dx(0)  -inner(vvi,-(omega/(beta*c0))**2 *nur*uvi)*dx(0)  \
                -inner(vvi,nui*ZAmat(usi,omega,beta))*dx(0)  \
                +inner(Bmat(vvi),nui*Bmat(uvr))*dx(0)  -inner(vvi,-(omega/(beta*c0))**2 *nui*uvr)*dx(0) \
                +inner(vvi,nur*ZAmat(usr,omega,beta))*dx(0)   \
                -inner(Amat(vsi),nui*Zmat(uvi,omega,beta))*dx(0)  \
                +inner(Amat(vsi),nur*Amat(usi))*dx(0)  \
                +inner(Amat(vsi),nur*Zmat(uvr,omega,beta))*dx(0)  \
                +inner(Amat(vsi),nui*Amat(usr))*dx(0)

                # testf, nu, ansatzf, tr/long: +rrrt
                # testf, nu, ansatzf, t/l: -rirl
                # testf, nu, ansatzf, t/l: -riit
                # testf, nu, ansatzf, t/l: -rril
                # testf, nu, ansatzf, t/l: -rirt
                # testf, nu, ansatzf, t/l: +rrrl
                # testf, nu, ansatzf, t/l: -rrit
                # testf, nu, ansatzf, t/l: -riil

                # testf, nu, ansatzf, t/l: +irrt
                # testf, nu, ansatzf, t/l: -iirl
                # testf, nu, ansatzf, t/l: +iiit
                # testf, nu, ansatzf, t/l: +iril
                # testf, nu, ansatzf, t/l: -iirt
                # testf, nu, ansatzf, t/l: +irrl
                # testf, nu, ansatzf, t/l: +irit
                # testf, nu, ansatzf, t/l: +iiil


    if SIBC:
        if twolayer:
            d1=0.00008
            kappa1=6e9
            kappa2=1.8e8
            mu1=mu0
            mu2=mu0
            delta1=(2.0/(omega*mu1*kappa1))**0.5
            kz1=(1-1j)/delta1
            R=(mu1*kappa2/(mu2*kappa1))**0.5
            M=1+R
            N=1-R
            Zs=(1+1j)/(kappa1*delta1)*(M*np.exp(1j*kz1*d1)+N*np.exp(-1j*kz1*d1))/(M*np.exp(1j*kz1*d1)-N*np.exp(-1j*kz1*d1))
            Ys=1.0/Zs
            Yr=Ys.real
            Yi=Ys.imag
            print ("Yr: ", Yr )
            print ("Yi: ", Yi )
            kr=omega*mu0*Yi
            ki=-omega*mu0*Yr
        else:
            delta=(2.0/(omega*mu0*kappa_s))**0.5
            kr=-1.0/delta
            ki=-1.0/delta

        print (ds)
        print ("kr: ", kr )
        print ("ki: ", ki )
        def t():
            #n = -FacetNormal(mesh)
            return as_vector((-n[1],n[0]))

        a_boundary= \
            0.0*(-nu0*inner(vsr,(omega/(beta*c0))*inner(uvi,n))*ds(0)) -kr*nu0*inner(vsr,usr)*ds(0) + ki*nu0*inner(vsr,usi)*ds(0) \
            +0.0*(nu0*inner(vsi,(omega/(beta*c0))*inner(uvr,n))*ds(0)) -ki*nu0*inner(vsi,usr)*ds(0) - kr*nu0*inner(vsi,usi)*ds(0) \
                +kr*nu0*inner(dot(vvr,t()),dot(uvr,t()))*ds(0) - ki*nu0*inner(dot(vvr,t()),dot(uvi,t()))*ds(0) \
                +ki*nu0*inner(dot(vvi,t()),dot(uvr,t()))*ds(0) + kr*nu0*inner(dot(vvi,t()),dot(uvi,t()))*ds(0)
            #S1
            #S2
            #S3
            #S4
    else:
        a_boundary=0


    RHS=RHSsr*vsi*dx(0)  + RHSsi*vsi*dx(0)  + inner(RHSvr,vvr)*dx(0)  + inner(RHSvi,vvi) *dx(0)

    equation = a_curlcurl +a_kappa + a_epsilon +a_boundary == RHS

    #assemble(equation,exterior_facet_domains=boundaries)


    ##################################################################################
    Zero = Expression(('0','0','0','0','0','0'), degree=6)

    def u0_boundary(x, on_boundary):    # returns boolean if x on boundary
        return on_boundary

    ElectricBC=DirichletBC(V, Zero, u0_boundary)
    ####################################################################################

    ################
    # set_log_level(PROGRESS)
    if SIBC:
        solve(equation, Ecurl,solver_parameters={"linear_solver": "lu","preconditioner": "none"})
    else:
        solve(equation, Ecurl,ElectricBC,solver_parameters={"linear_solver": "mumps","preconditioner": "none"})
        # solve(equation, Ecurl,ElectricBC,solver_parameters={"linear_solver": "lu","preconditioner": "none"})
        # solve(equation, Ecurl,ElectricBC,solver_parameters={"linear_solver": "gmres","preconditioner": "sor"})
    print("solver done")
    #################

    [Ecurltr,Ecurlti,Ecurllr,Ecurlli]=Ecurl.split(deepcopy=True)

    ####################################################################
    #SIBC testing
    if(SIBC):
        SIBCtest_l=assemble(Ecurllr*ds(0))
        SIBCtest_t=assemble(inner(Ecurltr,my_Tangential(n))*ds(0))
        print ("SIBCtest_l: ", SIBCtest_l)
        print ("SIBCtest_t: ", SIBCtest_t)
    #####################################################################

    if(plot3Dflag):
        plot(Ecurltr,title='Ecurltr',basename='Ecurltr')
        plot(Ecurlti,title='Ecurlti',basename='Ecurlti')
        plot(Ecurllr,title='Ecurllr',basename='Ecurllr')
        plot(Ecurlli,title='Ecurlli',basename='Ecurlli')
        interactive()


    return [Ecurltr,Ecurlti,Ecurllr,Ecurlli]
#################################################################################################



###################################################################
def CurlCurl(mesh,omega, beta, epsilon, kappa, nu, RHSsr, RHSsi, RHSvr,RHSvi):
    Vtr = FiniteElement("Nedelec 1st kind H(curl)", mesh.ufl_cell(), curl_order)
    Vlr = FiniteElement("Lagrange", mesh.ufl_cell(), curl_long_order)
    Vli = Vlr
    Vti = Vtr

    V= FunctionSpace(mesh, MixedElement([Vtr, Vti, Vlr, Vli])) # 3D complex vector space

    (uvr,uvi,usr,usi) = TrialFunctions(V)
    (vvr,vvi,vsr,vsi) = TestFunctions(V)

    Ecurl = Function(V)
    #Ecurltr=Function(Vtr)
    #Ecurlti=Function(Vti)
    #Ecurllr=Function(Vlr)
    #Ecurlli=Function(Vli)

    """
    $
    A_\varepsilon=-\omega^2\left(
    \int_\Omega{v_v^\Re\cdot\varepsilon u_v^\Re d\Omega}+
    \int_\Omega{v_v^\Im\cdot\varepsilon u_v^\Im d\Omega}+
    \int_\Omega{v_s^\Re\cdot\varepsilon u_s^\Re d\Omega}+
    \int_\Omega{v_s^\Im\cdot\varepsilon u_s^\Im d\Omega}
    \right)
    $
    $
    A_\kappa=\omega\left(
    \int_\Omega{v_v^\Im\cdot\kappa u_v^\Re d\Omega}-
    \int_\Omega{v_v^\Re\cdot\kappa u_v^\Im d\Omega}+
    \int_\Omega{v_s^\Im\cdot\kappa u_s^\Re d\Omega}-
    \int_\Omega{v_s^\Re\cdot\kappa u_s^\Im d\Omega}
    \right)
    $
    $
    A_{\nabla}=\left(
    \int_\Omega{v_v^\Im\cdot\kappa u_v^\Re d\Omega}-
    \int_\Omega{v_v^\Re\cdot\kappa u_v^\Im d\Omega}+
    \int_\Omega{v_s^\Im\cdot\kappa u_s^\Re d\Omega}-
    \int_\Omega{v_s^\Re\cdot\kappa u_s^\Im d\Omega}
    \right)
    $
    """
    a_epsilon=-omega*omega*(inner(vvr,epsilon*uvr)*dx + inner(vvi,epsilon*uvi)*dx + inner(vsr,epsilon*usr)*dx + inner(vsi,epsilon*usi)*dx)
    a_kappa=omega*(inner(vvi,kappa*uvr)*dx - inner(vvr,kappa*uvi)*dx + inner(vsi,kappa*usr)*dx - inner(vsr,kappa*usi)*dx)

    a_curlcurl= -inner(vvr,-(omega/(beta*c0))**2 *nu*uvr)*dx -inner(vvr,nu*ZAmat(usi,omega,beta))*dx +inner(Bmat(vvr),nu*Bmat(uvr))*dx  \
                    -inner(Amat(vsr),nu*Zmat(uvi,omega,beta))*dx + inner(Amat(vsr),nu*Amat(usr))*dx \
                 -inner(vvi,-(omega/(beta*c0))**2 *nu*uvi)*dx +inner(vvi,nu*ZAmat(usr,omega,beta))*dx +inner(Bmat(vvi),nu*Bmat(uvi))*dx  \
                    +inner(Amat(vsi),nu*Zmat(uvr,omega,beta))*dx + inner(Amat(vsi),nu*Amat(usi))*dx

    RHS=RHSsr*vsi*dx + RHSsi*vsi*dx + inner(RHSvr,vvr)*dx + inner(RHSvi,vvi) *dx

    equation = a_curlcurl +a_kappa + a_epsilon == RHS

    ##################################################################################
    Zero = Constant(('0','0','0','0','0','0'))

    def u0_boundary(x, on_boundary):    # returns boolean if x on boundary
        return on_boundary

    ElectricBC=DirichletBC(V, Zero, u0_boundary)
    ####################################################################################

    ################
    # set_log_level(PROGRESS)
    # solve(equation, Ecurl,ElectricBC,solver_parameters={"linear_solver": "lu","preconditioner": "none"})
    solve(equation, Ecurl,ElectricBC,solver_parameters={"linear_solver": "mumps","preconditioner": "none"})
    # solve(equation, Ecurl,ElectricBC,solver_parameters={"linear_solver": "gmres","preconditioner": "ilu"})
    print("solver done")
    #################

    [Ecurltr,Ecurlti,Ecurllr,Ecurlli]=Ecurl.split(deepcopy=True)

    if(plot3Dflag):
        plot(Ecurltr,title='Ecurltr')
        plot(Ecurlti,title='Ecurlti')
        plot(Ecurllr,title='Ecurllr')
        plot(Ecurlli,title='Ecurlli')
        interactive()


    return [Ecurltr,Ecurlti,Ecurllr,Ecurlli]
