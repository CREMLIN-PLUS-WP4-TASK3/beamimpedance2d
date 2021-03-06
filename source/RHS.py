"""
FEniCS 2D beam coupling impedance simulation in frequency domain

by Uwe Niedermayer 2014
"""

from dolfin import *

class Beam(SubDomain):
   def inside(self, x, a):
       return x[0]*x[0] + x[0]*x[0] < a+DOLFIN_EPS

# """
# class MyFunction(Expression):
#     def eval(self, values, x):
#         r2 = x[0]*x[0] + x[1]*x[1]
#         if r2 < 0.5*0.5:
#             values[0] = 3.91
#         else:
#             values[1] = 0.0
# """

class Step(Expression):
    def eval(self, values, x):
        if x>0:
            values[0]=1
        else:
            values[0]=0

class Ring(Expression):
    def __init__(self,r1,r2,c1,c2):
        self.r1_=r1
        self.r2_=r2
        self.c1_=c1
        self.c2_=c2

    def eval(self,values,x):
        r=((x[0]-self.c1_)**2+(x[1]-self.c2_)**2)**0.5
        if (r>self.r1_ and r<self.r2_):
            values[0]=1.0
        else:
            values[0]=0.0

class RingDipole(UserExpression):
    def __init__(self,a,eps,xd,yd,*args,**kwargs):
        self.a_=a
        self.eps_=eps
        self.xd_=xd
        self.yd_=yd
        super().__init__(*args,**kwargs)

    def eval(self,values,x):
        assert(self.a_ != None)
        assert(self.eps_ != None)
        assert(self.xd_ != None)
        assert(self.yd_ != None)
        r=((x[0]-self.xd_)**2+(x[1]-self.yd_)**2)**0.5
        if near(r, 0.0):
            r = DOLFIN_EPS
        cosPhi=(x[0]-self.xd_)/r
        sinPhi=(x[1]-self.yd_)/r
        if (r>self.a_-self.eps_ and r<self.a_+self.eps_):
            if horizontal:
                values[0]=cosPhi
                #if (r<self.a_):
                #    values[0]=(r-self.a_+self.eps_)*cosPhi
                #else:
                #    values[0]=(r-self.a_+self.eps_)*cosPhi
            else:
                values[0]=(r-self.a_+self.eps_)*sinPhi
        else:
            values[0]=0.0

# class RingQuadrupole(Expression):
#     def __init__(self,a,eps,xd,yd):
#         self.a_=a
#         self.eps_=eps
#         self.xd_=xd
#         self.yd_=yd

#     def eval(self,values,x):
#         r=((x[0]-self.xd_)**2+(x[1]-self.yd_)**2)**0.5
#         cos2Phi=((x[0]-self.xd_)**2-(x[1]-self.yd_)**2)/r**2
#         sin2Phi=(x[1]-self.yd_)/r###do
#         Quad=(x[0]-self.xd_)**2-(x[1]-self.yd_)**2
#         if (r>self.a_-self.eps_ and r<self.a_+self.eps_):
#             if horizontal:
#                 values[0]=Quad
#                 #values[0]=cos2Phi+sin2Phi
#                 #if (r<self.a_):
#                 #    values[0]=(r-self.a_+self.eps_)*cosPhi
#                 #else:
#                 #    values[0]=(r-self.a_+self.eps_)*cosPhi
#             else:
#                 values[0]=0.0 #Implement skew quad
#         else:
#             values[0]=0.0


"""
###########################################################
#obsolete
def ExCurrent(mesh,subdomains,q,a):

    #beam=Beam()
    #markers = CellFunction("bool", mesh)
    #markers.set_all(False)
    #beam.mark(markers, True)

    #restriction = Restriction(markers,0)
    limit=1e6
    Jz=Expression('q/(pi*a*a)*(0.5*(0.0+2/pi*(atan(limit*(sqrt(x[0]*x[0]+x[1]*x[1])+a))-atan(limit*(sqrt(x[0]*x[0]+x[1]*x[1])-a)))))' , q=q, limit=limit,a=a)    #rect
    #V_Ex=FunctionSpace(restriction,'CG',1)
    #V=FunctionSpace(mesh,'CG',div_long_order)
    #Jz=project(Jz,V)




    #for cell in cells(mesh)
    #    for vertex in vertices(cell):
    #    s#ubdomain_no = subdomains.array()[cell_number]
    #    if vertex  !=75:
    #        Jz.vector()[vertex]=0.0

    return Jz
##################################################################
"""

def ExCurrentShifted(mesh,subdomains,q,a,delta_x,delta_y):
    """
    @param q[Float] charge
    @param a[Float] beam radius
    $$\frac{1}{\pi}
        \arctan{
                \left(
                        C_{limit}
                        \left(
                                a-\sqrt{\left(x-\Delta x\right)^2 + \left(y-\Delta y\right)^2}
                        \right)
                \right)
        }
        + 1/2
    $$
    """
    Climit=1e6
    #Jz=Expression('q/(pi*a*a)*   (1.0/pi *atan(limit*(a-sqrt(pow(x[0]-delta_x,2.0)+pow(x[1]-delta_y,2.0) ))) +0.5)  ' , q=q, limit=limit,a=a,delta_x=delta_x,delta_y=delta_y)    #normalized current desity
    Jz=Expression('1.0/pi *atan(Climit*(a-sqrt(pow(x[0]-delta_x,2.0)+pow(x[1]-delta_y,2.0) ))) +0.5',
                  q=q,
                  Climit=Climit,
                  a=a,
                  delta_x=delta_x,
                  delta_y=delta_y,
                  degree=2)    #uniform current desity
    #######################################
    #calculate proper area for normalzation
    VArea=FunctionSpace(mesh,'Discontinuous Lagrange',0)
    AreaFunction=project(Jz,VArea)
    Area=assemble(AreaFunction*dx)
    print ("Excurrentshifted Area: ", Area)
    #######################################

    return Jz/Area


def ExCurrentShiftedDipole(mesh,subdomains,q,a,eps,xd,yd):
    #######################################
    #calculate proper dipole moment for normalzation
    Vdip=FunctionSpace(mesh,'Lagrange',1)


    Source=RingDipole(a,eps,xd,yd,element=Vdip.ufl_element())
    # Source=UserExpression(ring_dipole_code)


    #DipEx=Expression('x[0]-xd',xd=xd)
    if horizontal:
        DipTest=Expression('x[0]-xd',xd=xd, degree=2)
    else:
        DipTest=Expression('x[1]-yd',yd=yd, degree=2)

    #DipExFunction=project(Source,Vdip)
    #Dipolemoment=assemble(DipExFunction*DipTest*dx)
    DipTestFunction=project(Source,Vdip)
    DipolemomentTEST=assemble(DipTestFunction*DipTest*dx)
    if DipolemomentTEST==0:
        DipolemomentTEST = DOLFIN_EPS
    print("Dipole Moment: " , DipolemomentTEST)
    #######################################

    return Source/DipolemomentTEST


# def ExCurrentShiftedQuadrupole(mesh,subdomains,q,a,eps,xd,yd):
#     Source=RingQuadrupole(a,eps,xd,yd)


#     #DipEx=Expression('x[0]-xd',xd=xd)

#     QuadTest=Expression('(x[0]-xd)*(x[0]-xd)-(x[1]-yd)*(x[1]-yd)',xd=xd, yd=yd)


#     #######################################
#     #calculate proper dipole moment for normalzation
#     Vquad=FunctionSpace(mesh,'CG',1)

#     #DipExFunction=project(Source,Vdip)
#     #Dipolemoment=assemble(DipExFunction*DipTest*dx)
#     QuadTestFunction=project(Source,Vquad)
#     QuadrupolemomentTEST=assemble(QuadTestFunction*QuadTest*dx)
#     print("Quadrupole Moment: " , QuadrupolemomentTEST)
#     #######################################

#     return Source/QuadrupolemomentTEST
