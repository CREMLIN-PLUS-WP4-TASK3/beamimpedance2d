from dolfin import *
def CplxPoisson(mesh,omega, beta, epsilon, kappa, Jszr,Jszi):
    H1 = FunctionSpace(mesh, "CG",div_long_order)    #Space for the Potential
    Mix=H1*H1
    
    (ur,ui) = TrialFunctions(Mix)
    (vr,vi) = TestFunctions(Mix)

    Phi=Function(Mix)
    Phir=Function(H1)
    Phii=Function(H1)
    
    ar=inner(grad(vr),epsilon*grad(ur))*dx +inner(grad(vr),(kappa/omega)*grad(ui))*dx
    ai=inner(grad(vi),epsilon*grad(ui))*dx -inner(grad(vi),(kappa/omega)*grad(ur))*dx
    RHSr=1/(beta*c0) *vr*Jszr*dx
    RHSi=1/(beta*c0) *vi*Jszi*dx
    
    eq= ar+ai==RHSr+RHSi
    
    Zero = Expression(('0.0','0.0'))
    

    def u0_boundary(x, on_boundary):    # returns boolean if x on boundary
        return on_boundary

    BC=DirichletBC(Mix, Zero, u0_boundary)
    
    
    set_log_level(PROGRESS)
    #solve(eq, Phi,BC,solver_parameters={"linear_solver": "gmres","preconditioner": "sor"})
    solve(eq, Phi,BC,solver_parameters={"linear_solver": "lu","preconditioner": "none"})
    (Phir,Phii)=Phi.split(deepcopy=False) 
    
    if(plot3Dflag):
        plot(Phir,title='Phir')
        plot(Phii,title='Phii')
        interactive()
    
    return [Phir,Phii]





"""
#This is obsolete!

def CplxMixedPoisson(mesh,omega, beta, epsilon, kappa, Jszr,Jszi):
    Hdiv = FunctionSpace(mesh, "BDM", 1)
    H1 = FunctionSpace(mesh, "CG", 1)
    Mix = MixedFunctionSpace([H1, H1, Hdiv, Hdiv,Hdiv,Hdiv])
    (uPhir,uPhii, uDr,uDi,uJr,uJi) = TrialFunctions(Mix)
    (vPhir,vPhii, vDr,vDi,vJr,vJi) = TestFunctions(Mix)
    
    
    sol=Function(Mix)
    #Phir=Function(H1)
    #Phii=Function(H1)
    #Dr=Function(Hdiv)
    #Di=Function(Hdiv)
    #Jr=Function(Hdiv)
    #Ji=Function(Hdiv)
    
    
    aSrhsr=1/(beta*c0) *inner(vPhir,Jszr)*dx
    aSrhsi=1/(beta*c0) *inner(vPhii,Jszr)*dx
    
    ar=-inner(grad(vPhir),uDr)*dx  -(1.0/omega)*inner(grad(vPhir),uJi)*dx
    ai=-inner(grad(vPhii),uDi)*dx  +(1.0/omega)*inner(grad(vPhii),uJi)*dx
    
    aJr= kappa*dot(div(vJr),uPhir)*dx
    aJi= kappa*dot(div(vJi),uPhii)*dx

    aDr= -inner(vDr,epsilon*grad(uPhir))*dx
    aDi= -inner(vDi,epsilon*grad(uPhii))*dx

    bDr= inner(vDr,uDr)*dx
    bDi= inner(vDi,uDi)*dx
    bJr= inner(vJr,uJr)*dx
    bJi= inner(vJi,uJi)*dx

    
    eq= ar+ai +aDr+aDi-bDr-bDi  +aJr+aJi-bJr-bJi== aSrhsr+aSrhsi
    
    
    Zero = Expression('0')
    VecZero=Expression(('0.0','0.0'))
    def u0_boundary(x, on_boundary):    # returns boolean if x on boundary
        return on_boundary

    bcPhir=DirichletBC(Mix.sub(0), Zero, u0_boundary)
    bcPhii=DirichletBC(Mix.sub(1), Zero, u0_boundary)
    bct=DirichletBC(Mix.sub(2), VecZero, u0_boundary)
    bctt=DirichletBC(Mix.sub(3), VecZero, u0_boundary)
    bcttt=DirichletBC(Mix.sub(4), VecZero, u0_boundary)
    bctttt=DirichletBC(Mix.sub(5), VecZero, u0_boundary)
    BC=([bcPhir,bcPhii,bct,bctt,bcttt,bctttt])
    
    set_log_level(PROGRESS)
    solve(eq, sol,BC,solver_parameters={"linear_solver": "lu","preconditioner": "none"})
    (Phir,Phii,Dr,Di,Jr,Ji)=sol.split(deepcopy=True) 
    
    if(plot3Dflag):
        plot(Phir,title='Phir')
        plot(Phii,title='Phii')
        plot(Dr,title='Dr')
        plot(Di,title='Di')
        plot(Jr,title='Jr')
        plot(Ji,title='Ji')
        interactive()
    
    return [Phir,Phii]
"""
