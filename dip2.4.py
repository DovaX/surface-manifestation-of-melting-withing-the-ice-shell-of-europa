import os
from dolfin import *
import math
import numpy

"""2.3 - puvodni conditional strain-rate, 2.4 harmonicky prumer, strain-weakening - plastic strain"""


# Uflacs representation scales better with form complexity
# parameters['form_compiler']['representation'] = 'uflacs'
comm = mpi_comm_world()
rank = MPI.rank(comm)
# set_log_level(WARNING if rank==0 else WARNING+1)


parameters["ghost_mode"] = "shared_facet"
# CASES
"""
var2 = int(input("1) Blankenbach 2) Poorman ALE 3)Viscoplastic benchmark\n"))
if var2==1:
    var = int(input("1) Case 1, 2) Case 2a: ,3) Case 2b\n"))
elif var2==3:
    var = int(input("1) Case 1, 2) Case 2, 3) Case 3, 4) Case 4, 5) Case 5\n"))
else:
    var = 1
"""

#Parametry behu
var2 =4
var = 1

picard = True
results_number = 85 #Number of the directory

harmonic=True #pridani harmonickeho prumeru pro eff viskozitu, deformacni zmekceni plastic strain


disc = 400 # Number of elements in x-axis

angle = 3


if var2==4 or var2==5:    
    write_data_interval = 1  # 20 #will print data once per x steps
else:
    write_data_interval = 10

lmin = 0.0
hmin = 0.0
lmax = 1.0  # width
hmax = 1.0  # height
t_zero = True  # first step is t=0 - plastic part of viscosity is substituted by linear part
t = 0.0 # time variable

# Case 1
if (var == 1 or var == 2 or var == 3 or var == 4 or var == 5):
    lmax = 1.0  # width
    hmax = 1.0  # height
    if (var2 == 1 and var == 3):
        x = 2.5
        x = 1.0  # ????
    # Equation coefficients
    alpha = 2.5e-5  # rozmerna teplotni roztaznost
    Ra = 1e4  # Raighley number
    eta0 = 1.0  # viskozitni parametr
    # Time-stepping parameters
    if var2!=4 and var2!=5:
        t_end = 1.0
    else:
        t_end = 5.0
    dt = Constant(0.0001)  # Constant(0.0005)
    theta = Constant(1.0)  # Crank-Nicolson scheme
    ez = Constant((0.0, 1.0))
    k = 1.0  # difuzivita
    deltaT = 1000.0  # teplota na dolni hranici
    if var2 == 1 or var2 == 2:
        if (var == 1):
            def eta(T, x):
                return(eta0)
        if (var == 2):
            b = math.log(1000)
            c = 0.0

            def eta(T, x):
                return(eta0*exp(-b*T/deltaT+c*(1-x[1])/h))
        if (var == 3):
            b = math.log(16384)
            c = math.log(64)

            def eta(T, x):
                return(eta0*exp(-b*T/deltaT+c*(1-x[1])/h))

    if (var2 == 3):
        Ra = 1e2
        def dot_epsilon(v):

            return(sqrt(inner(sym(grad(v)), sym(grad(v)))))
            # return(sym(grad(v)))

        def eta_lin(T, x):
            gamma_T = ln(1e5)
            if var > 2:
                gamma_y = ln(10)  # case 3,4,5a,5b
            else:
                gamma_y = ln(1)  # case 1,2
            return(exp(-gamma_T*T+gamma_y*(1-x[1])))

        def eta_plast(v, x):  # depends on \dot{epsilon}
            if var == 5:
                yield_stress = 4.0
            else:
                yield_stress = 1.0
            eta_star = 0.001
            # yield_stress=1.0

            # epsilon_doubledot_epsilon = dot_epsilon(v)#sqrt(inner(dot_epsilon(v),dot_epsilon(v)))
            # print(project(dot_epsilon(v),Pspace)(0.5,0.5))
            #print(str(epsilon_doubledot_epsilon(1.0,1.0)), "EDE")
            # print(epsilon_doubledot_epsilon(0.0,1.0))
            return(eta_star+yield_stress/dot_epsilon(v))

        def eta_tzero(T, x, v):
            eta_l = eta_lin(T, x)
            eta_p = eta_lin(T, x)

            return(2*(eta_l*eta_p)/(eta_l+eta_p))

        def eta(T, x, v):
            eta_l = eta_lin(T, x)
            eta_p = eta_plast(v, x)
            return(2*(eta_l*eta_p)/(eta_l+eta_p))


if var2 == 4 or var2==5:

    eps_0 = 0.0
    eps_inf = 0.1
    C_0 = 400
    C_inf = 20
    rho = 2700
    g = 1.0
    angle_phi = math.pi/18*angle
    lmin = -2.0
    lmax = 2.0
    hmin = 0.0
    hmax = 1.0
    eta_bg = 100000
    eta_inclusion = 1
    lmin_inclusion = -0.02
    lmax_inclusion = 0.02
    hmin_inclusion = 0.0
    hmax_inclusion = 0.02

    def C_fc(eps_II):
        C1 = conditional(le(eps_II, eps_0), C_0, C_0 +
                         (C_inf-C_0)*(eps_II-eps_0)/(eps_inf-eps_0))
        C2 = conditional(le(eps_inf, eps_II), C_inf, C1)
        return(C2)

    def eta_viscous(x):
        eta_c1 = conditional(
            le(abs(x[1]), hmax_inclusion), eta_inclusion, eta_bg)
        eta_visc = conditional(le(abs(x[0]), lmax_inclusion), eta_c1, eta_bg)
        return(eta_visc)
    



    def sigma(x, v):
        sigma = eta_viscous(x)*(grad(v)+grad(v).T)
        return(sigma)

    def sigma_II(x, v):
        sigma_invariant = second_invariant(sigma(x, v))
        return(sigma_invariant)

    # plotting
    # sigma_min = 0.1, pridat (2.30 formula)
    def yield_str(p, eps_II):

        sigma_min = 0.1
        y_stress= p*math.sin(angle_phi)+C_fc(eps_II)*math.cos(angle_phi)
        result=conditional(le(sigma_min, y_stress), y_stress, sigma_min)
        return(result)

    def second_invariant(dot_eps):

        dot_eps_II = sqrt(inner(dot_eps, dot_eps)/2)

        # old - principal invariant - we use main invariants instead - because principal invariant can be negative
        # dot_eps_I = tr(dot_eps)
        # dot_eps_II = (dot_eps_I**2 - tr(dot_eps*dot_eps))/2

        return(dot_eps_II)

    def plastic_strain_rate(p, eps_II,v,x,t):
        dot_eps_plast = sigma(x,v)/2/eta_eff(p,eps_II,v,x,t,True)-sigma(x,v)/2/eta_viscous(x)
        return(dot_eps_plast)


    def eta_eff(p, eps_II, v, x, t,harmonic):
        yield_stress = yield_str(p, eps_II)
        sigma_invariant = sigma_II(x, v)
        # change this part
        # dot_eps_II = conditional(le(second_invariant(grad(v)+grad(v).T),0),0.1,second_invariant(grad(v)+grad(v).T))
        dot_eps_II = second_invariant(grad(v)+grad(v).T)
        if harmonic:
            eta_effective = 1/(1/eta_viscous(x)+1/(yield_stress/(2*dot_eps_II)))
        else:
            eta_effective = conditional(
            le(yield_stress, sigma_invariant), yield_stress/(2*dot_eps_II), eta_viscous(x))
        
    
        return(eta_effective)


# Create mesh and build function space
xdiscretization = disc
ydiscretization = int(xdiscretization*(hmax-hmin)/(lmax-lmin))

if var2 == 4 or var2 == 5:
    mesh = RectangleMesh(Point(lmin, hmin), Point(
        lmax, hmax), xdiscretization, ydiscretization, "crossed")
else:
    mesh = UnitSquareMesh(xdiscretization, ydiscretization, 'crossed')
x = SpatialCoordinate(mesh)  # create coordinates variable


# ========================== PARAMETERS FILE ==========================
if var2 == 1:
    name = "Blankenbach"
if var2 == 2:
    name = "ALE"
if var2 == 3:
    name = "Viscoplastic"
if var2 == 4:
    name = "Plastic"
if var2 == 5:
    name = "ALE_Plastic"


directory = name+"_results_c"+str(var)+"_"+str(results_number)
if picard:
    directory += "_picard"

try:
    if not os.path.exists(os.path.dirname(directory+"/setup.txt")):
        os.makedirs(os.path.dirname(directory+"/setup.txt"))
except:
    pass

file_setup = open(directory+"/setup.txt", "w")
file_setup.write(str(name))
file_setup.write(str("var"+str(var)))
try:
    file_setup.write(str("var2"+str(var2)))
except:
    pass
file_setup.write(str("lmax,hmax"+str(lmax)+"\t"+str(hmax)))
file_setup.write(str("alpha"+str(alpha)))
file_setup.write(str("Ra"+str(Ra)))
file_setup.write(str("t_end"+str(t_end)))
file_setup.write(str("eta0"+str(eta0)))
file_setup.write(str("theta"+str(theta)))
file_setup.write(str("deltaT"+str(deltaT)))
file_setup.write(str("initial dt"+str(dt)))
file_setup.flush()
file_setup.close()
# ========================== END PARAMETERS FILE ==========================


""" Function spaces, elements"""
Velement = VectorElement("CG", mesh.ufl_cell(),  2)  # .ufl_cell()
Pelement = FiniteElement("CG", mesh.ufl_cell(), 1)
Telement = FiniteElement("CG", mesh.ufl_cell(), 2)
Vspace = FunctionSpace(mesh,  Velement)
Pspace = FunctionSpace(mesh, Pelement)
Tspace = FunctionSpace(mesh, Telement)
VPspace = FunctionSpace(mesh, MixedElement([Velement, Pelement]))

rms_fc = Function(Pspace)

if var2 == 2 or var2 == 5:
    HVelement = VectorElement("CG", mesh.ufl_cell(), 1)
    HVspace = FunctionSpace(mesh, HVelement)

if var2 == 4 or var2==5:
    EPSelement = FiniteElement("DG", mesh.ufl_cell(), 0)
    EPSspace = FunctionSpace(mesh, EPSelement)


# Create boundary markers
boundary_parts = FacetFunction('size_t', mesh)
top = AutoSubDomain(lambda x: near(x[1], hmax))
bottom = AutoSubDomain(lambda x: near(x[1], hmin))
left = AutoSubDomain(lambda x: near(x[0], lmin))
right = AutoSubDomain(lambda x: near(x[0], lmax))

if var2 == 4 or var2 == 5:
    top.mark(boundary_parts, 1)
    bottom.mark(boundary_parts, 2)
    left.mark(boundary_parts, 3)
    right.mark(boundary_parts, 4)
else:
    top.mark(boundary_parts, 1)
    bottom.mark(boundary_parts, 2)
    left.mark(boundary_parts, 3)
    right.mark(boundary_parts, 3)


# left.mark(boundary_parts, 4)
# Initial condition and right-hand side


if var2 == 3:
    T_init = Expression("""1-x[1]+cos(pi*x[0])*sin(pi*x[1])*0.01""", degree=2)
else:
    # *0.01 slower beginning
    T_init = Expression('(1-x[1])+sin(2*pi*x[0])*sin(pi*x[1])*0.01', degree=2)

if var2 == 4 or var2 == 5:
    # eps_init = Expression('(x[0]+2)*0.01',degree=2)
    # eps_init = Expression('x[1]-x[0]*x[0]*0.01',degree=2)
    eps_init = Expression('0.0', degree=2)


h_init = Expression(("0.0", "0.0"), hmax=hmax, hmin=hmin, degree=1)

# boundary conditions
freeslip = Constant(0)
bc_top = DirichletBC(VPspace.sub(0).sub(1), freeslip, boundary_parts, 1)
bc_bottom = DirichletBC(VPspace.sub(0).sub(1), freeslip, boundary_parts, 2)
bc_sides = DirichletBC(VPspace.sub(0).sub(0), freeslip, boundary_parts, 3)
if (var2 == 1 or var2 == 3):  # zafixovany tlak pouze pro Blankenbacha
    bc_pressure = DirichletBC(VPspace.sub(1), Expression(
        "0.0", degree=1), "near(x[0],0.0) && near(x[1],1.0)", method="pointwise")
    bcs = [bc_top, bc_bottom, bc_sides, bc_pressure]
    print(1)
else:
    bcs = [bc_bottom, bc_sides]
    print(2)
    # bcs_h = []

if var2 == 4 or var2==5:
    # bc_pressure = DirichletBC(VPspace.sub(1), Expression("0.0",degree=1), "near(x[0],0.0) && near(x[1],1.0)", method="pointwise")
    bc_left = DirichletBC(VPspace.sub(0).sub(
        0), Constant(2.0), boundary_parts, 3)
    bc_right = DirichletBC(VPspace.sub(0).sub(
        0), Constant(-2.0), boundary_parts, 4)
    bcs = [bc_left, bc_right, bc_bottom]  # ,bc_pressure]

    bcEPS1 = DirichletBC(EPSspace, Constant(0.0), boundary_parts, 3)
    bcEPS2 = DirichletBC(EPSspace, Constant(0.0), boundary_parts, 4)
    bcEPS = [bcEPS1, bcEPS2]


bcT1 = DirichletBC(Tspace, Constant(0.0), boundary_parts, 1)
bcT2 = DirichletBC(Tspace, Constant(1.0), boundary_parts, 2)
bcT = [bcT1, bcT2]




class Omega(SubDomain):
    def inside(self, x, on_boundary):
        return True


omega = Omega()

# Define boundary measure on Neumann part of boundary
# dsN = Measure("ds", subdomain_id=3, subdomain_data=boundary_parts)
# Collect boundary conditions
n = FacetNormal(mesh)
I = Identity(mesh.geometry().dim())
ds = Measure("ds", subdomain_data=boundary_parts)
dS = Measure("dS")  # (subdomain_data=boundary_parts)


# Functions
T0 = Function(Tspace)
phi = TestFunction(Tspace)
T_ = TrialFunction(Tspace)
psi, ksi = TestFunctions(VPspace)
w = Function(VPspace)
w_ = TrialFunction(VPspace)
T = Function(Tspace)
v, p = split(w)
v_, p_ = split(w_)


# Mesh functions
if var2 == 2 or var2==5:
    h = Function(HVspace)
    h_k = Function(HVspace)
    h_ = TrialFunction(HVspace)
    htest = TestFunction(HVspace)
    dh = Function(HVspace)

if var2 == 3 or var2 == 4 or var2 == 5:
    eta_fc = Function(Pspace)


if var2 == 4 or var2==5:
    eps_II_prev = Function(EPSspace)
    eps_II = Function(EPSspace)
    eps_II_ = TrialFunction(EPSspace)
    test_eps = TestFunction(EPSspace)
    eps_II_k = Function(EPSspace)
    doteps_II = Function(EPSspace)
    v_diff = Function(Vspace)
    plastic_dot_eps = Function(EPSspace)
# ========================================================
# Free surface problem
# ========================================================
# delta h free surface bcs
if(var2 == 2 or var2==5):

    bc_h_bot = DirichletBC(HVspace.sub(1), Constant(0.0), boundary_parts, 2)
    bcs_hx = DirichletBC(HVspace.sub(0), Constant(0.0), omega)
    bcs_h = [bc_h_bot, bcs_hx]

    # delta_h predictor
    gamma_h = Constant(0.005/xdiscretization)
    F_h_pr = inner(grad(h_), grad(htest))*dx - inner(grad(h_[1]), n)*htest[1]*ds(1) \
        - (h_[1] - h_k[1] + dt*(v[0]*h_[1].dx(0) - v[1])) * \
        (inner(grad(htest[1]), n) - htest[1]/gamma_h)*ds(1)
    lhs_h_pr = lhs(F_h_pr)
    rhs_h_pr = rhs(F_h_pr)
# ========================================================


if (var2 == 1 or var2 == 3):
    F1 = (1.0/dt)*(T_-T0)*phi*dx + theta*dot(v, grad(T_))*phi*dx + (1.0-theta)*dot(v, grad(T0)) * \
        phi*dx + k*dot(grad(T_), grad(phi))*theta*dx + k * \
        dot(grad(T0), grad(phi))*(1-theta)*dx
if (var2 == 1):
    F2 = inner(eta(T, x)*(grad(v_)+transpose(grad(v_))), grad(psi)) * \
        dx-p_*div(psi)*dx-ksi*div(v_)*dx - Ra*T*dot(psi, ez)*dx

if (var2 == 3):
    # Only eta definition varies - doplnit
    # F1 = (1.0/dt)*(T-T0)*phi*dx + theta*dot(v, grad(T))*phi*dx + (1.0-theta)*dot(v, grad(T0))*phi*dx + k*dot(grad(T),grad(phi))*theta*dx +k*dot(grad(T0),grad(phi))*(1-theta)*dx

    # Otazka - ma to byt tu, nebo pred teplotni rovnici?
    F2_tzero = inner(eta_tzero(T, x, v)*(grad(v)+transpose(grad(v))),
                     grad(psi))*dx-p*div(psi)*dx-ksi*div(v)*dx - Ra*T*dot(psi, ez)*dx
    F2 = inner(eta(T, x, v)*(grad(v)+transpose(grad(v))), grad(psi)) * \
        dx-p*div(psi)*dx-ksi*div(v)*dx - Ra*T*dot(psi, ez)*dx
#   F3 = (eta_fc - eta(T,x,v))*dx
if var2 == 2:  # STABILIZACE
    v_cor = dh/dt
    lam = Constant(1.0)
    F1 = (1.0/dt)*(T_-T0)*phi*dx + theta*dot(v-v_cor, grad(T_))*phi*dx +  \
        (1.0-theta)*dot(v-v_cor, grad(T0))*phi*dx +  \
        k*dot(grad(T_), grad(phi))*theta*dx + \
        k*dot(grad(T0), grad(phi))*(1-theta)*dx
    F2 = inner(eta(T, x)*(grad(v_)+transpose(grad(v_))), grad(psi))*dx \
        - p_*div(psi)*Ra/(alpha*deltaT)*dx-ksi*div(v_)*dx \
        + Ra*(1/(alpha*deltaT)-T)*dot(psi, ez)*dx \
        + lam*dt*Ra*(1/(alpha*deltaT)-T)*dot(ez, psi)*inner(v_, n)*ds(1)


# ========================== PROBLEM SOLVERS ==========================


if var2 == 3:
    # DF1 = derivative(F1, T)
    # problem1 = NonlinearVariationalProblem(F1, T, bcT,DF1) #DF1 a DF2 - Jacobians
    # solver1  = NonlinearVariationalSolver(problem1)

    DF2 = derivative(F2, w)  # Jacobian
    problem2 = NonlinearVariationalProblem(F2, w, bcs, DF2)
    solver2 = NonlinearVariationalSolver(problem2)

    DF2_tzero = derivative(F2_tzero, w)  # Jacobian
    problem2_tzero = NonlinearVariationalProblem(F2_tzero, w, bcs, DF2_tzero)
    solver2_tzero = NonlinearVariationalSolver(problem2_tzero)

    solver2.parameters['newton_solver']['absolute_tolerance'] = 1e-6
    solver2.parameters['newton_solver']['maximum_iterations'] = 100
elif var2 == 2 or var2 == 1:

    problem2 = LinearVariationalProblem(lhs(F2), rhs(F2), w, bcs)
    solver2 = LinearVariationalSolver(problem2)

# Prepare solution function and solver

if var2 != 4 and var2!=5:
    problem1 = LinearVariationalProblem(lhs(F1), rhs(F1), T, bcT)
    solver1 = LinearVariationalSolver(problem1)

if(var2 == 2) or var2==5:
    problem_h_pr = LinearVariationalProblem(lhs_h_pr, rhs_h_pr, h, bcs_h)
    solver_h_pr = LinearVariationalSolver(problem_h_pr)
    solver_h_pr.parameters['linear_solver'] = "cg"


# ========================== PICARD ITERATIONS ==========================

# ###### Tosi et al. #####

if picard and var2 == 3:
    # assigner = FunctionAssigner(VPspace, Vspace*Pspace)
    # dot_eps_II = grad(v)+grad(v).T
    # v_k = Function(VPspace.sub(0))
    # w_const=Constant((0.0,0.0,0.0))
    # w_k=interpolate(w_const, VPspace)

    # w_k.assign(w)
    # v_k=w.split(deepcopy=True)[0]
    # v_k0=w_k.split(deepcopy=True)[0]

    w_k = Function(VPspace)
    v_k, p_k = split(w_k)

    F2_picard = inner(eta(T, x, v_k)*(grad(v_)+transpose(grad(v_))),
                      grad(psi))*dx-p_*div(psi)*dx-ksi*div(v_)*dx - Ra*T*dot(psi, ez)*dx
    problem2_picard = LinearVariationalProblem(
        lhs(F2_picard), rhs(F2_picard), w, bcs)  # ,DF2_picard)
    solver2_picard = LinearVariationalSolver(problem2_picard)

    F2_picard_tzero = inner(eta_tzero(T, x, v_k)*(grad(v_)+transpose(grad(v_))),
                            grad(psi))*dx-p_*div(psi)*dx-ksi*div(v_)*dx - Ra*T*dot(psi, ez)*dx
    problem2_picard_tzero = LinearVariationalProblem(
        lhs(F2_picard_tzero), rhs(F2_picard_tzero), w, bcs)  # ,DF2_picard)
    solver2_picard_tzero = LinearVariationalSolver(problem2_picard_tzero)


# ========================== Maierova Benchmark ==========================


if var2 == 4:

    w_k = Function(VPspace)
    v_k, p_k = split(w_k)

    vn = abs(dot(v("+"), n("+")))/2.0
    F_cell = -dot(grad(test_eps), v*eps_II_)*dx
    F_int = dot(jump(test_eps, n), jump(eps_II_, n))*vn*dS + \
        dot(v("+"), jump(test_eps, n))*avg(eps_II_)*dS
    F_outflow = dot(v,n)*eps_II_*test_eps*ds(1)#dot(v, jump(test_eps, n))*avg(eps_II_)*ds(1)#(1)

    F1_harmonic = (1.0/dt)*(eps_II_-eps_II_prev)*test_eps*dx + F_cell + \
    F_int - second_invariant(plastic_strain_rate(p, eps_II_prev,v,x,t))*test_eps*dx + F_outflow
    """nahradit eps_II_prev uvnitr plastic strain rate za epsilon"""
    
    
    #F1_harmonic = (1.0/dt)*(eps_II_-eps_II_prev)*test_eps*dx + F_cell + \
    #F_int - second_invariant((grad(v)+grad(v).T)/2)*test_eps*dx + F_outflow
    
    
    F1 = (1.0/dt)*(eps_II_-eps_II_prev)*test_eps*dx + F_cell + \
    F_int - second_invariant((grad(v)+grad(v).T)/2)*test_eps*dx + F_outflow
    
    # F1 = (1.0/dt)*(eps_II_-eps_II_prev)*test_eps*dx + dot(v,grad(eps_II_))*test_eps*dx - second_invariant(grad(v)+grad(v).T)*test_eps*dx
    F2_picard = inner(eta_eff(p, eps_II, v_k, x, t,False)*(grad(v_)+transpose(grad(v_))),
                      grad(psi))*dx-p_*div(psi)*dx-ksi*div(v_)*dx + rho*g*dot(psi, ez)*dx

    F2_picard_harmonic = inner(eta_eff(p, eps_II, v_k, x, t,True)*(grad(v_)+transpose(grad(v_))),
                      grad(psi))*dx-p_*div(psi)*dx-ksi*div(v_)*dx + rho*g*dot(psi, ez)*dx
  

    problem1 = LinearVariationalProblem(lhs(F1), rhs(F1), eps_II, bcEPS)
    solver1 = LinearVariationalSolver(problem1)
    
    problem2_picard = LinearVariationalProblem(
        lhs(F2_picard), rhs(F2_picard), w, bcs)  # ,DF2_picard)
    solver2_picard = LinearVariationalSolver(problem2_picard)
    
    solver1.parameters['linear_solver'] = "mumps"
    solver2_picard.parameters['linear_solver'] = "mumps"

# =========================== Shearbands with ALE ==============================

if var2 == 5:


    v_cor = dh/dt
    lam = Constant(1.0)
    #F1 = (1.0/dt)*(T_-T0)*phi*dx + theta*dot(v-v_cor, grad(T_))*phi*dx +  \
    #    (1.0-theta)*dot(v-v_cor, grad(T0))*phi*dx +  \
    #    k*dot(grad(T_), grad(phi))*theta*dx + \
    #    k*dot(grad(T0), grad(phi))*(1-theta)*dx
    #F2 = inner(eta(T, x)*(grad(v_)+transpose(grad(v_))), grad(psi))*dx \
    #    - p_*div(psi)*Ra/(alpha*deltaT)*dx-ksi*div(v_)*dx \
    #    + Ra*(1/(alpha*deltaT)-T)*dot(psi, ez)*dx \
    #    + lam*dt*Ra*(1/(alpha*deltaT)-T)*dot(ez, psi)*inner(v_, n)*ds(1)
    
    w_k = Function(VPspace)
    v_k, p_k = split(w_k)

    vn = abs(dot(v("+")-v_cor("+"), n("+")))/2.0
    F_cell = -dot(grad(test_eps), (v-v_cor)*eps_II_)*dx
    F_int = dot(jump(test_eps, n), jump(eps_II_, n))*vn*dS + \
        dot(v("+")-v_cor("+"), jump(test_eps, n))*avg(eps_II_)*dS
    F_outflow = dot(v-v_cor,n)*eps_II_*test_eps*ds(1)#dot(v, jump(test_eps, n))*avg(eps_II_)*ds(1)#(1)

    F1 = (1.0/dt)*(eps_II_-eps_II_prev)*test_eps*dx + F_cell + \
        F_int - second_invariant((grad(v)+grad(v).T)/2)*test_eps*dx + F_outflow
    
    # F1 = (1.0/dt)*(eps_II_-eps_II_prev)*test_eps*dx + dot(v,grad(eps_II_))*test_eps*dx - second_invariant(grad(v)+grad(v).T)*test_eps*dx
    F2_picard = inner(eta_eff(p, eps_II, v_k, x, t)*(grad(v_)+transpose(grad(v_))),
                      grad(psi))*dx-p_*div(psi)*dx-ksi*div(v_)*dx + rho*g*dot(psi, ez)*dx+lam*dt*rho*g*dot(ez, psi)*inner(v_, n)*ds(1)

    problem1 = LinearVariationalProblem(lhs(F1), rhs(F1), eps_II, bcEPS)
    solver1 = LinearVariationalSolver(problem1)
    
    problem2_picard = LinearVariationalProblem(
        lhs(F2_picard), rhs(F2_picard), w, bcs)  # ,DF2_picard)
    solver2_picard = LinearVariationalSolver(problem2_picard)
    
    solver1.parameters['linear_solver'] = "mumps"
    solver2_picard.parameters['linear_solver'] = "mumps"






# Prepare initial condition
T0.interpolate(T_init)

if var2 == 4 or var2 == 5:
    eps_II_prev.interpolate(eps_init)

if (var2 == 2) or var2 == 5:
    h.interpolate(h_init)
    h_k.interpolate(h_init)







def open_XDMF_files(var2, directory):
    file_h = 0
    file_eta = 0
    file_eps_II = 0
    file_doteps = 0

    # Create file for storing results
    fT = XDMFFile(directory+"/T.xdmf")
    fV = XDMFFile(directory+"/V.xdmf")
    fP = XDMFFile(directory+"/P.xdmf")
    fT.parameters["flush_output"] = True
    fV.parameters["flush_output"] = True
    fP.parameters["flush_output"] = True

    if var2 == 2 or var2 == 5:
        file_h = XDMFFile(comm, directory+"/h.xdmf")
        file_h.parameters["flush_output"] = True
        file_h.parameters["rewrite_function_mesh"] = True
        h.rename("h", "h")
        file_h.write(h, 0)

    if var2 == 3 or var2 == 4 or var2 == 5:
        file_eta = XDMFFile(comm, directory+"/eta.xdmf")
        file_eta.parameters["flush_output"] = True

        eta_fc.rename("eta_fc", "eta_fc")
        file_eta.write(eta_fc, 0)

    if var2 == 4 or var2 == 5:
        file_eps_II = XDMFFile(comm, directory+"/eps_II.xdmf")
        file_eps_II.parameters["flush_output"] = True

        eps_II.rename("eps_II", "eps_II")
        file_eps_II.write(eps_II, 0)

        file_doteps = XDMFFile(comm, directory+"/doteps_II.xdmf")
        file_doteps.parameters["flush_output"] = True

        doteps_II.rename("doteps_II", "doteps_II")
        file_doteps.write(doteps_II, 0)

    file_rms = XDMFFile(comm, directory+"/rmsvel.xdmf")
    file_rms.parameters["flush_output"] = True

    rms_fc.rename("rms_vel", "rms_vel")
    file_rms.write(rms_fc, 0)
    return(fT, fV, fP, file_h, file_eta, file_eps_II, file_rms, file_doteps)


fT, fV, fP, file_h, file_eta, file_eps_II, file_rms, file_doteps = open_XDMF_files(
    var2, directory)

# Time-stepping
if var2 == 4 or var2 == 5:
    eps_II.interpolate(eps_init)


T.interpolate(T_init)
T.rename("T", "temperature")

# save initial solution
fT.write(T)

if var2 == 4 or var2 == 5:
    file_eps_II.write(eps_II)
    file_doteps.write(doteps_II)

# computation of time step using CFL criterion


def Compute_dt(uvm):
    
    h_min = MPI.min(comm, mesh.hmin())
    vm_array = uvm.vector().get_local()  # .array()
    vm_max = numpy.abs(vm_array).max()
    vm_max2 = MPI.max(comm, vm_max)
    if var2==4 or var2==5:
    
        C_CFL = 0.1      # value of CFL - it was 0.5   #0.3,0.1 vyzkouset
    else:
        C_CFL = 0.5
    dt = (C_CFL*h_min/vm_max2)
    # dt=0.001
    return (dt)


def nusselt(T, hmax, n):
    int1 = assemble(dot(grad(T), n)*ds(1))
    int2 = assemble(T*ds(2))
    return(-hmax*int1/int2)


def nusselt_top(T, n):
    int1 = -assemble(dot(grad(T), n)*ds(1))
    return(int1)


def nusselt_bot(T, n):
    int1 = -assemble(dot(grad(T), -n)*ds(2))
    return(int1)


def mean_T(T):
    meant = assemble(T*dx)
    return(meant)


def depth_minmax(value_depth_list):
    """averages values with same y coordinate - input in format [[value1,y1],[value2,y2],...]"""
    """returns dict1 (lateral maxima) and dict2 (lateral minima)"""

    list2 = set([x[1] for x in value_depth_list])
    dict1 = {depth: 0 for depth in list2}

    # DEPTH MAXIMUM
    for i in range(len(value_depth_list)):
        if dict1[value_depth_list[i][1]] < value_depth_list[i][0]:
            dict1[value_depth_list[i][1]] = value_depth_list[i][0]

    overall_max = max(list(dict1.values()))
    dict2 = {depth: overall_max for depth in list2}  # initialize with maximum

    # DEPTH MINIMUM
    for i in range(len(value_depth_list)):
        if dict2[value_depth_list[i][1]] > value_depth_list[i][0]:
            dict2[value_depth_list[i][1]] = value_depth_list[i][0]

    return(dict1, dict2)


def minmax(value_depth_list):

    dict1, dict2 = depth_minmax(value_depth_list)
    overall_max = max(list(dict1.values()))
    overall_min = min(list(dict2.values()))
    return(overall_max, overall_min)


def depth_average(value_depth_list):
    """averages values with same y coordinate - input in format [[value1,y1],[value2,y2],...]"""

    list2 = set([x[1] for x in value_depth_list])
    dict1 = {depth: 0 for depth in list2}
    dict2 = {depth: 0 for depth in list2}
    dict3 = {depth: 0 for depth in list2}

    for i in range(len(value_depth_list)):
        dict1[value_depth_list[i][1]] += value_depth_list[i][0]
        dict2[value_depth_list[i][1]] += 1

    for i in range(len(value_depth_list)):
        dict3[value_depth_list[i][1]] = dict1[value_depth_list[i][1]] / \
            dict2[value_depth_list[i][1]]

    return(dict3)


def mean_T_depth(T, directory):
    depth_temperature = []
    for ver in vertices(mesh):
        x = ver.point().x()
        y = ver.point().y()
        T_value = T(x, y)
        depth_temperature.append([T_value, ver.point().y()])
    mean_T_dict = depth_average(depth_temperature)
    f_T_depth = open(directory+"/T_depth.dat", "w")
    for i in range(len(list(mean_T_dict.keys()))):
        f_T_depth.write(str(list(mean_T_dict.keys())[
                        i])+"\t"+str(list(mean_T_dict.values())[i])+"\n")
    f_T_depth.flush()
    f_T_depth.close()


def mean_eta_depth(eta_fc, directory):
    depth_eta = []
    for ver in vertices(mesh):
        x = ver.point().x()
        y = ver.point().y()
        eta_value = eta_fc(x, y)
        depth_eta.append([eta_value, ver.point().y()])
    mean_eta_dict = depth_average(depth_eta)  # Depth average of eta
    eta_max, eta_min = minmax(depth_eta)  # Minimum and maximum over all domain
    f_eta_depth = open(directory+"/eta_depth.dat", "w")
    for i in range(len(list(mean_eta_dict.keys()))):
        f_eta_depth.write(str(list(mean_eta_dict.keys())[
                          i])+"\t"+str(list(mean_eta_dict.values())[i])+"\n")
    f_eta_depth.flush()
    f_eta_depth.close()
    return(eta_max, eta_min)


def mean_rms_depth(rms_fc, directory):
    depth_rms = []
    depth_v_x = []

    for ver in vertices(mesh):
        x = ver.point().x()
        y = ver.point().y()
        rms_value = rms_fc(x, y)
        v_x_value = v(x, y)[0]

        depth_rms.append([rms_value, ver.point().y()])
        depth_v_x.append([v_x_value, ver.point().y()])
        dict1, dict2 = depth_minmax(depth_v_x)
        # print(dict1)
        try:
            # max lateral velocity in top coordinate
            max_top_rms_x_vel = dict1[1.0]
        except:
            max_top_rms_x_vel = 0
    mean_rms_dict = depth_average(depth_rms)
    f_rms_depth = open(directory+"/rms_depth.dat", "w")
    for i in range(len(list(mean_rms_dict.keys()))):
        f_rms_depth.write(str(list(mean_rms_dict.keys())[
                          i])+"\t"+str(list(mean_rms_dict.values())[i])+"\n")
    f_rms_depth.flush()
    f_rms_depth.close()
    return(max_top_rms_x_vel)


def avg_rms_vel(v, hmax, lmax):
    return(sqrt(hmax/lmax*assemble(dot(v, v)*dx)))


def rms_vel(v, hmax, lmax):
    return(sqrt(hmax/lmax*(dot(v, v))))


def mean_top_rms_x(v, hmax, lmax):  # v^surf_RMS
    return(sqrt(hmax/lmax*assemble(v[0]*v[0]*ds(1))))


def tau(v):
    dimension = v.geometric_dimension()
    sigma = eta(T, x)*(grad(v)+grad(v).T)
    tau = sigma-p*Identity(dimension)
    return(tau)


def tau2(v):  # nepouzivam
    return (eta(T, x)*dev(2*sym(grad(v))))


def dynamic_topography(v):
    tauu = tau(v)  # +eta(T,x)*grad(v_)
    tensor_space = TensorFunctionSpace(mesh, 'CG', 1)
    stress = project(tauu, tensor_space)
    #stress = stress.vector().array()
    tau_22_up_left = stress(0.0, 1.0)[3]  # tau_22
    h_up_left = -alpha*deltaT/Ra*(tau_22_up_left-p(0.0, 1.0))
    print("TAU", h_up_left)
    h = -alpha*deltaT/Ra*dot(dot(n, tauu), n)
    mean_h = assemble(h*ds(1))
    result = h_up_left-mean_h
    return(result)


def open_files(var2, directory):
    # Initialize files with 0
    f_top_v_x, fNusselt, fRmsvel, file_h1, file_h2, file_h3, file_h4, file_h5, fNusseltTop, fNusseltBot, fMeanT, f_eta_minmax = (
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    f_top_v_x = open(directory+"/top_v_x.dat", "a+")
    fNusselt = open(directory+"/Nusselt.dat", "a+")
    fRmsvel = open(directory+"/Rmsvel.dat", "a+")
    if var2 == 1 or var2 == 2:
        file_h1 = open(directory+"/h1.dat", "a+")
    if var2 == 2:
        file_h2 = open(directory+"/h2.dat", "a+")
        file_h3 = open(directory+"/h3.dat", "a+")
        file_h4 = open(directory+"/h4.dat", "a+")
        file_h5 = open(directory+"/h5.dat", "a+")

    if var2 == 3:
        fNusseltTop = open(directory+"/nus_top.dat", "a+")
        fNusseltBot = open(directory+"/nus_bot.dat", "a+")
        fMeanT = open(directory+"/mean_T.dat", "a+")
        f_eta_minmax = open(directory+"/eta_minmax.dat", "a+")
    return(f_top_v_x, fNusselt, fRmsvel, file_h1, file_h2, file_h3, file_h4, file_h5, fNusseltTop, fNusseltBot, fMeanT, f_eta_minmax)


def close_files(var2, f_top_v_x, fNusselt, fRmsvel, file_h1, file_h2, file_h3, file_h4, file_h5, fNusseltTop, fNusseltBot, fMeanT, f_eta_minmax):
    f_top_v_x.close()
    fNusselt.close()
    fRmsvel.close()

    if var2 == 1 or var2 == 2:
        file_h1.close()

    if var2 == 2:
        file_h2.close()
        file_h3.close()
        file_h4.close()
        file_h5.close()
    if var2 == 3:
        fNusseltTop.close()
        fNusseltBot.close()
        fMeanT.close()
        f_eta_minmax.close()


def close_XDMF_files(var2, fT, fV, fP, file_rms, file_eta, file_h, file_eps_II, file_doteps):
    fT.close()
    fV.close()
    fP.close()
    file_rms.close()
    if var2 == 3 or var2 == 4 or var2 == 5:
        file_eta.close()
    if var2 == 2 or var2 == 5:
        file_h.close()
    if var2 == 4 or var2 == 5:
        file_eps_II.close()
        file_doteps.close()


write_t = 0

if var2 == 3:
    eta_fc.assign(project(eta(T, x, v), Pspace))

if var2 == 4 or var2 == 5:
    #both quantities are projected in logarithmic scale
    eta_fc.assign(project(ln(eta_eff(p, eps_II, v_k, x, t,True))/ln(10), EPSspace))
    doteps_II.assign(
        project(ln(second_invariant(grad(v)+grad(v).T))/ln(10), EPSspace))

try:
    while (t < t_end):
        f_top_v_x, fNusselt, fRmsvel, file_h1, file_h2, file_h3, file_h4, file_h5, fNusseltTop, fNusseltBot, fMeanT, f_eta_minmax = open_files(
            var2, directory)
        # Solve the problem
        if var2 == 4 or var2 == 5:
            if t_zero == False and harmonic:
                problem1 = LinearVariationalProblem(lhs(F1_harmonic), rhs(F1_harmonic), eps_II, bcEPS)
                solver1 = LinearVariationalSolver(problem1)
                solver1.parameters['linear_solver'] = "mumps"
                problem2_picard = LinearVariationalProblem(
                        lhs(F2_picard), rhs(F2_picard), w, bcs)  # ,DF2_picard)
                solver2_picard = LinearVariationalSolver(problem2_picard)
                solver2.parameters['linear_solver'] = "mumps"
            
            
            solver1.solve()
            tol = 0.0001
            iter = 0  # iteration counter
            """nastavit at dokonverguje"""
            maxiter = 50  # max no of iterations allowed
            eps = 1.0
            while eps > tol and iter < maxiter:
                # print(w_k,v_k,p_k,type(w_k),type(v_k))
                iter += 1
                v_k0 = w.split(deepcopy=True)[0]
                solver2_picard.solve()
                v_k1 = w.split(deepcopy=True)[0]
                #diff = v_k1.vector().get_local() - v_k0.vector().get_local()
                v_diff.assign(project(v_k1-v_k0,Vspace))
                eps = norm(v_diff.vector(),'linf')
                
                #eps = numpy.linalg.norm(diff, ord=numpy.Inf)
                info('iter=%d: norm=%g' % (iter, eps))
                w_k.assign(w)
            if write_t>50:
                write_data_interval = 5
        else:
            solver1.solve()
            if picard:
                """OUTDATED"""
                tol = 0.0001
                iter = 0  # iteration counter
                maxiter = 30  # max no of iterations allowed
                eps = 1.0
                while eps > tol and iter < maxiter:
                    iter += 1
                    v_k0 = w.split(deepcopy=True)[0]
                    if t_zero:
                        solver2_picard_tzero.solve()
                    else:
                        solver2_picard.solve()
                    v_k1 = w.split(deepcopy=True)[0]
                    diff = v_k1.vector().get_local() - v_k0.vector().get_local()
                    eps = numpy.linalg.norm(diff, ord=numpy.Inf)
                    info('iter=%d: norm=%g' % (iter, eps))
                    w_k.assign(w)
            else:
                if t_zero:
                    solver2_tzero.solve()
                else:
                    solver2.solve()

        # Extract solutions
        v, p = w.split()
        if var2 == 3:
            eta_fc.assign(project(eta(T, x, v), Pspace))

        if var2 == 4 or var2 == 5:
            eta_fc.assign(project(ln(eta_eff(p, eps_II, v_k, x, t,True))/ln(10), EPSspace))
            doteps_II.assign(
                project(ln(second_invariant(grad(v)+grad(v).T))/ln(10), EPSspace))
        rms_fc.assign(project(rms_vel(v, hmax, lmax), Pspace))
        if var2 == 2 or var2 == 5:
            # mesh movement
            h_k.assign(h)
            solver_h_pr.solve()
            dh.assign(project(h-h_k, HVspace))
            ALE.move(mesh, dh)

        T.rename("T", "temperature")
        v.rename("v", "velocity")
        p.rename("p", "pressure")
        # mean_T_depth_var = mean_T_depth(T)
        # print(mean_T_depth_var)
        if var2 == 2:
            h.rename("h", "h")
            if t == 0:
                for ver in vertices(mesh):
                    x1 = ver.point().x()
                    y1 = ver.point().y()
                    if abs(x1) < 0.01 and abs(y1-1.0) < 0.01:
                        ver1 = ver

                    if abs(x1-0.25) < 0.01 and abs(y1-1.0) < 0.01:
                        ver2 = ver

                    if abs(x1-0.5) < 0.01 and abs(y1-1.0) < 0.01:
                        ver3 = ver

                    if abs(x1-0.75) < 0.01 and abs(y1-1.0) < 0.01:
                        ver4 = ver

                    if abs(x1-1.00) < 0.01 and abs(y1-1.0) < 0.01:
                        ver5 = ver

            x1 = ver1.point().x()
            y1 = ver1.point().y()
            x2 = ver2.point().x()
            y2 = ver2.point().y()
            x3 = ver3.point().x()
            y3 = ver3.point().y()
            x4 = ver4.point().x()
            y4 = ver4.point().y()
            x5 = ver5.point().x()
            y5 = ver5.point().y()
            assert xdiscretization % 4 == 0  # assume we have discretization divisible by 4
            file_h1.write(str(t)+"\t"+str(y1-1)+"\n")
            file_h2.write(str(t)+"\t"+str(y2-1)+"\n")
            file_h3.write(str(t)+"\t"+str(y3-1)+"\n")
            file_h4.write(str(t)+"\t"+str(y4-1)+"\n")
            file_h5.write(str(t)+"\t"+str(y5-1)+"\n")
            file_h1.flush()
            file_h2.flush()
            file_h3.flush()
            file_h4.flush()
            file_h5.flush()
  
        # Move to next time step
        if var2 == 4 or var2 == 5:
            eps_II_prev.assign(eps_II)

        T0.assign(T)
        wvm = w.split(deepcopy=True)[0]
        dt_value = (Compute_dt(wvm))

        if var2 == 1:
            h_dyn = dynamic_topography(v)
            print(h_dyn)
            file_h1.write(str(t)+"\t"+str(h_dyn)+"\n")
            file_h1.flush()
        # Store solution to file
        # second_invariant((grad(v)+grad(v).T))
        if write_t % write_data_interval == 0:
            fT.write(T, t)
            fV.write(v, t)
            fP.write(p, t)
            if var2 == 2 or var2 == 5:
                file_h.write(h, t)
            if var2 == 4 or var2 == 5:
                file_eps_II.write(eps_II, t)
                file_doteps.write(doteps_II, t)

        # Report nusselt and root mean square velocity
        nus = nusselt(T, hmax, n)
        rms = avg_rms_vel(v, hmax, lmax)
        info('t = %g, nus = %g, rmsvel = %g' % (t, nus, rms))

        fNusselt.write((2*'%-15s ' + '\n') % (t, nus))
        fRmsvel.write((2*'%-15s ' + '\n') % (t, rms))
        fNusselt.flush()
        fRmsvel.flush()
        max_top_v_x = mean_rms_depth(rms_fc, directory)
        avg_top_rms_x = mean_top_rms_x(v, hmax, lmax)
        f_top_v_x.write(str(t)+"\t"+str(avg_top_rms_x) +
                        "\t"+str(max_top_v_x)+"\n")
        f_top_v_x.flush()

        if var2 == 4 or var2 == 5:

            eta_fc.rename("eta", "eta")
            if write_t % write_data_interval == 0:
                file_eta.write(eta_fc, t)

        if var2 == 3:
            nus_top = nusselt_top(T, n)
            nus_bot = nusselt_bot(T, n)
            info('t = %g, nus_top = %g, nus_bot = %g' % (t, nus_top, nus_bot))
            meant = mean_T(T)
#            max_eta = max_eta(T,x,v)
            fNusseltTop.write((2*'%-15s ' + '\n') % (t, nus_top))
            fNusseltBot.write((2*'%-15s ' + '\n') % (t, nus_bot))
            fMeanT.write((2*'%-15s ' + '\n') % (t, meant))
            fNusseltTop.flush()
            fNusseltBot.flush()
            fMeanT.flush()
#            fMaxEta.write((2*'%-15s ' + '\n') % (t,max_eta))
            mean_T_depth(T, directory)
            eta_max, eta_min = mean_eta_depth(eta_fc, directory)
            # eta_fc.assign(project(eta,Pspace))
            eta_fc.rename("eta", "eta")
            if write_t % write_data_interval == 0:
                file_eta.write(eta_fc, t)
            f_eta_minmax.write(str(t)+"\t"+str(eta_min)+"\t"+str(eta_max)+"\n")
            f_eta_minmax.flush()
        if write_t % write_data_interval == 0:
            file_rms.write(rms_fc, t)
        dt.assign(dt_value)

        if t_zero:
            t_zero = False
        t += dt_value
        write_t += 1
        close_files(var2, f_top_v_x, fNusselt, fRmsvel, file_h1, file_h2, file_h3,
                    file_h4, file_h5, fNusseltTop, fNusseltBot, fMeanT, f_eta_minmax)

    close_XDMF_files(var2, fT, fV, fP, file_rms, file_eta,
                     file_h, file_eps_II, file_doteps)

except KeyboardInterrupt:
    close_XDMF_files(var2, fT, fV, fP, file_rms, file_eta,
                     file_h, file_eps_II, file_doteps)
    close_files(var2, f_top_v_x, fNusselt, fRmsvel, file_h1, file_h2, file_h3,
                file_h4, file_h5, fNusseltTop, fNusseltBot, fMeanT, f_eta_minmax)