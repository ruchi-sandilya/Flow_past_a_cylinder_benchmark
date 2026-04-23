import numpy as np
import meshio
from fenics import *
from dolfin import *
import matplotlib.pyplot as plt
import dolfin.common.plotting as fenicsplot 
parameters = " matplotlib "
import sys
import os
import h5py


r = float(sys.argv[1]) # radius of obstacle (cylinder)
Re = int(sys.argv[2]) # Reynolds Number

print(Re, r)

msh = meshio.read("generated_mesh/mesh2D_r_{}.msh".format(r))

tri_cells = []
for cell in msh.cells:
    if  cell.type == "triangle":
        if len(tri_cells) == 0:
            tri_cells = cell.data
        else:
            tri_cells = np.vstack([tri_cells, cell.data])

tri_mesh = meshio.Mesh(points=msh.points[:,:2], cells={"triangle": tri_cells})
meshio.write("generated_mesh/mesh2D_r_{}.xdmf".format(r), tri_mesh)
mesh = Mesh()
XDMFFile("generated_mesh/mesh2D_r_{}.xdmf".format(r)).read(mesh)

T = 1.5            # final time
num_steps = 3000   # number of time steps
dt = T / num_steps # time step size
#dt = 0.0005
mu = 0.001         # dynamic viscosity
rho = 1            # density

# Define function spaces
V = VectorFunctionSpace(mesh, 'P', 2)
Q = FunctionSpace(mesh, 'P', 1)


# Define boundaries
inflow   = 'near(x[0], 0)'
outflow  = 'near(x[0], 2.2)'
walls    = 'near(x[1], 0) || near(x[1], 0.41)'
cylinder = 'on_boundary && x[0]>0.1 && x[0]<0.3 && x[1]>0.1 && x[1]<0.3'

# Define inflow profile
_U = Constant(0.003*Re/(4*r)) 
inflow_profile = ('4.0*_U*x[1]*(0.41 - x[1]) / pow(0.41, 2)', '0')

# Define boundary conditions
bcu_inflow = DirichletBC(V, Expression(inflow_profile, degree=2, _U =_U), inflow)
bcu_walls = DirichletBC(V, Constant((0, 0)), walls)
bcu_cylinder = DirichletBC(V, Constant((0, 0)), cylinder)
bcp_outflow = DirichletBC(Q, Constant(0), outflow)
bcu = [bcu_inflow, bcu_walls, bcu_cylinder]
bcp = [bcp_outflow]

# Define trial and test functions
u = TrialFunction(V)
v = TestFunction(V)
p = TrialFunction(Q)
q = TestFunction(Q)

# Define functions for solutions at previous and current time steps
u_n = Function(V)
u_  = Function(V)
p_n = Function(Q)
p_  = Function(Q)

# Define expressions used in variational forms
U  = 0.5*(u_n + u)
n  = FacetNormal(mesh)
f  = Constant((0, 0))
k  = Constant(dt)
mu = Constant(mu)
rho = Constant(rho)

# Define symmetric gradient
def epsilon(u):
    return sym(nabla_grad(u))

# Define stress tensor
def sigma(u, p):
    return 2*mu*epsilon(u) - p*Identity(len(u))

# Define variational problem for step 1
F1 = rho*dot((u - u_n) / k, v)*dx + rho*dot(dot(u_n, nabla_grad(u_n)), v)*dx + inner(sigma(U, p_n), epsilon(v))*dx + dot(p_n*n, v)*ds - dot(mu*nabla_grad(U)*n, v)*ds - dot(f, v)*dx
a1 = lhs(F1)
L1 = rhs(F1)

# Define variational problem for step 2
a2 = dot(nabla_grad(p), nabla_grad(q))*dx
L2 = dot(nabla_grad(p_n), nabla_grad(q))*dx - (1/k)*div(u_)*q*dx

# Define variational problem for step 3
a3 = dot(u, v)*dx
L3 = dot(u_, v)*dx - k*dot(nabla_grad(p_ - p_n), v)*dx

# Assemble matrices
A1 = assemble(a1)
A2 = assemble(a2)
A3 = assemble(a3)

# Apply boundary conditions to matrices
[bc.apply(A1) for bc in bcu]
[bc.apply(A2) for bc in bcp]

# create folder if it doesn't exist
os.makedirs("navier_stokes_cylinder", exist_ok=True)


# Create velocity time series for saving data (velocity and pressure in .h5 format)
timeseries_u = TimeSeries('navier_stokes_cylinder/velocity_series_r_{}_Re_{}'.format(r,Re))

# Create progress bar
progress = Progress('Time-stepping')
set_log_level(LogLevel.PROGRESS)

# Time-stepping
t = 0

for n in range(num_steps):
    print('steps', n)
    # Update current time
    t += dt

    # Step 1: Tentative velocity step
    b1 = assemble(L1)
    [bc.apply(b1) for bc in bcu]
    solve(A1, u_.vector(), b1, 'bicgstab', 'hypre_amg')

    # Step 2: Pressure correction step
    b2 = assemble(L2)
    [bc.apply(b2) for bc in bcp]
    solve(A2, p_.vector(), b2, 'bicgstab', 'hypre_amg')

    # Step 3: Velocity correction step
    b3 = assemble(L3)
    solve(A3, u_.vector(), b3, 'cg', 'sor')
    
    #File('navier_stokes_cylinder/u_{}.pvd'.format(n+1)) << (u_, t)
    #File('navier_stokes_cylinder/p_{}.pvd'.format(n+1)) << (p_, t)
    
    # Save nodal values to file
    if (n+1) % 30 == 0:
        timeseries_u.store(u_.vector(), t)
        # timeseries_p.store(p_.vector(), t)

    # Update previous solution
    u_n.assign(u_)
    p_n.assign(p_)

    # Update progress bar
    #progress.update(t / T)
    progress = t/T
    
    set_log_level(LogLevel.ERROR)
    if (n+1)%100 == 0:
        print ("[Step %d/%d] [u max: %.4f] [p max: %.4f] [Time: %.4f]" % 
               (n+1, num_steps, u_.vector().get_local().max(), p_.vector().get_local().max(),t))
    
    # Save images for each timesteps    
    #plt.figure(figsize=(5,3))
    plot(sqrt(pow(u_[0],2)+pow(u_[1],2)), cmap = 'inferno')
    plt.xticks([])
    plt.yticks([])
    plt.xlim((0,1.5))
    plt.tight_layout()
    os.makedirs("Images", exist_ok=True)
    plt.savefig('Images/velmag_r_{}_Re_{}_step_{}.png'.format(r,Re,n))
    plt.close()


# v_f = 'navier_stokes_cylinder/velocity_series_r_{}_Re_{}.h5'.format(r,Re)

# v_df = h5py.File(v_f, 'r')

# arr_v = []
# for v_dset in v_df['Vector'].keys():
#     arr_v.append(v_df['Vector'][v_dset][:])  
# arr_v = np.asarray(arr_v)   

# for n in range(num_steps):
#     v.vector()[:] = arr_v[n,:]
#     print("[Step %d] [Re = : %d] [r =: %.4f]" % (n+1, Re, r))    
#     # Plot images
#     #plt.figure(figsize=(6,2))
#     v_m = sqrt(pow(v[0],2)+pow(v[1],2))
#     plot(v_m, cmap = 'inferno')
#     plt.xticks([])
#     plt.yticks([])
#     plt.xlim((0,1.5))
#     plt.tight_layout()
#     plt.savefig('Images/velmag_r_{}_Re_{}_step_{}.png'.format(r,Re,n))
#     plt.close()
