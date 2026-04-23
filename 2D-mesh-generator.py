import warnings
import os
warnings.filterwarnings("ignore")
import gmsh
gmsh.initialize()
import numpy as np
import sys

gmsh.model.add("DFG 2D")
L = 2.2
H = 0.41
c_x = c_y =0.2
r_ = 0.05
r = float(sys.argv[1]) # change for different radius
gdim = 2


# ---------------------------
# Create geometry
# ---------------------------
channel_tag = gmsh.model.occ.addRectangle(0.0, 0.0, 0.0, L, H, tag=1)
cylinder_tag = gmsh.model.occ.addDisk(c_x, c_y, 0.0, r, r)

# Subtract cylinder from channel
cut_result = gmsh.model.occ.cut([(gdim, channel_tag)], [(gdim, cylinder_tag)])
gmsh.model.occ.synchronize()

# Get resulting fluid surface
fluid_surfaces = gmsh.model.getEntities(dim=gdim)
assert len(fluid_surfaces) == 1, f"Expected 1 fluid surface, found {len(fluid_surfaces)}"
fluid_dim, fluid_tag = fluid_surfaces[0]

# ---------------------------
# Physical groups
# ---------------------------
fluid_marker = 1
inlet_marker = 2
outlet_marker = 3
wall_marker = 4
obstacle_marker = 5

gmsh.model.addPhysicalGroup(fluid_dim, [fluid_tag], fluid_marker)
gmsh.model.setPhysicalName(fluid_dim, fluid_marker, "Fluid")

inflow = []
outflow = []
walls = []
obstacle_edges = []

boundaries = gmsh.model.getBoundary([(fluid_dim, fluid_tag)], oriented=False)

tol = 1e-6
for dim, tag in boundaries:
    com = gmsh.model.occ.getCenterOfMass(dim, tag)

    if np.allclose(com, [0.0, H / 2.0, 0.0], atol=tol):
        inflow.append(tag)
    elif np.allclose(com, [L, H / 2.0, 0.0], atol=tol):
        outflow.append(tag)
    elif np.allclose(com, [L / 2.0, 0.0, 0.0], atol=tol) or np.allclose(com, [L / 2.0, H, 0.0], atol=tol):
        walls.append(tag)
    else:
        obstacle_edges.append(tag)

gmsh.model.addPhysicalGroup(1, inflow, inlet_marker)
gmsh.model.setPhysicalName(1, inlet_marker, "Inlet")

gmsh.model.addPhysicalGroup(1, outflow, outlet_marker)
gmsh.model.setPhysicalName(1, outlet_marker, "Outlet")

gmsh.model.addPhysicalGroup(1, walls, wall_marker)
gmsh.model.setPhysicalName(1, wall_marker, "Walls")

gmsh.model.addPhysicalGroup(1, obstacle_edges, obstacle_marker)
gmsh.model.setPhysicalName(1, obstacle_marker, "Obstacle")

# ---------------------------
# Mesh sizing controls
# ---------------------------
# Let the background fields control the mesh size
gmsh.option.setNumber("Mesh.CharacteristicLengthFromPoints", 0)
gmsh.option.setNumber("Mesh.CharacteristicLengthFromCurvature", 0)
gmsh.option.setNumber("Mesh.CharacteristicLengthExtendFromBoundary", 0)

# Global mesh size scale
lc_far = 0.02          # coarse size away from cylinder/wake
lc_cyl = r / 8.0       # finer around cylinder
lc_wake = r / 6.0      # refinement in wake

# 1) Distance-based refinement around cylinder
dist_field = gmsh.model.mesh.field.add("Distance")
gmsh.model.mesh.field.setNumbers(dist_field, "EdgesList", obstacle_edges)

thresh_field = gmsh.model.mesh.field.add("Threshold")
gmsh.model.mesh.field.setNumber(thresh_field, "IField", dist_field)
gmsh.model.mesh.field.setNumber(thresh_field, "LcMin", lc_cyl)
gmsh.model.mesh.field.setNumber(thresh_field, "LcMax", lc_far)
gmsh.model.mesh.field.setNumber(thresh_field, "DistMin", 0.5 * r)
gmsh.model.mesh.field.setNumber(thresh_field, "DistMax", 4.0 * r)

# 2) Wake refinement box downstream of cylinder
wake_field = gmsh.model.mesh.field.add("Box")
gmsh.model.mesh.field.setNumber(wake_field, "VIn", lc_wake)
gmsh.model.mesh.field.setNumber(wake_field, "VOut", lc_far)

# Wake region: from just behind cylinder toward outlet
gmsh.model.mesh.field.setNumber(wake_field, "XMin", c_x + 0.5 * r)
gmsh.model.mesh.field.setNumber(wake_field, "XMax", min(L, c_x + 12.0 * r))
gmsh.model.mesh.field.setNumber(wake_field, "YMin", max(0.0, c_y - 3.0 * r))
gmsh.model.mesh.field.setNumber(wake_field, "YMax", min(H, c_y + 3.0 * r))
gmsh.model.mesh.field.setNumber(wake_field, "Thickness", 0.03)

# Take the minimum of refinement fields
min_field = gmsh.model.mesh.field.add("Min")
gmsh.model.mesh.field.setNumbers(min_field, "FieldsList", [thresh_field, wake_field])
gmsh.model.mesh.field.setAsBackgroundMesh(min_field)

# ---------------------------
# Mesh generation options
# ---------------------------
# Keep triangles for easier FEniCS compatibility
gmsh.option.setNumber("Mesh.Algorithm", 6)   # Frontal-Delaunay for 2D
gmsh.model.mesh.generate(gdim)
gmsh.model.mesh.optimize("Netgen")


# create folder if it doesn't exist
os.makedirs("generated_mesh", exist_ok=True)

gmsh.write("generated_mesh/mesh2D_r_{}.msh".format(r))
