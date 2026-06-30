# Flow Past a Cylinder Benchmark

This repository contains a 2-D benchmark for simulating incompressible flow past a circular cylinder in a rectangular channel. The benchmark is designed to generate flow fields for different cylinder radii and Reynolds numbers, making it useful for studying laminar and periodic wake dynamics and for producing datasets for scientific machine learning or generative modeling experiments.

<table align="center">
  <tr>
    <td align="center" width="24%"><img src="gif_images/Re_20_r_0.01.gif" width="100%" alt="Flow past cylinder, Re=20, r=0.01" /><br><sub>Re = 20, r = 0.01</sub></td>
    <td align="center" width="24%"><img src="gif_images/Re_20_r_0.05.gif" width="100%" alt="Flow past cylinder, Re=20, r=0.05" /><br><sub>Re = 20, r = 0.05</sub></td>
    <td align="center" width="24%"><img src="gif_images/Re_120_r_0.05.gif" width="100%" alt="Flow past cylinder, Re=120, r=0.05" /><br><sub>Re = 120, r = 0.05</sub></td>
    <td align="center" width="24%"><img src="gif_images/Re_120_r_0.1.gif" width="100%" alt="Flow past cylinder, Re=120, r=0.1" /><br><sub>Re = 120, r = 0.10</sub></td>
  </tr>
</table>

## Overview
The code solves a canonical flow-past-a-cylinder problem in a 2-D channel with a circular obstacle. The workflow consists of two main steps:

1. **Mesh generation** using Gmsh for a rectangular channel with a circular cylinder.
2. **Flow simulation** using FEniCS/DOLFIN to solve the time-dependent Navier-Stokes equations and save velocity-magnitude snapshots.

The benchmark supports parameter sweeps over:

- Cylinder radius, `r`
- Reynolds number, `Re`

The generated outputs include mesh files, velocity time-series data, log files, and image frames that can be converted into animations.

## Repository structure

```text
Flow_past_a_cylinder_benchmark/
├── 2D-mesh-generator.py
├── flow_simulation_generate_images.py
├── script_mesh_generator.sh
├── script_simulation.sh
├── mesh-visualizations.ipynb
├── fenicsproject.yml
├── gif_images/
└── README.md
```
## Problem setup

The computational domain is a 2-D rectangular channel with a circular cylinder removed from the interior.

Default geometry in `2D-mesh-generator.py`:

| Quantity | Value |
|---|---:|
| Channel length, `L` | `2.2` |
| Channel height, `H` | `0.41` |
| Cylinder center | `(0.2, 0.2)` |
| Cylinder radius | User-specified, `r` |

The mesh is refined near the cylinder and downstream wake region. Physical boundary groups are assigned for the fluid, inlet, outlet, walls, and obstacle.

## Installation

Clone the repository:

```bash
git clone https://github.com/ruchi-sandilya/Flow_past_a_cylinder_benchmark.git
cd Flow_past_a_cylinder_benchmark
```

Create and activate the conda environment:

```bash
conda env create -f fenicsproject.yml
conda activate fenicsproject
```

The code uses Python with Gmsh, MeshIO, FEniCS/DOLFIN, NumPy, Matplotlib, and HDF5 support.

## Usage

### 1. Generate a mesh for one cylinder radius

```bash
python 2D-mesh-generator.py 0.05
```

This creates:

```text
generated_mesh/mesh2D_r_0.05.msh
```

### 2. Run one flow simulation

```bash
python flow_simulation_generate_images.py 0.05 20
```

Arguments:

```text
python flow_simulation_generate_images.py <radius> <Reynolds_number>
```

Example:

```bash
python flow_simulation_generate_images.py 0.05 120
```

This reads the generated mesh, converts it to XDMF format, solves the flow problem, stores velocity time-series data, and saves velocity-magnitude images.

### 3. Generate meshes over a radius sweep

```bash
bash script_mesh_generator.sh
```

The default sweep generates meshes for:

```text
r = 0.01, 0.012, 0.014, ..., 0.10
```

Mesh-generation logs are saved in:

```text
txt_output_files/
```

### 4. Run simulation sweeps

```bash
bash script_simulation.sh
```

The script includes example cases for laminar and periodic-flow regimes:

| Case | Radius range | Reynolds number |
|---|---:|---:|
| Case 1 | `0.01` to `0.048` | `20` |
| Case 2 | `0.062` to `0.10` | `120` |
| Case 3 | `0.05` | `20` to `40` |
| Case 4 | `0.05` | `100` to `120` |

Simulation logs are written to:

```text
txt_output_files/
```

## Outputs

The main generated folders are:

```text
generated_mesh/          # .msh and .xdmf mesh files
navier_stokes_cylinder/  # FEniCS velocity time-series data
Images/                  # Velocity-magnitude image snapshots
txt_output_files/        # Log files from batch runs
```

Velocity-magnitude frames are saved with names of the form:

```text
Images/velmag_r_<radius>_Re_<Reynolds_number>_step_<step>.png
```

## Notes

- Run mesh generation before running the flow simulation for a given radius.
- The simulation script assumes the corresponding mesh file already exists in `generated_mesh/`.
- Generated image and velocity files can become large. Consider adding generated output folders to `.gitignore` if you do not want to track them in Git.
- FEniCS/DOLFIN installation can be system-dependent; using the provided conda environment is recommended.

## Citation

If you use this benchmark or the generated data in your work, please cite:

```bibtex
@article{sandilya2026conda,
  title={Contrastive Diffusion Alignment: Learning Structured Latents for Controllable Generation},
  author={Sandilya, R. and Perez, S. and Lynch, C. and Victoria, L. and Zebley, B. and Buchanan, D. M. and Bhati, M. T. and Williams, N. and Spellman, T. and Gunning, F. M. and Liston, C. and Grosenick, L.},
  journal={Forty-Third International Conference on Machine Learning (ICML). arXiv:2510.14190},
  year={2026}
}
```
