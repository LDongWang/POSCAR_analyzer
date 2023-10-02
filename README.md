# POSCAR_analyzer
A tool reading POSCAR file and outputting crystal information.

Function:
1. Compute lattice information
2. Compute distance matrice for each atoms' pair


Usage: python POSCAR2info [-h] [--output OUTPUT] [--symmetry SYMMETRY] [--dist] [--info] filename
positional arguments:
  filename             input POSCAR name

options:
  -h, --help           show this help message and exit
  --output OUTPUT      output filename
  --symmetry SYMMETRY  specify symmetry point group
  --dist               output distance matrice
  --info               output lattice information
