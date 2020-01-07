# octahedron-rotation-for-vasp

the python file is used to rotate the octahedron to align the M-O bond with inner coordinate, i.e., x, y, z axis.

# Usage

`python octahedron_rotation.py POSCAR`

# Notice
By default, here the `POSCAR.vasp` and `POSCAR.xyz` are read in for processing.
Only cartesian coordinate is considred, but fractional one can be easily transfred by `VESTA`
The input argument is filename with no extension
The example here shows how to rotate the `CoO6` octahedron in the LiCoO2 material.