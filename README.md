# clamato_voids
Voids identified in the following manner:

(for the map)

load_and_smooth_map.py: reads in the map and smooths it with 2 h^-1 Mpc Gaussian kernel

make_underdensity_file.py: identifies all regions in the smoothed map with delta_F > 0.224.
Grows a sphere around each until the average delta_F is 0.167.

find_voids.py: removes overlapping voids, keeping the largest void in the case of a conflict.
writes voids to an output file.

(for the Nyx simulation)

make_underdensity_file_nyx_mpi.py: identify all regions in the redshift-space density field 
with Delta < 0.15
and grows spheres out to Delta < 0.3.  Code is very slow so I run it on Edison using mpi4py
identify_void_regions.py: removes overlapping voids, keeping the largest void in the case of a conflict.
Removes all voids with R < 2 h^-1 Mpc.  Writes an output file with void regions identified.


cross correlation files: in progress!
