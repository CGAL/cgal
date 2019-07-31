#!/bin/sed -f

s/build-check-deps/BUILD_DIR/g;

#
# LGPL packages
#
s/Algebraic_foundations/Foundation_packages__LGPL_/g;
s/Algebraic_kernel_d/Foundation_packages__LGPL_/g;
s/Arithmetic_kernel/Foundation_packages__LGPL_/g;
s/BGL/Foundation_packages__LGPL_/g;
s/Cartesian_kernel/Foundation_packages__LGPL_/g;
s/CGAL_Core/Foundation_packages__LGPL_/g;
s/CGAL_ImageIO/Foundation_packages__LGPL_/g;
s/Circulator/Foundation_packages__LGPL_/g;
s/Combinatorial_map/Foundation_packages__LGPL_/g;
s/Conic_2/Foundation_packages__LGPL_/g;
s/Distance_2/Foundation_packages__LGPL_/g;
s/Distance_3/Foundation_packages__LGPL_/g;
s/Filtered_kernel/Foundation_packages__LGPL_/g;
s/Generator/Foundation_packages__LGPL_/g;
s/Geomview/Foundation_packages__LGPL_/g;
s/HalfedgeDS/Foundation_packages__LGPL_/g;
s/Hash_map/Foundation_packages__LGPL_/g;
s/Homogeneous_kernel/Foundation_packages__LGPL_/g;
s/Installation/Foundation_packages__LGPL_/g;
s/Intersections_2/Foundation_packages__LGPL_/g;
s/Intersections_3/Foundation_packages__LGPL_/g;
s/Interval_support/Foundation_packages__LGPL_/g;
s/Inventor/Foundation_packages__LGPL_/g;
s/Kernel_23/Foundation_packages__LGPL_/g;
s/\bKernel_d\b/Foundation_packages__LGPL_/g;
s/LEDA/Foundation_packages__LGPL_/g;
s/Modifier/Foundation_packages__LGPL_/g;
s/Modular_arithmetic/Foundation_packages__LGPL_/g;
s/NewKernel_d/Foundation_packages__LGPL_/g;
s/Number_types/Foundation_packages__LGPL_/g;
s/OpenNL/Foundation_packages__LGPL_/g;
s/Optimisation_basic/Foundation_packages__LGPL_/g;
s/\bPolygon\b/Foundation_packages__LGPL_/g;
s/Polynomial/Foundation_packages__LGPL_/g;
s/Profiling_tools/Foundation_packages__LGPL_/g;
s/Property_map/Foundation_packages__LGPL_/g;
s/Random_numbers/Foundation_packages__LGPL_/g;
s/Solver_interface/Foundation_packages__LGPL_/g;
s/Spatial_sorting/Foundation_packages__LGPL_/g;
s/STL_Extension/Foundation_packages__LGPL_/g;
s/Stream_support/Foundation_packages__LGPL_/g;
s/Subdivision_method_3/Foundation_packages__LGPL_/g;
s/Testsuite/Foundation_packages__LGPL_/g;
s/Union_find/Foundation_packages__LGPL_/g;

#
# Special LGPL packages
#
s/CGAL_ipelets/Ipelets__LGPL__/g;
s/Linear_cell_complex/Linear_cell_complex__LGPL__/g;
s/Generalized_map/Generalized_map__LGPL__/g;
s/Kinetic_data_structures/Kinetic__LGPL__/g;
s/Kinetic_framework/Kinetic__LGPL__/g;

# renaming
s/Principal_component_analysis_LGPL/PCA__LGPL_files__/g;
