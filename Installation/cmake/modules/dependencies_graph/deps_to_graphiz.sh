#!/bin/bash

if [ -n "$1" ]; then
    set > $1
fi

printf "strict digraph CGAL_Dependencies {\n  rank=source\n  rankdir=LR\n"
echo <<EOF
subgraph Foundations {
label=Foundation
node [shape=box]
Algebraic_foundations
Algebraic_kernel_d
Arithmetic_kernel
BGL
Cartesian_kernel
CGAL_Core
CGAL_ImageIO
Circulator
Combinatorial_map
Distance_2
Distance_3
Filtered_kernel
Generator
Geomview
HalfedgeDS
Hash_map
Homogeneous_kernel
Installation
Intersections_2
Intersections_3
Interval_support
Inventor
Kernel_23
Kernel_d
LEDA
Modifier
Modular_arithmetic
NewKernel_d
Number_types
OpenNL
Optimisation_basic
Polygon
Polynomial
Profiling_tools
Property_map
Random_numbers
Solver_interface
Spatial_sorting
STL_Extension
Stream_support
Subdivision_method_3
Testsuite
Union_find    
} 
EOF
for f in package_info/*/dependencies; do
    pkg=${f#package_info/}
    pkg=${pkg%/dependencies}
    printf "  %s;\n" $pkg
    for dep in $(cat $f); do
        printf "  %s -> %s;\n" $pkg $dep
    done
done
printf "}\n"
