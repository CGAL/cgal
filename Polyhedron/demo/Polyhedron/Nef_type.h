#ifndef NEF_TYPE_H
#define NEF_TYPE_H

// CGAL
// kernel
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

// surface mesh
#include <CGAL/Polyhedron_3.h>

// nef
#include <CGAL/Nef_polyhedron_3.h> 

// Boolean operations work only with exact kernel
typedef CGAL::Exact_predicates_exact_constructions_kernel Exact_Kernel;
typedef CGAL::Polyhedron_3<Exact_Kernel> Exact_polyhedron;

typedef CGAL::Nef_polyhedron_3<Exact_Kernel> Nef_polyhedron; 

#endif // NEF_TYPE_H
