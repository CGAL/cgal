#ifndef NEF_TYPE_H
#define NEF_TYPE_H

// CGAL
// kernel
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

// surface mesh
#include <CGAL/Polyhedron_3.h>

// nef
#include <CGAL/Nef_polyhedron_3.h> 
#include <CGAL/Nef_3/SNC_indexed_items.h>

// Boolean operations work only with exact kernel
#ifdef USE_FORWARD_DECL
// struct Exact_Kernel : public CGAL::Exact_predicates_exact_constructions_kernel {};
struct Exact_Kernel : public CGAL::Simple_cartesian<CGAL::Gmpq> {};
#else
typedef CGAL::Exact_predicates_exact_constructions_kernel Exact_Kernel;
#endif
typedef CGAL::Polyhedron_3<Exact_Kernel> Exact_polyhedron;

typedef CGAL::Nef_polyhedron_3<Exact_Kernel,
			       CGAL::SNC_indexed_items,
			       bool> Nef_polyhedron; 

#endif // NEF_TYPE_H
