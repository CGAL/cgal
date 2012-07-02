// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#ifndef _cmdline_h
#define _cmdline_h

//#include <tchar.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <vector>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Periodic_implicit_mesh_domain_3.h>
#include <CGAL/make_periodic_mesh_3.h>
#include <CGAL/Mesh_3_periodic_triangulation_3.h>

#include <CGAL/Mesh_constant_domain_field_3.h>

#include <CGAL/Bbox_3.h>
#include <CGAL/bounding_box.h>

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT FT;
typedef K::Point_3 Point;
typedef K::Point_3 Point_3;
typedef FT (Function)(const Point&);
typedef CGAL::Periodic_implicit_mesh_domain_3<Function,K> Periodic_mesh_domain;

// Triangulation
typedef CGAL::Mesh_3_periodic_triangulation_3_generator<Periodic_mesh_domain>::type Mesh_3_periodic_triangulation_3;
typedef Mesh_3_periodic_triangulation_3 Tr;

typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

typedef Tr::Geom_traits GT;
typedef GT::Sphere_3 Sphere_3;
typedef GT::Vector_3 Vector_3;
typedef GT::FT FT;

typedef Periodic_mesh_domain::Iso_cuboid_3 Iso_cuboid_3;
typedef Periodic_mesh_domain::Bbox_3 Bbox_3;

/*
// TODO: reference additional headers your program requires here

#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>

#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h> 
#include <fstream>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/output_surface_facets_to_polyhedron.h>

#include <CGAL/Bbox_3.h>
#include <CGAL/bounding_box.h>

typedef CGAL::Bbox_3 Bbox_3;

// default triangulation for Surface_mesher
typedef CGAL::Surface_mesh_default_triangulation_3 Tr;

// c2t3
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;

typedef Tr::Geom_traits GT;
typedef GT::Sphere_3 Sphere_3;
typedef GT::Point_3 Point_3;
typedef GT::Vector_3 Vector_3;
typedef GT::FT FT;

typedef FT (*Function)(Point_3);

typedef CGAL::Implicit_surface_3<GT, Function> Surface_3;

typedef GT::Iso_cuboid_3 Iso_cuboid_3;

typedef CGAL::Polyhedron_3<GT>         Polyhedron;
*/

typedef enum {
	righthanded,
	lefthanded
} t_chirality;

typedef struct t_crystalsize {
	t_crystalsize() { x=0; y=0; z=0; }
	t_crystalsize(int dx, int dy, int dz) { x=dx; y=dy; z=dz; }
	int x;
	int y;
	int z;
} t_crystalsize;


#endif
