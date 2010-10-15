#ifndef C2T3_TYPE_H
#define C2T3_TYPE_H

#include "Polyhedron_type.h"

#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>


typedef CGAL::Surface_mesh_default_triangulation_3 Tr;
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;

#endif
