#ifndef TEXTURED_POLYHEDRON_TYPE_H
#define TEXTURED_POLYHEDRON_TYPE_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/textured_polyhedron.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Textured_polyhedron<Kernel,CGAL::Textured_items> Textured_polyhedron;

#endif // TEXTURED_POLYHEDRON_TYPE_H
