#ifndef TEXTURED_POLYHEDRON_TYPE_H 
#define TEXTURED_POLYHEDRON_TYPE_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/textured_polyhedron.h>

#ifdef USE_FORWARD_DECL
struct EPIC_Kernel : public CGAL::Exact_predicates_inexact_constructions_kernel {};
#else
typedef CGAL::Exact_predicates_inexact_constructions_kernel EPIC_Kernel;
#endif
typedef CGAL::Textured_polyhedron<EPIC_Kernel,CGAL::Textured_items> Textured_polyhedron;

#endif // TEXTURED_POLYHEDRON_TYPE_H
