#ifndef TEXTURED_POLYHEDRON_TYPE_H 
#define TEXTURED_POLYHEDRON_TYPE_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/textured_polyhedron.h>

#include "Textured_polyhedron_type_fwd.h"

#ifdef USE_FORWARD_DECL
struct EPIC_kernel : public CGAL::Exact_predicates_inexact_constructions_kernel {};
#else
typedef CGAL::Exact_predicates_inexact_constructions_kernel EPIC_kernel;
typedef CGAL::Textured_polyhedron<EPIC_kernel,CGAL::Textured_items> Textured_polyhedron;
#endif

#endif // TEXTURED_POLYHEDRON_TYPE_H
