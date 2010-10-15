#ifndef TEXTURED_POLYHEDRON_TYPE_FWD_H 
#define TEXTURED_POLYHEDRON_TYPE_FWD_H

#ifdef USE_FORWARD_DECL

struct EPIC_kernel;
namespace CGAL {
  struct Textured_items;
  template <class Kernel, class Items> class Textured_polyhedron;
}

typedef CGAL::Textured_polyhedron< ::EPIC_kernel,CGAL::Textured_items> Textured_polyhedron;

#else // USE_FORWARD_DECL
#  include "Textured_polyhedron_type.h"
#endif // USE_FORWARD_DECL

#endif // TEXTURED_POLYHEDRON_TYPE_FWD_H
