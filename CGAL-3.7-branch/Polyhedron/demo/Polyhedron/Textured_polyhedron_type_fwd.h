#ifndef TEXTURED_POLYHEDRON_TYPE_FWD_H 
#define TEXTURED_POLYHEDRON_TYPE_FWD_H

struct EPIC_kernel;
namespace CGAL {
  struct Textured_items;
  template <class Kernel, class Items> class Textured_polyhedron;
}

typedef CGAL::Textured_polyhedron< ::EPIC_kernel,CGAL::Textured_items> Textured_polyhedron;

#endif // TEXTURED_POLYHEDRON_TYPE_FWD_H
