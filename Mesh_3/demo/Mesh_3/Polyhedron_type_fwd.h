#ifndef POLYHEDRON_TYPE_FWD_H
#define POLYHEDRON_TYPE_FWD_H

#include <memory>  // for std::allocator
#include <utility> // for std::pair

#ifdef USE_FORWARD_DECL

namespace CGAL {

  template < typename FT_ >
  struct Simple_cartesian;

  class Epick;

  class Polyhedron_items_3;

  template < class T, class I, class A>
  class HalfedgeDS_default;

  template < class PolyhedronTraits_3,
             class PolyhedronItems_3,
             template < class T, class I, class A>
             class T_HDS, 
             class Alloc
             >
  class Polyhedron_3;

  namespace Mesh_3 {
    template <typename Patch_id>
    class Mesh_polyhedron_items;
  } // end namespace Mesh_3
} // end namespace CGAL

// kernel

typedef CGAL::Epick Kernel;

typedef std::pair<int, int> Patch_id;
// surface mesh
typedef CGAL::Polyhedron_3<Kernel,
                           CGAL::Mesh_3::Mesh_polyhedron_items<Patch_id>,
                           CGAL::HalfedgeDS_default,
                           std::allocator<int> > Polyhedron;

#else // USE_FORWARD_DECL

#include "Polyhedron_type.h"

#endif // USE_FORWARD_DECL

#endif // POLYHEDRON_TYPE_FWD_H
