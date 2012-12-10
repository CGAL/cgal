#ifndef POLYHEDRON_TYPE_FWD_H
#define POLYHEDRON_TYPE_FWD_H

#include <CGAL/Filtered_kernel_fwd.h>
#include <memory>

#ifdef USE_FORWARD_DECL

#include <CGAL/Filtered_kernel_fwd.h>

template <typename Patch_id>
class Polyhedron_demo_items;

namespace CGAL {

  namespace Mesh_3 {
    template <typename Kernel>
    struct Robust_intersection_traits_3;
  }

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
  
// changed since CGAL-3.8-Ic-8
#if CGAL_VERSION_NR > 1030800008 
  class Epick;
#endif
} // end namespace CGAL

// kernel
namespace polyhedron_type_fwd_h {
// changed since CGAL-3.8-Ic-8
#if CGAL_VERSION_NR > 1030800008 
  typedef CGAL::Epick K1;
#else
  typedef CGAL::Filtered_kernel< CGAL::Simple_cartesian<double>, true > K1;
#endif
  typedef CGAL::Mesh_3::Robust_intersection_traits_3<K1> Kernel;
}

#else // USE_FORWARD_DECL

#include "Polyhedron_type.h"

#endif // USE_FORWARD_DECL

// surface mesh
typedef CGAL::Polyhedron_3<polyhedron_type_fwd_h::Kernel,
                           Polyhedron_demo_items<int>,
                           // CGAL::Polyhedron_items_3,
                           CGAL::HalfedgeDS_default,
                           std::allocator<int> > Polyhedron;

#endif // POLYHEDRON_TYPE_FWD_H
