#ifndef POLYHEDRON_TYPE_FWD_H
#define POLYHEDRON_TYPE_FWD_H


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <memory>


typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
//struct Kernel;
namespace CGAL {
  class Polyhedron_items_3;

#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
  template < class T, class I, class A>
#endif
  class HalfedgeDS_default;

  template < class PolyhedronTraits_3,
             class PolyhedronItems_3,
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
             template < class T, class I, class A>
#endif
             class T_HDS, 
             class Alloc
             >
  class Polyhedron_3;
}

typedef CGAL::Polyhedron_3<Kernel,
                           CGAL::Polyhedron_items_3,
                           CGAL::HalfedgeDS_default,
                           std::allocator<int> > Polyhedron;

#endif // POLYHEDRON_TYPE_FWD_H
