#ifndef CGAL_TOPOLOGICAL_MAP_ITEMS_H
#define CGAL_TOPOLOGICAL_MAP_ITEMS_H 1

#ifndef CGAL_HALFEDGEDS_VERTEX_BASE_H
#include <CGAL/HalfedgeDS_vertex_base.h>
#endif
#ifndef CGAL_HALFEDGEDS_HALFEDGE_BASE_H
#include <CGAL/HalfedgeDS_halfedge_base.h>
#endif
#ifndef CGAL_HALFEDGEDS_FACE_BASE_H
#include <CGAL/Topological_map_face_base.h>
#endif

CGAL_BEGIN_NAMESPACE

class Topological_map_items {
public:
  template < class Refs, class Traits >
  struct Vertex_wrapper {
    typedef typename Traits::Point_2 Point;
    typedef HalfedgeDS_vertex_base< Refs, Tag_true, Point > Vertex;
  };
  template < class Refs, class Traits>
  struct Halfedge_wrapper {
    typedef typename Traits::X_curve X_curve;
    typedef HalfedgeDS_halfedge_base< Refs >                Halfedge;
  };
  template < class Refs, class Triats>
  struct Face_wrapper {
    typedef Topological_map_face_list_base< Refs >          Face;
  };
};

CGAL_END_NAMESPACE
#endif // CGAL_TOPOLOGICAL_MAP_ITEMS_H
// EOF //
