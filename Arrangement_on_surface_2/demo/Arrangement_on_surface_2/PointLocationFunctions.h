#ifndef ARRANGEMENT_DEMO_POINT_LOCATION_FUNCTIONS
#define ARRANGEMENT_DEMO_POINT_LOCATION_FUNCTIONS

#include <CGAL/Object.h>

class QPointF;

template <typename Arr_>
struct PointLocationFunctions
{
  using Arrangement = Arr_;
  using Traits = typename Arrangement::Geometry_traits_2;
  using Halfedge_const_handle = typename Arrangement::Halfedge_const_handle;
  using Face_const_handle = typename Arrangement::Face_const_handle;
  using Vertex_const_handle = typename Arrangement::Vertex_const_handle;
  using Halfedge_around_vertex_const_circulator =
    typename Arrangement::Halfedge_around_vertex_const_circulator;

  CGAL::Object locate(const Arrangement*, const QPointF&);
  Face_const_handle getFace(const Arrangement*, const QPointF&);
  CGAL::Object rayShootUp(const Arrangement*, const QPointF&);
  CGAL::Object rayShootDown(const Arrangement*, const QPointF&);
};

#endif
