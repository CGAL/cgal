// (LicenseStuffHere)
//
// $URL$
// $Id$
// 
//
// Author(s)     : Stephen Kiazyk

#ifndef CGAL_POLYHEDRON_SHORTEST_PATH_INTERNAL_MISC_H
#define CGAL_POLYHEDRON_SHORTEST_PATH_INTERNAL_MISC_H

#include <CGAL/boost/graph/properties.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/iterator.h>

#include <CGAL/Interval_nt.h>

namespace CGAL {

namespace internal {

template <class Triangle_3, class Polyhedron, class VertexPointMap>
Triangle_3 triangle_from_halfedge(typename boost::graph_traits<Polyhedron>::halfedge_descriptor edge, const Polyhedron& polyhedron, VertexPointMap vertexPointMap)
{
  typedef typename boost::graph_traits<Polyhedron>::halfedge_descriptor halfedge_descriptor;
  
  halfedge_descriptor e0 = edge;
  halfedge_descriptor e1 = CGAL::next(edge, polyhedron);

  return Triangle_3(vertexPointMap[boost::source(e0, polyhedron)], vertexPointMap[boost::target(e0, polyhedron)], vertexPointMap[boost::target(e1, polyhedron)]);
}

template <class Triangle_3, class Polyhedron>
Triangle_3 triangle_from_halfedge(typename boost::graph_traits<Polyhedron>::halfedge_descriptor edge, const Polyhedron& polyhedron)
{
  return triangle_from_halfedge<Triangle_3, Polyhedron, typename boost::property_map<Polyhedron, CGAL::vertex_point_t>::type>(edge, polyhedron, CGAL::get(CGAL::vertex_point, polyhedron));
}


template <class Polyhedron>
size_t edge_index(typename boost::graph_traits<Polyhedron>::halfedge_descriptor he, Polyhedron& p)
{
  typedef typename boost::graph_traits<Polyhedron> GraphTraits;
  typedef typename GraphTraits::face_descriptor face_descriptor;
  typedef typename GraphTraits::halfedge_descriptor halfedge_descriptor;
  
  face_descriptor f = CGAL::face(he, p);
  
  halfedge_descriptor start = CGAL::halfedge(f, p);
  halfedge_descriptor current = start;
  
  size_t count = 0;
  
  while (current != he)
  {
    current = CGAL::next(current, p);
    ++count;
  }

  return count;
}

template <class P, class FT>
P interpolate_points(const P& p0, const P& p1, FT t)
{
  FT t0 = FT(1.0) - t;
  
  return P(CGAL::ORIGIN) + (((p0 - P(CGAL::ORIGIN)) * t0) + ((p1 - P(CGAL::ORIGIN)) * t));
}

template <class V>
V shift_vector_3_left(const V& v, size_t by)
{
  return V(v[by], v[(by + 1) % 3], v[(by + 2) % 3]);
}

template <class V>
V shift_vector_3_right(const V& v, size_t by)
{
  by %= 3;
  return V(v[(3 - by) % 3], v[(4 - by) % 3], v[(5 - by) % 3]);
}

/*
enum Line_relation
{
  LINE_RELATION_UNKNOWN = 0,

  LINE_RELATION_INTERSECT,
  LINE_RELATION_NO_INTERSECT,
  LINE_RELATION_PARALLEL,
};

template <class FT>
struct Intersection_result
{
  Intersection_result()
  {
    relation = LINE_RELATION_UNKNOWN;
    t0 = 0.0;
    t1 = 0.0;
  }

  Intersection_result(Line_relation _relation, FT _t0, FT _t1)
  {
    relation = _relation;
    t0 = _t0;
    t1 = _t1;
  }

  Line_relation relation;
  FT t0;
  FT t1;
};

template <class FT, class V>
FT cross_product(const V& v0, const V& v1)
{
  return v0.x()*v1.y() - v0.y()*v1.x();
}

template <class FT, class R, class V>
Intersection_result<FT> intersect_rays(const R& r0, const R& r1)
{
  V r0D = r0.to_vector();
  V r1D = r1.to_vector();

  FT crossProd = cross_product<FT,V>(r0D, r1D);

  V diff = r1.start() - r0.start();

  FT r0Cross = cross_product<FT,V>(diff, r0D);
  FT r1Cross = cross_product<FT,V>(diff, r1D);

  if (CGAL::abs(crossProd) < 1e-10)
  {
    if (CGAL::abs(r0Cross) < 1e-10 || CGAL::abs(r1Cross) < 1e-10)
    {
      return Intersection_result<FT>(LINE_RELATION_PARALLEL, 0.0, 0.0);
    }
    else
    {
      return Intersection_result<FT>(LINE_RELATION_NO_INTERSECT, 0.0, 0.0);
    }
  }

  // Yes, this is correct
  FT t1 = r0Cross / crossProd;
  FT t0 = r1Cross / crossProd;

  return Intersection_result<FT>(LINE_RELATION_INTERSECT, t0, t1);
}
*/

} // namespace internal

} // namespace CGAL

#endif // CGAL_POLYHEDRON_SHORTEST_PATH_INTERNAL_MISC_H