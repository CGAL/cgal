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

template <class Triangle_3, class FaceGraph, class VertexPointMap>
Triangle_3 triangle_from_halfedge(typename boost::graph_traits<FaceGraph>::halfedge_descriptor edge, const FaceGraph& faceGraph, VertexPointMap vertexPointMap)
{
  typedef typename boost::graph_traits<FaceGraph>::halfedge_descriptor halfedge_descriptor;
  
  halfedge_descriptor e0 = edge;
  halfedge_descriptor e1 = next(edge, faceGraph);

  return Triangle_3(get(vertexPointMap, boost::source(e0, faceGraph)), get(vertexPointMap, boost::target(e0, faceGraph)), get(vertexPointMap, boost::target(e1, faceGraph)));
}

template <class Triangle_3, class FaceGraph>
Triangle_3 triangle_from_halfedge(typename boost::graph_traits<FaceGraph>::halfedge_descriptor edge, const FaceGraph& faceGraph)
{
  return triangle_from_halfedge<Triangle_3, FaceGraph, typename boost::property_map<FaceGraph, boost::vertex_point_t>::type>(edge, faceGraph, get(boost::vertex_point, faceGraph));
}


template <class FaceGraph>
size_t edge_index(typename boost::graph_traits<FaceGraph>::halfedge_descriptor he, FaceGraph& p)
{
  typedef typename boost::graph_traits<FaceGraph> GraphTraits;
  typedef typename GraphTraits::face_descriptor face_descriptor;
  typedef typename GraphTraits::halfedge_descriptor halfedge_descriptor;
  
  face_descriptor f = face(he, p);
  
  halfedge_descriptor start = halfedge(f, p);
  halfedge_descriptor current = start;
  
  size_t count = 0;
  
  while (current != he)
  {
    current = next(current, p);
    ++count;
  }

  return count;
}

template <class FT>
FT internal_sqrt(const FT& x, CGAL::Tag_true)
{
  return CGAL::sqrt(x);
}

template <class FT>
FT internal_sqrt(const FT& x, CGAL::Tag_false)
{
  return FT(std::sqrt(CGAL::to_double(x)));
}

template <class FT>
FT select_sqrt(const FT& x)
{
  typedef ::CGAL::Algebraic_structure_traits<FT> AST; 
  static const bool has_sqrt = ! ::boost::is_same< ::CGAL::Null_functor, typename AST::Sqrt >::value;
  return internal_sqrt(x, ::CGAL::Boolean_tag<has_sqrt>());
}

} // namespace internal

} // namespace CGAL

#endif // CGAL_POLYHEDRON_SHORTEST_PATH_INTERNAL_MISC_H