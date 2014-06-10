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

namespace CGAL {

namespace internal {

template <class GraphTraits>
size_t edge_index(typename GraphTraits::halfedge_descriptor he, GraphTraits& p)
{
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

template <class P, class V, class FT>
P interpolate_points(const P& p0, const P& p1, FT t)
{
  FT t0 = FT(1.0) - t;
  
  return P(CGAL::ORIGIN) + (((p0 - P(CGAL::ORIGIN)) * t0) + ((p1 - P(CGAL::ORIGIN)) * t));
}

template <class V>
V shift_vector_3(const V& v, size_t by)
{
  return V(v[by], v[(by + 1) % 3], v[(by + 2) % 3]);
}

} // namespace internal

} // namespace CGAL

#endif // CGAL_POLYHEDRON_SHORTEST_PATH_INTERNAL_MISC_H