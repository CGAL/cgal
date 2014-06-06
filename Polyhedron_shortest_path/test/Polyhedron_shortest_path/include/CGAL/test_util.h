#ifndef CGAL_TEST_UTIL_H
#define CGAL_TEST_UTIL_H

#include <algorithm>

namespace CGAL {

namespace test {

#include <CGAL/boost/graph/properties.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>

template<class Polyhedron>
struct Plane_from_facet {
  typedef typename Polyhedron::Plane_3 Plane_3;
  typedef typename Polyhedron::Facet Facet;
  typedef typename Polyhedron::Halfedge_handle Halfedge_handle;

  Plane_3 operator()(Facet& f) {
      Halfedge_handle h = f.halfedge();
      return Plane_3( h->vertex()->point(),
                                    h->next()->vertex()->point(),
                                    h->opposite()->vertex()->point());
  }
};

template <class Polyhedron>
void construct_polyhedron_planes(Polyhedron& out)
{
  std::transform( out.facets_begin(), out.facets_end(), out.planes_begin(), Plane_from_facet<Polyhedron>());
}

template <class Polyhedron>
typename Polyhedron::Halfedge_handle make_regular_tetrahedron(Polyhedron& out)
{
  typedef typename Polyhedron::Traits::FT FT;
  
  FT rsqrt2 = FT(1.0) / CGAL::sqrt(FT(2.0));
  out.clear();
  typename Polyhedron::Halfedge_handle result = out.make_tetrahedron(
    typename Polyhedron::Point_3(FT(1.0), FT(0.0), -rsqrt2),
    typename Polyhedron::Point_3(-FT(1.0), FT(0.0), -rsqrt2),
    typename Polyhedron::Point_3(FT(0.0), FT(1.0), rsqrt2),
    typename Polyhedron::Point_3(FT(0.0), -FT(1.0), rsqrt2));
  construct_polyhedron_planes(out);
  return result;
}

template <class Polyhedron>
size_t face_vertex_index(typename boost::graph_traits<Polyhedron>::face_descriptor face, typename boost::graph_traits<Polyhedron>::vertex_descriptor vertex, Polyhedron& P)
{
  size_t index = 0;
  
  typedef typename boost::graph_traits<Polyhedron>::halfedge_descriptor halfedge_descriptor;
  
  halfedge_descriptor currentEdge(CGAL::halfedge(face, P));
  halfedge_descriptor startEdge = currentEdge;
  
  do
  {
    if (CGAL::source(currentEdge, P) == vertex)
    {
      return index;
    }
    
    ++index;
    currentEdge = CGAL::next(currentEdge, P);
  }
  while (currentEdge != startEdge);
  
  return index;
}

} // namespace util

} // namespace CGAL

#endif // CGAL_TEST_UTIL_H