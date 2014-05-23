#ifndef CGAL_TEST_UTIL_H
#define CGAL_TEST_UTIL_H

#include <algorithm>

namespace CGAL {

namespace test {

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
void Construct_polyhedron_planes(Polyhedron& out)
{
  std::transform( out.facets_begin(), out.facets_end(), out.planes_begin(), Plane_from_facet<Polyhedron>());
}

template <class Polyhedron>
void Make_regular_tetrahedron(Polyhedron& out)
{
  typedef typename Polyhedron::Kernel::FT FT;
  
  FT rsqrt2 = FT(1.0) / CGAL::sqrt(FT(2.0));
  out.clear();
  out.make_tetrahedron(
    Polyhedron::Point_3(FT(1.0), FT(0.0), -rsqrt2),
    Polyhedron::Point_3(-FT(1.0), FT(0.0), -rsqrt2),
    Polyhedron::Point_3(FT(0.0), FT(1.0), rsqrt2),
    Polyhedron::Point_3(FT(0.0), FT(1.0), rsqrt2));
  Construct_polyhedron_planes(out);
}

} // namespace util

} // namespace CGAL

#endif // CGAL_TEST_UTIL_H