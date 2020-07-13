#include <cassert>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_traits_3.h>

int main() {

  using Kernel  = CGAL::Simple_cartesian<double>;
  using Point_3 = typename Kernel::Point_3;
  using Plane_3 = typename Kernel::Plane_3;
  using Traits  = CGAL::Polyhedron_traits_3<Kernel>;

  using Polyhedron = CGAL::Polyhedron_3<Traits>;
  using HDS = typename Polyhedron::HDS;
  using Vertex = typename Polyhedron::Vertex;

  using Vertex_iterator = typename Polyhedron::Vertex_iterator;
  using Halfedge_iterator = typename Polyhedron::Halfedge_iterator;
  using Facet_iterator = typename Polyhedron::Facet_iterator;

  using Vertex_handle = typename Polyhedron::Vertex_handle;
  using Halfedge_handle = typename Polyhedron::Halfedge_handle;
  using Facet_handle = typename Polyhedron::Facet_handle;

  using Point_iterator = typename Polyhedron::Point_iterator;
  using Plane_iterator = typename Polyhedron::Plane_iterator;

  Polyhedron polyhedron;
  auto h = polyhedron.make_tetrahedron(
    Point_3(1, 0, 0), Point_3(0, 1, 0),
    Point_3(0, 0, 1), Point_3(0, 0, 0));
  assert(polyhedron.is_valid());
  assert(polyhedron.is_tetrahedron(h));

  {
    Point_iterator begin(polyhedron.points_begin());
    Point_iterator end(polyhedron.points_end());
    Vertex_iterator vit = polyhedron.vertices_begin();
    for (auto const& vh : polyhedron.vertex_handles()){
      assert(vh->point() == *begin);
      ++vit; ++begin;
    }
    assert(vit == polyhedron.vertices_end());
  }

  {
    Plane_iterator begin(polyhedron.planes_begin());
    Plane_iterator end(polyhedron.planes_end());
    Facet_iterator fit = polyhedron.facets_begin();
    for (auto const& fh : polyhedron.facet_handles()){
      assert(fh->plane() == *begin);
      ++fit; ++begin;
    }
    assert(fit == polyhedron.facets_end());
  }

  {
    Halfedge_iterator hit = polyhedron.halfedges_begin();
    for (auto const& h : polyhedron.halfedge_handles()) {
      ++hit;
    }
    assert(hit == polyhedron.halfedges_end());
  }

  return EXIT_SUCCESS;
}
