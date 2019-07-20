#define CGAL_IDENTIFICATION_XY 2

// Testing the spherical gaussian map
#include <iostream>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Arr_spherical_gaussian_map_3/Arr_polyhedral_sgm_traits.h>
#include <CGAL/Arr_spherical_gaussian_map_3/Arr_polyhedral_sgm.h>
#include <CGAL/Arr_spherical_gaussian_map_3/Arr_polyhedral_sgm_polyhedron_3.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel       Kernel;
typedef Kernel::Point_3                                         Point_3;
typedef CGAL::Arr_polyhedral_sgm_traits<Kernel>                 Gm_traits;
typedef CGAL::Arr_polyhedral_sgm<Gm_traits>                     Gm;
typedef CGAL::Arr_polyhedral_sgm_polyhedron_3<Gm, Kernel>       Gm_polyhedron;
typedef CGAL::Arr_polyhedral_sgm_initializer<Gm, Gm_polyhedron> Gm_initializer;

int main()
{
  // Construct the Gaussian map of a tetrahedron
  Point_3 points[] = {
    Point_3(1.0, 0.0, 0.0),
    Point_3(0.0, 1.0, 0.0),
    Point_3(0.0, 0.0, 1.0),
    Point_3(0.0, 0.0, 0.0)
  };
  Gm_polyhedron P1;
  P1.make_tetrahedron(points[0], points[1], points[2], points[3]);
  Gm gm1;
  Gm_initializer gm_initializer1(gm1);
  gm_initializer1(P1);
  if (! gm1.is_valid()) return -1;

  // Construct the Gaussian map of the reflection of a tetrahedron
  Gm_polyhedron P2;
  for (Point_3* p = points; p != &points[4]; ++p) {
    Kernel::Vector_3 v = CGAL::ORIGIN - *p;
    *p = CGAL::ORIGIN + v;
  }
  P2.make_tetrahedron(points[1], points[0], points[2], points[3]);
  std::cout << P2 << std::endl;
  Gm gm2;
  Gm_initializer gm_initializer2(gm2);
  gm_initializer2(P2);
  if (! gm2.is_valid()) return -1;

  // Compute the Minowski sum of the Gaussian maps
  Gm gm;
  gm.minkowski_sum(gm1, gm2);
  // std::cout << gm << std::endl;
  if (! gm.is_valid()) return -1;

  Kernel::FT sw(16);
  Gm::Vertex_const_handle it;
  for (it = gm.vertices_begin(); it != gm.vertices_end(); ++it) {
    if (it->degree() < 3) continue;
    Gm::Halfedge_around_vertex_const_circulator hec3(it->incident_halfedges());
    Gm::Halfedge_around_vertex_const_circulator hec1 = hec3++;
    Gm::Halfedge_around_vertex_const_circulator hec2 = hec3++;
    std::cout << (*hec1).face()->point() << ", "
              << (*hec2).face()->point() << ", "
              << (*hec3).face()->point() << std::endl;
    Kernel::Plane_3 plane((*hec1).face()->point(), (*hec2).face()->point(),
                          (*hec3).face()->point());
    Kernel::Vector_3 v(CGAL::ORIGIN, plane.projection(CGAL::ORIGIN));
    Kernel::FT tmp = v.squared_length();
    if (tmp < sw) sw = tmp;
  }
  // std::cout << sw << std::endl;
  if ((3 * sw) != 1) return -1;
  return 0;
}
