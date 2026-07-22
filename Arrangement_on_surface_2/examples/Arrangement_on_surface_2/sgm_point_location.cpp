#include <list>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Arr_spherical_gaussian_map_3/Arr_polyhedral_sgm.h>
#include <CGAL/Arr_spherical_gaussian_map_3/Arr_polyhedral_sgm_polyhedron_3.h>

#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/Arr_landmarks_point_location.h>
#include <CGAL/Arr_walk_along_line_point_location.h>
#include <CGAL/Arr_trapezoid_ric_point_location.h>
#include <CGAL/Arr_batched_point_location.h>

#include "point_location_utils.h"

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Point_3 = Kernel::Point_3;

#if 0
using Gm_traits = CGAL::Arr_geodesic_arc_on_sphere_traits_2<Kernel, -8, 6>;
#elif 0
using Gm_traits = CGAL::Arr_geodesic_arc_on_sphere_traits_2<Kernel, -11, 7>;
#else
using Gm_traits = CGAL::Arr_geodesic_arc_on_sphere_traits_2<Kernel, -1, 0>;
#endif

using Gm = CGAL::Arr_polyhedral_sgm<Gm_traits>;
using Gm_polyhedron = CGAL::Arr_polyhedral_sgm_polyhedron_3<Gm, Kernel>;
using Gm_initializer = CGAL::Arr_polyhedral_sgm_initializer<Gm, Gm_polyhedron>;

using Naive_pl = CGAL::Arr_naive_point_location<Gm>;
using Walk_pl = CGAL::Arr_walk_along_line_point_location<Gm>;
using Landmarks_pl = CGAL::Arr_landmarks_point_location<Gm>;
using Trap_pl = CGAL::Arr_trapezoid_ric_point_location<Gm>;

using Geom_traits = Gm::Geometry_traits_2;
using Point_2 = Geom_traits::Point_2;

int main() {
  Gm_polyhedron p;
  p.make_tetrahedron(Point_3(1.0, 0.0, 0.0), Point_3(0.0, 1.0, 0.0),
                     Point_3(0.0, 0.0, 1.0), Point_3(0.0, 0.0, 0.0));
  Gm gm;

  Naive_pl naive_pl(gm);
  Landmarks_pl landmarks_pl(gm);
  Walk_pl walk_pl(gm);
  // Trap_pl trap_pl(gm);

  Gm_traits traits;
  Gm_initializer gm_initializer(gm);
  gm_initializer(p);
  if (! gm.is_valid()) return -1;

  auto ctr_point = traits.construct_point_2_object();
  Point_2 points[] = {
    ctr_point(-1, 0, 0),
    ctr_point(0, -1, 0),
    ctr_point(0, 0, -1)
  };

  locate_point(naive_pl, points[0]);
  locate_point(naive_pl, points[1]);
  locate_point(naive_pl, points[2]);

  // locate_point(walk_pl, points[0]);
  // locate_point(walk_pl, points[1]);
  // locate_point(walk_pl, points[2]);

  locate_point(landmarks_pl, points[0]);
  locate_point(landmarks_pl, points[1]);
  locate_point(landmarks_pl, points[2]);

  // locate_point(trap_pl, points[0]);
  // locate_point(trap_pl, points[1]);
  // locate_point(trap_pl, points[2]);

  return 0;
}
