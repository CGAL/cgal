#include <list>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Arr_spherical_gaussian_map_3/Arr_polyhedral_sgm_traits.h>
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
using Direction_3 = Kernel::Direction_3;

#if 0
using Gm_traits = CGAL::Arr_polyhedral_sgm_traits<Kernel, -8, 6>;
#elif 0
using Gm_traits = CGAL::Arr_polyhedral_sgm_traits<Kernel, -11, 7>;
#else
using Gm_traits = CGAL::Arr_polyhedral_sgm_traits<Kernel, -1, 0>;
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

using Point_location_result = CGAL::Arr_point_location_result<Gm>;
using Query_result = std::pair<Point_2, Point_location_result::Type>;

using Vertex_const_handle = Gm::Vertex_const_handle;
using Halfedge_const_handle = Gm::Halfedge_const_handle;
using Face_const_handle = Gm::Face_const_handle;

int main() {
  Gm_polyhedron p;
  p.make_tetrahedron(Point_3(1.0, 0.0, 0.0), Point_3(0.0, 1.0, 0.0),
                     Point_3(0.0, 0.0, 1.0), Point_3(0.0, 0.0, 0.0));
  Gm gm;

  Naive_pl naive_pl(gm);
  // Landmarks_pl landmarks_pl(gm);
  Walk_pl walk_pl(gm);
  // Trap_pl trap_pl(gm);
  /* Need to add the code below to both Arr_spherical_gaussian_map_3 and
   * Arr_polyhedral_sgm, and then work on the trap point location code...

private:
  friend class Arr_observer<Self>;
  friend class Arr_accessor<Self>;

protected:
  typedef Arr_observer<Self>                      Observer;

  void _register_observer(Observer *p_obs)
  { Base::_register_observer((typename Base::Observer*)p_obs); }

  bool _unregister_observer(Observer *p_obs)
  { return (Base::_unregister_observer ((typename Base::Observer*)p_obs)); }
  */

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

  // locate_point(trap_pl, points[0]);

  ////////
  std::list<Query_result> results;
  // The following cause an assertion failure.
  // CGAL::locate(gm, &points[0], &points[3], std::back_inserter(results));

  // Print the results.
  for (auto it = results.begin(); it != results.end(); ++it) {
    std::cout << "The point (" << it->first << ") is located ";
    if (const Face_const_handle* f =
        std::get_if<Face_const_handle>(&(it->second)))       // inside a face
      std::cout << "inside "
                << (((*f)->is_unbounded()) ? "the unbounded" : "a bounded")
                << " face.\n";
    else if (const Halfedge_const_handle* e =
             std::get_if<Halfedge_const_handle>(&(it->second))) // on an edge
      std::cout << "on an edge: " << (*e)->curve() << std::endl;
    else if (const Vertex_const_handle* v =
             std::get_if<Vertex_const_handle>(&(it->second)))  // on a vertex
      std::cout << "on "
                << (((*v)->is_isolated()) ? "an isolated" : "a")
                << " vertex: " << (*v)->point() << std::endl;
  }

  return 0;
}
