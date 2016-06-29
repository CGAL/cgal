//! \file examples/Arrangement_on_surface_2/spherical_degenerate_sweep.cpp
// Using the global aggregated insertion function.

//#define CGAL_SL_VERBOSE 0

//#define CGAL_ARRANGEMENT_ON_SURFACE_INSERT_VERBOSE 0

#include <vector>

#include <CGAL/config.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/Arr_geodesic_arc_on_sphere_traits_2.h>
#include <CGAL/Arr_spherical_topology_traits_2.h>

typedef CGAL::Exact_rational                               Number_type;
typedef CGAL::Cartesian<Number_type>                       Kernel;
typedef CGAL::Arr_geodesic_arc_on_sphere_traits_2<Kernel>  Geom_traits_2;
typedef Geom_traits_2::Point_2                             Point_2;
typedef Geom_traits_2::X_monotone_curve_2                  X_monotone_curve_2;
typedef CGAL::Arr_spherical_topology_traits_2<Geom_traits_2> Topol_traits_2;
typedef CGAL::Arrangement_on_surface_2<Geom_traits_2, Topol_traits_2>
                                                           Arrangement_2;
typedef Arrangement_2::Vertex_handle                       Vertex_handle;

int main()
{
  Geom_traits_2 traits;
  Geom_traits_2::Construct_point_2 ctr_p = traits.construct_point_2_object();
  Geom_traits_2::Construct_x_monotone_curve_2 ctr_xcv =
    traits.construct_x_monotone_curve_2_object();

  std::vector< Point_2 > points;
  std::vector< X_monotone_curve_2 > xcvs;

  CGAL::set_pretty_mode(std::cout);

  Point_2 sp = ctr_p(0, 0, -1);
  Point_2 np = ctr_p(0, 0, 1);
  points.push_back(sp);
  points.push_back(np);

  Point_2 p1 = ctr_p(-1,  0, 0);
  Point_2 p2 = ctr_p( 0, -1, 0);
  Point_2 p3 = ctr_p( 1,  0, 0);
  points.push_back(p1);
  points.push_back(p2);
  points.push_back(p3);

  X_monotone_curve_2 xcv_sp1 = ctr_xcv(sp, p1);
  X_monotone_curve_2 xcv_sp2 = ctr_xcv(sp, p2);
  X_monotone_curve_2 xcv_sp3 = ctr_xcv(sp, p3);
  CGAL_assertion(xcv_sp1.is_vertical());
  CGAL_assertion(xcv_sp2.is_vertical());
  CGAL_assertion(xcv_sp3.is_vertical());
  xcvs.push_back(xcv_sp1);
  xcvs.push_back(xcv_sp2);
  xcvs.push_back(xcv_sp3);

  X_monotone_curve_2 xcv_12 = ctr_xcv(p1, p2);
  X_monotone_curve_2 xcv_23 = ctr_xcv(p2, p3);
  CGAL_assertion(!xcv_12.is_vertical());
  CGAL_assertion(!xcv_23.is_vertical());
  xcvs.push_back(xcv_12);
  xcvs.push_back(xcv_23);

  X_monotone_curve_2 xcv_np1 = ctr_xcv(np, p1);
  X_monotone_curve_2 xcv_np2 = ctr_xcv(np, p2);
  X_monotone_curve_2 xcv_np3 = ctr_xcv(np, p3);
  CGAL_assertion(xcv_np1.is_vertical());
  CGAL_assertion(xcv_np2.is_vertical());
  CGAL_assertion(xcv_np3.is_vertical());
  xcvs.push_back(xcv_np1);
  xcvs.push_back(xcv_np2);
  xcvs.push_back(xcv_np3);

  unsigned subsetsp = (1 << points.size());
  std::cout << "#subsets points: " << subsetsp << std::endl;

  unsigned subsets = (1 << xcvs.size());
  std::cout << "#subsets curves: " << subsets << std::endl;

  std::cout << "total combinations: " << (subsetsp)*(subsets) << std::endl<< std::endl;;

  for (unsigned up = 0; up < subsetsp; up++) {
    std::vector< Point_2 > points_sub;
    for (unsigned ep = 0; ep <= points.size(); ep++) {
      if (up & (1 << ep)) {
        std::cout << "take pt: "  << ep << std::endl;
        points_sub.push_back(points[ep]);
      }
    }

    for (unsigned u = 0; u < subsets; u++) {
      std::vector< X_monotone_curve_2 > xcvs_sub;
      for (unsigned e = 0; e <= xcvs.size(); e++) {
        if (u & (1 << e)) {
          std::cout << "take xcv: "  << e << std::endl;
          xcvs_sub.push_back(xcvs[e]);
        }
      }

      std::cout << "subsetpoints #" << up << " has size: "  << points_sub.size() << std::endl;
      std::cout << "subsetcurves #" << u << " has size: "  << xcvs_sub.size() << std::endl;

#if 1
      Arrangement_2 arr;
      std::cout << "inserting "
                << xcvs_sub.size() << " x-monotone curves and "
                << points_sub.size() << " isolated points."
                << std::endl;

      // TODO why is this signature not available as "insert(...)"
      CGAL::insert_empty(arr, xcvs_sub.begin(), xcvs_sub.end(), points_sub.begin(), points_sub.end());

      // Print the size of the arrangement.
      std::cout << "The arrangement size:" << std::endl
                << "   V = " << arr.number_of_vertices()
                << ",  E = " << arr.number_of_edges()
                << ",  F = " << arr.number_of_faces() << std::endl;

      std::cout << "=======================================================" << std::endl << std::endl << std::endl;
#endif
      //std::cout << "arr: " << arr << std::endl;
      std::cout << std::endl;
    }
  }
  return 0;
}
