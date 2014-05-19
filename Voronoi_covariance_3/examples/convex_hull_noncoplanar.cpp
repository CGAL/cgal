#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/point_generators_3.h>
#include <CGAL/algorithm.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Voronoi_covariance_3/Convex_hull_traits_dual_3.h>

#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedron_3<K>                               Polyhedron_3;
typedef CGAL::Voronoi_covariance_3::Convex_hull_traits_dual_3<K>
                                                            Hull_traits_dual_3;
typedef CGAL::Polyhedron_3<Hull_traits_dual_3>              Polyhedron_dual_3;

int main ()
{
  // define polyhedron to hold convex hull
  Polyhedron_dual_3 dual_poly;
  
  // traits
  Hull_traits_dual_3 dual_traits;
  
  // dual convex hull
  std::cout << "Dual convex hull" << std::endl;
  std::list<K::Plane_3> planes;
  planes.push_back(K::Plane_3(1, 0, 0, -1));
  planes.push_back(K::Plane_3(1, 0, 0, 1));
  planes.push_back(K::Plane_3(0, 1, 0, -1));
  planes.push_back(K::Plane_3(0, 1, 0, 1));
  planes.push_back(K::Plane_3(0, 0, 1, -1));
  planes.push_back(K::Plane_3(0, 0, 1, 1));
  
  typedef CGAL::Triangulation_data_structure_2<
    CGAL::Triangulation_vertex_base_with_info_2<int, CGAL::GT3_for_CH3<Hull_traits_dual_3> >,
    CGAL::Convex_hull_face_base_2<int, Hull_traits_dual_3> >                           Tds;  
  Tds dual_tds;
  CGAL::internal::Convex_hull_3::non_coplanar_quickhull_3(planes, dual_tds, dual_traits);
  return 0;
}

