#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Delaunay_triangulation_on_sphere_traits_2.h>
#include <CGAL/Delaunay_triangulation_on_sphere_2.h>
#include <CGAL/Projection_on_sphere_traits_3.h>

#include <CGAL/algorithm.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/Timer.h>

#include <boost/iterator/transform_iterator.hpp>

#include <cmath>
#include <fstream>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel    K;
typedef CGAL::Surface_mesh<K>                                  Surface_mesh;

typedef K::Segment_3                                           Segment_3;
typedef CGAL::Delaunay_triangulation_3<K>                      Delaunay;

typedef CGAL::Delaunay_triangulation_on_sphere_traits_2<K>     Gt;
typedef CGAL::Projection_on_sphere_traits_3<K>                 Gt2;
typedef CGAL::Delaunay_triangulation_on_sphere_2<Gt>           DTOS;
typedef CGAL::Delaunay_triangulation_on_sphere_2<Gt2>          DTOS2;
typedef K::Point_3                                             Point;

typedef CGAL::Delaunay_triangulation_3<K, CGAL::Fast_location> Delaunay_fast;
typedef CGAL::Creator_uniform_3<double, Point>                 Creator;

int main(int, char**)
{
  CGAL::Timer time;

  const std::size_t nu_of_pts = 1e7;
  const double radius = 5184.152;

  CGAL::Random random;
  std::cout << "Seed is " << random.get_seed() << std::endl;

  CGAL::Random_points_on_sphere_3<Point, Creator> on_sphere(radius, random);

  std::vector<Point> points;
  points.reserve(nu_of_pts);

  for(std::size_t count=0; count<nu_of_pts; ++count)
    points.push_back(*on_sphere++);
  std::cout << points.size() << " points" << std::endl;

  // Delaunay_traits
  DTOS dtos;
  dtos.set_radius(radius);

  std::cout << " ***STARTING***" << std::endl;
  time.start();
  dtos.insert(points.begin(), points.end());
  time.stop();
  assert(dtos.number_of_vertices() == nu_of_pts);
  std::cout << "Triangulation sphere: "
            << dtos.number_of_vertices() << " vertices in " << time.time() << " sec" << std::endl;

  //Triangulation with points on the sphere (projection_traits)
  Gt2 traits(K::Point_3(0, 0, 0), radius);
  DTOS2 dtos2(traits);
  Gt2::Construct_point_on_sphere_2 cst = traits.construct_point_on_sphere_2_object();

  time.reset();
  time.start();
  dtos2.insert(boost::make_transform_iterator(points.begin(), cst),
               boost::make_transform_iterator(points.end(), cst));
  time.stop();
  std::cout << "Triangulation w/ sphere projection traits: "
            << dtos2.number_of_vertices() << " vertices in " << time.time() << " sec" << std::endl;

//  Surface_mesh sm;

//  time.reset();
//  time.start();
//  CGAL::convex_hull_3(points.begin(), points.end(), sm);
//  time.stop();
//  std::cout << "Convex hull 3D: " << time.time() << " " << std::endl;

  time.reset();
  time.start();
  Delaunay T;
  T.insert(Point(0, 0, 0));
  T.insert(points.begin(), points.end());
  time.stop();
  std::cout << "Delaunay 3D with origin: "
            << T.number_of_vertices() << " vertices in " << time.time() << " sec" << std::endl;

  time.reset();
  time.start();
  Delaunay_fast T_fast_on2;
  T_fast_on2.insert(Point(0, 0, 0));
  T_fast_on2.insert(points.begin(), points.end());
  time.stop();
  std::cout << "Delaunay 3D with origin, fast location: " << time.time() << std::endl;

  return EXIT_SUCCESS;
}
