#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>

#include <vector>
#include <CGAL/Convex_hull_traits_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Object.h>
#include <CGAL/point_generators_3.h>

typedef CGAL::Cartesian<double>                 R;
typedef CGAL::Convex_hull_traits_3<R>           Traits;
typedef Traits::Polyhedron_3                    Polyhedron_3;

typedef R::Point_3                              Point_3;
typedef R::Segment_3                            Segment_3;

typedef CGAL::Creator_uniform_3<double,Point_3>               Creator;
typedef CGAL::Random_points_in_sphere_3<Point_3,Creator>      Generator;

const unsigned int num = 40;

void test_small_hull()
{
  std::vector<Point_3> points;
  points.push_back(Point_3(0,0,0));
  points.push_back(Point_3(1,0,0));
  points.push_back(Point_3(0,1,0));
  points.push_back(Point_3(0,0,1));
  points.push_back(Point_3(0.5,0.5,0.5));
  points.push_back(Point_3(0.2,0.5,0.3));
  points.push_back(Point_3(0.1,0.3,0.2));

  Polyhedron_3 polyhedron1;
  CGAL::convex_hull_3(points.begin(), points.end(), polyhedron1, Traits());
  assert ( polyhedron1.size_of_vertices() == 5 && 
           polyhedron1.size_of_facets() == 6 );
  Polyhedron_3 polyhedron2;
  CGAL::convex_hull_3(points.begin(), points.end(), polyhedron2);
  assert (CGAL::is_strongly_convex_3(polyhedron2)); // test default traits class
  assert ( polyhedron2.size_of_vertices() == 5 && 
           polyhedron2.size_of_facets() == 6 );
}


int main()
{
  test_small_hull();

  std::vector<Point_3> points;
  Generator g(5000.0);
  CGAL::copy_n( g, num, std::back_inserter(points));

  assert(points.size() == num);

  CGAL::Object ch_object;
  CGAL::convex_hull_3(points.begin(), points.end(), ch_object, Traits());
  CGAL::convex_hull_3(points.begin(), points.end(), ch_object);

  Segment_3 segment;

  Polyhedron_3 polyhedron;

  assert( assign(segment, ch_object) || assign(polyhedron, ch_object) );
  return 0;
}
