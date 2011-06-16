#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/algorithm.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_3.h>
#include <vector>


typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef CGAL::Polyhedron_3<K>                     Polyhedron_3;
typedef K::Segment_3                              Segment_3;

// define point creator
typedef K::Point_3                                Point_3;
typedef CGAL::Creator_uniform_3<double, Point_3>  PointCreator;


int main()
{
  CGAL::Random_points_in_sphere_3<Point_3, PointCreator> gen(100.0);

  // generate 250 points randomly on a sphere of radius 100.0
  // and copy them to a vector
  std::vector<Point_3> points;
  CGAL::copy_n( gen, 250, std::back_inserter(points) );

  // define object to hold convex hull
  CGAL::Object ch_object;

  // compute convex hull
  CGAL::convex_hull_3(points.begin(), points.end(), ch_object);

  // determine what kind of object it is
  if ( const Segment_3* segment=CGAL::object_cast<Segment_3>(&ch_object) ){
     std::cout << "convex hull is the segment " << *segment << std::endl;
  }
  else if (const Polyhedron_3* poly = CGAL::object_cast<Polyhedron_3>(&ch_object) )
     std::cout << "convex hull is a polyhedron with " 
               << poly->size_of_vertices() << " vertices" << std::endl;
  else
     std::cout << "convex hull error!" << std::endl;

  return 0;
}
