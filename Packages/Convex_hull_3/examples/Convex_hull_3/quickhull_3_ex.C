//
// file: examples/Convex_hull_3/ch_quickhull_3_ex.C
//
#include <CGAL/Homogeneous.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/copy_n.h>
#include <CGAL/Convex_hull_traits_3.h>
#include <CGAL/convex_hull_3.h>
#include <vector>
#ifdef CGAL_USE_LEDA
#include <CGAL/leda_integer.h>
typedef leda_integer RT;
#else
#ifdef CGAL_USE_GMP
#include <CGAL/Gmpz.h>
typedef CGAL::Gmpz RT;
#else
// NOTE: the choice of double here for a number type may cause problems
//       for degenerate point sets
#include <CGAL/double.h>
typedef double RT;
#endif
#endif


typedef CGAL::Homogeneous<RT>                     K;
typedef CGAL::Convex_hull_traits_3<K>             Traits;
typedef Traits::Polyhedron_3                      Polyhedron_3;
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
  Segment_3 segment;
  Polyhedron_3 polyhedron;
  if ( assign(segment, ch_object) )
     cout << "convex hull is a segment " << endl;
  else if ( assign (polyhedron, ch_object) )
     cout << "convex hull is a polyhedron " << endl;
  else
     cout << "convex hull error!" << endl;

  return 0;
}
