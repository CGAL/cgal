#include <CGAL/Cartesian.h>
#include <CGAL/Point_3.h>
#include <CGAL/Segment_3.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/copy_n.h>
#include <CGAL/Convex_hull_traits_3.h>
#include <CGAL/convex_hull_3.h>
#include <vector>

// NOTE: the choice of double here for a number type may cause problems 
//       for degenerate point sets
typedef CGAL::Cartesian<double>                   R;
typedef CGAL::Convex_hull_traits_3<R>             Traits;
typedef Traits::Polyhedron_3                      Polyhedron_3;
typedef CGAL::Segment_3<R>                        Segment_3;

/* define point creator */
typedef CGAL::Point_3<R>                          Point_3;
typedef CGAL::Creator_uniform_3<double, Point_3>  PointCreator;


int main()
{
  /* generate 250 points randomly on a sphere of radius 100.0 */
  CGAL::Random_points_in_sphere_3<Point_3, PointCreator> gen(100.0);

  /* and copy them to a vector */
  std::vector<Point_3> points;
  CGAL::copy_n( gen, 250, std::back_inserter(points) );
  
  /* define polyhedron to hold convex hull */
  CGAL::Object ch_object;

  /* compute convex hull */
  CGAL::convex_hull_3(points.begin(), points.end(), ch_object);

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
