#include <CGAL/Cartesian.h>
#include <CGAL/Point_3.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/copy_n.h>
#include <CGAL/Halfedge_data_structure_polyhedron_default_3.h>
#include <CGAL/Polyhedron_default_traits_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_3_from_d.h>
#include <vector>

/* representation class */
typedef CGAL::Cartesian<double>                                R;

/* define polyhedron type */
typedef CGAL::Halfedge_data_structure_polyhedron_default_3<R>  HDS;
typedef CGAL::Polyhedron_default_traits_3<R>                   PolyTraits;
typedef CGAL::Polyhedron_3< PolyTraits, HDS>                   Polyhedron;

/* define point creator */
typedef CGAL::Point_3<R>                                       Point;
typedef CGAL::Creator_uniform_3<double,Point>                  PointCreator;

int 
main()
{
  /* generate 250 points randomly on a sphere of radius 100.0 */
  CGAL::Random_points_in_sphere_3<Point, PointCreator> gen(100.0);

  /* and copy them to a vector */
  std::vector<Point> V;
  CGAL::copy_n( gen, 250, std::back_inserter(V) );
  
  /* define polyhedron to hold convex hull */
  Polyhedron P;

  /* compute convex hull */
  CGAL::convex_hull_3_from_d( V.begin(), V.end(), P, true);

  return 0;
}
