//
// file: examples/Convex_hull_3/convex_hull_3_from_d_ex.C
//
#include <CGAL/Cartesian.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/copy_n.h>
#include <CGAL/Polyhedron_default_traits_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_3_from_d.h>
#include <vector>

// NOTE: the choice of double here for a number type may cause problems 
//       for degenerate point sets
typedef CGAL::Cartesian<double>                R;
typedef CGAL::Polyhedron_default_traits_3<R>   PolyTraits;
typedef CGAL::Polyhedron_3< PolyTraits >       Polyhedron;

int main()
{
  /* generate 250 points randomly on a sphere of radius 100.0 */
  CGAL::Random_points_in_sphere_3<Point> gen(100.0);

  std::vector<Point> V;
  CGAL::copy_n( gen, 250, std::back_inserter(V) ); /* copy them to a vector */
  
  Polyhedron P; /* define polyhedron to hold convex hull */

  /* compute convex hull */
  CGAL::convex_hull_3_from_d( V.begin(), V.end(), P, true);

  return 0;
}
