//
// file: examples/Convex_hull_3/incremental_hull_3_ex.C
//
#include <CGAL/Homogeneous.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/copy_n.h>
#include <CGAL/Polyhedron_default_traits_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_incremental_3.h>
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


typedef CGAL::Homogeneous<RT>                  K;
typedef K::Point_3                             Point_3;
typedef CGAL::Polyhedron_default_traits_3<K>   PolyTraits;
typedef CGAL::Polyhedron_3< PolyTraits >       Polyhedron;
typedef CGAL::Creator_uniform_3<int, Point_3>  Creator;

int main()
{
  CGAL::Random_points_in_sphere_3<Point_3, Creator> gen(100.0);

  std::vector<Point_3> V;
  // generate 250 points randomly on a sphere of radius 100.0 and copy 
  // them to a vector
  CGAL::copy_n( gen, 250, std::back_inserter(V) ); 
  
  Polyhedron P; // define polyhedron to hold convex hull 

  // compute convex hull 
  CGAL::convex_hull_incremental_3( V.begin(), V.end(), P, true);


  return 0;
}
