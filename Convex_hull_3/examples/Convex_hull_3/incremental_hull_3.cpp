#include <CGAL/Homogeneous.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/algorithm.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_incremental_3.h>
#include <vector>

#ifdef CGAL_USE_GMP
#include <CGAL/Gmpz.h>
typedef CGAL::Gmpz RT;
#else
#include <CGAL/MP_Float.h>
typedef CGAL::MP_Float RT;
#endif


typedef CGAL::Homogeneous<RT>                  K;
typedef K::Point_3                             Point_3;
typedef CGAL::Polyhedron_3< K>                 Polyhedron;
typedef CGAL::Creator_uniform_3<int, Point_3>  Creator;

int main()
{
  CGAL::Random_points_in_sphere_3<Point_3, Creator> gen(100.0);

  std::vector<Point_3> V;
  // generate 250 points randomly on a sphere of radius 100.0 and copy
  // them to a vector
  CGAL::cpp11::copy_n( gen, 250, std::back_inserter(V) );

  Polyhedron P; // define polyhedron to hold convex hull

  // compute convex hull
  CGAL::convex_hull_incremental_3( V.begin(), V.end(), P, true);


  return 0;
}
