// 
// file: demo/Convex_hull_3/incremental_3_demo.C
//
#include <CGAL/Cartesian.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Convex_hull_d.h>
#include <CGAL/Convex_hull_d_traits_3.h>
#include <CGAL/Convex_hull_d_to_polyhedron_3.h>
#include <CGAL/Polyhedron_default_traits_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Geomview_stream.h>
#include <CGAL/IO/Polyhedron_geomview_ostream.h>

// NOTE: the choice of double here for a number type may cause problems 
//       for degenerate point sets
typedef CGAL::Cartesian<double>                K;
typedef K::Point_3                             Point_3;
typedef CGAL::Polyhedron_default_traits_3<K>   PolyTraits;
typedef CGAL::Polyhedron_3< PolyTraits >       Polyhedron_3;

typedef CGAL::Convex_hull_d_traits_3<K>        Hull_traits_3;
typedef CGAL::Convex_hull_d< Hull_traits_3 >   Convex_hull_3;

int main ()
{
  Convex_hull_3 CH(3);  // create instance of the class with dimension == 3

  // generate 250 points randomly on a sphere of radius 100.0 
  // and insert them into the convex hull
  CGAL::Random_points_in_sphere_3<Point_3> gen(100.0);
  for (int i = 0; i < 250 ; i++, gen++)
     CH.insert(*gen);

  assert(CH.is_valid());

  // define polyhedron to hold convex hull and create it
  Polyhedron_3 P; 
  CGAL::convex_hull_d_to_polyhedron_3(CH,P);

  // display polyhedron in a geomview window
  CGAL::Geomview_stream geomview;
  geomview << CGAL::RED;
  geomview << P;

  std::cout << "Press any key to end the program: ";
  char wait;
  std::cin >> wait;
}

