// 
// file: demo/Convex_hull_3/incremental_3_demo.C
//
#include <CGAL/Homogeneous.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Convex_hull_d.h>
#include <CGAL/Convex_hull_d_traits_3.h>
#include <CGAL/Convex_hull_d_to_polyhedron_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/copy_n.h>
#include <CGAL/IO/Geomview_stream.h>
#include <CGAL/IO/Polyhedron_geomview_ostream.h>
#include <vector>
#include <cassert>

#if !defined(__BORLANDC__) && !defined(_MSC_VER)

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
typedef CGAL::Polyhedron_3< K>                 Polyhedron_3;

typedef CGAL::Convex_hull_d_traits_3<K>        Hull_traits_3;
typedef CGAL::Convex_hull_d< Hull_traits_3 >   Convex_hull_3;
typedef CGAL::Creator_uniform_3<double, Point_3>   Creator;

int main ()
{
  Convex_hull_3 CH(3);  // create instance of the class with dimension == 3

  // generate 250 points randomly on a sphere of radius 100 
  // and insert them into the convex hull
  CGAL::Random_points_in_sphere_3<Point_3, Creator> gen(100);

  for (int i = 0; i < 250 ; i++, ++gen)
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
  std::cout.flush();
  char ch;
  std::cin.get(ch);

  return 0;
}

#else // on windows:

int main() {
  std::cerr <<
  "This demo requires geomview, which is is not present on windows\n";
  return 0;
}

#endif

