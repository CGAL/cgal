#include <CGAL/basic.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Convex_hull_d_traits_3.h>
#include <CGAL/Convex_hull_d.h>
#include <CGAL/random_selection.h>

#include <CGAL/point_generators_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Convex_hull_d_to_polyhedron_3.h>
#include <CGAL/IO/Polyhedron_geomview_ostream.h>
#include <iostream>
#include <string>

#if !defined(__BORLANDC__) && !defined(_MSC_VER)

#ifdef CGAL_USE_LEDA
#include <CGAL/leda_integer.h>
typedef leda_integer RT;
#else
#ifdef CGAL_USE_GMP
#include <CGAL/Gmpz.h>
typedef CGAL::Gmpz RT;
#else
#include <CGAL/double.h>
typedef double RT;
#endif
#endif

typedef CGAL::Homogeneous<RT> Kernel_3;
typedef CGAL::Convex_hull_d_traits_3<Kernel_3> Kernel_d_3;
typedef CGAL::Convex_hull_d<Kernel_d_3> Convex_hull_d;

typedef  CGAL::Point_3<Kernel_3> Point_3;
typedef  CGAL::Polyhedron_3<Kernel_3> Polyhedron;

typedef CGAL::Creator_uniform_3<RT,Point_3> Creator;
typedef CGAL::Random_points_in_cube_3<Point_3,Creator> Point_source;

int main(int argc, char* argv[]) {
  int dimension = 3;  
  int n = 100; 
  if (argc > 1 && std::string(argv[1])=="-h") {
    std::cout << "usage: ch5-demo [#points]\n";
    exit(1);
  }
  if (argc > 1) n = atoi(argv[1]);
 
  int r = 2*n;
  CGAL::Geomview_stream gv(CGAL::Bbox_3(-r, -r, -r, r, r, r));
  gv.clear();


  Convex_hull_d T(dimension);  
  Point_source g(n);
  for(int i=0; i<n; ++i) {
    T.insert(*g++);
    if (i%10==0) std::cout << i << " points inserted" << std::endl;
  }
  T.is_valid(true); 

  Polyhedron P;
  CGAL::convex_hull_d_to_polyhedron_3(T,P);
  gv << P;
  std::cout << "Enter a key to finish" << std::endl;
  char ch;
  std::cin >> ch;

  return 0;
}

#else // on windows:

int main(int argc, char* argv[]) {
  std::cerr << 
  "This demo requires geomview, that is is not present on windows\n";
  return 0;
}

#endif
