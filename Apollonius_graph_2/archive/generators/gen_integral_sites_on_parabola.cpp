#include <CGAL/basic.h>
#include <CGAL/Apollonius_graph_2/random_integral_sites_on_parabola_2.h>
#include <CGAL/Random.h>

#include <iostream>
#include <sstream>
#include <iomanip>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Apollonius_site_2.h>

typedef CGAL::Simple_cartesian<double>  K;
typedef CGAL::Apollonius_site_2<K>      Site_2;
typedef CGAL::Random                    Random;

int usage (int argc, char *argv[])
{
  std::cerr << "usage: " << argv[0]
	    << " <number of points> <bit size of coordinates> "
	    << "[bit size of perturbation] [seed]" << std::endl;
  return 2;
}

int main(int argc, char *argv[])
{
  int num, seed = 17;
  unsigned int b, p = 0;

  if ( argc < 3 ) {
    return usage(argc, argv);
  }

  {
    std::istringstream is(argv[1]);
    if ( !(is >> num) ) { return usage(argc, argv); }
  }

  {
    std::istringstream is(argv[2]);
    if ( !(is >> b) ) { return usage(argc, argv); }
  }

  if ( argc > 3 ) {
    std::istringstream is(argv[3]);
    if ( !(is >> p) ) { return usage(argc, argv); }
  }

  if ( argc > 4 ) {
    std::istringstream is(argv[4]);
    if ( !(is >> seed) ) { return usage(argc, argv); }
  }

  CGAL::Random_integral_sites_on_parabola_2<Site_2,Random> g(b, p, seed);

  std::cout << std::setprecision(17);
  for (int i = 0; i < num; ++i) {
    std::cout << *g << std::endl;
  }

  return 0;
}
