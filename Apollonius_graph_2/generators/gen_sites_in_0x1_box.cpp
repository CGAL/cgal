#include <CGAL/basic.h>
#include <CGAL/Apollonius_graph_2/random_sites_in_0x1_box.h>

#include <iostream>
#include <sstream>
#include <iomanip>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Apollonius_site_2.h>

typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Apollonius_site_2<K>   Site_2;

int usage (int argc, char **argv)
{
    std::cerr << "usage: " << argv[0]
	      << " <number of points> <max radius> [seed]" << std::endl;
    return 2;
}


int main (int argc, char **argv)
{
    int num, seed = 42;
    double rmax;

    if (argc < 3)
        return usage (argc, argv);

    {
      std::istringstream is(argv[1]);
      if ( !(is >> num) ) { return usage(argc, argv); }
    }
    {
      std::istringstream is(argv[2]);
      if ( !(is >> rmax) ) { return usage(argc, argv); }
    }
    if (argc > 3) {
      std::istringstream is(argv[3]);
      if ( !(is >> seed) ) { return usage(argc, argv); }
    }

    CGAL::Random_sites_in_0x1_box<Site_2> g(rmax, seed);

    std::cout << std::setprecision(17);

    for (int i = 0; i < num; ++i) {
      std::cout << *g << std::endl;
    }

    return 0;
}
