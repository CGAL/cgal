#include <CGAL/Kinetic/basic.h>

#include <CGAL/Random.h>
#include <algorithm>
#include <CGAL/Polynomial/Polynomial.h>
#if CGAL_USE_BOOST_PROGRAM_OPTIONS
#include <boost/program_options.hpp>
#endif


template <class NT>
void generate(int n, int d, int num_coords)
{
  CGAL::Random rand;
  for (int i=0; i< n; ++i) {
    std::vector<std::vector<double> > coefs(num_coords);
    for (int j=0; j<= d; ++j) {
      for (int k=0; k < num_coords; ++k) {
	coefs[k].push_back((rand.get_double()*10-5)/(j+1));
      }
    }
    for (int j=0; j< num_coords; ++j) {
      CGAL::POLYNOMIAL::Polynomial<NT> p(coefs[j].begin(), coefs[j].end());
      std::cout << p;
      if (j != num_coords-1) {
	std::cout << ", ";
      }
      else {
	std::cout << std::endl;
      }
    }

  }
}


int main(int argc, char *argv[])
{
  int n=10;
  int d=2;
  bool threed=false;
  bool inexact=true;
  bool weighted=false;

#ifdef CGAL_USE_BOOST_PROGRAM_OPTIONS
  bool print_help=false;

  boost::program_options::options_description desc("Allowed options");
  desc.add_options()
    ("help", boost::program_options::bool_switch(&print_help), "produce help message")
    ("num-points,n", boost::program_options::value<int>(&n), "Number of points to use.")
    ("float,f", boost::program_options::bool_switch(&inexact), "Output the coordinates as floats rather than rationals.")
    ("three-dimensions,3", boost::program_options::bool_switch(&threed), "Write three dimensional points.")
    ("weighted,w", boost::program_options::bool_switch(&weighted), "Write weighted points.")
    ("degree,d", boost::program_options::value<int>(&d), "The degree of the motions to use.");

  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::command_line_parser(argc, argv).
				options(desc).run(), vm);
  boost::program_options::notify(vm);


  if (print_help) {
    std::cout << "This program generates a set of moving points and outputs it to a file.\n";
    std::cout << desc << "\n";
    return EXIT_FAILURE;
  }

#else

  bool boost_program_options_disabled;
#endif

  int num_coords=2;
  if (weighted) ++num_coords;
  if (threed) ++num_coords;
  if (inexact) {
    generate<double>(n,d,num_coords);
  }
  else {
    generate<CGAL::Kinetic::Default_field_nt>(n,d,num_coords);
  }

  return EXIT_SUCCESS;
}
