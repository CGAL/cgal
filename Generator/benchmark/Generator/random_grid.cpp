
#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <CGAL/Simple_cartesian.h>
#include <cassert>
#include <vector>
#include <algorithm>
#include <CGAL/point_generators_2.h>
#include <CGAL/algorithm.h>

using namespace CGAL;

typedef Simple_cartesian<int>         K;
typedef K::Point_2                    Point;
typedef Creator_uniform_2<int,Point>  Creator;


void
grid(int N, double eps)
{
  CGAL::Random rng;
  std::cout << N*N << std::endl;
  for(double i = 0; i < N; i++){
    for(double j = 0; j < N; j++){
      std::cout << i + rng.get_double(-eps,eps) << " "
                << j + rng.get_double(-eps,eps) << "\n";
    }
  }
}

int main(int argc, char* argv[])
{
  int N= 10;
  double eps = 0;
  try {

    po::options_description desc("Allowed options");
    desc.add_options()
      ("help", "Generator of perturbed points on a grid")
      ("N", po::value<int>(), "generate a grid with N x N points")
      ("eps", po::value<double>(), "perturb x and y of points by eps")
      ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
      std::cout << desc << "\n";
      return 1;
    }

    if (vm.count("N")) {
      N = vm["N"].as<int>();
    }

    if (vm.count("eps")) {
      eps = vm["eps"].as<double>();
    }
  }
  catch(std::exception& e) {
    std::cerr << "error: " << e.what() << "\n";
    return 1;
  }

  grid(N,eps);
  return 0;
}
