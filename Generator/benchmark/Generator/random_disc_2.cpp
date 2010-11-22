
#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <CGAL/Simple_cartesian.h>
#include <cassert>
#include <vector>
#include <algorithm>
#include <CGAL/point_generators_2.h>
#include <CGAL/algorithm.h>

using namespace CGAL;

typedef Simple_cartesian<double>         K;
typedef K::Point_2                    Point;
typedef Creator_uniform_2<int,Point>  Creator;


void
disk(int N, double radius)
{
  std::cout.precision(12);
  CGAL::Random_points_in_disc_2<Point> rng(radius);
  std::cout << N << std::endl;
  for(double i = 0; i < N; i++){
    std::cout << *rng << std::endl;
    ++rng;
  }
}

int main(int argc, char* argv[]) 
{
  int N= 10;
  double radius = 1;
  try {

    po::options_description desc("Allowed options");
    desc.add_options()
      ("help", "Generator of perturbed points in a disk")
      ("N", po::value<int>(), "generate N points")
      ("radius", po::value<double>(), "radius of the disc")
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

    if (vm.count("radius")) {
      radius = vm["radius"].as<double>();
    }
  }
  catch(std::exception& e) {
    std::cerr << "error: " << e.what() << "\n";
    return 1;
  }

  disk(N, radius);
  return 0;
}
