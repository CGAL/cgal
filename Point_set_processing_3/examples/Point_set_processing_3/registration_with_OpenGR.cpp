#include <CGAL/Simple_cartesian.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/IO/write_xyz_points.h>
#include <CGAL/property_map.h>

#include <CGAL/OpenGR/align.h>

#include <fstream>
#include <iostream>
#include <utility>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point_3;
typedef K::Vector_3 Vector_3;
typedef std::pair<Point_3, Vector_3> Pwn;
typedef CGAL::First_of_pair_property_map<Pwn> Point_map;
typedef CGAL::Second_of_pair_property_map<Pwn> Normal_map;

namespace params = CGAL::parameters;

int main(int argc, const char** argv)
{
  const char* fname1 = (argc>1)?argv[1]:"data/hippo1.xyz";
  const char* fname2 = (argc>2)?argv[2]:"data/hippo2.xyz";

  std::vector<Pwn> pwns1, pwns2;
  std::ifstream input(fname1);
  if (!input ||
      !CGAL::read_xyz_points(input, std::back_inserter(pwns1),
            CGAL::parameters::point_map (CGAL::First_of_pair_property_map<Pwn>()).
            normal_map (Normal_map())))
  {
    std::cerr << "Error: cannot read file " << fname1 << std::endl;
    return EXIT_FAILURE;
  }
  input.close();

  input.open(fname2);
  if (!input ||
      !CGAL::read_xyz_points(input, std::back_inserter(pwns2),
            CGAL::parameters::point_map (Point_map()).
            normal_map (Normal_map())))
  {
    std::cerr << "Error: cannot read file " << fname2 << std::endl;
    return EXIT_FAILURE;
  }
  input.close();

  // OpenGR options
  GlobalRegistration::Match4PCSOptions options;
  bool overlap_OK = options.configureOverlap(0.7);
  if(!overlap_OK)
  {
    std::cerr << "Invalid overlap configuration.\n";
    return EXIT_FAILURE;
  }
  options.sample_size = 200;
  options.max_time_seconds = 1000;
  options.delta = 0.01;

  // call the registration method Super4PCS from OpenGR
  double score =
    CGAL::OpenGR::align(pwns1, pwns2,
                        params::point_map(Point_map())
                               .normal_map(Normal_map())
                               .opengr_options(options),
                        params::point_map(Point_map())
                               .normal_map(Normal_map()));

  std::ofstream out("pwns2_aligned.xyz");
  if (!out ||
      !CGAL::write_xyz_points(
        out, pwns2,
        CGAL::parameters::point_map(Point_map()).
        normal_map(Normal_map())))
  {
    return EXIT_FAILURE;
  }

  std::cout << "Registration score: " << score << ".\n"
            << "Transformed version of " << fname2
            << " written to pwn2_aligned.xyz.\n";

  return EXIT_SUCCESS;
}