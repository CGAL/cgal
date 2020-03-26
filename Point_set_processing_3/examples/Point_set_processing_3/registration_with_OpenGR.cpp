#include <CGAL/Simple_cartesian.h>
#include <CGAL/IO/read_ply_points.h>
#include <CGAL/IO/write_ply_points.h>
#include <CGAL/property_map.h>

#include <CGAL/OpenGR/compute_registration_transformation.h>
#include <CGAL/OpenGR/register_point_sets.h>

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
  const char* fname1 = (argc>1)?argv[1]:"data/hippo1.ply";
  const char* fname2 = (argc>2)?argv[2]:"data/hippo2.ply";

  std::vector<Pwn> pwns1, pwns2;
  std::ifstream input(fname1);
  if (!input ||
      !CGAL::read_ply_points(input, std::back_inserter(pwns1),
            CGAL::parameters::point_map (CGAL::First_of_pair_property_map<Pwn>()).
            normal_map (Normal_map())))
  {
    std::cerr << "Error: cannot read file " << fname1 << std::endl;
    return EXIT_FAILURE;
  }
  input.close();

  input.open(fname2);
  if (!input ||
      !CGAL::read_ply_points(input, std::back_inserter(pwns2),
            CGAL::parameters::point_map (Point_map()).
            normal_map (Normal_map())))
  {
    std::cerr << "Error: cannot read file " << fname2 << std::endl;
    return EXIT_FAILURE;
  }
  input.close();

  // EITHER call the registration method Super4PCS from OpenGR to get the transformation to apply to pwns2
  // std::pair<K::Aff_transformation_3, double> res =
    CGAL::OpenGR::compute_registration_transformation(pwns1, pwns2,
                                                      params::point_map(Point_map())
                                                      .normal_map(Normal_map())
                                                      .number_of_samples(200)
                                                      .maximum_running_time(60)
                                                      .accuracy(0.01),
                                                      params::point_map(Point_map())
                                                      .normal_map(Normal_map()));

  // OR call the registration method Super4PCS from OpenGR and apply the transformation to pwn2
  double score =
    CGAL::OpenGR::register_point_sets(pwns1, pwns2,
                                      params::point_map(Point_map())
                                      .normal_map(Normal_map())
                                      .number_of_samples(200)
                                      .maximum_running_time(60)
                                      .accuracy(0.01),
                                      params::point_map(Point_map())
                                      .normal_map(Normal_map()));

  std::ofstream out("pwns2_aligned.ply");
  if (!out ||
      !CGAL::write_ply_points(
        out, pwns2,
        CGAL::parameters::point_map(Point_map()).
        normal_map(Normal_map())))
  {
    return EXIT_FAILURE;
  }

  std::cout << "Registration score: " << score << ".\n"
            << "Transformed version of " << fname2
            << " written to pwn2_aligned.ply.\n";

  return EXIT_SUCCESS;
}
