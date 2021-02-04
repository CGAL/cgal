#include <CGAL/Simple_cartesian.h>

#include <CGAL/IO/read_points.h>
#include <CGAL/IO/write_points.h>
#include <CGAL/property_map.h>
#include <CGAL/Aff_transformation_3.h>

#include <CGAL/pointmatcher/register_point_sets.h>

#include <CGAL/OpenGR/compute_registration_transformation.h>

#include <fstream>
#include <iostream>
#include <utility>
#include <vector>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point_3;
typedef K::Vector_3 Vector_3;
typedef std::pair<Point_3, Vector_3> Pwn;
typedef CGAL::First_of_pair_property_map<Pwn> Point_map;
typedef CGAL::Second_of_pair_property_map<Pwn> Normal_map;

namespace params = CGAL::parameters;

int main(int argc, const char** argv)
{
  const std::string fname1 = (argc>1) ? argv[1] : CGAL::data_file_path("points_3/hippo1.ply");
  const std::string fname2 = (argc>2) ? argv[2] : CGAL::data_file_path("points_3/hippo2.ply");

  std::vector<Pwn> pwns1, pwns2;
  if(!CGAL::IO::read_points(fname1, std::back_inserter(pwns1),
                            CGAL::parameters::point_map(CGAL::First_of_pair_property_map<Pwn>())
                                             .normal_map (Normal_map())))
  {
    std::cerr << "Error: cannot read file " << fname1 << std::endl;
    return EXIT_FAILURE;
  }

  if(!CGAL::IO::read_points(fname2, std::back_inserter(pwns2),
                            CGAL::parameters::point_map(Point_map())
                                             .normal_map(Normal_map())))
  {
    std::cerr << "Error: cannot read file " << fname2 << std::endl;
    return EXIT_FAILURE;
  }

  std::cerr << "Computing registration transformation using OpenGR Super4PCS.." << std::endl;
  // First, compute registration transformation using OpenGR Super4PCS
  K::Aff_transformation_3 res =
    std::get<0>( // get first of pair, which is the transformation
      CGAL::OpenGR::compute_registration_transformation
        (pwns1, pwns2,
         params::point_map(Point_map()).normal_map(Normal_map()),
         params::point_map(Point_map()).normal_map(Normal_map())));

  std::cerr << "Computing registration transformation using PointMatcher ICP, "
            << "taking transformation computed by OpenGR Super4PCS as initial transformation.." << std::endl;
  // Then, compute registration transformation using PointMatcher ICP, taking transformation computed
  // by OpenGR as initial transformation, and apply the transformation to pwns2
  // bool converged =
    CGAL::pointmatcher::register_point_sets
      (pwns1, pwns2,
       params::point_map(Point_map()).normal_map(Normal_map()),
       params::point_map(Point_map()).normal_map(Normal_map()).transformation(res));

  if(!CGAL::IO::write_points("pwns2_aligned.ply", pwns2,
                             CGAL::parameters::point_map(Point_map())
                                              .normal_map(Normal_map())))
  {
    return EXIT_FAILURE;
  }

  std::cerr << "Transformed version of " << fname2
            << " written to pwn2_aligned.ply.\n";

  return EXIT_SUCCESS;
}
