// TODO: Copyright info

#include <CGAL/Simple_cartesian.h>
#include <CGAL/IO/read_ply_points.h>
#include <CGAL/IO/write_ply_points.h>
#include <CGAL/property_map.h>

#include <CGAL/pointmatcher/compute_registration_transformation.h>

#include <fstream>
#include <iostream>
#include <utility>
#include <vector>

typedef CGAL::Simple_cartesian<float> K;
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

  // TODO: Another example with default settings, refering to default settings in doc
  
  //
  // Prepare ICP config
  //
  using CGAL::pointmatcher::ICP_config;

  // TODO: Change naming convention of filters? (for example, reference data points filters -> reference point set filters)
  // TODO: Get reference point set filters from np2?
  // TODO: Refer to https://github.com/ethz-asl/libpointmatcher/blob/master/doc/Configuration.md while doc

  // Prepare reading data points filters
  std::vector<ICP_config> reading_data_points_filters;
  reading_data_points_filters.push_back( ICP_config { .name="MinDistDataPointsFilter"       , .params={ {"minDist", "0.5" }}  } );
  reading_data_points_filters.push_back( ICP_config { .name="RandomSamplingDataPointsFilter", .params={ {"prob"   , "0.05"}}  } );

  // Prepare reference data points filters
  std::vector<ICP_config> reference_data_points_filters;
  reference_data_points_filters.push_back( ICP_config { .name="MinDistDataPointsFilter"       , .params={ {"minDist", "0.5" }}  } );
  reference_data_points_filters.push_back( ICP_config { .name="RandomSamplingDataPointsFilter", .params={ {"prob"   , "0.05"}}  } );

	// Prepare matcher function
  ICP_config matcher { .name="KDTreeMatcher", .params={ {"knn", "1"}, {"epsilon", "3.16"} } };

  // Prepare outlier filters
  std::vector<ICP_config> outlier_filters;
  outlier_filters.push_back( ICP_config { .name="TrimmedDistOutlierFilter", .params={ {"ratio", "0.75" }}  } );

  // Prepare error minimizer
  ICP_config error_minimizer { .name = "PointToPointErrorMinimizer"};

  // Prepare transformation checker
  std::vector<ICP_config> transformation_checkers;
  transformation_checkers.push_back( ICP_config { .name="CounterTransformationChecker", .params={ {"maxIterationCount", "150" }}  } );
  transformation_checkers.push_back( ICP_config { .name="DifferentialTransformationChecker", .params={ {"minDiffRotErr"  , "0.001" },
                                                                                                       {"minDiffTransErr", "0.01"  },
                                                                                                       {"smoothLength"   , "4"     } }  
                                                } );
  // Prepare inspector
  ICP_config inspector { .name="NullInspector" };

  // Prepare logger
  ICP_config logger { .name= "FileLogger" };


K::Aff_transformation_3 res = 
  CGAL::pointmatcher::compute_registration_transformation
    (pwns1, pwns2, 
     params::point_map(Point_map()).normal_map(Normal_map())
     .pm_reading_data_points_filters(reading_data_points_filters)
     .pm_reference_data_points_filters(reference_data_points_filters)
     .pm_matcher(matcher)
     .pm_outlier_filters(outlier_filters)
     .pm_error_minimizer(error_minimizer)
     .pm_transformation_checkers(transformation_checkers)
     .pm_inspector(inspector)
     .pm_logger(logger),
     params::point_map(Point_map()).normal_map(Normal_map()));

  std::ofstream out("pwns2_aligned.ply");
  if (!out ||
      !CGAL::write_ply_points(
        out, pwns2,
        CGAL::parameters::point_map(Point_map()).
        normal_map(Normal_map())))
  {
    return EXIT_FAILURE;
  }

  std::cout << "Transformed version of " << fname2
            << " written to pwn2_aligned.ply.\n";

  return EXIT_SUCCESS;
}