#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/property_map.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/IO/write_xyz_points.h>
#include <CGAL/bilateral_smooth_point_set.h>
#include <CGAL/tags.h>

#include <utility> // defines std::pair
#include <list>
#include <string>
#include <fstream>

// Types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;

// Point with normal vector stored in a std::pair.
typedef std::pair<Point, Vector> PointVectorPair;


int main(void)
{
  const std::string INPUT_FILENAME_WITHOUT_EXT = "data/fin90_with_PCA_normals";

  // Reads a .xyz point set file in points[].
  std::vector<PointVectorPair> points;
  std::ifstream stream(INPUT_FILENAME_WITHOUT_EXT + ".xyz");
  if (!stream ||
      !CGAL::read_xyz_points_and_normals(stream,
                     std::back_inserter(points),
                     CGAL::First_of_pair_property_map<PointVectorPair>(),
                     CGAL::Second_of_pair_property_map<PointVectorPair>()))
  {
     std::cerr << "Error: cannot read file " 
               << INPUT_FILENAME_WITHOUT_EXT << ".xyz" << std::endl;
     return EXIT_FAILURE;
  }

  // Algorithm parameters
  int k = 120;                 // size of neighborhood 
  double sharpness_angle = 25; // control sharpness of the result
                               // the bigger the smoother the result will be
  int iter_number = 3;         // number of times the projection is applied
  
  for (int i = 0; i < iter_number; ++i)
  {
    double error = 
    CGAL::bilateral_smooth_point_set <CGAL::Parallel_tag>(
          points.begin(), 
          points.end(),
          CGAL::First_of_pair_property_map<PointVectorPair>(),
          CGAL::Second_of_pair_property_map<PointVectorPair>(),
          k,
          sharpness_angle);
  }
  
  //// Save point set.
  std::ofstream out(INPUT_FILENAME_WITHOUT_EXT + "_bilateral DENOISED.xyz");   
  if (!out ||
      !CGAL::write_xyz_points_and_normals(
      out, points.begin(), points.end(), 
      CGAL::First_of_pair_property_map<PointVectorPair>(),
      CGAL::Second_of_pair_property_map<PointVectorPair>()))
  {
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

