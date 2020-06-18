#include <CGAL/Simple_cartesian.h>
#include <CGAL/property_map.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/IO/write_xyz_points.h>
#include <CGAL/bilateral_smooth_point_set.h>
#include <CGAL/tags.h>

#include <utility> // defines std::pair
#include <fstream>

// Types
typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;

// Point with normal vector stored in a std::pair.
typedef std::pair<Point, Vector> PointVectorPair;

// Concurrency
typedef CGAL::Parallel_if_available_tag Concurrency_tag;

int main(int argc, char*argv[])
{
  const char* input_filename =  (argc>1)?argv[1]:"data/fin90_with_PCA_normals.xyz";
  const char* output_filename = (argc>2)?argv[2]:"data/fin90_with_PCA_normals_bilateral_smoothed.xyz";

  // Reads a .xyz point set file in points[] * with normals *.
  std::vector<PointVectorPair> points;
  std::ifstream stream(input_filename);
  if (!stream ||
      !CGAL::read_xyz_points(stream,
                     std::back_inserter(points),
                     CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PointVectorPair>()).
                     normal_map(CGAL::Second_of_pair_property_map<PointVectorPair>())))
  {
     std::cerr << "Error: cannot read file " << input_filename << std::endl;
     return EXIT_FAILURE;
  }

  // Algorithm parameters
  int k = 120;                 // size of neighborhood. The bigger the smoother the result will be.
                               // This value should bigger than 1.
  double sharpness_angle = 25; // control sharpness of the result.
                               // The bigger the smoother the result will be
  int iter_number = 3;         // number of times the projection is applied

  for (int i = 0; i < iter_number; ++i)
  {
    /* double error = */
    CGAL::bilateral_smooth_point_set <Concurrency_tag>(
      points,
      k,
      CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PointVectorPair>()).
      normal_map(CGAL::Second_of_pair_property_map<PointVectorPair>()).
      sharpness_angle (sharpness_angle));
  }

  //// Save point set.
  std::ofstream out(output_filename);
  out.precision(17);
  if (!out ||
      !CGAL::write_xyz_points(
      out, points,
      CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PointVectorPair>()).
      normal_map(CGAL::Second_of_pair_property_map<PointVectorPair>())))
  {
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

