#include <CGAL/Simple_cartesian.h>

#include <CGAL/constraint_based_smooth_point_set2.h>
#include <CGAL/IO/read_points.h>
#include <CGAL/IO/write_points.h>
#include <CGAL/property_map.h>
#include <CGAL/tags.h>

#include <utility> // defines std::pair
#include <string>
#include <fstream>
#include <chrono>

// Types
typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef std::array<unsigned char, 3> Color;

typedef std::tuple<Point, Vector, Color> PNC;
typedef CGAL::Nth_of_tuple_property_map<0, PNC> Point_map;
typedef CGAL::Nth_of_tuple_property_map<1, PNC> Normal_map;
typedef CGAL::Nth_of_tuple_property_map<2, PNC> Color_map;

// Concurrency
typedef CGAL::Parallel_if_available_tag Concurrency_tag;

// Define how a color should be stored
namespace CGAL {

template< class F >
struct Output_rep< ::Color, F > {
  const ::Color& c;
  static const bool is_specialized = true;
  Output_rep (const ::Color& c) : c(c)
  { }
  std::ostream& operator() (std::ostream& out) const
  {
    if (IO::is_ascii(out))
      out << int(c[0]) << " " << int(c[1]) << " " << int(c[2]);
    else
      out.write(reinterpret_cast<const char*>(&c), sizeof(c));
    return out;
  }
};

} // namespace CGAL

int main(int argc, char*argv[])
{
  const std::string input_filename = (argc>1) ? argv[1] : CGAL::data_file_path("points_3/fin90_with_PCA_normals.xyz");
  const char* output_filename = (argc>2) ? argv[2] : "data/fin90_with_PCA_normals_constraint_based_smoothed.ply";

  // Reads a point set file in points[] * with normals *.
  std::vector<PNC> points;
  if(!CGAL::IO::read_points(input_filename, std::back_inserter(points),
                            CGAL::parameters::point_map(CGAL::Nth_of_tuple_property_map<0, PNC>())
                                             .normal_map(CGAL::Nth_of_tuple_property_map<1, PNC>())))
  {
     std::cerr << "Error: cannot read file " << input_filename << std::endl;
     return EXIT_FAILURE;
  }

  // Algorithm parameters
  const int iter_number = 150;

  auto start = std::chrono::high_resolution_clock::now();

  for(int i = 0; i < iter_number; ++i)
  {
    auto curr_start = std::chrono::high_resolution_clock::now();
    /* double error = */
    CGAL::constraint_based_smooth_point_set <Concurrency_tag>(
      points,
      CGAL::parameters::point_map(CGAL::Nth_of_tuple_property_map<0, PNC>())
                       .normal_map(CGAL::Nth_of_tuple_property_map<1, PNC>()));

    std::cout << i << std::endl;

    auto curr_stop = std::chrono::high_resolution_clock::now();
    std::cout << "Iteration time (ms): " << std::chrono::duration_cast<std::chrono::milliseconds>(curr_stop - curr_start).count() << std::endl;

    std::cout << "iterations/iteration_" + std::to_string(i) + ".ply" << std::endl;
    std::ofstream f("data/iterations/iteration_" + std::to_string(i) + ".ply", std::ios::binary);
    CGAL::IO::set_binary_mode(f); // The PLY file will be written in the binary format
    if(!CGAL::IO::write_PLY_with_properties(f, points,
                                        CGAL::make_ply_point_writer (Point_map()),
                                        CGAL::make_ply_normal_writer(Normal_map()),
                                        std::make_tuple(Color_map(),
                                                        CGAL::IO::PLY_property<unsigned char>("red"),
                                                        CGAL::IO::PLY_property<unsigned char>("green"),
                                                        CGAL::IO::PLY_property<unsigned char>("blue"))))
    {
      std::cout << "AHHHHHHHHHHHHHHHHHH" << std::endl;
      return EXIT_FAILURE;
    }

  }

  auto stop = std::chrono::high_resolution_clock::now();

  for(size_t i=0;i<points.size();++i){
    Color c = get<2>(points[i]);
    // std::cout << int(c[0]) << " " << int(c[1]) << " " << int(c[2]) << std::endl;
  }

  std::cout << "Average iteration time (ms): " << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count() / iter_number << std::endl;

  //// Save point set.
  // if(!CGAL::IO::write_XYZ(output_filename, points,
  //                         CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PointVectorPair>())
  //                                          .normal_map(CGAL::Second_of_pair_property_map<PointVectorPair>())
  //                                          .stream_precision(17)))
  //   return EXIT_FAILURE;

  std::ofstream f(output_filename, std::ios::binary);
  CGAL::IO::set_binary_mode(f); // The PLY file will be written in the binary format
  if(!CGAL::IO::write_PLY_with_properties(f, points,
                                      CGAL::make_ply_point_writer (Point_map()),
                                      CGAL::make_ply_normal_writer(Normal_map()),
                                      std::make_tuple(Color_map(),
                                                      CGAL::IO::PLY_property<unsigned char>("red"),
                                                      CGAL::IO::PLY_property<unsigned char>("green"),
                                                      CGAL::IO::PLY_property<unsigned char>("blue"))))
  {
    std::cout << "AHHHHHHHHHHHHHHHHHH" << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}