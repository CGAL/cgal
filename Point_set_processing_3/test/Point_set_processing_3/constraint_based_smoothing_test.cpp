#include <CGAL/Simple_cartesian.h>

#include <CGAL/constraint_based_smooth_point_set.h>
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

namespace internal {

  template <typename Kernel>
  typename Kernel::FT measure_sphere_normal_deviation(
    std::vector<PNC> points
  )
  {
    // basic geometric types
    typedef typename Kernel::FT FT;
    typedef typename Kernel::Vector_3 Vector;
    typedef typename Kernel::Point_3 Point;

    FT cum_angle_diff = 0;

    for (PNC point : points)
    {
      Point p = get<0>(point);
      Vector vp(p[0], p[1], p[2]);

      Vector n = get<1>(point);

      cum_angle_diff += CGAL::approximate_angle(n, vp);
    }

    return cum_angle_diff / points.size();
  }

  template <typename Kernel>
  typename Kernel::FT measure_sphere_point_deviation(
    std::vector<PNC> points
  )
  {
    // basic geometric types
    typedef typename Kernel::FT FT;
    typedef typename Kernel::Vector_3 Vector;
    typedef typename Kernel::Point_3 Point;

    FT cum_point_diff = 0;

    for (PNC point : points)
    {
      Point p = get<0>(point);
      
      cum_point_diff += std::abs(p[0]*p[0] + p[1]*p[1] + p[2]*p[2] - 1); 
    }

    return cum_point_diff / points.size();
  }

} // namespace CGAL

int main(int argc, char*argv[])
{
  const std::string input_filename = CGAL::data_file_path("points_3/sphere_20k_normal.xyz");
  const char* output_filename = (argc>2) ? argv[2] : "data/sphere_20k_normal_constraint_based_smoothed.ply";

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
  const int iter_number = 5;

  // saved parameters for sphere
  // FT neighbor_radius = 0;
  // FT normal_threshold = 0.9 * (180/M_PI);
  // FT damping_factor = 3;
  // FT eigenvalue_threshold_nvt = 0.5;
  // FT eigenvalue_threshold_covm = .25;
  // FT update_threshold = 0;

  Kernel::FT original_angle_diff = internal::measure_sphere_normal_deviation<Kernel>(points);
  std::cout << "Original avg. angle deviation (deg): " << original_angle_diff << std::endl;

  Kernel::FT original_point_diff = internal::measure_sphere_point_deviation<Kernel>(points);
  std::cout << "Original avg. point deviation (deg): " << original_point_diff << std::endl;

  for(int i = 0; i < iter_number; ++i)
  {
    CGAL::constraint_based_smooth_point_set <Concurrency_tag>(
      points,
      CGAL::parameters::point_map(CGAL::Nth_of_tuple_property_map<0, PNC>())
                       .normal_map(CGAL::Nth_of_tuple_property_map<1, PNC>()));

    Kernel::FT curr_angle_diff = internal::measure_sphere_normal_deviation<Kernel>(points);
    std::cout << "Avg. angle deviation (deg): " << curr_angle_diff << std::endl;

    Kernel::FT curr_point_diff = internal::measure_sphere_point_deviation<Kernel>(points);
    std::cout << "Avg. point deviation: " << curr_point_diff << std::endl;

    if (curr_angle_diff > original_angle_diff)
    {
      std::cout << "TEST FAILURE: Normal deviation grew." << std::endl;
      // return EXIT_FAILURE;
    }
    if (curr_point_diff > original_point_diff)
    {
      std::cout << "TEST FAILURE: Point deviation grew." << std::endl;
      // return EXIT_FAILURE;
    }
  }

  for(size_t i=0;i<points.size();++i){
    Color c = get<2>(points[i]);
  }

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
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}