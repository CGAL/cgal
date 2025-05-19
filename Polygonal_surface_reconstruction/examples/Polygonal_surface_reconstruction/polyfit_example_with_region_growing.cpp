#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/read_points.h>
#include <CGAL/property_map.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Point_set.h>
#include <CGAL/Polygonal_surface_reconstruction.h>

#ifdef CGAL_USE_SCIP  // defined (or not) by CMake scripts, do not define by hand

#include <CGAL/SCIP_mixed_integer_program_traits.h>
typedef CGAL::SCIP_mixed_integer_program_traits<double> MIP_Solver;

#elif defined(CGAL_USE_GLPK)  // defined (or not) by CMake scripts, do not define by hand

#include <CGAL/GLPK_mixed_integer_program_traits.h>
typedef CGAL::GLPK_mixed_integer_program_traits<double>        MIP_Solver;

#endif

#if defined(CGAL_USE_GLPK) || defined(CGAL_USE_SCIP)

#include <fstream>
#include <CGAL/Timer.h>
#include <boost/range/irange.hpp>

typedef CGAL::Exact_predicates_inexact_constructions_kernel        Kernel;

typedef Kernel::FT       FT;
typedef Kernel::Point_3         Point;
typedef Kernel::Vector_3 Vector;

// Point with normal, and plane index.
typedef boost::tuple<Point, Vector, int> PNI;
typedef std::vector<PNI> Point_vector;

typedef CGAL::Nth_of_tuple_property_map<0, PNI>        Point_map;
typedef CGAL::Nth_of_tuple_property_map<1, PNI>        Normal_map;
typedef CGAL::Nth_of_tuple_property_map<2, PNI>        Plane_index_map;

using Point_map_region_growing = CGAL::Compose_property_map<CGAL::Random_access_property_map<Point_vector>, Point_map >;
using Normal_map_region_growing = CGAL::Compose_property_map<CGAL::Random_access_property_map<Point_vector>, Normal_map >;

using Region_type = CGAL::Shape_detection::Point_set::Least_squares_plane_fit_region<Kernel, std::size_t, Point_map_region_growing, Normal_map_region_growing>;
using Neighbor_query = CGAL::Shape_detection::Point_set::Sphere_neighbor_query<Kernel, std::size_t, Point_map_region_growing>;
using Region_growing = CGAL::Shape_detection::Region_growing<Neighbor_query, Region_type>;

typedef CGAL::Surface_mesh<Point>        Surface_mesh;
typedef CGAL::Polygonal_surface_reconstruction<Kernel> Polygonal_surface_reconstruction;

/*
* This example first extracts planes from the input point cloud
* (using region growing) and then reconstructs
* the surface model from the planes.
*/

int main(int argc, char* argv[])
{
  Point_vector points;

  // Load point set from a file.
  const std::string input_file = (argc > 1) ? argv[1] : CGAL::data_file_path("points_3/cube.pwn");
  std::ifstream input_stream(input_file.c_str());
  if (input_stream.fail()) {
    std::cerr << "Failed open file \'" << input_file << "\'" << std::endl;
    return EXIT_FAILURE;
  }
  input_stream.close();
  std::cout << "Loading point cloud: " << input_file << "...";

  CGAL::Timer t;
  t.start();
  if (!CGAL::IO::read_points(input_file.c_str(), std::back_inserter(points),
    CGAL::parameters::point_map(Point_map()).normal_map(Normal_map()))) {

    std::cerr << "Error: cannot read file " << input_file << std::endl;
    return EXIT_FAILURE;
  }
  else
    std::cout << " Done. " << points.size() << " points. Time: "
    << t.time() << " sec." << std::endl;

  //////////////////////////////////////////////////////////////////////////

  // Shape detection.

  // Default parameter values for the data file cube.pwn.
  const FT          search_sphere_radius = FT(2) / FT(100);
  const FT          max_distance_to_plane = FT(2) / FT(1000);
  const FT          max_accepted_angle = FT(25);
  const std::size_t min_region_size = 200;

  Point_map_region_growing point_map_rg(CGAL::make_random_access_property_map(points));
  Normal_map_region_growing normal_map_rg(CGAL::make_random_access_property_map(points));

  // Create instances of the classes Neighbor_query and Region_type.
  Neighbor_query neighbor_query(
    boost::irange<std::size_t>(0, points.size()), CGAL::parameters::sphere_radius(search_sphere_radius).point_map(point_map_rg));

  Region_type region_type(
    CGAL::parameters::
    maximum_distance(max_distance_to_plane).
    maximum_angle(max_accepted_angle).
    minimum_region_size(min_region_size).
    point_map(point_map_rg).
    normal_map(normal_map_rg));

  // Create an instance of the region growing class.
  Region_growing region_growing(
    boost::irange<std::size_t>(0, points.size()), neighbor_query, region_type);

  std::cout << "Extracting planes...";
  std::vector<typename Region_growing::Primitive_and_region> regions;
  t.reset();
  region_growing.detect(std::back_inserter(regions));
  std::cout << " Done. " << regions.size() << " planes extracted. Time: "
    << t.time() << " sec." << std::endl;

  // Stores the plane index of each point as the third element of the tuple.
  for (std::size_t i = 0; i < points.size(); ++i)
    // Uses the get function from the property map that accesses the 3rd element of the tuple.
    points[i].get<2>() = static_cast<int>(get(region_growing.region_map(), i));

  //////////////////////////////////////////////////////////////////////////

  // Reconstruction.

  std::cout << "Generating candidate faces...";
  t.reset();
  Polygonal_surface_reconstruction algo(
    points,
    Point_map(),
    Normal_map(),
    Plane_index_map()
  );
  std::cout << " Done. Time: " << t.time() << " sec." << std::endl;

  Surface_mesh model;
  std::cout << "Reconstructing...";
  t.reset();
  if (!algo.reconstruct<MIP_Solver>(model)) {
    std::cerr << "Failed: " << algo.error_message() << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << " Done. Time: " << t.time() << " sec." << std::endl;

  std::cout << "Saving...";
  t.reset();
  const std::string& output_file("with_region_growing_result.off");
  if (CGAL::IO::write_OFF(output_file, model))
    std::cout << " Done. Saved to " << output_file << ". Time: " << t.time() << " sec." << std::endl;
  else {
    std::cerr << " Failed saving file." << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

#else

int main(int, char**)
{
  std::cerr << "This test requires either GLPK or SCIP.\n";
  return EXIT_SUCCESS;
}

#endif  // defined(CGAL_USE_GLPK) || defined(CGAL_USE_SCIP)
