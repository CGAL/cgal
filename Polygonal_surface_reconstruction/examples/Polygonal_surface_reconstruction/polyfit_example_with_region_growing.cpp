#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/IO/Writer_OFF.h>
#include <CGAL/property_map.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_point_set.h>
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

typedef CGAL::Shape_detection::Point_set::
Sphere_neighbor_query<Kernel, Point_vector, Point_map> Neighbor_query;
typedef CGAL::Shape_detection::Point_set::
Least_squares_plane_fit_region<Kernel, Point_vector, Point_map, Normal_map> Region_type;
typedef CGAL::Shape_detection::
Region_growing<Point_vector, Neighbor_query, Region_type> Region_growing;

typedef CGAL::Surface_mesh<Point>        Surface_mesh;
typedef        CGAL::Polygonal_surface_reconstruction<Kernel> Polygonal_surface_reconstruction;

class Index_map {

public:
  using key_type = std::size_t;
  using value_type = int;
  using reference = value_type;
  using category = boost::readable_property_map_tag;

  Index_map() { }
  template<typename PointRange>
  Index_map(
    const PointRange& points,
    const std::vector< std::vector<std::size_t> >& regions) :
  m_indices(new std::vector<int>(points.size(), -1)) {

    for (std::size_t i = 0; i < regions.size(); ++i)
      for (const std::size_t idx : regions[i])
        (*m_indices)[idx] = static_cast<int>(i);
  }

  inline friend value_type get(
    const Index_map& index_map,
    const key_type key) {

    const auto& indices = *(index_map.m_indices);
    return indices[key];
  }

private:
  std::shared_ptr< std::vector<int> > m_indices;
};

/*
* This example first extracts planes from the input point cloud
* (using region growing) and then reconstructs
* the surface model from the planes.
*/

int main()
{
  Point_vector points;

  // Load point set from a file.
  const std::string input_file("data/cube.pwn");
    std::ifstream input_stream(input_file.c_str());
  if (input_stream.fail()) {
    std::cerr << "Failed open file \'" << input_file << "\'" << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "Loading point cloud: " << input_file << "...";

  CGAL::Timer t;
  t.start();
    if (!input_stream ||
    !CGAL::read_xyz_points(input_stream,
      std::back_inserter(points),
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
  const FT          search_sphere_radius  = FT(2) / FT(100);
  const FT          max_distance_to_plane = FT(2) / FT(1000);
  const FT          max_accepted_angle    = FT(25);
  const std::size_t min_region_size       = 200;

  // Create instances of the classes Neighbor_query and Region_type.
  Neighbor_query neighbor_query(
    points,
    search_sphere_radius);

  Region_type region_type(
    points,
    max_distance_to_plane, max_accepted_angle, min_region_size);

  // Create an instance of the region growing class.
  Region_growing region_growing(
    points, neighbor_query, region_type);

  std::cout << "Extracting planes...";
  std::vector< std::vector<std::size_t> > regions;
  t.reset();
  region_growing.detect(std::back_inserter(regions));
  std::cout << " Done. " << regions.size() << " planes extracted. Time: "
  << t.time() << " sec." << std::endl;

  // Stores the plane index of each point as the third element of the tuple.
  Index_map index_map(points, regions);
  for (std::size_t i = 0; i < points.size(); ++i) {
    // Uses the get function from the property map that accesses the 3rd element of the tuple.
    const int plane_index = get(index_map, i);
    points[i].get<2>() = plane_index;
  }

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
  const std::string& output_file("data/cube_result.off");
  std::ofstream output_stream(output_file.c_str());
  if (output_stream && CGAL::write_off(output_stream, model)) {
    // flush the buffer
    output_stream << std::flush;
    std::cout << " Done. Saved to " << output_file << ". Time: " << t.time() << " sec." << std::endl;
  }
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
