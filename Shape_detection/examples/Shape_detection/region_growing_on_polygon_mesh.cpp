// STL includes.
#include <string>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>

// CGAL includes.
#include <CGAL/memory.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Iterator_range.h>
#include <CGAL/HalfedgeDS_vector.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_polygon_mesh.h>

// Type declarations.
using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;

using FT      = typename Kernel::FT;
using Point_3 = typename Kernel::Point_3;

using Color = CGAL::Color;

// Choose the type of a container for a polygon mesh.
#define USE_SURFACE_MESH

#if defined(USE_SURFACE_MESH)

    using Polygon_mesh = CGAL::Surface_mesh<Point_3>;
    using Face_range   = typename Polygon_mesh::Face_range;

    using Neighbor_query = CGAL::Shape_detection::Polygon_mesh::One_ring_neighbor_query<Polygon_mesh>;
    using Region_type    = CGAL::Shape_detection::Polygon_mesh::Least_squares_plane_fit_region<Kernel, Polygon_mesh>;
    using Sorting        = CGAL::Shape_detection::Polygon_mesh::Least_squares_plane_fit_sorting<Kernel, Polygon_mesh, Neighbor_query>;

#else

    using Polygon_mesh = CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_3, CGAL::HalfedgeDS_vector>;
    using Face_range   = typename CGAL::Iterator_range<typename boost::graph_traits<Polygon_mesh>::face_iterator>;

    using Neighbor_query = CGAL::Shape_detection::Polygon_mesh::One_ring_neighbor_query<Polygon_mesh, Face_range>;
    using Region_type    = CGAL::Shape_detection::Polygon_mesh::Least_squares_plane_fit_region<Kernel, Polygon_mesh, Face_range>;
    using Sorting        = CGAL::Shape_detection::Polygon_mesh::Least_squares_plane_fit_sorting<Kernel, Polygon_mesh, Neighbor_query, Face_range>;

#endif

using Region  = std::vector<std::size_t>;
using Regions = std::vector<Region>;

using Vertex_to_point_map = typename Region_type::Vertex_to_point_map;

using Region_growing = CGAL::Shape_detection::Region_growing<Face_range, Neighbor_query, Region_type, typename Sorting::Seed_map>;

int main(int argc, char *argv[]) {

  std::cout << std::endl <<
    "region_growing_on_polygon_mesh example started"
  << std::endl << std::endl;

  // Load off data either from a local folder or a user-provided file.
  std::ifstream in(argc > 1 ? argv[1] : "data/polygon_mesh.off");
  CGAL::set_ascii_mode(in);

  if (!in) {
    std::cout <<
    "Error: cannot read the file polygon_mesh.off!" << std::endl;
    std::cout <<
    "You can either create a symlink to the data folder or provide this file by hand."
    << std::endl << std::endl;
    return EXIT_FAILURE;
  }

  Polygon_mesh polygon_mesh;
  in >> polygon_mesh;

  in.close();
  const Face_range face_range = faces(polygon_mesh);

  std::cout <<
    "* polygon mesh with "
  << face_range.size() <<
    " faces is loaded"
  << std::endl;

  // Default parameter values for the data file polygon_mesh.off.
  const FT          max_distance_to_plane = FT(1);
  const FT          max_accepted_angle    = FT(45);
  const std::size_t min_region_size       = 5;

  // Create instances of the classes Neighbor_query and Region_type.
  Neighbor_query neighbor_query(polygon_mesh);

  const Vertex_to_point_map vertex_to_point_map(
    get(CGAL::vertex_point, polygon_mesh));

  Region_type region_type(
    polygon_mesh,
    max_distance_to_plane, max_accepted_angle, min_region_size,
    vertex_to_point_map);

  // Sort face indices.
  Sorting sorting(
    polygon_mesh, neighbor_query,
    vertex_to_point_map);
  sorting.sort();

  // Create an instance of the region growing class.
  Region_growing region_growing(
    face_range, neighbor_query, region_type,
    sorting.seed_map());

  // Run the algorithm.
  Regions regions;
  region_growing.detect(std::back_inserter(regions));

  // Print the number of found regions.
  std::cout << "* " << regions.size() <<
    " regions have been found"
  << std::endl;

  // Save the result in a file only if it is stored in CGAL::Surface_mesh.
  #if defined(USE_SURFACE_MESH)

    using Face_index = typename Polygon_mesh::Face_index;

    // Save the result to a file in the user-provided path if any.
    srand(static_cast<unsigned int>(time(NULL)));
    if (argc > 2) {

      bool created;
      typename Polygon_mesh::template Property_map<Face_index, Color> face_color;
      boost::tie(face_color, created) =
        polygon_mesh.template add_property_map<Face_index, Color>(
          "f:color", Color(0, 0, 0));

      if (!created) {

        std::cout << std::endl <<
          "region_growing_on_polygon_mesh example finished"
        << std::endl << std::endl;

        return EXIT_FAILURE;
      }

      const std::string path     = argv[2];
      const std::string fullpath = path + "regions_polygon_mesh.off";

      std::ofstream out(fullpath);

      // Iterate through all regions.
      for (const auto& region : regions) {

        // Generate a random color.
        const Color color(
          static_cast<unsigned char>(rand() % 256),
          static_cast<unsigned char>(rand() % 256),
          static_cast<unsigned char>(rand() % 256));

        // Iterate through all region items.
        using size_type = typename Polygon_mesh::size_type;
        for (const auto index : region)
          face_color[Face_index(static_cast<size_type>(index))] = color;
      }

      out << polygon_mesh;
      out.close();

      std::cout <<
        "* polygon mesh is saved in "
      << fullpath << std::endl;
    }

  #endif

  std::cout << std::endl <<
    "region_growing_on_polygon_mesh example finished"
  << std::endl << std::endl;

  return EXIT_SUCCESS;
}
