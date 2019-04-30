#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/PLY_reader.h>
#include <CGAL/IO/PLY_writer.h>

#define CGAL_KSR_VERBOSE_LEVEL 3
#include <CGAL/Kinetic_shape_reconstruction_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point_3;
typedef std::vector<std::size_t> Polygon;

typedef CGAL::Kinetic_shape_reconstruction_3<Kernel> Reconstruction;

struct My_polygon_map
{
  typedef std::vector<std::size_t> key_type;
  typedef std::vector<Point_3> value_type;
  typedef value_type reference;
  typedef boost::readable_property_map_tag category;
  
  const std::vector<Point_3>* points;

  My_polygon_map (const std::vector<Point_3>& points) : points (&points) { }

  friend reference get (const My_polygon_map& map, const key_type& k)
  {
    reference out;
    out.reserve (k.size());
    std::transform (k.begin(), k.end(), std::back_inserter (out),
                    [&](const std::size_t& idx) -> Point_3 { return (*(map.points))[idx]; });
    return out;
  }
};

int main (int argc, char** argv)
{
  std::string input_shapes_filename = (argc > 1 ? argv[1] : "data/simple_planes.ply");
  std::ifstream input_shapes_file (input_shapes_filename);

  std::vector<Point_3> vertices;
  std::vector<Polygon> facets;

  if (!CGAL::read_PLY (input_shapes_file, vertices, facets))
  {
    std::cerr << "Error: can't read " << input_shapes_filename << std::endl;
    return EXIT_FAILURE;
  }

  Reconstruction reconstruction;

  reconstruction.partition (facets, My_polygon_map (vertices));

  vertices.clear();
  facets.clear();

  reconstruction.output_partition_facets_to_polygon_soup (std::back_inserter (vertices),
                                                          std::back_inserter (facets));

  std::ofstream output_shapes_file ("out.ply");
//  CGAL::set_binary_mode (output_shapes_file);
  CGAL::write_PLY (output_shapes_file, vertices, facets);

  return EXIT_SUCCESS;
}
