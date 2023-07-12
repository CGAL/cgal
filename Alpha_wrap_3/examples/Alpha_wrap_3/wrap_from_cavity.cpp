#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/alpha_wrap_3.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

#include <iostream>
#include <string>

namespace PMP = CGAL::Polygon_mesh_processing;

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = K::Point_3;

using Mesh = CGAL::Surface_mesh<Point_3>;

int main(int argc, char** argv)
{
  // Read the input
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/armadillo.off");
  std::cout << "Reading " << filename << "..." << std::endl;

  Mesh input;
  if(!PMP::IO::read_polygon_mesh(filename, input) ||
     is_empty(input) || !is_triangle_mesh(input))
  {
    std::cerr << "Invalid input." << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Input: " << num_vertices(input) << " vertices, " << num_faces(input) << " faces" << std::endl;

  const double relative_alpha = (argc > 2) ? std::stod(argv[2]) : 30.;
  const double relative_offset = (argc > 3) ? std::stod(argv[3]) : 600.;

  // Compute the alpha and offset values
  CGAL::Bbox_3 bbox = CGAL::Polygon_mesh_processing::bbox(input);
  const double diag_length = std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                       CGAL::square(bbox.ymax() - bbox.ymin()) +
                                       CGAL::square(bbox.zmax() - bbox.zmin()));

  const double alpha = diag_length / relative_alpha;
  const double offset = diag_length / relative_offset;

  // Construct the wrap
  CGAL::Real_timer t;
  t.start();

  // There is no limit on how many seeds can be used.
  // However, the algorithm automatically determines whether a seed can be used
  // to initialize the refinement based on a few conditions (distance to the offset, value of alpha, etc.)
  // See internal function Alpha_wrap_3::initialize_from_cavities() for more information
  std::vector<Point_3> seeds =
  {
    Point_3(0, 50, 0) // a point within the armadillo surface
  };

  Mesh wrap;
  alpha_wrap_3(input, alpha, offset, wrap, CGAL::parameters::seed_points(std::ref(seeds)));

  t.stop();
  std::cout << "Result: " << num_vertices(wrap) << " vertices, " << num_faces(wrap) << " faces" << std::endl;
  std::cout << "Took " << t.time() << " s." << std::endl;

  // Save the result
  std::string input_name = std::string(filename);
  input_name = input_name.substr(input_name.find_last_of("/") + 1, input_name.length() - 1);
  input_name = input_name.substr(0, input_name.find_last_of("."));
  std::string output_name = input_name + "_cavity_" + std::to_string(static_cast<int>(relative_alpha))
                            + "_" + std::to_string(static_cast<int>(relative_offset)) + ".off";
  std::cout << "Writing to " << output_name << std::endl;
  CGAL::IO::write_polygon_mesh(output_name, wrap, CGAL::parameters::stream_precision(17));

  return EXIT_SUCCESS;
}
