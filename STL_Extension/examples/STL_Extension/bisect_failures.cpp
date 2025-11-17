#include <CGAL/bisect_failures.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/clip.h>
#include <CGAL/Surface_mesh/Surface_mesh.h>

#include <iostream>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point = Kernel::Point_3;
using Mesh = CGAL::Surface_mesh<Point>;


int main(int argc, char* argv[]) {
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/non_manifold.off");

  Mesh mesh_a;
  if(!CGAL::IO::read_polygon_mesh(filename, mesh_a)) {
    std::cerr << "Error: cannot read file " << filename << std::endl
    return EXIT_FAILURE;
  }

  std::cout << "Loaded mesh with " << mesh_a.number_of_vertices() << " vertices and "
            << mesh_a.number_of_faces() << " faces" << std::endl;

  const std::string clip_name = (argc > 2) ? argv[2] : CGAL::data_file_path("meshes/cube.off");
  Mesh mesh_b;
  if(!CGAL::IO::read_polygon_mesh(clip_name, mesh_b)) {
    std::cerr << "Error: cannot read file " << clip_name << std::endl
    return EXIT_FAILURE;
  }

  std::cout << "Loaded clipping mesh mesh with " << mesh_b.number_of_vertices() << " vertices and "
            << mesh_b.number_of_faces() << " faces" << std::endl

  //! [bisect_failures_snippet]
  // Define the callbacks for bisect_failures

  auto get_size = [](const Mesh& m) -> std::size_t {
    return m.number_of_faces();
  };

  auto simplify = [](Mesh& m, int start, int end) -> bool {
    for(auto i = end - 1; i >= start; --i) {
      const auto f = m.faces().begin() + i;
      CGAL::Euler::remove_face(halfedge(*f, m), m);
    }

    return m.is_valid();
  };

  auto run = [&mesh_b](Mesh& mesh) -> int {
    return CGAL::Polygon_mesh_processing::clip(mesh, mesh_b) ? EXIT_SUCCESS : EXIT_FAILURE;
  };

  auto save = [](const Mesh& m, const std::string& prefix) {
    std::string out_filename = prefix + ".off";
    if(!CGAL::IO::write_polygon_mesh(out_filename, m)) {
      std::cerr << "Warning: Could not save mesh to " << out_filename << std::endl
    } else {
      std::cout << "Saved mesh with " << m.number_of_faces()
                << " faces to " << out_filename << std::endl
    }
  };

  // Run bisection to find minimal failing case
  std::cout << "\n=== Starting bisection to find minimal failing case ===\n" << std::endl

  int result = CGAL::bisect_failures(mesh_a, get_size, simplify, run, save);
  //! [bisect_failures_snippet]

  if(result == EXIT_SUCCESS) {
    std::cout << "\nNo failure detected during bisection." << std::endl
  } else {
    std::cout << "\nFailure detected during bisection. Result code: " << result << std::endl
  }

  return EXIT_SUCCESS;
}
