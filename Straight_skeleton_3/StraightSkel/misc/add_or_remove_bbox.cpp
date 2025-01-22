#include "algo/3d/MeshOffset.h"

// For now, there are issues when EPECK is embedded back into doubles,
// so dump files with EPECK and read them again with EPECK
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Bbox_3.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/generators.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

#include <boost/property_map/function_property_map.hpp>

#include <iostream>
#include <fstream>
#include <type_traits>
#include <unordered_map>

using EPECK = CGAL::Exact_predicates_exact_constructions_kernel;
using K = EPECK;
using FT = K::FT;
using Point = K::Point_3;
using Vector = K::Vector_3;
using Iso_cuboid = K::Iso_cuboid_3;

using Mesh = CGAL::Surface_mesh<Point>;
using vertex_descriptor = boost::graph_traits<Mesh>::vertex_descriptor;
using face_descriptor = boost::graph_traits<Mesh>::face_descriptor;

namespace PMP = ::CGAL::Polygon_mesh_processing;

int main(int argc, char** argv)
{
  if(argc < 2)
  {
    std::cerr << "Usage: "
              << argv[0] << "\n"
              << "\tinput_filename\n"
              << "\t[invert]\n"
              << "\t[output_filename]\n"
              << std::endl;
    std::cerr << "Input format: any format readable with CGAL::IO::read_polygon_mesh()" << std::endl;
    std::cerr << "'invert' format: 'true'/'false'" << std::endl;
    std::cerr << "Output format: PLY (adding bbox); OBJ (removing bbox)" << std::endl;
    return EXIT_FAILURE;
  }

  const char* input_filename = argv[1];
  const std::string do_add_str = (argc > 2) ? argv[2] : "add";
  const char* output_filename = (argc > 3) ? argv[3] : "out.ply";
  const bool do_add = (do_add_str == "add");

  std::cout << "in: " << input_filename << std::endl;
  std::cout << "add: " << do_add << std::endl;
  std::cout << "out: " << output_filename << std::endl;

  Mesh sm;

  // @tmp, while https://github.com/CGAL/cgal/pull/7874 is not integrated
  // because we don't want to lose the f:weight property map while reading/writing PLY files
  const std::string ext = CGAL::IO::internal::get_file_extension(input_filename);
  if(ext == "ply") {
    std::ifstream in(input_filename);
    if(!in || !CGAL::IO::read_PLY(in, sm)) {
      std::cerr << "Error: failed to read input (PLY)" << std::endl;
      return EXIT_FAILURE;
    }
  } else {
    CGAL_assertion(!do_add); // expecting to be here for removal only
    // PMP in case the mesh is non-manifold
    if(!PMP::IO::read_polygon_mesh(input_filename, sm, CGAL::parameters::verbose(true))) {
      std::cerr << "Error: failed to read input" << std::endl;
      return EXIT_FAILURE;
    }
  }

  std::cout << "Input mesh: " << num_vertices(sm) << " NV " << num_faces(sm) << " NF" << std::endl;

  if(CGAL::is_empty(sm) || !CGAL::is_closed(sm)) {
    std::cerr << "Error: empty or open input" << std::endl;
    return EXIT_FAILURE;
  }

  if(do_add) {
    if(!algo::_3d::MeshOffset::invert_and_add_bbox(sm)) {
        return EXIT_FAILURE;
    }

    std::ofstream out(output_filename);
    if(!out) {
      std::cerr << "Error: failed to create output file: " << output_filename << std::endl;
      return EXIT_FAILURE;
    }

    if(!CGAL::IO::write_PLY(out, sm,
                            CGAL::parameters::use_binary_mode(false)
                                             .stream_precision(17))) {
      std::cerr << "Error: failed to write in " << output_filename << std::endl;
      return EXIT_FAILURE;
    }
  } else {
    if (!algo::_3d::MeshOffset::remove_bbox_and_invert(sm)) {
      return EXIT_FAILURE;
    }

    std::ofstream out(output_filename);
    if(!out) {
      std::cerr << "Error: failed to create output file: " << output_filename << std::endl;
      return EXIT_FAILURE;
    }

    if constexpr (std::is_same<K, EPECK>::value) {
#if 0 // this writes 'doubles' so not so generic, huh...
      auto ef = [&sm](vertex_descriptor v) { return get(CGAL::vertex_point, sm, v).exact(); };
      auto evpm = boost::make_function_property_map<vertex_descriptor>(ef);
      if(!CGAL::IO::write_OBJ(out, sm, CGAL::parameters::vertex_point_map(evpm))) {
        std::cerr << "Error: failed to write in " << output_filename << std::endl;
        return EXIT_FAILURE;
      }
#else
      sm.collect_garbage();

      std::ofstream out(output_filename);
      if(!out) {
        std::cerr << "Error: failed to create output file: " << output_filename << std::endl;
        return EXIT_FAILURE;
      }

      for(vertex_descriptor v : vertices(sm)) {
        const Point& p = get(CGAL::vertex_point, sm, v);
        out << "v " << p.x().exact() << " " << p.y().exact() << " " << p.z().exact() << std::endl;
      }

      for(face_descriptor f : faces(sm)) {
        out << "f";
        for(vertex_descriptor v : vertices_around_face(halfedge(f, sm), sm)) {
          out << " " << v.id() + 1;
        }
        out << "\n";
      }
#endif
    } else {
      if(!CGAL::IO::write_OBJ(out, sm, CGAL::parameters::stream_precision(17))) {
        std::cerr << "Error: failed to write in " << output_filename << std::endl;
        return EXIT_FAILURE;
      }
    }
  }

  std::cout << "Done!" << std::endl;
}