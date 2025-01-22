#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include "algo/3d/MeshOffset.h"

#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>

#include <iostream>
#include <fstream>

namespace PMP = ::CGAL::Polygon_mesh_processing;

using K = CGAL::Exact_predicates_exact_constructions_kernel;
using FT = K::FT;
using Point = K::Point_3;
using Vector = K::Vector_3;

using Mesh = CGAL::Surface_mesh<Point>;
using face_descriptor = boost::graph_traits<Mesh>::face_descriptor;

int main(int argc, char** argv)
{
  if(argc < 2)
  {
    std::cerr << "Usage: "
              << argv[0] << "\n"
              << "\tinput_filename\n"
              << "\t[output_filename] (PLY)\n"
              << "\t[weights_filename.txt]\n"
              << std::endl;
    std::cerr << "Input format: any format readable with CGAL::IO::read_polygon_mesh()" << std::endl;
    std::cerr << "Output format: PLY" << std::endl;
    std::cerr << "Weight format:" << std::endl;
    std::cerr << "  x1: val" << std::endl;
    std::cerr << "  x2: val" << std::endl;
    std::cerr << "  y1: val" << std::endl;
    std::cerr << "  y2: val" << std::endl;
    std::cerr << "  bottom: val (optional line)" << std::endl;
    std::cerr << "  top: val (optional line)" << std::endl;
    return EXIT_FAILURE;
  }

  const char* input_filename = argv[1];
  const char* output_filename = (argc > 2) ? argv[2] : "out.ply";
  const char* weights_filename = (argc > 3) ? argv[3] : nullptr;

  std::cout << "in: " << input_filename << std::endl;
  std::cout << "out: " << output_filename << std::endl;
  std::cout << "weight: " << weights_filename << std::endl;

  Mesh sm;
  if(!CGAL::IO::read_polygon_mesh(input_filename, sm)) {
    std::cerr << "Error: failed to read input " << input_filename << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Input mesh: " << num_vertices(sm) << " NV " << num_faces(sm) << " NF" << std::endl;

  if(!algo::_3d::MeshOffset::assign_weights(sm, weights_filename))
    return EXIT_FAILURE;

  std::ofstream out(output_filename);
  if(!out) {
    std::cerr << "Error: failed to create output file" << std::endl;
    return EXIT_FAILURE;
  }

  CGAL::IO::write_PLY(out, sm,
                      CGAL::parameters::use_binary_mode(false)
                                       .stream_precision(17));


  return EXIT_SUCCESS;
}
