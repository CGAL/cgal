#define CGAL_MESH_3_VERBOSE 1
#define CGAL_CONCURRENT_MESH_3_VERBOSE 1
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/File_binary_mesh_3.h>
#include <CGAL/IO/output_to_vtu.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/make_mesh_3.h>
#include <fstream>
#include <iostream>

// Domain
using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Polyhedron = CGAL::Surface_mesh<K::Point_3>;
using Mesh_domain = CGAL::Polyhedral_mesh_domain_with_features_3<K, Polyhedron>;

// Triangulation
using Tr = CGAL::Mesh_triangulation_3<Mesh_domain, CGAL::Default,
                                      CGAL::Parallel_if_available_tag>::type;
using C3t3 = CGAL::Mesh_complex_3_in_triangulation_3<Tr>;

// Criteria
using Mesh_criteria = CGAL::Mesh_criteria_3<Tr>;

namespace params = CGAL::parameters;

void check_stream(const std::ios& stream,
                  const std::string& filename,
                  const std::string& operation,
                  bool ok = true) {
  if(!stream || !ok) {
    std::cerr << "Stream error after " << operation << " file: " << filename << std::endl;
    std::cerr << "Stream state: ";
    if(stream.rdstate() == std::ios::goodbit) {
      std::cerr << "no error";
    } else {
      if(stream.rdstate() & std::ios::eofbit) {
        std::cerr << "eofbit ";
      }
      if(stream.rdstate() & std::ios::failbit) {
        std::cerr << "failbit ";
      }
      if(stream.rdstate() & std::ios::badbit) {
        std::cerr << "badbit ";
      }
    }
    std::cerr << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

int main(int argc, char* argv[]) {
  const std::string fname = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/cube.off");

  // Create input polyhedron
  Polyhedron polyhedron;
  std::ifstream input(fname);
  check_stream(input, fname, "opening");
  input >> polyhedron;
  check_stream(input, fname, "reading polyhedron from");

  // Create domain
  Mesh_domain domain(polyhedron);
  domain.detect_features();

  // Mesh criteria (no cell_size set)
  Mesh_criteria criteria(params::facet_angle(25).facet_size(0.15).facet_distance(0.05).cell_radius_edge_ratio(3));

  // Mesh generation
  const C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, params::no_perturb().no_exude());

  const auto nb_vertices = c3t3.triangulation().number_of_vertices();
  const auto nb_far_vertices = c3t3.number_of_far_points();
  const auto nb_cells = c3t3.number_of_cells();

  std::cout << "Created a c3t3 with " << c3t3.triangulation().number_of_vertices() << " vertices, and "
            << c3t3.number_of_cells() << " cells" << std::endl;

  // Output
  {
    const std::string filename = "ascii.vtu";
    std::ofstream out(filename);
    check_stream(out, filename, "opening");
    CGAL::IO::output_to_vtu(out, c3t3, CGAL::IO::ASCII);
    check_stream(out, filename, "writing (ASCII)");
  }

  {
    const std::string filename = "binary.vtu";
    std::ofstream out(filename, std::ios::binary);
    check_stream(out, filename, "opening");
    CGAL::IO::output_to_vtu(out, c3t3, CGAL::IO::BINARY);
    check_stream(out, filename, "writing (BINARY)");
  }

  const std::string filename = "mesh.binary.cgal";
  {
    std::ofstream out(filename, std::ios::binary);
    check_stream(out, filename, "opening");
    bool ok = CGAL::IO::save_binary_file(out, c3t3);
    check_stream(out, filename, "writing (binary)", ok);
  }

  // Input
  C3t3 bis;
  {
    std::ifstream in(filename, std::ios::binary);
    check_stream(in, filename, "opening (binary)");
    bool ok = CGAL::IO::load_binary_file(in, bis);
    check_stream(in, filename, "reading binary file", ok);
  }

  {
    const std::string bis_filename = "bis.vtu";
    std::ofstream bis_os(bis_filename);
    check_stream(bis_os, bis_filename, "opening (bis.vtu)");
    CGAL::IO::output_to_vtu(bis_os, bis, CGAL::IO::ASCII);
    check_stream(bis_os, bis_filename, "writing (ASCII)");
  }

  assert(bis.number_of_cells() == nb_cells);
  assert(bis.triangulation().number_of_vertices() == nb_vertices);
  assert(bis.number_of_far_points() == nb_far_vertices);

  std::cout << "Mesh generation and output completed successfully." << std::endl;

  return EXIT_SUCCESS;
}
