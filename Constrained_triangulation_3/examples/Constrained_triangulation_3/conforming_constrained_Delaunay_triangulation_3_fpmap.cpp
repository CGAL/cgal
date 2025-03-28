#include <CGAL/make_conforming_constrained_Delaunay_triangulation_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/detect_features.h>

#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/IO/File_medit.h>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;

using Mesh = CGAL::Surface_mesh<K::Point_3>;
using face_descriptor = boost::graph_traits<Mesh>::face_descriptor;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char* argv[])
{
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/cross_quad.off");

  Mesh mesh;
  if(!PMP::IO::read_polygon_mesh(filename, mesh)) {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }

  using EIFMap = boost::property_map<Mesh, CGAL::edge_is_feature_t>::type;
  using FPMap = boost::property_map<Mesh, CGAL::face_patch_id_t<int>>::type;

  EIFMap eif = get(CGAL::edge_is_feature, mesh);
  FPMap fpmap = get(CGAL::face_patch_id_t<int>(), mesh);

  std::size_t number_of_patches = PMP::sharp_edges_segmentation(mesh, 80, eif, fpmap);

  std::cout << "Read " << mesh.number_of_vertices() << " vertices and "
                       << mesh.number_of_faces() << " faces" << std::endl;

  auto ccdt = CGAL::make_conforming_constrained_Delaunay_triangulation_3(mesh,
                CGAL::parameters::face_patch_map(fpmap));

  std::cout << "Number of vertices in the CDT: "
            << ccdt.triangulation().number_of_vertices() << '\n'
            << "Number of constrained facets in the CDT: "
            << ccdt.number_of_constrained_facets() << '\n';

  std::ofstream out("ccdt.mesh");
  CGAL::IO::write_MEDIT(out, ccdt.triangulation());
  out.close();

  return EXIT_SUCCESS;
}
