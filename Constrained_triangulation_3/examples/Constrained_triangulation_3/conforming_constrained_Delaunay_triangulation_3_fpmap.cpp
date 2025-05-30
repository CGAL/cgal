#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/write_MEDIT.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/detect_features.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/make_conforming_constrained_Delaunay_triangulation_3.h>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;

using Mesh = CGAL::Surface_mesh<K::Point_3>;
using face_descriptor = boost::graph_traits<Mesh>::face_descriptor;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char* argv[])
{
  std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/cross_quad.off");

  Mesh mesh;
  if(!PMP::IO::read_polygon_mesh(filename, mesh)) {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }

  std::cout << "Read " << mesh.number_of_vertices() << " vertices and "
                       << mesh.number_of_faces() << " facets\n";

  auto edge_is_feature_map = get(CGAL::edge_is_feature, mesh);
  auto face_patch_map = get(CGAL::face_patch_id_t<int>(), mesh);

  std::size_t number_of_patches = PMP::sharp_edges_segmentation(mesh, 10, edge_is_feature_map, face_patch_map);

  std::cout << "Number of patches: " << number_of_patches << std::endl;

  filename = argc > 2 ? argv[2] : "mesh.ply";
  CGAL::IO::write_polygon_mesh(filename, mesh, CGAL::parameters::stream_precision(17));
  std::cout << "Wrote segmented mesh to " << filename << "\n";

  auto ccdt = CGAL::make_conforming_constrained_Delaunay_triangulation_3(mesh,
                CGAL::parameters::plc_face_id(face_patch_map));

  std::cout << "Number of vertices in the CDT: "
            << ccdt.triangulation().number_of_vertices() << '\n'
            << "Number of constrained facets in the CDT: "
            << ccdt.number_of_constrained_facets() << '\n';

  filename = argc > 3 ? argv[3] : "out.mesh";
  std::ofstream out(filename);
  out.precision(17);
  CGAL::IO::write_MEDIT(out, ccdt, CGAL::parameters::with_plc_face_id(true));
  std::cout << "Wrote CDT to " << filename << "\n";

  return EXIT_SUCCESS;
}
