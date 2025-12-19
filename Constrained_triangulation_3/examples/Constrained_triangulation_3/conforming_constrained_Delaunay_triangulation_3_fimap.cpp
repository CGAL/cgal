#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/remesh_planar_patches.h>
#include <CGAL/make_conforming_constrained_Delaunay_triangulation_3.h>

#include <CGAL/IO/write_MEDIT.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;

using Mesh = CGAL::Surface_mesh<K::Point_3>;
using face_descriptor = boost::graph_traits<Mesh>::face_descriptor;


int main(int argc, char* argv[])
{
  std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/cross.off");

  Mesh input;
  if(!CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(filename, input)) {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }

  std::cout << "Read " << input.number_of_vertices() << " vertices and "
                       << input.number_of_faces() << " facets\n";

  Mesh mesh;
  auto plc_facet_map = get(CGAL::face_patch_id_t<int>(), mesh);

  // Remesh planar patches and segment the mesh into planar patches
  CGAL::Polygon_mesh_processing::remesh_planar_patches(input, mesh,
                             CGAL::parameters::default_values(),
                             CGAL::parameters::face_patch_map(plc_facet_map)
                                              .do_not_triangulate_faces(true));


  filename = argc > 2 ? argv[2] : "mesh.ply";
  CGAL::IO::write_polygon_mesh(filename, mesh,
      CGAL::parameters::stream_precision(17));
  std::cout << "Wrote segmented mesh to " << filename << "\n";

  // Create a conforming constrained Delaunay triangulation from the mesh
  auto ccdt = CGAL::make_conforming_constrained_Delaunay_triangulation_3(mesh,
                CGAL::parameters::plc_face_id(plc_facet_map));

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
