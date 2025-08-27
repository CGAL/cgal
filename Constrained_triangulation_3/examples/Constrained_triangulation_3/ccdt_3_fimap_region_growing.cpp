#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/write_MEDIT.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/region_growing.h>
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

  auto [face_patch_map, _] =mesh.add_property_map<face_descriptor, std::size_t>("f:patch_id");

  const auto bbox = CGAL::Polygon_mesh_processing::bbox(mesh);
  const double bbox_max_span = (std::max)({bbox.x_span(), bbox.y_span(), bbox.z_span()});

  std::cout << "Merging facets into coplanar patches..." << std::endl;

  auto number_of_patches = CGAL::Polygon_mesh_processing::region_growing_of_planes_on_faces(
      mesh,
      face_patch_map,
      CGAL::parameters::maximum_distance(bbox_max_span * 1.e-6)
                       .maximum_angle(5.));

  for(auto f: faces(mesh))
  {
    // if region growing did not assign a patch id, assign one
    if(get(face_patch_map, f) == static_cast<std::size_t>(-1)) {
      put(face_patch_map, f, number_of_patches++);
    }
  }

  std::cout << "Number of patches: " << number_of_patches << std::endl;

  filename = argc > 2 ? argv[2] : "mesh.ply";
  CGAL::IO::write_polygon_mesh(filename, mesh,
                               CGAL::parameters::stream_precision(17)
                                                .use_binary_mode(false)
                                                .face_patch_map(face_patch_map));
  std::cout << "-- Wrote segmented mesh to \"" << filename << "\"\n";

  std::cout << "Creating a conforming constrained Delaunay triangulation...\n";
  auto ccdt = CGAL::make_conforming_constrained_Delaunay_triangulation_3(mesh,
                CGAL::parameters::plc_face_id(face_patch_map));

  std::cout << "Number of vertices in the CDT: "
            << ccdt.triangulation().number_of_vertices() << '\n'
            << "Number of constrained facets in the CDT: "
            << ccdt.number_of_constrained_facets() << '\n';

  // Write the CDT to a file, with the PLC face ids
  filename = argc > 3 ? argv[3] : "out.mesh";
  std::ofstream out(filename);
  out.precision(17);
  CGAL::IO::write_MEDIT(out, ccdt, CGAL::parameters::with_plc_face_id(true));
  std::cout << "-- Wrote CDT to \"" << filename << "\"\n";

  return EXIT_SUCCESS;
}
