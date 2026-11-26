#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/generate_hexahedral_mesh_using_two_refinement.h>
#include <CGAL/draw_linear_cell_complex.h>
#include <CGAL/Linear_cell_complex/IO/VTK.h>
#include <CGAL/Linear_cell_complex/lcc_scaled_jacobian.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Surface_mesh.h>
#include <string>

int main(int argc, char** argv)
{
  auto filename=(argc<2?CGAL::data_file_path("meshes/bunny00.off"):argv[1]);
  int cube_cells_per_dim=(argc<3?10:std::atoi(argv[2]));
  int nb_levels=(argc<4?2:std::atoi(argv[3]));
  bool trim=true;
  bool smooth=true;
  if((argc>4 && std::string("-notrim")==argv[4]) ||
     (argc>5 && std::string("-notrim")==argv[5])) { trim=false; }
  if((argc>4 && std::string("-nosmooth")==argv[4]) ||
     (argc>5 && std::string("-nosmooth")==argv[5])) { smooth=false; }

  CGAL::Surface_mesh<CGAL::Exact_predicates_inexact_constructions_kernel::Point_3> poly;
  CGAL::IO::read_polygon_mesh(filename, poly);
  CGAL::Polygon_mesh_processing::triangulate_faces(poly);

  auto lcc=CGAL::generate_hexahedral_mesh_using_two_refinement
    (poly, cube_cells_per_dim, nb_levels, CGAL::parameters::use_trimming(trim).use_smoothing(smooth));

  // Each hexahedron carries its scaled Jacobian: in the viewer, press D for a
  // palette, then Shift+D until "scaled Jacobian (per cell)".
  using LCC=decltype(lcc);
  using Dart=typename LCC::Dart_const_descriptor;
  CGAL::Graphics_scene_options<LCC, Dart, Dart, Dart, Dart> gso;
  gso.is_volume_valued=[](const LCC& alcc, Dart d)
  { return alcc.is_volume_combinatorial_hexahedron(d); };
  gso.volume_value=[](const LCC& alcc, Dart d)
  { return float(CGAL::scaled_jacobian_of_hexahedron(alcc, d)); };
  gso.volume_value_name="scaled Jacobian";

  // The faces of the hexahedra are oriented toward their inside: reverse the
  // normals so that the mesh is lit (N in the viewer swaps them back).
  CGAL::Graphics_scene scene;
  CGAL::add_to_graphics_scene(lcc, scene, gso);
  scene.reverse_all_normals();
  CGAL::draw_graphics_scene(scene, "Hexahedral mesh");
  CGAL::IO::write_VTK("hexmesh.vtk", lcc);

  return EXIT_SUCCESS;
}
