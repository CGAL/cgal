#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/generate_hexahedral_mesh_using_two_refinement.h>
#include <CGAL/draw_linear_cell_complex.h>
#include <CGAL/Linear_cell_complex/IO/VTK.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Surface_mesh.h>
#include <string>

int main(int argc, char** argv)
{
  std::string filename=(argc<2?CGAL::data_file_path("meshes/bunny00.off"):argv[1]);

  CGAL::Surface_mesh<CGAL::Exact_predicates_inexact_constructions_kernel::Point_3> poly;
  CGAL::IO::read_polygon_mesh(filename, poly);
  CGAL::Polygon_mesh_processing::triangulate_faces(poly);

  auto lcc=CGAL::generate_hexahedral_mesh_using_two_refinement
    (poly, 18, 2, CGAL::parameters::use_trimming(true).use_smoothing(true));

  CGAL::draw(lcc);
  CGAL::IO::write_VTK("hexmesh.vtk", lcc);

  return EXIT_SUCCESS;
}
