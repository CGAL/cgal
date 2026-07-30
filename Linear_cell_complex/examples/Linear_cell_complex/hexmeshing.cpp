#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Hexmeshing_generate_two_refinement_mesh.h>
#include <CGAL/draw_linear_cell_complex.h>
#include <CGAL/Linear_cell_complex/IO/VTK.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polyhedron_3.h>
#include <string>

int main(int argc, char** argv)
{
  std::string filename=(argc<2?CGAL::data_file_path("meshes/bunny00.off"):argv[1]);
  int cube_cells_per_dim=(argc<3?18:std::atoi(argv[2]));
  int nb_levels=(argc<4?2:std::atoi(argv[3]));
  bool trim=(argc<5||std::string("-notrim")!=argv[4])?true:false;

  std::ifstream off_file(filename);
  if(!off_file.good())
  {
    std::cout<<"Input mesh couldn't be read: "<<filename<<std::endl;
    return EXIT_FAILURE;
  }
  
  CGAL::Polyhedron_3<CGAL::Exact_predicates_inexact_constructions_kernel> poly;
  CGAL::IO::read_polygon_mesh(filename, poly);
  CGAL::Polygon_mesh_processing::triangulate_faces(poly);
  
  auto lcc=CGAL::generate_hexahedral_mesh_using_two_refinement
    (poly, cube_cells_per_dim, nb_levels, trim);

  CGAL::draw(lcc);
  CGAL::IO::write_VTK("hexmesh.vtk", lcc);

  return EXIT_SUCCESS;
}
