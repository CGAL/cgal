#include <CGAL/Hexmeshing_generate_two_refinement_mesh.h>
#include <CGAL/draw_linear_cell_complex.h>
#include <string>

int main(int argc, char** argv)
{
  std::string filename=(argc<2?CGAL::data_file_path("meshes/cube.off"):argv[1]);

  auto lcc=CGAL::generate_two_refinement_mesh(filename, 20, 2, true);
  CGAL::draw(lcc);

  return EXIT_SUCCESS;
}
