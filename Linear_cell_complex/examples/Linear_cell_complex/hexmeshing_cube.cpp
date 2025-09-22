#include <CGAL/Hexmeshing_generate_two_refinement_mesh.h>
#include <CGAL/draw_linear_cell_complex.h>
#include <string>


void render_two_refinement(const std::string& file, int cube_cells_per_dim, int nb_levels = 1){
  auto lcc = CGAL::generate_two_refinement_mesh(CGAL::data_file_path("meshes/" + file), cube_cells_per_dim, nb_levels, true);

  CGAL::draw(lcc);
}


int main(){
  render_two_refinement("cube.off", 20, 1);
}
