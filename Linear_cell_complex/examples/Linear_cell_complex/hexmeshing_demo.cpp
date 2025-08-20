
#include <CGAL/Hexmeshing_for_linear_cell_complex_sequential.h>
#include <CGAL/Hexmeshing_mesh_data_for_hexmeshing.h>
#include <CGAL/Hexmeshing_render_results.h>

#include <CGAL/Combinatorial_map_save_load.h>
#include <CGAL/config.h>
#include <CGAL/Graphics_scene_options.h>
#include <string>
#include <cstdlib>
#include <filesystem>
#include <iostream>


void render_two_refinement(const std::string& file, int cube_cells_per_dim, int nb_levels = 1){
  CGAL::MeshDataForHexmeshing mesh;

  mesh.load_surface(CGAL::data_file_path("meshes/" + file));
  mesh.cubic_grid_from_aabb(cube_cells_per_dim);

  CGAL::HexMeshingData hdata;

  hdata.two_refinement(
    mesh,
    nb_levels,
    true
  );

  CGAL::render_two_refinement_result(hdata);
}


int main(){
  render_two_refinement("anc101.off", 20, 0);
  render_two_refinement("bunny.off", 20, 0);
  render_two_refinement("dragon_res3.off", 20, 0);

  render_two_refinement("anc101.off", 20, 1);
  render_two_refinement("bunny00.off", 20, 1);
  render_two_refinement("dragon_res3.off", 20, 1);
  // std::string s = "mesh2.off";
  // render_two_refinement(s, 20, 0);
  // render_two_refinement(s, 20, 1);

  render_two_refinement("anc101.off", 20, 2);
  render_two_refinement("bunny00.off", 20, 2);
  render_two_refinement("dragon_res2.off", 20, 2);
  
}
