#include <CGAL/Hexmeshing_generate_two_refinement_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/config.h>
#include <string>


void render_two_refinement(const std::string& file, int cube_cells_per_dim, int nb_levels = 1)
{
  CGAL::Polyhedron_3<CGAL::Exact_predicates_inexact_constructions_kernel> poly;
  CGAL::IO::read_polygon_mesh(CGAL::data_file_path("meshes/" + file), poly);

  CGAL::generate_hexahedral_mesh_using_two_refinement
    (poly, cube_cells_per_dim, nb_levels, true);

  // CGAL::render_two_refinement_result(hdata);
  std::cout << "file: " << file << std::endl;
  std::cout << "nb_levels: " << nb_levels << std::endl;
  std::cout << "test passed\n" << std::endl;
}

int main()
{
  // render_two_refinement("anc101.off", 20, 0);
  render_two_refinement("bunny00.off", 16, 0);
  render_two_refinement("dragon_res3.off", 16, 0);

  // render_two_refinement("anc101.off", 20, 1);
  render_two_refinement("bunny00.off", 16, 1);
  render_two_refinement("dragon_res3.off", 16, 1);
  // std::string s = "mesh2.off";
  // render_two_refinement(s, 20, 0);
  // render_two_refinement(s, 20, 1);

  // render_two_refinement("anc101.off", 20, 2);
  render_two_refinement("bunny00.off", 16, 2);
  render_two_refinement("dragon_res2.off", 16, 2);

  return EXIT_SUCCESS;
}
