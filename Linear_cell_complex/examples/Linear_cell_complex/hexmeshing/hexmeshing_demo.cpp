
#include "hexmeshing.h"
#include "utils.h"

#include <CGAL/Combinatorial_map_save_load.h>
#include <CGAL/config.h>
#include <CGAL/draw_polyhedron.h>
#include <CGAL/Graphics_scene_options.h>
#include <cstdlib>
#include <filesystem>
#include <iostream>

using namespace CGAL::HexRefinement;

Polyhedron load_surface(const std::string& file) {
  std::ifstream off_file(file);
  CGAL_precondition_msg(off_file.good(), ("Input .off couldn't be read : " + file).c_str());

  Polyhedron surface;
  off_file>>surface;
  return surface;
}

Tree get_surface_aabb(Polyhedron& poly) {
  // Triangulate before AABB
  Tree tree;
  CGAL::Polygon_mesh_processing::triangulate_faces(poly);
  // Compute AABB tree
  tree.insert(faces(poly).first, faces(poly).second, poly);
  tree.accelerate_distance_queries();
  tree.bbox();

  return tree;
}

TwoRefinement::Grid cubic_grid_from_aabb(Tree& aabb, int cube_cells_per_dim){
  auto bbox = aabb.bbox();

  Point center = {bbox.xmin() + (bbox.x_span()/2),
                  bbox.ymin() + (bbox.y_span()/2),
                  bbox.zmin() + (bbox.z_span()/2)};

  double max_size = std::max(std::max(bbox.x_span(), bbox.y_span()), bbox.z_span());
  return TwoRefinement::Grid::make_centered_cube(center, max_size / cube_cells_per_dim, cube_cells_per_dim);
}

void render_two_refinement(const std::string& file, int cube_cells_per_dim, int nb_levels = 1){
  Polyhedron poly = load_surface(CGAL::data_file_path("hexmeshing/mesh/" + file));
  Tree aabb = get_surface_aabb(poly);
  TwoRefinement::Grid grid = cubic_grid_from_aabb(aabb, cube_cells_per_dim);

  LCC lcc = two_refinement(
    grid,
    TwoRefinement::is_volume_intersecting_poly(aabb),
    TwoRefinement::is_volume_intersecting_poly(aabb),
    // 2,
    nb_levels
  );

  render_two_refinement_result(lcc, aabb, false);
  render_two_refinement_result(lcc, aabb);
}


int main(){
  render_two_refinement("mesh1.off", 20, 1);
  render_two_refinement("mesh1.off", 20, 2);
  render_two_refinement("mesh2.off", 20, 3);
}
