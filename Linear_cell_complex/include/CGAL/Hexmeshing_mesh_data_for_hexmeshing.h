#ifndef HEXMESHING_MESH_DATA_FOR_HEXMESHING_H
#define HEXMESHING_MESH_DATA_FOR_HEXMESHING_H

#include <CGAL/hexmeshing/Hexmeshing_outer_alias.h>
#include <CGAL/hexmeshing/Hexmeshing_grid.h>

#include <CGAL/Combinatorial_map_save_load.h>
#include <CGAL/config.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <cstdlib>
#include <filesystem>

namespace CGAL {
  struct MeshDataForHexmeshing {
    MeshDataForHexmeshing() {}
    MeshDataForHexmeshing(Hexmeshing::Polyhedron& poly_out) : poly(poly_out) {
      // Triangulate before AABB
      CGAL::Polygon_mesh_processing::triangulate_faces(poly);
      // Compute AABB tree
      tree.insert(faces(poly).first, faces(poly).second, poly);
      tree.accelerate_distance_queries();
      tree.bbox();    
    }
    MeshDataForHexmeshing(Hexmeshing::Polyhedron poly_out, Hexmeshing::Grid grid_out) : poly(poly_out), grid(grid_out) {
      // Triangulate before AABB
      CGAL::Polygon_mesh_processing::triangulate_faces(poly);
      // Compute AABB tree
      tree.insert(faces(poly).first, faces(poly).second, poly);
      tree.accelerate_distance_queries();
      tree.bbox();    
    }

    void load_surface(const std::string& file) {
      std::ifstream off_file(file);
      CGAL_precondition_msg(off_file.good(), ("Input .off couldn't be read : " + file).c_str());
  
      off_file>>poly;
      
      // Triangulate before AABB
      CGAL::Polygon_mesh_processing::triangulate_faces(poly);
      // Compute AABB tree
      tree.insert(faces(poly).first, faces(poly).second, poly);
      tree.accelerate_distance_queries();
      tree.bbox();
    }
  
    void cubic_grid_from_aabb(int cube_cells_per_dim){
      assert(cube_cells_per_dim > 2);
      auto bbox = tree.bbox();
  
      Hexmeshing::Point center = {bbox.xmin() + (bbox.x_span()/2),
                      bbox.ymin() + (bbox.y_span()/2),
                      bbox.zmin() + (bbox.z_span()/2)};
  
      double max_size = std::max(std::max(bbox.x_span(), bbox.y_span()), bbox.z_span());
      grid = Hexmeshing::Grid::make_centered_cube(center, max_size / (cube_cells_per_dim-2), cube_cells_per_dim);
    }

    Hexmeshing::Grid* get_grid_pointer() {
      return &grid;
    }

    Hexmeshing::Tree* get_tree_pointer() {
      return &tree;
    }

  private:
    Hexmeshing::Polyhedron poly;
    Hexmeshing::Tree tree;
    Hexmeshing::Grid grid;
  };
}


#endif