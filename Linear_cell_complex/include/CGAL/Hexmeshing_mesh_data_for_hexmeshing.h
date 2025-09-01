// Copyright (c) 2025 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
// Contributor(s): Soichiro Yamazaki <soichiro19998@gmail.com>
//
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
  class Mesh_data_for_hexmeshing {
  public:
    Mesh_data_for_hexmeshing() {}
    Mesh_data_for_hexmeshing(internal::Hexmeshing::Polyhedron& poly_out) : poly(poly_out) {
      construct_tree_from_poly();
    }
    Mesh_data_for_hexmeshing(internal::Hexmeshing::Polyhedron poly_out, internal::Hexmeshing::Grid grid_out) : poly(poly_out), grid(grid_out) {
      construct_tree_from_poly();
    }

    void load_surface(const std::string& file) {
      std::ifstream off_file(file);
      CGAL_precondition_msg(off_file.good(), ("Input .off couldn't be read : " + file).c_str());

      off_file>>poly;

      construct_tree_from_poly();
    }

    void cubic_grid_from_aabb(int cube_cells_per_dim){
      assert(cube_cells_per_dim > 2);
      auto bbox = tree.bbox();

      internal::Hexmeshing::Point center = {bbox.xmin() + (bbox.x_span()/2),
                      bbox.ymin() + (bbox.y_span()/2),
                      bbox.zmin() + (bbox.z_span()/2)};

      double max_size = std::max(std::max(bbox.x_span(), bbox.y_span()), bbox.z_span());
      grid = internal::Hexmeshing::Grid::make_centered_cube(center, max_size / (cube_cells_per_dim-2), cube_cells_per_dim);
    }

    internal::Hexmeshing::Grid* get_grid_pointer() {
      return &grid;
    }

    internal::Hexmeshing::Tree* get_tree_pointer() {
      return &tree;
    }

  private:
    void construct_tree_from_poly() {
      // Triangulate before AABB
      CGAL::Polygon_mesh_processing::triangulate_faces(poly);
      // Compute AABB tree
      tree.insert(faces(poly).first, faces(poly).second, poly);
      tree.accelerate_distance_queries();
      tree.bbox();
    }
    internal::Hexmeshing::Polyhedron poly;
    internal::Hexmeshing::Tree tree;
    internal::Hexmeshing::Grid grid;
  };
}


#endif