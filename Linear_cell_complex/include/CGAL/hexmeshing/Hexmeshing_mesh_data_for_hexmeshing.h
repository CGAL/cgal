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
#ifndef CGAL_HEXMESHING_MESH_DATA_FOR_HEXMESHING_H
#define CGAL_HEXMESHING_MESH_DATA_FOR_HEXMESHING_H

#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/hexmeshing/Hexmeshing_grid.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <cstdlib>

namespace CGAL::internal
{
  template<typename TriangleMesh>
  class Mesh_data_for_hexmeshing
  {
  public:
    using Point=typename TriangleMesh::Point;
    using Kernel=typename CGAL::Kernel_traits<Point>::Kernel;
    using FT=typename Kernel::FT;
    using Vector=typename Kernel::Vector_3;
    using Triangle=typename Kernel::Triangle_3;
    using Segment=typename Kernel::Segment_3;

    using Primitive=CGAL::AABB_face_graph_triangle_primitive<TriangleMesh>;
    using AABB_Traits=CGAL::AABB_traits_3<Kernel, Primitive>;
    using Tree=CGAL::AABB_tree<AABB_Traits>;
    using Primitive_id=typename Tree::Primitive_id;
    using Side_of_mesh=CGAL::Side_of_triangle_mesh<TriangleMesh, Kernel>;

    Mesh_data_for_hexmeshing(const TriangleMesh& poly_out, int cube_cells_per_dim) :
        poly(poly_out)
    {
      construct_tree_from_poly();
      cubic_grid_from_aabb(cube_cells_per_dim);
    }
    Mesh_data_for_hexmeshing(TriangleMesh poly_out, Hexmeshing::Grid grid_out) :
        poly(poly_out), grid(grid_out)
    {
      construct_tree_from_poly();
    }

    Hexmeshing::Grid* get_grid_pointer() {
      return &grid;
    }

    Tree* get_tree_pointer()
    { return &tree; }

  private:
    void construct_tree_from_poly()
    {
      // Compute AABB tree
      tree.insert(faces(poly).first, faces(poly).second, poly);
      tree.accelerate_distance_queries();
      tree.bbox();
    }

    void cubic_grid_from_aabb(int cube_cells_per_dim)
    {
      CGAL_assertion(cube_cells_per_dim > 2);
      auto bbox = tree.bbox();

      Hexmeshing::Point center = {bbox.xmin() + (bbox.x_span()/2),
                      bbox.ymin() + (bbox.y_span()/2),
                      bbox.zmin() + (bbox.z_span()/2)};

      double max_size=std::max(std::max(bbox.x_span(), bbox.y_span()), bbox.z_span());
      grid = Hexmeshing::Grid::make_centered_cube
          (center, max_size / (cube_cells_per_dim-2), cube_cells_per_dim);
    }

    TriangleMesh poly;
    Tree tree;
    Hexmeshing::Grid grid;
  };
}

#endif
