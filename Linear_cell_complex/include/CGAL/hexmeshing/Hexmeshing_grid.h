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
// Contributor(s): Soichiro Yamazaki <soichiro19998@gmail.com>, Théo Bénard <benard320@gmail.com>
//
#ifndef HEXMESHING_GRID_H
#define HEXMESHING_GRID_H

#include <CGAL/hexmeshing/LCC_items_for_hexmeshing.h>
#include <CGAL/hexmeshing/Hexmeshing_generic_point.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Point_3.h>
#include <CGAL/Vector_3.h>
#include <array>
#include <cassert>


namespace CGAL::internal::Hexmeshing {
  /**
   * @brief A structure representing a 3D grid for hexahedral mesh generation
   *
   * This structure defines a regular grid in 3D space, specified by its position,
   * cell size, and dimensions. It provides methods for creating both centered and
   * non-centered grids, as well as specialized cube grid creation.
   */
  struct Grid {
    using Kernel = Exact_predicates_inexact_constructions_kernel;
    using Point = Kernel::Point_3;
    using PointInt = GenericPointForHexmeshing<int>;
    using LCCTraits = Linear_cell_complex_traits<3,Kernel>;
    using LCC = Linear_cell_complex_for_combinatorial_map<3,3, LCCTraits, LCCItemsForHexmeshing>;

    Point pos;      ///< Starting position (origin) of the grid
    Point size;     ///< Size of each cell in the grid
    PointInt dims;  ///< Number of cells in each dimension (x, y, z)
    PointInt dims_id_start;  ///< Starting indices for each dimension used only for parallelization

    /// @brief Default constructor
    Grid() {}

    /**
     * @brief Constructs a grid with specified parameters
     * @param from Starting position of the grid
     * @param cell_size Size of each cell
     * @param dims Number of cells in each dimension
     */
    Grid(Point from, Point cell_size, PointInt dims)
    : pos(from), size(cell_size), dims(dims) {}

    /**
     * @brief Creates a grid centered around a specified point
     * @param center Center point of the grid
     * @param cell_size Size of each cell
     * @param dims Number of cells in each dimension
     * @return A Grid object centered at the specified point
     */
    static Grid make_centered_grid(Point center, Point cell_size, PointInt dims) {
      Kernel::Vector_3 offset = {
        dims.x / 2 * cell_size.x(),
        dims.y / 2 * cell_size.y(),
        dims.z / 2 * cell_size.z(),
      };

      Point from = center - offset;
      return Grid(from, cell_size, dims);
    }

    /**
     * @brief Creates a grid starting from a specified point
     * @param from Starting position of the grid
     * @param cell_size Size of each cell
     * @param dims Number of cells in each dimension
     * @return A Grid object starting at the specified point
     */
    static Grid make_grid(Point from, Point cell_size, PointInt dims){
      return Grid(from, cell_size, dims);
    }

    /**
     * @brief Creates a cubic grid centered around a point
     * @param center Center point of the cube
     * @param cell_size Size of each cubic cell
     * @param dim Number of cells in each dimension (same for x, y, z)
     * @return A Grid object representing a cubic grid
     */
    static Grid make_centered_cube(Point center, double cell_size, int dim){
      Kernel::Vector_3 offset = {
        dim / 2 * cell_size,
        dim / 2 * cell_size,
        dim / 2 * cell_size,
      };

      Point from = center - offset;
      auto& s = cell_size;
      return Grid(from, {s,s,s}, {dim,dim,dim});
    }

    /**
     * @brief Creates a cubic grid starting from a point
     * @param from Starting position of the cube
     * @param cell_size Size of each cubic cell
     * @param dim Number of cells in each dimension (same for x, y, z)
     * @return A Grid object representing a cubic grid
     */
    static Grid make_cube(Point from, double cell_size, int dim){
      auto& s = cell_size;
      return Grid(from, {s,s,s}, {dim,dim,dim});
    }

    /**
     * @brief Generates a regular hexahedral grid in the Linear Cell Complex
     *
     * This function creates a regular 3D grid of hexahedral cells based on the provided
     * grid configuration. It generates the initial mesh structure that will be used
     * as the starting point for the hexahedral refinement algorithm.
     *
     * The function performs the following operations:
     *
     * 1. **Grid Iteration**: Iterates through all grid positions in three dimensions
     *    (x, y, z) based on the grid dimensions specified in `grid.dims`
     *
     * 2. **Coordinate Calculation**: For each grid position, calculates the coordinates
     *    of the eight vertices that define a hexahedral cell:
     *    - Bottom face vertices: (x1,y1,z1), (x2,y1,z1), (x2,y2,z1), (x1,y2,z1)
     *    - Top face vertices: (x1,y1,z2), (x2,y1,z2), (x2,y2,z2), (x1,y2,z2)
     *    Where x1, y1, z1 are the coordinates of the current grid position, and
     *    x2, y2, z2 are the coordinates of the next grid position
     *
     * 3. **Hexahedron Creation**: Creates a hexahedral cell at each grid position using
     *    `lcc.make_hexahedron()` with the calculated vertex coordinates
     *
     * 4. **Face Sewing**: After creating all hexahedra, calls `lcc.sew3_same_facets()`
     *    to properly connect adjacent faces and establish the correct topological
     *    relationships between neighboring cells
     *
     * The resulting mesh is a regular grid where each cell is a hexahedron with
     * faces properly connected to its neighbors. This provides the foundation for
     * the subsequent refinement operations.
     *
     * @param lcc The Linear Cell Complex where the grid will be created
     * @param grid Grid configuration containing position, size, and dimensions
     */
    static void generate_grid(LCC& lcc, Grid& grid) {
      for (int x = 0; x < grid.dims.x; x++) {
        for (int y = 0; y < grid.dims.y; y++) {
          for (int z = 0; z < grid.dims.z; z++) {
            double x1 = grid.pos.x() + x * grid.size.x();
            double y1 = grid.pos.y() + y * grid.size.y();
            double z1 = grid.pos.z() + z * grid.size.z();

            double x2 = grid.pos.x() + (x+1)*grid.size.x();
            double y2 = grid.pos.y() + (y+1)*grid.size.y();
            double z2 = grid.pos.z() + (z+1)*grid.size.z();

            lcc.make_hexahedron(Point(x1,y1,z1), Point(x2,y1,z1),
                                Point(x2,y2,z1), Point(x1,y2,z1),
                                Point(x1,y2,z2), Point(x1,y1,z2),
                                Point(x2,y1,z2), Point(x2,y2,z2));
          }
        }
      }

      lcc.sew3_same_facets();
    }
  };


}


#endif