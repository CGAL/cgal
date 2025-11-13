// Copyright (c) 2025 LIS Marseille (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Alexandra Bac <alexandra.bac@univ-amu.fr>

#ifndef CGAL_HDVF_TRIANGULATION_3_IO_H
#define CGAL_HDVF_TRIANGULATION_3_IO_H

#include <CGAL/license/HDVF.h>

#include <iostream>
#include <cassert>
#include <vector>
#include <map>
#include <cmath>
#include <CGAL/boost/graph/IO/polygon_mesh_io.h>
#include <CGAL/HDVF/Mesh_object_io.h>

namespace CGAL {
namespace Homological_discrete_vector_field {


/*!
 \ingroup PkgHDVFAlgorithmClasses

 The class `Triangulation_3_io` is an intermediate IO class, used to load a `Triangulation_3` and produce simplicial complexes. The class loads the Vertices and the Cells (ie. tetrahedra) of the `Triangulation_3` into a `Mesh_object_io`.
\tparam Triangulation3 a model of  `CGAL::Triangulation_3`.
\tparam Traits a geometric traits class model of the `HDVFTraits` concept.
 */

template <typename Triangulation3, typename Traits>
class Triangulation_3_io : public Mesh_object_io<Traits>
{
public:
    typedef Triangulation3 Triangulation_3;
    typedef typename Triangulation3::Vertex_handle vertex_descriptor ;
    typedef typename Triangulation3::Cell_handle cell_descriptor ;
    typedef typename Triangulation3::Point Point ;
private:

    /** \brief Storage of vertices Io_cell_type <-> Cell_handle permutation.
     *
     * Vertices maps:
     *- `_io_cell_to_vertex_handle` maps a vertex Io_cell_type to its associated Vertex_handle in the triangulation
     *- `_vertex_handle_to_io_cell` stores  the map between Vertex_handle and Io_cell_type of vertices (of dimension 0).
     *
     *Cells maps (dimension 3):
     *- `_io_cell_to_cell_handle` maps a Io_cell_type to its associated Cell_handle in the triangulation
     *
     *- `_cell_handle_to_io_cell` stores  the map between Cell_handle and Io_cell_type of Cells (of dimension 3).
     */
    std::map<size_t, vertex_descriptor> _io_cell_to_vertex_handle;
    std::map<vertex_descriptor, size_t> _vertex_handle_to_io_cell;
    
    std::map<Io_cell_type, cell_descriptor> _io_cell_to_cell_handle;
    std::map<cell_descriptor, Io_cell_type> _cell_handle_to_io_cell;

public:
    /** \brief Constructor from a `Triangulation_3`.
     *
     */
    Triangulation_3_io(const Triangulation_3& triangulation) : Mesh_object_io<Traits>(3) {
        typedef typename Traits::Point Point;

        this->nvertices = triangulation.number_of_vertices();
        this->ncells = triangulation.number_of_vertices() + triangulation.number_of_finite_cells();

        // Load nodes
        for (Point p : triangulation.points()) {
            this->nodes.push_back(p);
        }
        // Load vertices
        size_t tmp_index(0);
        for (vertex_descriptor v : triangulation.finite_vertex_handles()) {
            // Corresponding Io_cell_type
            Io_cell_type tmp_io_cell({tmp_index});
            // Add to cells
            this->cells.push_back(tmp_io_cell);
            // Store the permutation
            _io_cell_to_vertex_handle[tmp_index] = v;
            _vertex_handle_to_io_cell[v] = tmp_index;
            ++tmp_index;
        }
        // Load cells
        for (cell_descriptor c : triangulation.finite_cell_handles()) {
            // Compute corresponding Io_cell
            std::vector<size_t> tmp_cell;
            // Visit vertices around the face to get the indices of vertices
            std::array<vertex_descriptor, 4> verts(triangulation.vertices(c));
            for (int i=0; i<4; ++i){
                tmp_cell.push_back(_vertex_handle_to_io_cell[verts[i]]);
            }
            // Sort this vector
            std::sort(tmp_cell.begin(), tmp_cell.end());
            // Add to cells
            this->cells.push_back(tmp_cell);
            // Set the permutation
            _io_cell_to_cell_handle[tmp_cell] = c;
            _cell_handle_to_io_cell[c] = tmp_cell;
        }
    }
} ;

} /* end namespace Homological_discrete_vector_field */
} /* end namespace CGAL */


#endif // CGAL_HDVF_TRIANGULATION_3_IO_H
