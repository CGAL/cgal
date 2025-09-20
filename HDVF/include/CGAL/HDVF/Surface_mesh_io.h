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

#ifndef CGAL_HDVF_SURFACE_MESH_H
#define CGAL_HDVF_SURFACE_MESH_H

#include <CGAL/license/HDVF.h>

#include <iostream>
#include <cassert>
#include <vector>
#include <map>
#include <cmath>
#include <CGAL/HDVF/Mesh_object_io.h>
#include <CGAL/Surface_mesh.h>

namespace CGAL {
namespace Homological_discrete_vector_field {


/*!
 \ingroup PkgHDVFAlgorithmClasses

 The class `Surface_mesh_io` is an intermediate IO class, used to load `CGAL::Surface_mesh` and produce simplicial complexes.

 */

template <typename SurfaceMesh>
class Surface_mesh_io : public Mesh_object_io
{
public:
    typedef SurfaceMesh Surface_mesh;
private:
    
    /** \brief Storage of vertices Io_cell_type <-> Vertex_index permutation.
     *
     * `_io_cell_to_he_index` maps a Io_cell_type to its associated Halfedge_index in the mesh
     *
     * `_he_index_to_io_cell.at(q)` stores, for each dimension `q`,  the map between Halfedge_index and Io_cell_type of cells of dimension `q`
     *
     * All cells, whatever their dimension, are identified by halfedges:
     * - a vertex is the `target` of the halfedge
     * - an edge is the `edge` of the halfedge
     * - a face is the `face` of the halfedge
     */
    std::map<Io_cell_type, typename Surface_mesh::halfedge_index> _io_cell_to_he_index;
    std::vector<std::map<typename Surface_mesh::halfedge_index, Io_cell_type> > _he_index_to_io_cell;
    
public:
    /** \brief Constructor from a `CGAL::Surface_mesh` encoding a triangular mesh.
     *
     * Build a Surface_mesh_io from a `CGAL::Surface_mesh`.
     */
    Surface_mesh_io(const Surface_mesh& mesh) : Mesh_object_io(2) {
        typedef typename Surface_mesh::Point Point;
        
        this->nvertices = mesh.vertices().size();
        this->ncells = mesh.vertices().size() + mesh.edges().size() + mesh.faces().size();
        _he_index_to_io_cell.resize(3);
        
        // Load nodes
        typename Surface_mesh::Vertex_range vr = mesh.vertices();
        for (typename Surface_mesh::Vertex_range::iterator it = vr.begin(); it != vr.end(); ++it) {
            Point p(mesh.point(*it));
            Io_node_type node;
            for (size_t q=0; q<p.dimension(); ++q)
                node.push_back(p[q]);
            this->nodes.push_back(node);
        }
        // Load vertices
        size_t tmp_index(0);
        for (typename Surface_mesh::Vertex_range::iterator it = vr.begin(); it != vr.end(); ++it) {
            // Corresponding Io_cell_type
            Io_cell_type tmp_io_cell({tmp_index});
            // Add to cells
            this->cells.push_back(tmp_io_cell);
            // Associated he
            typename Surface_mesh::Halfedge_index vertex_he(mesh.halfedge(*it));
            // Store the permutation
            _io_cell_to_he_index[tmp_io_cell] = vertex_he;
            _he_index_to_io_cell.at(0)[vertex_he] = tmp_io_cell;
            ++tmp_index;
        }
        // Load edges
        typename Surface_mesh::Edge_range er = mesh.edges();
        for (typename Surface_mesh::Edge_range::iterator it = er.begin(); it != er.end(); ++it) {
            // Associated he
            typename SurfaceMesh::Halfedge_index edge_he(mesh.halfedge(*it));
            // Compute corresponding Io_cell
            std::vector<size_t> tmp_cell;
            // Vertex1
            {
                // Get the halfedge "encoding" the target vertex
                typename Surface_mesh::Halfedge_index vert_he_ind(mesh.halfedge(mesh.target(edge_he)));
                // Get corresponding io_cell
                Io_cell_type vert(_he_index_to_io_cell.at(0).at(vert_he_ind));
                assert(vert.size()==1); // vertex
                // Push_back the index
                tmp_cell.push_back(vert.at(0));
            }
            // Vertex2
            {
                // Get the halfedge "encoding" the source vertex
                typename Surface_mesh::Halfedge_index vert_he_ind(mesh.halfedge(mesh.target(mesh.opposite(edge_he))));
                // Get the corresponding io_cell
                Io_cell_type vert(_he_index_to_io_cell.at(0).at(vert_he_ind));
                assert(vert.size()==1); // vertex
                // Push_back the index
                tmp_cell.push_back(vert.at(0));
            }
            std::sort(tmp_cell.begin(), tmp_cell.end());
            // Add to cells
            this->cells.push_back(tmp_cell);
            // Set the permutation
            _io_cell_to_he_index[tmp_cell] = edge_he;
            _he_index_to_io_cell.at(1)[edge_he] = tmp_cell;
        }
        // Load faces
        typename Surface_mesh::Face_range fr = mesh.faces();
        for (typename Surface_mesh::Face_range::iterator it = fr.begin(); it != fr.end(); ++it) {
            // Associated he
            typename SurfaceMesh::Halfedge_index face_he(mesh.halfedge(*it));
            // Compute corresponding Io_cell
            std::vector<size_t> tmp_cell;
            // Visit vertices around the face
            CGAL::Vertex_around_face_circulator<Surface_mesh> hebegin(mesh.halfedge(*it), mesh), hedone(hebegin);
            size_t cpt_verts(0);
            do {
                ++cpt_verts;
                // Get the halfedge stored in the vertex
                typename Surface_mesh::Halfedge_index vert_he_ind(mesh.halfedge(*hebegin));
                // Get the corresponding io_cell
                Io_cell_type vert(_he_index_to_io_cell.at(0).at(vert_he_ind));
                assert(vert.size()==1); // vertex
                // Push_back the index
                tmp_cell.push_back(vert.at(0));
                ++hebegin;
            } while (hebegin != hedone);
            std::cout << cpt_verts << std::endl;
            std::sort(tmp_cell.begin(), tmp_cell.end());
            // Add to cells
            this->cells.push_back(tmp_cell);
            // Set the permutation
            _io_cell_to_he_index[tmp_cell] = face_he;
            _he_index_to_io_cell.at(2)[face_he] = tmp_cell;
        }
    }
} ;

} /* end namespace Homological_discrete_vector_field */
} /* end namespace CGAL */


#endif // CGAL_HDVF_SURFACE_MESH_H
