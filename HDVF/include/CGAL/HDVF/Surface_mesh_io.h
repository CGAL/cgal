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

#ifndef CGAL_HDVF_SURFACE_MESH_IO_H
#define CGAL_HDVF_SURFACE_MESH_IO_H

#include <CGAL/license/HDVF.h>

#include <iostream>
#include <cassert>
#include <vector>
#include <map>
#include <cmath>
#include <CGAL/HDVF/Mesh_object_io.h>

namespace CGAL {
namespace Homological_discrete_vector_field {


/*!
 \ingroup PkgHDVFAlgorithmClasses

 The class `Surface_mesh_io` is an intermediate IO class, used to load a triangle mesh and produce simplicial complexes.
\tparam TriangleMesh a model of `FaceGraph` and `HalfedgeGraph` concepts, e.g., a `CGAL::Surface_mesh`.
\tparam Traits a geometric traits class model of the `HDVFTraits` concept.
 */

template <typename TriangleMesh, typename Traits>
class Surface_mesh_io : public Mesh_object_io<Traits>
{
public:
    typedef TriangleMesh Surface_mesh;
    typedef typename boost::graph_traits<Surface_mesh>::vertex_descriptor vertex_descriptor ;
    typedef typename boost::graph_traits<Surface_mesh>::halfedge_descriptor halfedge_descriptor ;
    typedef typename boost::graph_traits<Surface_mesh>::edge_descriptor edge_descriptor ;
    typedef typename boost::graph_traits<Surface_mesh>::face_descriptor face_descriptor ;
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
    std::map<Io_cell_type, halfedge_descriptor> _io_cell_to_he_index;
    std::vector<std::map<halfedge_descriptor, Io_cell_type> > _he_index_to_io_cell;

public:
    /** \brief Constructor from atriangule mesh.
     *
     */
    Surface_mesh_io(const TriangleMesh& mesh) : Mesh_object_io(2) {
        typedef typename Traits::Point Point;

        this->nvertices = num_vertices(mesh);
        this->ncells = num_vertices(mesh) + num_edges(mesh) + num_faces(mesh) ;
        _he_index_to_io_cell.resize(3);

        auto vpm = get(CGAL::vertex_point, mesh);
        // Load nodes
        for (vertex_descriptor v : vertices(mesh)) {
            Point p(get(vpm,v));
            this->nodes.push_back(p);
        }
        // Load vertices
        size_t tmp_index(0);
        for (vertex_descriptor v : vertices(mesh)) {
            // Corresponding Io_cell_type
            Io_cell_type tmp_io_cell({tmp_index});
            // Add to cells
            this->cells.push_back(tmp_io_cell);
            // Associated he
            halfedge_descriptor vertex_he(halfedge(v,mesh));
            // Store the permutation
            _io_cell_to_he_index[tmp_io_cell] = vertex_he;
            _he_index_to_io_cell.at(0)[vertex_he] = tmp_io_cell;
            ++tmp_index;
        }
        // Load edges
        for (edge_descriptor e : edges(mesh)) {
            // Associated he
            halfedge_descriptor edge_he(halfedge(e,mesh));
            // Compute corresponding Io_cell
            std::vector<size_t> tmp_cell;
            // Vertex1
            {
                // Get the halfedge "encoding" the target vertex
                halfedge_descriptor vert_he_ind(halfedge(target(edge_he, mesh), mesh));
                // Get corresponding io_cell
                Io_cell_type vert(_he_index_to_io_cell.at(0).at(vert_he_ind));
                assert(vert.size()==1); // vertex
                // Push_back the index
                tmp_cell.push_back(vert.at(0));
            }
            // Vertex2
            {
                // Get the halfedge "encoding" the source vertex
                halfedge_descriptor vert_he_ind(halfedge(target(opposite(edge_he, mesh), mesh), mesh));
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
        for (face_descriptor f : faces(mesh)) {
            // Associated he
            halfedge_descriptor face_he(halfedge(f,mesh));
            // Compute corresponding Io_cell
            std::vector<size_t> tmp_cell;
            // Visit vertices around the face
            size_t cpt_verts(0);
            for(vertex_descriptor v : vertices_around_face(halfedge(f,mesh), mesh)){
                ++cpt_verts;
                // Get the halfedge stored in the vertex
                halfedge_descriptor vert_he_ind(halfedge(v,mesh));
                // Get the corresponding io_cell
                Io_cell_type vert(_he_index_to_io_cell.at(0).at(vert_he_ind));
                assert(vert.size()==1); // vertex
                // Push_back the index
                tmp_cell.push_back(vert.at(0));
            }
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


#endif // CGAL_HDVF_SURFACE_MESH_IO_H
