// Copyright (c) 2005  INRIA (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Laurent Saboret, Pierre Alliez, Bruno Levy


/// ParameterizationPatchableMesh_3 inherits from concept ParameterizationMesh_3,
/// thus is a concept of a 3D surface mesh.
///
/// ParameterizationPatchableMesh_3 adds the ability to support patches and virtual seams.
/// "Patches" are a subset of a 3D mesh. "Virtual seams" are the ability
/// to behave exactly as if the surface was cut following a certain path.
///
/// This mainly means that:
/// - vertices can be tagged as inside or outside the patch to parameterize.
/// - the fields specific to parameterizations (index, u, v, is_parameterized)
///   can be set per "corner" (aka half-edge).
///
/// The main purpose of this feature is to allow the Surface_mesh_parameterization package
/// to parameterize any 3D surface by decomposing it as a list of topological disks.
///
/// @heading Design Pattern:
/// ParameterizationPatchableMesh_3 is an Adaptor [GHJV95]: it changes the
/// interface of a 3D mesh to match the interface expected by class Parameterization_mesh_patch_3.
///
/// @heading Has Models:
/// Adaptator for Polyhedron_3 is provided.

class ParameterizationPatchableMesh_3 : public ParameterizationMesh_3
{
// Public types
public:

    // Same sub-types as ParameterizationMesh_3


// Public operations
public:

    // Construction and destruction are undefined.

    // VERTEX INTERFACE

    /// Get/set vertex seaming flag. Default value is undefined.
    int  get_vertex_seaming(Vertex_const_handle vertex) const;
    void set_vertex_seaming(Vertex_handle vertex, int seaming);


    // EDGE INTERFACE

    /// Get/set oriented edge's seaming flag, i.e. position of the oriented edge
    /// w.r.t. to the UNIQUE main border.
    int  get_halfedge_seaming(Vertex_const_handle source, Vertex_const_handle target) const;
    void set_halfedge_seaming(Vertex_handle source, Vertex_handle target, int seaming);


    // CORNER INTERFACE

    /// Get/set the 2D position (= (u,v) pair) of corners at the "right"
    /// of the prev_vertex -> vertex -> next_vertex line.
    /// Default value is undefined.
    Point_2 get_corners_uv(Vertex_const_handle vertex,
                           Vertex_const_handle prev_vertex,
                           Vertex_const_handle next_vertex) const;
    void set_corners_uv(Vertex_handle vertex,
                        Vertex_const_handle prev_vertex,
                        Vertex_const_handle next_vertex,
                        const Point_2& uv);

    /// Get/set "is parameterized" field of corners at the "right"
    /// of the prev_vertex -> vertex -> next_vertex line.
    /// Default value is undefined.
    bool are_corners_parameterized(Vertex_const_handle vertex,
                                   Vertex_const_handle prev_vertex,
                                   Vertex_const_handle next_vertex) const;
    void set_corners_parameterized(Vertex_handle vertex,
                                   Vertex_const_handle prev_vertex,
                                   Vertex_const_handle next_vertex,
                                   bool parameterized);

    /// Get/set index of corners at the "right"
    /// of the prev_vertex -> vertex -> next_vertex line.
    /// Default value is undefined.
    int get_corners_index(Vertex_const_handle vertex,
                          Vertex_const_handle prev_vertex,
                          Vertex_const_handle next_vertex) const;
    void set_corners_index(Vertex_handle vertex,
                           Vertex_const_handle prev_vertex,
                           Vertex_const_handle next_vertex,
                           int index);

    /// Get/set all purpose tag of corners at the "right"
    /// of the prev_vertex -> vertex -> next_vertex line.
    /// Default value is undefined.
    int get_corners_tag(Vertex_const_handle vertex,
                        Vertex_const_handle prev_vertex,
                        Vertex_const_handle next_vertex) const;
    void set_corners_tag(Vertex_handle vertex,
                         Vertex_const_handle prev_vertex,
                         Vertex_const_handle next_vertex,
                         int tag);

}; /// ParameterizationPatchableMesh_3

