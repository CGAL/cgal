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


/// \ingroup PkgSurfaceParameterizationConcepts
/// \cgalconcept
/// ParameterizationPatchableMesh_3 inherits from concept ParameterizationMesh_3,
/// thus is a concept of a 3D surface mesh.
///
/// ParameterizationPatchableMesh_3 adds support for patches and virtual seams.
/// "Patches" are a subset of a 3D mesh. "Virtual seams" are the ability
/// to behave as if the surface was cut following a certain path.
///
/// This mainly means that:
/// - vertices can be tagged as inside or outside the patch to parameterize.
/// - the fields specific to parameterizations (`index`, `u`, `v`, `is_parameterized`)
///   can be set per "corner" (aka half-edge).
///
/// This allows to parameterize any 3D surface by decomposing it as a list of topological disks.
///

/// \refines ParameterizationMesh_3
/// \hasModel `CGAL::Parameterization_polyhedron_adaptor_3`
class ParameterizationPatchableMesh_3
{
// Public types
public:

    // Same sub-types as ParameterizationMesh_3


// Public operations
public:

    // Construction and destruction are undefined.

    // VERTEX INTERFACE

    /// Get vertex seaming flag. Default value is undefined.
    int  get_vertex_seaming(Vertex_const_handle vertex) const;

    /// Set vertex seaming flag. Default value is undefined.
    void set_vertex_seaming(Vertex_handle vertex, int seaming);


    // EDGE INTERFACE

    /// Get oriented edge's seaming flag, i.e.\ position of the oriented edge
    /// w.r.t.\ to the unique main border.
    int  get_halfedge_seaming(Vertex_const_handle source, Vertex_const_handle target) const;


    /// Set oriented edge's seaming flag, i.e.\ position of the oriented edge
    /// w.r.t.\ to the unique main border.
    void set_halfedge_seaming(Vertex_handle source, Vertex_handle target, int seaming);


    // CORNER INTERFACE

    /// Get the 2D position (= `(u,v)` pair) of corners at the "right"
    /// of the `prev_vertex` -> `vertex` -> `next_vertex` line.
    /// Default value is undefined.
    Point_2 get_corners_uv(Vertex_const_handle vertex,
                           Vertex_const_handle prev_vertex,
                           Vertex_const_handle next_vertex) const;


    /// Set the 2D position (= `(u,v)` pair) of corners at the "right"
    /// of the `prev_vertex` -> `vertex` -> `next_vertex` line.
    /// Default value is undefined.
    void set_corners_uv(Vertex_handle vertex,
                        Vertex_const_handle prev_vertex,
                        Vertex_const_handle next_vertex,
                        const Point_2& uv);

    /// Get `is_parameterized` field of corners at the "right"
    /// of the `prev_vertex` -> `vertex` -> `next_vertex` line.
    /// Default value is undefined.
    bool are_corners_parameterized(Vertex_const_handle vertex,
                                   Vertex_const_handle prev_vertex,
                                   Vertex_const_handle next_vertex) const;

    /// Set `is_parameterized` field of corners at the "right"
    /// of the `prev_vertex` -> `vertex` -> `next_vertex` line.
    /// Default value is undefined.
    void set_corners_parameterized(Vertex_handle vertex,
                                   Vertex_const_handle prev_vertex,
                                   Vertex_const_handle next_vertex,
                                   bool parameterized);

    /// Get `index` of corners at the "right"
    /// of the `prev_vertex` -> `vertex` -> `next_vertex line`.
    /// Default value is undefined.
    int get_corners_index(Vertex_const_handle vertex,
                          Vertex_const_handle prev_vertex,
                          Vertex_const_handle next_vertex) const;

    /// Set `is_parameterized` field of corners at the "right"
    /// of the `prev_vertex` -> `vertex` -> `next_vertex` line.
    /// Default value is undefined.
    void set_corners_index(Vertex_handle vertex,
                           Vertex_const_handle prev_vertex,
                           Vertex_const_handle next_vertex,
                           int index);

    /// Get all purpose tag of corners at the "right"
    /// of the `prev_vertex` -> `vertex` -> `next_vertex` line.
    /// Default value is undefined.
    int get_corners_tag(Vertex_const_handle vertex,
                        Vertex_const_handle prev_vertex,
                        Vertex_const_handle next_vertex) const;

    /// Set all purpose tag of corners at the "right"
    /// of the `prev_vertex` -> `vertex` -> `next_vertex` line.
    /// Default value is undefined.
    void set_corners_tag(Vertex_handle vertex,
                         Vertex_const_handle prev_vertex,
                         Vertex_const_handle next_vertex,
                         int tag);

}; /// ParameterizationPatchableMesh_3

