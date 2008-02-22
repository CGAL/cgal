// Copyright (c) 2005  INRIA (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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


/// ParameterizationMesh_3 is a concept for a 3D surface mesh.
/// Its main purpose is to allow the parameterization methods to access meshes in a uniform manner.
///
/// A ParameterizationMesh_3 surface consists of vertices, facets and an incidence relation on them.
/// No notion of edge is requested.
/// Vertices represent points in 3d-space. Facets are planar polygons without holes
/// defined by the circular sequence of vertices along their border.
/// The surface itself can have holes. The vertices
/// along the border of a hole are called "border vertices".
/// A surface is "closed" if it contains no border vertices.
///
/// The surface must be an oriented 2-manifold with border vertices, i.e.
/// the neighborhood of each point on the surface is either
/// homeomorphic to a disc or to a half disc, except for vertices where
/// many holes and surfaces with border can join.
///
/// ParameterizationMesh_3 defines the types, data and methods that a mesh must implement
/// to allow surface parameterization.
/// Among other things, this concept defines accessors to fields specific
/// to parameterizations methods: index, u, v, is_parameterized.
///
/// ParameterizationMesh_3 meshes can have any genus, arity or number of components. On the other hand,
/// as parameterization methods deal only with topological disks, ParameterizationMesh_3
/// defines an interface oriented towards topological disks.
///
/// @heading Has Models:
/// We provide 2 models of this concept:
/// - Parameterization_polyhedron_adaptor_3<Polyhedron_3_>
/// - Parameterization_mesh_patch_3<ParameterizationPatchableMesh_3>
///
/// @heading Design Pattern:
/// ParameterizationMesh_3 is an Adaptor [GHJV95]: it changes the
/// interface of a 3D mesh to match the interface expected by the parameterization methods.

class ParameterizationMesh_3
{
// Public types
public:

    /// Number type to represent coordinates.
    typedef xxx NT;

    /// 2D point that represents (u,v) coordinates computed
    /// by parameterization methods. Must provide X() and Y() methods.
    typedef xxx Point_2;
    /// 3D point that represents vertices coordinates. Must provide X() and Y() methods.
    typedef xxx Point_3;
    /// 2D vector. Must provide X() and Y() methods.
    typedef xxx Vector_2;
    /// 3D vector. Must provide X() and Y() methods.
    typedef xxx Vector_3;

    /// Opaque type representing a facet of the 3D mesh. No methods are expected.
    typedef xxx Facet;
    /// Handle to a facet. Model of the Handle concept.
    typedef xxx Facet_handle;
    typedef xxx Facet_const_handle;
    /// Iterator over all mesh facets. Model of the ForwardIterator concept.
    typedef xxx Facet_iterator;
    typedef xxx Facet_const_iterator;

    /// Opaque type representing a vertex of the 3D mesh. No methods are expected.
    typedef xxx Vertex;
    /// Handle to a vertex. Model of the Handle concept.
    typedef xxx Vertex_handle;
    typedef xxx Vertex_const_handle;
    /// Iterator over all vertices of a mesh. Model of the ForwardIterator concept.
    typedef xxx Vertex_iterator;
    typedef xxx Vertex_const_iterator;
    /// Iterator over vertices of the mesh "main border".
    /// Model of the ForwardIterator concept.
    typedef xxx Border_vertex_iterator;
    typedef xxx Border_vertex_const_iterator;
    /// Counter-clockwise circulator over a facet's vertices.
    /// Model of the BidirectionalCirculator concept.
    typedef xxx Vertex_around_facet_circulator;
    typedef xxx Vertex_around_facet_const_circulator;
    /// Clockwise circulator over the vertices incident to a vertex.
    /// Model of the BidirectionalCirculator concept.
    typedef xxx Vertex_around_vertex_circulator;
    typedef xxx Vertex_around_vertex_const_circulator;


// Public operations
public:

    // Construction and destruction are undefined.

    // MESH INTERFACE

    /// Indicate if the mesh matches the ParameterizationMesh_3 concept.
    bool is_valid() const;

    /// Get iterator over first vertex of mesh.
    Vertex_iterator  mesh_vertices_begin ();
    Vertex_const_iterator  mesh_vertices_begin () const;

    /// Get iterator over past-the-end vertex of mesh.
    Vertex_iterator  mesh_vertices_end ();
    Vertex_const_iterator  mesh_vertices_end () const;

    /// Count the number of vertices of the mesh.
    int  count_mesh_vertices () const;

    /// Index vertices of the mesh from 0 to count_mesh_vertices()-1.
    void  index_mesh_vertices ();

    /// Get iterator over first vertex of mesh's "main border".
    Border_vertex_iterator  mesh_main_border_vertices_begin ();
    Border_vertex_const_iterator  mesh_main_border_vertices_begin () const;

    /// Get iterator over past-the-end vertex of mesh's "main border".
    Border_vertex_iterator  mesh_main_border_vertices_end ();
    Border_vertex_const_iterator  mesh_main_border_vertices_end () const;

    /// Return the border containing seed_vertex.
    /// Return an empty list if not found.
    std::list<Vertex_handle> get_border(Vertex_handle seed_vertex);

    /// Get iterator over first facet of mesh
    Facet_iterator  mesh_facets_begin ();
    Facet_const_iterator  mesh_facets_begin () const;

    /// Get iterator over past-the-end facet of mesh.
    Facet_iterator  mesh_facets_end ();
    Facet_const_iterator  mesh_facets_end () const;

    /// Count the number of facets of the mesh.
    int  count_mesh_facets () const;

    /// Return true of all mesh's facets are triangles.
    bool  is_mesh_triangular () const;

    /// Count the number of halfedges of the mesh.
    int  count_mesh_halfedges() const;

    // FACET INTERFACE

    /// Get circulator over facet's vertices.
    Vertex_around_facet_circulator facet_vertices_begin(Facet_handle facet);
    Vertex_around_facet_const_circulator  facet_vertices_begin(Facet_const_handle facet) const;

    /// Count the number of vertices of a facet.
    int  count_facet_vertices(Facet_const_handle facet) const;

    // VERTEX INTERFACE

    /// Get the 3D position of a vertex.
    Point_3  get_vertex_position (Vertex_const_handle vertex) const;

    /// Get/set the 2D position (u/v pair) of a vertex. Default value is undefined.
    Point_2  get_vertex_uv (Vertex_const_handle vertex) const;
    void  set_vertex_uv (Vertex_handle vertex, const Point_2& uv);

    /// Get/set "is parameterized" field of vertex. Default value is undefined.
    bool  is_vertex_parameterized (Vertex_const_handle vertex) const;
    void  set_vertex_parameterized (Vertex_handle vertex, bool parameterized);

    /// Get/set vertex index. Default value is undefined.
    int  get_vertex_index (Vertex_const_handle vertex) const;
    void  set_vertex_index (Vertex_handle vertex, int index);

    /// Get/set vertex' all purpose tag. Default value is undefined.
    int  get_vertex_tag(Vertex_const_handle vertex) const;
    void set_vertex_tag(Vertex_handle vertex, int tag);

    /// Return true if a vertex belongs to ANY mesh's border.
    bool  is_vertex_on_border(Vertex_const_handle vertex) const;

    /// Return true if a vertex belongs to the UNIQUE mesh's main border.
    bool  is_vertex_on_main_border(Vertex_const_handle vertex) const;

    /// Get circulator over the vertices incident to 'vertex'.
    /// 'start_position' defines the optional initial position of the circulator.
    Vertex_around_vertex_circulator vertices_around_vertex_begin(
                            Vertex_handle vertex,
                            Vertex_handle start_position = Vertex_handle());
    Vertex_around_vertex_const_circulator vertices_around_vertex_begin(
                            Vertex_const_handle vertex,
                            Vertex_const_handle start_position = Vertex_const_handle()) const;

}; /// ParameterizationMesh_3

