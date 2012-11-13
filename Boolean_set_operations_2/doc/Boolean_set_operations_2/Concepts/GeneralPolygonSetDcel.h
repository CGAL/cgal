
/*!
\ingroup PkgBooleanSetOperations2Concepts
\cgalConcept

A doubly-connected edge-list (<span class="textsc">Dcel</span> for short) data-structure. It consists 
of three containers of records: vertices \f$ V\f$, halfedges \f$ E\f$, and faces \f$ F\f$. 
It maintains the incidence relation among them. The halfedges are ordered 
in pairs sometimes referred to as twins, such that each halfedge pair 
represent an edge. 

A model of the `GeneralPolygonSetDcel` concept must provide the following types and 
operations. (In addition to the requirements here, the local types 
`Vertex`,`Halfedge`, `Face`, 
`Hole`, and `Isolated_vertex` 
must be models of the concepts 
`ArrangementDcelVertex`, 
`ArrangementDcelHalfedge`, 
`GeneralPolygonSetDcelFace` 
, 
`ArrangementDcelHole`, and 
`ArrangementDcelIsolatedVertex` 
respectively.) 
Notice that this concept differs from the concept `ArrangemenDcel` 
only in the type `Face`. 

\cgalHasModel `CGAL::Arr_dcel_base<V,H,F>` 
\cgalHasModel `CGAL::Arr_default_dcel<Traits>` 
\cgalHasModel `CGAL::Arr_face_extended_dcel<Traits,FData,V,H,F>` 
\cgalHasModel `CGAL::Arr_extended_dcel<Traits,VData,HData,FData,V,H,F>` 

\sa `ArrangementDcelVertex` 
\sa `ArrangementDcelHalfedge` 
\sa `GeneralPolygonSetDcelFace` 
\sa `ArrangementDcelHole` 
\sa `ArrangementDcelIsolatedVertex` 

*/

class GeneralPolygonSetDcel {
public:

/// \name Types 
/// @{

/*! 
the vertex type. 
*/ 
typedef Hidden_type Vertex; 

/*! 
the halfedge type. 
*/ 
typedef Hidden_type Halfedge; 

/*! 
the face type. 
*/ 
typedef Hidden_type Face; 

/*! 
the hole type. 
*/ 
typedef Hidden_type Hole; 

/*! 
the isolated vertex type. 
*/ 
typedef Hidden_type Isolated_vertex; 

/*! 
used to represent size values (e.g., `size_t`). 
*/ 
typedef Hidden_type Size; 

/// @}

/// \name Iterators
/// The non-mutable iterators `Vertex_const_iterator`,
/// `Halfedge_const_iterator` and `Face_const_iterator` are also
/// defined.
/// @{

/*! 
a bidirectional iterator over the vertices. Its value-type is 
`Vertex`. 
*/ 
typedef Hidden_type Vertex_iterator; 

/*! 
a bidirectional iterator over the halfedges. Its value-type is 
`Halfedge`. 
*/ 
typedef Hidden_type Halfedge_iterator; 

/*! 
a bidirectional iterator over the faces. Its value-type is `Face`. 
*/ 
typedef Hidden_type Face_iterator; 

/// @} 

/// \name Creation 
/// @{

/*! 
constructs an empty <span class="textsc">Dcel</span> with one unbounded face. 
*/ 
Arr_dcel(); 

/*! 
assigns the contents of the `other` <span class="textsc">Dcel</span> whose unbounded face 
is given by `uf`, to `dcel`. The function returns a pointer to 
the unbounded face of `dcel` after the assignment. 
*/ 
Face* assign (const Self& other, const Face *uf); 

/// @} 

/// \name Access Functions 
/// @{

/*! 
returns the number of vertices.
*/ 
Size size_of_vertices() const; 

/*! 
returns the number of halfedges (always even). 
*/ 
Size size_of_halfedges() const; 

/*! 
returns the number of faces. 
*/ 
Size size_of_faces() const; 

/*! 
returns the number of holes (the number of connected components). 
*/ 
Size size_of_holes() const; 

/*! 
returns the number of isolated vertices. 
*/ 
Size size_of_isolated_vertices() const; 

/// @}

/// \name Iterator Access
/// The following operations have an equivalent `const` operations
/// that return the corresponding non-mutable iterators:
/// @{

/*!
returns a begin-iterator of the vertices in `dcel`.
*/ 
Vertex_iterator vertices_begin(); 

/*! 
returns a past-the-end iterator of the vertices in `dcel`. 
*/ 
Vertex_iterator vertices_end(); 

/*! 
returns a begin-iterator of the halfedges in `dcel`. 
*/ 
Halfedge_iterator halfedges_begin(); 

/*! 
returns a past-the-end iterator of the halfedges in `dcel`. 
*/ 
Halfedge_iterator halfedges_end(); 

/*! 
returns a begin-iterator of the faces in `dcel`. 
*/ 
Vertex_iterator faces_begin(); 

/*! 
returns a past-the-end iterator of the faces in `dcel`. 
*/ 
Vertex_iterator faces_end(); 

/// @} 

/// \name Modifiers 
/// The following operations allocate a new element of the respective
/// type. Halfedges are always allocated in pairs of opposite
/// halfedges. The halfedges and their opposite pointers are
/// automatically set.
/// @{

/*! 
creates a new vertex. 
*/ 
Vertex* new_vertex(); 

/*! 
creates a new pair of twin halfedges. 
*/ 
Halfedge* new_edge(); 

/*! 
creates a new face. 
*/ 
Face* new_face(); 

/*! 
creates a new hole record. 
*/ 
Hole* new_hole(); 

/*! 
creates a new isolated vertex record. 
*/ 
Isolated_vertex* new_isolated_vertex(); 

/*! 
deletes the vertex `v`. 
*/ 
void delete_vertex(Vertex* v); 

/*! 
deletes the halfedge `e` as well as its twin. 
*/ 
void delete_edge(Halfedge* e); 

/*! 
deletes the face `f`. 
*/ 
void delete_face(Face* f); 

/*! 
deletes the hole `ho`. 
*/ 
void delete_hole(Hole* ho); 

/*! 
deletes the isolated vertex `iv`. 
*/ 
void delete_isolated_vertex(Isolated_vertex* iv); 

/// @}

}; /* end GeneralPolygonSetDcel */

