
/*!
\ingroup PkgArrangementOnSurface2ConceptsDCEL
\cgalConcept

A doubly-connected edge-list (\dcel for short) data-structure. It consists
of three containers of records: vertices \f$ V\f$, halfedges \f$ E\f$, and faces \f$ F\f$.
It maintains the incidence relation among them. The halfedges are ordered
in pairs sometimes referred to as twins, such that each halfedge pair
represent an edge.

A model of the `ArrangementDcel` concept must provide the following types and
operations. (In addition to the requirements here, the local types
`Vertex`,`Halfedge`, `Face`
`Hole` and `Isolated_vertex`
must be models of the concepts
`ArrangementDcelVertex`,
`ArrangementDcelHalfedge`,
`ArrangementDcelFace`,
`ArrangementDcelHole`, and
`ArrangementDcelIsolatedVertex` respectively.)

\cgalHasModel `CGAL::Arr_dcel_base<V,H,F>`
\cgalHasModel `CGAL::Arr_default_dcel<Traits>`
\cgalHasModel `CGAL::Arr_face_extended_dcel<Traits,FData,V,H,F>`
\cgalHasModel `CGAL::Arr_extended_dcel<Traits,VData,HData,FData,V,H,F>`

\sa `ArrangementDcelVertex`
\sa `ArrangementDcelHalfedge`
\sa `ArrangementDcelFace`
\sa `ArrangementDcelHole`
\sa `ArrangementDcelIsolatedVertex`

*/

class ArrangementDcel {
public:

/// \name Types
/// @{

/*!
the vertex type.
*/
typedef unspecified_type Vertex;

/*!
the halfedge type.
*/
typedef unspecified_type Halfedge;

/*!
the face type.
*/
typedef unspecified_type Face;

/*!
the hole type.
*/
typedef unspecified_type Hole;

/*!
the isolated vertex type.
*/
typedef unspecified_type Isolated_vertex;

/*!
used to represent size values (e.g., `size_t`).
*/
typedef unspecified_type Size;

/*!
a bidirectional iterator over the vertices. Its value-type is
`Vertex`.
*/
typedef unspecified_type Vertex_iterator;

/*!
a bidirectional iterator over the vertices. Its value-type is
`Vertex`.
*/
typedef unspecified_type Vertex_const_iterator;


/*!
a bidirectional iterator over the halfedges. Its value-type is
`Halfedge`.
*/
typedef unspecified_type Halfedge_iterator;

/*!
a bidirectional iterator over the halfedges. Its value-type is
`Halfedge`.
*/
typedef unspecified_type Halfedge_const_iterator;


/*!
a bidirectional iterator over the faces. Its value-type is `Face`.
*/
typedef unspecified_type Face_iterator;

/*!
a bidirectional iterator over the faces. Its value-type is `Face`.
*/
typedef unspecified_type Face_const_iterator;


/// @}

/// \name Creation
/// @{

/*!
constructs an empty \dcel with one unbounded face.
*/
Arr_dcel();

/*!
assigns the contents of the `other` \dcel whose unbounded face
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

/*!
returns a begin-iterator of the vertices in `dcel`.
*/
Vertex_iterator vertices_begin();

/*!
returns a past-the-end iterator of the vertices in `dcel`.
*/
Vertex_iterator vertices_end();

/*!
returns a range over handles of the vertices in `dcel`.
*/
unspecified_type vertex_handles();

/*!
returns a begin-iterator of the vertices in `dcel`.
*/
Vertex_const_iterator vertices_begin() const;

/*!
returns a past-the-end iterator of the vertices in `dcel`.
*/
Vertex_const_iterator vertices_end() const;

/*!
returns a const range (model of `ConstRange`) over handles of the vertices in `dcel`.
*/
unspecified_type vertex_handles() const;

/*!
returns a begin-iterator of the halfedges in `dcel`.
*/
Halfedge_iterator halfedges_begin();

/*!
returns a past-the-end iterator of the halfedges in `dcel`.
*/
Halfedge_iterator halfedges_end();

/*!
returns a range over handles of the halfedges in `dcel`.
*/
unspecified_type halfedge_handles();

/*!
returns a begin-iterator of the halfedges in `dcel`.
*/
Halfedge_const_iterator halfedges_begin() const;

/*!
returns a past-the-end iterator of the halfedges in `dcel`.
*/
Halfedge_const_iterator halfedges_end() const;

/*!
returns a const range (model of `ConstRange`) over handles of the halfedges in `dcel`.
*/
unspecified_type halfedge_handles() const;

/*!
returns a begin-iterator of the faces in `dcel`.
*/
Face_iterator faces_begin();

/*!
returns a past-the-end iterator of the faces in `dcel`.
*/
Face_iterator faces_end();

/*!
returns a range over handles of the faces in `dcel`.
*/
unspecified_type face_handles();

/*!
returns a begin-iterator of the faces in `dcel`.
*/
Face_const_iterator faces_begin() const;

/*!
returns a past-the-end iterator of the faces in `dcel`.
*/
Face_const_iterator faces_end() const;

/*!
returns a const range (model of `ConstRange`) over handles of the faces in `dcel`.
*/
unspecified_type face_handles() const;

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

}; /* end ArrangementDcel */

