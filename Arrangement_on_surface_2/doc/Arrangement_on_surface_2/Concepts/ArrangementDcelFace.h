
/*!
\ingroup PkgArrangementOnSurface2ConceptsDCEL
\cgalConcept

A face record in a \dcel data structure. A face may either be unbounded,
otherwise it has an incident halfedge along the chain defining its outer
boundary. A face may also contain holes and isolated vertices in its
interior.

\sa `ArrangementDcel`
\sa `ArrangementDcelVertex`
\sa `ArrangementDcelHalfedge`

*/

class ArrangementDcelFace {
public:

/// \name Types
/// The non-mutable iterators `Hole_const_iterator`, and `Isolated_vertex_const_iterator` are also defined.
/// @{

/*!
the corresponding \dcel vertex type.
*/
typedef unspecified_type Vertex;

/*!
the corresponding \dcel halfedge type.
*/
typedef unspecified_type Halfedge;

/*!
a bidirectional iterator over the holes in
inside the face. Its value-type is `Halfedge*`.
*/
typedef unspecified_type Hole_iterator;

/*!
a bidirectional iterator over the
isolated vertices in inside the face. Its value-type is `Vertex*`.
*/
typedef unspecified_type Isolated_vertex_iterator;

/// @}

/// \name Creation
/// @{

/*!
default constructor.
*/
Arr_dcel_face();

/*!
assigns `f` with the contents of the `other` face.
*/
void assign (const Self& other);

/// @}

/// \name Access Functions
/// All functions below also have `const` counterparts, returning
/// non-mutable pointers or iterators:
/// @{

/*!
returns whether the face is unbounded.
*/
bool is_unbounded() const;

/*!
returns an incident halfedge along the outer boundary of the face.
If `f` has no outer boundary, the function returns `nullptr`.
*/
Halfedge* halfedge();

/*!
returns the number of holes inside `f`.
*/
size_t number_of_holes() const;

/*!
returns a begin-iterator for the holes inside `f`.
*/
Hole_iterator holes_begin();

/*!
returns a past-the-end iterator for the holes inside `f`.
*/
Hole_iterator holes_end();

/*!
returns the number of isolated vertices inside `f`.
*/
size_t number_of_isolated_vertices() const;

/*!
returns a begin-iterator for the isolated vertices inside `f`.
*/
Isolated_vertex_iterator isolated_vertices_begin();

/*!
returns a past-the-end iterator for the isolated vertices inside
`f`.
*/
Isolated_vertex_iterator isolated_vertices_end();

/// @}

/// \name Modifiers
/// @{

/*!
sets the face as unbounded (if `flag` is `true`), or as a bounded
face (if it is `false`).
*/
void set_unbounded (bool flag);

/*!
sets the incident halfedge.
*/
void set_halfedge (Halfedge* e);

/*!
adds `e` as a hole inside `f`.
*/
void add_hole (Halfedge* e);

/*!
removes the hole that `it` points to from inside `f`.
*/
void erase_hole (Hole_iterator it);

/*!
adds `v` as an isolated vertex inside `f`.
*/
void add_isolated_vertex (Vertex* v);

/*!
removes the isolated vertex that `it` points to from inside `f`.
*/
void erase_isolated_vertex (Isolated_vertex_iterator it);

/// @}

}; /* end ArrangementDcelFace */

