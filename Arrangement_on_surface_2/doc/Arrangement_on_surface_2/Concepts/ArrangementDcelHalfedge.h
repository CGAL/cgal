
/*!
\ingroup PkgArrangementOnSurface2ConceptsDCEL
\cgalConcept

A halfedge record in a \dcel data structure. Two halfedges with opposite
directions always form an edge (a halfedge pair). The halfedges form together
chains, defining the boundaries of connected components, such that all
halfedges along a chain have the same incident face. Note that the chain the
halfedge belongs to may form the outer boundary of a bounded face (an outer
CCB) or the boundary of a hole inside a face (an inner CCB).

An edge is always associated with a curve, but the halfedge records only
store a pointer to the associated curve, and the actual curve objects
are stored elsewhere. Two opposite halfedges are always associated with
the same curve.

\sa `ArrangementDcel`
\sa `ArrangementDcelVertex`
\sa `ArrangementDcelFace`
\sa `ArrangementDcelHole`

*/

class ArrangementDcelHalfedge {
public:

/// \name Types
/// @{

/*!
the corresponding \dcel vertex type.
*/
typedef unspecified_type Vertex;

/*!
the corresponding \dcel face type.
*/
typedef unspecified_type Face;

/*!
the corresponding \dcel hole type.
*/
typedef unspecified_type Hole;

/*!
the curve type associated with the edge.
*/
typedef unspecified_type X_monotone_curve;

/// @}

/// \name Creation
/// @{

/*!
default constructor.
*/
Arr_dcel_halfedge();

/*!
assigns `e` with the contents of the `other` halfedge.
*/
void assign (const Self& other);

/// @}

/// \name Access Functions
/// @{

/*!
returns `ARR_LEFT_TO_RIGHT` if `e`'s source vertex is
lexicographically smaller than it target, and
`ARR_RIGHT_TO_LEFT` if it is lexicographically larger than
the target.
*/
Arr_halfedge_direction direction() const;

/*!
determines whether the `e` lies on an outer CCB of a bounded face,
or on an inner CCB (a hole inside a face). The function returns `true`
if `e` lies on a hole.
*/
bool is_on_hole() const;

/// @}

/// \name
/// All functions below also have `const` counterparts, returning non-mutable pointers or references:
/// @{

/*!
returns the twin halfedge.
*/
Halfedge* opposite();

/*!
returns the previous halfedge along the chain.
*/
Halfedge* prev();

/*!
returns the next halfedge along the chain.
*/
Halfedge* next();

/*!
returns the target vertex.
*/
Vertex* vertex();

/*!
returns the incident face.
\pre `e` lies on the outer boundary of this face.
*/
Face* face();

/*!
returns the hole (inner CCB) `e` belongs to.
\pre `e` lies on a hole inside its incident face.
*/
Hole* hole();

/*!
returns whether the vertex is not associated with a valid curve.
*/
bool has_null_curve() const;

/*!
returns the associated curve.
\pre `e` is associated with a valid curve.
*/
X_monotone_curve& curve();

/// @}

/// \name Modifiers
/// @{

/*!
sets the opposite halfedge.
*/
void set_opposite (Halfedge* opp);

/*!
sets the lexicographical order between `e`'s source and target
vertices to be `dir`.
The direction of the opposite halfedge is also set to the
opposite direction.
*/
void set_direction (Arr_halfedge_direction dir);

/*!
sets the previous halfedge of `e` along the chain,
and updates the cross-pointer `prev->next()`.
*/
void set_prev (Halfedge* prev);

/*!
sets the next halfedge of `e` along the chain,
and updates the cross-pointer `next->prev()`.
*/
void set_next (Halfedge* next);

/*!
sets the target vertex.
*/
void set_vertex (Vertex* v);

/*!
sets the incident face, marking that `e` lies on the outer CCB
of the face `f`.
*/
void set_face (Face* f);

/*!
sets the incident hole, marking that `e` lies on an inner CCB.
*/
void set_hole (Hole* ho);

/*!
sets the associated curve of `e` and its opposite halfedge.
*/
void set_curve (X_monotone_curve* c);

/// @}

}; /* end ArrangementDcelHalfedge */

