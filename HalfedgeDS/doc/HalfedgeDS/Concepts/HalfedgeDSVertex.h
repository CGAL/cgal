/*!
\ingroup PkgHalfedgeDSConcepts
\cgalConcept

The concept `HalfedgeDSVertex` defines the requirements for the local `Vertex`
type in the `HalfedgeDS` concept. It is also required in
the `Vertex_wrapper<Refs,Traits>` member class template of an
items class, see the `HalfedgeDSItems` concept.

A vertex optionally stores a reference to an incident halfedge that
points to the vertex. A type tag indicates whether the related member
functions are supported.
Figure \ref figureHalfedgeDSOptionalMethods
depicts the relationship between a halfedge and its incident
halfedges, vertices, and faces.

For the protection of the integrity of the data structure classes such as
`CGAL::Polyhedron_3` are allowed to redefine the modifying member
functions to be private. In order to make them accessible for the
halfedge data structure they must be derived from a base class `Base`
where the modifying member functions are still public. (The protection
could be bypassed by an user, but not by accident.)

\cgalRefines{CopyConstructible,DefaultConstructible}

\cgalHasModelsBegin
\cgalHasModels{CGAL::HalfedgeDS_vertex_base<Refs>}
\cgalHasModels{CGAL::HalfedgeDS_vertex_min_base<Refs>}
\cgalHasModelsEnd

\sa `HalfedgeDS<Traits,Items,Alloc>`
\sa `HalfedgeDSItems`
\sa `HalfedgeDSHalfedge`
\sa `HalfedgeDSFace`

*/

class HalfedgeDSVertex {
public:

/// \name Types
/// @{

/*!
instantiated `HalfedgeDS` ( \f$ \equiv\f$ `Refs`).
*/
typedef unspecified_type HalfedgeDS;

/*!
base class that allows modifications.
*/
typedef unspecified_type Base;

/*!
model of `HalfedgeDSHalfedge`.
*/
typedef unspecified_type Halfedge;

/*!
model of `HalfedgeDSFace`.
*/
typedef unspecified_type Face;

/*!
handle to vertex.
*/
typedef unspecified_type Vertex_handle;

/*!
handle to halfedge.
*/
typedef unspecified_type Halfedge_handle;

/*!
handle to face.
*/
typedef unspecified_type Face_handle;

/*!

*/
typedef unspecified_type Vertex_const_handle;

/*!

*/
typedef unspecified_type Halfedge_const_handle;

/*!

*/
typedef unspecified_type Face_const_handle;

/*!
`CGAL::Tag_true` or
`CGAL::Tag_false`.
*/
typedef unspecified_type Supports_vertex_halfedge;

/// @}

/// \name Creation
/// @{

/*!
default constructor.
*/
Vertex();

/// @}

/// \name Operations available if Supports_vertex_halfedge == CGAL::Tag_true
/// @{

/*!

*/
Halfedge_handle halfedge();
/*!

incident halfedge that points to the vertex.
*/
Halfedge_const_handle halfedge() const;
/*!

sets incident halfedge to `h`.
*/
void set_halfedge( Halfedge_handle h);

/// @}

}; /* end HalfedgeDSVertex */
