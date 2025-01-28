/*!
\ingroup PkgStraightSkeleton2Concepts
\cgalConcept

\cgalRefines{HalfedgeDSVertex}

The concept `StraightSkeletonVertex_2` describes the requirements for the vertex type of the
`StraightSkeleton_2` concept. It is a refinement of the `HalfedgeDSVertex` concept
with support for storage of the incident halfedge. The `StraightSkeletonVertex_2` concept
requires the geometric embedding to be a 2D point.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Straight_skeleton_vertex_base_2}
\cgalHasModelsEnd

\sa `StraightSkeletonHalfedge_2`
\sa `StraightSkeletonFace_2`
*/
class StraightSkeletonVertex_2 {
public:

/// \name Types
/// @{

/*!
The type of the 2D point being the geometric embedding of the vertex
*/
typedef unspecified_type Point_2;

/*!
A model of the `FieldWithSqrt` concept representing the time of a vertex (a Euclidean distance)
*/
typedef unspecified_type FT;

/*!

*/
typedef unspecified_type Halfedge_around_vertex_const_circulator;

/*!
The circulator type used to visit all the incident halfedges around a vertex
*/
typedef unspecified_type Halfedge_around_vertex_circulator;

/*!

*/
typedef unspecified_type Defining_contour_halfedge_const_circulator;

/*!
The circulator type used to visit all the defining contour halfedges of a vertex
*/
typedef unspecified_type Defining_contour_halfedge_circulator;

/// @}

/// \name Creation
/// @{

/*!
%Default constructor
*/
StraightSkeletonVertex_2();

/*!
constructs a contour vertex with ID number `id`, at the point `p`
*/
StraightSkeletonVertex_2(int id, const Point_2& p);

/*!
constructs a skeleton vertex with ID number `id`, at point `p` and time `time`.
*/
StraightSkeletonVertex_2(int id, const Point_2& p, FT time );

/// @}

/// \name Access Functions
/// @{

/*!
The ID of the vertex
*/
int id() const;

/*!
The vertex point
*/
const Point_2& point() const;

/*!
The time of the vertex: the distance from the vertex point to the lines supporting the defining contour edges
*/
FT time() const;

/*!

*/
Halfedge_handle primary_bisector();

/*!
returns the skeleton halfedge incident upon the vertex (called the <I>primary</I> bisector).
*/
Halfedge_const_handle primary_bisector() const;

/*!

*/
Halfedge_around_vertex_circulator halfedge_around_vertex_begin();

/*!
returns a bi-directional circulator pointing to one of the incident halfedges (which one is unspecified).

There will always be as many incident halfedges as the degree of the vertex.

If this is a <I>contour</I> vertex, its degree is exactly 3, and from the halfedges pointed to by the circulator, 2 are contour and 1 is a bisector.

If this is an `skeleton` vertex, its degree is at least 3 and all of the halfedges pointed to by the circulator are bisectors.

Each halfedge pointed to by this circulator is the one which is oriented towards the vertex (according to the geometric embedding).
*/
Halfedge_around_vertex_const_circulator halfedge_around_vertex_begin() const;

/*!

*/
Defining_contour_halfedge_circulator defining_contour_halfedges_begin();

/*!
returns a bi-directional circulator pointing to one of the defining contour halfedges of the vertex (which one is unspecified).

There will always be as many incident defining contour halfedges as the degree of the vertex.

Each halfedge pointed to by this circulator is the one having its left side facing inwards
(which happens to be the contour halfedge for which `is_border()` is `false`).
*/
Defining_contour_halfedge_const_circulator defining_contour_halfedges_begin() const;

/// @}

/// \name Queries
/// @{

/*!
returns `true` iff this is a contour vertex.
*/
bool is_contour() const;

/*!
returns `true` iff this is a skeleton vertex.
*/
bool is_skeleton() const;

/// @}

}; /* end StraightSkeletonVertex_2 */
