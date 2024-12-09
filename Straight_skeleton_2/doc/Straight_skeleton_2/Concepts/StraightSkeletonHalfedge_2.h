/*!
\ingroup PkgStraightSkeleton2Concepts
\cgalConcept

\cgalRefines{HalfedgeDSHalfedge}

The concept `StraightSkeletonHalfedge_2` describes the requirements for the halfedge type of the
`StraightSkeleton_2` concept. It is a refinement of the `HalfedgeDSHalfedge` concept.

The `StraightSkeletonHalfedge_2` concept requires no geometric embedding at all.
The only geometric embedding used by the Straight Skeleton Data Structure are the 2D points
in the contour and skeleton vertices. However, for any halfedge, there is a 2D segment implicitly
given by its `source` and `target` vertices.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Straight_skeleton_halfedge_base_2}
\cgalHasModelsEnd

\sa `StraightSkeleton_2`
\sa `StraightSkeletonHalfedge_2`
\sa `CGAL::Straight_skeleton_vertex_base_2<Refs,Point,FT>`
\sa `CGAL::Straight_skeleton_halfedge_base_2<Refs>`
*/
class StraightSkeletonHalfedge_2 {
public:

/// \name Creation
/// @{

/*!
%Default Constructor.
*/
StraightSkeletonHalfedge_2();

/*!
constructs a halfedge with ID `id`.

It is the links to other halfedges that determines if this is a contour edge,
a contour-skeleton edge, or an inner-skeleton edge.
*/
StraightSkeletonHalfedge_2( int id );

/// @}

/// \name Access Functions
/// @{

/*!
returns the ID of the halfedge
*/
int id() const;

/*!
resets the ID of the halfedge to `aID`
*/
void reset_id ( int aID );

/*!
returns the weight of the halfedge
*/
FT weight() const;

/*!
sets the weight of the halfedge to `aWeight`
*/
void set_weight( FT aWeight );

/*!

*/
Halfedge_handle defining_contour_edge();

/*!
If this is a bisector halfedge, returns a handle to the inward-facing (non-border) contour halfedge
corresponding to the defining contour edge which is to its left; if this is a contour halfedge,
returns a handle to itself if `is_border()` is `false`, or to its opposite if it is `true`.
*/
Halfedge_const_handle defining_contour_edge() const;

/// @}

/// \name Queries
/// @{

/*!
returns `true` iff this is a bisector (or skeleton) halfedge (i.e.\ is not a contour halfedge).
*/
bool is_bisector() const;

/*!
returns `true` iff this is a bisector and is inner (i.e.\ is not incident upon a contour vertex).
*/
bool is_inner_bisector() const;

/// @}

}; /* end StraightSkeletonHalfedge_2 */
