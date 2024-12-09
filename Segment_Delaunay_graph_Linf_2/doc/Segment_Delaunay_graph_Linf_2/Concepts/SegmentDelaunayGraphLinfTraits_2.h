/*!
\ingroup PkgSegmentDelaunayGraphLinf2Concepts
\cgalConcept

The concept `SegmentDelaunayGraphLinfTraits_2` provides traits
for constructing the segment Delaunay graph under the
\f$ L_{\infty} \f$ distance.
The segment Delaunay graph is the dual of the segment Voronoi diagram.
We stress that we consider the 1-dimensionalization of \f$ L_{\infty} \f$
bisectors between two sites which is explained in
Section \ref subsecbis1dim of the User Manual,
and this reflects on the constructed graph (and
its dual diagram).
These traits should be used in the Gt template parameter of the
`CGAL::Segment_Delaunay_graph_Linf_2<Gt,DS>` and
`CGAL::Segment_Delaunay_graph_Linf_hierarchy_2<Gt,STag,DS>` class templates.
The concept is a refinement of `SegmentDelaunayGraphTraits_2`.
In particular, it provides a type `Site_2`, which must be a model of
the concept `SegmentDelaunayGraphSite_2`. It also provides
constructions for sites and several function object
types for the predicates.

In contrast with \f$ L_{2} \f$, the concept also contains drawing
methods for the edges of the \f$ L_{\infty} \f$ segment Voronoi diagram
(class templates
`Construct_sdg_bisector_2`,
`Construct_sdg_bisector_ray_2`, and
`Construct_sdg_bisector_segment_2`).
These methods are used when, additionally, the tag
type `Has_bisector_constructions_type` is defined in the concept.

We only describe the refined and additional requirements
of the `SegmentDelaunayGraphLinfTraits_2` concept
with respect to the
`SegmentDelaunayGraphTraits_2` concept.

\cgalRefines{SegmentDelaunayGraphTraits_2}

\cgalHasModelsBegin
\cgalHasModels{CGAL::Segment_Delaunay_graph_Linf_traits_2<K,MTag>}
\cgalHasModels{CGAL::Segment_Delaunay_graph_Linf_traits_without_intersections_2<K,MTag>}
\cgalHasModels{CGAL::Segment_Delaunay_graph_Linf_filtered_traits_2<CK,CM,EK,EM,FK,FM>}
\cgalHasModels{CGAL::Segment_Delaunay_graph_Linf_filtered_traits_without_intersections_2<CK,CM,EK,EM,FK,FM>}
\cgalHasModelsEnd

\sa `SegmentDelaunayGraphSite_2`
\sa `CGAL::Segment_Delaunay_graph_Linf_2<Gt,DS>`
\sa `CGAL::Segment_Delaunay_graph_Linf_hierarchy_2<Gt,STag,DS>`

*/

class SegmentDelaunayGraphLinfTraits_2 {
public:

/// \name Bisector construction tag
/// @{

/*!
If this tag type is defined in the concept, then
the concept should also contain three bisector drawing
class templates, `Construct_sdg_bisector_2`,
`Construct_sdg_bisector_ray_2`, and
`Construct_sdg_bisector_segment_2`, that will be used when drawing the
dual of the \f$L_{\infty}\f$ segment Delaunay graph.
This is the only way to draw correctly the \f$L_{\infty}\f$
segment Voronoi diagram. If this type is omitted, then the default
drawing methods are used, which are only good for the \f$L_2\f$
segment Voronoi diagram.
*/
  typedef char Has_bisector_constructions_type;

/// @}


/// \name Types
/// @{

/*!
A constructor for the point \f$ v \f$ at which the three \f$ L_{\infty} \f$
bisectors between the three given sites `s1`, `s2` and `s3` intersect
and the regions of `s1`, `s2` and `s3` appear in the counter-clockwise order
`s1`, `s2`, `s3` around \f$ v\f$.

Point \f$ v \f$ is equidistant
from the three sites, under the \f$ L_{\infty} \f$ distance.

It must provide
`Point_2 operator()(Site_2 s1, Site_2 s2, Site_2 s3)`, which
constructs this point \f$ v \f$ which is
equidistant (under the \f$ L_{\infty} \f$ distance)
from the sites `s1`, `s2` and `s3`.
*/
typedef Hidden_type Construct_svd_vertex_2;

/*!
A predicate object type.
It must provide `Oriented_side operator()(Site_2 s1, Site_2 s2, Point_2 p)`,
which returns
the oriented side of the \f$ L_{\infty} \f$ bisector of `s1` and `s2` that
contains `p`. It returns `ON_POSITIVE_SIDE` if `p` lies in
the nearest region of `s1` (i.e., `p` is closer to `s1` than
`s2`); returns `ON_NEGATIVE_SIDE` if `p` lies in the nearest
region of `s2`; returns `ON_ORIENTED_BOUNDARY` if `p`
lies on the \f$ L_{\infty} \f$ bisector of `s1` and `s2`.
*/
typedef Hidden_type Oriented_side_of_bisector_2;

/*!
A predicate object type.

It must provide `Sign operator()(Site_2 s1, Site_2 s2, Site_2 s3, Site_2 q)`,
which
returns the sign of the distance of `q` from the
\f$ L_{\infty} \f$ Voronoi square
of `s1`, `s2`, `s3`. The \f$ L_{\infty} \f$ Voronoi square
of three sites
`s1`, `s2`, `s3` is an axis-parallel square which is passing through all three
sites and touches them in the `s1`, `s2`, `s3`
order as we walk on the square in the counter-clockwise sense.
The center of the square is at the intersection of the three
\f$ L_{\infty} \f$ bisectors of the three sites.

\pre the \f$ L_{\infty} \f$ Voronoi square of `s1`, `s2`, `s3` must exist.

It must also provide `Sign operator()(Site_2 s1, Site_2 s2, Site_2 q)`,
which returns the sign of the distance of
`q` from the
\f$ L_{\infty} \f$ Voronoi square of sites `s1`, `s2`, \f$s_\infty\f$,
where \f$s_\infty\f$ is the dummy site at infinity.
This is a degenerate
\f$ L_{\infty} \f$ Voronoi square, with its center at infinity, which
is either a line or a right angle wedge.
*/
typedef Hidden_type Vertex_conflict_2;

/*!
A predicate object type.

It must provide
`bool operator()(Site_2 s1, Site_2 s2, Site_2 s3, Site_2 s4, Site_2 q, Sign sgn)`.
The sites `s1`, `s2`,
`s3` and `s4` define a Voronoi edge that lies on the
bisector of `s1` and `s2` and has as endpoints the \f$L_{\infty}\f$ Voronoi
vertices
\f$ v_{123} \f$ and \f$ v_{142} \f$
defined by the ordered triplets `s1`, `s2`, `s3` and
`s1`, `s4`, `s2`, respectively. The sign `sgn` is the common sign
of the distance of the site `q` from the \f$L_{\infty}\f$ Voronoi square of the
triplets `s1`, `s2`, `s3` and `s1`, `s4`,
`s2`. In case that `sgn` is equal to `NEGATIVE`, the
predicate returns `true` if and only if the entire Voronoi edge is
in conflict with `q`. If `sgn` is equal to `POSITIVE` or
`ZERO`, the predicate returns `false` if and only if `q`
is not in conflict with the Voronoi edge.
\pre
the \f$L_{\infty}\f$ Voronoi vertices
\f$ v_{123} \f$ and \f$ v_{142} \f$
must exist
and the sign of the distance of `q` from these two vertices
must be common and equal to `sgn`.

It must also provide
`bool operator()(Site_2 s1, Site_2 s2, Site_2 s3, Site_2 q, Sign sgn)`.
The
sites `s1`, `s2`, `s3` and the dummy site at infinity
\f$ s_\infty\f$ define a Voronoi edge that lies on the \f$L_{\infty}\f$
bisector of
`s1` and `s2` and has as endpoints the Voronoi vertices
\f$ v_{123}\f$ and \f$ v_{{1}\infty{2}}\f$ defined by the triplets `s1`,
`s2`, `s3` and `s1`, \f$ s_\infty\f$,  `s2` (the second
vertex is at infinity). The sign `sgn` is the common sign
of the distance of the site `q` from the two \f$L_{\infty}\f$ Voronoi squares
centered at the Voronoi vertices \f$ v_{123}\f$ and \f$ v_{1\infty{2}}\f$.
In case that `sgn` is `NEGATIVE`, the predicate
returns `true` if and only if the entire Voronoi edge is in
conflict with `q`. If `sgn` is `POSITIVE` or `ZERO`,
the predicate returns `false` if and only if `q` is not in
conflict with the Voronoi edge.
\pre the \f$L_{\infty}\f$
Voronoi vertex \f$ v_{123}\f$ of `s1`, `s2`, `s3` must exist
and the sign of the distance of `q` from the vertices \f$ v_{123}\f$
and \f$ v_{1{\infty}2}\f$
must be common and equal to `sgn`.

It must finally provide
`bool operator()(Site_2 s1, Site_2 s2, Site_2 q, Sign sgn)`.
The
sites `s1`, `s2` and the dummy site at infinity
\f$ s_\infty\f$ define a Voronoi edge that is equal to the
\f$L_{\infty}\f$ bisector of `s1` and `s2`.
The endpoints of this edge are the \f$L_{\infty}\f$ Voronoi
vertices
\f$ v_{12\infty}\f$ and \f$ v_{1\infty{}2}\f$, defined
by the triplets `s1`, `s2`, \f$ s_\infty\f$ and `s1`,
\f$ s_\infty\f$, `s2` (both vertices are at
infinity).
The sign `sgn` denotes the common sign of the distance
of the site `q` from the \f$L_{\infty}\f$ Voronoi squares centered at
\f$ v_{12\infty}\f$ and \f$ v_{1\infty{}2}\f$.
If `sgn` is `NEGATIVE`, the predicate
returns `true` if and only if the entire Voronoi edge is in
conflict with `q`. If `sgn` is `POSITIVE` or `ZERO`,
the predicate returns `false` if and only if `q` is not in
conflict with the Voronoi edge.
\pre
the sign of the distance of `q` from the \f$L_{\infty}\f$ Voronoi vertices
\f$ v_{12\infty}\f$ and \f$ v_{1\infty{}2}\f$
must be common and equal to `sgn`.
*/
typedef Hidden_type Finite_edge_interior_conflict_2;

/*!
A predicate object type.

It must provide
`bool operator()(Site_2 s1, Site_2 s2, Site_2 s3, Site_2 q, Sign sgn)`.
Let \f$s_\infty\f$ be the dummy site at infinity.
The
sites \f$ s_\infty\f$, `s1`, `s2` and `s3` define a
Voronoi edge that lies on the bisector of \f$ s_\infty\f$ and `s1`
and has as endpoints the \f$L_{\infty}\f$
Voronoi vertices \f$ v_{\infty{}12}\f$ and
\f$ v_{\infty{}31}\f$ defined by the triplets
\f$ s_\infty\f$, `s1`, `s2` and \f$ s_\infty\f$, `s3`,
`s1`, respectively.
The sign `sgn` is the common sign of the distances of
`q` from the \f$L_{\infty}\f$ Voronoi squares centered at the vertices
\f$ v_{\infty{}{1}{2}}\f$ and \f$ v_{\infty{}31}\f$. If `sgn` is `NEGATIVE`,
the predicate returns `true` if and only if the entire Voronoi
edge is in conflict with `q`. If `sgn` is `POSITIVE` or
`ZERO`, the predicate returns `false` if and only if `q`
is not in conflict with the Voronoi edge.
\pre
the sign of the distance of `q` from the \f$L_{\infty}\f$ Voronoi vertices
\f$ v_{{\infty}{1}{2}}\f$ and \f$ v_{{\infty}{3}{1}}\f$
must be common and equal to `sgn`.
*/
typedef Hidden_type Infinite_edge_interior_conflict_2;

/*!
A predicate object type.

First, we define the notion of (non-oriented) \f$L_{\infty}\f$-perpendicular lines to
a given non-trivial (non-oriented) segment \f$ s \f$.
If \f$ s \f$ is horizontal, then the perpendicular lines are the
vertical lines.
If \f$ s \f$ is vertical, then the perpendicular lines are the
horizontal lines.
If \f$ s \f$ has positive slope, then the perpendicular lines are the
lines of slope -1.
If \f$ s \f$ has negative slope, then the perpendicular lines are the
lines of slope +1.

Since CGAL segments have also an orientation, we also orient
\f$L_{\infty}\f$-perpendicular lines, as follows.
For an <I>oriented</I> segment \f$ s \f$, we orient its
\f$L_{\infty}\f$-perpendicular lines so that the lines'
orientation is closest to the following orientation:
the orientation of \f$ s \f$ rotated <I>counter-clockwise</I> by
\f$ \pi/2 \f$.

Let `s` be a segment and `p` a point contained in its interior.
Let \f$ \ell\f$ be the line which is \f$L_{\infty}\f$-perpendicular
to segment `s` and is passing through point `p`.

The predicate object type must provide
`Oriented_side operator()(Site_1 s1, Site_2 s2, Site_2 s3, Site_2 s, Site_2 p)`.
This
determines the oriented side of the line \f$ \ell\f$
in which the \f$L_{\infty}\f$ Voronoi vertex \f$v_{123}\f$ of the sites `s1`,
`s2`, `s3` is contained.
\pre `s` must be a segment, `p` must be a point, and `p` must
be contained in the interior of `s`.

The predicate object type must also provide
`bool operator()(Site_2 s1, Site_2 s2, Site_2 s, Site_2 p)`.
This
determines the oriented side of the line \f$ \ell\f$
in which the \f$L_{\infty}\f$ Voronoi vertex \f$v_{12\infty}\f$
of the sites `s1`,
`s2`, \f$s_{\infty}\f$ is contained.
\pre `s` must be a segment, `p` must be a point, and `p` must
be contained in the interior of `s`.
*/
typedef Hidden_type Oriented_side_2;

/// @}

/// \name Access to predicate objects
/// @{


/*!

*/
Oriented_side_of_bisector_2
oriented_side_of_bisector_test_2_object();

/*!

*/
Vertex_conflict_2 vertex_conflict_2_object();

/*!

*/
Finite_edge_interior_conflict_2
finite_edge_interior_conflict_2_object();

/*!

*/
Infinite_edge_interior_conflict_2
infinite_edge_interior_conflict_2_object();

/*!

*/
Oriented_side_2 oriented_side_2_object();


/// @}

/// \name Access to constructor objects
/// @{

/*!

*/
Construct_svd_vertex_2
construct_svd_vertex_2_object();

/// @}



/// \name Bisector construction class templates
/// @{


/*!
The class template drawing the
\f$L_{\infty}\f$ bisector between two sites.

The class should define a
`Bisector_type operator()(Site_2 p, Site_2 q)`
that returns the  \f$L_{\infty}\f$ bisector between
sites `p` and `q`.
*/
template<class Gt, class M>
class Construct_sdg_bisector_2
{};

/*!
The class template drawing the \f$L_{\infty}\f$
edge between two sites, that is bounded by another site and
the dummy site \f$s_{\infty}\f$ (at infinity).

The class should define a
`Bisector_ray_type operator()(Site_2 p, Site_2 q, Site_2 r)`
that returns the edge between sites `p` and `q` that is bounded
by the \f$L_{\infty}\f$ Voronoi vertices \f$v_{pqr}\f$ and
\f$v_{qp{\infty}}\f$.
*/
template<class Gt, class M>
class Construct_sdg_bisector_ray_2
{};

/*!
The class template drawing the \f$L_{\infty}\f$
edge between two sites, that is bounded by two other sites.

The class should define a
`Bisector_segment_type operator()(Site_2 p, Site_2 q, Site_2 r, Site_2 s)`
that returns the edge between sites `p` and `q` that is bounded
by the \f$L_{\infty}\f$ Voronoi vertices \f$v_{pqr}\f$ and
\f$v_{qps}\f$.
*/
template<class Gt, class M>
class Construct_sdg_bisector_segment_2
{};

/// @}


}; /* end SegmentDelaunayGraphTraits_2 */


