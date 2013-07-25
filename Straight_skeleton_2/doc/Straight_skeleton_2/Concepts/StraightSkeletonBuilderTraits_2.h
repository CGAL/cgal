/*!
\ingroup PkgStraightSkeleton2Concepts
\cgalConcept

The concept `StraightSkeletonBuilderTraits_2` describes the requirements for the geometric traits class required by the algorithm class `Straight_skeleton_builder_2<Gt,Ss>`. 

\cgalHasModel CGAL::Straight_skeleton_builder_traits_2

\sa `CGAL::Straight_skeleton_builder_2<Gt,Ss>` 
\sa `CGAL::Straight_skeleton_builder_traits_2<K>` 

*/

class StraightSkeletonBuilderTraits_2 {
public:

/// \name Types 
/// @{

/*!
A model of the `Kernel` concept. 
*/ 
typedef unspecified_type Kernel; 

/*!
A model of the `FieldWithSqrt` concept provided by the kernel. This type is used to represent the 
coordinates of the input points, of the skeleton nodes and as the event time stored in the skeleton nodes. 
*/ 
typedef unspecified_type FT; 

/*!
A pair of (x,y) coordinates representing a 2D Cartesian point. 
*/ 
boost::tuple<FT,FT> Vertex; 

/*!
A pair of vertices representing an edge 
*/ 
boost::tuple<Vertex,Vertex> Edge; 

/*!
A triple of edges representing an event 
*/ 
boost::tuple<Edge,Edge,Edge> EdgeTriple; 

/*!
A predicate object type being a model of the `Kernel::Equal_2` function object concept. 
*/ 
typedef unspecified_type Equal_2; 

/*!
A predicate object type being a model of the `Kernel::LeftTurn_2` function object concept. 
*/ 
typedef unspecified_type Left_turn_2; 

/*!
A predicate object type being a model of the `Kernel::Collinear_2` function object concept. 
*/ 
typedef unspecified_type Collinear_2; 

/*!
A predicate object type. 

Must provide `bool operator()( EdgeTriple const& et) const`, which determines if, given the 3 <I>oriented</I> lines defined by the 3 input edges (3 pair of points), there exist an Euclidean distance `t >= 0` for which the corresponding 3 <I>offset lines at `t`</I> (parallel lines at an Euclidean distance of `t`) intersect in a single point. 

\pre each edge in the triple must properly define an oriented line, that is, such points cannot be coincident. 
*/ 
typedef unspecified_type Do_ss_event_exist_2; 

/*!
A predicate object type. 

Must provide `Comparison_result operator()( EdgeTriple const& x, EdgeTriple const& y) const`, which compares the `times` for the events determined by the edge-triples `x` and `y`. 

The time of an event given by an edge triple (which defines 3 oriented lines) is the Euclidean distance `t` at which the corresponding offset lines at `t` intersect in a single point. 

\pre `x` and `y` must be edge-triples corresponding to events that actually exist (as determined by the predicate `Exist_sls_event_2`). 
*/ 
typedef unspecified_type Compare_ss_event_times_2; 

/*!
A predicate object type. 

Must provide `Comparison_result operator()( Point_2 const& p, EdgeTriple const& x, EdgeTriple const& y)`, which compares the Euclidean distance of the points of event for `x` and `y` to the point `p`. 

The point of an event given by an edge triple (which defines 3 oriented lines) is the intersection point of the 3 corresponding offset lines at the time `t` of the event. 

It must also provide `Comparison_result operator( EdgeTriple const& seed, EdgeTriple const& x, EdgeTriple const& y)()` which makes the same comparison as the first overload but where the seed point is given implicitly as the point of event for `seed`. 

\pre `seed`, `x` and `y` must be edge-triples corresponding to events that actually exist (as determined by the predicate `Exist_sls_event_2`). 
*/ 
typedef unspecified_type Compare_ss_event_distance_to_seed_2; 

/*!
A predicate object type. 

Must provide `bool operator()( EdgeTriple const& e, EdgeTriple const& zone)`, which determines if the point of event for `e` is inside the <I>offset zone</I> defined by the 3 <I>oriented</I> lines given by `zone`. 

An offset zone given by 3 <I>oriented</I> lines is the intersection of the halfplanes to the left of each oriented line. 

\pre `e` must be an edge-triple corresponding to an event that actually exist (as determined by the predicate `Exist_sls_event_2`), and the 3 oriented lines given by `zone` must be well defined (no point-pair can have coincident points). 
*/ 
typedef unspecified_type Is_ss_event_inside_offset_zone_2; 

/*!
A predicate object type. 

Must provide `bool operator()( EdgeTriple const& x, EdgeTriple const& y)`, which determines if the events given by `x` and `y` are coincident in time and space; that is, both triples of offset lines intersect at the same point and at the same Euclidean distance from their sources. 

\pre `x` and `y` must be edge-triples corresponding to events that actually exist (as determined by the predicate `Exist_sls_event_2`). 
*/ 
typedef unspecified_type Are_ss_events_simultaneous_2; 

/*!
A construction object type. 

Must provide `boost::tuple< boost::optional<FT>, boost::optional<Point_2> > operator()( EdgeTriple const& e)`, which given the 3 <I>oriented</I> lines defined by the 3 input edges (3 pair of points), returns the Euclidean distance `t >= 0` and intersection point at which the corresponding 3 <I>offset lines at `t`</I> intersect. 

If the values cannot be computed, not even approximately (because of overflow for instance), an empty optional must be returned. 

\pre `e` must be an edge-triple corresponding to an event that actually exist (as determined by the predicate `Exist_sls_event_2`). 
*/ 
typedef unspecified_type Construct_ss_event_time_and_point_2; 

/*!
A construction object type. 

Must provide `Vertex operator()( Point_2 const& p)`, which given a `Point_2` `p` returns a Vertex encapsulating the corresponding (x,y) pair of <I>Cartesian</I> coordinates. 
*/ 
typedef unspecified_type Construct_ss_vertex_2; 

/*!
A construction object type. 

Must provide `Edge operator()( Point_2 const& s, Point_2 const& t)`, which given source and target points `s` and `t` returns an Edge encapsulating the corresponding input segment (in <I>Cartesian</I> coordinates.) 
*/ 
typedef unspecified_type Construct_ss_edge_2; 

/*!
A construction object type. 

Must provide `Triedge operator()( Edge const& e0, Edge const& e1, Edge const& e2)`, which given the 3 edges that define an event, `e0`, `e1` and `e2`, returns a Triedge encapsulating them. 
*/ 
typedef unspecified_type Construct_ss_triedge_2; 

/// @}

}; /* end StraightSkeletonBuilderTraits_2 */
