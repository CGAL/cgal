/*!
\ingroup PkgStraightSkeleton2Concepts
\cgalConcept

\cgalRefines `DefaultConstructible`
\cgalRefines `CopyConstructible`

The concept `StraightSkeletonBuilderTraits_2` describes the requirements
for the geometric traits class required by the algorithm class `CGAL::Straight_skeleton_builder_2`.

\cgalHasModel `CGAL::Straight_skeleton_builder_traits_2<K>`
*/
class StraightSkeletonBuilderTraits_2 {
public:

/// \name Kernel Types
/// @{

/*!
A model of the `Kernel` concept.
*/
typedef unspecified_type Kernel;

/*!
A model of the `FieldWithSqrt` concept provided by the kernel.
This type is used to represent the event time stored in the skeleton nodes.
*/
typedef unspecified_type FT;

/*!
A model of the `Kernel::Point_2` concept
*/
typedef unspecified_type Point_2;

/*!
A model of the `Kernel::Segment_2` concept
*/
typedef unspecified_type Segment_2;

/*!
A model of the `Kernel::Vector_2` concept
*/
typedef unspecified_type Vector_2;

/*!
A model of the `Kernel::Direction_2` concept
*/
typedef unspecified_type Direction_2;

/*!
*/
typedef CGAL::Trisegment_2<Kernel> Trisegment_2;

/*!
*/
typedef boost::intrusive_ptr<Trisegment_2> Trisegment_2_ptr;

/// @}

/// \name Predicates and Constructors Types
/// @{

/*!
A predicate object type.

Must provide `bool operator()( const Trisegment_2_ptr& tri_segment, boost::optional<FT> max_time ) const`,
which determines if, given the three <I>oriented</I> lines defined by the three input edges,
there exists an Euclidean distance `t >= 0` and `t <= max_time` for which the corresponding three
<I>offset lines at `t`</I> (parallel lines at an Euclidean distance of `t`) intersect in a single point.

\pre Each edge in the triple must properly define an oriented line, that is, its points cannot be coincident.
*/
typedef unspecified_type Do_ss_event_exist_2;

/*!
A predicate object type.

Must provide `CGAL::Comparison_result operator()(const Trisegment_2_ptr& l, const Trisegment_2_ptr& r) const`,
which compares the times for the events determined by the edge-triples `l` and `r`.

The time of an event given by an edge triple (which defines three oriented lines) is the Euclidean distance `t`
at which the corresponding offset lines at `t` intersect in a single point.

\pre `l` and `r` must be edge-triples corresponding to events that actually exist (as determined
by the predicate `Do_ss_event_exist_2`).
*/
typedef unspecified_type Compare_ss_event_times_2;

/*!
A predicate object type.

Must provide `bool operator()( const Trisegment_2_ptr& l, const Trisegment_2_ptr& r)`, which determines
if the events given by `l` and `r` are coincident in time and space; that is, both triples
of offset lines intersect at the same point and at the same Euclidean distance from their sources.

\pre `l` and `r` must be edge-triples corresponding to events that actually exist (as determined
by the predicate `Do_ss_event_exist_2`).
*/
typedef unspecified_type Are_ss_events_simultaneous_2;

/*!
A construction object type.

Must provide `boost::optional< boost::tuple<FT, Point_2> > operator()( const Trisegment_2_ptr& e)`,
which returns the Euclidean distance `t >= 0` and the intersection point at which the corresponding
three <I>offset lines at `t`</I> intersect if they do.

If the values cannot be computed, not even approximately (because of overflow for instance),
an empty optional must be returned.
*/
typedef unspecified_type Construct_ss_event_time_and_point_2;

/*!
A construction object type.

Must provide `Trisegment_2_ptr operator()( const Segment_2& e0, const Segment_2& e1, const Segment_2& e2 ) const`,
which returns an intrusive pointer to an object containing the three <I>oriented</I> lines
defined by the three input edges.
*/
typedef unspecified_type Construct_ss_trisegment_2;

/*!
A predicate object type.

Must provide
`bool operator()( const Point_2& contour_node, const Segment_2& opposite_edge ) const`,
and
`bool operator()( const Trisegment_2_ptr& skeleton_node, const Segment_2& opposite_edge ) const`,
which return whether the point (either defined directly as a `Point_2` if it is a contour node,
or indirectly by its edge triple if it is a skeleton node) lies in the halfplane to the left
of the directed edge `opposite_edge`.

\pre `skeleton_node` must be an edge-triple corresponding to an event that actually exists (as determined
by the predicate `Do_ss_event_exist_2`).
*/
typedef unspecified_type Is_edge_facing_ss_node_2;

/*!
A predicate object type.

Must provide `CGAL::Oriented_side operator()( const Trisegment_2_ptr& e, const Segment_2& e0, const Segment_2& e1, const Trisegment_2_ptr& e01_event, bool e0_is_primary) const`,
which returns the oriented side of the event point described by the edge triple `e`
w.r.t the (positive) bisector `[e0,e1]`.

\pre `skeleton_node` must be an edge-triple corresponding to an event that actually exists (as determined
by the predicate `Do_ss_event_exist_2`).
*/
typedef unspecified_type Oriented_side_of_event_point_wrt_bisector_2;

/// @}

/// \name Predicates and Constructors Functors
/// @{

/*! */
Are_ss_events_simultaneous_2 are_ss_events_simultaneous_2_object() const;
/*! */
Compare_ss_event_times_2 compare_ss_event_times_2_object() const;
/*! */
Do_ss_event_exist_2 do_ss_event_exist_2_object() const;
/*! */
Is_edge_facing_ss_node_2 is_edge_facing_ss_node_2_object() const;
/*! */
Oriented_side_of_event_point_wrt_bisector_2 oriented_side_of_event_point_wrt_bisector_2_object() const;

/*! */
Construct_ss_event_time_and_point_2 construct_ss_event_time_and_point_2_object() const;
/*! */
Construct_ss_trisegment_2 construct_ss_trisegment_2_object() const;

/// @}

}; /* end StraightSkeletonBuilderTraits_2 */
