/*!
\ingroup PkgStraightSkeleton2Concepts
\cgalConcept

\cgalRefines{StraightSkeletonBuilderTraits_2}

The concept `PolygonOffsetBuilderTraits_2` describes the requirements for the geometric traits class
required by the algorithm class `CGAL::Polygon_offset_builder_2`.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Polygon_offset_builder_traits_2}
\cgalHasModelsEnd

\sa `CGAL::Polygon_offset_builder_2<Ss,Traits,Container>`
*/
class PolygonOffsetBuilderTraits_2 {
public:

/// \name Types
/// @{

/*!
A predicate object type.

Must provide

`Comparison_result operator()( FT d, const Trisegment_2_ptr& et) const`,

which compares the Euclidean distance `d` with the event time for `et`. Such event time
is the Euclidean distance at which the <I>offset lines</I> intersect in a single point.
The source of such offset lines is given by the three <I>oriented</I> lines defined by the edge-triple `et`.

\pre `et` must be an edge-triple corresponding to an event that actually exist (that is,
there must exist an offset distance `t > 0` at which the offset lines do intersect at a single point).
*/
typedef unspecified_type Compare_offset_against_event_time_2;

/*!
A construction object type.

Must provide

`std::optional<Point_2> operator()(const FT& t, const Segment_2& e0, const Segment_2& e1, const Trisegment_2_ptr& et) const`,

which constructs the point of intersection of the lines obtained by offsetting
the oriented lines given by `e0` and `e0` a Euclidean distance `t`.

If `e0` and `e1` are collinear, if `et` is not specified (`nullptr`), then the midpoint should be returned,
otherwise, the event point of `et` should be returned.

If the point cannot be computed, not even approximately (because of overflow for instance),
an empty optional must be returned.

\pre `et` must be an edge-triple corresponding to an event that actually exist (that is,
there must exist an offset distance `t > 0` at which the offset lines do intersect at a single point).
*/
typedef unspecified_type Construct_offset_point_2;

/*! */
Compare_offset_against_event_time_2 compare_offset_against_event_time_2_object() const;
/*! */
Construct_offset_point_2 construct_offset_point_2_object() const;

/// @}

}; /* end PolygonOffsetBuilderTraits_2 */
