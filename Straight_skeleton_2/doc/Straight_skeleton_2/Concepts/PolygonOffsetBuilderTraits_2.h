/*!
\ingroup PkgStraightSkeleton2Concepts
\cgalConcept

The concept `PolygonOffsetBuilderTraits_2` describes the requirements for the geometric traits class required by the algorithm class `Polygon_offset_builder_2<Ss,Gt,Polygon_2>`. 

\cgalHasModel CGAL::Polygon_offset_builder_traits_2

\sa `CGAL::Polygon_offset_builder_2<Ss,Gt,Polygon_2>` 
\sa `CGAL::Polygon_offset_builder_traits_2<K>` 

*/

class PolygonOffsetBuilderTraits_2 {
public:

/// \name Types 
/// @{

/*!
A model of the `Kernel` concept. 
*/ 
typedef unspecified_type Kernel; 

/*!
A model of the `FieldWithSqrt` concept provided by the kernel. This type is used to represent the 
coordinates of the input points and to specify the desired offset distance. 
*/ 
typedef unspecified_type FT; 

/*!
A 2D point type 
*/ 
typedef unspecified_type Point_2; 

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
A predicate object type. 

Must provide `Comparison_result operator()( FT d, EdgeTriple const& et) const`, which compares the Euclidean distance `d` with the event time for `et`. Such event time is the Euclidean distance at which the <I>offset lines</I> intersect in a single point. The source of such offset lines is given by the 3 <I>oriented</I> lines defined by the edge-triple `et` 
\pre `et` must be an edge-triple corresponding to an event that actually exist (that is, there must exist an offset distance `t > 0` at which the offset lines do intersect at a single point. 
*/ 
typedef unspecified_type Compare_offset_against_event_time_2; 

/*!
A construction object type. 

Must provide `boost::optional<Point_2> operator()( FT t, Edge const& x, Edge const& y) const`, which constructs the point of intersection of the lines obtained by offsetting the oriented lines given by `x` and `y` an Euclidean distance `t`. 
If the point cannot be computed, not even approximately (because of overflow for instance), an empty optional must be returned. 

\pre `x` and `y` must intersect in a single point 
*/ 
typedef unspecified_type Construct_offset_point_2; 

/*!
A construction object type. 

Must provide `Vertex operator()( Point_2 const& p)`, which given a `Point_2` `p` returns a Vertex encapsulating the corresponding (x,y) pair of Cartesian coordinates. 
*/ 
typedef unspecified_type Construct_ss_vertex_2; 

/*!
A construction object type. 

Must provide `Edge operator()( Point_2 const& s, Point_2 const& t)`, which given source and target points `s` and `t` returns an Edge encapsulating the corresponding input segment (in Cartesian coordinates.) 
*/ 
typedef unspecified_type Construct_ss_edge_2; 

/*!
A construction object type. 

Must provide `Triedge operator()( Edge const& e0, Edge const& e1, Edge const& e2)`, which given the 3 edges that define an event, `e0`, `e1` and `e2`, returns a Triedge encapsulating them. 
*/ 
typedef unspecified_type Construct_ss_triedge_2; 

/// @}

}; /* end PolygonOffsetBuilderTraits_2 */
