
/*!
\ingroup PkgBooleanSetOperations2Concepts
\cgalConcept

This concept defines the minimal set of geometric predicates needed to 
perform the Boolean-set operations. It refines the directional \f$ x\f$-monotone 
arrangement-traits concept. In addition to the `Point_2` and 
`X_monotone_curve_2` types defined in the generalized concept, it defines 
a type that represents a general polygon and another one that represents 
general polygon with holes. It also requires operations that operate on these 
types. 

\cgalRefines `ArrangementDirectionalXMonotoneTraits_2` 

\cgalHasModel `CGAL::Gps_segment_traits_2<Kernel,Container,ArrSegmentTraits>`
\cgalHasModel `CGAL::Gps_circle_segment_traits_2<Kernel>`
\cgalHasModel `CGAL::Gps_traits_2<ArrTraits,GeneralPolygon>` 

\sa `ArrangementDirectionalXMonotoneTraits_2` 

*/

class GeneralPolygonSetTraits_2 {
public:

/// \name Types 
/// @{

/*!
represents a general polygon. Must be a model of the `GpsTraitsGeneralPolygon_2` concept. 
*/ 
typedef unspecified_type Polygon_2; 

/*!
represents a general polygon with holes. Must be a model of the `GpsTraitsGeneralPolygonWithHoles_2` concept. 
*/ 
typedef unspecified_type Polygon_with_holes_2; 

/*!
A const iterator of curves. Its value type is const 
`X_monotone_curve_2`. 
*/ 
typedef unspecified_type Curve_const_iterator; 

/// @} 

/// \name Functor Types 
/// @{

/*!
a functor that constructs a general polygon from a range of 
\f$ x\f$-monotone curves. It uses the operator 

`void operator() (InputIterator begin, Input iterator end, Polygon_2 & pgn)`, 

parameterized by the `InputIterator` type. 
*/ 
typedef unspecified_type Construct_polygon_2; 

/*!
a functor that constructs a general polygon with holes from a general polygon and, optionally, a range of holes. It uses the operator 

`void operator() (const General_polygon_2& pgn_boundary)` or 

`void operator() (const General_polygon_2& pgn_boundary, HolesInputIterator h_begin, HolesInputIterator h_end)` parameterized by the `HolesInputIterator` type. 
*/ 
typedef unspecified_type Construct_general_polygon_with_holes_2; 

/*!
A functor for constructing the outer boundary of a polygon with holes. It uses the operator 

`General_polygon_2 operator()(const General_polygon_with_holes_2& pol_wh)`. 
*/ 
typedef unspecified_type Construct_outer_boundary; 

/*!
A functor for constructing the container of holes of a polygon with holes. It returns the begin/end iterators for the holes It uses the operator 

`std::pair<Hole_const_iterator, Hole_const_iterator> operator()(const General_polygon_with_holes_2& pol_wh)`. 
*/ 
typedef unspecified_type Construct_holes; 

/*!
A functor for checking if polygon with holes has no outer boundary. It uses the operator 

`bool operator()(const General_polygon_with_holes_2& pol_wh)`. 
*/ 
typedef unspecified_type Is_unbounded; 

/// @} 

/// \name Creation 
/// @{

/*!
default constructor. 
*/ 
GeneralPolygonSetTraits_2(); 

/*!
copy constructor 
*/ 
GeneralPolygonSetTraits_2(GeneralPolygonSetTraits_2 other); 

/*!
assignment operator. 
*/ 
GeneralPolygonSetTraits_2 operator=(other); 

/// @} 

/// \name Accessing Functor Objects 
/// @{

/*!
returns a functor that constructs a polygon. 
*/ 
Construct_polygon_2 construct_polygon_2_object(); 

/*!
returns a functor that obtains the curves of a polygon. 
*/ 
Construct_curves_2 construct_curves_2_object(); 

/*!
returns a functor that constructs a polygon with holes. 
*/ 
Construct_general_polygon_with_holes_2 construct_polygon_with_holes_2_object() const; 

/*!
returns a functor that obtains the outer boundary of a polygon with holes. 
*/ 
Construct_outer_boundary construct_outer_boundary_object() const; 

/*!
returns a functor that obtains begin/end iterators over a container of holes. 
*/ 
Construct_holes construct_holes_object() const; 

/*!
returns a functor that determines if the polygon with holes is unbounded 
*/ 
Is_unbounded construct_is_unbounded_object(); 

/// @}

}; /* end GeneralPolygonSetTraits_2 */

