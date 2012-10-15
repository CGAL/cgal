
/*!
\ingroup PkgTriangulation2Concepts
\cgalconcept

The concept `RegularTriangulationTraits_2` describe the requirements 
for the traits class of regular triangulations. It refines the 
concept `TriangulationTraits_2` requiring the type 
`CGAL::Weighted_point` and the *power test* predicate on those 
weighted points. 
A weighted point is basically 
a point augmented with a scalar weight. It can be seen as a circle 
when the weight is interpreted as a square radius. 
The power test on weighted points 
is the fundamental test to build regular triangulations 
as the `side_of_oriented_circle` test is the fundamental test 
of Delaunay triangulations. 

\refines ::TriangulationTraits_2 

\hasModel `CGAL::Regular_triangulation_euclidean_traits_2<K,Weight>` 
\hasModel `CGAL::Regular_triangulation_filtered_traits_2<FK>` 

\sa `TriangulationTraits_2` 
\sa `CGAL::Regular_triangulation_2<Traits,Tds>` 

*/

class RegularTriangulationTraits_2 {
public:

/// \name Types 
/// @{

/*! 
Another name for the point type. 
*/ 
typedef Hidden_type Bare_point; 

/*! 
The weighted point type, it has to be 
a model of the concept `WeightedPoint`. 
*/ 
typedef Hidden_type Weighted_point_2; 

/*! 
A predicate object type. Must provide 
the operators: 

- `Oriented_side operator() ( Weighted_point_2 p, Weighted_point_2 q, Weighted_point_2 r, Weighted_point_2 s) ` 
which is the power test for points `p`, `q`, `r` and 
`s`. \pre the bare points corresponding to `p`, `q`, `r` are not collinear. 

- `Oriented_side operator() ( Weighted_point_2 p, Weighted_point_2 q, Weighted_point_2 r) ` 
which is the degenerated power test for collinear points 
`p`, `q`, `r`. 
\pre the bare points corresponding to `p`, `q`, `r` are collinear and `p != q`. 

- `Oriented_side operator() ( Weighted_point_2 p, Weighted_point_2 q) ` 
which is the degenerated power test for weighted points 
`p` and `q` whose corresponding bare-points are identical. 
\pre the bare points corresponding to `p` and `q` are identical. 

*/ 
typedef Hidden_type Power_test_2; 

/// @}

/// \name Optional Types
/// The following type/predicate is required if a call to `nearest_power_vertex` is issued:
/// @{

/*! 
A predicate object type. Must 
provide the operator: 

`Comparison_result operator()(Bare_point p, Weighted_point_2 q, Weighted_point_2 r)`, 
which compares the power distance between `p` and `q` to the 
power distance between `p` and `r`. 
*/ 
typedef Hidden_type Compare_power_distance_2; 

/*! 
A constructor 
object which constructs the weighted circumcenter of three 
weighted points. Provides the operator 

`Bare_point operator() ( Weighted_point_2 p, Weighted_point_2 q, Weighted_point_2 r);` 
*/ 
typedef Hidden_type Construct_weighted_circumcenter_2; 

/*! 
A constructor type which 
constructs 
the radical axis of two weighted points. Provides the operator : 

`Line_2 operator() ( Weighted_point_2 p, Weighted_point_2 q);` 
*/ 
typedef Hidden_type Construct_radical_axis_2; 

/// @} 

/// \name Creation 
/// @{

/*! 
default constructor. 
*/ 
RegularTriangulationTraits_2(); 

/*! 
copy constructor. 
*/ 
RegularTriangulationTraits_2 ( 
const RegularTriangulatioTraits_2& ); 

/*! 
assignment operator 
*/ 
RegularTriangulationTraits_2& operator= 
(const RegularTriangulationTraits_2& ); 

/// @} 

/// \name Access to predicate and constructors objects 
/// @{

/*! 

*/ 
Power_test_2 power_test_2_object(); 

/*! 

*/ 
Compare_power_distance_2 compare_power_distance_2_object(); 

/*! 

*/ 
Construct_weighted_circumcenter_2 
construct_weighted_circumcenter_2_object(); 

/*! 

*/ 
Construct_radical_axis_2 
construct_radical_axis_2_object(); 

/// @}

}; /* end RegularTriangulationTraits_2 */

