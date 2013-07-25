
/*!
\ingroup PkgTriangulation2Concepts
\cgalConcept

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

\cgalRefines `TriangulationTraits_2` 

\cgalHasModel `CGAL::Regular_triangulation_euclidean_traits_2<K,Weight>` 
\cgalHasModel `CGAL::Regular_triangulation_filtered_traits_2<FK>` 

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
typedef unspecified_type Bare_point; 

/*!
The weighted point type. It has to be 
a model of the concept `WeightedPoint`. 
*/ 
typedef unspecified_type Weighted_point_2; 

/*!
A function object which must provide operators for the power test applied to two, three, and four points. 

Must provide the operators: 

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
typedef unspecified_type Power_test_2; 

/// @}

/// \name Optional Types
/// The following type/predicate is required for a call to `nearest_power_vertex`:
/// @{

/*!
A function object for computing two power distances. 
Must provide the operator: 

`Comparison_result operator()(Bare_point p, Weighted_point_2 q, Weighted_point_2 r)`, 
which compares the power distance between `p` and `q` to the 
power distance between `p` and `r`. 
*/ 
typedef unspecified_type Compare_power_distance_2; 

/*!
A function
object which constructs the weighted circumcenter of three 
weighted points. 

Must provide the operator:

`Bare_point operator() ( Weighted_point_2 p, Weighted_point_2 q, Weighted_point_2 r);` 
*/ 
typedef unspecified_type Construct_weighted_circumcenter_2; 

/*!
A function object constructs 
the radical axis of two weighted points. 

Must provide the operator: 

`Line_2 operator() ( Weighted_point_2 p, Weighted_point_2 q);` 
*/ 
typedef unspecified_type Construct_radical_axis_2; 

/// @} 

/// \name Creation 
/// @{

/*!
%Default constructor. 
*/ 
RegularTriangulationTraits_2(); 

/*!
Copy constructor. 
*/ 
RegularTriangulationTraits_2 ( 
const RegularTriangulatioTraits_2& ); 

/*!
Assignment operator.
*/ 
RegularTriangulationTraits_2& operator= 
(const RegularTriangulationTraits_2& ); 

/// @} 

/// \name Access to Predicate and Constructors Objects 
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

