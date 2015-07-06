
namespace CGAL {

/*!
\ingroup PkgTriangulation3TraitsClasses

The class `Regular_triangulation_euclidean_traits_3` is designed as a default traits class for the 
class `Regular_triangulation_3<RegularTriangulationTraits_3,TriangulationDataStructure_3>`. 
It provides `Weighted_point_3`, a class for weighted points, which derives 
from the three dimensional point class `K::Point_3`. 

\tparam K must be a model of the `Kernel` concept. 

\tparam Weight is optional. If is it not provided, `K::RT` will be used. 

The class is a model of the concept `RegularTriangulationTraits_3` 
but it also contains predicates and constructors on weighted points 
that are not required in the 
concept `RegularTriangulationTraits_3`. 

Note that filtered predicates are automatically used if the 
Boolean `Has_filtered_predicates` in the kernel provided as template parameter 
of that class is set to `true`. 

\cgalModels `RegularTriangulationTraits_3`

\cgalHeading{Operations}

The following functions give access to the predicate and constructor 
functors. 

*/
template< typename K, typename Weight >
class Regular_triangulation_euclidean_traits_3 : public K {
public:

/// \name Types 
/// @{

/*!
The type for point \f$ p\f$ of a weighted point \f$ {p}^{(w)}=(p,w_p)\f$. 
*/ 
typedef K::Point_3 Bare_point; 

/*!
The type for weighted points. 
*/ 
typedef Weighted_point <Bare_point, Weight> 
Weighted_point_3; 

/*!
The type for points. 
*/ 
typedef Weighted_point <Bare_point, Weight> 
Point_3; 

/// @} 

/// \name Types for Predicate Functors 
/// @{

/*!
A predicate type for power test. Belongs to the 
`RegularTriangulationTraits_3` concept. 
*/ 
typedef unspecified_type Power_test_3; 

/*!
A predicate type to compare power distance. Belongs to the 
`RegularTriangulationTraits_3` concept. 
*/ 
typedef unspecified_type Compare_power_distance_3; 

/*!
A predicate type. The operator() takes weighted point(s) as arguments, 
together with one weight. It compares the weight of the smallest sphere 
orthogonal to the weighted points with the input weight. 

`Comparison_result operator()( Weighted_point_3 p, Weighted_point_3 q, Weighted_point_3 r, Weighted_point_3 s, FT w) ;` 

`Comparison_result operator()( Weighted_point_3 p, Weighted_point_3 q, Weighted_point_3 r, FT w) ;` 

`Comparison_result operator()( Weighted_point_3 p, Weighted_point_3 q, FT w) ;` 

`Comparison_result operator()( Weighted_point_3 p, FT w) ;` 

*/ 
typedef unspecified_type Compare_weighted_squared_radius_3; 

/*!
A predicate type. The operator() takes weighted points as arguments 
and returns the sign of the power distance of the last one 
with respect to the smallest sphere orthogonal to the others. 

`Sign operator()( Weighted_point_3 p, 		 Weighted_point_3 q, 		 Weighted_point_3 r, 		 Weighted_point_3 s, 		 Weighted_point_3 t) ;` 

`Sign operator()( Weighted_point_3 p, 		 Weighted_point_3 q, 		 Weighted_point_3 r, 		 Weighted_point_3 s) ;` 

`Sign operator()( Weighted_point_3 p, 		 Weighted_point_3 q, 		 Weighted_point_3 r) ;` 

`Sign operator()( Weighted_point_3 p, 		 Weighted_point_3 q) ;` 

*/ 
typedef unspecified_type In_smallest_orthogonal_sphere_3; 

/*!
A predicate type. The operator() is similar to the operator() of 
`In_smallest_orthogonal_sphere_3` 
except that the returned type is not a `Sign` 
but belongs to the enum `Bounded_side` 
(A `NEGATIVE`, `ZERO` and `POSITIVE`) 
corresponding respectively to 
`ON_BOUNDED_SIDE`, `ON_BOUNDARY` and `ON_UNBOUNDED_SIDE`)). 

`Bounded_side operator() ( Weighted_point_3 p, 			 Weighted_point_3 q, 			 Weighted_point_3 r, 			 Weighted_point_3 s, 			 Weighted_point_3 t) ; ` 

`Bounded_side operator() ( Weighted_point_3 p, 			 Weighted_point_3 q, 			 Weighted_point_3 r, 			 Weighted_point_3 s) ; ` 

`Bounded_side operator() ( Weighted_point_3 p, 			 Weighted_point_3 q, 			 Weighted_point_3 r) ; ` 

*/ 
typedef unspecified_type Side_of_bounded_orthogonal_sphere_3; 

/*!
A predicate type. The operator() 
takes weighted points as arguments, 
considers the subspace of points with equal power distance 
with respect to its arguments 
and the intersection of this subspace with the affine hull of the bare 
points associated to the arguments. 
The operator() returns `ON_BOUNDED_SIDE`, 
`ON_BOUNDARY` or `ON_UNBOUNDED_SIDE` 
according to the position of this intersection 
with respect to the simplex formed by the bare 
points. 
This predicate is useful for flow computations. 

`Bounded_side operator()( Weighted_point_3 p, 		 Weighted_point_3 q, 		 Weighted_point_3 r, 		 Weighted_point_3 s) ; ` 

`Bounded_side operator()( Weighted_point_3 p, 		 Weighted_point_3 q, 		 Weighted_point_3 r) ; ` 

`Bounded_side operator()( Weighted_point_3 p, 		 Weighted_point_3 q) ; ` 

*/ 
typedef unspecified_type Does_simplex_intersect_dual_support_3; 

/// @} 

/// \name Types for Constructor Functors 
/// @{

/*!
A constructor type. The operator() constructs the bare point 
which is the center of the smallest orthogonal sphere to the input 
weighted points. 

`Bare_point operator() ( Weighted_point_3 p, 			 Weighted_point_3 q, 			 Weighted_point_3 r, 			 Weighted_point_3 s);` 

`Bare_point operator() ( Weighted_point_3 p, 			 Weighted_point_3 q, 			 Weighted_point_3 r);` 

`Bare_point operator() ( Weighted_point_3 p, 			 Weighted_point_3 q);` 

*/ 
typedef unspecified_type Construct_weighted_circumcenter_3; 

/*!
A functor type. The operator() computes the power distance between its 
arguments. 

`FT operator() ( Weighted_point_3 p, 		 Weighted_point_3 q) ; ` 

*/ 
typedef unspecified_type Compute_power_product_3; 

/*!
A functor type. The operator() computes the squared radius of the 
smallest sphere orthogonal to the argument(s). 

`FT operator() ( Weighted_point_3 p, 		 Weighted_point_3 q, 		 Weighted_point_3 r, 		 Weighted_point_3 s);` 

`FT operator() ( Weighted_point_3 p, 		 Weighted_point_3 q, 		 Weighted_point_3 r);` 

`FT operator() ( Weighted_point_3 p, 		 Weighted_point_3 q);` 

`FT operator() ( Weighted_point_3 p);` 

*/ 
typedef unspecified_type Compute_squared_radius_smallest_orthogonal_sphere_3; 

/*!
A functor type. The operator() takes weighted points as arguments 
and computes the squared radius 
of the sphere centered in the last point and orthogonal 
to the other weighted points. The last argument is a weighted point 
but its weight does not matter. 
This construction is ad hoc for pumping slivers. 
For robustness issue, a predicate to compare critical squared radii 
for a given last point should be needed. 

`FT operator() ( Weighted_point_3 p, 			 Weighted_point_3 q, 			 Weighted_point_3 r, 			 Weighted_point_3 s, 			 Weighted_point_3 t);` 

*/ 
typedef unspecified_type Compute_critical_squared_radius_3; 

/// @} 

/// \name Operations 
/// @{

/*!

*/ 
Power_test_3 power_test_3_object(); 

/*!

*/ 
Compare_power_distance_3 
compare_power_distance_3_object(); 

/*!

*/ 
Compare_weighted_squared_radius_3 
compare_weighted_squared_radius_3_object(); 

/*!

*/ 
In_smallest_orthogonal_sphere_3 
in_smallest_orthogonal_sphere_3_object(); 

/*!

*/ 
Side_of_bounded_orthogonal_sphere_3 
side_of_bounded_orthogonal_sphere_3_object(); 

/*!

*/ 
Does_simplex_intersect_dual_support_3 
does_simplex_intersect_dual_support_3_object(); 

/*!

*/ 
Construct_weighted_circumcenter_3 
construct_weighted_circumcenter_3_object(); 

/*!

*/ 
Compute_power_product_3 
compute_power_product_3_object() ; 

/*!

*/ 
Compute_squared_radius_smallest_orthogonal_sphere_3 
compute_squared_radius_smallest_orthogonal_sphere_3_object() ; 

/*!

*/ 
Compute_critical_squared_radius_3 
compute_critical_squared_radius_3_object(); 

/// @}

}; /* end Regular_triangulation_euclidean_traits_3 */
} /* end namespace CGAL */
