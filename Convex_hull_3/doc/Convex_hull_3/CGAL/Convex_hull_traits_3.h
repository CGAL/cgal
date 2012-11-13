namespace CGAL {

/*!
\ingroup PkgConvexHull3Traits

The class `Convex_hull_traits_3` serves as a traits class for the function 
`convex_hull_3()`. This is the default traits class for this 
function when `R` is a kernel with exact predicates but inexact constructions 
(note that the type `Plane_3` is a triple of `Point_3` and not `R::Plane_3`). 

\cgalModels ::ConvexHullTraits_3 
\cgalModels ::IsStronglyConvexTraits_3 

*/
template< typename R >
class Convex_hull_traits_3 {
public:

/// \name Types 
/// @{

/*! 

*/ 
typedef R::Point_3 Point_3; 

/*! 

*/ 
typedef R::Segment_3 Segment_3; 

/*! 

*/ 
typedef R::Triangle_3 Triangle_3; 

/*! 
A triple of points, in order to avoid the need for exact constructions. 
*/ 
typedef Hidden_type Plane_3;; 

/*! 

*/ 
typedef R::Vector_3 Vector_3; 

/*! 

*/ 
typedef Polyhedron_default_traits_3<R> Poly_traits; 

/*! 

*/ 
typedef Halfedge_data_structure_polyhedron_default_3<R> HDS; 

/*! 

*/ 
typedef Polyhedron_3<Poly_traits, HDS> Polyhedron_3; 

/*! 

*/ 
typedef R::Construct_segment_3 Construct_segment_3; 

/*! 

*/ 
typedef R::Construct_ray_3 Construct_ray_3; 

/*! 

*/ 
typedef R::Construct_plane_3 Construct_plane_3; 

/*! 

*/ 
typedef R::Construct_vector_3 Construct_vector_3; 

/*! 

*/ 
typedef R::Construct_triangle_3 Construct_triangle_3; 

/*! 

*/ 
typedef R::Construct_centroid_3 Construct_centroid_3; 

/*! 

*/ 
typedef R::Construct_orthogonal_vector_3 
Construct_orthogonal_vector_3; 

/*! 

*/ 
R::Equal_3 Equal_3; 

/*! 

*/ 
R::Collinear_3 Collinear_3; 

/*! 

*/ 
R::Coplanar_3 Coplanar_3; 

/*! 

*/ 
R::Less_distance_to_point_3 Less_distance_to_point_3; 

/*! 

*/ 
R::Has_on_positive_side_3 Has_on_positive_side_3; 

/*! 

*/ 
R::Less_signed_dist_to_plane_3 
Less_signed_distance_to_plane_3; 

/*! 

*/ 
Projection_traits_xy_3<R> Traits_xy; 

/*! 

*/ 
Projection_traits_xz_3<R> Traits_xz; 

/*! 

*/ 
Projection_traits_yz_3<R> Traits_yz; 

/*! 

*/ 
R::Ray_3 Ray_3; 

/*! 

*/ 
R::Has_on_3 Has_on_3; 

/*! 

*/ 
R::Oriented_side_3 Oriented_side_3; 

/*! 

*/ 
R::Do_intersect_3 Do_intersect_3; 

/// @} 

/// \name Creation 
/// @{

/*! 
copy constructor. 
*/ 
Convex_hull_traits_3(Convex_hull_traits_3& t); 

/// @} 

/// \name Operations 
/// @{

/*! 

*/ 
Construct_segment_3 
construct_segment_3_object() const; 

/*! 

*/ 
Construct_ray_3 
construct_ray_3_object() const; 

/*! 

*/ 
Construct_plane_3 
construct_plane_3_object() const; 

/*! 

*/ 
Construct_triangle_3 
construct_triangle_3_object() const; 

/*! 

*/ 
Construct_vector_3 
construct_vector_3_object() const; 

/*! 

*/ 
Construct_centroid_3 
construct_centroid_3_object() const; 

/*! 

*/ 
Construct_orthogonal_vector_3 
construct_orthogonal_vector_3_object() const; 

/*! 

*/ 
Equal_3 
equal_3_object() const; 

/*! 

*/ 
Collinear_3 
collinear_3_object() const; 

/*! 

*/ 
Coplanar_3 
coplanar_3_object() const; 

/*! 

*/ 
Has_on_3 
has_on_3_object() const; 

/*! 

*/ 
Less_distance_to_point_3 
less_distance_to_point_3_object() const; 

/*! 

*/ 
Has_on_positive_side_3 
has_on_positive_side_3_object() const; 

/*! 

*/ 
Oriented_side_3 
oriented_side_3_object() const; 

/*! 

*/ 
Do_intersect_3 
do_intersect_3_object() const; 

/*! 

*/ 
Less_signed_distance_to_plane_3 
less_signed_distance_to_plane_3_object() const; 

/// @}

}; /* end Convex_hull_traits_3 */
} /* end namespace CGAL */
