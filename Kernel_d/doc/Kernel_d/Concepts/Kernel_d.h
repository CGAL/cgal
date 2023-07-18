
/*!
\ingroup PkgKernelDKernelConcept
\cgalConcept

The concept of a <I>kernel</I> is defined by a set of requirements on
the provision of certain types and access member functions to create
objects of these types. The types are function object classes to be used
within the algorithms and data structures in the basic library of \cgal.
This allows you to use any model of a kernel as a traits class in
the \cgal algorithms and data structures, unless they require types
beyond those provided by a kernel.

`Kernel_d` subsumes the concept of a <I>\f$ d\f$-dimensional kernel</I>.

A kernel provides types, construction objects, and generalized
predicates. The former replace constructors of the kernel classes and
constructive procedures in the kernel. There are also function objects
replacing operators, especially for equality testing.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Cartesian_d<FieldNumberType>}
\cgalHasModels{CGAL::Homogeneous_d<RingNumberType>}
\cgalHasModels{CGAL::Epick_d<DimensionTag>}
\cgalHasModels{CGAL::Epeck_d<DimensionTag>}
\cgalHasModelsEnd
*/
class Kernel_d {
public:

/// \name Types
/// @{

/*!
a number type that is a model for `FieldNumberType`
*/
typedef unspecified_type FT;

/*!
a number type that is a model for `RingNumberType`
*/
typedef unspecified_type RT;

/*!
the dimension of the ambient space. It must be either `Dimension_tag<d>` for some integer d or `Dynamic_dimension_tag`.
*/
typedef unspecified_type Dimension;

/// @}

/// \name Coordinate Access
/// @{

/*!
a type that allows to iterate over
the %Cartesian coordinates
*/
typedef unspecified_type Cartesian_const_iterator_d;

/// @}

/// \name Geometric Objects
/// @{

/*!

*/
typedef unspecified_type Point_d;

/*!

*/
typedef unspecified_type Vector_d;

/*!

*/
typedef unspecified_type Direction_d;

/*!

*/
typedef unspecified_type Hyperplane_d;

/*!

*/
typedef unspecified_type Line_d;

/*!

*/
typedef unspecified_type Ray_d;

/*!

*/
typedef unspecified_type Segment_d;

/*!

*/
typedef unspecified_type Iso_box_d;

/*!

*/
typedef unspecified_type Sphere_d;

/*!

*/
typedef unspecified_type Aff_transformation_d;

/// @}

/// \name Constructions
/// @{

/*!

*/
typedef unspecified_type Barycentric_coordinates_d;

/*!
a model of `Kernel_d::Center_of_sphere_d`
*/
typedef unspecified_type Center_of_sphere_d;

/*!
a model of `Kernel_d::Compute_coordinate_d`
*/
typedef unspecified_type Compute_coordinate_d;

/*!

*/
typedef unspecified_type Construct_point_d;

/*!

*/
typedef unspecified_type Construct_vector_d;

/*!

*/
typedef unspecified_type Construct_direction_d;

/*!

*/
typedef unspecified_type Construct_hyperplane_d;

/*!

*/
typedef unspecified_type Construct_segment_d;

/*!

*/
typedef unspecified_type Construct_iso_box_d;

/*!

*/
typedef unspecified_type Construct_line_d;

/*!

*/
typedef unspecified_type Construct_ray_d;

/*!

*/
typedef unspecified_type Construct_sphere_d;

/*!

*/
typedef unspecified_type Construct_aff_transformation_d;

/*!
a model of `Kernel_d::ConstructCartesianConstIterator_d`
*/
typedef unspecified_type Construct_cartesian_const_iterator_d;

/*!
a model of `Kernel_d::Intersect_d`
*/
typedef unspecified_type Intersect_d;

/*!
a model of `Kernel_d::Linear_base_d`
*/
typedef unspecified_type Linear_base_d;

/*!
a model of `Kernel_d::Midpoint_d`
*/
typedef unspecified_type Midpoint_d;

/*!
a model of `Kernel_d::Orthogonal_vector_d`
*/
typedef unspecified_type Orthogonal_vector_d;

/*!
a model of `Kernel_d::Point_of_sphere_d`
*/
typedef unspecified_type Point_of_sphere_d;

/*!
a model of `Kernel_d::Point_to_vector_d`
*/
typedef unspecified_type Point_to_vector_d;

/*!
a model of `Kernel_d::Squared_distance_d`
*/
typedef unspecified_type Squared_distance_d;

/*!
a model of `Kernel_d::Value_at_d`
*/
typedef unspecified_type Value_at_d;

/*!
a model of `Kernel_d::Vector_to_point_d`
*/
typedef unspecified_type Vector_to_point_d;

/// @}

/// \name Generalized Predicates
/// @{

/*!
a model of `Kernel_d::Affine_rank_d`
*/
typedef unspecified_type Affine_rank_d;

/*!
a model of `Kernel_d::Affinely_independent_d`
*/
typedef unspecified_type Affinely_independent_d;

/*!
a model of `Kernel_d::Compare_lexicographically_d`
*/
typedef unspecified_type Compare_lexicographically_d;

/*!
a model of `Kernel_d::Component_accessor_d`
*/
typedef unspecified_type Component_accessor_d;

/*!
a model of `Kernel_d::Contained_in_affine_hull_d`
*/
typedef unspecified_type Contained_in_affine_hull_d;

/*!
a model of `Kernel_d::Contained_in_linear_hull_d`
*/
typedef unspecified_type Contained_in_linear_hull_d;

/*!
a model of `Kernel_d::Contained_in_simplex_d`
*/
typedef unspecified_type Contained_in_simplex_d;

/*!
a model of `Kernel_d::Equal_d`
*/
typedef unspecified_type Equal_d;

/*!
a model of `Kernel_d::Has_on_positive_side_d`
*/
typedef unspecified_type Has_on_positive_side_d;

/*!
a model of `Kernel_d::Less_coordinate_d`
*/
typedef unspecified_type Less_coordinate_d;

/*!
a model of `Kernel_d::Less_lexicographically_d`
*/
typedef unspecified_type Less_lexicographically_d;

/*!
a model of `Kernel_d::Less_or_equal_lexicographically_d`
*/
typedef unspecified_type Less_or_equal_lexicographically_d;

/*!
a model of `Kernel_d::Linear_rank_d`
*/
typedef unspecified_type Linear_rank_d;

/*!
a model of `Kernel_d::Linearly_independent_d`
*/
typedef unspecified_type Linearly_independent_d;

/*!
a model of `Kernel_d::Orientation_d`
*/
typedef unspecified_type Orientation_d;

/*!
a model of `Kernel_d::Oriented_side_d`
*/
typedef unspecified_type Oriented_side_d;

/*!
a model of `Kernel_d::Point_dimension_d`
*/
typedef unspecified_type Point_dimension_d;

/*!

*/
typedef unspecified_type Position_on_line_d;

/*!
a model of `Kernel_d::Side_of_bounded_sphere_d`
*/
typedef unspecified_type Side_of_bounded_sphere_d;

/*!
a model of `Kernel_d::Side_of_oriented_sphere_d`
*/
typedef unspecified_type Side_of_oriented_sphere_d;

/// @}

/// \name Operations
/// The following member functions return function objects of the
/// types listed above. The name of the access function is the name of
/// the type returned with an `_object` suffix and no capital letter
/// at the beginning. We only give two examples to show the
/// scheme. For the functors `Construct_point_d` and `Orientation_d`
/// the corresponding functions are:
/// @{

/*!

*/
Kernel_d::Construct_point_d construct_point_d_object() const;

/*!

*/
Kernel_d::Orientation_d orientation_d_object() const;

/// @}

}; /* end Kernel_d */

