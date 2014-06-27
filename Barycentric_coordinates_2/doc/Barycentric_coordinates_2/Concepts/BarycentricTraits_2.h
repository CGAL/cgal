/*!
\ingroup PkgBarycentric_coordinates_2Concepts
\cgalConcept

Requirements of the template parameter `Traits` for all the classes with two-dimensional barycentric coordinates from the namespace `CGAL::Barycentric_coordinates`. 

\cgalHasModel `CGAL::Barycentric_coordinates::Barycentric_traits_2`
\cgalHasModel All models of `Kernel`

*/

class BarycentricTraits_2 {

public:

/// \name Types 
/// @{

/*!
	a model of `FieldNumberType`
*/ 
typedef unspecified_type FT;

/// @}

/// \name Two-dimensional Geometric Objects
/// @{

/*!
	a model of `Kernel::Point_2`
*/ 
typedef unspecified_type Point_2; 

/*!
	a model of `Kernel::Vector_2`
*/ 
typedef unspecified_type Vector_2; 

/// @} 

/// \name Two-dimensional Constructions
/// @{

/*!
	a model of `Kernel::ComputeArea_2`
*/ 
typedef unspecified_type Compute_area_2; 

/*!
    a model of `Kernel::ComputeSquaredDistance_2` 
*/ 
typedef unspecified_type Compute_squared_distance_2; 

/*!
    a model of `Kernel::ComputeSquaredLength_2` 
*/ 
typedef unspecified_type Compute_squared_length_2;

/*!
    a model of `Kernel::ComputeScalarProduct_2`
*/
typedef unspecified_type Compute_scalar_product_2;

/*!
    a model of `RealEmbeddableTraits_::ToDouble` and `AlgebraicStructureTraits_::Sqrt`
*/
typedef unspecified_type Sqrt;

/// @}

/// \name Two-dimensional Generalized Predicates
/// @{

/*!
	a model of `Kernel::Equal_2`
*/ 
typedef unspecified_type Equal_2; 

/*!
	a model of `Kernel::Collinear_2` 
*/ 
typedef unspecified_type Collinear_2; 

/*!
	a model of `Kernel::CollinearAreOrderedAlongLine_2` 
*/ 
typedef unspecified_type Collinear_are_ordered_along_line_2; 

/// @} 

}; /* end BarycentricTraits_2 */