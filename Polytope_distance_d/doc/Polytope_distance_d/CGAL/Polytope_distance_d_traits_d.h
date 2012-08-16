
namespace CGAL {

/*!
\ingroup PkgOptimalDistances

The class `Polytope_distance_d_traits_d` is a traits class for the \f$ d\f$-dimensional 
optimisation algorithms using the \f$ d\f$-dimensional \cgal kernel. 

Requirements 
-------------- 

The template parameter `K` is a model for `Kernel`. Template 
parameters `ET` and `NT` are models for `RingNumberType`. 

The second and third template parameter have default type `K::RT`. 

\models ::PolytopeDistanceDTraits 

\sa `CGAL::Polytope_distance_d<Traits>` 
\sa `CGAL::Polytope_distance_d_traits_2<K,ET,NT>` 
\sa `CGAL::Polytope_distance_d_traits_3<K,ET,NT>` 
\sa `PolytopeDistanceDTraits` 

*/
template< typename K, typename ET, typename NT >
class Polytope_distance_d_traits_d {
public:

/// \name Types 
/// @{

/*! 
typedef to `K::Point_d`. 
*/ 
typedef Hidden_type Point_d; 

/*! 
typedef to `K::Rep_tag`. 
*/ 
typedef Hidden_type Rep_tag; 

/*! 
typedef to `K::RT`. 
*/ 
typedef Hidden_type RT; 

/*! 
typedef to `K::FT`. 
*/ 
typedef Hidden_type FT; 

/*! 
typedef to `K::Access_dimension_d`. 
*/ 
typedef Hidden_type Access_dimension_d; 

/*! 
typedef to `K::Access_coordinates_begin_d`. 
*/ 
typedef Hidden_type Access_coordinates_begin_d; 

/*! 
typedef to `K::Construct_point_d`. 
*/ 
typedef Hidden_type Construct_point_d; 

/*! 
second template parameter (default is `K::RT`). 
*/ 
typedef Hidden_type ET; 

/*! 
third template parameter (default is `K::RT`). 
*/ 
typedef Hidden_type NT; 

/// @} 

/// \name Creation 
/// @{

/*! 
default constructor. 
*/ 
Polytope_distance_d_traits_d( ); 

/*! 
copy constructor. 
*/ 
Polytope_distance_d_traits_d( 
const Polytope_distance_d_traits_d<K,ET,NT>&); 

/// @} 

/// \name Operations 
/// The following functions just return the corresponding function class object.
/// @{

/*! 

*/ 
Access_dimension_d 
access_dimension_d_object() const; 

/*! 

*/ 
Access_coordinates_begin_d 
access_coordinates_begin_d_object() const; 

/*! 

*/ 
Construct_point_d 
construct_point_d_object() const; 

/// @}

}; /* end Polytope_distance_d_traits_d */
} /* end namespace CGAL */
