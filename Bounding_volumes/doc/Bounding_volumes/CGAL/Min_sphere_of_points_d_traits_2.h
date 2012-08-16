namespace CGAL {

/*!
\ingroup PkgBoundingVolumes

The class 
`Min_sphere_of_points_d_traits_2<K,FT,UseSqrt,Algorithm>` is a 
model for concept `MinSphereOfSpheresTraits`. It uses the \cgal 
type `Point_2` to represent circles. 

\models ::MinSphereOfSpheresTraits 

\name Parameters 
The last two template parameters, `UseSqrt` and `Algorithm`, have
default arguments, namely `CGAL::Tag_false` and
`CGAL::Default_algorithm`, respectively. The template parameters
of class `Min_sphere_of_points_d_traits_2<K,FT,UseSqrt,Algorithm>`
must fulfill the following requirements:

\tparam K is a model for `Kernel`.
\tparam FT is a number type, which fulfills the requirements of 
type `FT` of concept `MinSphereOfSpheresTraits`: It must be 
either `double` or `float`, or an exact number type. 
\tparam UseSqrt fulfills the 
requirements of type `Use_square_roots` of concept 
`MinSphereOfSpheresTraits`: It must be either `Tag_true` or `Tag_false`. 
\tparam Algorithm fulfills the requirements of type `Algorithm` of
concept `MinSphereOfSpheresTraits`: It must be either
`Default_algorithm`, `LP_algorithm` or `Farthest_first_heuristic`.

*/
template< typename K, typename FT, typename UseSqrt, typename Algorithm >
class Min_sphere_of_points_d_traits_2 {
public:

/// \name Constants 
/// @{

/*! 
is the constant 2, i.e. the dimension of \f$ \R^2\f$. 
*/ 
typedef Hidden_type D; 

/// @} 

/// \name Types 
/// In addition to the types required by the concept
/// `MinSphereOfSpheresTraits`, this model also defines the types
/// `Radius` and `Point`. Here's the complete list of defined types:
/// @{

/*! 

*/ 
typedef Hidden_type FT; 

/*! 

*/ 
typedef Hidden_type Use_square_roots; 

/*! 

*/ 
typedef Hidden_type Algorithm; 

/*! 
is a typedef to the template parameter `FT` 
*/ 
typedef Hidden_type Radius; 

/*! 
is a typedef to `K::Point_2`. 
*/ 
typedef Hidden_type Point; 

/*! 
is a typedef to 
`Point`. 
*/ 
typedef Hidden_type Sphere; 

/*! 
is a typedef to 
`K::Cartesian_const_iterator_2`. 
*/ 
typedef Hidden_type Cartesian_const_iterator; 

/// @} 

/// \name Access Functions 
/// The class provides the access functions required by the concept 
/// `MinSphereOfSpheresTraits`; they simply map to the corresponding 
/// routines of class `K::Point_2`: 
/// @{

/*! 
returns `0`. 
*/ 
FT radius(const 
Sphere& s); 

/*! 
maps to `s.first.cartesian_begin()`. 
*/ 
Cartesian_const_iterator center_cartesian_begin(const 
Sphere& s); 

/// @}

}; /* end Min_sphere_of_points_d_traits_2 */
} /* end namespace CGAL */

