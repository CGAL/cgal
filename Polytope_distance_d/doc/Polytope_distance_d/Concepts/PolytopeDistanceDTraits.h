
/*!
\ingroup PkgOptimalDistancesConcepts
\cgalConcept

This concept defines the requirements for traits classes of \f$ d\f$-dimensional 
optimisation algorithms. 

\cgalHasModel `CGAL::Polytope_distance_d_traits_2<K,ET,NT>`
\cgalHasModel `CGAL::Polytope_distance_d_traits_3<K,ET,NT>`
\cgalHasModel `CGAL::Polytope_distance_d_traits_d<K,ET,NT>`

\sa `CGAL::Polytope_distance_d<Traits>` 

*/

class PolytopeDistanceDTraits {
public:

/// \name Types 
/// @{

/*!
point type used to represent the input points. 
*/ 
typedef unspecified_type Point_d; 

/*!
compile time tag to distinguish between %Cartesian and homogeneous 
representation of the input points. `Rep_tag` has to be either 
`CGAL::``Cartesian_tag` or 
`CGAL::``Homogeneous_tag`. 
*/ 
typedef unspecified_type Rep_tag; 

/*!
number type used to represent the coordinates of the input points. 
It has to be a model for `RingNumberType`. 
*/ 
typedef unspecified_type RT; 

/*!
number type used to return either the squared radius of the smallest 
enclosing sphere or annulus, or the squared distance of the polytopes. 
`FT` has to be either `RT` or `CGAL::``Quotient<RT>` if 
the input points have %Cartesian or homogeneous representation, 
respectively (cf. `Rep_tag`). 
*/ 
typedef unspecified_type FT; 

/*!
data accessor object used to access the dimension of the input points. 
*/ 
typedef unspecified_type Access_dimension_d; 

/*!
data accessor object used to access the coordinates of the input points. 
*/ 
typedef unspecified_type Access_coordinates_begin_d; 

/*!
constructor object used to construct either the center of the smallest 
enclosing sphere or annulus, or the points realizing the distance between 
the two polytopes. 
*/ 
typedef unspecified_type Construct_point_d; 

/// @}

/// \name Special types
/// The following two number types are only needed for
/// `CGAL::Min_annulus_d<Traits>` and
/// `CGAL::Polytope_distance_d<Traits>`.
/// @{

/*!
exact number type used to do the exact computations in the 
underlying solver for linear programs. It has to to be a model for 
`RingNumberType`. There must be an implicit conversion from 
`RT` to `ET` available. 
*/ 
typedef unspecified_type ET; 

/*!
fast (possibly inexact) number type used to speed up the pricing step in 
the underlying solver for linear programs. It has to be a model for 
`RingNumberType`. There must be implicit conversions from `RT` to 
`NT` and from `NT` to `ET` available. 
*/ 
typedef unspecified_type NT; 

/// @} 

/// \name Creation 
/// Only default and copy constructor are required.
/// @{

/*!

*/ 
PolytopeDistanceDTraits( ); 

/*!

*/ 
PolytopeDistanceDTraits( const PolytopeDistanceDTraits&); 

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

}; /* end PolytopeDistanceDTraits */

