
/*!
\ingroup PkgAlgebraicFoundationsFractionsConcepts
\cgalconcept

A model of `FractionTraits` is associated with a type `Type`. 

In case the associated type is a `Fraction`, a model of `FractionTraits` provides the relevant functionality for decomposing and re-composing as well 
as the numerator and denominator type. 

\hasModel `CGAL::Fraction_traits<T>`

\sa `FractionTraits::Decompose` 
\sa `FractionTraits::Compose` 
\sa `FractionTraits::CommonFactor` 

*/
class FractionTraits {
public:

/// \name Types 
/// @{

/*! 
the associated type 
*/ 
typedef Hidden_type Type; 

/*! 

Tag indicating whether the associated type is a fraction and can be 
decomposed into a numerator and denominator. 

This is either `CGAL::Tag_true` or `CGAL::Tag_false`. 
*/ 
typedef Hidden_type Is_fraction; 

/*! 
The type to represent the numerator. 
This is undefined in case the associated type is not a fraction. 
*/ 
typedef Hidden_type Numerator_type ; 

/*! 
The (simpler) type to represent the denominator. 
This is undefined in case the associated type is not a fraction. 
*/ 
typedef Hidden_type Denominator_type; 

/// @} 

/// \name Functors 
/// In case `Type` is not a `Fraction` all functors are `Null_functor`.
/// @{


/*!
\ingroup PkgAlgebraicFoundationsFractionsConcepts
\cgalconcept

Functor decomposing a `Fraction` into its numerator and denominator. 

\sa `Fraction` 
\sa `FractionTraits` 
\sa `FractionTraits::Compose` 
\sa `FractionTraits::CommonFactor` 

*/

class Decompose {
public:

/// \name Operations 
/// @{

/*! 
decompose \f$ f\f$ into numerator \f$ n\f$ and denominator \f$ d\f$. 
*/ 
void operator()( FractionTraits::Type f, 
FractionTraits::Numerator_type & n, 
FractionTraits::Denominator_type & d); 

/// @}

}; /* end FractionTraits::Decompose */

/*!
\ingroup PkgAlgebraicFoundationsFractionsConcepts
\cgalconcept

`AdaptableBinaryFunction`, returns the fraction of its arguments. 

\refines ::AdaptableBinaryFunction 

\sa `Fraction` 
\sa `FractionTraits` 
\sa `FractionTraits::Decompose` 
\sa `FractionTraits::CommonFactor` 

*/

class Compose {
public:

/// \name Types 
/// @{

/*! 

*/ 
typedef FractionTraits::Type result_type; 

/*! 

*/ 
typedef FractionTraits::Numerator_type first_argument_type; 

/*! 

*/ 
typedef FractionTraits::Denominator_type second_argument_type; 

/// @} 

/// \name Operations 
/// @{

/*! 
return the fraction \f$ n/d\f$. 
*/ 
result_type operator()(first_argument_type n, second_argument_type d); 

/// @}

}; /* end FractionTraits::Compose */


/*!
\ingroup PkgAlgebraicFoundationsFractionsConcepts
\cgalconcept

`AdaptableBinaryFunction`, finds great common factor of denominators. 

This can be considered as a relaxed version of `AlgebraicStructureTraits::Gcd`, 
this is needed because it is not guaranteed that `FractionTraits::Denominator_type` is a model of 
`UniqueFactorizationDomain`. 

\refines ::AdaptableBinaryFunction 

\sa `Fraction` 
\sa `FractionTraits` 
\sa `FractionTraits::Decompose` 
\sa `FractionTraits::Compose` 
\sa `AlgebraicStructureTraits::Gcd` 

*/

class CommonFactor {
public:

/// \name Types 
/// @{

/*! 

*/ 
typedef FractionTraits::Denominator_type result_type; 

/*! 

*/ 
typedef FractionTraits::Denominator_type first_argument_type; 

/*! 

*/ 
typedef FractionTraits::Denominator_type second_argument_type; 

/// @} 

/// \name Operations 
/// @{

/*! 
return a great common factor of \f$ d1\f$ and \f$ d2\f$. 

Note: <TT>operator()(0,0) = 0</TT> 
*/ 
result_type operator()(first_argument_type d1, second_argument_type d2); 

/// @}

}; /* end FractionTraits::CommonFactor */

  /*!
    A model of FractionTraits::Compose.
   */
typedef Hidden_type Compose;


  /*!
    A model of FractionTraits::Decompose.
   */
typedef Hidden_type Decompose;


  /*!
    A model of FractionTraits::CommonFactor.
   */
typedef Hidden_type Common_factor;


/// @}

}; /* end FractionTraits */
