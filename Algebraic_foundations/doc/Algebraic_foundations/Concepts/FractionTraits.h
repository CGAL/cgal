
/*!
\ingroup PkgAlgebraicFoundationsAlgebraicStructuresConcepts
\cgalConcept

A model of `FractionTraits` is associated with a type `Type`. 

In case the associated type is a `Fraction`, a model of `FractionTraits` provides the relevant functionality for decomposing and re-composing as well 
as the numerator and denominator type. 

\cgalHasModel `CGAL::Fraction_traits<T>`

\sa `FractionTraits_::Decompose` 
\sa `FractionTraits_::Compose` 
\sa `FractionTraits_::CommonFactor` 

*/
class FractionTraits {
public:

/// \name Types 
/// @{

/*!
The associated type 
*/ 
typedef unspecified_type Type; 

/*!

Tag indicating whether the associated type is a fraction and can be 
decomposed into a numerator and denominator. 

This is either `CGAL::Tag_true` or `CGAL::Tag_false`. 
*/ 
typedef unspecified_type Is_fraction; 

/*!
The type to represent the numerator. 
This is undefined in case the associated type is not a fraction. 
*/ 
typedef unspecified_type Numerator_type ; 

/*!
The (simpler) type to represent the denominator. 
This is undefined in case the associated type is not a fraction. 
*/ 
typedef unspecified_type Denominator_type; 

/// @} 

/// \name Functors 
/// In case `Type` is not a `Fraction` all functors are `Null_functor`.
/// @{

  /*!
    A model of FractionTraits_::Compose.
   */
typedef unspecified_type Compose;


  /*!
    A model of FractionTraits_::Decompose.
   */
typedef unspecified_type Decompose;


  /*!
    A model of FractionTraits_::CommonFactor.
   */
typedef unspecified_type Common_factor;


/// @}

}; /* end FractionTraits */

namespace FractionTraits_ {

/*!
\ingroup PkgAlgebraicFoundationsAlgebraicStructuresConcepts
\cgalConcept

Functor decomposing a `Fraction` into its numerator and denominator. 

\sa `Fraction` 
\sa `FractionTraits` 
\sa `FractionTraits_::Compose` 
\sa `FractionTraits_::CommonFactor` 

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

}; /* end Decompose */

/*!
\ingroup PkgAlgebraicFoundationsAlgebraicStructuresConcepts
\cgalConcept

`AdaptableBinaryFunction`, returns the fraction of its arguments. 

\cgalRefines `AdaptableBinaryFunction` 

\sa `Fraction` 
\sa `FractionTraits` 
\sa `FractionTraits_::Decompose` 
\sa `FractionTraits_::CommonFactor` 

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

}; /* end Compose */


/*!
\ingroup PkgAlgebraicFoundationsAlgebraicStructuresConcepts
\cgalConcept

`AdaptableBinaryFunction`, finds great common factor of denominators. 

This can be considered as a relaxed version of `AlgebraicStructureTraits_::Gcd`, 
this is needed because it is not guaranteed that `FractionTraits::Denominator_type` is a model of 
`UniqueFactorizationDomain`. 

\cgalRefines `AdaptableBinaryFunction` 

\sa `Fraction` 
\sa `FractionTraits` 
\sa `FractionTraits_::Decompose` 
\sa `FractionTraits_::Compose` 
\sa `AlgebraicStructureTraits_::Gcd` 

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

}; /* end CommonFactor */

} /* end of namespace FractionTraits_ */
