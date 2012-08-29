/// \addtogroup PkgAlgebraicFoundations Algebraic Foundations
/// @{
/// \addtogroup PkgAlgebraicFoundationsConcepts Concepts
/// @{

 
///  
///  A model of `FractionTraits` is associated with a type `Type`.
///  In case the associated type is a `Fraction`, a model of `FractionTraits` provides the relevant functionality for decomposing and re-composing as well
///  as the numerator and denominator type.
///  
///  
///  
///  Functors 
///  -------------- 
///  
///  In case `Type` is not a `Fraction` all functors are `Null_functor`.
///  
///  \hasModel CGAL::Fraction_traits
///
///  \sa `FractionTraits::Decompose`
///  \sa `FractionTraits::Compose`
///  \sa `FractionTraits::CommonFactor`
class FractionTraits {
public:

/// \name Types
/// @{
/*!
 the associated type
*/
typedef Hidden_type Type;
/// @}

/// \name Types
/// @{
/*!
  
        Tag indicating whether the associated type is a fraction and can be 
        decomposed into a numerator and denominator. 
        This is either `CGAL::Tag_true` or `CGAL::Tag_false`.
*/
typedef Hidden_type Is_fraction;
/// @}

/// \name Types
/// @{
/*!
 The type to represent the numerator. 
        This is undefined in case the associated type is not a fraction. 
*/
typedef Hidden_type Numerator_type  ;
/// @}

/// \name Types
/// @{
/*!
 The (simpler) type to represent the denominator.
        This is undefined in case the associated type is not a fraction. 
*/
typedef Hidden_type Denominator_type;
/// @}

/// \name Functors
/// @{
/*!
 A model of `FractionTraits::Decompose`.
*/
typedef Hidden_type Decompose;
/// @}

/// \name Functors
/// @{
/*!
 A model of `FractionTraits::Compose`.
*/
typedef Hidden_type Compose;
/// @}

/// \name Functors
/// @{
/*!
 A model of `FractionTraits::CommonFactor`.
*/
typedef Hidden_type Common_factor;
/// @}

}; /* concept FractionTraits */
/// @}
/// @}

