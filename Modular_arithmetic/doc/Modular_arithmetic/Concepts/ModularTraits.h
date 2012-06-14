/// \ingroup PkgModularArithmeticConcepts
/// 
/// A model of `ModularTraits` is associated to a specific `Type`.
/// In case this associated type is `Modularizable`, this is indicated by the 
/// boolean tag `Is_modularizable`. The mapping into the `Residue_type` is
/// provided by the functor `Modular_image`.
/// 
/// A model of `ModularTraits` is supposed to provide:
/// 
/// \hasModel CGAL::Modular_traits
///
/// \sa `CGAL::Residue`
/// \sa `Modularizable`
///  
class ModularTraits {
public:

/// \ingroup PkgModularArithmeticConcepts
/// 
/// This `AdaptableUnaryFunction` computes the modular image of the given value 
/// with respect to a homomorphism \f$\phi\f$ from the 
/// `ModularTraits::Type` into the `ModularTraits::Residue_type`.
/// The homomorphism preserves the mapping of `int` into both types
/// , i.e., \f$\phi(Type(i)) == Residue_type(i)\f$.
/// 
/// \refines ::AdaptableUnaryFunction
///
/// \sa `ModularTraits`
class ModularImage {
public:

/// \name Types
/// @{
/*!
 
*/
typedef ModularTraits::Residue_type result_type;
/// @}

/// \name Types
/// @{
/*!
 
*/
typedef ModularTraits::Type         argument_type;
/// @}

/// \name Types
/// @{
/*!
 computes \f$\phi(x)\f$. 
*/

result_type
operator()(const  argument_type & x);
/// @}

}; /* concept ModularTraits::ModularImage */

/// \ingroup PkgModularArithmeticConcepts
/// 
/// This `AdaptableUnaryFunction` returns a representative in the
/// original type of a given modular image. More precisely, it
/// implements the \f$right inverse\f$ of a proper restriction of the
/// homomorphism \f$\phi\f$, which is implemented by
/// `ModularTraits::ModularImage`.
/// 
/// 
/// \refines ::AdaptableUnaryFunction
/// \sa `ModularTraits`
///  
class ModularImageRepresentative {
public:

/// \name Types
/// @{
/*!
 
*/
typedef ModularTraits::Type         result_type;
/// @}

/// \name Types
/// @{
/*!
 
*/
typedef ModularTraits::Residue_type argument_type;
/// @}

/// \name Types
/// @{
/*!
  computes \f$\phi^{-1}(x)\f$. 
*/
result_type
operator()(const  argument_type &x);
/// @}

}; /* concept ModularTraits::ModularImageRepresentative */

/// \name Types
/// @{
/*!
 The associated type.
*/
typedef Hidden_type Type        ;
/// @}

/// \name Types
/// @{
/*!
  
  Tag indicating whether the associated type is modularizable. 
  This is either `CGAL::Tag_true` or `CGAL::Tag_false`. 
*/
typedef Hidden_type Is_modularizable;
/// @}

/// \name Types
/// @{
/*!
 
   The type of the modular image. 
   In case the type is not `Modularizable` this is undefined. 

*/
typedef Hidden_type Residue_type;
/// @}

/// \name Functors
/// In case the associated type is `Modularizable` all functors are provided.
/// In case a functor is not provided, it is set to `CGAL::Null_functor`.
/// @{

/*!
 A model of `ModularTraits::ModularImage` 
*/
typedef Hidden_type Modular_image;

/*!
 A model of `ModularTraits::ModularImageRepresentative` 
*/
typedef Hidden_type Modular_image_representative;
/// @}

}; /* concept ModularTraits */
///

                   
  

