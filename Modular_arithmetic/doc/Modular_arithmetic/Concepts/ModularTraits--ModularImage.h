
/*!
\ingroup PkgModularArithmeticConcepts
\cgalConcept

This `AdaptableUnaryFunction` computes the modular image of the given value 
with respect to a homomorphism \f$ \varphi\f$ from the 
`ModularTraits::Type` into the `ModularTraits::Residue_type`. 

The homomorphism preserves the mapping of `int` into both types 
, i.e., \f$ \varphi(\mathrm{Type}(i)) == \mathrm{Residue\_type}(i)\f$. 

\cgalRefines `AdaptableUnaryFunction` 

\sa `ModularTraits` 

*/

class ModularTraits::ModularImage {
public:

/// \name Types 
/// @{

/*!

*/ 
typedef ModularTraits::Residue_type result_type; 

/*!

*/ 
typedef ModularTraits::Type argument_type; 

/*!

computes \f$ \varphi(x)\f$. 

*/ 
result_type 
operator()(const argument_type & x); 

/// @}

}; /* end ModularTraits::ModularImage */

