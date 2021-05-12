
/*!
\ingroup PkgModularArithmeticConcepts
\cgalConcept

This `AdaptableUnaryFunction` returns a representative in the original type of a given modular image. More precisely, it implements the *right inverse* of a proper restriction of the homomorphism \f$ \varphi\f$, which is implemented by `ModularTraits::ModularImage`.

\cgalRefines `AdaptableUnaryFunction`

\sa `ModularTraits`

*/

class ModularTraits::ModularImageRepresentative {
public:

/// \name Types
/// @{

/*!

*/
typedef ModularTraits::Type result_type;

/*!

*/
typedef ModularTraits::Residue_type argument_type;

/*!

computes \f$ \varphi^{-1}(x)\f$.

*/
result_type
operator()(const argument_type &x);

/// @}

}; /* end ModularTraits::ModularImageRepresentative */

