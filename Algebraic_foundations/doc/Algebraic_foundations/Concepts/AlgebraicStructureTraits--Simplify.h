
namespace AlgebraicStructureTraits_{

/*!
\ingroup PkgAlgebraicFoundationsAlgebraicStructuresConcepts
\cgalConcept

This `AdaptableUnaryFunction` may simplify a given object. 

\cgalRefines `AdaptableUnaryFunction` 

\sa `AlgebraicStructureTraits`

*/

class Simplify {
public:

/// \name Types 
/// @{

/*!
Is void. 
*/ 
typedef unspecified_type result_type; 

/*!
Is `AlgebraicStructureTraits::Type`. 
*/ 
typedef unspecified_type argument_type; 

/// @} 

/// \name Operations 
/// @{

/*!
may simplify \f$ x\f$. 
*/ 
result_type operator()(argument_type x); 

/// @}

}; /* end Simplify */

} /* end of namespace AlgebraicStructureTraits_ */
