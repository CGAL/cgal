
namespace AlgebraicStructureTraits_{

/*!
\ingroup PkgAlgebraicFoundationsAlgebraicStructuresConcepts
\cgalConcept

`AdaptableUnaryFunction`, computing the square of the argument. 

\cgalRefines `AdaptableUnaryFunction` 

\sa `AlgebraicStructureTraits`

*/

class Square {
public:

/// \name Types 
/// @{

/*!
Is `AlgebraicStructureTraits::Type`. 
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
returns the square of \f$ x\f$. 
*/ 
result_type operator()(argument_type x); 

/// @}

}; /* end Square */

} /* end of namespace AlgebraicStructureTraits_ */
