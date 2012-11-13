
/*!
\ingroup PkgAlgebraicFoundationsAlgebraicStructuresConcepts
\cgalConcept

`AdaptableUnaryFunction` providing the square root. 

\cgalRefines `AdaptableUnaryFunction` 

\sa ::AlgebraicStructureTraits 

*/

class AlgebraicStructureTraits::Sqrt {
public:

/// \name Types 
/// @{

/*! 
Is `AlgebraicStructureTraits::Type`. 
*/ 
typedef Hidden_type result_type; 

/*! 
Is `AlgebraicStructureTraits::Type`. 
*/ 
typedef Hidden_type argument_type; 

/// @} 

/// \name Operations 
/// @{

/*! 
returns \f$ \sqrt{x}\f$. 
*/ 
result_type operator()(argument_type x) const; 

/// @}

}; /* end AlgebraicStructureTraits::Sqrt */

