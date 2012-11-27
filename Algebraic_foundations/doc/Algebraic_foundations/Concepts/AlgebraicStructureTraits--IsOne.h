
/*!
\ingroup PkgAlgebraicFoundationsAlgebraicStructuresConcepts
\cgalConcept

`AdaptableUnaryFunction`, 
returns true in case the argument is the one of the ring. 

\cgalRefines `AdaptableUnaryFunction` 

\sa `AlgebraicStructureTraits`

*/

class AlgebraicStructureTraits::IsOne {
public:

/// \name Types 
/// @{

/*! 
Is `AlgebraicStructureTraits::Boolean`. 
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

returns true in case \f$ x\f$ is the one of the ring. 
*/ 
result_type operator()(argument_type x); 

/// @}

}; /* end AlgebraicStructureTraits::IsOne */

