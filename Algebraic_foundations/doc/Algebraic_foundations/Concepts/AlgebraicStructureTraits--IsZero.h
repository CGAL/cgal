
/*!
\ingroup PkgAlgebraicFoundationsAlgebraicStructuresConcepts
\cgalconcept

`AdaptableUnaryFunction`, returns true in case the argument is the zero element of the ring. 

\refines ::AdaptableUnaryFunction 

\sa ::AlgebraicStructureTraits 
\sa ::RealEmbeddableTraits::IsZero 

*/

class AlgebraicStructureTraits::IsZero {
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

returns true in case \f$ x\f$ is the zero element of the ring. 
*/ 
result_type operator()(argument_type x) const; 

/// @}

}; /* end AlgebraicStructureTraits::IsZero */

