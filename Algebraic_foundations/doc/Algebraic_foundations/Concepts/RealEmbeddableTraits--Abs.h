
/*!
\ingroup PkgAlgebraicFoundationsRealEmbeddableConcepts
\cgalConcept

`AdaptableUnaryFunction` computes the absolute value of a number. 

\cgalRefines `AdaptableUnaryFunction` 

\sa ::RealEmbeddableTraits 

*/

class RealEmbeddableTraits::Abs {
public:

/// \name Types 
/// @{

/*! 
Is `RealEmbeddableTraits::Type`. 
*/ 
typedef Hidden_type result_type; 

/*! 
Is `RealEmbeddableTraits::Type`. 
*/ 
typedef Hidden_type argument_type; 

/// @} 

/// \name Operations 
/// @{

/*! 
computes the absolute value of \f$ x\f$. 
*/ 
result_type operator()(argument_type x); 

/// @}

}; /* end RealEmbeddableTraits::Abs */

