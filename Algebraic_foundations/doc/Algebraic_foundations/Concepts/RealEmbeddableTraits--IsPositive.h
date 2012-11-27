
/*!
\ingroup PkgAlgebraicFoundationsConcepts
\cgalConcept

`AdaptableUnaryFunction`, returns true in case the argument is positive. 

\cgalRefines `AdaptableUnaryFunction` 

\sa `RealEmbeddableTraits`

*/

class RealEmbeddableTraits::IsPositive {
public:

/// \name Types 
/// @{

/*! 
Type convertible to `bool`. 
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
returns true in case \f$ x\f$ is positive. 
*/ 
result_type operator()(argument_type x); 

/// @}

}; /* end RealEmbeddableTraits::IsPositive */

