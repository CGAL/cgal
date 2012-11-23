
/*!
\ingroup PkgAlgebraicFoundationsConcepts
\cgalConcept

This `AdaptableUnaryFunction` computes the sign of a real embeddable number. 

\cgalRefines `AdaptableUnaryFunction` 

\sa `RealEmbeddableTraits`

*/

class RealEmbeddableTraits::Sgn {
public:

/// \name Types 
/// @{

/*! 
Type convertible to `CGAL::Sign`. 
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
Computes the sign of \f$ x\f$. 
*/ 
result_type operator()(argument_type x); 

/// @}

}; /* end RealEmbeddableTraits::Sgn */

