
namespace RealEmbeddableTraits_ {

/*!
\ingroup PkgAlgebraicFoundationsAlgebraicStructuresConcepts
\cgalConcept

`AdaptableUnaryFunction` computes a double approximation of a real 
embeddable number. 

Remark: In order to control the quality of approximation one has to resort 
to methods that are specific to NT. There are no general guarantees whatsoever. 

\cgalRefines `AdaptableUnaryFunction` 

\sa `RealEmbeddableTraits`

*/

class ToDouble {
public:

/// \name Types 
/// @{

/*!
The result type. 
*/ 
typedef double result_type;

/*!
Is `RealEmbeddableTraits::Type`. 
*/ 
typedef unspecified_type argument_type; 

/// @} 

/// \name Operations 
/// @{

/*!
computes a double approximation of a real embeddable number. 
*/ 
result_type operator()(argument_type x); 

/// @}

}; /* end ToDouble */

} /* end of namespace RealEmbeddableTraits_ */
