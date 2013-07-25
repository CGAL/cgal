
namespace RealEmbeddableTraits_ {

/*!
\ingroup PkgAlgebraicFoundationsAlgebraicStructuresConcepts
\cgalConcept

This `AdaptableUnaryFunction` computes the sign of a real embeddable number. 

\cgalRefines `AdaptableUnaryFunction` 

\sa `RealEmbeddableTraits`

*/

class Sgn {
public:

/// \name Types 
/// @{

/*!
Type convertible to `CGAL::Sign`. 
*/ 
typedef unspecified_type result_type; 

/*!
Is `RealEmbeddableTraits::Type`. 
*/ 
typedef unspecified_type argument_type; 

/// @} 

/// \name Operations 
/// @{

/*!
Computes the sign of \f$ x\f$. 
*/ 
result_type operator()(argument_type x); 

/// @}

}; /* end Sgn */

} /* end of namespace RealEmbeddableTraits_ */
