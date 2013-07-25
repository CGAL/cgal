
namespace AlgebraicStructureTraits_{

/*!
\ingroup PkgAlgebraicFoundationsAlgebraicStructuresConcepts
\cgalConcept

`AdaptableUnaryFunction` providing the inverse element with 
respect to multiplication of a `Field`. 

\cgalRefines `AdaptableUnaryFunction` 

\sa `AlgebraicStructureTraits`

*/

class Inverse {
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
returns the inverse element of \f$ x\f$ with respect to multiplication. 
\pre \f$ x \neq0\f$. 

*/ 
result_type operator()(argument_type x) const; 

/// @}

}; /* end Inverse */

} /* end of namespace AlgebraicStructureTraits_ */
