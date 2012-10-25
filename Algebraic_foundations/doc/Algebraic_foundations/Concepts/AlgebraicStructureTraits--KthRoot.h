
/*!
\ingroup PkgAlgebraicFoundationsAlgebraicStructuresConcepts
\cgalconcept

`AdaptableBinaryFunction` providing the k-th root. 

\refines `AdaptableBinaryFunction` 

\sa ::FieldWithRootOf 
\sa ::AlgebraicStructureTraits 

*/

class AlgebraicStructureTraits::KthRoot {
public:

/// \name Types 
/// @{

/*! 
Is `AlgebraicStructureTraits::Type`. 
*/ 
typedef Hidden_type result_type; 

/*! 
Is int. 
*/ 
typedef Hidden_type first_argument; 

/*! 
Is `AlgebraicStructureTraits::Type`. 
*/ 
typedef Hidden_type second_argument; 

/// @} 

/// \name Operations 
/// @{

/*! 
returns the \f$ k\f$-th root of \f$ x\f$. 
\pre \f$ k \geq1\f$ 

*/ 
result_type operator()(int k, second_argument_type x); 

/// @}

}; /* end AlgebraicStructureTraits::KthRoot */

