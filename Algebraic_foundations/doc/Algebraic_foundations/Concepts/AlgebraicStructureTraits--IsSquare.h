
/*!
\ingroup PkgAlgebraicFoundationsAlgebraicStructuresConcepts
\cgalConcept

`AdaptableBinaryFunction` that computes whether the first argument is a square. 
If the first argument is a square the second argument, which is taken by reference, contains the square root. 
Otherwise, the content of the second argument is undefined. 

A ring element \f$ x\f$ is said to be a square iff there exists a ring element \f$ y\f$ such 
that \f$ x= y*y\f$. In case the ring is a `UniqueFactorizationDomain`, 
\f$ y\f$ is uniquely defined up to multiplication by units. 

\cgalRefines `AdaptableBinaryFunction` 

\sa `AlgebraicStructureTraits`

*/

class AlgebraicStructureTraits::IsSquare {
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
typedef Hidden_type first_argument; 

/*! 
Is `AlgebraicStructureTraits::Type&`. 
*/ 
typedef Hidden_type second_argument; 

/// @} 

/// \name Operations 
/// @{

/*! 
returns <TT>true</TT> in case \f$ x\f$ is a square, i.e.\ \f$ x = y*y\f$. 
\post \f$ unit\_part(y) == 1\f$. 

*/ 
result_type operator()(first_argument_type x, 
second_argument_type y); 

/*! 
returns <TT>true</TT> in case \f$ x\f$ is a square. 

*/ 
result_type operator()(first_argument_type x); 

/// @}

}; /* end AlgebraicStructureTraits::IsSquare */

