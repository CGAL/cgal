
namespace AlgebraicStructureTraits_{

/*!
\ingroup PkgAlgebraicFoundationsAlgebraicStructuresConcepts
\cgalConcept

This `AdaptableUnaryFunction` computes the unit part of a given ring 
element. 

The mathematical definition of unit part is as follows: Two ring elements \f$ a\f$ 
and \f$ b\f$ are said to be associate if there exists an invertible ring element 
(i.e.\ a unit) \f$ u\f$ such that \f$ a = ub\f$. This defines an equivalence relation. 
We can distinguish exactly one element of every equivalence class as being 
unit normal. Then each element of a ring possesses a factorization into a unit 
(called its unit part) and a unit-normal ring element 
(called its unit normal associate). 

For the integers, the non-negative numbers are by convention unit normal, 
hence the unit-part of a non-zero integer is its sign. For a `Field`, every 
non-zero element is a unit and is its own unit part, its unit normal 
associate being one. The unit part of zero is, by convention, one. 

\cgalRefines `AdaptableUnaryFunction` 

\sa `AlgebraicStructureTraits`

*/


class UnitPart {
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
returns the unit part of \f$ x\f$. 
*/ 
result_type operator()(argument_type x); 

/// @}

}; /* end UnitPart */

} /* end of namespace AlgebraicStructureTraits_ */
