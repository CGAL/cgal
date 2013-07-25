
/*!
\ingroup PkgAlgebraicKerneldConceptsUni
\cgalConcept

Compares `AlgebraicKernel_d_1::Algebraic_real_1` values. 

\cgalRefines `AdaptableBinaryFunction` 

*/
class AlgebraicKernel_d_1::Compare_1 {
public:

/// \name Types 
/// @{

/*!
Type convertible to `CGAL::Comparison_result` 
*/ 
typedef unspecified_type result_type; 

/*!

*/ 
typedef AlgebraicKernel_d_1::Algebraic_real_1 first_argument_type; 

/*!

*/ 
typedef AlgebraicKernel_d_1::Algebraic_real_1 second_argument_type; 

/// @} 

/// \name Operations 
/// The following operators and their symmetric counterparts are
/// required:
/// @{

/*!
Compares `a` and `b`. 
*/ 
result_type operator()(AlgebraicKernel_d_1::Algebraic_real_1 a, AlgebraicKernel_d_1::Algebraic_real_1 b); 

/*!
Compares `a` and `b`. 
*/ 
result_type operator()(AlgebraicKernel_d_1::Algebraic_real_1 a, int b); 

/*!
Compares `a` and `b`. 
*/ 
result_type operator()(AlgebraicKernel_d_1::Algebraic_real_1 a, AlgebraicKernel_d_1::Bound b); 

/*!
Compares `a` and `b`. 
*/ 
result_type operator()(AlgebraicKernel_d_1::Algebraic_real_1 a, AlgebraicKernel_d_1::Coefficient b); 

/// @}

}; /* end AlgebraicKernel_d_1::Compare_1 */

