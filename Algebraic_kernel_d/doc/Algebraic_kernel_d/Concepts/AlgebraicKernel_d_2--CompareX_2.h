
/*!
\ingroup PkgAlgebraicKerneldConceptsBi
\cgalConcept

Compares the first coordinates of `AlgebraicKernel_d_2::Algebraic_real_2`s. 

\cgalRefines `AdaptableBinaryFunction` 

\sa `AlgebraicKernel_d_2::CompareY_2`
\sa `AlgebraicKernel_d_2::CompareXY_2`

*/

class AlgebraicKernel_d_2::CompareX_2 {
public:

/// \name Types 
/// @{

/*!
Type convertible to `CGAL::Comparison_result` 
*/ 
typedef unspecified_type result_type; 

/*!

*/ 
typedef AlgebraicKernel_d_2::Algebraic_real_2 first_argument_type; 

/*!

*/ 
typedef AlgebraicKernel_d_2::Algebraic_real_2 second_argument_type; 

/// @} 

/// \name Operations 
/// The following operators and their symmetric counterparts are required:
/// @{

/*!
Compares the first coordinates of \f$ a\f$ and \f$ b\f$. 
*/ 
result_type operator()(const first_argument_type & a, const second_argument_type & b); 

/*!
Compares the first coordinate of \f$ a\f$ with \f$ x\f$. 
*/ 
result_type operator()(AlgebraicKernel_d_2::Algebraic_real_2 a, int x); 

/*!
Compares the first coordinate of \f$ a\f$ with \f$ x\f$. 
*/ 
result_type operator()(AlgebraicKernel_d_2::Algebraic_real_2 a, AlgebraicKernel_d_2::Bound x); 

/*!
Compares the first coordinate of \f$ a\f$ with \f$ x\f$. 
*/ 
result_type operator()(AlgebraicKernel_d_2::Algebraic_real_2 a, AlgebraicKernel_d_2::Coefficient x); 

/*!
Compares the first coordinate of \f$ a\f$ with \f$ x\f$. 
*/ 
result_type operator()(AlgebraicKernel_d_2::Algebraic_real_2 a, AlgebraicKernel_d_2::Algebraic_real_1 x); 

/// @}

}; /* end AlgebraicKernel_d_2::CompareX_2 */

