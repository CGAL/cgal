
/*!
\ingroup PkgAlgebraicKerneldConceptsBi
\cgalConcept

Compares the second coordinated of `AlgebraicKernel_d_2::Algebraic_real_2`s. 

\cgalRefines `AdaptableBinaryFunction` 

\sa `AlgebraicKernel_d_2::CompareX_2`
\sa `AlgebraicKernel_d_2::CompareXY_2`

*/

class AlgebraicKernel_d_2::CompareY_2 {
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
Compares the second coordinates of \f$ a\f$ and \f$ b\f$. 
*/ 
result_type operator()(const first_argument_type & a, const second_argument_type & b); 

/*!
Compares the second coordinate of \f$ a\f$ with \f$ y\f$. 
*/ 
result_type operator()(AlgebraicKernel_d_2::Algebraic_real_2 a, int y); 

/*!
Compares the second coordinate of \f$ a\f$ with \f$ y\f$. 
*/ 
result_type operator()(AlgebraicKernel_d_2::Algebraic_real_2 a, AlgebraicKernel_d_2::Bound y); 

/*!
Compares the second coordinate of \f$ a\f$ with \f$ y\f$. 
*/ 
result_type operator()(AlgebraicKernel_d_2::Algebraic_real_2 a, AlgebraicKernel_d_2::Coefficient y); 

/*!
Compares the second coordinate of \f$ a\f$ with \f$ y\f$. 
*/ 
result_type operator()(AlgebraicKernel_d_2::Algebraic_real_2 a, AlgebraicKernel_d_2::Algebraic_real_1 y); 

/// @}

}; /* end AlgebraicKernel_d_2::CompareY_2 */

