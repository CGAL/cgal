
namespace CGAL {

/*!
\ingroup PkgAlgebraicKerneldModels

\anchor Algebraic_kernel_rs_gmpz_d_1 

This univariate algebraic kernel uses the Rs library to perform 
integer univariate polynomial root isolation. It is a model of the 
`AlgebraicKernel_d_1` concept. 

\cgalModels `AlgebraicKernel_d_1`

\sa `Algebraic_kernel_rs_gmpz_d_1` 

*/

class Algebraic_kernel_rs_gmpz_d_1 {
public:

/// \name Types 
/// @{

/*!
It is a typedef to `CGAL::Gmpz`. 
*/ 
typedef unspecified_type Coefficient; 

/*!
It is defined as `CGAL::Polynomial<CGAL::Gmpz>`. 
*/ 
typedef unspecified_type Polynomial_1; 

/*!
Type that represents the real roots of 
integer univariate polynomials, containing a pointer to the polynomial of 
which the represented algebraic number is root and and a `CGAL::Gmpfi` 
isolating interval. 
*/ 
typedef unspecified_type Algebraic_real_1; 

/*!
Since the isolating intervals of the roots have type 
`CGAL::Gmpfi`, the bounds have type `CGAL::Gmpfr`. 
*/ 
typedef unspecified_type Bound; 

/*!
The multiplicity is an `int`. 
*/ 
typedef unspecified_type Multiplicity_type; 

/// @}

}; /* end Algebraic_kernel_rs_gmpz_d_1 */
} /* end namespace CGAL */
