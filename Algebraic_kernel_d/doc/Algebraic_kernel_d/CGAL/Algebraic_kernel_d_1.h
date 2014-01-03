
namespace CGAL {

/*!
\ingroup PkgAlgebraicKerneldModels

The class represents an algebraic real root by a square free polynomial and an 
isolating interval that uniquely defines the root. 
The template argument `Coeff` determines the coefficient type of the 
kernel, which is also the coefficient type of the supported polynomials. 

Currently, the following coefficient types are supported: 

- `Gmpz`, `Gmpq`, (requires configuration with external libraries GMP, MPFR and MPFI) 

- `CORE::BigInt`, `CORE::BigRat`, (requires configuration with external library GMP) 

- `leda_integer`, `leda_rational`. (requires configuration with external library LEDA) 

\cgalAdvancedBegin
The template argument type can also be set to `Sqrt_extension<NT,ROOT>`, where `NT` is 
one of the types listed above. `ROOT` should be one of the integer types. 
See also the documentation of `Sqrt_extension<NT,ROOT>`. 
\cgalAdvancedEnd

The current method to isolate roots is the bitstream Descartes method
presented in \cgalCite{eigenwillig-phd-08}.  The used method to refine the
approximation of an algebraic real root is a slightly modified
(filtered) version of the one presented in \cgalCite{abbott-qir-06}. The
method has quadratic convergence.

\cgalModels `AlgebraicKernel_d_1`

\sa `AlgebraicKernel_d_1` 
\sa `Polynomial_d` 
\sa `CGAL::Algebraic_kernel_d_2<Coeff>`

*/
template< typename Coeff >
class Algebraic_kernel_d_1 {
public:

/// \name Types 
/// @{

/*!
Same type as the template argument `Coeff`. 
*/ 
typedef unspecified_type Coefficient; 

/*!
A model of `::AlgebraicKernel_d_1::Polynomial_1`.
*/ 
typedef unspecified_type Polynomial_1; 

/*!
A model of `::AlgebraicKernel_d_1::Algebraic_real_1`.
*/ 
typedef unspecified_type Algebraic_real_1; 

/*!
The choice of `Coeff` also determines the provided bound type.
In case of `Coeff` is: 

- `Gmpz` or `Gmpq`, this is `Gmpq`,

- `CORE::BigInt` or `CORE::BigRat`, this is `CORE::BigRat`,

- `leda_integer` or `leda_rational`, this is `leda_rational`.
*/ 
typedef unspecified_type Bound; 

/*!
The multiplicity type is `int`. 
*/ 
typedef unspecified_type Multiplicity_type; 

/// @}

}; /* end Algebraic_kernel_d_1 */
} /* end namespace CGAL */
