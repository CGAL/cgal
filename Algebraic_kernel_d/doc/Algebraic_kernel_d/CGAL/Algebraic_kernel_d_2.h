
namespace CGAL {

/*!
\ingroup PkgAlgebraicKerneldModels

This class gathers necessary tools for solving and handling bivariate 
polynomial systems of general degree \f$ d\f$.

This class is based on an algorithm computing a 
geometric-topological analysis of a single curve \cite ekw-fast-07 and of a 
pair of curves \cite ek-exact-08. 
The main idea behind both analyses is to compute the critical 
x-coordinates of curves and curve pairs by projection (resultants), and compute 
additional information about the critical fibers using subresultants 
and Sturm-Habicht sequences \cite grlr-sturm-habicht-98. 
With that information, the fiber at 
critical x-coordinates is computed by a variant of the Bitstream 
Descartes method. 
See also \cite kerber-phd-09 for a comprehensive description of 
these techniques. 

A point \f$ p\f$ of type `Algebraic_real_2` is represented 
by its \f$ x\f$-coordinate \f$ x_0\f$ (as described in the `Algebraic_kernel_d_1` 
paragraph above), an algebraic curve where \f$ p\f$ lies on, and an 
integer \f$ i\f$, denoting that \f$ p\f$ is the \f$ i\f$th point in the fiber at \f$ x_0\f$, 
counted from the bottom (ignoring a possible vertical line at \f$ x_0\f$). 
This determines the point uniquely, but the \f$ y\f$-coordinate 
is not stored internally in terms of an `Algebraic_real_1` object. 
Querying such a representation by calling `Compute_y_2` is a 
time-consuming step, and should be avoided for efficiency reasons if possible. 
Note that this representation is not exposed in the interface. 

The template argument `Coeff` determines the coefficient type of the 
kernel, which is also the innermost coefficient type of the supported polynomials. 

Currently, the following coefficient types are supported: 

- `Gmpz`, `Gmpq`, (requires configuration with external libraries GMP, MPFR and MPFI) 
- `CORE::BigInt`, `CORE::BigRat`, (requires configuration with external library GMP) 
- `leda_integer`, `leda_rational`. (requires configuration with external library LEDA) 

\cgalAdvancedBegin
The template argument type can also be set to
`Sqrt_extension<NT,ROOT>`, where `NT` is one of the types listed
above. `ROOT` should be one of the integer types.  See also the
documentation of `Sqrt_extension<NT,ROOT>`.
\cgalAdvancedEnd

\cgalModels `AlgebraicKernel_d_2`

\sa `AlgebraicKernel_d_1` 
\sa `AlgebraicKernel_d_2` 
\sa `Polynomial_d` 
\sa `CGAL::Algebraic_kernel_d_2<Coeff>`

*/
template< typename Coeff >
class Algebraic_kernel_d_2 {
public:

/// \name Types 
/// @{

/*! 
Same type as the template argument `Coeff`. 
*/ 
typedef unspecified_type Coefficient; 

/*! 
A model of `AlgebraicKernel_d_2::Polynomial_2` 
*/ 
typedef unspecified_type Polynomial_2; 

/*! 
A model of `AlgebraicKernel_d_2::AlgebraicReal_2` 
*/ 
typedef unspecified_type Algebraic_real_2; 

/*! 
The choice of `Coeff` also determines the provided bound type.
In case of `Coeff` is 
- `Gmpz` or `Gmpq`, this is `Gmpq`

- `CORE::BigInt` or `CORE::BigRat`, this is `CORE::BigRat`

- `leda_integer` or `leda_rational`, this is `leda_rational`

*/ 
typedef unspecified_type Bound; 

/*! 
The multiplicity type is `int`. 
*/ 
typedef unspecified_type Multiplicity_type; 

/// @}

}; /* end Algebraic_kernel_d_2 */
} /* end namespace CGAL */
