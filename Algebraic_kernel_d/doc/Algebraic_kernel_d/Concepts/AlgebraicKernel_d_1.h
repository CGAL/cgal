
/*!
\ingroup PkgAlgebraicKerneldConceptsUni
\cgalConcept

A model of the `AlgebraicKernel_d_1` concept is meant to provide the 
algebraic functionalities on univariate polynomials of general degree \f$ d\f$. 

\cgalRefines `CopyConstructible` 
\cgalRefines `Assignable` 

\cgalHasModel `CGAL::Algebraic_kernel_rs_gmpz_d_1`
\cgalHasModel `CGAL::Algebraic_kernel_rs_gmpq_d_1`

\sa `AlgebraicKernel_d_2`

A model of `AlgebraicKernel_d_1` must provide: 
*/

class AlgebraicKernel_d_1 {
public:

/// \name Types 
/// @{

/*!

A model of `IntegralDomain` and `RealEmbeddable`. 

`ExplicitInteroperable` with `AlgebraicKernel_d_1::Bound`. 

*/ 
typedef unspecified_type Coefficient; 

/*!

A univariate polynomial that is a model of `Polynomial_d`, 
where \link PolynomialTraits_d::Innermost_coefficient_type `CGAL::Polynomial_traits_d<Polynomial_1>::Innermost_coefficient_type` \endlink
is `AlgebraicKernel_d_1::Coefficient`. 

*/ 
typedef unspecified_type Polynomial_1; 

/*!

A type that is used to represent real roots of univariate polynomials. 
The type must be a model of `DefaultConstructible`, 
`CopyConstructible`, `Assignable` and `RealEmbeddable`. 

*/ 
typedef unspecified_type Algebraic_real_1; 

/*!

A type to represent upper and lower bounds of `AlgebraicKernel_d_1::Algebraic_real_1`. 

The type is `ExplicitInteroperable` with 
`AlgebraicKernel_d_1::Coefficient` and must be a model 
`IntegralDomain`, `RealEmbeddable` and dense in \f$ \mathbb{R}\f$. 

*/ 
typedef unspecified_type Bound; 

/*!
Size type (unsigned integral type). 
*/ 
typedef unspecified_type size_type; 

/*!
Multiplicity type (unsigned integral type). 
*/ 
typedef unspecified_type Multiplicity_type; 

/// @} 

/// \name Functors 
/// @{

/*!
A model of `AlgebraicKernel_d_1::ConstructAlgebraicReal_1`. 
*/ 
typedef unspecified_type Construct_algebraic_real_1; 

/*!
A model of `AlgebraicKernel_d_1::ComputePolynomial_1`. 
*/ 
typedef unspecified_type Compute_polynomial_1; 

/*!
A model of `AlgebraicKernel_d_1::Isolate_1`. 
*/ 
typedef unspecified_type Isolate_1; 

/*!
A model of `AlgebraicKernel_d_1::IsSquareFree_1`. 
*/ 
typedef unspecified_type Is_square_free_1; 

/*!
A model of `AlgebraicKernel_d_1::MakeSquareFree_1`. 
*/ 
typedef unspecified_type Make_square_free_1; 

/*!
A model of `AlgebraicKernel_d_1::SquareFreeFactorize_1`. 
*/ 
typedef unspecified_type Square_free_factorize_1; 

/*!
A model of `AlgebraicKernel_d_1::IsCoprime_1`. 
*/ 
typedef unspecified_type Is_coprime_1; 

/*!
A model of `AlgebraicKernel_d_1::MakeCoprime_1`. 
*/ 
typedef unspecified_type Make_coprime_1; 

/*!
A model of `AlgebraicKernel_d_1::Solve_1`. 
*/ 
typedef unspecified_type Solve_1; 

/*!
A model of `AlgebraicKernel_d_1::NumberOfSolutions_1`. 
*/ 
typedef unspecified_type Number_of_solutions_1; 

/*!
A model of `AlgebraicKernel_d_1::SignAt_1`. 
*/ 
typedef unspecified_type Sign_at_1; 

/*!
A model of `AlgebraicKernel_d_1::Compare_1`. 
*/ 
typedef unspecified_type Compare_1; 

/*!
A model of `AlgebraicKernel_d_1::BoundBetween_1`. 
*/ 
typedef unspecified_type Bound_between_1; 

/*!
A model of `AlgebraicKernel_d_1::ApproximateAbsolute_1`. 
*/ 
typedef unspecified_type Approximate_absolute_1; 

/*!
A model of `AlgebraicKernel_d_1::ApproximateRelative_1`. 
*/ 
typedef unspecified_type Approximate_relative_1; 

/// @} 

/// \name Operations 
/// For each of the function objects above, there must exist a member
/// function that requires no arguments and returns an instance of
/// that function object. The name of the member function is the
/// uncapitalized name of the type returned with the suffix `_object`
/// appended. For example, for the function object
/// `AlgebraicKernel_d_1::Bound_between_1` the following member
/// function must exist:
/// @{

/*!

*/ 
AlgebraicKernel_d_1::Bound_between_1 bound_between_1_object() const; 

/// @}

}; /* end AlgebraicKernel_d_1 */

