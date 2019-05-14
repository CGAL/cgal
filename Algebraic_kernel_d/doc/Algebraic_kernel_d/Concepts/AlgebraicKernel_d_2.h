
/*!
\ingroup PkgAlgebraicKerneldConceptsBi
\cgalConcept

A model of the `AlgebraicKernel_d_2` concept gathers necessary tools 
for solving and handling bivariate polynomial systems of general degree \f$ d\f$. 

\cgalRefines `AlgebraicKernel_d_1` 
\cgalRefines `CopyConstructible` 
\cgalRefines `Assignable` 

\sa `AlgebraicKernel_d_1`

*/

class AlgebraicKernel_d_2 {
public:

/// \name Types 
/// @{

/*!

A bivariate polynomial that is a model of `Polynomial_d`, 
where \link PolynomialTraits_d::Innermost_coefficient_type `CGAL::Polynomial_traits_d<Polynomial_2>::Innermost_coefficient_type` \endlink
is `AlgebraicKernel_d_1::Coefficient`. 

*/ 
typedef unspecified_type Polynomial_2; 

/*!

A type that is used to represent real solutions of bivariate zero dimensional polynomial systems. 
A model of `DefaultConstructible`, `CopyConstructible` and 
`Assignable`. 

*/ 
typedef unspecified_type Algebraic_real_2; 

/// @} 

/// \name Functors 
/// @{

/*!
A model of `AlgebraicKernel_d_2::ConstructAlgebraicReal_2`. 
*/ 
typedef unspecified_type Construct_algebraic_real_2; 

/*!
A model of `AlgebraicKernel_d_2::ComputePolynomialX_2`. 
*/ 
typedef unspecified_type Compute_polynomial_x_2; 

/*!
A model of `AlgebraicKernel_d_2::ComputePolynomialY_2`. 
*/ 
typedef unspecified_type Compute_polynomial_y_2; 

/*!
A model of `AlgebraicKernel_d_2::Isolate_2`. 
*/ 
typedef unspecified_type Isolate_2; 

/*!
A model of `AlgebraicKernel_d_2::IsolateX_2`. 
*/ 
typedef unspecified_type Isolate_x_2; 

/*!
A model of `AlgebraicKernel_d_2::IsolateY_2`. 
*/ 
typedef unspecified_type Isolate_y_2; 

/*!
A model of `AlgebraicKernel_d_2::IsSquareFree_2`. 
*/ 
typedef unspecified_type Is_square_free_2; 

/*!
A model of `AlgebraicKernel_d_2::MakeSquareFree_2`. 
*/ 
typedef unspecified_type Make_square_free_2; 

/*!
A model of `AlgebraicKernel_d_2::SquareFreeFactorize_2`. 
*/ 
typedef unspecified_type Square_free_factorize_2; 

/*!
A model of `AlgebraicKernel_d_2::IsCoprime_2`. 
*/ 
typedef unspecified_type Is_coprime_2; 

/*!
A model of `AlgebraicKernel_d_2::MakeCoprime_2`. 
*/ 
typedef unspecified_type Make_coprime_2; 

/*!
A model of `AlgebraicKernel_d_2::Solve_2`. 
*/ 
typedef unspecified_type Solve_2; 

/*!
A model of `AlgebraicKernel_d_2::NumberOfSolutions_2`. 
*/ 
typedef unspecified_type Number_of_solutions_2; 

/*!
A model of `AlgebraicKernel_d_2::SignAt_2`. 
*/ 
typedef unspecified_type Sign_at_2; 

/*!
A model of `AlgebraicKernel_d_2::CompareX_2`. 
*/ 
typedef unspecified_type Compare_x_2; 

/*!
A model of `AlgebraicKernel_d_2::CompareY_2`. 
*/ 
typedef unspecified_type Compare_y_2; 

/*!
A model of `AlgebraicKernel_d_2::CompareXY_2`. 
*/ 
typedef unspecified_type Compare_xy_2; 

/*!
A model of `AlgebraicKernel_d_2::BoundBetweenX_2`. 
*/ 
typedef unspecified_type Bound_between_x_2; 

/*!
A model of `AlgebraicKernel_d_2::BoundBetweenY_2`. 
*/ 
typedef unspecified_type Bound_between_y_2; 

/*!
A model of `AlgebraicKernel_d_2::ApproximateAbsoluteX_2`. 
*/ 
typedef unspecified_type Approximate_absolute_x_2; 

/*!
A model of `AlgebraicKernel_d_2::ApproximateAbsoluteY_2`. 
*/ 
typedef unspecified_type Approximate_absolute_y_2; 

/*!
A model of `AlgebraicKernel_d_2::ApproximateRelativeX_2`. 
*/ 
typedef unspecified_type Approximate_relative_x_2; 

/*!
A model of `AlgebraicKernel_d_2::ApproximateRelativeY_2`. 
*/ 
typedef unspecified_type Approximate_relative_y_2; 

/// @} 

/// \name Operations 
/// For each of the function objects above, there must exist a member
/// function that requires no arguments and returns an instance of
/// that function object. The name of the member function is the
/// uncapitalized name of the type returned with the suffix `_object`
/// appended. For example, for the function object
/// `AlgebraicKernel_d_2::Bound_betweenX_2` the following member
/// function must exist:
/// @{

/*!

*/ 
AlgebraicKernel_d_2::Bound_between_x_2 bound_between_x_2_object() const; 

/// @}

}; /* end AlgebraicKernel_d_2 */

