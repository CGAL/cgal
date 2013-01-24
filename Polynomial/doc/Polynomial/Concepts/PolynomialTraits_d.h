
/*!
\ingroup PkgPolynomialConcepts
\cgalConcept

A model of `PolynomialTraits_d` is associated with a type 
`Polynomial_d`. 
The type `Polynomial_d` represents a multivariate polynomial. 
The number of variables is denoted as the dimension \f$ d\f$ of the polynomial, 
it is arbitrary but fixed for a certain model of this concept. 
Note that univariate polynomials are not excluded by this concept. In this case 
\f$ d\f$ is just set to one. 

`PolynomialTraits_d` provides two different views on the 
multivariate polynomial. 

<UL> 
<LI>The recursive view: In this view, the polynomial is considered as 
an element of \f$ R[x_0,\dots,x_{d-2}][x_{d-1}]\f$. That is, the polynomial 
is treated as a univariate polynomial over the ring \f$ R[x_0,\dots,x_{d-2}]\f$. 
<LI>The symmetric or multivariate view: This view is symmetric 
with respect to all variables, 
considering the polynomial as an element of \f$ R [x_0,\dots,x_{d-1}]\f$. 
</UL> 

Many functors consider the polynomial as a univariate polynomial in one variable. 
By default this is the outermost variable \f$ x_{d-1}\f$. However, in general it 
is possible to select a certain variable. 

\cgalRefines `AlgebraicStructureTraits` 

\sa `Polynomial_d`

\cgalHasModel `CGAL::Polynomial_traits_d<Polynomial_d>`

*/

class PolynomialTraits_d {
public:

/// \name Constants 
/// @{

/*! 
The dimension and the number of variables respectively. 
*/ 
static const int d; 

/// @} 

/// \name Types 
/// @{

/*! 
Type representing \f$ R[x_0,\dots,x_{d-1}]\f$. 
*/ 
typedef Hidden_type Polynomial_d; 

/*! 
Type representing \f$ R[x_0,\dots,x_{d-2}]\f$. 
*/ 
typedef Hidden_type Coefficient_type ; 

/*! 
Type representing the base ring \f$ R\f$. 
*/ 
typedef Hidden_type Innermost_coefficient_type; 

/*! 
Const iterator used to iterate through all coefficients of the polynomial. 
*/ 
typedef Hidden_type Coefficient_const_iterator; 

/*! 
Const iterator used to iterate through all innermost coefficients of the polynomial. 
*/ 
typedef Hidden_type Innermost_coefficient_const_iterator; 

/*! 
This template class has to define a type `Rebind<T,d>::%Other` which is a model
of the concept `PolynomialTraits_d`, where `d` is the number of 
variables and `T` the `Innermost_coefficient_type`. 
*/ 
template <typename T, int d>
using Rebind = Hidden_type;

/// @} 

/// \name Functors 
/// In case a functor is not provided it is set to `CGAL::Null_functor`.
/// @{

/*! 
A model of `PolynomialTraits_d::ConstructPolynomial`. 
*/ 
typedef Hidden_type Construct_polynomial; 

/*! 
A model of `PolynomialTraits_d::GetCoefficient`. 
*/ 
typedef Hidden_type Get_coefficient; 

/*! 
A model of `PolynomialTraits_d::GetInnermostCoefficient`. 
*/ 
typedef Hidden_type Get_innermost_coefficient; 

/*! 
A model of 
`PolynomialTraits_d::ConstructCoefficientConstIteratorRange`. 
*/ 
typedef Hidden_type Construct_coefficient_const_iterator_range; 

/*! 
A model of 
`PolynomialTraits_d::ConstructInnermostCoefficientConstIteratorRange`. 
*/ 
typedef Hidden_type Construct_innermost_coefficient_const_iterator_range; 

/*! 
A model of `PolynomialTraits_d::Swap`. 
*/ 
typedef Hidden_type Swap; 

/*! 
A model of `PolynomialTraits_d::Move`. 
*/ 
typedef Hidden_type Move; 

/*! 
A model of `PolynomialTraits_d::Degree`. 
*/ 
typedef Hidden_type Degree; 

/*! 
A model of `PolynomialTraits_d::TotalDegree`. 
*/ 
typedef Hidden_type Total_degree; 

/*! 
A model of `PolynomialTraits_d::DegreeVector`. 
*/ 
typedef Hidden_type Degree_vector; 

/*! 
A model of `PolynomialTraits_d::LeadingCoefficient`. 
*/ 
typedef Hidden_type Leading_coefficient; 

/*! 
A model of `PolynomialTraits_d::InnermostLeadingCoefficient`. 
*/ 
typedef Hidden_type Innermost_leading_coefficient; 

/*! 
A model of `PolynomialTraits_d::Canonicalize`. 
*/ 
typedef Hidden_type Canonicalize; 

/*! 
A model of `PolynomialTraits_d::Differentiate`. 
*/ 
typedef Hidden_type Differentiate; 

/*! 
A model of `PolynomialTraits_d::Evaluate`. 
*/ 
typedef Hidden_type Evaluate; 

/*! 
A model of `PolynomialTraits_d::EvaluateHomogeneous`. 
*/ 
typedef Hidden_type Evaluate_homogeneous; 

/*! 
A model of `PolynomialTraits_d::Substitute`. 
*/ 
typedef Hidden_type Substitute; 

/*! 
A model of `PolynomialTraits_d::SubstituteHomogeneous`. 
*/ 
typedef Hidden_type Substitute_homogeneous; 

/*! 
A model of `PolynomialTraits_d::IsZeroAt`. 
*/ 
typedef Hidden_type Is_zero_at; 

/*! 
A model of `PolynomialTraits_d::IsZeroAtHomogeneous`. 
*/ 
typedef Hidden_type Is_zero_at_homogeneous; 

/*! 

A model of `PolynomialTraits_d::SignAt`. 

In case `Innermost_coefficient_type` is not `RealEmbeddable` this 
is `CGAL::Null_functor`. 
*/ 
typedef Hidden_type Sign_at; 

/*! 

A model of `PolynomialTraits_d::SignAtHomogeneous`. 

In case `Innermost_coefficient_type` is not `RealEmbeddable` this 
is `CGAL::Null_functor`. 
*/ 
typedef Hidden_type Sign_at_homogeneous; 

/*! 

A model of `PolynomialTraits_d::Compare`. 

In case `Innermost_coefficient_type` is not `RealEmbeddable` this 
is `CGAL::Null_functor`. 
*/ 
typedef Hidden_type Compare; 

/*! 

In case `PolynomialTraits_d::Coefficient_type` is <B>not</B> a model of 
`UniqueFactorizationDomain`, this is `CGAL::Null_functor`, otherwise this is 
a model of `PolynomialTraits_d::UnivariateContent`. 
*/ 
typedef Hidden_type Univariate_content; 

/*! 

In case `PolynomialTraits_d::Innermost_coefficient_type` is <B>not</B> 
a model of `UniqueFactorizationDomain`, this is `CGAL::Null_functor`, 
otherwise this is a model of 
`PolynomialTraits_d::MultivariateContent`. 
*/ 
typedef Hidden_type Multivariate_content; 

/*! 
A model of `PolynomialTraits_d::Shift`. 
*/ 
typedef Hidden_type Shift; 

/*! 
A model of `PolynomialTraits_d::Negate`. 
*/ 
typedef Hidden_type Negate; 

/*! 
A model of `PolynomialTraits_d::Invert`. 
*/ 
typedef Hidden_type Invert; 

/*! 
A model of `PolynomialTraits_d::Translate`. 
*/ 
typedef Hidden_type Translate; 

/*! 
A model of `PolynomialTraits_d::TranslateHomogeneous`. 
*/ 
typedef Hidden_type Translate_homogeneous; 

/*! 
A model of `PolynomialTraits_d::Scale`. 
*/ 
typedef Hidden_type Scale; 

/*! 
A model of `PolynomialTraits_d::ScaleHomogeneous`. 
*/ 
typedef Hidden_type Scale_homogeneous; 

/*! 
A model of `PolynomialTraits_d::MakeSquareFree`. 
*/ 
typedef Hidden_type Make_square_free; 

/*! 
In case `PolynomialTraits::Polynomial_d` 
is not a model of `UniqueFactorizationDomain`, this is of type `CGAL::Null_functor`, 
otherwise this is a model of `PolynomialTraits_d::SquareFreeFactorize`. 
*/ 
typedef Hidden_type Square_free_factorize; 

/*! 
A model of `PolynomialTraits_d::PseudoDivision`. 
*/ 
typedef Hidden_type Pseudo_division ; 

/*! 
A model of `PolynomialTraits_d::PseudoDivisionRemainder`. 
*/ 
typedef Hidden_type Pseudo_division_remainder; 

/*! 
A model of `PolynomialTraits_d::PseudoDivisionQuotient`. 
*/ 
typedef Hidden_type Pseudo_division_quotient ; 

/*! 
A model of `PolynomialTraits_d::GcdUpToConstantFactor`. 
*/ 
typedef Hidden_type Gcd_up_to_constant_factor; 

/*! 
A model of `PolynomialTraits_d::IntegralDivisionUpToConstantFactor`. 
*/ 
typedef Hidden_type Integral_division_up_to_constant_factor; 

/*! 
A model of `PolynomialTraits_d::UnivariateContentUpToConstantFactor`. 
*/ 
typedef Hidden_type Content_up_to_constant_factor; 

/*! 
A model of `PolynomialTraits_d::SquareFreeFactorizeUpToConstantFactor`. 
*/ 
typedef Hidden_type Square_free_factorize_up_to_constant_factor; 

/*! 
A model of `PolynomialTraits_d::Resultant`. 
*/ 
typedef Hidden_type Resultant; 

/*! 
Either `CGAL::Null_functor` or a model of `PolynomialTraits_d::PolynomialSubresultants`. 
*/ 
typedef Hidden_type Polynomial_subresultants; 

/*! 
Either `CGAL::Null_functor` or a model of `PolynomialTraits_d::PolynomialSubresultants_with_cofactors`. 
*/ 
typedef Hidden_type Polynomial_subresultants_with_cofactors; 

/*! 
Either `CGAL::Null_functor` or a model of `PolynomialTraits_d::PrincipalSubresultants`. 
*/ 
typedef Hidden_type Principal_subresultants; 

/*! 
Either `CGAL::Null_functor` or a model of `PolynomialTraits_d::SturmHabichtSequence`. 
*/ 
typedef Hidden_type Sturm_habicht_sequence; 

/*! 
Either `CGAL::Null_functor` or a model of `PolynomialTraits_d::SturmHabichtSequenceWithCofactors`. 
*/ 
typedef Hidden_type Sturm_habicht_sequence_with_cofactors; 

/*! 
Either `CGAL::Null_functor` or a model of `PolynomialTraits_d::PrincipalSturmHabichtSequence`. 
*/ 
typedef Hidden_type Principal_sturm_habicht_sequence; 

/// @}

}; /* end PolynomialTraits_d */

