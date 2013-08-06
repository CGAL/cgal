
/*!
\ingroup PkgPolynomialConcepts
\cgalConcept

This `AdaptableBinaryFunction` computes the <I>resultant</I> of two polynomials 
\f$ f\f$ and \f$ g\f$ of type `PolynomialTraits_d::Polynomial_d` with respect to a 
certain variable. 

Note that this functor operates on the polynomial in the univariate view, 
that is, the polynomial is considered as a univariate polynomial in one 
specific variable. 

Let \f$ f\f$ and \f$ g\f$ be two univariate polynomials over some commutative ring \f$ A\f$, 
where 
\f[ f = f_mx^m + \dots + f_0 \f] and 
\f[ g = g_nx^n + \dots + g_0. \f] 
The resultant of \f$ f\f$ and \f$ g\f$ is defined as the determinant of the <I>Sylvester matrix</I>: 

\image html sylvester_matrix.png
\image latex sylvester_matrix.png

Note that this is a \f$ (n+m)\times(n+m)\f$ matrix as there are \f$ n\f$ rows for \f$ f\f$ 
and \f$ m\f$ rows that are used for \f$ g\f$. The blank spaces are supposed to be 
filled with zeros. 

\cgalAdvancedBegin
Let \f$ L\f$ be the algebraic closure of \f$ A\f$, and write \f$ f\f$ and \f$ g\f$ as 
\f[ f := f_m \ccProd{i=1}{m}{(x-\alpha_i)},\ \alpha_i \in L \f] 
and 
\f[ g := g_n \ccProd{j=1}{n}{(x-\beta_j)},\ \beta_i \in L, \f] then 
the resultant of \f$ f\f$ and \f$ g\f$ is (up to leading coefficients) 
the product of all pairwise differences of the roots of \f$ f\f$ and \f$ g\f$, namely 
\f[ res(f,g) = f_m^n g_n^m \ccProd{i=1}{m}{\ccProd{j=1}{n}{(\alpha_i-\beta_j)}}. \f] 
In particular, \f$ res(f,g) \neq 0\f$ iff \f$ f\f$ and \f$ g\f$ have a common factor with a 
positive degree in \f$ X\f$. 
\cgalAdvancedEnd

There are various ways to compute the resultant. 
Naive options are the computation of the resultant as the determinant of 
the Sylvester Matrix or the Bezout 
Matrix as well as the so called subresultant algorithm, 
which is a variant of the Euclidean Algorithm. 
More sophisticated methods may use modular arithmetic and interpolation. 
For more information we refer to, e.g., \cgalCite{gg-mca-99}. 

\cgalRefines `AdaptableBinaryFunction` 
\cgalRefines `CopyConstructible` 
\cgalRefines `DefaultConstructible` 

\sa `Polynomial_d`
\sa `PolynomialTraits_d`
\sa `PolynomialTraits_d::UnivariateContent`
\sa `PolynomialTraits_d::PolynomialSubresultants`
\sa `PolynomialTraits_d::PrincipalSubresultants`

*/

class PolynomialTraits_d::Resultant {
public:

/// \name Types 
/// @{

/*!

*/ 
typedef PolynomialTraits_d::Coefficient_type result_type; 

/*!

*/ 
typedef PolynomialTraits_d::Polynomial_d first_argument_type; 

/*!

*/ 
typedef PolynomialTraits_d::Polynomial_d second_argument_type; 

/// @} 

/// \name Operations 
/// @{

/*!
Computes the resultant of \f$ f\f$ and \f$ g\f$, 
with respect to the outermost variable. 
*/ 
result_type operator()(first_argument_type f, 
second_argument_type g); 

/// @}

}; /* end PolynomialTraits_d::Resultant */

