
/*!
\ingroup PkgPolynomialConcepts
\cgalConcept

A model of `Polynomial_d` is representing a multivariate
polynomial in \f$ d \geq 1\f$ variables over some basic ring \f$ R\f$.
This type is denoted as the innermost coefficient.
A model of `Polynomial_d` must be accompanied by a traits class
`CGAL::Polynomial_traits_d<Polynomial_d>`, which is a model of
`PolynomialTraits_d`.
Please have a look at the concept `PolynomialTraits_d`, since nearly
all functionality related to polynomials is provided by the traits.

\cgalRefines `IntegralDomainWithoutDivision`

The algebraic structure of `Polynomial_d` depends on the
algebraic structure of `PolynomialTraits_d::Innermost_coefficient_type`:

Innermost_coefficient_type       | %Polynomial_d
---------------------------------|--------------------------------
::IntegralDomainWithoutDivision  | ::IntegralDomainWithoutDivision
::IntegralDomain                 | ::IntegralDomain
::UniqueFactorizationDomain      | ::UniqueFactorizationDomain
::EuclideanRing                  | ::UniqueFactorizationDomain
::Field                          | ::UniqueFactorizationDomain


\note In case the polynomial is univariate and the innermost
coefficient is a `Field` the polynomial is model of `EuclideanRing`.

\sa `AlgebraicStructureTraits`
\sa `PolynomialTraits_d`

\cgalHasModel `CGAL::Polynomial<Coeff>`

*/

class Polynomial_d {

}; /* end Polynomial_d */

