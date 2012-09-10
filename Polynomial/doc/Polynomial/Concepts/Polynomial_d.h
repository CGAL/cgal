
/*!
\ingroup PkgPolynomialConcepts
\cgalconcept

A model of `Polynomial_d` is representing a multivariate 
polynomial in \f$ d \geq 1\f$ variables over some basic ring \f$ R\f$. 
This type is denoted as the innermost coefficient. 
A model of `Polynomial_d` accompanied by a traits class 
`CGAL::Polynomial_traits_d<Polynomial_d>`, which is a model of 
`PolynomialTraits_d`. 
Please have a look at the concept `PolynomialTraits_d`, since nearly 
all functionality related to polynomials is provided by the traits. 

\refines ::IntegralDomainWithoutDivision 
The algebraic structure of \refines ::Polynomial_d depends on the 
algebraic structure of \refines ::Innermost_coefficient_type: 

<TABLE CELLSPACING=5 > 
<TR><TD ALIGN=LEFT NOWRAP COLSPAN=2><HR> 
<TR> 
<TD ALIGN=LEFT NOWRAP> 
\refines ::Innermost_coefficient_type 
<TD ALIGN=LEFT NOWRAP> 
\refines ::Polynomial_d 
<TR><TD ALIGN=LEFT NOWRAP COLSPAN=2><HR> 
<TR> 
<TD ALIGN=LEFT NOWRAP> 
\refines ::IntegralDomainWithoutDivision 
<TD ALIGN=LEFT NOWRAP> 
\refines ::IntegralDomainWithoutDivision 
<TR> 
<TD ALIGN=LEFT NOWRAP> 
\refines ::IntegralDomain 
<TD ALIGN=LEFT NOWRAP> 
\refines ::IntegralDomain 
<TR> 
<TD ALIGN=LEFT NOWRAP> 
\refines ::UniqueFactorizationDomain 
<TD ALIGN=LEFT NOWRAP> 
\refines ::UniqueFactorizationDomain 
<TR> 
<TD ALIGN=LEFT NOWRAP> 
\refines ::EuclideanRing 
<TD ALIGN=LEFT NOWRAP> 
\refines ::UniqueFactorizationDomain 
<TR> 
<TD ALIGN=LEFT NOWRAP> 
\refines ::Field 
<TD ALIGN=LEFT NOWRAP> 
\refines ::UniqueFactorizationDomain 
<TR><TD ALIGN=LEFT NOWRAP COLSPAN=2><HR> 
</TABLE> 

\note In case the polynomial is univariate and the innermost
coefficient is a Field the polynomial is model of EuclideanRing.

\sa \ref ::AlgebraicStructureTraits 
\sa \ref ::PolynomialTraits_d 

\hasModel CGAL::Polynomial<Coeff>

*/

class Polynomial_d {
public:

/// @}

}; /* end Polynomial_d */

