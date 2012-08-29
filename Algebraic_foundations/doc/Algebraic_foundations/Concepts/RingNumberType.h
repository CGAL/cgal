
/*!
\ingroup PkgAlgebraicFoundationsRealNumberTypesConcepts
\cgalconcept

The concept `RingNumberType` combines the requirements of the concepts 
`IntegralDomainWithoutDivision` and `RealEmbeddable`. 
A model of `RingNumberType` can be used as a template parameter 
for Homogeneous kernels. 

\refines ::IntegralDomainWithoutDivision 
\refines ::RealEmbeddable 

\hasModel \cpp built-in number types 
\hasModel CGAL::Gmpq 
\hasModel CGAL::Gmpz 
\hasModel CGAL::Interval_nt 
\hasModel CGAL::Interval_nt_advanced 
\hasModel CGAL::Lazy_exact_nt<RingNumberType> 
\hasModel CGAL::MP_Float 
\hasModel CGAL::Gmpzf 
\hasModel CGAL::Quotient<RingNumberType> 
\hasModel leda_integer 
\hasModel leda_rational 
\hasModel leda_bigfloat 
\hasModel leda_real 

\sa `FieldNumberType` 

*/

class RingNumberType {
public:

/// @}

}; /* end RingNumberType */

