
/*!
\ingroup PkgAlgebraicFoundationsAlgebraicStructuresConcepts
\cgalConcept

The concept `RingNumberType` combines the requirements of the concepts 
`IntegralDomainWithoutDivision` and `RealEmbeddable`. 
A model of `RingNumberType` can be used as a template parameter 
for Homogeneous kernels. 

\cgalRefines `IntegralDomainWithoutDivision` 
\cgalRefines `RealEmbeddable` 

\cgalHasModel \cpp built-in number types 
\cgalHasModel `CGAL::Gmpq` 
\cgalHasModel `CGAL::Gmpz` 
\cgalHasModel` CGAL::Interval_nt` 
\cgalHasModel \ref CGAL::Interval_nt_advanced
\cgalHasModel `CGAL::Lazy_exact_nt<RingNumberType>` 
\cgalHasModel `CGAL::MP_Float` 
\cgalHasModel `CGAL::Gmpzf`
\cgalHasModel `CGAL::Quotient<RingNumberType>` 
\cgalHasModel `leda_integer` 
\cgalHasModel `leda_rational` 
\cgalHasModel `leda_bigfloat` 
\cgalHasModel `leda_real` 

\sa `FieldNumberType` 

*/

class RingNumberType {
public:

/// @}

}; /* end RingNumberType */

