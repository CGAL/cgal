/*!
\ingroup PkgAlgebraicFoundationsRealNumberTypesConcepts
\cgalConcept

The concept `FieldNumberType` combines the requirements of the concepts 
`Field` and `RealEmbeddable`. 
A model of `FieldNumberType` can be used as a template parameter 
for Cartesian kernels. 

\cgalRefines `Field` 
\cgalRefines `RealEmbeddable` 

\cgalHasModel float 
\cgalHasModel double 
\cgalHasModel `CGAL::Gmpq` 
\cgalHasModel `CGAL::Interval_nt` 
\cgalHasModel \ref CGAL::Interval_nt_advanced
\cgalHasModel `CGAL::Lazy_exact_nt<FieldNumberType>` 
\cgalHasModel `CGAL::Quotient<RingNumberType>` 
\cgalHasModel `CGAL::leda_rational` 
\cgalHasModel `CGAL::leda_bigfloat` 
\cgalHasModel `CGAL::leda_real` 

\sa `RingNumberType` 
\sa `Kernel` 

*/

class FieldNumberType {
public:

/// @}

}; /* end FieldNumberType */

