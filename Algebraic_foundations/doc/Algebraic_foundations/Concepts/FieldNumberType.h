typedef int X;


/*!
\ingroup PkgAlgebraicFoundationsRealNumberTypesConcepts
\cgalconcept

The concept `FieldNumberType` combines the requirements of the concepts 
`Field` and `RealEmbeddable`. 
A model of `FieldNumberType` can be used as a template parameter 
for Cartesian kernels. 

\refines ::Field 
\refines ::RealEmbeddable 

\hasModel float 
\hasModel double 
\hasModel `CGAL::Gmpq` 
\hasModel `CGAL::Interval_nt` 
\hasModel \ref CGAL::Interval_nt_advanced
\hasModel `CGAL::Lazy_exact_nt<FieldNumberType>` 
\hasModel `CGAL::Quotient<RingNumberType>` 
\hasModel `CGAL::leda_rational` 
\hasModel `CGAL::leda_bigfloat` 
\hasModel `CGAL::leda_real` 

\sa `RingNumberType` 
\sa `Kernel` 

*/

class FieldNumberType {
public:

/// @}

}; /* end FieldNumberType */

