/*!
\ingroup PkgAlgebraicFoundationsAlgebraicStructuresConcepts
\cgalConcept

The concept `FieldNumberType` combines the requirements of the concepts
`Field` and `RealEmbeddable`.
A model of `FieldNumberType` can be used as a template parameter
for Cartesian kernels.

\cgalRefines{Field,RealEmbeddable}

\cgalHasModelsBegin
\cgalModels{float}
\cgalModels{double}
\cgalModels{CGAL::Gmpq}
\cgalModels{CGAL::Interval_nt}
\cgalModels{CGAL::Interval_nt_advanced}
\cgalModels{CGAL::Lazy_exact_nt<FieldNumberType>}
\cgalModels{CGAL::Quotient<RingNumberType>}
\cgalModels{leda_rational}
\cgalModels{leda_bigfloat}
\cgalModels{leda_real}
\cgalHasModelsEnd

\sa `RingNumberType`
\sa `Kernel`

*/

class FieldNumberType {
public:

/// @{
/// @}

}; /* end FieldNumberType */

