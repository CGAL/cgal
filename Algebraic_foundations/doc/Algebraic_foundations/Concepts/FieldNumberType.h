/*!
\ingroup PkgAlgebraicFoundationsAlgebraicStructuresConcepts
\cgalConcept

The concept `FieldNumberType` combines the requirements of the concepts
`Field` and `RealEmbeddable`.
A model of `FieldNumberType` can be used as a template parameter
for %Cartesian kernels.

\cgalRefines{Field,RealEmbeddable}

\cgalHasModelsBegin
\cgalHasModels{float}
\cgalHasModels{double}
\cgalHasModels{CGAL::Gmpq}
\cgalHasModels{CGAL::Interval_nt}
\cgalHasModels{CGAL::Interval_nt_advanced}
\cgalHasModels{CGAL::Lazy_exact_nt<FieldNumberType>}
\cgalHasModels{CGAL::Quotient<RingNumberType>}
\cgalHasModels{leda_rational}
\cgalHasModels{leda_bigfloat}
\cgalHasModels{leda_real}
\cgalHasModelsEnd

\sa `RingNumberType`
\sa `Kernel`

*/

class FieldNumberType {
public:

/// @{
/// @}

}; /* end FieldNumberType */
