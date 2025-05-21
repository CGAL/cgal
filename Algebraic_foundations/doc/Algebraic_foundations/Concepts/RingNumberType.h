
/*!
\ingroup PkgAlgebraicFoundationsAlgebraicStructuresConcepts
\cgalConcept

The concept `RingNumberType` combines the requirements of the concepts
`IntegralDomainWithoutDivision` and `RealEmbeddable`.
A model of `RingNumberType` can be used as a template parameter
for homogeneous kernels.

\cgalRefines{IntegralDomainWithoutDivision,RealEmbeddable}

\cgalHasModelsBegin
\cgalHasModelsBare{\cpp built-in number types}
\cgalHasModels{CGAL::Gmpq}
\cgalHasModels{CGAL::Gmpz}
\cgalHasModels{CGAL::Interval_nt}
\cgalHasModels{CGAL::Interval_nt_advanced}
\cgalHasModels{CGAL::Lazy_exact_nt<RingNumberType>}
\cgalHasModels{CGAL::MP_Float}
\cgalHasModels{CGAL::Gmpzf}
\cgalHasModels{CGAL::Quotient<RingNumberType>}
\cgalHasModels{leda_integer}
\cgalHasModels{leda_rational}
\cgalHasModels{leda_bigfloat}
\cgalHasModels{leda_real}
\cgalHasModelsEnd

\sa `FieldNumberType`

*/

class RingNumberType {
public:

}; /* end RingNumberType */
