/*!
\ingroup PkgHDVFConcepts
\cgalConcept

The concept `GeometricChainComplex` refines the concept `AbstractChainComplex` and describes the requirements for (topological) chain complexes associated to geometric complexes used in the concept `HDVF`.
 It adds to `AbstractChainComplex` methods to get vertex coordinates.

 \cgalRefines{AbstractChainComplex}

\cgalHasModelsBegin
\cgalHasModelsBare{`CGAL::Homological_discrete_vector_field::Simplicial_chain_complex<CoefficientRing,Traits>`}
\cgalHasModelsBare{`CGAL::Homological_discrete_vector_field::Cubical_chain_complex<CoefficientRing,Traits>`}
\cgalHasModelsEnd

*/

class GeometricChainComplex
{
public:
/// \name Types
/// @{

/*!
 Type of coordinates (the vector size is free, hence coordinates can be any dimension).
 */
typedef unspecified_type Point ;
/// @}


/// \name Access functions
/// @{

/*!
Returns the vector of vertex coordinates.
 */
const std::vector<Point>& points() const;

/*!
 Returns the coordinates of the cell of index `i` and dimension 0.
 */
Point point(size_t i) const;

/// @}

};
