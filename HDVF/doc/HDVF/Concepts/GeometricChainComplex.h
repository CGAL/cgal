/*!
\ingroup PkgHDVFConcepts
\cgalConcept

The concept `GeometricChainComplex` refines the concept `AbstractChainComplex` and describes the requirements for (topological) chain complexes associated to geometric complexes used in the concept `HDVF`.
 It adds to `AbstractChainComplex` methods to get vertex coordinates.

 \cgalRefines{AbstractChainComplex}

\cgalHasModelsBegin
\cgalHasModelsBare{`CGAL::HDVF::Simplicial_chain_complex<CoefficientType>`}
\cgalHasModelsBare{`CGAL::HDVF::Cubical_chain_complex<CoefficientType>`}
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
typedef std::vector<double> Point ;
/// @}


/// \name Access functions
/// @{

/*!
Returns the vector of vertex coordinates.
 */
const std::vector<Point>& get_vertices_coords() const;

/*!
 Returns the coordinates of the cell of index `i` and dimension 0.
 */
Point get_vertex_coords (size_t i) const;

/// @}

};
