namespace CGAL {
namespace Weights {

/*!
\ingroup PkgWeightsRefConcepts
\cgalConcept

A concept that describes the set of methods required in all classes for computing
2D generalized barycentric weights with respect to polygons.

\cgalHasModel
- `CGAL::Weights::Wachspress_weights_2`
- `CGAL::Weights::Mean_value_weights_2`
- `CGAL::Weights::Discrete_harmonic_weights_2`
*/
class BarycentricWeights_2 {

public:

  /*!
    fills a destination range with 2D generalized barycentric weights
    computed at the `query` point with respect to the vertices of the input polygon.

    The number of computed weights equals to the number of polygon vertices.
  */
  template<typename OutputIterator>
  OutputIterator operator()(
    const Point_2& query, OutputIterator w_begin)
  { }
};

} // namespace Weights
} // namespace CGAL
