/*!
\ingroup PkgWeightsRefConcepts
\cgalConcept

A concept that describes the set of methods required in all classes used in
the computation of 2D generalized barycentric weights.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Weights::Wachspress_weights_2}
\cgalHasModels{CGAL::Weights::Mean_value_weights_2}
\cgalHasModels{CGAL::Weights::Discrete_harmonic_weights_2}
\cgalHasModelsEnd
*/
class BarycentricWeights_2 {

public:

  /*!
    fills a destination range with 2D generalized barycentric weights
    computed at the `query` point with respect to the vertices of the input polygon.

    \tparam OutIterator
    a model of `OutputIterator` whose value type is `FieldNumberType`

    The number of computed weights is equal to the number of polygon vertices.
  */
  template<typename OutIterator>
  OutIterator operator()(
    const Point_2& query, OutIterator w_begin)
  { }
};
