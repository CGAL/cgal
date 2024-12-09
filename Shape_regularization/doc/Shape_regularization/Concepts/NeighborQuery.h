namespace CGAL {
namespace Shape_regularization {

/*!
\ingroup PkgShapeRegularizationRefConcepts
\cgalConcept

A concept that describes the set of methods used by the class
`QP_regularization` to access neighbors of
a geometric object being regularized.

\cgalHasModelsBegin
\cgalHasModels{Segments::Delaunay_neighbor_query_2}
\cgalHasModelsEnd
*/
class NeighborQuery {

public:

  /*!
    fills in `neighbors` with indices of all geometric objects, which are
    direct neighbors of the object with the index `query_index`.

    `QP_regularization` calls this method
    once for each object from the input range.
  */
  void operator()(
    const std::size_t query_index,
    std::vector<std::size_t>& neighbors) {

  }
};

} // namespace Shape_regularization
} // namespace CGAL
