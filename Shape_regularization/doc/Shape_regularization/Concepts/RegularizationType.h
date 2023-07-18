namespace CGAL {
namespace Shape_regularization {

/*!
\ingroup PkgShapeRegularizationRefConcepts
\cgalConcept

A concept that describes the set of methods used by the class
`QP_regularization` to access various data
required for setting up the the global regularization problem.

\cgalHasModelsBegin
\cgalHasModels{Segments::Angle_regularization_2}
\cgalHasModels{Segments::Offset_regularization_2}
\cgalHasModelsEnd
*/
class RegularizationType {

public:

  /*!
    returns the maximum bound on a regularization characteristic (angle-orientation/
    distance-offset/etc.) with respect to which a geometric object with the index
    `query_index` is being regularized.

    `QP_regularization` calls this method
    once for each object from the input range.
  */
  FieldNumberType bound(
    const std::size_t query_index) const {

  }

  /*!
    returns an objective function value between two geometric objects, which are
    direct neighbors, that is they form a neighbor pair `i <-> j`. Neighbors are
    provided by the concept `NeighborQuery`.

    `QP_regularization` calls this method
    once for each neighbor pair being regularized.
  */
  FieldNumberType target(
    const std::size_t i, const std::size_t j) {

  }

  /*!
    updates regularization characteristics (angle-orientation/distance-offset/etc.)
    of the geometric objects being regularized using values from `solution`, one
    value per one regularized object. These values depend on what is being regularized,
    they could be angle or offset differences for example. The solution vector is
    computed by the `QuadraticProgramTraits`.

    Number of values in `solution` equals to the number n of geometric objects being
    regularized + the number m of neighbor pairs between these objects. The first
    n values are the values that should be used.

    `QP_regularization` calls this method
    once after the global regularization QP problem has been solved.
  */
  void update(
    const std::vector<FieldNumberType>& solution) {

  }
};

} // namespace Shape_regularization
} // namespace CGAL
