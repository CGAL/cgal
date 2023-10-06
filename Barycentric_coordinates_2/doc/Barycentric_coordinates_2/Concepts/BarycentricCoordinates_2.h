namespace CGAL {
namespace Barycentric_coordinates {

/*!
\ingroup PkgBarycentricCoordinates2RefConcepts
\cgalConcept

A concept that describes the set of methods that should be defined for all coordinate
models used to parameterize the class `Generalized_barycentric_coordinates_2`.

\cgalHasModelsBegin
\cgalHasModels{Wachspress_2}
\cgalHasModels{Mean_value_2}
\cgalHasModels{Discrete_harmonic_2}
\cgalHasModelsEnd

\deprecated This part of the package is deprecated since the version 5.4 of \cgal.
*/
class BarycentricCoordinates_2 {

public:

  /// \name Creation
  /// @{

  /*!
    Creates a class that implements generalized barycentric coordinates for any query
    point that does not belong to the polygon's boundary. The polygon is given by a
    range of vertices of the type `Traits::Point_2` stored in a container of the
    type <a href="https://en.cppreference.com/w/cpp/container/vector">`std::vector`</a>.
  */
  BarycentricCoordinates_2(
    const std::vector<Traits::Point_2>& vertices, const Traits& barycentric_traits) {

  }

  /// @}

  /// \name Functions
  /// @{

  /*!
    A function that computes generalized barycentric coordinates without normalization
    that are called generalized barycentric weights (as fast as possible algorithm is used).
    Weights are computed with respect to a query point of the type `Traits::Point_2` and
    stored in the output iterator `output`. The function returns a pointer to the last stored element.
  */
  std::optional<OutputIterator>
  weights(
    const Traits::Point_2& query_point, OutputIterator& output) {

  }

  /*!
    A function that computes generalized barycentric coordinates on the bounded side
    of a polygon with one of two possible algorithms: one is precise and one is fast.
    The algorithm type is specified by the parameter `type_of_algorithm`. Coordinates
    are computed with respect to a query point of the type `Traits::Point_2` and stored
    in the output iterator `output`. The function returns a pointer to the last stored element.
  */
  std::optional<OutputIterator>
  coordinates_on_bounded_side(
    const Traits::Point_2& query_point,
    OutputIterator& output,
    const Type_of_algorithm type_of_algorithm) {

  }

  /*!
    A function that computes generalized barycentric coordinates on the unbounded side
    of a polygon with one of two possible algorithms: one is precise and one is fast.
    The algorithm type is specified by the parameter `type_of_algorithm`. Coordinates
    are computed with respect to a query point of the type `Traits::Point_2` and stored
    in the output iterator `output`. The function returns a pointer to the last stored element.
  */
  std::optional<OutputIterator>
  coordinates_on_unbounded_side(
    const Traits::Point_2& query_point,
    OutputIterator& output,
    const Type_of_algorithm type_of_algorithm) {

  }

/// @}

};

} // namespace Barycentric_coordinates
} // namespace CGAL
