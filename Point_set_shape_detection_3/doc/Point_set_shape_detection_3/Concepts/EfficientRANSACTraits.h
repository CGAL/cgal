/*!
\ingroup PkgPointSetShapeDetection3Concepts
\cgalConcept

Concept describing the set of types required by the class `CGAL::Shape_detection_3::Efficient_RANSAC` and all shape classes.

To avoid copying potentially large input data, the shape detection class
`CGAL::Shape_detection_3::Efficient_RANSAC` will work on the input
data directly and no internal copy will be made. For this reason, the
input data has to be provided in form of a random access iterator.
Point and normal property maps have to
be provided to extract the points and the normals from the input.

\cgalHasModel `CGAL::Shape_detection_3::Efficient_RANSAC_traits`

*/
class EfficientRANSACTraits{
public:

  /// the point type
  typedef unspecified_type Point_3;

  /// the vector type
  typedef unspecified_type Vector_3;

  /// Model of the concept `Range` with random access iterators, providing input points and normals
  /// through the two property maps `Point_map` and `Normal_map`.
  typedef unspecified_type Input_range;
  /// a model of `ReadablePropertyMap` with `std::iterator_traits<Input_range::iterator>::%value_type` as key type and `Point_3` as value type.
  typedef unspecified_type Point_map;
  /// a model of `ReadablePropertyMap` with `std::iterator_traits<Input_range::iterator>::%value_type` as key type and `Vector_3` as value type.
  typedef unspecified_type Normal_map;
};
