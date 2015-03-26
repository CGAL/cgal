/*!
\ingroup PkgPointSetShapeDetection3Concepts
\cgalConcept

Concept describing the set of require type needed by the class `CGAL::Shape_detection_3::Efficient_RANSAC` and all shape classes.

To avoid copying of potentially large input data, the shape detection class
`CGAL::Efficient_RANSAC` will work on the input
data directly and no internal copy will be done. For this reason the
input data has to be provided in form of a random access iterator.
Point and normal property maps have to
be provided to extract from the input the points and the normals respectively.

\cgalHasModel `CGAL::Shape_detection_3::Efficient_RANSAC_traits`

*/
class EfficientRANSACTraits{
public:
  /// Geometric traits.
  /// It must provide `Geom_traits::FT`, `Geom_traits::Point_3` and `Geom_traits::Vector_3`.
  /// `Geom_traits::FT` must be a floating point number type like `double` or `float`.
  typedef unspecified_type Geom_traits;
  /// Random access iterator used to get the input points and normals.
  typedef InputIt Input_iterator;
  /// property map: a model of `ReadablePropertyMap` with `Input_iterator` as key type and `Geom_traits::Point_3` as value type
  typedef Ppmap Point_pmap;
  /// property map: a model of `ReadablePropertyMap` with `Input_iterator` as key type and `Geom_traits::Vector_3` as value type
  typedef Npmap Normal_pmap;
};
