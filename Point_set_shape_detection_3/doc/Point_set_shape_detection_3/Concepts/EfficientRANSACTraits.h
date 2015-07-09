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
  
/// \name Types
/// @{

  /// The point type
  typedef unspecified_type Point_3;
  /// The vector type
  typedef unspecified_type Vector_3;
  /// The sphere type, only required if you want to detect spheres
  typedef unspecified_type Sphere_3;
  /// The line type, only required if you want to detect cylinders
  typedef unspecified_type Line_3;
  /// The plane type, only required if you want to detect planes
  typedef unspecified_type Plane_3;
  /// The 2D point type, only required if you want to detect tori
  typedef unspecified_type Point_2;
  /// The circle type, only required if you want to detect tori
  typedef unspecified_type Circle_2;

  /// The number type of the Cartesian coordinates of types Point_3
  typedef unspecified_type FT;

  /// A model of the concept `Range` with random access iterators, providing input points and normals
  /// through the two property maps `Point_map` and `Normal_map`.
  typedef unspecified_type Input_range;
  /// A model of `ReadablePropertyMap` with 
  /// `std::iterator_traits<Input_range::iterator>::%value_type` 
  /// as key type and `Point_3` as value type.
  typedef unspecified_type Point_map;
  /// A model of `ReadablePropertyMap` with 
  /// `std::iterator_traits<Input_range::iterator>::%value_type` 
  /// as key type and `Vector_3` as value type.
  typedef unspecified_type Normal_map;

  /// a model of `SearchTraits`
  /// where `SearchTraits::point_d` is `Point_3`,
  /// `SearchTraits::Dimension` is ` CGAL::Dimension_tag<3>`,
  /// and `SearchTraits::FT` is ` FT`,
  typedef unspecified_type Search_traits;
  
  /*!
   * Function object type that provides
   * `Point_3 operator()(FT x, FT y, FT z)`
   * returning the point with `x`, `y` and `z` as Cartesian coordinates.
   */
  typedef unspecified_type Construct_point_3;

  /*!
   * Function object type that provides `Point_3 operator()(Point_3 p1, Point_3 p2)`
   * returning the vector `p1p2`, `Vector_3 operator()(NULL_VECTOR)` returning 
   * the null vector, and `Vector_3 operator()(Line_3 l)` returning 
   * a vector having the same direction as `l`.
   */
  typedef unspecified_type Construct_vector_3;
  
  /*!
   * Function object type that provides `Sphere_3 operator()(Point_3 c, FT r)`
   * returning the sphere of center `p` and radius `r`.
   */
  typedef unspecified_type Construct_sphere_3;
  
  /*!
   * Function object type that provides `Line_3 operator()(Point_3 p, Vector_3 d)`
   * returning the line going through  `p` in the direction of `d`.
   */
  typedef unspecified_type Construct_line_3;
  
  /*!
   * Function object type that provides
   * `Point_3 operator()(const Line_3& l)`
   * returning an arbitrary point on `l`. It holds `point(i) == point(j)`,
   * iff `i == j`. Furthermore, it is directed from `point(i)` to `point(j)`,
   * for all `i < j`.
   */
  typedef unspecified_type Construct_point_on_3;

  /*!
   * Function object type that provides
   * `FT operator()(Point_3 p)` and `FT operator()(Vector_3 v)`
   * returning the `x` coordinate of a point and a vector respectively.
   */
  typedef unspecified_type Compute_x_3;

  /*!
   * Function object type that provides
   * `FT operator()(Point_3 p)` and `FT operator()(Vector_3 v)`
   * returning the `y` coordinate of a point and a vector respectively.
   */
  typedef unspecified_type Compute_y_3;

  /*!
   * Function object type that provides
   * `FT operator()(Point_3 p)` and `FT operator()(Vector_3 v)`
   * returning the `z` coordinate of a point and a vector respectively.
   */
  typedef unspecified_type Compute_z_3;
  
  /*!
   * Function object type that provides
   * `FT operator()(Vector_3 v)`
   * returning the squared length of `v`.
   */
  typedef unspecified_type Compute_squared_length_3;

  /*!
   * Function object type that provides
   * `Vector_3 operator()(Vector_3 v, FT t)`
   * returning the vector `t * v`.
   */
  typedef unspecified_type Construct_scaled_vector_3;
  
  /*!
   * Function object type that provides
   * `Vector_3 operator()(Vector_3 v1, Vector_3 v2)`
   * returning the `v1+v2`.
   */
  typedef unspecified_type Construct_sum_of_vectors_3;

  /*!
  * Function object type that provides:
  * `Point_3 operator()(Point_3 p, Vector_3 v)`
  * returning  the point obtained by translating `p` by the vector `v`. 
  */ 
  typedef unspecified_type Construct_translated_point_3; 

  /*!
   * Function object type that provides
   * `FT operator()(Vector_3 v1, Vector_3 v2)`
   * returning the scalar product of `v1` and `v2`.
   */
  typedef unspecified_type Compute_scalar_product_3;
 
  /*!
   * Function object type that provides
   * `Vector_3 operator()(Vector_3 v1, Vector_3 v2)`
   * returning the cross-product vector of `v1` and `v2`.
   */
  typedef unspecified_type Construct_cross_product_vector_3;
  
  /*!
   * Function object type that provides
   * `Point_3 operator()(const Sphere_3& s)`
   * returning the center of the sphere `s`.
   */
  typedef unspecified_type Construct_center_3;

  /*!
   * Function object type that provides
   * `FT operator()(const Sphere_3& s)`
   * returning the squared radius of the sphere `s`.
   */
  typedef unspecified_type Compute_squared_radius_3;

/// @}

/// \name Access to Function Objects
/// @{
  
  Construct_point_3
  construct_point_3_object();

  Construct_vector_3
  construct_vector_3_object();

  Construct_sphere_3
  construct_sphere_3_object();
  
  Construct_line_3
  construct_line_3_object();
  
  Construct_point_on_3
  construct_point_on_3_object();

  Compute_x_3
  compute_x_3_object();

  Compute_y_3
  compute_y_3_object();

  Compute_z_3
  compute_z_3_object();
  
  Compute_squared_length_3
  compute_squared_length_3_object();

  Construct_scaled_vector_3
  construct_scaled_vector_3_object();

  Construct_sum_of_vectors_3
  construct_sum_of_vectors_3_object();

  Compute_scalar_product_3
  compute_scalar_product_3_object();
  
  Construct_cross_product_vector_3
  construct_cross_product_vector_3_object();

  Construct_translated_point_3 
  construct_translated_point_3_object();

  Construct_center_3
  construct_center_3_object();

  Compute_squared_radius_3
  compute_squared_radius_3_object();

/// @}
};
