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
  /// The sphere type, only required if you want to detect spheres or tori
  typedef unspecified_type Sphere_3;
  /// The line type, only required if you want to detect cylinders
  typedef unspecified_type Line_3;
  /// The plane type, only required if you want to detect planes
  typedef unspecified_type Plane_3;
  /// The 2D point type, only required if you want to detect tori
  typedef unspecified_type Point_2;
  /// The circle type, only required if you want to detect tori
  typedef unspecified_type Circle_2;
  /// The 2D vector type, only required if you want to detect tori
  typedef unspecified_type Vector_3;

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
   * `Point_3 operator()(Origin p)`
   * returning the point with 0, 0, 0 as Cartesian coordinates
   * and `Point_3 operator()(FT x, FT y, FT z)`
   * returning the point with `x`, `y` and `z` as Cartesian coordinates.
   */
  typedef unspecified_type Construct_point_3;

  /*!
   * Function object type that provides 
   * `Vector_3 operator()(Point_3 p1, Point_3 p2)` and 
   * `Vector_3 operator()(Origin p1, Point_3 p2)`
   * returning the vector `p1p2`, `Vector_3 operator()(NULL_VECTOR)` returning 
   * the null vector, and `Vector_3 operator()(Line_3 l)` returning 
   * a vector having the same direction as `l` 
   * (this last one is only required if you want to detect cylinders).
   */
  typedef unspecified_type Construct_vector_3;
  
  /*!
   * Function object type that provides `Sphere_3 operator()(Point_3 c, FT r)`
   * returning the sphere of center `p` and radius `r`.
   * Only required if you want to detect spheres or tori.
   */
  typedef unspecified_type Construct_sphere_3;
  
  /*!
   * Function object type that provides `Line_3 operator()(Point_3 p, Vector_3 d)`
   * returning the line going through  `p` in the direction of `d`.
   * Only required if you want to detect cylinders.
   */
  typedef unspecified_type Construct_line_3;
  
  /*!
   * Function object type that provides
   * `Point_3 operator()(Line_3 l, int i)`
   * returning an arbitrary point on `l`. `i` is not used and can vbe of 
   * any value.
   * Only required if you want to detect cylinders.
   */
  typedef unspecified_type Construct_point_on_3;
  
  /*!
   * Function object type that provides
   * `Point_2 operator()(FT x, FT y)`
   * returning the 2D point with `x` and `y` as Cartesian coordinates.
   * Only required if you want to detect tori.
   */
  typedef unspecified_type Construct_point_2;
  
  /*!
   * Function object type that provides 
   * `Vector_2 operator()(Point_2 p1, Point_2 p2)`
   * returning the vector `p1p2`, `Vector_2 operator()(NULL_VECTOR)` returning 
   * the null vector.
   * Only required if you want to detect tori.
   */
  typedef unspecified_type Construct_vector_2;
  
  /*!
   * Function object type that provides 
   * `Circle_2 operator()(Point_2 p1, Point_2 p2, Point_2 p3)`
   * returning the circle going through `p1`, `p2` and `p3`.
   * Only required if you want to detect tori.
   */
  typedef unspecified_type Construct_circle_2;

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
   * `FT operator()(Point_2 p)` and `FT operator()(Vector_2 v)`
   * returning the `x` coordinate of a point and a vector respectively.
   * Only required if you want to detect tori.
   */
  typedef unspecified_type Compute_x_2;

  /*!
   * Function object type that provides
   * `FT operator()(Point_2 p)` and `FT operator()(Vector_2 v)`
   * returning the `y` coordinate of a point and a vector respectively.
   * Only required if you want to detect tori.
   */
  typedef unspecified_type Compute_y_2;

  /*!
   * Function object type that provides
   * `FT operator()(Vector_3 v)`
   * returning the squared length of `v`.
   */
  typedef unspecified_type Compute_squared_length_3;
  
  /*!
   * Function object type that provides
   * `FT operator()(Vector_2 v)`
   * returning the squared length of `v`.
   * Only required if you want to detect tori.
   */
  typedef unspecified_type Compute_squared_length_2;

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
   * `Point_3 operator()(Sphere_3 s)`
   * returning the center of the sphere `s`.
   * Only required if you want to detect spheres or tori.
   */
  typedef unspecified_type Construct_center_3;
  
  /*!
   * Function object type that provides
   * `Point_2 operator()(Circle_2 c)`
   * returning the center of the circle `c`.
   * Only required if you want to detect tori.
   */
  typedef unspecified_type Construct_center_2;

  /*!
   * Function object type that provides
   * `FT operator()(Sphere_3 s)`
   * returning the squared radius of the sphere `s`.
   * Only required if you want to detect spheres or tori.
   */
  typedef unspecified_type Compute_squared_radius_3;
  
  /*!
   * Function object type that provides
   * `FT operator()(Circle_2 c)`
   * returning the squared radius of the circle `c`.
   * Only required if you want to detect tori.
   */
  typedef unspecified_type Compute_squared_radius_2;
  
  /*!
   * Function object type that provides
   * `bool operator(Point_2 p, Point_2 q, Point_2 r)`
   * returning true if the points `p`, `q`, and `r` are collinear and
   * false otherwise.
   * Only required if you want to detect tori.
  */ 
  typedef unspecified_type Collinear_2; 

/// @}

/// \name Access to Function Objects
/// @{
  
  Construct_point_3
  construct_point_3_object();

  Construct_vector_3
  construct_vector_3_object();

  /// Only required if you want to detect spheres or tori.
  Construct_sphere_3
  construct_sphere_3_object();
  
  /// Only required if you want to detect cylinders.
  Construct_line_3
  construct_line_3_object();
  
  /// Only required if you want to detect cylinders.
  Construct_point_on_3
  construct_point_on_3_object();
  
  /// Only required if you want to detect tori.
  Construct_point_2 
  construct_point_2_object();

  /// Only required if you want to detect tori.
  Construct_vector_2 
  construct_vector_2_object();

  /// Only required if you want to detect tori.
  Construct_circle_2 
  construct_circle_2_object();

  Compute_x_3
  compute_x_3_object();

  Compute_y_3
  compute_y_3_object();

  Compute_z_3
  compute_z_3_object();
  
  /// Only required if you want to detect tori.
  Compute_x_2
  compute_x_2_object();

  /// Only required if you want to detect tori.
  Compute_y_2
  compute_y_2_object();

  Compute_squared_length_3
  compute_squared_length_3_object();

  /// Only required if you want to detect tori.
  Compute_squared_length_2 
  compute_squared_length_2_object();

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
  
  /// Only required if you want to detect spheres.
  Construct_center_3
  construct_center_3_object();

  /// Only required if you want to detect tori.
  Construct_center_2 
  construct_center_2_object();
  
  /// Only required if you want to detect spheres or tori.
  Compute_squared_radius_3
  compute_squared_radius_3_object();

  /// Only required if you want to detect tori.
  Compute_squared_radius_2 
  compute_squared_radius_2_object();
  
  /// Only required if you want to detect tori.
  Collinear_2 
  collinear_2_object();

/// @}
};
