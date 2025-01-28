
/*!
\ingroup PkgMesh3Concepts
\cgalConcept

The concept `MeshTriangulationTraits_3` describes the requirements for
the geometric traits class of the underlying regular triangulation used during
a mesh generation process.

\cgalRefines{RegularTriangulationTraits_3}

\cgalHasModelsBegin
\cgalHasModelsBare{All models of the \cgal concept `Kernel`}
\cgalHasModelsEnd

In addition to the requirements described for the traits class
RegularTriangulationTraits_3, the geometric traits class of a
regular triangulation used in mesh generation must fulfill the following
requirements.
*/
class MeshTriangulationTraits_3
{
public:

  /// \name Types
  /// @{

  /*!
  Numerical type. Must be a model of `::FieldNumberType` and
  `::FieldWithSqrt`, and constructible from a `double`.
  */
  typedef unspecified_type FT;

  /*!
  The vector type.
  */
  typedef unspecified_type Vector_3;

  /*!
  The sphere type.
  */
  typedef unspecified_type Sphere_3;

  /*!
  The isocuboid type.
  */
  typedef unspecified_type Iso_cuboid_3;

  /*!
  The affine transformation type.
  */
  typedef unspecified_type Aff_transformation_3;

  /*!
  A predicate object that must provide the function operator:

  `bool operator()(Sphere_3 s1, Sphere_3 s2)`

  which returns `true` iff the spheres `s1` and `s2` intersect.
  */
  typedef unspecified_type Do_intersect_3;

  /*!
  A predicate object that must provide the function operator:

  `bool operator()(Point_3 p, Point_3 q)`

  which returns `true` iff `p` and `q` are equal.
  */
  typedef unspecified_type Equal_3;

  /*!
  A predicate object that must provide the function operators:

  `bool operator()(Segment_3 s)`

  `bool operator()(Ray_3 s)`

  which return `true` iff the object is degenerate.
  */
  typedef unspecified_type Is_degenerate_3;

  /*!
  A predicate object that must provide the function operator:

  `bool operator()(Point_3 p, Point_3 q, Point_3 r)`

  which returns `true` iff `p`, `q`, and `r` are collinear.
  */
  typedef unspecified_type Collinear_3;

  /*!
  A constructor object that must provide the function operators:

  `bool operator()(Point_3 p, FT w)`

  `bool operator()(Point_3 p, Point_3 q, FT w)`

  `bool operator()(Point_3 p, Point_3 q, Point_3 r, FT w)`

  `bool operator()(Point_3 p, Point_3 q, Point_3 r, Point_3 s, FT w)`

  which compare the weight of the smallest sphere orthogonal to the input weighted point(s)
  with the input weight.
  */
  typedef unspecified_type Compare_weighted_squared_radius_3;

  /*!
  A constructor object that must provide the function operator:

  `FT operator()(Point_3 p, Point_3 q, Point_3 r, Point_3 s)`

  which returns an approximation of the signed dihedral angle in the tetrahedron `pqrs` of edge `pq`.
  The sign is negative if `orientation(p,q,r,s)` is `CGAL::NEGATIVE` and positive otherwise.
  The angle is given in degrees.
  */
  typedef unspecified_type Compute_approximate_dihedral_angle_3;

  /*!
  A constructor object that must provide the function operator:

  `FT operator()(Point_3 p, Point_3 q, Point_3 r)`

  which returns the area of the triangle `p`, `q`, `r`.
  */
  typedef unspecified_type Compute_area_3;

  /*!
  A constructor object that must provide the function operator:

  `FT operator()(Point_3 p, Point_3 q, Point_3 r, Point_3 s, Point_3 t)`

  which returns the squared radius of the sphere centered in `t` and orthogonal
  to the sphere orthogonal to `p`, `q`, `r`, and `s`.
  */
  typedef unspecified_type Compute_power_distance_to_power_sphere_3;

  /*!
  A constructor object that must provide the function operator:

  `FT operator()(Triangle_3 t)`

  which returns the squared area of the triangle `t`.
  */
  typedef unspecified_type Compute_squared_area_3;

  /*!
  A constructor object that must provide the function operator:

  `FT operator()(Point_3 p, Point_3 q)`

  which returns the squared distance between two points `p` and `q`.
  */
  typedef unspecified_type Compute_squared_distance_3;

  /*!
  A constructor object that must provide the function operator:

  `FT operator()(Vector_3 v)`

  which returns the squared length of `v`.
  */
  typedef unspecified_type Compute_squared_length_3;

  /*!
  A constructor object that must provide the function operator:

  `FT operator()(Sphere_3 s)`

   which returns the squared radius of the sphere `s`.
  */
  typedef unspecified_type Compute_squared_radius_3;

  /*!
  A constructor object that must provide the function operators:

  `FT operator()(Point_3 p, Point_3 q, Point_3 r, Point_3 s)`

  `FT operator()(Point_3 p, Point_3 q, Point_3 r)`

  `FT operator()(Point_3 p, Point_3 q)`

   which returns the squared radius of the smallest sphere orthogonal to the argument(s).
  */
  typedef unspecified_type Compute_squared_radius_smallest_orthogonal_sphere_3;

  /*!
  A constructor object that must provide the function operators:

  `FT operator()(Point_3 p, Point_3 q, Point_3 r, Point_3 s)`

  `FT operator()(Tetrahedron_3 t)`

  which respectively return the volume of the tetrahedron formed by the points
  `p`, `q`, `r`, and `s` and of the tetrahedron `t`.
  */
  typedef unspecified_type Compute_volume_3;

  /*!
  A constructor object that must provide the function operators:

  `FT operator()(Weighted_point_3 p)`

  which returns the weight of the weighted point `p`.
  */
  typedef unspecified_type Compute_weight_3;

  /*!
  A constructor object that must provide the function operator:

  `Vector_3 operator()(Plane_3 p, int index)`,

  which returns:
  - when `index == 1`: a vector `b1` that is orthogonal to the normal `n` to plane `h`;
  - when `index == 2`: a vector `b2` that is orthogonal to `n` and `b1` and such that
  for an arbitrary point `p` on the plane `h`, the orientation of `p`, `p + b1`, `p + b2`,
  and `p + n` is positive.
  */
  typedef unspecified_type Construct_base_vector_3;

  /*!
  A constructor object that must provide the function operator:

  `Vector_3 operator()(Plane_3 p)`,

  which returns a vector that is orthogonal to the plane `p` and directed to the positive side of `p`.
  */
  typedef unspecified_type Construct_orthogonal_vector_3;

  /*!
  A constructor object that must provide the function operator:

  `Point_3 operator()(Sphere_3 q)`,

  which returns the center of the sphere `s`.
  */
  typedef unspecified_type Construct_center_3;

  /*!
  A constructor object that must provide the function operators:

  `Point_3 operator()(Weighted_point_3 p, Weighted_point_3 q, Weighted_point_3 r, Weighted_point_3 s)`,

  `Point_3 operator()(Weighted_point_3 p, Weighted_point_3 q, Weighted_point_3 r, Weighted_point_3 s)`,

  `Point_3 operator()(Weighted_point_3 p, Weighted_point_3 q, Weighted_point_3 r, Weighted_point_3 s)`,

  which return the center of the smallest orthogonal sphere to the input weighted points.
  */
  typedef unspecified_type Construct_weighted_circumcenter_3;

  /*!
  A constructor object that must provide the function operators:

  `Point_3 operator()(Point_3 p, Point_3 q, Point_3 r)`,

  `Point_3 operator()(Tetrahedron_3 t)`,

  which respectively return the centroid of the points `p`, `q`, and `r`,
  and the centroid of the tetrahedron `t`.
  */
  typedef unspecified_type Construct_centroid_3;

  /*!
  A constructor object that must provide the function operator:

  `Point_3 operator()(Point_3 p, Point_3 q)`,

  which returns the midpoint of the segment pq.
  */
  typedef unspecified_type Construct_midpoint_3;

  /*!
  A constructor object that must provide the function operator:

  `Vector_3 operator()(Point_3 p, Point_3 q, Point_3 r)`,

  which returns the normal of the vectors `q-p` and `r-p`.
  */
  typedef unspecified_type Construct_normal_3;

  /*!
  A constructor object that must provide the function operator:

  `Sphere_3 operator()(Point_3 p, Point_3 q)`,

  which returns the smallest sphere which passes through the points p and q.
  */
  typedef unspecified_type Construct_sphere_3;

  /*!
  A constructor object that must provide the function operator:

  `Tetrahedron_3 operator()(Point_3 p, Point_3 q, Point_3 r, Point_3 s)`,

  which returns the tetrahedron with vertices `p`, `q`, `r`, and `s`.
  */
  typedef unspecified_type Construct_tetrahedron_3;

  /*!
  A constructor object that must provide the function operator:

  `Point_3 operator()(Point_3 p, Vector_3 v)`,

  which returns the point obtained by translating `p` by the vector `v`.
  */
  typedef unspecified_type Construct_translated_point_3;

  /*!
  A constructor object that must provide the function operator:

  `Vector_3 operator()(Point_3 p, Point_3 q)`,

  which returns a vector from two points.
  */
  typedef unspecified_type Construct_vector_3;

  /*!
  A constructor object that must provide the function operator:

  `Vector_3 operator()(Vector_3 v)`

  which returns `-v`.
  */
  typedef unspecified_type Construct_opposite_vector_3;

  /*!
  A constructor object that must provide the function operator:

  `Vector_3 operator()(Vector_3 v, FT scale)`

  which returns the vector `v` scaled by a factor `scale`.
  */
  typedef unspecified_type Construct_scaled_vector_3;

  /*!
  A constructor object that must provide the function operators:

  `std::optional< std::variant< T... > > operator()(Segment_3 s, Plane_3 p)`

  `std::optional< std::variant< T... > > operator()(Ray_3 r, Iso_cuboid i)`

  `std::optional< std::variant< T... > > operator()(Segment_3 s, Iso_cuboid i)`

  which returns the intersection region of two geometrical objects.
  */
  typedef unspecified_type Intersect_3;

  /// @}

  /// \name Operations
  /// The following functions give access to the predicate and construction objects:
  /// @{

  /*!

  */
  Do_intersect_3 do_intersect_3_object();

  /*!

  */
  Equal_3 equal_3_object();

  /*!

  */
  Is_degenerate_3 equal_3_object();

  /*!

  */
  Compare_weighted_squared_radius_3 compare_weighted_squared_radius_3_object();

  /*!

  */
  Compute_approximate_dihedral_angle_3 compute_approximate_dihedral_angle_3_object();

  /*!

  */
  Compute_area_3 compute_area_3_object();

  /*!

  */
  Compute_power_distance_to_power_sphere_3 compute_power_distance_to_power_sphere_3_object();

  /*!

  */
  Compute_squared_area_3 compute_squared_area_3_object();

  /*!

  */
  Compute_squared_distance_3 compute_squared_distance_3_object();

  /*!

  */
  Compute_squared_length_3 compute_squared_length_3_object();

  /*!

  */
  Compute_squared_radius_3 compute_squared_radius_3_object();

  /*!

  */
  Compute_squared_radius_smallest_orthogonal_sphere_3 compute_squared_radius_smallest_orthogonal_sphere_3_object();

  /*!

  */
  Compute_volume_3 compute_volume_3_object();

  /*!

  */
  Compute_weight_3 compute_weight_3_object();

  /*!

  */
  Construct_base_vector_3 construct_base_vector_3_object();

  /*!

  */
  Construct_orthogonal_vector_3 construct_orthogonal_vector_3_object();

  /*!

  */
  Construct_center_3 construct_center_3_object();

  /*!

  */
  Construct_midpoint_3 construct_midpoint_3_object();

  /*!

  */
  Construct_normal_3 construct_normal_3_object();

  /*!

  */
  Construct_sphere_3 construct_sphere_3_object();

  /*!

  */
  Construct_tetrahedron_3 construct_tetrahedron_3_object();

  /*!

  */
  Construct_translated_point_3 construct_translated_point_3_object();

  /*!

  */
  Construct_vector_3 construct_vector_3_object();

  /*!

  */
  Construct_scaled_vector_3 construct_scaled_vector_3_object();

  /*!

  */
  Construct_opposite_vector_3 construct_opposite_vector_3_object();

  /*!

  */
  Intersect_3 intersect_3_object();

  /// @}

}; /* end MeshTriangulationTraits_3 */

