
/*!
\ingroup PkgCircularKernel3GeometricConcepts
\cgalConcept

\cgalRefines{Kernel}

\cgalHasModelsBegin
\cgalHasModels{CGAL::Spherical_kernel_3<Kernel,AlgebraicKernelForSpheres>}
\cgalHasModels{CGAL::Exact_spherical_kernel_3}
\cgalHasModelsEnd

\sa `Kernel`

*/

class SphericalKernel {
public:

/// \name Types
/// A model of `SphericalKernel` must provide some basic
/// types
/// @{

/*!
Model of `Kernel`.
*/
typedef unspecified_type Linear_kernel;

/*!
Model of `AlgebraicKernelForSpheres`.
*/
typedef unspecified_type Algebraic_kernel;

/*!
Model of `RootOf_2`.
*/
typedef unspecified_type Root_of_2;

/*!
Model of `AlgebraicKernelForSpheres::RootForSpheres_2_3`.
*/
typedef unspecified_type Root_for_spheres_2_3;

/*!
Model of `AlgebraicKernelForSpheres::Polynomial_1_3`.
*/
typedef unspecified_type Polynomial_1_3;

/*!
Model of `AlgebraicKernelForSpheres::PolynomialsForLines_3`.
*/
typedef unspecified_type Polynomials_for_line_3;

/*!
Model of `AlgebraicKernelForSpheres::PolynomialForSpheres_2_3`.
*/
typedef unspecified_type Polynomial_for_spheres_2_3;

/*!
Model of `AlgebraicKernelForSpheres::PolynomialsForCircles_3`.
*/
typedef unspecified_type Polynomials_for_circle_3;

/// @}

/// \name
/// and to define the following geometric objects
/// @{

/*!
Model of `Kernel::Point_3`.
*/
typedef unspecified_type Point_3;

/*!
Model of `Kernel::Vector_3`.
*/
typedef unspecified_type Vector_3;

/*!
Model of `Kernel::Line_3`.
*/
typedef unspecified_type Line_3;

/*!
Model of `Kernel::Plane_3`.
*/
typedef unspecified_type Plane_3;

/*!
Model of `Kernel::Sphere_3`.
*/
typedef unspecified_type Sphere_3;

/*!
Model of `Kernel::Circle_3`.
*/
typedef unspecified_type Circle_3;

/*!
Model of `SphericalKernel::LineArc_3`.
*/
typedef unspecified_type Line_arc_3;

/*!
Model of `SphericalKernel::CircularArc_3`.
*/
typedef unspecified_type Circular_arc_3;

/*!
Model of `SphericalKernel::CircularArcPoint_3`.
*/
typedef unspecified_type Circular_arc_point_3;

/// @}


/// \name Predicates
/// Moreover, a model of `SphericalKernel` must provide predicates,
/// constructions and other functionalities.
/// @{

/*!
Model of `SphericalKernel::CompareX_3`.
*/
typedef unspecified_type Compare_x_3;

/*!
Model of `SphericalKernel::CompareY_3`.
*/
typedef unspecified_type Compare_y_3;

/*!
Model of `SphericalKernel::CompareZ_3`.
*/
typedef unspecified_type Compare_z_3;

/*!
Model of `SphericalKernel::CompareXY_3`.
*/
typedef unspecified_type Compare_xy_3;

/*!
Model of `SphericalKernel::CompareXYZ_3`.
*/
typedef unspecified_type Compare_xyz_3;

/*!
Model of `SphericalKernel::CompareTheta_3`.
*/
typedef unspecified_type Compare_theta_3;

/*!
Model of `SphericalKernel::CompareThetaZ_3`.
*/
typedef unspecified_type Compare_theta_z_3;

/*!
Model of `SphericalKernel::CompareZAtTheta_3`.
*/
typedef unspecified_type Compare_z_at_theta_3;

/*!
Model of `SphericalKernel::CompareZToRight_3`.
*/
typedef unspecified_type Compare_z_to_right_3;

/*!
Model of `SphericalKernel::Equal_3`.
*/
typedef unspecified_type Equal_3;

/*!
Model of `SphericalKernel::HasOn_3`.
*/
typedef unspecified_type Has_on_3;

/*!
Model of `SphericalKernel::DoOverlap_3`.
*/
typedef unspecified_type Do_overlap_3;

/*!
Model of `SphericalKernel::DoIntersect_3`.
*/
typedef unspecified_type Do_intersect_3;

/*!
Model of `SphericalKernel::BoundedSide_3`.
*/
typedef unspecified_type Bounded_side_3;

/*!
Model of `SphericalKernel::HasOnBoundedSide_3`.
*/
typedef unspecified_type Has_on_bounded_side_3;

/*!
Model of `SphericalKernel::HasOnUnboundedSide_3`.
*/
typedef unspecified_type Has_on_unbounded_side_3;

/*!
Model of `SphericalKernel::IsThetaMonotone_3`.
*/
typedef unspecified_type Is_theta_monotone_3;

/// @}

/// \name Constructions
/// @{

/*!
Model of `SphericalKernel::ConstructLine_3`.
*/
typedef unspecified_type Construct_line_3;

/*!
Model of `SphericalKernel::ConstructPlane_3`.
*/
typedef unspecified_type Construct_plane_3;

/*!
Model of `SphericalKernel::ConstructSphere_3`.
*/
typedef unspecified_type Construct_sphere_3;

/*!
Model of `SphericalKernel::ConstructCircle_3`.
*/
typedef unspecified_type Construct_circle_3;

/*!
Model of `SphericalKernel::ConstructLineArc_3`.
*/
typedef unspecified_type Construct_line_arc_3;

/*!
Model of `SphericalKernel::ConstructCircularArc_3`.
*/
typedef unspecified_type Construct_circular_arc_3;

/*!
Model of `SphericalKernel::ConstructCircularArcPoint_3`.
*/
typedef unspecified_type Construct_circular_arc_point_3;

/*!
Model of `SphericalKernel::ConstructCircularMinVertex_3`.
*/
typedef unspecified_type Construct_circular_min_vertex_3;

/*!
Model of `SphericalKernel::ConstructCircularMaxVertex_3`.
*/
typedef unspecified_type Construct_circular_max_vertex_3;

/*!
Model of `SphericalKernel::ConstructCircularSourceVertex_3`.
*/
typedef unspecified_type Construct_circular_source_vertex_3;

/*!
Model of `SphericalKernel::ConstructCircularTargetVertex_3`.
*/
typedef unspecified_type Construct_circular_target_vertex_3;

/*!
Model of `SphericalKernel::ConstructBbox_3`.
*/
typedef unspecified_type Construct_bbox_3;

/*!
Model of `SphericalKernel::Intersect_3`.
*/
typedef unspecified_type Intersect_3;

/*!
Model of `SphericalKernel::Split_3`.
*/
typedef unspecified_type Split_3;

/*!
Model of `SphericalKernel::MakeThetaMonotone_3`.
*/
typedef unspecified_type Make_theta_monotone_3;

/// @}

/// \name Computations
/// @{

/*!
Model of `SphericalKernel::ComputeCircularX_3`.
*/
typedef unspecified_type Compute_circular_x_3;

/*!
Model of `SphericalKernel::ComputeCircularY_3`.
*/
typedef unspecified_type Compute_circular_y_3;

/*!
Model of `SphericalKernel::ComputeCircularZ_3`.
*/
typedef unspecified_type Compute_circular_z_3;

/*!
Model of `SphericalKernel::ComputeApproximateSquaredLength_3`.
*/
typedef unspecified_type Compute_approximate_squared_length_3;

/*!
Model of `SphericalKernel::ComputeApproximateAngle_3`.
*/
typedef unspecified_type Compute_approximate_angle_3;

/// @}

/// \name Link with the algebraic kernel
/// @{

/*!
Model of `SphericalKernel::GetEquation`.
*/
typedef unspecified_type Get_equation;

/// @}

/// \name Operations
/// As in the `Kernel` concept, for each of the function objects
/// above, there must exist a member function that requires no
/// arguments and returns an instance of that function object. The
/// name of the member function is the uncapitalized name of the type
/// returned with the suffix `_object` appended. For example, for the
/// function object `SphericalKernel::Construct_circular_arc_3` the
/// following member function must exist:
/// @{

/*!

*/
Construct_circular_arc_3 construct_circular_arc_3_object() const;

/// @}

/// \name Operations on a given sphere,
/// A <I>context</I> sphere must be provided to the following functions:
/// @{

/*!

*/
Compare_theta_3 compare_theta_3_object(const Sphere_3& sphere) const;

/*!

*/
Compare_theta_z_3 compare_theta_z_3_object(const Sphere_3& sphere) const;

/*!

*/
Compare_z_at_theta_3 compare_z_at_theta_3_object(const Sphere_3& sphere) const;

/*!

*/
Compare_z_to_right_3 compare_z_to_right_3_object(const Sphere_3& sphere) const;

/*!

*/
Make_theta_monotone_3 make_theta_monotone_3_object(const Sphere_3& sphere) const;

/*!

*/
Is_theta_monotone_3 is_theta_monotone_3_object(const Sphere_3& sphere) const;

/// @}

}; /* end SphericalKernel */

