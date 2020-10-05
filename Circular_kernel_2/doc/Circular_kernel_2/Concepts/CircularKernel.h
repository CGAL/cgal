
/*!
\ingroup PkgCircularKernel2GeometricConcepts
\cgalConcept

\cgalRefines `Kernel`

\cgalHasModel `CGAL::Circular_kernel_2<Kernel,AlgebraicKernelForCircles>`
\cgalHasModel `CGAL::Exact_circular_kernel_2`

\sa `Kernel`

*/

class CircularKernel {
public:

/// \name Types
/// A model of `CircularKernel` is supposed to provide some basic
/// types and to define the following geometric objects Moreover, a
/// model of `CircularKernel` must provide predicates, constructions
/// and other functionalities.
/// @{

/*!
Model of `LinearKernel`.
*/
typedef unspecified_type Linear_kernel;

/*!
Model of `AlgebraicKernelForCircles`.
*/
typedef unspecified_type Algebraic_kernel;

/*!
Model of `RingNumberType`.
*/
typedef unspecified_type RT;

/*!
Model of `FieldNumberType`.
*/
typedef unspecified_type FT;

/*!
Model of `RootOf_2`.
*/
typedef unspecified_type Root_of_2;

/*!
Model of `AlgebraicKernelForCircles::RootForCircles_2_2`.
*/
typedef unspecified_type Root_for_circles_2_2;

/*!
Model of `AlgebraicKernelForCircles::Polynomial_1_2`.
*/
typedef unspecified_type Polynomial_1_2;

/*!
Model of `AlgebraicKernelForCircles::PolynomialForCircles_2_2`.
*/
typedef unspecified_type Polynomial_for_circles_2_2;

/*!
Model of `Kernel::Point_2`.
*/
typedef unspecified_type Point_2;

/*!
Model of `Kernel::Line_2`.
*/
typedef unspecified_type Line_2;

/*!
Model of `Kernel::Circle_2`.
*/
typedef unspecified_type Circle_2;

/*!
Model of `CircularKernel::LineArc_2`.
*/
typedef unspecified_type Line_arc_2;

/*!
Model of `CircularKernel::CircularArc_2`.
*/
typedef unspecified_type Circular_arc_2;

/*!
Model of `CircularKernel::CircularArcPoint_2`.
*/
typedef unspecified_type Circular_arc_point_2;

/// @}

/// \name Predicates
/// @{

/*!
Model of `CircularKernel::CompareX_2`.
*/
typedef unspecified_type Compare_x_2;

/*!
Model of `CircularKernel::CompareY_2`.
*/
typedef unspecified_type Compare_y_2;

/*!
Model of `CircularKernel::CompareXY_2`.
*/
typedef unspecified_type Compare_xy_2;

/*!
Model of `CircularKernel::Equal_2`.
*/
typedef unspecified_type Equal_2;

/*!
Model of `CircularKernel::CompareYatX_2`.
*/
typedef unspecified_type Compare_y_at_x_2;

/*!
Model of `CircularKernel::CompareYtoRight_2`.
*/
typedef unspecified_type Compare_y_to_right_2;

/*!
Model of `CircularKernel::HasOn_2`.
*/
typedef unspecified_type Has_on_2;

/*!
Model of `CircularKernel::DoOverlap_2`.
*/
typedef unspecified_type Do_overlap_2;

/*!
Model of `CircularKernel::DoIntersect_2`.
*/
typedef unspecified_type Do_intersect_2;

/*!
Model of `CircularKernel::BoundedSide_2`.
*/
typedef unspecified_type Bounded_side_2;

/*!
Model of `CircularKernel::HasOnBoundedSide_2`.
*/
typedef unspecified_type Has_on_bounded_side_2;

/*!
Model of `CircularKernel::HasOnUnboundedSide_2`.
*/
typedef unspecified_type Has_on_unbounded_side_2;

/*!
Model of `CircularKernel::InXRange_2`.
*/
typedef unspecified_type In_x_range_2;

/*!
Model of `CircularKernel::IsVertical_2`.
*/
typedef unspecified_type Is_vertical_2;

/*!
Model of `CircularKernel::IsXMonotone_2`.
*/
typedef unspecified_type Is_x_monotone_2;

/*!
Model of `CircularKernel::IsYMonotone_2`.
*/
typedef unspecified_type Is_y_monotone_2;

/// @}

/// \name Constructions
/// @{

/*!
Model of `CircularKernel::ConstructLine_2`.
*/
typedef unspecified_type Construct_line_2;

/*!
Model of `CircularKernel::ConstructCircle_2`.
*/
typedef unspecified_type Construct_circle_2;

/*!
Model of `CircularKernel::ConstructBbox_2`.
*/
typedef unspecified_type Construct_bbox_2;

/*!
Model of `CircularKernel::ConstructCircularArcPoint_2`.
*/
typedef unspecified_type Construct_circular_arc_point_2;

/*!
Model of `CircularKernel::ConstructLineArc_2`.
*/
typedef unspecified_type Construct_line_arc_2;

/*!
Model of `CircularKernel::ConstructCircularArc_2`.
*/
typedef unspecified_type Construct_circular_arc_2;

/*!
Model of `CircularKernel::ComputeCircularX_2`
*/
typedef unspecified_type Compute_circular_x_2;

/*!
Model of `CircularKernel::ComputeCircularY_2`
*/
typedef unspecified_type Compute_circular_y_2;

/*!
Model of `CircularKernel::ConstructCircularMinVertex_2`.
*/
typedef unspecified_type Construct_circular_min_vertex_2;

/*!
Model of `CircularKernel::ConstructCircularMaxVertex_2`.
*/
typedef unspecified_type Construct_circular_max_vertex_2;

/*!
Model of `CircularKernel::ConstructCircularSourceVertex_2`.
*/
typedef unspecified_type Construct_circular_source_vertex_2;

/*!
Model of `CircularKernel::ConstructCircularTargetVertex_2`.
*/
typedef unspecified_type Construct_circular_target_vertex_2;

/*!
Model of `CircularKernel::Intersect_2`.
*/
typedef unspecified_type Intersect_2;

/*!
Model of `CircularKernel::Split_2`.
*/
typedef unspecified_type Split_2;

/*!
Model of `CircularKernel::MakeXMonotone_2`.
*/
typedef unspecified_type Make_x_monotone_2;

/*!
Model of `CircularKernel::MakeXYMonotone_2`.
*/
typedef unspecified_type Make_xy_monotone_2;

/// @}

/// \name Link with the algebraic kernel
/// @{

/*!
Model of `CircularKernel::GetEquation`.
*/
typedef unspecified_type Get_equation;

/// @}

/// \name Operations
/// As in the `Kernel` concept, for each of the function objects
/// above, there must exist a member function that requires no
/// arguments and returns an instance of that function object. The
/// name of the member function is the uncapitalized name of the type
/// returned with the suffix `_object` appended. For example, for the
/// function object `CircularKernel::Construct_circular_arc_2` the
/// following member function must exist:
/// @{

/*!

*/
Construct_circular_arc_2 construct_circular_arc_2_object() const;

/// @}

}; /* end CircularKernel */

