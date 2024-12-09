
/*!
\ingroup PkgKernelDKernelConcept
\cgalConcept

The concept of a <I>kernel with lifting</I> is a small refinement of the
general kernel concept. It adds 2 functors, the meaning of which would be
unclear in kernels of fixed dimension.

\cgalRefines{Kernel_d}
\cgalHasModelsBegin
\cgalHasModels{CGAL::Cartesian_d<FieldNumberType>}
\cgalHasModels{CGAL::Homogeneous_d<RingNumberType>}
\cgalHasModelsEnd
*/
class KernelWithLifting_d {
public:

/// \name Constructions
/// @{

/*!
a model of `KernelWithLifting_d::Lift_to_paraboloid_d`
*/
typedef unspecified_type Lift_to_paraboloid_d;

/*!
a model of `KernelWithLifting_d::Project_along_d_axis_d`
*/
typedef unspecified_type Project_along_d_axis_d;

/// @}

/// \name Operations
/// The following member functions return function objects of the
/// types listed above.
/// @{

/*!

*/
KernelWithLifting_d::Lift_to_paraboloid_d lift_to_paraboloid_d_object() const;

/*!

*/
KernelWithLifting_d::Project_along_d_axis_d project_along_d_axis_d_object() const;

/// @}

}; /* end KernelWithLifting_d */

