
/*!
 * \ingroup PkgMeanCurvatureSkeleton3Concepts
 * \cgalConcept
 *
 * Traits class concept defining the requirement of the class `Mean_curvature_flow_skeletonization`.
 *
 * \cgalHasModel any \cgal Kernel
 *
 * \todo finish to fill this concept
 *
 */
class MeanCurvatureSkeletonizationTraits {
public:

/// \name Types
/// @{

/// The point type.
typedef unspecified_type Point_3;


/// The vector type.
typedef unspecified_type Vector_3;


/*!
 * Function object type that provides
 * `double operator()(Point_3 p, Point_3 q)`
 * returning the squared distance between `p` and `q`.
 */
typedef unspecified_type Compute_squared_distance_3;

/// @}

/// \name Access to Function Objects
/// @{

Compute_squared_distance_3
compute_squared_distance_3_object();

/// @}


}; /* end DelaunayTriangulationTraits_2 */

