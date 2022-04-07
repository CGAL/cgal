/*!
\ingroup PkgSurfaceMeshSegmentationConcepts
\cgalConcept

The concept `SegmentationGeomTraits` describes the set of requirements
of the geometric traits needed by the segmentation functions.

\cgalHasModel All the \cgal Kernels

\cgalRefines AABBGeomTraits

*/

class SegmentationGeomTraits {
public:

/// \name Number Type
/// @{
/*!
A number type model of `FieldNumberType` for which an overload of the function `to_double()` is defined
*/
typedef unspecified_type FT;
/// @}

/// \name Geometric Types
/// @{

/*!
The point type.
*/
typedef unspecified_type Point_3;

/*!
The vector type.
*/
typedef unspecified_type Vector_3;

/*!
The segment type.
*/
typedef unspecified_type Segment_3;

/*!
The plane type.
*/
typedef unspecified_type Plane_3;

/*!
The ray type.
*/
typedef unspecified_type Ray_3;

/// @}

/// \name Function Object Types
/// @{

/*!
Function object model of Kernel::Angle_3 for the aforementioned geometric types.
*/
typedef unspecified_type Angle_3;

/*!
Function object model of Kernel::ComputeScalarProduct_3 for the aforementioned geometric types.
*/
typedef unspecified_type Compute_scalar_product_3;

/*!
Function object model of Kernel::ComputeSquaredDistance_3 for the aforementioned geometric types.
*/
typedef unspecified_type Compute_squared_distance_3;

/*!
Function object model of Kernel::ConstructCrossProductVector_3 for the aforementioned geometric types.
*/
typedef unspecified_type Construct_cross_product_vector_3;

/*!
Function object model of Kernel::ConstructNormal_3 for the aforementioned geometric types.
*/
typedef unspecified_type Construct_normal_3;

/*!
Function object model of Kernel::ConstructCentroid_3 for the aforementioned geometric types.
*/
typedef unspecified_type Construct_centroid_3;

/*!
Function object model of Kernel::ConstructTranslatedPoint_3 for the aforementioned geometric types.
*/
typedef unspecified_type Construct_translated_point_3;

/*!
Function object model of Kernel::ConstructScaledVector_3 for the aforementioned geometric types.
*/
typedef unspecified_type Construct_scaled_vector_3;

/*!
Function object model of Kernel::ConstructSumOfVectors_3 for the aforementioned geometric types.
*/
typedef unspecified_type Construct_sum_of_vectors_3;

/*!
Function object model of Kernel::Collinear_3 for the aforementioned geometric types.
*/
class Collinear_3
{
public:
 //!determines if points `p`, `q` and `r` are collinear.
  bool operator()(const Point_3 &p,
                  const Point_3 &q,
                  const Point_3 &r);
};

/// @}

/// \name Operations
/// @{
/*! For each of the above function and predicate object types, `Func_obj_type`,
 * a function must exist with the name `func_obj_type_object` that creates
 * an instance of the function or predicate object type.

For example:
*/
Collinear_3
collinear_3_object()const
{return Collinear_3();}

/// @}
}; /* end SegmentationGeomTraits */

