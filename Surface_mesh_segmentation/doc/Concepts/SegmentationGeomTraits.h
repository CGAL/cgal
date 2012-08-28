/*!
\ingroup PkgSurfaceMeshSegmentationConcepts
\cgalconcept

The concept `SegmentationGeomTraits` describes the set of requirements
of the geometric traits needed by the segmentation fonctions.

\hasModel All the \cgal Kernels 

\refines AABBGeomTraits

*/

class SegmentationGeomTraits {
public:

/// \name Number type
/// @{
/*! 
A number type model of FieldWithSqrt (double or float is recommanded)
*/ 
typedef Hidden_type NT;
/// @}

/// \name Geometric Types 
/// @{

/*! 
The vector type. 
*/ 
typedef Hidden_type Point_3; 

/*! 
The vector type. 
*/ 
typedef Hidden_type Vector_3; 

/*! 
The segment type. 
*/ 
typedef Hidden_type Segment_3; 

/*! 
The plane type. 
*/ 
typedef Hidden_type Plane_3; 

/*! 
The ray type. 
*/ 
typedef Hidden_type Ray_3; 

/// @}

/// \name Function Object Types
/// @{ 

/*! 
Function object model of Kernel::Angle_3 for the aforementioned geometric types.
*/
typedef Hidden_type Angle_3;

/*! 
Function object model of Kernel::ConstructNormal_3 for the aforementioned geometric types.
*/
typedef Hidden_type Construct_normal_3;

/*! 
Function object model of Kernel::ConstructCentroid_3 for the aforementioned geometric types.
*/
typedef Hidden_type Construct_centroid_3;

/*! 
Function object model of Kernel::ConstructUnitNormal_3 for the aforementioned geometric types.
*/
typedef Hidden_type Construct_unit_normal_3;

/*! 
Function object model of Kernel::ConstructTranslatedPoint_3 for the aforementioned geometric types.
*/
typedef Hidden_type Construct_translated_point_3;

/*! 
Function object model of Kernel::ConstructScaledVector_3 for the aforementioned geometric types.
*/
typedef Hidden_type Construct_scaled_vector_3;

/*! 
Function object model of Kernel::ConstructSumOfVectors_3 for the aforementioned geometric types.
*/
typedef Hidden_type Construct_sum_of_vectors_3;

/// @}

}; /* end SegmentationGeomTraits */

