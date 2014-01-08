namespace Kernel {

/*!
  \ingroup PkgKernel23ConceptsGeomObject
  \cgalConcept

  A type representing circles in two dimensions. 

  \cgalRefines CopyConstructible
  \cgalRefines Assignable
  \cgalRefines DefaultConstructible

  \sa `CGAL::Circle_2<Kernel>` 
  \sa `Kernel::BoundedSide_2` 
  \sa `Kernel::ComputeSquaredRadius_2` 
  \sa `Kernel::ConstructCenter_2` 
  \sa `Kernel::ConstructCircle_2` 
  \sa `Kernel::ConstructOppositeCircle_2` 
  \sa `Kernel::Equal_2` 
  \sa `Kernel::HasOnBoundary_2` 
  \sa `Kernel::HasOnBoundedSide_2` 
  \sa `Kernel::HasOnNegativeSide_2` 
  \sa `Kernel::HasOnPositiveSide_2` 
  \sa `Kernel::HasOnUnboundedSide_2` 
  \sa `Kernel::IsDegenerate_2` 
  \sa `Kernel::OrientedSide_2` 

*/
class Circle_2 {
public:

}; /* end Kernel::Circle_2 */

/*!
  \ingroup PkgKernel23ConceptsGeomObject
  \cgalConcept

  A type representing circles in three dimensions. 

  \cgalRefines CopyConstructible
  \cgalRefines Assignable
  \cgalRefines DefaultConstructible

  \sa `CGAL::Circle_3<Kernel>` 
  \sa `Kernel::ComputeApproximateArea_3` 
  \sa `Kernel::ComputeApproximateSquaredLength_3` 
  \sa `Kernel::ComputeAreaDividedByPi_3` 
  \sa `Kernel::ComputeSquaredLengthDividedByPiSquare_3` 
  \sa `Kernel::ComputeSquaredRadius_3` 
  \sa `Kernel::ConstructBbox_3` 
  \sa `Kernel::ConstructCenter_3` 
  \sa `Kernel::ConstructCircle_3` 
  \sa `Kernel::ConstructSphere_3` 
  \sa `Kernel::ConstructPlane_3` 
  \sa `Kernel::Equal_3` 
  \sa `Kernel::HasOn_3` 
  \sa `Kernel::IsDegenerate_3` 

*/
class Circle_3 {
public:
}; /* end Kernel::Circle_3 */

/*!
  \ingroup PkgKernel23ConceptsGeomObject
  \cgalConcept

  A type representing directions in two dimensions. 

  \cgalRefines CopyConstructible
  \cgalRefines Assignable
  \cgalRefines DefaultConstructible

  \sa `CGAL::Direction_2<Kernel>` 
  \sa `Kernel::CompareAngleWithXAxis_2` 
  \sa `Kernel::ComputeDx_2` 
  \sa `Kernel::ComputeDy_2` 
  \sa `Kernel::ConstructDirection_2` 
  \sa `Kernel::ConstructOppositeDirection_2` 
  \sa `Kernel::ConstructPerpendicularDirection_2` 
  \sa `Kernel::CounterclockwiseInBetween_2` 
  \sa `Kernel::Equal_2` 

*/
class Direction_2 {
public:
}; /* end Kernel::Direction_2 */

/*!
  \ingroup PkgKernel23ConceptsGeomObject
  \cgalConcept

  A type representing directions in three dimensions. 

  \cgalRefines CopyConstructible
  \cgalRefines Assignable
  \cgalRefines DefaultConstructible

  \sa `CGAL::Direction_3<Kernel>` 
  \sa `Kernel::ConstructDirection_3` 
  \sa `Kernel::ConstructOppositeDirection_3` 
  \sa `Kernel::Equal_2` 

*/
class Direction_3 {
public:

}; /* end Kernel::Direction_3 */

/*!
  \ingroup PkgKernel23ConceptsGeomObject
  \cgalConcept

A type representing isocuboids in three dimensions. 

\cgalRefines CopyConstructible
\cgalRefines Assignable
\cgalRefines DefaultConstructible

\sa `CGAL::Iso_cuboid_3<Kernel>` 
\sa `Kernel::BoundedSide_3` 
\sa `Kernel::ComputeVolume_3` 
\sa `Kernel::ConstructIsoCuboid_3` 
\sa `Kernel::ConstructVertex_3` 
\sa `Kernel::Equal_2` 
\sa `Kernel::HasOnBoundary_3` 
\sa `Kernel::HasOnBoundedSide_3` 
\sa `Kernel::HasOnUnboundedSide_3` 
\sa `Kernel::IsDegenerate_3` 

*/
class IsoCuboid_3 {
public:
}; /* end Kernel::IsoCuboid_3 */

/*!
  \ingroup PkgKernel23ConceptsGeomObject
  \cgalConcept

  A type representing iso-rectangles in two dimensions. 

  \cgalRefines CopyConstructible
  \cgalRefines Assignable
  \cgalRefines DefaultConstructible

  \sa `CGAL::Iso_rectangle_2<Kernel>` 
  \sa `Kernel::ConstructIsoRectangle_2` 
  \sa `Kernel::ComputeXmin_2` 
  \sa `Kernel::ComputeXmax_2` 
  \sa `Kernel::ComputeYmin_2` 
  \sa `Kernel::ComputeYmax_2` 
  \sa `Kernel::BoundedSide_2` 
  \sa `Kernel::ComputeArea_2` 
  \sa `Kernel::ConstructIsoRectangle_2` 
  \sa `Kernel::ConstructVertex_2` 
  \sa `Kernel::DoIntersect_2` 
  \sa `Kernel::Equal_2` 
  \sa `Kernel::HasOnBoundary_2` 
  \sa `Kernel::HasOnBoundedSide_2` 
  \sa `Kernel::HasOnUnboundedSide_2` 
  \sa `Kernel::Intersect_2` 
  \sa `Kernel::IsDegenerate_2` 

*/
class IsoRectangle_2 {
public:
}; /* end Kernel::IsoRectangle_2 */

/*!
  \ingroup PkgKernel23ConceptsGeomObject
  \cgalConcept

  A type representing straight lines (and halfspaces) in two dimensions. 

  \cgalRefines CopyConstructible
  \cgalRefines Assignable
  \cgalRefines DefaultConstructible

  \sa `CGAL::Line_2<Kernel>` 
  \sa `Kernel::CompareXAtY_2` 
  \sa `Kernel::ComputeSquaredDistance_2` 
  \sa `Kernel::CompareYAtX_2` 
  \sa `Kernel::ConstructBisector_2` 
  \sa `Kernel::ConstructDirection_2` 
  \sa `Kernel::ConstructLine_2` 
  \sa `Kernel::ConstructOppositeLine_2` 
  \sa `Kernel::ConstructPerpendicularLine_2` 
  \sa `Kernel::ConstructPointOn_2` 
  \sa `Kernel::ConstructProjectedPoint_2` 
  \sa `Kernel::DoIntersect_2` 
  \sa `Kernel::Equal_2` 
  \sa `Kernel::HasOnNegativeSide_2` 
  \sa `Kernel::HasOnPositiveSide_2` 
  \sa `Kernel::HasOn_2` 
  \sa `Kernel::Intersect_2` 
  \sa `Kernel::IsDegenerate_2` 
  \sa `Kernel::IsHorizontal_2` 
  \sa `Kernel::IsVertical_2` 
  \sa `Kernel::OrientedSide_2` 

*/
class Line_2 {
public:
}; /* end Kernel::Line_2 */

/*!
  \ingroup PkgKernel23ConceptsGeomObject
  \cgalConcept

  A type representing straight lines in three dimensions. 

  \cgalRefines CopyConstructible
  \cgalRefines Assignable
  \cgalRefines DefaultConstructible

  \sa `CGAL::Line_3<Kernel>` 
  \sa `Kernel::ComputeSquaredDistance_3` 
  \sa `Kernel::ConstructDirection_3` 
  \sa `Kernel::ConstructLine_3` 
  \sa `Kernel::ConstructOppositeLine_3` 
  \sa `Kernel::ConstructPerpendicularLine_3` 
  \sa `Kernel::ConstructPlane_3` 
  \sa `Kernel::ConstructPointOn_3` 
  \sa `Kernel::ConstructProjectedPoint_3` 
  \sa `Kernel::DoIntersect_3` 
  \sa `Kernel::Equal_3` 
  \sa `Kernel::HasOn_3` 
  \sa `Kernel::Intersect_3` 
  \sa `Kernel::IsDegenerate_3` 

*/
class Line_3 {
public:
}; /* end Kernel::Line_3 */

/*!
  \ingroup PkgKernel23ConceptsGeomObject
  \cgalConcept

  A type representing different types of objects in two dimensions. 

  \deprecated This class is deprecated since \cgal 4.3 and type safe ways should be preferred. 

  \cgalRefines CopyConstructible
  \cgalRefines Assignable
  \cgalRefines DefaultConstructible

  \sa `CGAL::Object` 
  \sa `Kernel::Assign_2` 
  \sa `Kernel::ConstructObject_2` 
  \sa `Kernel::Intersect_2` 

*/
class Object_2 {
public:
}; /* end Kernel::Object_2 */

/*!
  \ingroup PkgKernel23ConceptsGeomObject
  \cgalConcept

  A type representing different types of objects in three dimensions. 

  \deprecated This class is deprecated since \cgal 4.3 and type safe ways should be preferred. 

  \cgalRefines CopyConstructible
  \cgalRefines Assignable
  \cgalRefines DefaultConstructible

  \sa `CGAL::Object` 
  \sa `Kernel::Assign_3` 
  \sa `Kernel::ConstructObject_3` 
  \sa `Kernel::Intersect_3` 

*/
class Object_3 {
public:
}; /* end Kernel::Object_3 */

/*!
  \ingroup PkgKernel23ConceptsGeomObject
  \cgalConcept

  A type representing planes (and half-spaces) in three dimensions. 

  \cgalRefines CopyConstructible
  \cgalRefines Assignable
  \cgalRefines DefaultConstructible

  \sa `CGAL::Plane_3<Kernel>` 
  \sa `Kernel::ComputeSquaredDistance_3` 
  \sa `Kernel::ConstructBaseVector_3` 
  \sa `Kernel::ConstructBisector_3` 
  \sa `Kernel::ConstructLiftedPoint_3` 
  \sa `Kernel::ConstructOppositePlane_3` 
  \sa `Kernel::ConstructOrthogonalVector_3` 
  \sa `Kernel::ConstructPerpendicularLine_3` 
  \sa `Kernel::ConstructPerpendicularPlane_3` 
  \sa `Kernel::ConstructPlane_3` 
  \sa `Kernel::ConstructPointOn_3` 
  \sa `Kernel::ConstructProjectedPoint_3` 
  \sa `Kernel::ConstructProjectedXYPoint_2` 
  \sa `Kernel::DoIntersect_3` 
  \sa `Kernel::Equal_3` 
  \sa `Kernel::HasOnNegativeSide_3` 
  \sa `Kernel::HasOnPositiveSide_3` 
  \sa `Kernel::HasOn_3` 
  \sa `Kernel::Intersect_3` 
  \sa `Kernel::IsDegenerate_3` 
  \sa `Kernel::LessSignedDistanceToPlane_3` 
  \sa `Kernel::OrientedSide_3` 

*/
class Plane_3 {
public:
}; /* end Kernel::Plane_3 */

/*!
  \ingroup PkgKernel23ConceptsGeomObject
  \cgalConcept

  A type representing points in two dimensions. 

  \cgalRefines CopyConstructible
  \cgalRefines Assignable
  \cgalRefines DefaultConstructible 

  \sa `Kernel::Angle_2` 
  \sa `Kernel::AreOrderedAlongLine_2` 
  \sa `Kernel::AreStrictlyOrderedAlongLine_2` 
  \sa `Kernel::Collinear_2` 
  \sa `Kernel::CollinearAreOrderedAlongLine_2` 
  \sa `Kernel::CollinearAreStrictlyOrderedAlongLine_2` 
  \sa `Kernel::CompareDistance_2` 
  \sa `Kernel::CompareXAtY_2` 
  \sa `Kernel::CompareXY_2` 
  \sa `Kernel::CompareX_2` 
  \sa `Kernel::CompareYAtX_2` 
  \sa `Kernel::CompareY_2` 
  \sa `Kernel::CompareYX_2` 
  \sa `Kernel::ComputeSquaredDistance_2` 
  \sa `Kernel::ComputeSquaredRadius_2` 
  \sa `Kernel::ComputeX_2` 
  \sa `Kernel::ComputeY_2` 
  \sa `Kernel::ComputeHx_2` 
  \sa `Kernel::ComputeHy_2` 
  \sa `Kernel::ConstructBisector_2` 
  \sa `Kernel::ConstructCircumcenter_2` 
  \sa `Kernel::ConstructLiftedPoint_3` 
  \sa `Kernel::ConstructMidpoint_2` 
  \sa `Kernel::ConstructPointOn_2` 
  \sa `Kernel::ConstructPoint_2` 
  \sa `Kernel::ConstructProjectedPoint_2` 
  \sa `Kernel::ConstructProjectedXYPoint_2` 
  \sa `Kernel::ConstructTranslatedPoint_2` 
  \sa `Kernel::DoIntersect_2` 
  \sa `Kernel::Equal_2` 
  \sa `Kernel::EqualX_2` 
  \sa `Kernel::EqualY_2` 
  \sa `Kernel::LeftTurn_2` 
  \sa `Kernel::LessDistanceToPoint_2` 
  \sa `Kernel::LessRotateCCW_2` 
  \sa `Kernel::LessSignedDistanceToLine_2` 
  \sa `Kernel::LessX_2` 
  \sa `Kernel::LessXY_2` 
  \sa `Kernel::LessY_2` 
  \sa `Kernel::LessYX_2` 
  \sa `Kernel::Orientation_2` 
  \sa `Kernel::SideOfBoundedCircle_2` 
  \sa `Kernel::SideOfOrientedCircle_2` 

*/
class Point_2 {
public:
}; /* end Kernel::Point_2 */

/*!
  \ingroup PkgKernel23ConceptsGeomObject
  \cgalConcept

  A type representing points in three dimensions. 

  \cgalRefines CopyConstructible
  \cgalRefines Assignable
  \cgalRefines DefaultConstructible

  \sa `Kernel::Angle_3` 
  \sa `Kernel::AreOrderedAlongLine_3` 
  \sa `Kernel::AreStrictlyOrderedAlongLine_3` 
  \sa `Kernel::Collinear_3` 
  \sa `Kernel::CollinearAreOrderedAlongLine_3` 
  \sa `Kernel::CollinearAreStrictlyOrderedAlongLine_3` 
  \sa `Kernel::CompareDistance_3` 
  \sa `Kernel::CompareXYZ_3` 
  \sa `Kernel::CompareXY_3` 
  \sa `Kernel::CompareX_3` 
  \sa `Kernel::CompareY_3` 
  \sa `Kernel::CompareZ_3` 
  \sa `Kernel::ComputeSquaredDistance_3` 
  \sa `Kernel::ComputeSquaredRadius_3` 
  \sa `Kernel::ComputeX_3` 
  \sa `Kernel::ComputeY_3` 
  \sa `Kernel::ComputeZ_3` 
  \sa `Kernel::ConstructBisector_3` 
  \sa `Kernel::ConstructCentroid_3` 
  \sa `Kernel::ConstructCircumcenter_3` 
  \sa `Kernel::ConstructLiftedPoint_3` 
  \sa `Kernel::ConstructMidpoint_3` 
  \sa `Kernel::ConstructPointOn_3` 
  \sa `Kernel::ConstructPoint_3` 
  \sa `Kernel::ConstructProjectedPoint_3` 
  \sa `Kernel::ConstructTranslatedPoint_3` 
  \sa `Kernel::CoplanarOrientation_3` 
  \sa `Kernel::CoplanarSideOfBoundedCircle_3` 
  \sa `Kernel::Coplanar_3` 
  \sa `Kernel::EqualXY_3` 
  \sa `Kernel::EqualX_3` 
  \sa `Kernel::EqualY_3` 
  \sa `Kernel::EqualZ_3` 
  \sa `Kernel::Equal_3` 
  \sa `Kernel::LessDistanceToPoint_3` 
  \sa `Kernel::LessSignedDistanceToPlane_3` 
  \sa `Kernel::LessXYZ_3` 
  \sa `Kernel::LessXY_3` 
  \sa `Kernel::LessX_3` 
  \sa `Kernel::LessY_3` 
  \sa `Kernel::LessZ_3` 
  \sa `Kernel::Orientation_3` 
  \sa `Kernel::SideOfBoundedSphere_3` 
  \sa `Kernel::SideOfOrientedSphere_3` 

*/
class Point_3 {
public:
}; /* end Kernel::Point_3 */

/*!
  \ingroup PkgKernel23ConceptsGeomObject
\cgalConcept

A type representing rays in two dimensions. 

\cgalRefines CopyConstructible
\cgalRefines Assignable
\cgalRefines DefaultConstructible

\sa `CGAL::Ray_2<Kernel>` 
\sa `Kernel::CollinearHasOn_2` 
\sa `Kernel::ComputeSquaredDistance_2` 
\sa `Kernel::ConstructDirection_2` 
\sa `Kernel::ConstructLine_2` 
\sa `Kernel::ConstructOppositeRay_2` 
\sa `Kernel::ConstructPointOn_2` 
\sa `Kernel::ConstructRay_2` 
\sa `Kernel::ConstructSource_2` 
\sa `Kernel::ConstructSecondPoint_2` 
\sa `Kernel::DoIntersect_2` 
\sa `Kernel::Equal_2` 
\sa `Kernel::HasOn_2` 
\sa `Kernel::Intersect_2` 
\sa `Kernel::IsDegenerate_2` 
\sa `Kernel::IsHorizontal_2` 
\sa `Kernel::IsVertical_2` 

*/
class Ray_2 {
public:
}; /* end Kernel::Ray_2 */

/*!
  \ingroup PkgKernel23ConceptsGeomObject
  \cgalConcept

  A type representing rays in three dimensions. 

  \cgalRefines CopyConstructible
  \cgalRefines Assignable
  \cgalRefines DefaultConstructible

  \sa `CGAL::Ray_3<Kernel>` 
  \sa `Kernel::ComputeSquaredDistance_3` 
  \sa `Kernel::ConstructDirection_3` 
  \sa `Kernel::ConstructLine_3` 
  \sa `Kernel::ConstructOppositeRay_3` 
  \sa `Kernel::ConstructPlane_3` 
  \sa `Kernel::ConstructPointOn_3` 
  \sa `Kernel::ConstructRay_3` 
  \sa `Kernel::DoIntersect_3` 
  \sa `Kernel::Equal_3` 
  \sa `Kernel::HasOn_3` 
  \sa `Kernel::Intersect_3` 
  \sa `Kernel::IsDegenerate_3` 

*/
class Ray_3 {
public:
}; /* end Kernel::Ray_3 */

/*!
  \ingroup PkgKernel23ConceptsGeomObject
  \cgalConcept

  A type representing segments in two dimensions. 

  \cgalRefines CopyConstructible
  \cgalRefines Assignable
  \cgalRefines DefaultConstructible

  \sa `CGAL::Segment_2<Kernel>` 
  \sa `Kernel::CollinearHasOn_2` 
  \sa `Kernel::ComputeSquaredDistance_2` 
  \sa `Kernel::ComputeSquaredLength_2` 
  \sa `Kernel::ConstructDirection_2` 
  \sa `Kernel::ConstructLine_2` 
  \sa `Kernel::ConstructOppositeSegment_2` 
  \sa `Kernel::ConstructPointOn_2` 
  \sa `Kernel::ConstructSegment_2` 
  \sa `Kernel::ConstructSource_2` 
  \sa `Kernel::ConstructTarget_2` 
  \sa `Kernel::ConstructVertex_2` 
  \sa `Kernel::DoIntersect_2` 
  \sa `Kernel::Equal_2` 
  \sa `Kernel::HasOn_2` 
  \sa `Kernel::Intersect_2` 
  \sa `Kernel::IsDegenerate_2` 
  \sa `Kernel::IsHorizontal_2` 
  \sa `Kernel::IsVertical_2` 

*/
class Segment_2 {
public:
}; /* end Kernel::Segment_2 */

/*!
  \ingroup PkgKernel23ConceptsGeomObject
  \cgalConcept

  A type representing segments in three dimensions. 

  \cgalRefines CopyConstructible
  \cgalRefines Assignable
  \cgalRefines DefaultConstructible

  \sa `CGAL::Segment_3<Kernel>` 
  \sa `Kernel::ComputeSquaredDistance_3` 
  \sa `Kernel::ComputeSquaredLength_3` 
  \sa `Kernel::ConstructDirection_3` 
  \sa `Kernel::ConstructLine_3` 
  \sa `Kernel::ConstructOppositeSegment_3` 
  \sa `Kernel::ConstructPlane_3` 
  \sa `Kernel::ConstructPointOn_3` 
  \sa `Kernel::ConstructSegment_3` 
  \sa `Kernel::ConstructVertex_3` 
  \sa `Kernel::DoIntersect_3` 
  \sa `Kernel::Equal_3` 
  \sa `Kernel::HasOn_3` 
  \sa `Kernel::Intersect_3` 
  \sa `Kernel::IsDegenerate_3` 

*/
class Segment_3 {
public:
}; /* end Kernel::Segment_3 */

/*!
  \ingroup PkgKernel23ConceptsGeomObject
  \cgalConcept

  A type representing spheres in three dimensions. 

  \cgalRefines CopyConstructible
  \cgalRefines Assignable
  \cgalRefines DefaultConstructible

  \sa `CGAL::Sphere_3<Kernel>` 
  \sa `Kernel::BoundedSide_3` 
  \sa `Kernel::ComputeSquaredRadius_3` 
  \sa `Kernel::ConstructCenter_3` 
  \sa `Kernel::ConstructOppositeSphere_3` 
  \sa `Kernel::ConstructRadicalPlane_3` 
  \sa `Kernel::ConstructSphere_3` 
  \sa `Kernel::Equal_2` 
  \sa `Kernel::HasOnBoundary_3` 
  \sa `Kernel::HasOnBoundedSide_3` 
  \sa `Kernel::HasOnNegativeSide_3` 
  \sa `Kernel::HasOnPositiveSide_3` 
  \sa `Kernel::HasOnUnboundedSide_3` 
  \sa `Kernel::IsDegenerate_3` 
  \sa `Kernel::OrientedSide_3` 

*/
class Sphere_3 {
public:
}; /* end Kernel::Sphere_3 */

/*!
  \ingroup PkgKernel23ConceptsGeomObject
  \cgalConcept

  A type representing tetrahedra in three dimensions. 

  \cgalRefines CopyConstructible
  \cgalRefines Assignable
  \cgalRefines DefaultConstructible

  \sa `CGAL::Tetrahedron_3<Kernel>` 
  \sa `Kernel::BoundedSide_3` 
  \sa `Kernel::ComputeVolume_3` 
  \sa `Kernel::ConstructCentroid_3` 
  \sa `Kernel::ConstructTetrahedron_3` 
  \sa `Kernel::ConstructVertex_3` 
  \sa `Kernel::Equal_2` 
  \sa `Kernel::HasOnBoundary_3` 
  \sa `Kernel::HasOnBoundedSide_3` 
  \sa `Kernel::HasOnNegativeSide_3` 
  \sa `Kernel::HasOnPositiveSide_3` 
  \sa `Kernel::HasOnUnboundedSide_3` 
  \sa `Kernel::IsDegenerate_3` 
  \sa `Kernel::OrientedSide_3` 

*/
class Tetrahedron_3 {
public:

}; /* end Kernel::Tetrahedron_3 */


/*!
  \ingroup PkgKernel23ConceptsGeomObject
  \cgalConcept

  A type representing triangles in two dimensions. 

  \cgalRefines CopyConstructible
  \cgalRefines Assignable
  \cgalRefines DefaultConstructible

  \sa `CGAL::Triangle_2<Kernel>` 
  \sa `Kernel::BoundedSide_2` 
  \sa `Kernel::ComputeArea_2` 
  \sa `Kernel::ComputeSquaredDistance_2` 
  \sa `Kernel::ConstructCentroid_2` 
  \sa `Kernel::ConstructOppositeTriangle_2` 
  \sa `Kernel::ConstructTriangle_2` 
  \sa `Kernel::ConstructVertex_2` 
  \sa `Kernel::DoIntersect_2` 
  \sa `Kernel::Equal_2` 
  \sa `Kernel::HasOnBoundary_2` 
  \sa `Kernel::HasOnBoundedSide_2` 
  \sa `Kernel::HasOnNegativeSide_2` 
  \sa `Kernel::HasOnPositiveSide_2` 
  \sa `Kernel::HasOnUnboundedSide_2` 
  \sa `Kernel::Intersect_2` 
  \sa `Kernel::IsDegenerate_2` 
  \sa `Kernel::OrientedSide_2` 

*/
class Triangle_2 {
public:
}; /* end Kernel::Triangle_2 */

/*!
  \ingroup PkgKernel23ConceptsGeomObject
  \cgalConcept

  A type representing triangles in three dimensions. 

  \cgalRefines CopyConstructible
  \cgalRefines Assignable
  \cgalRefines DefaultConstructible

  \sa `CGAL::Triangle_3<Kernel>` 
  \sa `Kernel::ComputeSquaredArea_3` 
  \sa `Kernel::ConstructCentroid_3` 
  \sa `Kernel::ConstructSupportingPlane_3` 
  \sa `Kernel::ConstructTriangle_3` 
  \sa `Kernel::ConstructVertex_3` 
  \sa `Kernel::DoIntersect_3` 
  \sa `Kernel::Equal_3` 
  \sa `Kernel::HasOn_3` 
  \sa `Kernel::IsDegenerate_3` 

*/
class Triangle_3 {
public:
}; /* end Kernel::Triangle_3 */

/*!
  \ingroup PkgKernel23ConceptsGeomObject
  \cgalConcept

  A type representing vectors in two dimensions. 

  \cgalRefines CopyConstructible
  \cgalRefines Assignable
  \cgalRefines DefaultConstructible

  \sa `CGAL::Vector_2<Kernel>` 
  \sa `Kernel::ComputeDeterminant_2` 
  \sa `Kernel::ComputeX_2` 
  \sa `Kernel::ComputeY_2` 
  \sa `Kernel::ComputeHx_2` 
  \sa `Kernel::ComputeHy_2` 
  \sa `Kernel::ConstructDirection_2` 
  \sa `Kernel::ConstructOppositeVector_2` 
  \sa `Kernel::ConstructPerpendicularVector_2` 
  \sa `Kernel::ConstructScaledVector_2` 
  \sa `Kernel::ConstructDividedVector_2` 
  \sa `Kernel::ConstructSumOfVectors_2` 
  \sa `Kernel::ConstructDifferenceOfVectors_2` 
  \sa `Kernel::ConstructVector_2` 
  \sa `Kernel::Equal_2` 
  \sa `Kernel::Orientation_2` 

*/
class Vector_2 {
public:
}; /* end Kernel::Vector_2 */

/*!
  \ingroup PkgKernel23ConceptsGeomObject
\cgalConcept

A type representing vectors in three dimensions. 

\cgalRefines CopyConstructible
\cgalRefines Assignable
\cgalRefines DefaultConstructible

\sa `CGAL::Vector_3<Kernel>` 
\sa `Kernel::ComputeDeterminant_3` 
\sa `Kernel::ComputeX_3` 
\sa `Kernel::ComputeY_3` 
\sa `Kernel::ComputeZ_3` 
\sa `Kernel::ConstructCrossProductVector_3` 
\sa `Kernel::ConstructDirection_3` 
\sa `Kernel::ConstructOppositeVector_3` 
\sa `Kernel::ConstructOrthogonalVector_3` 
\sa `Kernel::ConstructScaledVector_3` 
\sa `Kernel::ConstructDividedVector_3` 
\sa `Kernel::ConstructSumOfVectors_3` 
\sa `Kernel::ConstructDifferenceOfVectors_3` 
\sa `Kernel::ConstructVector_3` 
\sa `Kernel::Equal_3` 
\sa `Kernel::Orientation_3` 

*/
class Vector_3 {
public:

}; /* end Kernel::Vector_3 */

} /* end namespace Kernel */
