
namespace CGAL {

/*!
\ingroup kernel_predef

A typedef to a kernel that has the following properties:

<UL> 
<LI>It uses %Cartesian representation. 
<LI>It supports constructions of points from `double` %Cartesian 
coordinates. 
<LI>It provides exact geometric predicates, but geometric 
constructions are not exact in general. 
</UL>

<h4>Trivial Constructions</h4>  

Some geometric constructions, however, are exact because they only
copy parameters and do not involve any computations that could lead to
roundoff errors. We call such a construction <em>trivial</em>. For
instance, all copy constructors, such as `Point_2::Point_2(const
Point_2&)`, are trivial constructions. In addition, the following
constructions in `CGAL::Exact_predicates_inexact_constructions_kernel`
are trivial:

<ul> 
<li> `Point_2::Point_2(double,double)`
<li> `Point_2::Point_2(const Weighted_point_2&)`
<li> `Point_2::Point_2(Origin)`
<li> `Weighted_point_2::Weighted_point_2(double,double)`
<li> `Weighted_point_2::Weighted_point_2(const Point_2&, double = 0)`
<li> `Weighted_point_2::Weighted_point_2(Origin)`
<li> `Vector_2::Vector_2(double,double)`
<li> `Direction_2::Direction_2(double,double)`
<li> `Direction_2::Direction_2(const Vector_2&)`
<li> `Line_2::Line_2(double,double,double)`
<li> `Ray_2::Ray_2(const Point_2&, const Point_2&)`
<li> `Segment_2::Segment_2(const Point_2&, const Point_2&)`
<li> `Triangle_2::Triangle_2(const Point_2&, const Point_2&, const Point_2&)`
<li> `Iso_rectangle_2::Iso_rectangle_2(const Point_2&, const Point_2&)`
<li> `Iso_rectangle_2::Iso_rectangle_2(const Point_2&, const Point_2&, int)`
<li> `Iso_rectangle_2::Iso_rectangle_2(const Point_2&, const Point_2&, const Point_2&, const Point_2&)`
<li> `Iso_rectangle_2::Iso_rectangle_2(const Bbox_2&)`
<li> `Circle_2::Circle_2(const Point_2&, double, Orientation = COUNTERCLOCKWISE)`
<li> `Kernel::ConstructPoint_2()(const Weighted_point_2&)`
<li> `Kernel::ConstructPoint_2()(Origin)`
<li> `Kernel::ConstructWeightedPoint_2()(const Point_2&, double = 0)`
<li> `Kernel::ConstructWeightedPoint_2()(Origin)`
<li> `Kernel::ConstructVector_2()(Origin, const Point_2&)`
<li> `Kernel::ConstructVector_2()(const Point_2&, Origin)`
<li> `Kernel::ConstructVector_2()(Null_vector)`
<li> `Kernel::ConstructDirection_2()(const Vector_2&)`
<li> `Kernel::ConstructSegment_2()(const Point_2&, const Point_2&)`
<li> `Kernel::ConstructRay_2()(const Point_2&, const Point_2&)`
<li> `Kernel::ConstructCircle_2()(const Point_2&, double, Orientation = COUNTERCLOCKWISE)`
<li> `Kernel::ConstructCircle_2()(const Point_2&, Orientation = COUNTERCLOCKWISE)`
<li> `Kernel::ConstructTriangle_2()(const Point_2&, const Point_2&, const Point_2&)`
<li> `Kernel::ConstructIsoRectangle_2()(const Point_2&, const Point_2&)`
<li> `Kernel::ConstructIsoRectangle_2()(const Point_2&, const Point_2&, int)`
<li> `Kernel::ConstructIsoRectangle_2()(const Point_2&, const Point_2&, const Point_2&, const Point_2&)`
<li> `Kernel::ConstructVertex_2()(const Segment_2&, int)`
<li> `Kernel::ConstructVertex_2()(const Triangle_2&, int)`
<li> `Kernel::ConstructVertex_2()(const Iso_rectangle_2&, int)`
<li> `Kernel::ConstructBbox_2()(const Point_2&)`
<li> `Kernel::ConstructBbox_2()(const Segment_2&)`
<li> `Kernel::ConstructBbox_2()(const Triangle_2&)`
<li> `Kernel::ConstructBbox_2()(const Iso_rectangle_2&)`
<li> `Kernel::ConstructCenter_2()(const Circle_2&)`
<li> `Kernel::ConstructOppositeVector_2()(const Vector_2&)`
<li> `Kernel::ConstructOppositeDirection_2()(const Direction_2&)`
<li> `Kernel::ConstructOppositeSegment_2()(const Segment_2&)`
<li> `Kernel::ConstructOppositeRay_2()(const Ray_2&)`
<li> `Kernel::ConstructOppositeLine_2()(const Line_2&)`
<li> `Kernel::ConstructOppositeTriangle_2()(const Triangle_2&)`
<li> `Kernel::ConstructOppositeCircle_2()(const Circle_2&)`
<li> `Kernel::ConstructSource_2()(const Segment_2&)`
<li> `Kernel::ConstructSource_2()(const Ray_2&)`
<li> `Kernel::ConstructSecondPoint_2()(const Ray_2&)`
<li> `Kernel::ConstructTarget_2()(const Target_2&)`
<li> `Kernel::ConstructMinVertex_2()(const Segment_2&)`
<li> `Kernel::ConstructMinVertex_2()(const Iso_rectangle_2&)`
<li> `Kernel::ConstructMaxVertex_2()(const Segment_2&)`
<li> `Kernel::ConstructMaxVertex_2()(const Iso_rectangle_2&)`
<li> `Kernel::ComputeWeight_2()(const Weighted_point_2&)`
<li> `Kernel::ComputeX_2()(const Point_2&)`
<li> `Kernel::ComputeX_2()(const Vector_2&)`
<li> `Kernel::ComputeY_2()(const Point_2&)`
<li> `Kernel::ComputeY_2()(const Vector_2&)`
<li> `Kernel::ComputeA_2()(const Line_2&)`
<li> `Kernel::ComputeB_2()(const Line_2&)`
<li> `Kernel::ComputeC_2()(const Line_2&)`
<li> `Kernel::ComputeDx_2()(const Direction_2&)`
<li> `Kernel::ComputeDy_2()(const Direction_2&)`
<li> `Kernel::ComputeXmin_2()(const Iso_rectangle_2&)`
<li> `Kernel::ComputeXmax_2()(const Iso_rectangle_2&)`
<li> `Kernel::ComputeYmin_2()(const Iso_rectangle_2&)`
<li> `Kernel::ComputeYmax_2()(const Iso_rectangle_2&)`
<li> `Kernel::ComputeSquaredRadius_2()(const Circle_2&)`
<li> `Kernel::ComputeSquaredRadius_2()(const Point_2&)`
<li> `Kernel::CartesianConstIterator_2()(const Point_2&)`
<li> `Kernel::CartesianConstIterator_2()(const Point_2&, int)`
<li> `Kernel::CartesianConstIterator_2()(const Vector_2&)`
<li> `Kernel::CartesianConstIterator_2()(const Vector_2&, int)`
<li> `Point_3::Point_3(double,double,double)`
<li> `Point_3::Point_3(const Weighted_point_3&)`
<li> `Point_3::Point_3(Origin)`
<li> `Weighted_point_3::Weighted_point_3(double,double)`
<li> `Weighted_point_3::Weighted_point_3(const Point_3&, double = 0)`
<li> `Weighted_point_3::Weighted_point_3(Origin)`
<li> `Vector_3::Vector_3(double,double,double)`
<li> `Direction_3::Direction_3(double,double,double)`
<li> `Direction_3::Direction_3(const Vector_3&)`
<li> `Line_3::Line_3(const Point_3&, const Vector_3&)`
<li> `Line_3::Line_3(const Point_3&, const Direction_3&)`
<li> `Plane_3::Plane_3(double,double,double,double)`
<li> `Circle_3::Circle_3(const Point_3&, double, const Plane_3&)`
<li> `Sphere_3::Sphere_3(const Point_3&, double, Orientation = COUNTERCLOCKWISE)`
<li> `Ray_3::Ray_3(const Point_3&, const Point_3&)`
<li> `Segment_3::Segment_3(const Point_3&, const Point_3&)`
<li> `Triangle_3::Triangle_3(const Point_3&, const Point_3&, const Point_3&)`
<li> `Tetrahedron_3::Tetrahedron_3(const Point_3&, const Point_3&, const Point_3&, const Point_3&)`
<li> `Iso_cuboid_3::Iso_cuboid_3(const Point_3&, const Point_3&)`
<li> `Iso_cuboid_3::Iso_cuboid_3(const Point_3&, const Point_3&, int)`
<li> `Iso_cuboid_3::Iso_cuboid_3(const Point_3&, const Point_3&, const Point_3&, const Point_3&, const Point_3&, const Point_3&)`
<li> `Iso_cuboid_3::Iso_cuboid_3(const Bbox_3&)`
<li> `Kernel::ConstructPoint_3()(const Weighted_point_3&)`
<li> `Kernel::ConstructPoint_3()(Origin)`
<li> `Kernel::ConstructWeightedPoint_3()(const Point_3&, double = 0)`
<li> `Kernel::ConstructWeightedPoint_3()(Origin)`
<li> `Kernel::ConstructVector_3()(Origin, const Point_3&)`
<li> `Kernel::ConstructVector_3()(const Point_3&, Origin)`
<li> `Kernel::ConstructVector_3()(Null_vector)`
<li> `Kernel::ConstructDirection_3()(const Vector_3&)`
<li> `Kernel::ConstructPlane_3()(double,double,double,double)`
<li> `Kernel::ConstructIsoCuboid_3()(const Point_3&, const Point_3&)`
<li> `Kernel::ConstructIsoCuboid_3()(const Point_3&, const Point_3&, int)`
<li> `Kernel::ConstructIsoCuboid_3()(const Point_3&, const Point_3&, const Point_3&, const Point_3&, const Point_3&, const Point_3&)`
<li> `Kernel::ConstructLine_3()(const Point_3&, const Vector_3&)`
<li> `Kernel::ConstructLine_3()(const Point_3&, const Direction_3&)`
<li> `Kernel::ConstructRay_3()(const Point_3&, const Point_3&)`
<li> `Kernel::ConstructSphere_3()(const Point_3&, double, Orientation = COUNTERCLOCKWISE)`
<li> `Kernel::ConstructSphere_3()(const Point_3&, Orientation = COUNTERCLOCKWISE)`
<li> `Kernel::ConstructSegment_3()(const Point_3&, const Point_3&)`
<li> `Kernel::ConstructTriangle_3()(const Point_3&, const Point_3&, const Point_3&)`
<li> `Kernel::ConstructTetrahedron_3()(const Point_3&, const Point_3&, const Point_3&, const Point_3&)`
<li> `Kernel::ConstructCircle_3()(const Point_3&, double, const Plane_3&)`
<li> `Kernel::ConstructCircle_3()(
<li> `Kernel::ConstructVertex_3()(const Segment_3&, int)`
<li> `Kernel::ConstructVertex_3()(const Triangle_3&, int)`
<li> `Kernel::ConstructVertex_3()(const Tetrahedron_3&, int)`
<li> `Kernel::ConstructVertex_3()(const Iso_cuboid_3&, int)`
<li> `Kernel::ConstructSource_3()(const Segment_3&)`
<li> `Kernel::ConstructSource_3()(const Ray_3&)`
<li> `Kernel::ConstructTarget_3()(const Target_3&)`
<li> `Kernel::ConstructMinVertex_3()(const Segment_3&)`
<li> `Kernel::ConstructMinVertex_3()(const Iso_cuboid_3&)`
<li> `Kernel::ConstructMaxVertex_3()(const Segment_3&)`
<li> `Kernel::ConstructMaxVertex_3()(const Iso_cuboid_3&)`
<li> `Kernel::ConstructBbox_3()(const Point_3&)`
<li> `Kernel::ConstructBbox_3()(const Segment_3&)`
<li> `Kernel::ConstructBbox_3()(const Triangle_3&)`
<li> `Kernel::ConstructBbox_3()(const Tetrahedron_3&)`
<li> `Kernel::ConstructBbox_3()(const Iso_cuboid_3&)`
<li> `Kernel::ConstructCenter_3()(const Sphere_3&)`
<li> `Kernel::ConstructCenter_3()(const Circle_3&)`
<li> `Kernel::ConstructSecondPoint_3()(const Ray_3&)`
<li> `Kernel::ConstructOppositeVector_3()(const Vector_3&)`
<li> `Kernel::ConstructOppositeDirection_3()(const Direction_3&)`
<li> `Kernel::ConstructOppositeSegment_3()(const Segment_3&)`
<li> `Kernel::ConstructOppositeRay_3()(const Ray_3&)`
<li> `Kernel::ConstructOppositeLine_3()(const Line_3&)`
<li> `Kernel::ConstructOppositeTriangle_3()(const Triangle_3&)`
<li> `Kernel::ConstructOppositePlane_3()(const Plane_3&)`
<li> `Kernel::ConstructOppositeSphere_3()(const Sphere_3&)`
<li> `Kernel::ConstructOppositeCircle_3()(const Circle_3&)`
<li> `Kernel::ComputeWeight_3()(const Weighted_point_3&)`
<li> `Kernel::ComputeX_3()(const Point_3&)`
<li> `Kernel::ComputeX_3()(const Vector_3&)`
<li> `Kernel::ComputeY_3()(const Point_3&)`
<li> `Kernel::ComputeY_3()(const Vector_3&)`
<li> `Kernel::ComputeZ_3()(const Point_3&)`
<li> `Kernel::ComputeZ_3()(const Vector_3&)`
<li> `Kernel::ComputeA_3()(const Line_3&)`
<li> `Kernel::ComputeB_3()(const Line_3&)`
<li> `Kernel::ComputeC_3()(const Line_3&)`
<li> `Kernel::ComputeD_3()(const Line_3&)`
<li> `Kernel::ComputeDx_3()(const Direction_3&)`
<li> `Kernel::ComputeDy_3()(const Direction_3&)`
<li> `Kernel::ComputeDz_3()(const Direction_3&)`
<li> `Kernel::ComputeXmin_3()(const Iso_cuboid_3&)`
<li> `Kernel::ComputeXmax_3()(const Iso_cuboid_3&)`
<li> `Kernel::ComputeYmin_3()(const Iso_cuboid_3&)`
<li> `Kernel::ComputeYmax_3()(const Iso_cuboid_3&)`
<li> `Kernel::ComputeZmin_3()(const Iso_cuboid_3&)`
<li> `Kernel::ComputeZmax_3()(const Iso_cuboid_3&)`
<li> `Kernel::ComputeSquaredRadius_3()(const Sphere_3&)`
<li> `Kernel::ComputeSquaredRadius_3()(const Circle_3&)`
<li> `Kernel::ComputeSquaredRadius_3()(const Point_3&)`
<li> `Kernel::CartesianConstIterator_3()(const Point_3&)`
<li> `Kernel::CartesianConstIterator_3()(const Point_3&, int)`
<li> `Kernel::CartesianConstIterator_3()(const Vector_3&)`
<li> `Kernel::CartesianConstIterator_3()(const Vector_3&, int)`
</ul>

\cgalModels `Kernel`

\sa `CGAL::Exact_predicates_exact_constructions_kernel`
\sa `CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt`
\sa `CGAL::Exact_predicates_exact_constructions_kernel_with_kth_root`
\sa `CGAL::Exact_predicates_exact_constructions_kernel_with_root_of`
\sa `CGAL::Cartesian`

*/

class Exact_predicates_inexact_constructions_kernel {
public:

}; /* end Exact_predicates_inexact_constructions_kernel */
} /* end namespace CGAL */
