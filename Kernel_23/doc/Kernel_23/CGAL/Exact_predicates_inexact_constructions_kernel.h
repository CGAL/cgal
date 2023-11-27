
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

<dl class="section"><dt>Trivial Constructions</dt></dl>
<p>Some geometric constructions, however, are exact because they only
copy parameters and do not involve any computations that could lead to
roundoff errors. We call such a construction <em>trivial</em>. For
instance, all copy constructors, such as `Point_2::Point_2(const
Point_2&)`, are trivial constructions. In addition, the following
constructions in `CGAL::Exact_predicates_inexact_constructions_kernel`
are trivial:</p>

<ul>
<li> `Point_2::Point_2(FT, FT)`
<li> `Point_2::Point_2(Weighted_point_2)`
<li> `Point_2::Point_2(Origin)`
<li> `FT Point_2::x()`
<li> `FT Point_2::y()`
<li> `FT Point_2::cartesian(int)`
<li> `FT Point_2::operator[](int)`
<li> `Bbox_2 Point_2::bbox()`
<li> `Vector_2 operator-(Point_2, Origin)`
<li> `Vector_2 operator-(Origin, Point_2)`
<li> `Weighted_point_2::Weighted_point_2(FT, FT)`
<li> `Weighted_point_2::Weighted_point_2(Point_2, FT)`
<li> `Weighted_point_2::Weighted_point_2(Point_2)`
<li> `Weighted_point_2::Weighted_point_2(Origin)`
<li> `Point_2 Weighted_point_2::point()`
<li> `FT Weighted_point_2::weight()`
<li> `FT Weighted_point_2::x()`
<li> `FT Weighted_point_2::y()`
<li> `FT Weighted_point_2::cartesian(int)`
<li> `FT Weighted_point_2::operator[](int)`
<li> `Bbox_2 Weighted_point_2::bbox()`
<li> `Vector_2::Vector_2(FT, FT)`
<li> `FT Vector_2::x()`
<li> `FT Vector_2::y()`
<li> `FT Vector_2::cartesian(int)`
<li> `FT Vector_2::operator[](int)`
<li> `Direction_2 Vector_2::direction()`
<li> `Vector_2 Vector_2::operator-()`
<li> `Direction_2::Direction_2(FT, FT)`
<li> `Direction_2::Direction_2(Vector_2)`
<li> `RT Direction_2::delta(int)`
<li> `RT Direction_2::dx()`
<li> `RT Direction_2::dy()`
<li> `Direction_2 Direction_2::operator-()`
<li> `Vector_2 Direction_2::vector()`
<li> `Line_2::Line_2(RT, RT, RT)`
<li> `RT Line_2::a()`
<li> `RT Line_2::b()`
<li> `Line_2 Line_2::opposite()`
<li> `Ray_2::Ray_2(Point_2, Point_2)`
<li> `Point_2 Ray_2::source()`
<li> `Segment_2::Segment_2(Point_2, Point_2)`
<li> `Point_2 Segment_2::source()`
<li> `Point_2 Segment_2::target()`
<li> `Point_2 Segment_2::min()`
<li> `Point_2 Segment_2::max()`
<li> `Point_2 Segment_2::vertex(int)`
<li> `Point_2 Segment_2::point(int)`
<li> `Point_2 Segment_2::operator[](int)`
<li> `Segment_2 Segment_2::opposite()`
<li> `Bbox_2 Segment_2::bbox()`
<li> `Triangle_2::Triangle_2(Point_2, Point_2, Point_2)`
<li> `Point_2 Triangle_2::vertex(int)`
<li> `Point_2 Triangle_2::operator[](int)`
<li> `Triangle_2 Triangle_2::opposite()`
<li> `Bbox_2 Triangle_2::bbox()`
<li> `Iso_rectangle_2::Iso_rectangle_2(Point_2, Point_2)`
<li> `Iso_rectangle_2::Iso_rectangle_2(Point_2, Point_2, int)`
<li> `Iso_rectangle_2::Iso_rectangle_2(Point_2, Point_2, Point_2, Point_2)`
<li> `Iso_rectangle_2::Iso_rectangle_2(Bbox_2)`
<li> `Point_2 Iso_rectangle_2::vertex(int)`
<li> `Point_2 Iso_rectangle_2::operator[](int)`
<li> `Point_2 Iso_rectangle_2::min()`
<li> `Point_2 Iso_rectangle_2::max()`
<li> `FT Iso_rectangle_2::xmin()`
<li> `FT Iso_rectangle_2::ymin()`
<li> `FT Iso_rectangle_2::xmax()`
<li> `FT Iso_rectangle_2::ymax()`
<li> `FT Iso_rectangle_2::min_coord(int)`
<li> `FT Iso_rectangle_2::max_coord(int)`
<li> `Bbox_2 Iso_rectangle_2::bbox()`
<li> `Circle_2::Circle_2(Point_2, FT, Orientation)`
<li> `Circle_2::Circle_2(Point_2, Orientation)`
<li> `Point_2 Circle_2::center()`
<li> `FT Circle_2::squared_radius()`
<li> `Orientation Circle_2::orientation()`
<li> `Circle_2 Circle_2::opposite()`
<li> `Point_3::Point_3(FT, FT, FT)`
<li> `Point_3::Point_3(Weighted_point_3)`
<li> `Point_3::Point_3(Origin)`
<li> `FT Point_3::x()`
<li> `FT Point_3::y()`
<li> `FT Point_3::z()`
<li> `FT Point_3::cartesian(int)`
<li> `FT Point_3::operator[](int)`
<li> `Bbox_3 Point_3::bbox()`
<li> `Vector_3 operator-(Point_3, Origin)`
<li> `Vector_3 operator-(Origin, Point_3)`
<li> `Weighted_point_3::Weighted_point_3(FT, FT, FT)`
<li> `Weighted_point_3::Weighted_point_3(Point_3, FT)`
<li> `Weighted_point_3::Weighted_point_3(Point_3)`
<li> `Weighted_point_3::Weighted_point_3(Origin)`
<li> `Point_3 Weighted_point_3::point()`
<li> `FT Weighted_point_3::weight()`
<li> `FT Weighted_point_3::x()`
<li> `FT Weighted_point_3::y()`
<li> `FT Weighted_point_3::z()`
<li> `FT Weighted_point_3::cartesian(int)`
<li> `FT Weighted_point_3::operator[](int)`
<li> `Bbox_3 Weighted_point_3::bbox()`
<li> `Vector_3::Vector_3(FT, FT, FT)`
<li> `Vector_3::Vector_3(Null_vector)`
<li> `FT Vector_3::x()`
<li> `FT Vector_3::y()`
<li> `FT Vector_3::z()`
<li> `FT Vector_3::cartesian(int)`
<li> `FT Vector_3::operator[](int)`
<li> `Direction_3 Vector_3::direction()`
<li> `Vector_3 Vector_3::opposite()`
<li> `Direction_3::Direction_3(RT, RT, RT)`
<li> `Direction_3::Direction_3(Vector_3)`
<li> `RT Direction_3::delta(int)`
<li> `RT Direction_3::dx()`
<li> `RT Direction_3::dy()`
<li> `RT Direction_3::dz()`
<li> `Direction_3 Direction_3::opposite()`
<li> `Vector_3 Direction_3::vector()`
<li> `Line_3::Line_3(Point_3, Vector_3)`
<li> `Line_3::Line_3(Point_3, Direction_3)`
<li> `Line_3 Line_3::opposite()`
<li> `Vector_3 Line_3::to_vector()`
<li> `Direction_3 Line_3::direction()`
<li> `Plane_3::Plane_3(FT, FT, FT, FT)`
<li> `Plane_3::Plane_3(Circle_3)`
<li> `FT Plane_3::a()`
<li> `FT Plane_3::b()`
<li> `FT Plane_3::c()`
<li> `FT Plane_3::d()`
<li> `Plane_3 Plane_3::opposite()`
<li> `Vector_3 Plane_3::orthogonal_vector()`
<li> `Direction_3 Plane_3::orthogonal_direction()`
<li> `Circle_3::Circle_3(Point_3, FT, Plane_3)`
<li> `Point_3 Circle_3::center()`
<li> `FT Circle_3::squared_radius()`
<li> `FT Circle_3::supporting_plane()`
<li> `Sphere_3 Circle_3::diametral_sphere()`
<li> `Sphere_3::Sphere_3(Point_3, FT, Orientation)`
<li> `Sphere_3::Sphere_3(Point_3, Orientation)`
<li> `Sphere_3::Sphere_3(Circle_3)`
<li> `Point_3 Sphere_3::center()`
<li> `FT Sphere_3::squared_radius()`
<li> `Orientation Sphere_3::orientation()`
<li> `Sphere_3 Sphere_3::opposite()`
<li> `Ray_3::Ray_3(Point_3, Point_3)`
<li> `Point_3 Ray_3::source()`
<li> `Segment_3::Segment_3(Point_3, Point_3)`
<li> `Point_3 Segment_3::source()`
<li> `Point_3 Segment_3::target()`
<li> `Point_3 Segment_3::min()`
<li> `Point_3 Segment_3::max()`
<li> `Point_3 Segment_3::vertex(int)`
<li> `Point_3 Segment_3::point(int)`
<li> `Point_3 Segment_3::operator[](int)`
<li> `Segment_3 Segment_3::opposite()`
<li> `Bbox_3 Segment_3::bbox()`
<li> `Triangle_3::Triangle_3(Point_3, Point_3, Point_3)`
<li> `Point_3 Triangle_3::vertex(int)`
<li> `Point_3 Triangle_3::operator[](int)`
<li> `Bbox_3 Triangle_3::bbox()`
<li> `Tetrahedron_3::Tetrahedron_3(Point_3, Point_3, Point_3, Point_3)`
<li> `Point_3 Tetrahedron_3::vertex(int)`
<li> `Point_3 Tetrahedron_3::operator[](int)`
<li> `Bbox_3 Tetrahedron_3::bbox()`
<li> `Iso_cuboid_3::Iso_cuboid_3(Point_3, Point_3)`
<li> `Iso_cuboid_3::Iso_cuboid_3(Point_3, Point_3, int)`
<li> `Iso_cuboid_3::Iso_cuboid_3(Point_3, Point_3, Point_3, Point_3, Point_3, Point_3)`
<li> `Iso_cuboid_3::Iso_cuboid_3(Bbox_3)`
<li> `Point_3 Iso_cuboid_3::vertex(int)`
<li> `Point_3 Iso_cuboid_3::operator[](int)`
<li> `Point_3 Iso_cuboid_3::min()`
<li> `Point_3 Iso_cuboid_3::max()`
<li> `FT Iso_cuboid_3::xmin()`
<li> `FT Iso_cuboid_3::ymin()`
<li> `FT Iso_cuboid_3::zmin()`
<li> `FT Iso_cuboid_3::xmax()`
<li> `FT Iso_cuboid_3::ymax()`
<li> `FT Iso_cuboid_3::zmax()`
<li> `FT Iso_cuboid_3::min_coord(int)`
<li> `FT Iso_cuboid_3::max_coord(int)`
<li> `Bbox_3 Iso_cuboid_3::bbox()`
<li> `Point_2 min_vertex(Iso_rectangle_2)`
<li> `Point_2 min_vertex(Iso_cuboid_3)`
<li> `Point_2 max_vertex(Iso_rectangle_2)`
<li> `Point_2 max_vertex(Iso_cuboid_3)`
</ul>

\cgalModels{Kernel}

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
