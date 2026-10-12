#define CGAL_AW2_TIMER

#include <CGAL/alpha_wrap_2.h>
#include <CGAL/Alpha_wrap_2/internal/validation.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/bounding_box.h>
#include <CGAL/IO/WKT.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/General_polygon_with_holes_2.h>
#include <CGAL/Multipolygon_with_holes_2.h>
#include <CGAL/Random.h>

#include <fstream>
#include <iostream>
#include <vector>

using namespace CGAL::Alpha_wraps_2::internal;

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;

using FT = Kernel::FT;
using Point_2 = Kernel::Point_2;
using Segment_2 = Kernel::Segment_2;
using Triangle_2 = Kernel::Triangle_2;

using Points = std::vector<Point_2>;

using Polygon_2 = CGAL::Polygon_2<Kernel>;
using Polygon_with_holes = CGAL::Polygon_with_holes_2<Kernel>;
using Multipolygon = CGAL::Multipolygon_with_holes_2<Kernel>;

// This is just to test the API, polygons are converted into polylines anyway

template <typename K>
struct Raw_traits
{
  // Should be in some concept like TriangulationTraits_2 *at least* (bug)
  using Compute_x_2 = typename K::Compute_x_2;
  using Compute_y_2 = typename K::Compute_y_2;

  Compute_x_2 compute_x_2_object() const { return Compute_x_2(); }
  Compute_y_2 compute_y_2_object() const { return Compute_y_2(); }

  // SpatialSortingTraits_2
  using Point_2 = typename K::Point_2;

  using Compare_x_2 = typename K::Compare_x_2;
  using Compare_y_2 = typename K::Compare_y_2;
  using Less_x_2 = typename K::Less_x_2;
  using Less_y_2 = typename K::Less_y_2;

  Less_x_2 less_x_2_object() const { return Less_x_2(); }
  Less_y_2 less_y_2_object() const { return Less_y_2(); }
  Compare_x_2 compare_x_2_object() const { return Compare_x_2(); }
  Compare_y_2 compare_y_2_object() const { return Compare_y_2(); }

  // TriangulationTaits_2
  using Segment_2 = typename K::Segment_2;
  using Triangle_2 = typename K::Triangle_2;

  using Construct_point_2 = typename K::Construct_point_2;
  using Construct_segment_2 = typename K::Construct_segment_2;
  using Construct_triangle_2 = typename K::Construct_triangle_2;
  using Compare_xy_2 = typename K::Compare_xy_2;
  using Orientation_2 = typename K::Orientation_2;
  using Side_of_oriented_circle_2 = typename K::Side_of_oriented_circle_2;
  using Construct_circumcenter_2 = typename K::Construct_circumcenter_2;

  Construct_point_2 construct_point_2_object() const { return Construct_point_2(); }
  Construct_segment_2 construct_segment_2_object() const { return Construct_segment_2(); }
  Construct_triangle_2 construct_triangle_2_object() const { return Construct_triangle_2(); }
  Compare_xy_2 compare_xy_2_object() const { return Compare_xy_2(); }
  Orientation_2 orientation_2_object() const { return Orientation_2(); }
  Side_of_oriented_circle_2 side_of_oriented_circle_2_object() const { return Side_of_oriented_circle_2(); }
  Construct_circumcenter_2 construct_circumcenter_2_object() const { return Construct_circumcenter_2(); }

  // DelaunayTriangulationTraits_2
  using Line_2 = typename K::Line_2;
  using Ray_2 = typename K::Ray_2;

  using Compare_distance_2 = typename K::Compare_distance_2;
  using Construct_bisector_2 = typename K::Construct_bisector_2;
  using Construct_ray_2 = typename K::Construct_ray_2;

  Compare_distance_2 compare_distance_2_object() const { return Compare_distance_2(); }
  Construct_bisector_2 construct_bisector_2_object() const { return Construct_bisector_2(); }
  Construct_ray_2 construct_ray_2_object() const { return Construct_ray_2(); }

  // SearchGeomTraits_2
  using Iso_rectangle_2 = typename K::Iso_rectangle_2;
  using Circle_2 = typename K::Circle_2;

  using Construct_min_vertex_2 = typename K::Construct_min_vertex_2;
  using Construct_max_vertex_2 = typename K::Construct_max_vertex_2;
  using Construct_center_2 = typename K::Construct_center_2;
  using Compute_squared_radius_2 = typename K::Compute_squared_radius_2;
  using Construct_iso_rectangle_2 = typename K::Construct_iso_rectangle_2;
  using Cartesian_const_iterator_2 = typename K::Cartesian_const_iterator_2;
  using Construct_cartesian_const_iterator_2 = typename K::Construct_cartesian_const_iterator_2;

    // no _object() for some reason...

  // AABBGeomTraits_2
  using Boolean = typename K::Boolean;

  using Do_intersect_2 = typename K::Do_intersect_2;
  using Intersect_2 = typename K::Intersect_2;
  using Construct_circle_2 = typename K::Construct_circle_2;
  using Construct_projected_point_2 = typename K::Construct_projected_point_2;
  using Compute_squared_distance_2 = typename K::Compute_squared_distance_2;
  using Equal_2 = typename K::Equal_2;

  Do_intersect_2 do_intersect_2_object() const { return Do_intersect_2(); }
  Intersect_2 intersect_2_object() const { return Intersect_2(); }
  Construct_circle_2 construct_circle_2_object() const { return Construct_circle_2(); }
  Construct_projected_point_2 construct_projected_point_2_object() const { return Construct_projected_point_2(); }
  Compute_squared_distance_2 compute_squared_distance_2_object() const { return Compute_squared_distance_2(); }
  Equal_2 equal_2_object() const { return Equal_2(); }

  // AABBRayIntersectionGeomTraits_2
  using Vector_2 = typename K::Vector_2;

  using Construct_vector_2 = typename K::Construct_vector_2;
  using Construct_source_2 = typename K::Construct_source_2;

  Construct_vector_2 construct_vector_2_object() const { return Construct_vector_2(); }
  Construct_source_2 construct_source_2_object() const { return Construct_source_2(); }

  // PolygonTraits_2
  using FT = typename K::FT;

  using Less_xy_2 = typename K::Less_xy_2;
  using Less_yx_2 = typename K::Less_yx_2;
  using Compute_area_2 = typename K::Compute_area_2;

  Less_xy_2 less_xy_2_object() const { return Less_xy_2(); };
  Less_yx_2 less_yx_2_object() const { return Less_yx_2(); };
  Compute_area_2 compute_area_2_object() const { return Compute_area_2(); };

  // AlphaWrapTraits_2

  using Angle_2 = typename K::Angle_2;
  using Construct_bbox_2 = typename K::Construct_bbox_2;
  using Construct_scaled_vector_2 = typename K::Construct_scaled_vector_2;
  using Construct_translated_point_2 = typename K::Construct_translated_point_2;
  using Construct_vertex_2 = typename K::Construct_vertex_2;
  using Has_on_bounded_side_2 = typename K::Has_on_bounded_side_2;
  using Has_on_unbounded_side_2 = typename K::Has_on_unbounded_side_2;
  using Is_degenerate_2 = typename K::Is_degenerate_2;
  using Side_of_bounded_circle_2 = typename K::Side_of_bounded_circle_2;

  Angle_2 angle_2_object() const { return Angle_2(); }
  Construct_bbox_2 construct_bbox_2_object() const { return Construct_bbox_2(); }
  Construct_scaled_vector_2 construct_scaled_vector_2_object() const { return Construct_scaled_vector_2(); }
  Construct_translated_point_2 construct_translated_point_2_object() const { return Construct_translated_point_2(); }
  Construct_vertex_2 construct_vertex_2_object() const { return Construct_vertex_2(); }
  Has_on_bounded_side_2 has_on_bounded_side_2_object() const { return Has_on_bounded_side_2(); }
  Has_on_unbounded_side_2 has_on_unbounded_side_2_object() const { return Has_on_unbounded_side_2(); }
  Is_degenerate_2 is_degenerate_2_object() const { return Is_degenerate_2(); }
  Side_of_bounded_circle_2 side_of_bounded_circle_2_object() const { return Side_of_bounded_circle_2(); }
};

int main(int, char**)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  Multipolygon wrap;

  Polygon_2 poly;
  CGAL::alpha_wrap_2(poly, 1, wrap);
  CGAL::alpha_wrap_2(poly, 1, wrap, CGAL::parameters::geom_traits(Raw_traits<Kernel>()));
  CGAL::alpha_wrap_2(poly, 1, 2, wrap);
  CGAL::alpha_wrap_2(poly, 1, 2, wrap, CGAL::parameters::geom_traits(Raw_traits<Kernel>()));

  Polygon_with_holes pwh;
  CGAL::alpha_wrap_2(pwh, 1, wrap);
  CGAL::alpha_wrap_2(pwh, 1, wrap, CGAL::parameters::geom_traits(Raw_traits<Kernel>()));
  CGAL::alpha_wrap_2(pwh, 1, 2, wrap);
  CGAL::alpha_wrap_2(pwh, 1, 2, wrap, CGAL::parameters::geom_traits(Raw_traits<Kernel>()));

  Polygon_with_holes mp;
  CGAL::alpha_wrap_2(mp, 1, wrap);
  CGAL::alpha_wrap_2(mp, 1, wrap, CGAL::parameters::geom_traits(Raw_traits<Kernel>()));
  CGAL::alpha_wrap_2(mp, 1, 2, wrap);
  CGAL::alpha_wrap_2(mp, 1, 2, wrap, CGAL::parameters::geom_traits(Raw_traits<Kernel>()));

  Points pts;
  CGAL::alpha_wrap_2(pts, 1, wrap);
  CGAL::alpha_wrap_2(pts, 1, wrap, CGAL::parameters::geom_traits(Raw_traits<Kernel>()));
  CGAL::alpha_wrap_2(pts, 1, 2, wrap);
  CGAL::alpha_wrap_2(pts, 1, 2, wrap, CGAL::parameters::geom_traits(Raw_traits<Kernel>()));

  std::deque<Segment_2> sss;
  CGAL::alpha_wrap_2(sss, 1, wrap);
  CGAL::alpha_wrap_2(sss, 1, wrap, CGAL::parameters::geom_traits(Raw_traits<Kernel>()));
  CGAL::alpha_wrap_2(sss, 1, 2, wrap);
  CGAL::alpha_wrap_2(sss, 1, 2, wrap, CGAL::parameters::geom_traits(Raw_traits<Kernel>()));

  std::vector<std::deque<Point_2> > mls;
  CGAL::alpha_wrap_2(mls, 1, wrap);
  CGAL::alpha_wrap_2(mls, 1, wrap, CGAL::parameters::geom_traits(Raw_traits<Kernel>()));
  CGAL::alpha_wrap_2(mls, 1, 2, wrap);
  CGAL::alpha_wrap_2(mls, 1, 2, wrap, CGAL::parameters::geom_traits(Raw_traits<Kernel>()));

  std::list<Triangle_2> trs;
  CGAL::alpha_wrap_2(trs, 1, wrap);
  CGAL::alpha_wrap_2(trs, 1, wrap, CGAL::parameters::geom_traits(Raw_traits<Kernel>()));
  CGAL::alpha_wrap_2(trs, 1, 2, wrap);
  CGAL::alpha_wrap_2(trs, 1, 2, wrap, CGAL::parameters::geom_traits(Raw_traits<Kernel>()));

  std::vector<Point_2> pos;
  std::vector<std::array<std::size_t, 3> > faces;
  CGAL::alpha_wrap_2(pos, faces, 1, wrap);
  CGAL::alpha_wrap_2(pos, faces, 1, wrap, CGAL::parameters::geom_traits(Raw_traits<Kernel>()));
  CGAL::alpha_wrap_2(pos, faces, 1, 2, wrap);
  CGAL::alpha_wrap_2(pos, faces, 1, 2, wrap, CGAL::parameters::geom_traits(Raw_traits<Kernel>()));

  std::cout << "Done!" << std::endl;
  return EXIT_SUCCESS;
}
