#include <CGAL/Arr_circle_segment_traits_2.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Gmpq.h>

#include <cassert>
#include <variant>
#include <vector>

void check_bounds(const CGAL::Bbox_2& box, double xmin, double ymin, double xmax, double ymax)
{
  assert(box.xmin() <= xmin && box.ymin() <= ymin);
  assert(box.xmax() >= xmax && box.ymax() >= ymax);
  assert(xmin - box.xmin() < 1e-12 && ymin - box.ymin() < 1e-12);
  assert(box.xmax() - xmax < 1e-12 && box.ymax() - ymax < 1e-12);
}

template <bool Filter>
void check_bbox()
{
  using Kernel = CGAL::Cartesian<CGAL::Gmpq>;
  using Traits = CGAL::Arr_circle_segment_traits_2<Kernel, Filter>;
  using Point = typename Traits::Point_2;
  using Curve = typename Traits::Curve_2;
  using X_curve = typename Traits::X_monotone_curve_2;
  using Coordinate = typename Traits::CoordNT;
  const Traits traits;
  const auto bbox = traits.construct_bbox_2_object();

  check_bounds(bbox(Point(2, 3)), 2, 3, 2, 3);
  check_bounds(bbox(Curve(Kernel::Point_2(-2, 3), Kernel::Point_2(4, -1))), -2, -1, 4, 3);
  check_bounds(bbox(Curve(Kernel::Point_2(2, 3), CGAL::Gmpq(5))), -3, -2, 7, 8);
  check_bounds(bbox(Curve(Kernel::Point_2(2, 3), CGAL::Gmpq(0))), 2, 3, 2, 3);

  const Kernel::Circle_2 ccw(Kernel::Point_2(0, 0), 25, CGAL::COUNTERCLOCKWISE);
  const Kernel::Circle_2 cw(Kernel::Point_2(0, 0), 25, CGAL::CLOCKWISE);
  check_bounds(bbox(Curve(ccw, Point(5, 0), Point(-5, 0))), -5, 0, 5, 5);
  check_bounds(bbox(Curve(cw, Point(5, 0), Point(-5, 0))), -5, -5, 5, 0);
  check_bounds(bbox(Curve(ccw, Point(5, 0), Point(3, 4))), 3, 0, 5, 4);
  check_bounds(bbox(Curve(cw, Point(3, 4), Point(5, 0))), 3, 0, 5, 4);
  check_bounds(bbox(Curve(ccw, Point(3, 4), Point(3, -4))), -5, -5, 3, 5);

  const Coordinate root(CGAL::Gmpq(0), CGAL::Gmpq(1), CGAL::Gmpq(2));
  const Kernel::Circle_2 circle(Kernel::Point_2(0, 0), 2);
  const Curve arc(circle, Point(root, 0), Point(0, root));
  std::vector<std::variant<Point, X_curve>> pieces;
  traits.make_x_monotone_2_object()(arc, std::back_inserter(pieces));
  assert(pieces.size() == 1);
  const auto& piece = std::get<X_curve>(pieces.front());
  for(const auto& box : {bbox(arc), bbox(piece), bbox(piece.construct_opposite())})
    for(const auto& point : {piece.source(), piece.target(), Point(1, 1)})
    {
      assert(Coordinate(CGAL::Gmpq(box.xmin())) <= point.x());
      assert(Coordinate(CGAL::Gmpq(box.ymin())) <= point.y());
      assert(Coordinate(CGAL::Gmpq(box.xmax())) >= point.x());
      assert(Coordinate(CGAL::Gmpq(box.ymax())) >= point.y());
    }
}

int main()
{
  check_bbox<true>();
  check_bbox<false>();
}
