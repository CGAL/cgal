// Regression test for issue #7124:
// Bbox_3 intersection with Segment_3/Ray_3/Line_3 lost precision due to
// intermediate to_double() conversions in the old intersection_bl() helper.
// Also tests degenerate segment handling.

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Intersections_3/Bbox_3_Segment_3.h>
#include <CGAL/Intersections_3/Bbox_3_Ray_3.h>
#include <CGAL/Intersections_3/Bbox_3_Line_3.h>

#include <iostream>
#include <cassert>

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef K::FT FT;
typedef K::Point_3 Point;
typedef K::Segment_3 Segment;
typedef K::Ray_3 Ray;
typedef K::Line_3 Line;

int main()
{
  // Original reproducer from issue #7124:
  // A segment with large coordinates that should intersect a small bbox.
  // src.y = 1/25 is inside the box; trg.y is far outside.
  {
    std::cout << "Test 1: Original issue #7124 reproducer (Segment_3)" << std::endl;
    CGAL::Bbox_3 box(0, 0, 0, 10, 10, 10);
    Point p1(FT(3) + FT(1)/FT(2475), FT(1)/FT(25), FT(1) + FT(1)/FT(2013));
    Point p2(FT(3) + FT(1)/FT(2475), FT(-65734357381440816LL), FT(1) + FT(1)/FT(2013));
    Segment seg(p1, p2);
    auto result = CGAL::intersection(seg, box);
    assert(result);
    const Segment* s = std::get_if<Segment>(&*result);
    assert(s != nullptr);
    std::cout << "  OK: intersection is a segment" << std::endl;
  }

  // Segment fully inside bbox
  {
    std::cout << "Test 2: Segment fully inside bbox" << std::endl;
    CGAL::Bbox_3 box(0, 0, 0, 10, 10, 10);
    Segment seg(Point(2, 3, 4), Point(5, 6, 7));
    auto result = CGAL::intersection(seg, box);
    assert(result);
    const Segment* s = std::get_if<Segment>(&*result);
    assert(s != nullptr);
    std::cout << "  OK" << std::endl;
  }

  // Segment crossing bbox boundary
  {
    std::cout << "Test 3: Segment crossing bbox boundary" << std::endl;
    CGAL::Bbox_3 box(0, 0, 0, 10, 10, 10);
    Segment seg(Point(5, 5, 5), Point(15, 5, 5));
    auto result = CGAL::intersection(seg, box);
    assert(result);
    const Segment* s = std::get_if<Segment>(&*result);
    assert(s != nullptr);
    assert(s->source() == Point(5, 5, 5));
    assert(s->target() == Point(10, 5, 5));
    std::cout << "  OK" << std::endl;
  }

  // Segment missing bbox entirely
  {
    std::cout << "Test 4: Segment missing bbox" << std::endl;
    CGAL::Bbox_3 box(0, 0, 0, 10, 10, 10);
    Segment seg(Point(20, 20, 20), Point(30, 30, 30));
    auto result = CGAL::intersection(seg, box);
    assert(!result);
    std::cout << "  OK" << std::endl;
  }

  // Degenerate segment (source == target) inside bbox
  {
    std::cout << "Test 5: Degenerate segment inside bbox" << std::endl;
    CGAL::Bbox_3 box(0, 0, 0, 10, 10, 10);
    Point p(5, 5, 5);
    Segment seg(p, p);
    auto result = CGAL::intersection(seg, box);
    assert(result);
    const Point* pt = std::get_if<Point>(&*result);
    assert(pt != nullptr);
    assert(*pt == p);
    std::cout << "  OK" << std::endl;
  }

  // Degenerate segment (source == target) outside bbox
  {
    std::cout << "Test 6: Degenerate segment outside bbox" << std::endl;
    CGAL::Bbox_3 box(0, 0, 0, 10, 10, 10);
    Point p(20, 20, 20);
    Segment seg(p, p);
    auto result = CGAL::intersection(seg, box);
    assert(!result);
    std::cout << "  OK" << std::endl;
  }

  // Degenerate segment on bbox boundary
  {
    std::cout << "Test 7: Degenerate segment on bbox boundary" << std::endl;
    CGAL::Bbox_3 box(0, 0, 0, 10, 10, 10);
    Point p(10, 5, 5);
    Segment seg(p, p);
    auto result = CGAL::intersection(seg, box);
    assert(result);
    const Point* pt = std::get_if<Point>(&*result);
    assert(pt != nullptr);
    assert(*pt == p);
    std::cout << "  OK" << std::endl;
  }

  // Ray through bbox with large coordinates
  {
    std::cout << "Test 8: Ray through bbox (large coords)" << std::endl;
    CGAL::Bbox_3 box(0, 0, 0, 10, 10, 10);
    Point origin(FT(3) + FT(1)/FT(2475), FT(-65734357381440816LL), FT(1) + FT(1)/FT(2013));
    Point target(FT(3) + FT(1)/FT(2475), FT(0), FT(1) + FT(1)/FT(2013));
    Ray ray(origin, target);
    auto result = CGAL::intersection(ray, box);
    assert(result);
    std::cout << "  OK" << std::endl;
  }

  // Line through bbox with large coordinates
  {
    std::cout << "Test 9: Line through bbox (large coords)" << std::endl;
    CGAL::Bbox_3 box(0, 0, 0, 10, 10, 10);
    Point p1(FT(3) + FT(1)/FT(2475), FT(-65734357381440816LL), FT(1) + FT(1)/FT(2013));
    Point p2(FT(3) + FT(1)/FT(2475), FT(0), FT(1) + FT(1)/FT(2013));
    Line line(p1, p2);
    auto result = CGAL::intersection(line, box);
    assert(result);
    std::cout << "  OK" << std::endl;
  }

  // Lower-precision kernel: Simple_cartesian<float>. The bbox coordinates are
  // double; computing in the coercion type (double here) preserves them,
  // instead of degrading the bbox to float as an Iso_cuboid<float> would
  // (the precision-loss case raised in review).
  {
    std::cout << "Test 10: Simple_cartesian<float> kernel" << std::endl;
    typedef CGAL::Simple_cartesian<float> Kf;
    typedef Kf::Point_3 Pf;
    typedef Kf::Segment_3 Sf;
    typedef Kf::Ray_3 Rf;
    typedef Kf::Line_3 Lf;
    CGAL::Bbox_3 box(0, 0, 0, 10, 10, 10);

    // Segment crossing the boundary: exact clip in float.
    Sf seg(Pf(5, 5, 5), Pf(15, 5, 5));
    auto rseg = CGAL::intersection(seg, box);
    assert(rseg);
    const Sf* s = std::get_if<Sf>(&*rseg);
    assert(s != nullptr);
    assert(s->source() == Pf(5, 5, 5));
    assert(s->target() == Pf(10, 5, 5));

    // Ray and line through the box.
    Rf ray(Pf(5, 5, 5), Pf(15, 5, 5));
    assert(CGAL::intersection(ray, box));
    Lf line(Pf(-5, 5, 5), Pf(15, 5, 5));
    assert(CGAL::intersection(line, box));

    // Degenerate segment inside the box returns a point.
    Sf dseg(Pf(3, 3, 3), Pf(3, 3, 3));
    auto rd = CGAL::intersection(dseg, box);
    assert(rd);
    assert(std::get_if<Pf>(&*rd) != nullptr);

    std::cout << "  OK" << std::endl;
  }

  std::cout << "All tests passed." << std::endl;
  return 0;
}
