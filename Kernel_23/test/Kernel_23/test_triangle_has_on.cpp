#undef NDEBUG

#include <CGAL/Cartesian.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Simple_homogeneous.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Gmpz.h>
#include <CGAL/assertions_behaviour.h>
#include <CGAL/exceptions.h>
#include <CGAL/intersections.h>

#include <cassert>

template<typename Function>
void check_precondition(Function function)
{
  bool rejected = false;
  try {
    function();
  }
  catch(const CGAL::Precondition_exception&) {
    rejected = true;
  }
  assert(rejected);
}

template<typename K>
void test()
{
  using Point = typename K::Point_3;
  using Triangle = typename K::Triangle_3;
  const Point a(0, 0, 0), b(2, 0, 0), c(0, 2, 0);
  const Point edge(1, 0, 0), inside(1, 1, 0), outside(3, 0, 0), off_plane(1, 1, 1);
  auto has_on = K().has_on_3_object();

  for(const Triangle& triangle : {Triangle(a, b, c), Triangle(a, c, b)})
  {
    assert(triangle.has_on(a));
    assert(triangle.has_on(edge));
    assert(triangle.has_on(inside));
    assert(!triangle.has_on(outside));
    assert(!triangle.has_on(off_plane));
    assert(has_on(triangle, inside));
    assert(!has_on(triangle, off_plane));
    assert(CGAL::do_intersect(triangle, inside));
    assert(!CGAL::do_intersect(off_plane, triangle));
  }

  for(const Triangle& triangle : {Triangle(a, a, b), Triangle(a, b, Point(4, 0, 0)), Triangle(a, a, a)})
  {
    assert(triangle.is_degenerate());
    check_precondition([&] { triangle.has_on(a); });
    check_precondition([&] { has_on(triangle, off_plane); });
    check_precondition([&] { CGAL::do_intersect(triangle, a); });
    check_precondition([&] { CGAL::do_intersect(a, triangle); });
  }
}

int main()
{
  CGAL::set_error_behaviour(CGAL::THROW_EXCEPTION);
  test<CGAL::Cartesian<CGAL::Gmpq>>();
  test<CGAL::Simple_cartesian<CGAL::Gmpq>>();
  test<CGAL::Homogeneous<CGAL::Gmpz>>();
  test<CGAL::Simple_homogeneous<CGAL::Gmpz>>();
  test<CGAL::Exact_predicates_inexact_constructions_kernel>();
  test<CGAL::Exact_predicates_exact_constructions_kernel>();
}
