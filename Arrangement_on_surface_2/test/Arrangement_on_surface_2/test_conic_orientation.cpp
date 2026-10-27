#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Cartesian.h>
#include <CGAL/CORE_algebraic_number_traits.h>

#include <array>
#include <cassert>

using Nt_traits = CGAL::CORE_algebraic_number_traits;
using Rational = Nt_traits::Rational;
using Rat_kernel = CGAL::Cartesian<Rational>;
using Alg_kernel = CGAL::Cartesian<Nt_traits::Algebraic>;
using Traits = CGAL::Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>;

Rational evaluate(const Traits::Curve_2& curve, const Rational& x, const Rational& y)
{
  return curve.r()*x*x + curve.s()*y*y + curve.t()*x*y +
         curve.u()*x + curve.v()*y + curve.w();
}

void test_ellipse(const Rat_kernel::Point_2& center,
                  const Rational& x_radius, const Rational& y_radius)
{
  Traits traits;
  auto construct = traits.construct_curve_2_object();
  for(auto orientation : {CGAL::CLOCKWISE, CGAL::COUNTERCLOCKWISE})
  {
    auto curve = construct(center, x_radius, y_radius, orientation);
    assert(curve.is_valid());
    assert(curve.is_full_conic());
    assert(curve.orientation() == orientation);
    assert(curve.t() == 0);
    assert(curve.r() != 0);
    assert(evaluate(curve, center.x() + x_radius, center.y()) == 0);
    assert(evaluate(curve, center.x() - x_radius, center.y()) == 0);
    assert(evaluate(curve, center.x(), center.y() + y_radius) == 0);
    assert(evaluate(curve, center.x(), center.y() - y_radius) == 0);
    assert(evaluate(curve, center.x(), center.y()) != 0);
    assert(evaluate(curve, center.x() + 2*x_radius, center.y()) != 0);
    assert(evaluate(curve, center.x(), center.y() + 2*y_radius) != 0);

    CGAL::Arrangement_2<Traits> arrangement(&traits);
    CGAL::insert(arrangement, curve);
    assert(arrangement.is_valid());
    assert(arrangement.number_of_vertices() == 2);
    assert(arrangement.number_of_edges() == 2);
    assert(arrangement.number_of_faces() == 2);
  }
}

void test_coefficient_orientation(const std::array<Rational, 6>& c)
{
  Traits traits;
  auto construct = traits.construct_curve_2_object();
  auto clockwise = construct(c[0], c[1], c[2], c[3], c[4], c[5]);
  auto counterclockwise = construct(-c[0], -c[1], -c[2], -c[3], -c[4], -c[5]);
  assert(clockwise.is_valid());
  assert(counterclockwise.is_valid());
  assert(clockwise.is_full_conic());
  assert(counterclockwise.is_full_conic());
  assert(clockwise.orientation() == CGAL::CLOCKWISE);
  assert(counterclockwise.orientation() == CGAL::COUNTERCLOCKWISE);
}

int main()
{
  test_ellipse({0, 0}, 1, 1);
  test_ellipse({2, -3}, 3, 2);
  test_ellipse({-4, 5}, 2, 3);
  test_ellipse({Rational(1, 3), Rational(-2, 5)}, Rational(3, 2), Rational(2, 3));
  test_ellipse({0, 0}, Rational(1, 1000), 1000);

  test_coefficient_orientation({1, 1, 0, 0, 0, -1});
  test_coefficient_orientation({2, 3, 1, -5, 4, -10});
  test_coefficient_orientation({Rational(1, 2), Rational(1, 3), Rational(1, 4),
                                Rational(-1, 5), Rational(1, 6), Rational(-1, 7)});
}
