#ifndef CGAL_TESTSUITE_APPROX_EQUAL_H
#define CGAL_TESTSUITE_APPROX_EQUAL_H
#include <boost/math/special_functions/next.hpp>

namespace CGAL {
namespace testsuite {

template <typename FT>
bool approx_equal(FT a, FT b) { return a == b; }

bool approx_equal(double a, double b) {
  return std::abs(boost::math::float_distance(a, b)) <= 1;
}

struct Xyz_tag {};
struct Xy_tag {};

template <typename Object>
bool approx_equal(Object a, Object b, CGAL::testsuite::Xyz_tag)
{
  return approx_equal(a.x(), b.x()) &&
    approx_equal(a.y(), b.y()) &&
    approx_equal(a.z(), b.z());
}

template <typename Object>
bool approx_equal(Object a, Object b, CGAL::testsuite::Xy_tag)
{
  return approx_equal(a.x(), b.x()) && approx_equal(a.y(), b.y());
}

struct Direction_2_tag {};
template <typename Object>
bool approx_equal(Object a, Object b, CGAL::testsuite::Direction_2_tag)
{
  return CGAL_NTS sign(a.dx()) == CGAL_NTS sign(a.dx())
    && CGAL_NTS sign(a.dy()) == CGAL_NTS sign(b.dy())
    && approx_equal(a.dx() * b.dy(), a.dy() * b.dx());
}

} // end namespace testsuite
} // end namespace CGAL

#endif // CGAL_TESTSUITE_APPROX_EQUAL_H
