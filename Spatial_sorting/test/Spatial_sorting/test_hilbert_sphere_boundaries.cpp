#include <CGAL/Simple_cartesian.h>
#include <CGAL/hilbert_sort_on_sphere.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <type_traits>
#include <vector>

using K = CGAL::Simple_cartesian<double>;
using Point = K::Point_3;

template <class Policy>
void check_face_order(Policy, double radius, const Point& center)
{
  const double t = std::sqrt(1.0 / 3.0) * radius;
  const auto point = [&](double x, double y, double z) {
    return Point(center.x() + x, center.y() + y, center.z() + z);
  };
  std::array<std::vector<Point>, 6> faces;
  faces[0].push_back(point(radius, 0, 0));
  faces[1].push_back(point(0, radius, 0));
  faces[2].push_back(point(-radius, 0, 0));
  faces[3].push_back(point(0, 0, radius));
  faces[4].push_back(point(0, -radius, 0));
  faces[5].push_back(point(0, 0, -radius));
  if(center == Point(0, 0, 0)) {
    for(int sy : {-1, 1}) {
      for(int sz : {-1, 1}) {
        faces[0].push_back(point(t, sy * t, sz * t));
        faces[sy > 0 ? 1 : 2].push_back(point(-t, sy * t, sz * t));
      }
    }
    const double below = std::nextafter(t, 0.0);
    faces[0].push_back(point(below, below, below));
    faces[0].push_back(point(below, below, -below));
    faces[2].push_back(point(-below, below, below));
    faces[2].push_back(point(-below, below, -below));
    const double smaller = std::nextafter(below, 0.0);
    faces[1].push_back(point(smaller, below, smaller));
    faces[3].push_back(point(smaller, smaller, below));
    faces[4].push_back(point(smaller, -below, smaller));
    faces[5].push_back(point(smaller, smaller, -below));
  }
  std::vector<Point> points, expected;
  for(const auto& face : faces)
    expected.insert(expected.end(), face.begin(), face.end());
  for(auto face = faces.rbegin(); face != faces.rend(); ++face)
    points.insert(points.end(), face->begin(), face->end());
  CGAL::Hilbert_sort_on_sphere_3<K, Policy> sort(K(), radius * radius,
                                               center, points.size());
  sort(points.begin(), points.end());
  assert(points == expected);
}

template <class Policy>
void check_public_sort(Policy policy, double radius, const Point& center)
{
  const double t = std::sqrt(1.0 / 3.0) * radius;
  std::vector<Point> points;
  for(int i = 0; i < 128; ++i)
    points.emplace_back(center.x() + t, center.y() + t, center.z() + (i % 2 ? t : -t));
  const auto original = points;
  CGAL::hilbert_sort_on_sphere(points.begin(), points.end(), policy, radius * radius, center);
  assert(std::is_permutation(points.begin(), points.end(), original.begin()));
  if constexpr(std::is_same<Policy, CGAL::Hilbert_sort_middle_policy>::value) {
    unsigned transitions = 0;
    for(std::size_t i = 0; i < points.size(); ++i)
      transitions += points[i] != points[(i + 1) % points.size()];
    assert(transitions == 2);
  }
  points.clear();
  CGAL::hilbert_sort_on_sphere(points.begin(), points.end(), policy, radius * radius, center);
  const Point singleton(center.x() + radius, center.y(), center.z());
  points.push_back(singleton);
  CGAL::hilbert_sort_on_sphere(points.begin(), points.end(), policy, radius * radius, center);
  assert(points.front() == singleton);
}

template <class Policy>
void check(Policy policy)
{
  for(double radius : {0.125, 1.0, 8.0}) {
    check_face_order(policy, radius, Point(0, 0, 0));
    check_face_order(policy, radius, Point(4, -2, 1));
    check_public_sort(policy, radius, Point(0, 0, 0));
    check_public_sort(policy, radius, Point(4, -2, 1));
  }
}

int main()
{
  check(CGAL::Hilbert_sort_median_policy());
  check(CGAL::Hilbert_sort_middle_policy());
  std::cout << "Spherical boundary regression passed.\n";
}
