#include <CGAL/basic.h>
#include <CGAL/Arrangement_2.h>

#include <vector>

#include "test_configuration.h"
#include "test_traits.h"
#include "Traits_test.h"

// Arrangement types:
typedef CGAL::Arr_default_dcel<Traits>                  Dcel;
typedef CGAL::Arrangement_2<Traits, Dcel>               Arr;

// Traits types:
typedef Traits::Point_2                                 Point_2;
typedef Traits::X_monotone_curve_2                      X_monotone_curve_2;
typedef Traits::Curve_2                                 Curve_2;

int main (int argc, char * argv[])
{
  Traits_test<Traits> test(argc, argv);
  return (test.start()) ? 0 /* SUCCESS */ : -1; /* FAILURE */
}

#if TEST_TRAITS == POLYLINE_TRAITS || TEST_TRAITS == NON_CACHING_POLYLINE_TRAITS

template <>
template <class stream>
bool
Traits_test<CGAL::Arr_polyline_traits_2<Segment_traits> >::
read_xcurve(stream & is,
            CGAL::
            Arr_polyline_traits_2<Segment_traits>::X_monotone_curve_2 & cv)
{
  unsigned int num_points;
  is >> num_points;
  std::vector<Point_2> points;
  points.clear();
  for (unsigned int j = 0; j < num_points; j++) {
    Basic_number_type x, y;
    is >> x >> y;
    Point_2 p(x, y);
    points.push_back(p);
  }
  cv =
    CGAL::
    Arr_polyline_traits_2<Segment_traits>::X_monotone_curve_2(points.begin(),
                                                              points.end());
  return true;
}

template <>
template <class stream>
bool
Traits_test<CGAL::Arr_polyline_traits_2<Segment_traits> >::
read_curve(stream & is,
           CGAL::Arr_polyline_traits_2<Segment_traits>::Curve_2 & cv)
{
  unsigned int num_points;
  is >> num_points;
  std::vector<Point_2> points;
  points.clear();
  for (unsigned int j = 0; j < num_points; j++) {
    Basic_number_type x, y;
    is >> x >> y;
    Point_2 p(x, y);
    points.push_back(p);
  }
  cv =
    CGAL::
    Arr_polyline_traits_2<Segment_traits>::Curve_2(points.begin(),
                                                   points.end());
  return true;
}

#elif TEST_TRAITS == CORE_CONIC_TRAITS

template <>
template <class stream>
bool
Traits_test<CGAL::Arr_conic_traits_2<Rat_kernel,Alg_kernel,Nt_traits> >::
read_point(stream & is,
           CGAL::Arr_conic_traits_2<Rat_kernel,Alg_kernel,Nt_traits>::
           Point_2 & p)
{
  Basic_number_type x, y;
  is >> x >> y;
  p = CGAL::Arr_conic_traits_2<Rat_kernel,Alg_kernel,Nt_traits>::Point_2(x, y);
  return true;
}

template <>
template <class stream>
bool
Traits_test<CGAL::Arr_conic_traits_2<Rat_kernel,Alg_kernel,Nt_traits> >::
read_xcurve(stream & is,
            CGAL::Arr_conic_traits_2<Rat_kernel,Alg_kernel,Nt_traits>::
            X_monotone_curve_2 & cv)
{
  return false;
}

template <>
template <class stream>
bool
Traits_test<CGAL::Arr_conic_traits_2<Rat_kernel,Alg_kernel,Nt_traits> >::
read_curve(stream & is,
           CGAL::Arr_conic_traits_2<Rat_kernel,Alg_kernel,Nt_traits>::
           Curve_2 & cv)
{
  return false;
}

#endif
