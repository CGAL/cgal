#ifndef TEST_ARC_TRAITS_H
#define TEST_ARC_TRAITS_H

#if TEST_TRAITS == LINE_ARC_TRAITS || \
  TEST_TRAITS == CIRCULAR_ARC_TRAITS || \
  TEST_TRAITS == CIRCULAR_LINE_ARC_TRAITS

/*! Read an arc point */
template <typename T_Traits, typename stream>
bool read_arc_point(stream & is, typename T_Traits::Point_2 & p)
{
  Basic_number_type x, y;
  is >> x >> y;
  Circular_kernel::Point_2 lp(x, y);
  p = typename T_Traits::Point_2(lp);
  // std::cout << "lp: " << lp << ", p: " << p << std::endl;
  return true;
}

#endif

#if TEST_TRAITS == LINE_ARC_TRAITS

/*! Read a line arc point */
template <>
template <class stream>
bool
Traits_test<CGAL::Arr_line_arc_traits_2<Circular_kernel> >::
read_point(stream & is,
           CGAL::Arr_line_arc_traits_2<Circular_kernel>::Point_2 & p)
{
  return
    read_arc_point<CGAL::Arr_line_arc_traits_2<Circular_kernel>, stream>(is, p);
}

/*! Read an x-monotone line arc curve */
template <>
template <class stream>
bool
Traits_test<CGAL::Arr_line_arc_traits_2<Circular_kernel> >::
read_xcurve(stream & is,
            CGAL::Arr_line_arc_traits_2<Circular_kernel>::X_monotone_curve_2 &
              xc)
{
  Basic_number_type x1, y1, x2, y2;
  is >> x1 >> y1 >> x2 >> y2;
  Circular_kernel::Point_2 lp1(x1, y1);
  Circular_kernel::Point_2 lp2(x2, y2);
  CGAL_assertion(lp1 != lp2);
  Circular_kernel::Line_arc_2 l(lp1, lp2);
  xc = CGAL::Arr_line_arc_traits_2<Circular_kernel>::X_monotone_curve_2(l);
  return true;
}

/*! Read a general line arc curve */
template <>
template <class stream>
bool
Traits_test<CGAL::Arr_line_arc_traits_2<Circular_kernel> >::
read_curve(stream & is,
           CGAL::Arr_line_arc_traits_2<Circular_kernel>::Curve_2 & cv)
{
  Basic_number_type x1, y1, x2, y2;
  is >> x1 >> y1 >> x2 >> y2;
  Circular_kernel::Point_2 lp1(x1, y1);
  Circular_kernel::Point_2 lp2(x2, y2);
  CGAL_assertion(lp1 != lp2);
  Circular_kernel::Line_arc_2 l(lp1, lp2);
  cv = CGAL::Arr_line_arc_traits_2<Circular_kernel>::Curve_2(l);
  return true;
}

#elif TEST_TRAITS == CIRCULAR_ARC_TRAITS

/*! Read a circular arc point */
template <>
template <class stream>
bool
Traits_test<CGAL::Arr_circular_arc_traits_2<Circular_kernel> >::
read_point(stream & is,
           CGAL::Arr_circular_arc_traits_2<Circular_kernel>::Point_2 & p)
{
  typedef CGAL::Arr_circular_arc_traits_2<Circular_kernel>          Traits;
  return read_arc_point<Traits, stream>(is, p);
}

/*! Read an x-monotone circular arc curve */
template <>
template <class stream>
bool
Traits_test<CGAL::Arr_circular_arc_traits_2<Circular_kernel> >::
read_xcurve(stream & is,
            CGAL::
            Arr_circular_arc_traits_2<Circular_kernel>::X_monotone_curve_2 & xc)
{
  return true;
}

/*! Read a general circular curve */
template <>
template <class stream>
bool
Traits_test<CGAL::Arr_circular_arc_traits_2<Circular_kernel> >::
read_curve(stream & is,
           CGAL::Arr_circular_arc_traits_2<Circular_kernel>::Curve_2 & cv)
{
  return true;
}

#elif TEST_TRAITS == CIRCULAR_LINE_ARC_TRAITS

/*! Read a circular-line arc point */
template <>
template <class stream>
bool
Traits_test<CGAL::Arr_circular_line_arc_traits_2<Circular_kernel> >::
read_point(stream & is,
           CGAL::Arr_circular_line_arc_traits_2<Circular_kernel>::Point_2 & p)
{
  typedef CGAL::Arr_circular_line_arc_traits_2<Circular_kernel>      Traits;
  return read_arc_point<Traits, stream>(is, p);
}

/*! Read an x-monotone circular-line arc curve */
template <>
template <class stream>
bool
Traits_test<CGAL::Arr_circular_line_arc_traits_2<Circular_kernel> >::
read_xcurve(stream & is,
            CGAL::
            Arr_circular_line_arc_traits_2<Circular_kernel>::
              X_monotone_curve_2 & xc)
{
  return true;
}

/*! Read a general circular-line curve */
template <>
template <class stream>
bool
Traits_test<CGAL::Arr_circular_line_arc_traits_2<Circular_kernel> >::
read_curve(stream & is,
           CGAL::Arr_circular_line_arc_traits_2<Circular_kernel>::Curve_2 & cv)
{
  return true;
}

#endif

#endif
