#include <CGAL/Segment_2.h>
#include <CGAL/Vector_2.h>

template <class R>
inline
bool CGAL::lexicographically_xy_smaller(const MyPoint<R>& p, const MyPoint<R>& q)
{
  return CGAL::lexicographically_xy_smaller((const CGAL::Point_2<R>&) p,
                                           (const CGAL::Point_2<R>&) q );
}

template <class R>
inline
bool CGAL::lexicographically_yx_smaller(const MyPoint<R>& p, const MyPoint<R>& q)
{
  return CGAL::lexicographically_yx_smaller((const CGAL::Point_2<R>&) p,
                                           (const CGAL::Point_2<R>&) q );
}

template <class R>
inline
bool CGAL::lexicographically_yx_smaller_or_equal(const MyPoint<R>& p, const MyPoint<R>& q)
{
  return CGAL::lexicographically_yx_smaller_or_equal((const CGAL::Point_2<R>&) p,
                                                    (const CGAL::Point_2<R>&) q );
}

template <class R>
CGAL::Comparison_result CGAL::compare_x(const MyPoint<R>& p, const MyPoint<R>& q)
{
  return CGAL::compare_x((const CGAL::Point_2<R>&) p,
                        (const CGAL::Point_2<R>&) q );
}

template <class R>
CGAL::Comparison_result CGAL::compare_y(const MyPoint<R>& p, const MyPoint<R>& q)
{
  return CGAL::compare_y((const CGAL::Point_2<R>&) p,
                        (const CGAL::Point_2<R>&) q );
}

template <class R>
ostream& operator<<(ostream& to, const MyPoint<R>& p)
{
  return to << (const CGAL::Point_2<R>&) p;
}

template <class R>
CGAL::Vector_2<R> operator-(const MyPoint<R>& p, const MyPoint<R>& q)
{
  return ((const CGAL::Point_2<R>&) p) - ((const CGAL::Point_2<R>&) q );
}

template < class R >
CGAL::Orientation CGAL::orientation(const MyPoint<R>& p,
                                  const MyPoint<R>& q,
                                  const MyPoint<R>& r)
{
  return CGAL::orientation((const CGAL::Point_2<R>&) p,
                          (const CGAL::Point_2<R>&) q,
                          (const CGAL::Point_2<R>&) r );
}

