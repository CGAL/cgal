#include <CGAL/Segment_2.h>
#include <CGAL/Vector_2.h>

template <class R>
inline
bool CGAL_lexicographically_xy_smaller(const MyPoint<R>& p, const MyPoint<R>& q)
{
  return CGAL_lexicographically_xy_smaller((const CGAL_Point_2<R>&) p,
                                           (const CGAL_Point_2<R>&) q );
}

template <class R>
inline
bool CGAL_lexicographically_yx_smaller(const MyPoint<R>& p, const MyPoint<R>& q)
{
  return CGAL_lexicographically_yx_smaller((const CGAL_Point_2<R>&) p,
                                           (const CGAL_Point_2<R>&) q );
}

template <class R>
inline
bool CGAL_lexicographically_yx_smaller_or_equal(const MyPoint<R>& p, const MyPoint<R>& q)
{
  return CGAL_lexicographically_yx_smaller_or_equal((const CGAL_Point_2<R>&) p,
                                                    (const CGAL_Point_2<R>&) q );
}

template <class R>
CGAL_Comparison_result CGAL_compare_x(const MyPoint<R>& p, const MyPoint<R>& q)
{
  return CGAL_compare_x((const CGAL_Point_2<R>&) p,
                        (const CGAL_Point_2<R>&) q );
}

template <class R>
CGAL_Comparison_result CGAL_compare_y(const MyPoint<R>& p, const MyPoint<R>& q)
{
  return CGAL_compare_y((const CGAL_Point_2<R>&) p,
                        (const CGAL_Point_2<R>&) q );
}

template <class R>
ostream& operator<<(ostream& to, const MyPoint<R>& p)
{
  return to << (const CGAL_Point_2<R>&) p;
}

template <class R>
CGAL_Vector_2<R> operator-(const MyPoint<R>& p, const MyPoint<R>& q)
{
  return ((const CGAL_Point_2<R>&) p) - ((const CGAL_Point_2<R>&) q );
}

template < class R >
CGAL_Orientation CGAL_orientation(const MyPoint<R>& p,
                                  const MyPoint<R>& q,
                                  const MyPoint<R>& r)
{
  return CGAL_orientation((const CGAL_Point_2<R>&) p,
                          (const CGAL_Point_2<R>&) q,
                          (const CGAL_Point_2<R>&) r );
}

