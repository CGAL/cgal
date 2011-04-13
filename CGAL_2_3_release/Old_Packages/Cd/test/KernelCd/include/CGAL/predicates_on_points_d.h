// revision      : 
// revision_date : 
// author(s)     : Hervé Brönnimann

#ifndef CGAL_PREDICATES_ON_POINTS_D_H
#define CGAL_PREDICATES_ON_POINTS_D_H

#ifndef CGAL_REP_CLASS_DEFINED
#error  no representation class defined
#endif  // CGAL_REP_CLASS_DEFINED

#ifdef CGAL_HOMOGENEOUS_H
#include <CGAL/predicates_on_pointsHd.h>
#endif // CGAL_HOMOGENEOUS_H

#ifdef CGAL_CARTESIAN_H
#include <CGAL/Cartesian/predicates_on_points_d.h>
#endif // CGAL_CARTESIAN_H

#ifndef CGAL_POINT_D_H
#include <CGAL/Point_d.h>
#endif // CGAL_POINT_D_H

CGAL_BEGIN_NAMESPACE

template < class R >
inline
Comparison_result
compare_lexicographically_d( const Point_d<R> &p,
                             const Point_d<R> &q)
{
  typedef typename  R::Point_d_base  RPoint_d;
  return compare_lexicographically_d((const RPoint_d& )p,
                                     (const RPoint_d& )q);
}

template < class R >
inline
bool
lexicographically_d_smaller_or_equal(const Point_d<R> &p,
                                       const Point_d<R> &q)
{
  typedef typename  R::Point_d_base  RPoint_d;
  return ( !( compare_lexicographically_d((const RPoint_d& )p,
                                          (const RPoint_d& )q)
              == LARGER ) );
}

template < class R >
inline
bool
lexicographically_d_smaller(const Point_d<R> &p,
                            const Point_d<R> &q)
{
  typedef typename  R::Point_d_base  RPoint_d;
  return ( compare_lexicographically_d((const RPoint_d& )p,
                                       (const RPoint_d& )q)
           == SMALLER );
}

template < class R >
inline
bool
x_equal(const Point_d<R> &p,
        const Point_d<R> &q,
	int i = 0)
{
  typedef typename  R::Point_d_base  RPoint_d;
  return x_equal((const RPoint_d& )p, (const RPoint_d& )q, i);
}

template < class R >
inline
Comparison_result
compare_x(const Point_d<R> &p, const Point_d<R> &q, int i = 0)
{
  typedef typename  R::Point_d_base  RPoint_d;
  return compare_x((const RPoint_d& )p, (const RPoint_d& )q, i);
}

template < class InputIterator >
inline
bool
coplanar(const InputIterator &first, const InputIterator &last)
{
  typedef typename std::iterator_traits<InputIterator>::value_type Point_d;
  typedef typename Point_d::R::Rep_tag Rep_tag;
  return coplanar(first, last, Rep_tag());
}

template < class R >
inline
bool
collinear(const Point_d<R> &p,
          const Point_d<R> &q,
          const Point_d<R> &r)
{
  typedef typename  R::Point_d_base  RPoint_d;
  return (collinear((const RPoint_d& )p,
                    (const RPoint_d& )q,
                    (const RPoint_d& )r));
}

template < class R >
inline
bool
are_ordered_along_line(const Point_d<R> &p,
                       const Point_d<R> &q,
                       const Point_d<R> &r)
{
  typedef typename  R::Point_d_base  RPoint_d;
  return (are_ordered_along_line((const RPoint_d& )p,
                                 (const RPoint_d& )q,
                                 (const RPoint_d& )r));
}

template < class R >
inline
bool
collinear_are_ordered_along_line(const Point_d<R> &p,
                                 const Point_d<R> &q,
                                 const Point_d<R> &r)
{
  typedef typename  R::Point_d_base  RPoint_d;
  return collinear_are_ordered_along_line((const RPoint_d& )p,
                                          (const RPoint_d& )q,
                                          (const RPoint_d& )r
                                              );
}

template < class R >
inline
bool
are_strictly_ordered_along_line(const Point_d<R> &p,
                                const Point_d<R> &q,
                                const Point_d<R> &r)
{
  typedef typename  R::Point_d_base  RPoint_d;
  return (are_strictly_ordered_along_line((const RPoint_d& )p,
                                          (const RPoint_d& )q,
                                          (const RPoint_d& )r));
}

template < class R >
inline
bool
collinear_are_strictly_ordered_along_line(const Point_d<R> &p,
                                          const Point_d<R> &q,
                                          const Point_d<R> &r)
{
  typedef typename  R::Point_d_base  RPoint_d;
  return collinear_are_strictly_ordered_along_line((const RPoint_d& )p,
                                                   (const RPoint_d& )q,
                                                   (const RPoint_d& )r
                                              );
}

template < class InputIterator >
inline
Orientation
orientation(const InputIterator &first, const InputIterator &last)
{
  typedef typename std::iterator_traits<InputIterator>::value_type Point_d;
  typedef typename Point_d::R::Rep_tag Rep_tag;
  return orientation(first, last, Rep_tag());
}

template <class R, class PointIterator >
inline
Bounded_side
side_of_bounded_sphere( const PointIterator &first,
                        const PointIterator &last,
                        const Point_d<R> &test)
{
  typedef typename  R::Point_d_base  RPoint_d;
  return side_of_bounded_sphere(first, last,
                                (const RPoint_d& )test);
}

template <class R, class PointIterator >
inline
Oriented_side
side_of_oriented_sphere( const PointIterator &first,
                         const PointIterator &last,
                         const Point_d<R> &test)
{
  typedef typename  R::Point_d_base  RPoint_d;
  return side_of_oriented_sphere(first, last,
                                 (const RPoint_d& )test);
}

CGAL_END_NAMESPACE

#endif // CGAL_PREDICATES_ON_POINTS_d_H
