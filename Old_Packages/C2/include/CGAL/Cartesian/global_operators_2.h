#ifndef CGAL_CARTESIAN_GLOBAL_OPERATORS_2_C
#define CGAL_CARTESIAN_GLOBAL_OPERATORS_2_C

#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_2_H
#include <CGAL/Cartesian/redefine_names_2.h>
#endif

CGAL_BEGIN_NAMESPACE

template < class R >
inline
PointC2<R CGAL_CTAG>
operator+(const PointC2<R CGAL_CTAG> &p, const VectorC2<R CGAL_CTAG> &v)
{
  return PointC2<R CGAL_CTAG>(p.x() + v.x(), p.y() + v.y()) ;
}

template < class R >
inline
PointC2<R CGAL_CTAG>
operator-(const PointC2<R CGAL_CTAG> &p, const VectorC2<R CGAL_CTAG> &v)
{
  return PointC2<R CGAL_CTAG>(p.x() - v.x(), p.y() - v.y()) ;
}

template < class R >
inline
PointC2<R CGAL_CTAG>
operator+(const Origin &, const VectorC2<R CGAL_CTAG> &v)
{
  return PointC2<R CGAL_CTAG>(v);
}

template < class R >
inline
PointC2<R CGAL_CTAG>
operator-(const Origin &, const VectorC2<R CGAL_CTAG> &v)
{
  return PointC2<R CGAL_CTAG>(-v);
}

template < class R >
inline
VectorC2<R CGAL_CTAG>
operator-(const PointC2<R CGAL_CTAG> &p, const PointC2<R CGAL_CTAG> &q)
{
  return VectorC2<R CGAL_CTAG>(p.x() - q.x(), p.y() - q.y()) ;
}

template < class R >
inline
VectorC2<R CGAL_CTAG>
operator-(const PointC2<R CGAL_CTAG> &p, const Origin &)
{
  return VectorC2<R CGAL_CTAG>(p) ;
}

template < class R >
inline
VectorC2<R CGAL_CTAG>
operator-(const Origin &, const PointC2<R CGAL_CTAG> &p)
{
  return VectorC2<R CGAL_CTAG>(-p.x(), -p.y()) ;
}

template < class R >
CGAL_KERNEL_INLINE
VectorC2<R CGAL_CTAG>
operator*(const typename R::FT &c, const VectorC2<R CGAL_CTAG> &w)
{
   return VectorC2<R CGAL_CTAG>( c* w.x(), c * w.y()) ;
}

template < class R >
CGAL_KERNEL_INLINE
VectorC2<R CGAL_CTAG>
operator*(const VectorC2<R CGAL_CTAG> &w, const typename R::FT &c)
{
   return VectorC2<R CGAL_CTAG>( c* w.x(), c * w.y()) ;
}

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_GLOBAL_OPERATORS_2_C
