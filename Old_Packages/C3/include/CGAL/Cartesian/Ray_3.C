// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri

#ifndef CGAL_CARTESIAN_RAY_3_C
#define CGAL_CARTESIAN_RAY_3_C

#include <CGAL/Cartesian/distance_computations_3.h>

#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_3_H
#define CGAL_CTAG
#endif

#ifdef CGAL_CFG_TYPENAME_BUG
#define typename
#endif

CGAL_BEGIN_NAMESPACE

template < class R >
_Twotuple< typename RayC3<R CGAL_CTAG>::Point_3 > *
RayC3<R CGAL_CTAG>::ptr() const
{
  return (_Twotuple< Point_3 >*)PTR;
}

template < class R >
RayC3<R CGAL_CTAG>::RayC3()
{
  PTR = new _Twotuple< Point_3 >;
}

template < class R >
RayC3<R CGAL_CTAG>::
RayC3(const RayC3<R CGAL_CTAG>  &r)
  : Handle((Handle&)r)
{}

template < class R >
RayC3<R CGAL_CTAG>::
RayC3(const typename RayC3<R CGAL_CTAG>::Point_3 &sp,
      const typename RayC3<R CGAL_CTAG>::Point_3 &secondp)
{
  PTR = new _Twotuple< Point_3 >(sp,secondp);
}

template < class R >
RayC3<R CGAL_CTAG>::
RayC3(const typename RayC3<R CGAL_CTAG>::Point_3 &sp,
      const typename RayC3<R CGAL_CTAG>::Direction_3 &d)
{
  PTR = new _Twotuple< Point_3 >(sp, sp + d.to_vector());
}

template < class R >
inline RayC3<R CGAL_CTAG>::~RayC3()
{}

template < class R >
RayC3<R CGAL_CTAG> &
RayC3<R CGAL_CTAG>::operator=(const RayC3<R CGAL_CTAG> &r)
{
  Handle::operator=(r);
  return *this;
}

template < class R >
inline bool RayC3<R CGAL_CTAG>::operator==(const RayC3<R CGAL_CTAG> &r) const
{
  return (source() == r.source()) && (direction() == r.direction());
}

template < class R >
inline bool RayC3<R CGAL_CTAG>::operator!=(const RayC3<R CGAL_CTAG> &r) const
{
  return !(*this == r);
}

template < class R >
inline
long RayC3<R CGAL_CTAG>::id() const
{
  return (long) PTR;
}

template < class R >
inline
typename RayC3<R CGAL_CTAG>::Point_3
RayC3<R CGAL_CTAG>::start() const
{
  return ptr()->e0;
}

template < class R >
inline
typename RayC3<R CGAL_CTAG>::Point_3
RayC3<R CGAL_CTAG>::source() const
{
  return ptr()->e0;
}

template < class R >
inline
typename RayC3<R CGAL_CTAG>::Point_3
RayC3<R CGAL_CTAG>::second_point() const
{
  return ptr()->e1;
}

template < class R >
CGAL_KERNEL_INLINE
typename RayC3<R CGAL_CTAG>::Point_3
RayC3<R CGAL_CTAG>::point(int i) const
{
  CGAL_kernel_precondition( i >= 0 );
  if (i == 0)
    return ptr()->e0;

  if (i == 1)
    return ptr()->e1;

  return source() + FT(i) * (second_point() - source());
}

template < class R >
inline
typename RayC3<R CGAL_CTAG>::Direction_3
RayC3<R CGAL_CTAG>::direction() const
{
  return Direction_3( second_point() - source() );
}

template < class R >
inline
typename RayC3<R CGAL_CTAG>::Line_3
RayC3<R CGAL_CTAG>::supporting_line() const
{
  return Line_3(*this);
}

template < class R >
inline
RayC3<R CGAL_CTAG>
RayC3<R CGAL_CTAG>::opposite() const
{
  return RayC3<R CGAL_CTAG>( source(), - direction() );
}

template < class R >
inline
RayC3<R CGAL_CTAG>
RayC3<R CGAL_CTAG>::
transform(const typename RayC3<R CGAL_CTAG>::Aff_transformation_3 &t) const
{
  return RayC3<R CGAL_CTAG>(t.transform(source()),
                            t.transform(second_point()));
}

template < class R >
bool
RayC3<R CGAL_CTAG>::
has_on(const typename RayC3<R CGAL_CTAG>::Point_3 &p) const
{
  return (p == source()) ||
         ( collinear(source(), p, second_point())
           && ( Direction_3(p - source()) == direction() ));
}

template < class R >
inline
bool
RayC3<R CGAL_CTAG>::is_degenerate() const
{
  return source() == second_point();
}

template < class R >
inline
bool
RayC3<R CGAL_CTAG>::
collinear_has_on(const typename RayC3<R CGAL_CTAG>::Point_3 &p) const
{
  CGAL_kernel_exactness_precondition( collinear(source(), p, second_point()) );

  Comparison_result cx = compare_x(source(), second_point());
  if (cx != EQUAL)
    return cx != compare_x(p, source());

  Comparison_result cy = compare_y(source(), second_point());
  if (cy != EQUAL)
    return cy != compare_y(p, source());

  Comparison_result cz = compare_z(source(), second_point());
  if (cz != EQUAL)
    return cz != compare_z(p, source());

  return true; // p == source()
}

#ifndef CGAL_NO_OSTREAM_INSERT_RAYC3
template < class R >
std::ostream &
operator<<(std::ostream &os, const RayC3<R CGAL_CTAG> &r)
{
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << r.start() << ' ' << r.direction();
    case IO::BINARY :
        return os<< r.start() << r.direction();
    default:
        return os << "RayC3(" << r.start() <<  ", " << r.direction() << ")";
    }
}
#endif // CGAL_NO_OSTREAM_INSERT_RAYC3

#ifndef CGAL_NO_ISTREAM_EXTRACT_RAYC3
template < class R >
std::istream &
operator>>(std::istream &is, RayC3<R CGAL_CTAG> &r)
{
    typename RayC3<R CGAL_CTAG>::Point_3 p;
    typename RayC3<R CGAL_CTAG>::Direction_3 d;

    is >> p >> d;

    r = RayC3<R CGAL_CTAG>(p, d);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_RAYC3

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_TYPENAME_BUG
#undef typename
#endif

#endif // CGAL_CARTESIAN_RAY_3_C
