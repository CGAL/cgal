// ======================================================================
//
// Copyright (c) 2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Cartesian/Direction_2.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri, Herve Bronnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_DIRECTION_2_H
#define CGAL_CARTESIAN_DIRECTION_2_H

#include <CGAL/Twotuple.h>

CGAL_BEGIN_NAMESPACE

template < class R_ >
class DirectionC2
  : public R_::template Handle<Twotuple<typename R_::FT> >::type
{
CGAL_VC7_BUG_PROTECTED
  typedef typename R_::FT                   FT;
  typedef typename R_::Point_2              Point_2;
  typedef typename R_::Vector_2             Vector_2;
  typedef typename R_::Line_2               Line_2;
  typedef typename R_::Ray_2                Ray_2;
  typedef typename R_::Segment_2            Segment_2;
  typedef typename R_::Direction_2          Direction_2;
  typedef typename R_::Aff_transformation_2 Aff_transformation_2;

  typedef Twotuple<FT>	                           rep;
  typedef typename R_::template Handle<rep>::type  base;

public:
  typedef R_                                     R;

  DirectionC2()
    : base(rep()) {}

  DirectionC2(const Vector_2 &v)
    : base(v) {}

  DirectionC2(const Line_2 &l)
    : base(l.direction()) {}

  DirectionC2(const Ray_2 &r)
    : base(r.direction()) {}

  DirectionC2(const Segment_2 &s)
    : base(s.direction()) {}

  DirectionC2(const FT &x, const FT &y)
    : base(rep(x, y)) {}

  bool operator==(const DirectionC2 &d) const;
  bool operator!=(const DirectionC2 &d) const;
  bool operator>=(const DirectionC2 &d) const;
  bool operator<=(const DirectionC2 &d) const;
  bool operator>(const DirectionC2 &d) const;
  bool operator<(const DirectionC2 &d) const;
  bool counterclockwise_in_between( const DirectionC2 &d1,
	                            const DirectionC2 &d2) const;
  
  Vector_2 to_vector() const;
  Vector_2 vector() const { return to_vector(); }

  Direction_2 perpendicular(const Orientation &o) const;
  Direction_2 transform(const Aff_transformation_2 &t) const
  {
    return t.transform(*this);
  }

  Direction_2 operator-() const;

  const FT & delta(int i) const;
  const FT & dx() const
  {
      return Ptr()->e0;
  }
  const FT & dy() const
  {
      return Ptr()->e1;
  }
};

#ifdef CGAL_CFG_TYPENAME_BUG
#define typename
#endif

template < class R >
inline
bool
DirectionC2<R>::operator==(const DirectionC2<R> &d) const
{
  if (identical(d))
      return true;
  return equal_direction(*this, d);
}

template < class R >
inline
bool
DirectionC2<R>::operator!=(const DirectionC2<R> &d) const
{
  return !( *this == d );
}

template < class R >
CGAL_KERNEL_MEDIUM_INLINE
bool
DirectionC2<R>::operator<(const DirectionC2<R> &d) const
{
  return compare_angle_with_x_axis(*this, d) == SMALLER;
}

template < class R >
CGAL_KERNEL_INLINE
bool
DirectionC2<R>::operator>(const DirectionC2<R> &d) const
{
  return d < *this;
}

template < class R >
CGAL_KERNEL_INLINE
bool
DirectionC2<R>::operator>=(const DirectionC2<R> &d) const
{
  return compare_angle_with_x_axis(*this, d) != SMALLER;
}

template < class R >
CGAL_KERNEL_INLINE
bool
DirectionC2<R>::operator<=(const DirectionC2<R> &d) const
{
  return compare_angle_with_x_axis(*this, d) != LARGER;
}

template < class R >
CGAL_KERNEL_INLINE
bool
DirectionC2<R>::
counterclockwise_in_between(const DirectionC2<R> &d1,
                            const DirectionC2<R> &d2) const
// returns true, iff \ccVar\ is not equal to \ccc{d1}, and 
// while rotating counterclockwise starting at \ccc{d1}, 
// \ccVar\ is reached strictly before \ccc{d2} is reached.
// Note that true is returned if \ccc{d1} == \ccc{d2}, unless
//  also \ccVar\ == \ccc{d1}.
{
  if ( d1 < *this)
    return ( *this < d2 )||( d2 <= d1 );
  else
    return ( *this < d2 )&&( d2 <= d1 );
}

template < class R >
inline
typename DirectionC2<R>::Vector_2
DirectionC2<R>::to_vector() const
{
  return Vector_2(dx(), dy());
}

template < class R >
CGAL_KERNEL_MEDIUM_INLINE
typename DirectionC2<R>::Direction_2
DirectionC2<R>::perpendicular(const Orientation &o) const
{
  CGAL_kernel_precondition(o != COLLINEAR);
  if (o == COUNTERCLOCKWISE)
    return DirectionC2<R>(-dy(), dx());
  else
    return DirectionC2<R>(dy(), -dx());
}

template < class R >
inline
typename DirectionC2<R>::Direction_2
DirectionC2<R>::operator-() const
{
  return DirectionC2<R>(-dx(), -dy());
}

template < class R >
CGAL_KERNEL_INLINE
const typename DirectionC2<R>::FT &
DirectionC2<R>::delta(int i) const
{
  CGAL_kernel_precondition( ( i == 0 ) || ( i == 1 ) );
  return (i==0) ? dx() : dy();
}

#ifndef CGAL_NO_OSTREAM_INSERT_DIRECTIONC2
template < class R >
std::ostream&
operator<<(std::ostream &os, const DirectionC2<R> &d)
{
    typename R::Vector_2 v = d.to_vector();
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << v.x() << ' ' << v.y();
    case IO::BINARY :
        write(os, v.x());
        write(os, v.y());
        return os;
    default:
        return os << "DirectionC2(" << v.x() << ", " << v.y() << ')';
    }
}
#endif // CGAL_NO_OSTREAM_INSERT_DIRECTIONC2

#ifndef CGAL_NO_ISTREAM_EXTRACT_DIRECTIONC2
template < class R >
std::istream&
operator>>(std::istream &is, DirectionC2<R> &p)
{
    typename R::FT x, y;
    switch(is.iword(IO::mode)) {
    case IO::ASCII :
        is >> x >> y;
        break;
    case IO::BINARY :
        read(is, x);
        read(is, y);
        break;
    default:
        std::cerr << std::endl << "Stream must be in ascii or binary mode"
	          << std::endl;
        break;
    }
    if (is)
	p = DirectionC2<R>(x, y);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_DIRECTIONC2

#ifdef CGAL_CFG_TYPENAME_BUG
#undef typename
#endif

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_DIRECTION_2_H
