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
// file          : include/CGAL/Cartesian/Direction_3.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_DIRECTION_3_H
#define CGAL_CARTESIAN_DIRECTION_3_H

#include <CGAL/Threetuple.h>

CGAL_BEGIN_NAMESPACE

template < class R_ >
class DirectionC3
  : public R_::template Handle<Threetuple<typename R_::FT> >::type
{
CGAL_VC7_BUG_PROTECTED
  typedef typename R_::FT                   FT;
  typedef typename R_::Vector_3             Vector_3;
  typedef typename R_::Line_3               Line_3;
  typedef typename R_::Ray_3                Ray_3;
  typedef typename R_::Segment_3            Segment_3;
  typedef typename R_::Direction_3          Direction_3;
  typedef typename R_::Aff_transformation_3 Aff_transformation_3;

  typedef Threetuple<FT>                           rep;
  typedef typename R_::template Handle<rep>::type  base;

public:
  typedef R_                                R;

  DirectionC3()
    : base(rep()) {}

  DirectionC3(const Vector_3 &v)
    : base(v) {}

  DirectionC3(const Line_3 &l)
    : base(l.direction()) {}

  DirectionC3(const Ray_3 &r)
    : base(r.direction()) {}

  DirectionC3(const Segment_3 &s)
    : base(s.direction()) {}

  DirectionC3(const FT &x, const FT &y, const FT &z)
    : base(rep(x, y, z)) {}

  bool           operator==(const DirectionC3 &d) const;
  bool           operator!=(const DirectionC3 &d) const;

  Vector_3       to_vector() const;
  Vector_3       vector() const { return to_vector(); }

  Direction_3    transform(const Aff_transformation_3 &t) const
  {
    return t.transform(*this);
  }

  Direction_3    operator-() const;

  const FT & delta(int i) const;
  const FT & dx() const
  {
      return Ptr()->e0;
  }
  const FT & dy() const
  {
      return Ptr()->e1;
  }
  const FT & dz() const
  {
      return Ptr()->e2;
  }

  const FT & hdx() const
  {
      return dx();
  }
  const FT & hdy() const
  {
      return dy();
  }
  const FT & hdz() const
  {
      return dz();
  }
  FT hw() const
  {
      return FT(1);
  }
};

template < class R >
inline
bool
DirectionC3<R>::operator==(const DirectionC3<R> &d) const
{
  if (identical(d))
      return true;
  return equal_directionC3(dx(), dy(), dz(), d.dx(), d.dy(), d.dz());
}

template < class R >
inline
bool
DirectionC3<R>::operator!=(const DirectionC3<R> &d) const
{
  return !(*this == d);
}

template < class R >
inline
typename DirectionC3<R>::Vector_3
DirectionC3<R>::to_vector() const
{
  return Vector_3(dx(), dy(), dz());
}

template < class R >
inline
typename DirectionC3<R>::Direction_3
DirectionC3<R>::operator-() const
{
  return DirectionC3<R>(-dx(), -dy(), -dz());
}

template < class R >
const typename DirectionC3<R>::FT &
DirectionC3<R>::delta(int i) const
{
  CGAL_kernel_precondition( i >= 0 && i <= 2 );
  if (i==0) return dx();
  if (i==1) return dy();
  return dz();
}

#ifndef CGAL_NO_OSTREAM_INSERT_DIRECTIONC3
template < class R >
std::ostream &
operator<<(std::ostream &os, const DirectionC3<R> &d)
{
  typename R::Vector_3 v = d.to_vector();
  switch(os.iword(IO::mode)) {
    case IO::ASCII :
      return os << v.x() << ' ' << v.y()  << ' ' << v.z();
    case IO::BINARY :
      write(os, v.x());
      write(os, v.y());
      write(os, v.z());
      return os;
    default:
      os << "DirectionC3(" << v.x() << ", " << v.y() << ", " << v.z() << ")";
      return os;
  }
}
#endif // CGAL_NO_OSTREAM_INSERT_DIRECTIONC3

#ifndef CGAL_NO_ISTREAM_EXTRACT_DIRECTIONC3
template < class R >
std::istream &
operator>>(std::istream &is, DirectionC3<R> &d)
{
  typename R::FT x, y, z;
  switch(is.iword(IO::mode)) {
    case IO::ASCII :
      is >> x >> y >> z;
      break;
    case IO::BINARY :
      read(is, x);
      read(is, y);
      read(is, z);
      break;
    default:
      std::cerr << "" << std::endl;
      std::cerr << "Stream must be in ascii or binary mode" << std::endl;
      break;
  }
  if (is)
      d = DirectionC3<R>(x, y, z);
  return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_DIRECTIONC3

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_DIRECTION_3_H
