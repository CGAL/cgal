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
// file          : include/CGAL/Cartesian/Vector_2.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri, Herve Bronnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_VECTOR_2_H
#define CGAL_CARTESIAN_VECTOR_2_H

#include <CGAL/Origin.h>
#include <CGAL/Twotuple.h>

CGAL_BEGIN_NAMESPACE

template < class R_ >
class VectorC2
  : public R_::template Handle<Twotuple<typename R_::FT> >::type
{
CGAL_VC7_BUG_PROTECTED
  typedef typename R_::FT                   FT;
  typedef typename R_::Point_2              Point_2;
  typedef typename R_::Direction_2          Direction_2;
  typedef typename R_::Vector_2             Vector_2;
  typedef typename R_::Aff_transformation_2 Aff_transformation_2;

  typedef Twotuple<FT>	                           rep;
  typedef typename R_::template Handle<rep>::type  base;

public:
  typedef R_                                     R;

  VectorC2()
    : base(rep()) {}

  VectorC2(const Null_vector &)
    : base(rep(FT(0), FT(0))) {}

  VectorC2(const Point_2 &p)
    : base(p) {}

  VectorC2(const Point_2 &a, const Point_2 &b)
    : base(b-a) {}

  VectorC2(const Direction_2 &d)
    : base(d) {}

  VectorC2(const FT &x, const FT &y)
    : base(rep(x, y)) {}

  VectorC2(const FT &hx, const FT &hy, const FT &hw)
  {
    if (hw != FT(1))
      initialize_with(rep(hx/hw, hy/hw));
    else
      initialize_with(rep(hx, hy));
  }

  bool operator==(const VectorC2 &v) const;
  bool operator!=(const VectorC2 &v) const;
  bool operator==(const Null_vector &) const;
  bool operator!=(const Null_vector &p) const;

  const FT & x() const
  {
      return Ptr()->e0;
  }
  const FT & y() const
  {
      return Ptr()->e1;
  }

  const FT & hx() const
  {
      return x();
  }
  const FT & hy() const
  {
      return y();
  }
  FT hw() const
  {
      return FT(1);
  }

  const FT & cartesian(int i) const;
  const FT & operator[](int i) const;
  FT homogeneous(int i) const;

  int dimension() const
  {
      return 2;
  }

  Vector_2 operator+(const VectorC2 &w) const;
  Vector_2 operator-(const VectorC2 &w) const;
  Vector_2 operator-() const;
  FT operator*(const VectorC2 &w) const;
  FT squared_length() const;
  Vector_2 operator/(const FT &c) const;
  Direction_2 direction() const;

  Vector_2 perpendicular(const Orientation &o) const;
  Vector_2 transform(const Aff_transformation_2 &t) const
  {
    return t.transform(*this);
  }
};

template < class R >
CGAL_KERNEL_INLINE
bool
VectorC2<R>::operator==(const VectorC2<R> &v) const
{
  if (identical(v))
      return true;
  return x() == v.x() && y() == v.y();
}

template < class R >
inline
bool
VectorC2<R>::operator!=(const VectorC2<R> &v) const
{
  return !(*this == v);
}

template < class R >
inline
bool
VectorC2<R>::operator==(const Null_vector &) const
{
  return CGAL_NTS is_zero(x()) && CGAL_NTS is_zero(y());
}

template < class R >
inline
bool
VectorC2<R>::operator!=(const Null_vector &v) const
{
  return !(*this == v);
}

template < class R >
CGAL_KERNEL_INLINE
const typename VectorC2<R>::FT &
VectorC2<R>::cartesian(int i) const
{
  CGAL_kernel_precondition( (i == 0) || (i == 1) );
  return (i == 0) ? x() : y();
}

template < class R >
inline
const typename VectorC2<R>::FT &
VectorC2<R>::operator[](int i) const
{
  return cartesian(i);
}

template < class R >
CGAL_KERNEL_INLINE
typename VectorC2<R>::FT 
VectorC2<R>::homogeneous(int i) const
{
  return (i == 2) ? FT(1) : cartesian(i);
}

template < class R >
CGAL_KERNEL_INLINE
typename VectorC2<R>::Vector_2
VectorC2<R>::operator+(const VectorC2<R> &w) const
{
  return VectorC2<R>(x() + w.x(), y() + w.y());
}

template < class R >
CGAL_KERNEL_INLINE
typename VectorC2<R>::Vector_2
VectorC2<R>::operator-(const VectorC2<R> &w) const
{
  return VectorC2<R>(x() - w.x(), y() - w.y());
}

template < class R >
CGAL_KERNEL_INLINE
typename VectorC2<R>::Vector_2
VectorC2<R>::operator-() const
{
  return VectorC2<R>(-x(), -y());
}

template < class R >
CGAL_KERNEL_INLINE
typename VectorC2<R>::FT
VectorC2<R>::operator*(const VectorC2<R> &w) const
{
  return x() * w.x() + y() * w.y();
}

template < class R >
CGAL_KERNEL_INLINE
typename VectorC2<R>::FT
VectorC2<R>::squared_length() const
{
  return CGAL_NTS square(x()) + CGAL_NTS square(y());
}

template < class R >
CGAL_KERNEL_INLINE
typename VectorC2<R>::Vector_2
VectorC2<R>::
operator/(const typename VectorC2<R>::FT &c) const
{
  return VectorC2<R>( x()/c, y()/c);
}

template < class R >
inline
typename VectorC2<R>::Direction_2
VectorC2<R>::direction() const
{
  return Direction_2(*this);
}

template < class R >
CGAL_KERNEL_MEDIUM_INLINE
typename VectorC2<R>::Vector_2
VectorC2<R>::perpendicular(const Orientation &o) const
{
  CGAL_kernel_precondition( o != COLLINEAR );
  if (o == COUNTERCLOCKWISE)
    return VectorC2<R>(-y(), x());
  else
    return VectorC2<R>(y(), -x());
}

#ifndef CGAL_NO_OSTREAM_INSERT_VECTORC2
template < class R >
std::ostream &
operator<<(std::ostream &os, const VectorC2<R> &v)
{
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << v.x() << ' ' << v.y();
    case IO::BINARY :
        write(os, v.x());
        write(os, v.y());
        return os;
    default:
        return os << "VectorC2(" << v.x() << ", " << v.y() << ')';
    }
}
#endif // CGAL_NO_OSTREAM_INSERT_VECTORC2

#ifndef CGAL_NO_ISTREAM_EXTRACT_VECTORC2
template < class R >
std::istream &
operator>>(std::istream &is, VectorC2<R> &p)
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
        std::cerr << "" << std::endl;
        std::cerr << "Stream must be in ascii or binary mode" << std::endl;
        break;
    }
    if (is)
	p = VectorC2<R>(x, y);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_VECTORC2

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_VECTOR_2_H
