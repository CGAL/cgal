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

#include <CGAL/Cartesian/redefine_names_2.h>

CGAL_BEGIN_NAMESPACE

template < class R_ >
class VectorC2 CGAL_ADVANCED_KERNEL_PARTIAL_SPEC
  : public R_::Vector_handle_2
{
public:
  typedef R_                                    R;
  typedef typename R::FT                        FT;
  typedef typename R::RT                        RT;

  typedef typename R::Vector_handle_2		Vector_handle_2_;
  typedef typename Vector_handle_2_::element_type	Vector_ref_2;

#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
  typedef VectorC2<R,Cartesian_tag>             Self;
  typedef typename R::Point_2                   Point_2;
  typedef typename R::Direction_2               Direction_2;
  typedef typename R::Line_2                    Line_2;
  typedef typename R::Ray_2                     Ray_2;
  typedef typename R::Triangle_2                Triangle_2;
  typedef typename R::Segment_2                 Segment_2;
  typedef typename R::Iso_rectangle_2           Iso_rectangle_2;
  typedef typename R::Aff_transformation_2      Aff_transformation_2;
  typedef typename R::Circle_2                  Circle_2;
#else
  typedef VectorC2<R>                           Self;
  typedef typename R::Point_2_base              Point_2;
  typedef typename R::Direction_2_base          Direction_2;
  typedef typename R::Line_2_base               Line_2;
  typedef typename R::Ray_2_base                Ray_2;
  typedef typename R::Triangle_2_base           Triangle_2;
  typedef typename R::Segment_2_base            Segment_2;
  typedef typename R::Iso_rectangle_2_base      Iso_rectangle_2;
  typedef typename R::Aff_transformation_2_base Aff_transformation_2;
  typedef typename R::Circle_2_base             Circle_2;
#endif

  VectorC2()
    : Vector_handle_2_(Vector_ref_2()) {}

  VectorC2(const Null_vector &)
    : Vector_handle_2_(Vector_ref_2(FT(0), FT(0))) {}

  VectorC2(const Point_2 &p)
    : Vector_handle_2_(p) {}

  VectorC2(const Point_2 &a, const Point_2 &b)
    : Vector_handle_2_(b-a) {}

  VectorC2(const Direction_2 &d)
    : Vector_handle_2_(d) {}

  VectorC2(const FT &x, const FT &y)
    : Vector_handle_2_(Vector_ref_2(x, y)) {}

  VectorC2(const FT &hx, const FT &hy, const FT &hw)
  {
    if (hw != FT(1))
      initialize_with(Vector_ref_2(hx/hw, hy/hw));
    else
      initialize_with(Vector_ref_2(hx, hy));
  }

  bool operator==(const Self &v) const;
  bool operator!=(const Self &v) const;
  bool operator==(const Null_vector &) const;
  bool operator!=(const Null_vector &p) const;

  FT x() const
  {
      return Ptr()->e0;
  }
  FT y() const
  {
      return Ptr()->e1;
  }

  FT hx() const
  {
      return x();
  }
  FT hy() const
  {
      return y();
  }
  FT hw() const
  {
      return FT(1);
  }

  FT cartesian(int i) const;
  FT operator[](int i) const;
  FT homogeneous(int i) const;

  int dimension() const
  {
      return 2;
  }

  Self operator+(const Self &w) const;
  Self operator-(const Self &w) const;
  Self operator-() const;
  FT operator*(const Self &w) const;
  FT squared_length() const;
  Self operator/(const FT &c) const;
  Direction_2 direction() const;

  Self perpendicular(const Orientation &o) const;
  Self transform(const Aff_transformation_2 &t) const
  {
    return t.transform(*this);
  }
};

#ifdef CGAL_CFG_TYPENAME_BUG
#define typename
#endif

template < class R >
CGAL_KERNEL_INLINE
bool
VectorC2<R CGAL_CTAG>::operator==(const VectorC2<R CGAL_CTAG> &v) const
{
  if (identical(v))
      return true;
  return x() == v.x() && y() == v.y();
}

template < class R >
inline
bool
VectorC2<R CGAL_CTAG>::operator!=(const VectorC2<R CGAL_CTAG> &v) const
{
  return !(*this == v);
}

template < class R >
inline
bool
VectorC2<R CGAL_CTAG>::operator==(const Null_vector &) const
{
  return CGAL_NTS is_zero(x()) && CGAL_NTS is_zero(y());
}

template < class R >
inline
bool
VectorC2<R CGAL_CTAG>::operator!=(const Null_vector &v) const
{
  return !(*this == v);
}

template < class R >
CGAL_KERNEL_INLINE
typename VectorC2<R CGAL_CTAG>::FT 
VectorC2<R CGAL_CTAG>::cartesian(int i) const
{
  CGAL_kernel_precondition( (i == 0) || (i == 1) );
  return (i == 0) ? x() : y();
}

template < class R >
inline
typename VectorC2<R CGAL_CTAG>::FT
VectorC2<R CGAL_CTAG>::operator[](int i) const
{
  return cartesian(i);
}

template < class R >
CGAL_KERNEL_INLINE
typename VectorC2<R CGAL_CTAG>::FT 
VectorC2<R CGAL_CTAG>::homogeneous(int i) const
{
  return (i == 2) ? FT(1) : cartesian(i);
}

template < class R >
CGAL_KERNEL_INLINE
VectorC2<R CGAL_CTAG>
VectorC2<R CGAL_CTAG>::operator+(const VectorC2<R CGAL_CTAG> &w) const
{
  return VectorC2<R CGAL_CTAG>(x() + w.x(), y() + w.y());
}

template < class R >
CGAL_KERNEL_INLINE
VectorC2<R CGAL_CTAG>
VectorC2<R CGAL_CTAG>::operator-(const VectorC2<R CGAL_CTAG> &w) const
{
  return VectorC2<R CGAL_CTAG>(x() - w.x(), y() - w.y());
}

template < class R >
CGAL_KERNEL_INLINE
VectorC2<R CGAL_CTAG>
VectorC2<R CGAL_CTAG>::operator-() const
{
  return VectorC2<R CGAL_CTAG>(-x(), -y());
}

template < class R >
CGAL_KERNEL_INLINE
typename VectorC2<R CGAL_CTAG>::FT
VectorC2<R CGAL_CTAG>::operator*(const VectorC2<R CGAL_CTAG> &w) const
{
  return x() * w.x() + y() * w.y();
}

template < class R >
CGAL_KERNEL_INLINE
typename VectorC2<R CGAL_CTAG>::FT
VectorC2<R CGAL_CTAG>::squared_length() const
{
  return CGAL_NTS square(x()) + CGAL_NTS square(y());
}

template < class R >
CGAL_KERNEL_INLINE
VectorC2<R CGAL_CTAG>
VectorC2<R CGAL_CTAG>::
operator/(const typename VectorC2<R CGAL_CTAG>::FT &c) const
{
  return VectorC2<R CGAL_CTAG>( x()/c, y()/c);
}

template < class R >
inline
typename VectorC2<R CGAL_CTAG>::Direction_2
VectorC2<R CGAL_CTAG>::direction() const
{
  return Direction_2(*this);
}

template < class R >
CGAL_KERNEL_MEDIUM_INLINE
VectorC2<R CGAL_CTAG>
VectorC2<R CGAL_CTAG>::perpendicular(const Orientation &o) const
{
  CGAL_kernel_precondition( o != COLLINEAR );
  if (o == COUNTERCLOCKWISE)
    return VectorC2<R CGAL_CTAG>(-y(), x());
  else
    return VectorC2<R CGAL_CTAG>(y(), -x());
}

#ifndef CGAL_NO_OSTREAM_INSERT_VECTORC2
template < class R >
std::ostream &
operator<<(std::ostream &os, const VectorC2<R CGAL_CTAG> &v)
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
operator>>(std::istream &is, VectorC2<R CGAL_CTAG> &p)
{
    typename VectorC2<R CGAL_CTAG>::FT x, y;
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
	p = VectorC2<R CGAL_CTAG>(x, y);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_VECTORC2

#ifdef CGAL_CFG_TYPENAME_BUG
#undef typename
#endif

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_VECTOR_2_H
