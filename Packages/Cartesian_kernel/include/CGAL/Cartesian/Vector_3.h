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
// file          : include/CGAL/Cartesian/Vector_3.h
// revision      : $Revision$
// revision_date : $Date$
// author        : Andreas Fabri
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_VECTOR_3_H
#define CGAL_CARTESIAN_VECTOR_3_H

#include <CGAL/Cartesian/redefine_names_3.h>
#include <CGAL/Threetuple.h>
#include <CGAL/Cartesian/Direction_3.h>

CGAL_BEGIN_NAMESPACE

template < class R_ >
class VectorC3 CGAL_ADVANCED_KERNEL_PARTIAL_SPEC
  : public Handle_for<Threetuple<typename R_::FT> >
{
public:
  typedef R_                               R;
  typedef typename R::FT                   FT;
  typedef typename R::RT                   RT;
#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
  typedef VectorC3<R CGAL_CTAG>            Self;
  typedef typename R::Point_3              Point_3;
  typedef typename R::Direction_3          Direction_3;
  typedef typename R::Aff_transformation_3 Aff_transformation_3;
#else
  typedef VectorC3<R>                           Self;
  typedef typename R::Point_3_base              Point_3;
  typedef typename R::Direction_3_base          Direction_3;
  typedef typename R::Aff_transformation_3_base Aff_transformation_3;
#endif

  VectorC3();
  VectorC3(const Self &v);
  VectorC3(const Null_vector &);
  VectorC3(const Point_3 &p);
  VectorC3(const Direction_3 &p);
  VectorC3(const FT &x, const FT &y, const FT &z);
  VectorC3(const FT &x, const FT &y, const FT &z, const FT &w);

  bool operator==(const Self &p) const;
  bool operator!=(const Self &p) const;

  bool operator==(const Null_vector &) const;
  bool operator!=(const Null_vector &) const;

  FT x() const
  {
      return ptr->e0;
  }
  FT y() const
  {
      return ptr->e1;
  }
  FT z() const
  {
      return ptr->e2;
  }

  FT hx() const
  {
      return x();
  }
  FT hy() const
  {
      return y();
  }
  FT hz() const
  {
      return z();
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
      return 3;
  }

  Self operator+(const Self &w) const;
  Self operator-(const Self &w) const;
  Self operator-() const;
  Self operator/(const FT &c) const;
  FT operator*(const Self &w) const;
  FT squared_length() const;
  Direction_3 direction() const;
  Self transform(const Aff_transformation_3 &) const;
};

template < class R >
CGAL_KERNEL_CTOR_INLINE
VectorC3<R CGAL_CTAG>::VectorC3()
  : Handle_for<Threetuple<typename R::FT> >( Threetuple<FT>()) {}

template < class R >
CGAL_KERNEL_CTOR_INLINE
VectorC3<R CGAL_CTAG>::VectorC3(const VectorC3<R CGAL_CTAG> &v)
  : Handle_for<Threetuple<typename R::FT> >(v) {}

template < class R >
CGAL_KERNEL_CTOR_INLINE
VectorC3<R CGAL_CTAG>::VectorC3(const Null_vector &)
  : Handle_for<Threetuple<typename R::FT> >( Threetuple<FT>(FT(0), FT(0), FT(0))) {}

template < class R >
CGAL_KERNEL_CTOR_INLINE
VectorC3<R CGAL_CTAG>::
VectorC3(const typename VectorC3<R CGAL_CTAG>::FT &x,
         const typename VectorC3<R CGAL_CTAG>::FT &y,
         const typename VectorC3<R CGAL_CTAG>::FT &z)
  : Handle_for<Threetuple<typename R::FT> >( Threetuple<FT>(x, y, z)) {}

template < class R >
CGAL_KERNEL_CTOR_INLINE
VectorC3<R CGAL_CTAG>::
VectorC3(const typename VectorC3<R CGAL_CTAG>::FT &x,
         const typename VectorC3<R CGAL_CTAG>::FT &y,
         const typename VectorC3<R CGAL_CTAG>::FT &z,
         const typename VectorC3<R CGAL_CTAG>::FT &w)
{
  if (w != FT(1))
    initialize_with( Threetuple<FT>(x/w, y/w, z/w) );
  else
    initialize_with( Threetuple<FT>(x, y, z) );
}

template < class R >
CGAL_KERNEL_CTOR_INLINE
VectorC3<R CGAL_CTAG>::
VectorC3(const typename VectorC3<R CGAL_CTAG>::Point_3 &p)
  : Handle_for<Threetuple<typename R::FT> >(p) {}

template < class R >
CGAL_KERNEL_CTOR_INLINE
VectorC3<R CGAL_CTAG>::
VectorC3(const typename VectorC3<R CGAL_CTAG>::Direction_3 &d)
  : Handle_for<Threetuple<typename R::FT> >(d) {}

template < class R >
inline
bool
VectorC3<R CGAL_CTAG>::operator==(const VectorC3<R CGAL_CTAG> &v) const
{
  if ( identical(v) ) return true;
  return x() == v.x() && y() == v.y() && z() == v.z();
}

template < class R >
inline
bool
VectorC3<R CGAL_CTAG>::operator!=(const VectorC3<R CGAL_CTAG> &v) const
{
  return !(*this == v);
}

template < class R >
inline
bool
VectorC3<R CGAL_CTAG>::operator==(const Null_vector &) const
{
  return CGAL_NTS is_zero(x()) && CGAL_NTS is_zero(y()) &&
         CGAL_NTS is_zero(z());
}

template < class R >
inline
bool
VectorC3<R CGAL_CTAG>::operator!=(const Null_vector &v) const
{
  return !(*this == v);
}

template < class R >
inline
typename VectorC3<R CGAL_CTAG>::FT
VectorC3<R CGAL_CTAG>::cartesian(int i) const
{
  CGAL_kernel_precondition( (i>=0) && (i<3) );
  if (i==0) return x();
  if (i==1) return y();
  return z();
}

template < class R >
inline
typename VectorC3<R CGAL_CTAG>::FT
VectorC3<R CGAL_CTAG>::operator[](int i) const
{
  return cartesian(i);
}

template < class R >
typename VectorC3<R CGAL_CTAG>::FT
VectorC3<R CGAL_CTAG>::homogeneous(int i) const
{
  if (i==3) return FT(1);
  return cartesian(i);
}

template < class R >
inline
VectorC3<R CGAL_CTAG>
VectorC3<R CGAL_CTAG>::
operator+(const VectorC3<R CGAL_CTAG> &w) const
{
  return VectorC3<R CGAL_CTAG>(x() + w.x(), y() + w.y(), z() + w.z());
}

template < class R >
inline
VectorC3<R CGAL_CTAG>
VectorC3<R CGAL_CTAG>::operator-(const VectorC3<R CGAL_CTAG> &w) const
{
  return VectorC3<R CGAL_CTAG>(x() - w.x(), y() - w.y(), z() - w.z());
}

template < class R >
inline
VectorC3<R CGAL_CTAG>
VectorC3<R CGAL_CTAG>::operator-() const
{
  return VectorC3<R CGAL_CTAG>(-x(), -y(), -z());
}

template < class R >
inline
typename VectorC3<R CGAL_CTAG>::FT
VectorC3<R CGAL_CTAG>::operator*(const VectorC3<R CGAL_CTAG> &w) const
{
  return x() * w.x() + y() * w.y() + z() * w.z();
}

template < class R >
inline
typename VectorC3<R CGAL_CTAG>::FT
VectorC3<R CGAL_CTAG>::squared_length() const
{
  return CGAL_NTS square(x()) + CGAL_NTS square(y()) + CGAL_NTS square(z());
}

template < class R >
inline
VectorC3<R CGAL_CTAG>
VectorC3<R CGAL_CTAG>::
operator/(const typename VectorC3<R CGAL_CTAG>::FT &c) const
{
  return VectorC3<R CGAL_CTAG>( x()/c, y()/c, z()/c);
}

template < class R >
inline
typename VectorC3<R CGAL_CTAG>::Direction_3
VectorC3<R CGAL_CTAG>::direction() const
{
  return Direction_3(*this);
}

template < class R >
VectorC3<R CGAL_CTAG>
VectorC3<R CGAL_CTAG>::
transform(const typename VectorC3<R CGAL_CTAG>::Aff_transformation_3 &t) const
{
  return t.transform(*this);
}

#ifndef CGAL_CARTESIAN_NO_OSTREAM_INSERT_VECTORC3
template < class R >
std::ostream &
operator<<(std::ostream &os, const VectorC3<R CGAL_CTAG> &v)
{
  switch(os.iword(IO::mode)) {
    case IO::ASCII :
      return os << v.x() << ' ' << v.y()  << ' ' << v.z();
    case IO::BINARY :
      write(os, v.x());
      write(os, v.y());
      write(os, v.z());
      return os;
    default:
      os << "VectorC3(" << v.x() << ", " << v.y() <<  ", " << v.z() << ")";
      return os;
  }
}
#endif // CGAL_CARTESIAN_NO_OSTREAM_INSERT_VECTORC3

#ifndef CGAL_CARTESIAN_NO_ISTREAM_EXTRACT_VECTORC3
template < class R >
std::istream &
operator>>(std::istream &is, VectorC3<R CGAL_CTAG> &p)
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
  p = VectorC3<R CGAL_CTAG>(x, y, z);
  return is;
}
#endif // CGAL_CARTESIAN_NO_ISTREAM_EXTRACT_VECTORC3

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_VECTOR_3_H
