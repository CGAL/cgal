// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// release       : 
// release_date  : 2000, August 16
//
// source        : webS2/S2.lw
// file          : include/CGAL/SimpleCartesian/VectorS2.h
// package       : S2 (1.7)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// revision      : 1.6
// revision_date : 27 Jun 2000
// author(s)     : Stefan Schirra
//                 based on code by
//                 Andreas Fabri and
//                 Herve Brönnimann
//
// coordinator   : MPI, Saarbrücken
// ======================================================================

#ifndef CGAL_VECTORS2_H
#define CGAL_VECTORS2_H

#include <CGAL/SimpleCartesian/PointS2.h>

CGAL_BEGIN_NAMESPACE

template < class FT >
class VectorS2
{
friend class DirectionS2<FT>;

public:
                  VectorS2();
                  VectorS2(const Null_vector& );
                  VectorS2(const FT& hx, const FT &hy, const FT &hw);
                  VectorS2(const FT& x, const FT &y);


  bool            operator==(const VectorS2<FT>& v) const;
  bool            operator!=(const VectorS2<FT>& v) const;

  bool            operator==(const Null_vector& ) const;
  bool            operator!=(const Null_vector& p) const;

  int             id() const;

  FT              hx() const;
  FT              hy() const;
  FT              hw() const;

  FT              x() const;
  FT              y() const;
  FT              cartesian(int i) const;
  FT              operator[](int i) const;

  FT              homogeneous(int i) const;

  int             dimension() const;

  VectorS2<FT>    operator+(const VectorS2<FT>& w) const;
  VectorS2<FT>    operator-(const VectorS2<FT>& w) const;
  VectorS2<FT>    operator-() const;
  FT               operator*(const VectorS2<FT>& w) const;
  VectorS2<FT>    operator/(const FT& c) const;
  DirectionS2<FT> direction() const;

  VectorS2<FT>    perpendicular(const Orientation& o) const;
  VectorS2<FT>    transform(const Aff_transformationS2<FT>& ) const;

protected:
                  VectorS2(const PointS2<FT>& p);
                  VectorS2(const DirectionS2<FT>& d);
public:
  FT  e0;
  FT  e1;
};


CGAL_END_NAMESPACE

#include <CGAL/SimpleCartesian/DirectionS2.h>

CGAL_BEGIN_NAMESPACE

template < class FT >
VectorS2<FT>::VectorS2() {}

template < class FT >
CGAL_KERNEL_CTOR_INLINE
VectorS2<FT>::VectorS2(const Null_vector& ) : e0(FT(0)), e1(FT(0)) {}

template < class FT >
CGAL_KERNEL_CTOR_MEDIUM_INLINE
VectorS2<FT>::VectorS2(const FT& hx, const FT &hy, const FT &hw)
{
  if( hw != FT(1))
  { e0 = hx/hw; e1 = hy/hw; }
  else
  { e0 = hx; e1 = hy; }
}

template < class FT >
CGAL_KERNEL_CTOR_INLINE
VectorS2<FT>::VectorS2(const FT& x, const FT &y) : e0(x), e1(y) {}

template < class FT >
CGAL_KERNEL_CTOR_INLINE
VectorS2<FT>::VectorS2(const DirectionS2<FT>& d)
{ e0 = d.e0; e1 = d.e1; }

template < class FT >
CGAL_KERNEL_CTOR_INLINE
VectorS2<FT>::VectorS2(const PointS2<FT>& d)
{ e0 = d.e0; e1 = d.e1; }

template < class FT >
CGAL_KERNEL_INLINE
bool VectorS2<FT>::operator==(const VectorS2<FT>& v) const
{ return (x() == v.x()) && (y() == v.y()); }

template < class FT >
inline
bool VectorS2<FT>::operator!=(const VectorS2<FT>& v) const
{ return !(*this == v); }

template < class FT >
inline
bool VectorS2<FT>::operator==(const Null_vector& ) const
{ return (x() == FT(0)) && (y() == FT(0)); }

template < class FT >
inline
bool VectorS2<FT>::operator!=(const Null_vector& v) const
{ return !(*this == v); }

template < class FT >
inline
FT VectorS2<FT>::x()  const
{
  return e0;
}

template < class FT >
inline
FT VectorS2<FT>::y()  const
{
  return e1;
}

template < class FT >
CGAL_KERNEL_INLINE
FT  VectorS2<FT>::cartesian(int i) const
{
  CGAL_kernel_precondition( (i == 0) || (i == 1) );
  return (i == 0) ? x() : y();
}

template < class FT >
inline
FT  VectorS2<FT>::operator[](int i) const
{ return cartesian(i); }

template < class FT >
inline
int VectorS2<FT>::dimension() const
{ return 2; }

template < class FT >
inline
FT VectorS2<FT>::hx()  const
{ return e0; }

template < class FT >
inline
FT VectorS2<FT>::hy()  const
{ return e1; }

template < class FT >
inline
FT VectorS2<FT>::hw()  const
{ return FT(1); }

template < class FT >
CGAL_KERNEL_INLINE
FT  VectorS2<FT>::homogeneous(int i) const
{ return (i == 2) ? FT(1) : cartesian(i); }

template < class FT >
CGAL_KERNEL_INLINE
VectorS2<FT> VectorS2<FT>::operator+(const VectorS2<FT>& w) const
{ return VectorS2<FT>(x() + w.x(), y() + w.y()) ; }

template < class FT >
CGAL_KERNEL_INLINE
VectorS2<FT> VectorS2<FT>::operator-(const VectorS2<FT>& w) const
{ return VectorS2<FT>(x() - w.x(), y() - w.y()) ; }

template < class FT >
CGAL_KERNEL_INLINE
VectorS2<FT> VectorS2<FT>::operator-() const
{ return VectorS2<FT>(-x(), -y()) ; }

template < class FT >
CGAL_KERNEL_INLINE
FT VectorS2<FT>::operator*(const VectorS2<FT>& w) const
{ return x() * w.x() + y() * w.y() ; }

template < class FT >
CGAL_KERNEL_INLINE
VectorS2<FT> operator*(const FT& c, const VectorS2<FT> &w)
{ return VectorS2<FT>( c* w.x(), c * w.y()) ; }

template < class FT >
CGAL_KERNEL_INLINE
VectorS2<FT> VectorS2<FT>::operator/(const FT& c) const
{ return VectorS2<FT>( x()/c, y()/c) ; }

template < class FT >
inline
DirectionS2<FT>   VectorS2<FT>::direction() const
{ return DirectionS2<FT>(*this) ; }

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
VectorS2<FT> VectorS2<FT>::perpendicular(const Orientation& o) const
{
  CGAL_kernel_precondition( o != COLLINEAR );
  if (o == COUNTERCLOCKWISE)
    return VectorS2<FT>(-y(), x());
  else
    return VectorS2<FT>(y(), -x());
}

template < class FT >
inline
VectorS2<FT> VectorS2<FT>::transform(const Aff_transformationS2<FT>& t) const
{ return t.transform(*this); }



#ifndef CGAL_NO_OSTREAM_INSERT_VECTORS2
template < class FT >
std::ostream& operator<<(std::ostream &os, const VectorS2<FT> &v)
{
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << v.x() << ' ' << v.y();
    case IO::BINARY :
        write(os, v.x());
        write(os, v.y());
        return os;
    default:
        return os << "VectorS2(" << v.x() << ", " << v.y() << ')';
    }
}
#endif // CGAL_NO_OSTREAM_INSERT_VECTORS2

#ifndef CGAL_NO_ISTREAM_EXTRACT_VECTORS2
template < class FT >
std::istream& operator>>(std::istream &is, VectorS2<FT> &p)
{
    FT x, y;
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
    p = VectorS2<FT>(x, y);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_VECTORS2



CGAL_END_NAMESPACE

#endif
