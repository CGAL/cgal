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
// file          : include/CGAL/SimpleCartesian/PointS2.h
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


#ifndef CGAL_POINTS2_H
#define CGAL_POINTS2_H

CGAL_BEGIN_NAMESPACE

template < class FT >
class PointS2
{
 public:
              PointS2() {}
              PointS2(const Origin& ) : e0(FT(0)), e1(FT(0)) {}
              PointS2(const FT& hx, const FT& hy, const FT& hw);
              PointS2(const FT& x,  const FT& y);

  bool        operator==(const PointS2<FT>& p) const;
  bool        operator!=(const PointS2<FT>& p) const;

  const FT&   hx() const;
  const FT&   hy() const;
  FT          hw() const;
  const FT&   x()  const;
  const FT&   y()  const;
  FT          cartesian(int i) const;
  FT          operator[](int i) const;
  FT          homogeneous(int i) const;

  int         dimension() const;
  Bbox_2      bbox() const;


  PointS2<FT> transform(const Aff_transformationS2<FT>& ) const;

// protected:
              PointS2(const VectorS2<FT>& v);

  FT          e0;
  FT          e1;
};


CGAL_END_NAMESPACE

#include <CGAL/Origin.h>
#include <CGAL/SimpleCartesian/VectorS2.h>
#include <CGAL/SimpleCartesian/Aff_transformationS2.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/number_utils.h>

CGAL_BEGIN_NAMESPACE

template < class FT >
CGAL_KERNEL_CTOR_INLINE
PointS2<FT>::PointS2(const VectorS2<FT>& v)
{ e0 = v.e0; e1 = v.e1; }

template < class FT >
CGAL_KERNEL_CTOR_INLINE
PointS2<FT>::PointS2(const FT& hx, const FT& hy, const FT& hw)
{
  if( hw != FT(1))
  { e0 = hx/hw; e1 = hy/hw; }
  else
  { e0 = hx; e1 = hy; }
}

template < class FT >
CGAL_KERNEL_CTOR_INLINE
PointS2<FT>::PointS2(const FT& x, const FT& y)
{ e0 = x; e1 = y; }

template < class FT >
inline
bool 
PointS2<FT>::operator==(const PointS2<FT>& p) const
{ return ((x() == p.x()) && (y() == p.y())) ; }

template < class FT >
inline
bool 
PointS2<FT>::operator!=(const PointS2<FT>& p) const
{ return !(*this == p); }

template < class FT >
inline
const FT& 
PointS2<FT>::x()  const
{ return e0; }

template < class FT >
inline
const FT& 
PointS2<FT>::y()  const
{ return  e1 ; }

template < class FT >
CGAL_KERNEL_INLINE
FT 
PointS2<FT>::cartesian(int i) const
{
  CGAL_kernel_precondition( (i == 0) || (i == 1) );
  return (i == 0) ? x() : y();
}

template < class FT >
inline
FT  
PointS2<FT>::operator[](int i) const
{ return cartesian(i); }

template < class FT >
inline
int 
PointS2<FT>::dimension() const
{ return 2; }

template < class FT >
inline
const FT& 
PointS2<FT>::hx()  const
{ return e0; }

template < class FT >
inline
const FT& 
PointS2<FT>::hy()  const
{ return e1; }

template < class FT >
inline
FT 
PointS2<FT>::hw()  const
{ return FT(1); }

template < class FT >
inline
FT  
PointS2<FT>::homogeneous(int i) const
{
  CGAL_kernel_precondition( (i>=0) && (i<=2) );
  if (i<2) return cartesian(i);
  return FT(1);
}

template < class FT >
inline
PointS2<FT> 
operator+(const PointS2<FT>& p, const VectorS2<FT>& v)
{ return PointS2<FT>(p.x() + v.x(), p.y() + v.y()) ; }

template < class FT >
inline
PointS2<FT> 
operator-(const PointS2<FT>& p, const VectorS2<FT>& v)
{ return PointS2<FT>(p.x() - v.x(), p.y() - v.y()) ; }

template < class FT >
inline
PointS2<FT> 
operator+(const Origin& , const VectorS2<FT>& v)
{ return PointS2<FT>(v) ; }

template < class FT >
inline
PointS2<FT> operator-(const Origin& , const VectorS2<FT>& v)
{ return PointS2<FT>(-v) ; }

template < class FT >
inline
VectorS2<FT> 
operator-(const PointS2<FT>& p, const PointS2<FT>& q)
{ return VectorS2<FT>(p.x() - q.x(), p.y() - q.y()) ; }

template < class FT >
inline
VectorS2<FT> 
operator-(const PointS2<FT>& p, const Origin& )
{ return VectorS2<FT>(p) ; }

template < class FT >
inline
VectorS2<FT> 
operator-(const Origin& , const PointS2<FT>& p)
{ return VectorS2<FT>(-p.x(), -p.y()) ; }

template < class FT >
CGAL_KERNEL_INLINE
PointS2<FT> 
PointS2<FT>::transform( const Aff_transformationS2<FT>& t) const
{ return t.transform(*this); }

template < class FT >
CGAL_KERNEL_INLINE
Bbox_2 
PointS2<FT>::bbox() const
{
  double bx = CGAL::to_double(x());
  double by = CGAL::to_double(y());
  return Bbox_2(bx,by, bx,by);
}

#ifndef CGAL_NO_OSTREAM_INSERT_POINTS2
template < class FT >
std::ostream& 
operator<<(std::ostream& os, const PointS2<FT>& p)
{
    switch(os.iword(IO::mode)) 
    {
    case IO::ASCII :
        return os << p.x() << ' ' << p.y();
    case IO::BINARY :
        write(os, p.x());
        write(os, p.y());
        return os;
    default:
        return os << "PointS2(" << p.x() << ", " << p.y() << ')';
    }
}
#endif // CGAL_NO_OSTREAM_INSERT_POINTS2

#ifndef CGAL_NO_ISTREAM_EXTRACT_POINTS2
template < class FT >
std::istream& 
operator>>(std::istream& is, PointS2<FT>& p)
{
    FT x, y;
    switch(is.iword(IO::mode)) 
    {
    case IO::ASCII :
        is >> x >> y;
        break;
    case IO::BINARY :
        read(is, x);
        read(is, y);
        break;
    default:
        CGAL_kernel_assertion_msg(false,"Stream must be in ascii or binary mode"); 
        // throw ios_base::failure("Stream must be in ascii or binary mode");
    }
    p = PointS2<FT>(x, y);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_POINTS2

CGAL_END_NAMESPACE

#endif // CGAL_POINTS2_H
