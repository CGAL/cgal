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
// release_date  : 2000, October 15
//
// source        : webS3/S3.lw
// file          : include/CGAL/SimpleCartesian/VectorS3.h
// package       : S3 (1.7)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// revision      : 1.7
// revision_date : 15 Oct 2000
// author(s)     : Stefan Schirra <Stefan.Schirra@@mpi-sb.mpg.de>
//                 based on code by
//                 Andreas Fabri and
//                 Herve Brönnimann
//
// coordinator   : MPI, Saarbrücken
// ======================================================================

#ifndef CGAL_VECTORS3_H
#define CGAL_VECTORS3_H

#include <CGAL/SimpleCartesian/PointS3.h>

CGAL_BEGIN_NAMESPACE

template < class FT >
inline VectorS3<FT>
operator-(const PointS3<FT>& p, const Origin& o);

template < class FT >
class VectorS3
{

friend class DirectionS3<FT>;

public:
                  VectorS3() {}
                  VectorS3(const Null_vector& )
                   : e0(FT(0)), e1(FT(0)), e2(FT(0)) {}
                  VectorS3(const FT& x, const FT& y, const FT& z)
                    : e0(x), e1(y), e2(z) {}
                  VectorS3(const FT& x, const FT& y, const FT& z, const FT& w);


  bool            operator==(const VectorS3<FT>& p) const;
  bool            operator!=(const VectorS3<FT>& p) const;

  bool            operator==(const Null_vector& ) const;
  bool            operator!=(const Null_vector& ) const;

  const FT&       x() const;
  const FT&       y() const;
  const FT&       z() const;
  const FT&       cartesian(int i) const;
  const FT&       operator[](int i) const;

  const FT&       hx() const;
  const FT&       hy() const;
  const FT&       hz() const;
  FT              hw() const;
  FT              homogeneous(int i) const;

  int             dimension() const;

  VectorS3<FT>    operator+(const VectorS3<FT>& w) const;
  VectorS3<FT>    operator-(const VectorS3<FT>& w) const;
  VectorS3<FT>    operator-() const;
  FT              operator*(const VectorS3<FT>& w) const;
  VectorS3<FT>    operator/(const FT& c) const;
  DirectionS3<FT> direction() const;
  VectorS3<FT> transform(const Aff_transformationS3<FT>& ) const;

// protected:
                  VectorS3(const PointS3<FT>& p);
                  VectorS3(const DirectionS3<FT>& p);

// private:
  FT   e0;
  FT   e1;
  FT   e2;
};


CGAL_END_NAMESPACE

#include <CGAL/SimpleCartesian/DirectionS3.h>

CGAL_BEGIN_NAMESPACE


template < class FT >
VectorS3<FT>::VectorS3(const FT& x, const FT& y, const FT& z, const FT& w)
{
  if (w != FT(1))
  {
    e0 = x/w; 
    e1 = y/w; 
    e2 = z/w;
  }
  else
  {
    e0 = x; 
    e1 = y; 
    e2 = z;
  }
}

template < class FT >
inline VectorS3<FT>::VectorS3(const PointS3<FT>& p)
{
  e0 = p.e0; 
  e1 = p.e1; 
  e2 = p.e2;
}

template < class FT >
inline VectorS3<FT>::VectorS3(const DirectionS3<FT>& d)
{
  e0 = d.e0; 
  e1 = d.e1; 
  e2 = d.e2;
}

template < class FT >
bool 
VectorS3<FT>::operator==(const VectorS3<FT>& v) const
{ return (x() == v.x()) && (y() == v.y()) && (z() == v.z()) ; }

template < class FT >
inline 
bool 
VectorS3<FT>::operator!=(const VectorS3<FT>& v) const
{ return !(*this == v); }


template < class FT >
bool 
VectorS3<FT>::operator==(const Null_vector& ) const
{ return (x() == FT(0)) && (y() == FT(0)) && (z() == FT(0)) ; }

template < class FT >
inline 
bool 
VectorS3<FT>::operator!=(const Null_vector& v) const
{ return !(*this == v); }

template < class FT >
inline
const FT&
VectorS3<FT>::x()  const
{ return e0; }

template < class FT >
inline
const FT&
VectorS3<FT>::y()  const
{ return  e1; }

template < class FT >
inline
const FT&
VectorS3<FT>::z()  const
{ return  e2; }

template < class FT >
inline
const FT&
VectorS3<FT>::cartesian(int i) const
{
  CGAL_kernel_precondition( (i>=0) && (i<3) );
  return (i==0) ? x() :
         (i==1) ? y() : z();
}

template < class FT >
inline
const FT&
VectorS3<FT>::operator[](int i) const
{ return cartesian(i); }

template < class FT >
inline
int
VectorS3<FT>::dimension() const
{ return 3; }

template < class FT >
inline
const FT&
VectorS3<FT>::hx()  const
{ return e0; }

template < class FT >
inline
const FT&
VectorS3<FT>::hy()  const
{ return e1; }

template < class FT >
inline
const FT&
VectorS3<FT>::hz()  const
{ return e2; }

template < class FT >
FT
VectorS3<FT>::hw()  const
{ return FT(1); }

template < class FT >
FT
VectorS3<FT>::homogeneous(int i) const
{ return (i==3) ? FT(1) : cartesian(i); }

template < class FT >
inline 
VectorS3<FT> 
VectorS3<FT>::operator+(const VectorS3<FT>& w) const
{ return VectorS3<FT>(x() + w.x(), y() + w.y(), z() + w.z()) ; }

template < class FT >
inline 
VectorS3<FT> 
VectorS3<FT>::operator-(const VectorS3<FT>& w) const
{ return VectorS3<FT>(x() - w.x(), y() - w.y(), z() - w.z()) ; }

template < class FT >
inline 
VectorS3<FT> 
VectorS3<FT>::operator-() const
{ return VectorS3<FT>(-x(), -y(), -z()) ; }

template < class FT >
inline 
FT 
VectorS3<FT>::operator*(const VectorS3<FT>& w) const
{ return x() * w.x() + y() * w.y() + z() * w.z() ; }

template < class FT >
inline 
VectorS3<FT> 
operator*(const FT& c, const VectorS3<FT>& w)
{ return VectorS3<FT>( c* w.x(), c * w.y(), c * w.z()) ; }

template < class FT >
inline 
VectorS3<FT> 
VectorS3<FT>::operator/(const FT& c) const
{ return VectorS3<FT>( x()/c, y()/c, z()/c) ; }

template < class FT >
VectorS3<FT> 
cross_product(const VectorS3<FT>& v, const VectorS3<FT>& w)
{
    return VectorS3<FT>( v.y() * w.z() - v.z() * w.y() ,
                         v.z() * w.x() - v.x() * w.z() ,
                         v.x() * w.y() - v.y() * w.x() );
}

template < class FT >
inline
DirectionS3<FT>
VectorS3<FT>::direction() const
{ return DirectionS3<FT>(*this); }


template < class FT >
VectorS3<FT>
VectorS3<FT>::transform(const Aff_transformationS3<FT>& t) const
{ return t.transform(*this); }


#ifndef CGAL_NO_OSTREAM_INSERT_VECTORS3
template < class FT >
std::ostream& operator<<(std::ostream& os, const VectorS3<FT>& v)
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
        os << "VectorS3(" << v.x() << ", " << v.y() <<  ", " << v.z() << ")";
        return os;
    }
}
#endif // CGAL_NO_OSTREAM_INSERT_VECTORS3

#ifndef CGAL_NO_ISTREAM_EXTRACT_VECTORS3
template < class FT >
std::istream& operator>>(std::istream& is, VectorS3<FT>& p)
{
    FT x, y, z;
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
        CGAL_kernel_assertion_msg(false,"Stream must be in ascii or binary mode"); 
        // throw ios_base::failure("Stream must be in ascii or binary mode");
        break;
    }
    p = VectorS3<FT>(x, y, z);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_VECTORS3


CGAL_END_NAMESPACE

#endif
