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
// file          : include/CGAL/SimpleCartesian/PointS3.h
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

#ifndef CGAL_POINTS3_H
#define CGAL_POINTS3_H

CGAL_BEGIN_NAMESPACE

template < class FT >
inline
PointS3<FT>
operator+(const Origin& o, const VectorS3<FT>& v);

template < class FT >
inline
PointS3<FT>
operator-(const Origin& o, const VectorS3<FT>& v);


template < class FT >
class PointS3
{
friend CGAL_FRIEND_INLINE 
       PointS3<FT> 
       operator+ CGAL_NULL_TMPL_ARGS( const Origin& o, const VectorS3<FT>& v);

friend CGAL_FRIEND_INLINE 
       PointS3<FT> 
       operator- CGAL_NULL_TMPL_ARGS( const Origin& o, const VectorS3<FT>& v);
public:
              PointS3() {}
              PointS3(const Origin& o) 
               : e0(FT(0)), e1(FT(0)), e2(FT(0)) {}
              PointS3(const FT& x, const FT& y, const FT& z)
               : e0(x), e1(y), e2(z) {}
              PointS3(const FT& x, const FT& y, const FT& z, const FT& hw);

  bool        operator==(const PointS3<FT>& p) const;
  bool        operator!=(const PointS3<FT>& p) const;

  const FT&   x() const;
  const FT&   y() const;
  const FT&   z() const;

  const FT&   hx() const;
  const FT&   hy() const;
  const FT&   hz() const;
  FT          hw() const;

  const FT&   cartesian(int i) const;
  const FT&   operator[](int i) const;

  FT          homogeneous(int i) const;

  int         dimension() const;
  Bbox_3      bbox() const;

  PointS3<FT> transform( const Aff_transformationS3<FT>& ) const;


// protected:
              PointS3(const VectorS3<FT>& v);
// private:
  FT   e0;
  FT   e1;
  FT   e2;
};

CGAL_END_NAMESPACE

#include <CGAL/Origin.h>
#include <CGAL/SimpleCartesian/VectorS3.h>
#include <CGAL/SimpleCartesian/Aff_transformationS3.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/number_utils.h>

CGAL_BEGIN_NAMESPACE


template < class FT >
PointS3<FT>::PointS3(const FT& x, const FT& y, const FT& z, const FT& w)
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
PointS3<FT>::PointS3(const VectorS3<FT>& v)
{
  e0 = v.e0; 
  e1 = v.e1; 
  e2 = v.e2;
}

template < class FT >
inline
bool
PointS3<FT>::operator==(const PointS3<FT>& p) const
{ return (x() == p.x()) && (y() == p.y()) && (z() == p.z()); }

template < class FT >
inline
bool
PointS3<FT>::operator!=(const PointS3<FT>& p) const
{ return !(*this == p); }


template < class FT >
inline
const FT&
PointS3<FT>::x()  const
{ return e0; }


template < class FT >
inline
const FT&
PointS3<FT>::y()  const
{ return  e1; }


template < class FT >
inline
const FT&
PointS3<FT>::z()  const
{ return  e2; }


template < class FT >
inline
int
PointS3<FT>::dimension() const
{ return 3; }


template < class FT >
inline
const FT&
PointS3<FT>::cartesian(int i) const
{
  CGAL_kernel_precondition( (i>=0) && (i<=2) );
  return (i==0) ? x() :
         (i==1) ? y() : z();
}


template < class FT >
inline
const FT&
PointS3<FT>::operator[](int i) const
{
  return cartesian(i);
}


template < class FT >
inline
const FT&
PointS3<FT>::hx()  const
{ return e0; }


template < class FT >
inline
const FT&
PointS3<FT>::hy()  const
{ return e1; }


template < class FT >
inline
const FT&
PointS3<FT>::hz()  const
{ return e2; }


template < class FT >
inline
FT
PointS3<FT>::hw()  const
{ return FT(1); }

template < class FT >
FT
PointS3<FT>::homogeneous(int i) const
{
  CGAL_kernel_precondition((i>=0) && (i<=3));
  return (i<3) ? cartesian(i) : FT(1);
}

template < class FT >
inline
PointS3<FT>
operator+(const PointS3<FT>& p, const VectorS3<FT>& v)
{ return PointS3<FT>(p.x() + v.x(), p.y() + v.y(), p.z() + v.z()); }

template < class FT >
inline
PointS3<FT>
operator-(const PointS3<FT>& p, const VectorS3<FT>& v)
{ return PointS3<FT>(p.x() - v.x(), p.y() - v.y(), p.z() - v.z()); }

template < class FT >
inline
PointS3<FT>
operator+(const Origin& , const VectorS3<FT>& v)
{ return PointS3<FT>(v); }

template < class FT >
inline
PointS3<FT>
operator-(const Origin& , const VectorS3<FT>& v)
{ return PointS3<FT>(-v); }

template < class FT >
inline
VectorS3<FT>
operator-(const PointS3<FT>& p, const PointS3<FT>& q)
{ return VectorS3<FT>(p.x() - q.x(), p.y() - q.y(), p.z() - q.z()); }


template < class FT >
inline
VectorS3<FT>
operator-(const PointS3<FT>& p, const Origin& )
{ return VectorS3<FT>(p); }


template < class FT >
inline
VectorS3<FT>
operator-(const Origin& , const PointS3<FT>& p)
{ return VectorS3<FT>(-p.x(), -p.y(), -p.z()); }


template < class FT >
inline
PointS3<FT>
PointS3<FT>::transform(const Aff_transformationS3<FT>& t) const
{ return t.transform(*this); }


template < class FT >
Bbox_3 PointS3<FT>::bbox() const
{
  // Not robust...
  double bx = CGAL::to_double(x());
  double by = CGAL::to_double(y());
  double bz = CGAL::to_double(z());
  return Bbox_3(bx, by, bz, bx, by, bz);
}


#ifndef CGAL_NO_OSTREAM_INSERT_POINTS3
template < class FT >
std::ostream& operator<<(std::ostream& os, const PointS3<FT>& p)
{
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << p.x() << ' ' << p.y()  << ' ' << p.z();
    case IO::BINARY :
        write(os, p.x());
        write(os, p.y());
        write(os, p.z());
        return os;
    default:
        os << "PointS3(" << p.x() << ", " << p.y()  << ", " << p.z() << ")";
        return os;
    }
}
#endif // CGAL_NO_OSTREAM_INSERT_POINTS3

#ifndef CGAL_NO_ISTREAM_EXTRACT_POINTS3
template < class FT >
std::istream& operator>>(std::istream& is, PointS3<FT>& p)
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
        CGAL_kernel_assertion_msg(false,"Stream must be in ascii or binary mode"
); 
        // throw ios_base::failure("Stream must be in ascii or binary mode");
        break;
    }
    p = PointS3<FT>(x, y, z);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_POINTS3


CGAL_END_NAMESPACE

#endif // CGAL_POINTS3_H
