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
// file          : include/CGAL/SimpleCartesian/DirectionS3.h
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

#ifndef CGAL_DIRECTIONS3_H
#define CGAL_DIRECTIONS3_H

#include <CGAL/SimpleCartesian/VectorS3.h>

CGAL_BEGIN_NAMESPACE

template < class FT >
class DirectionS3
{
public:
                 DirectionS3() {}
                 DirectionS3(const VectorS3<FT>& v);
                 DirectionS3(const FT& x, const FT& y, const FT& z)
                  : e0(x), e1(y), e2(z) {}

  bool           operator==(const DirectionS3<FT>& d) const;
  bool           operator!=(const DirectionS3<FT>& d) const;

  VectorS3<FT>   to_vector() const;
  VectorS3<FT>   vector() const { return to_vector(); } 


  DirectionS3    transform(const Aff_transformationS3<FT>& t) const;

  DirectionS3    operator-() const;

  const FT&      delta(int i) const;
  const FT&      dx() const;
  const FT&      dy() const;
  const FT&      dz() const;

  const FT&      hdx() const;
  const FT&      hdy() const;
  const FT&      hdz() const;
  FT             hw() const;


// private:
  FT   e0;
  FT   e1;
  FT   e2;
};

template < class FT >
DirectionS3<FT>::DirectionS3(const VectorS3<FT>& v)
{
  e0 = v.e0; 
  e1 = v.e1; 
  e2 = v.e2;
}

template < class FT >
bool 
DirectionS3<FT>::operator==(const DirectionS3<FT>& d) const
{
  return  ( dx()*d.dy() == dy()*d.dx() )
        &&( dx()*d.dz() == dz()*d.dx() )
        &&( dy()*d.dz() == dz()*d.dy() )
        &&( CGAL_NTS sign( dx() ) == CGAL_NTS sign( d.dx() ) )
        &&( CGAL_NTS sign( dy() ) == CGAL_NTS sign( d.dy() ) )
        &&( CGAL_NTS sign( dz() ) == CGAL_NTS sign( d.dz() ) );
}

template < class FT >
inline 
bool  
DirectionS3<FT>::operator!=(const DirectionS3<FT>& d) const
{ return !(*this == d); }

template < class FT >
inline 
VectorS3<FT> 
DirectionS3<FT>::to_vector() const
{ return VectorS3<FT>(*this); }


template < class FT >
inline
DirectionS3<FT>
DirectionS3<FT>::transform(const Aff_transformationS3<FT>& t) const
{ return t.transform(*this); }


template < class FT >
inline 
DirectionS3<FT> 
DirectionS3<FT>::operator-() const
{ return DirectionS3<FT>(-dx(), -dy(), -dz()); }


template < class FT >
const FT&
DirectionS3<FT>::delta(int i) const
{
  CGAL_kernel_precondition( i >= 0 && i <= 2 );
  return (i==0) ? dx() :
         (i==1) ? dy() : dz() ;
}


template < class FT >
inline 
const FT&
DirectionS3<FT>::dx() const
{ return e0; }


template < class FT >
inline 
const FT&
DirectionS3<FT>::dy() const
{ return e1; }


template < class FT >
inline 
const FT&
DirectionS3<FT>::dz() const
{ return e2; } 

template < class FT >
inline 
const FT&
DirectionS3<FT>::hdx() const
{ return e0; }


template < class FT >
inline 
const FT&
DirectionS3<FT>::hdy() const
{ return e1; }


template < class FT >
inline 
const FT&
DirectionS3<FT>::hdz() const
{ return e2; }

template < class FT >
inline 
FT 
DirectionS3<FT>::hw() const
{ return FT(1); }



#ifndef CGAL_NO_OSTREAM_INSERT_DIRECTIONS3
template < class FT >
std::ostream& operator<<(std::ostream& os, const DirectionS3<FT>& d)
{
  VectorS3<FT> v = d.vector();
  switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << v.x() << ' ' << v.y()  << ' ' << v.z();
    case IO::BINARY :
        write(os, v.x());
        write(os, v.y());
        write(os, v.z());
        return os;
    default:
        os << "DirectionS3(" << v.x() << ", " << v.y() << ", " << v.z() << ")";
        return os;
    }
}
#endif // CGAL_NO_OSTREAM_INSERT_DIRECTIONS3

#ifndef CGAL_NO_ISTREAM_EXTRACT_DIRECTIONS3
template < class FT >
std::istream& operator>>(std::istream& is, DirectionS3<FT>& p)
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
    p = DirectionS3<FT>(x, y, z);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_DIRECTIONS3



CGAL_END_NAMESPACE

#endif
