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
// file          : include/CGAL/SimpleCartesian/DirectionS2.h
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

#ifndef CGAL_DIRECTIONS2_H
#define CGAL_DIRECTIONS2_H

#include <CGAL/SimpleCartesian/VectorS2.h>

CGAL_BEGIN_NAMESPACE

template < class FT >
class DirectionS2
{
public:
                  DirectionS2();
                  DirectionS2(const VectorS2<FT>& v);
                  DirectionS2(const FT& x, const FT &y);

  bool            operator==(const DirectionS2<FT>& d) const;
  bool            operator!=(const DirectionS2<FT>& d) const;
  bool            operator>=(const DirectionS2<FT>& d) const;
  bool            operator<=(const DirectionS2<FT>& d) const;
  bool            operator>(const DirectionS2<FT>& d) const;
  bool            operator<(const DirectionS2<FT>& d) const;
  bool            counterclockwise_in_between( const DirectionS2<FT>& d1,
                                               const DirectionS2<FT>& d2) const;

  VectorS2<FT>    vector() const;

  DirectionS2<FT> perpendicular(const Orientation& o) const;
  DirectionS2<FT> transform(const Aff_transformationS2<FT>& t) const;

  DirectionS2<FT> operator-() const;

  FT              delta(int i) const;
  FT              dx() const;
  FT              dy() const;

// private:
  FT  e0;
  FT  e1;
};


template < class FT >
CGAL_KERNEL_CTOR_INLINE
DirectionS2<FT>::DirectionS2() {}

template < class FT >
CGAL_KERNEL_CTOR_INLINE
DirectionS2<FT>::DirectionS2(const VectorS2<FT>& v) 
{ e0 = v.e0; e1 = v.e1; }

template < class FT >
CGAL_KERNEL_CTOR_INLINE
DirectionS2<FT>::DirectionS2(const FT& x, const FT &y)
{ e0 = x; e1 = y; }

template < class FT >
bool
DirectionS2<FT>::operator==(const DirectionS2<FT>& d) const
{
// Use a S2 predicate for that ?
    return (CGAL_NTS sign(dx()) == CGAL_NTS sign(d.dx()))
        && (CGAL_NTS sign(dy()) == CGAL_NTS sign(d.dy()))
        && (dy()*d.dx() == d.dy()*dx());
}

template < class FT >
inline
bool
DirectionS2<FT>::operator!=(const DirectionS2<FT>& d) const
{ return !( *this == d ); }

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
bool
DirectionS2<FT>::operator<(const DirectionS2<FT>& d) const
{
  int quadrant_this = (dx() >= FT(0)) ? ((dy() >= FT(0))?1:4)
                                      : ((dy() >= FT(0))?2:3);
  int quadrant_d    = (d.dx() >= FT(0)) ? ((d.dy() >= FT(0))?1:4)
                                        : ((d.dy() >= FT(0))?2:3);

  if(quadrant_this < quadrant_d)
    return true;
  else if (quadrant_this > quadrant_d)
    return false;
  else
    return dy() * d.dx() < d.dy() * dx();
}

template < class FT >
CGAL_KERNEL_INLINE
bool
DirectionS2<FT>::operator>(const DirectionS2<FT>& d) const
{ return d < *this ; }

template < class FT >
CGAL_KERNEL_INLINE
bool
DirectionS2<FT>::operator>=(const DirectionS2<FT>& d) const
{ return (d < *this) || (d == *this) ; }

template < class FT >
CGAL_KERNEL_INLINE
bool
DirectionS2<FT>::operator<=(const DirectionS2<FT>& d) const
{ return (*this < d) || (d == *this) ; }

template < class FT >
CGAL_KERNEL_INLINE
bool
DirectionS2<FT>::counterclockwise_in_between(const DirectionS2<FT>& d1,
                                             const DirectionS2<FT>& d2) const
{
  if ( d1 < *this)
  {
      return ( *this < d2 )||( d2 <= d1 );
  }
  else
  {
      return ( *this < d2 )&&( d2 <= d1 );
  }
}

template < class FT >
inline
VectorS2<FT>
DirectionS2<FT>::vector() const
{ return VectorS2<FT>(*this); }

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
DirectionS2<FT>
DirectionS2<FT>::perpendicular(const Orientation& o) const
{
  CGAL_kernel_precondition(o != COLLINEAR);
  if (o == COUNTERCLOCKWISE)
    return DirectionS2<FT>(-dy(), dx());
  else
    return DirectionS2<FT>(dy(), -dx());
}

template < class FT >
CGAL_KERNEL_INLINE
DirectionS2<FT> DirectionS2<FT>::transform( const Aff_transformationS2<FT>& t) const
{ return t.transform(*this); }

template < class FT >
inline
DirectionS2<FT> DirectionS2<FT>::operator-() const
{ return DirectionS2<FT>(-dx(), -dy()); }



template < class FT >
CGAL_KERNEL_INLINE
FT
DirectionS2<FT>::delta(int i) const
{
  CGAL_kernel_precondition( ( i == 0 ) || ( i == 1 ) );
  return (i==0) ? dx() : dy();
}


template < class FT >
inline
FT
DirectionS2<FT>::dx() const
{ return e0; }

template < class FT >
inline
FT
DirectionS2<FT>::dy() const
{ return e1; }


#ifndef CGAL_NO_OSTREAM_INSERT_DIRECTIONS2
template < class FT >
std::ostream
&operator<<(std::ostream& os, const DirectionS2<FT> &d)
{
    VectorS2<FT> v = d.vector();
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << v.x() << ' ' << v.y();
    case IO::BINARY :
        write(os, v.x());
        write(os, v.y());
        return os;
    default:
        return os << "DirectionS2(" << v.x() << ", " << v.y() << ')';
    }
}
#endif // CGAL_NO_OSTREAM_INSERT_DIRECTIONS2

#ifndef CGAL_NO_ISTREAM_EXTRACT_DIRECTIONS2

template < class FT >
std::istream
&operator>>(std::istream& is, DirectionS2<FT> &p)
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
        std::cerr << std::endl << "Stream must be in ascii or binary mode" << std::endl;
        break;
    }
    p = DirectionS2<FT>(x, y);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_DIRECTIONS2



CGAL_END_NAMESPACE

#endif
