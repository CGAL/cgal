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
// file          : include/CGAL/SimpleCartesian/RayS2.h
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

#ifndef CGAL_RAYS2_H
#define CGAL_RAYS2_H

#include <CGAL/SimpleCartesian/PointS2.h>
#include <CGAL/SimpleCartesian/LineS2.h>

CGAL_BEGIN_NAMESPACE

template < class FT >
class RayS2
{
public:
                  RayS2() {};
                  RayS2(const PointS2<FT>& sp,
                        const PointS2<FT>& secondp) : s(sp), t(secondp) {};
                  RayS2(const PointS2<FT>& sp,
                        const DirectionS2<FT>& d);


  bool            operator==(const RayS2<FT>& r) const;
  bool            operator!=(const RayS2<FT>& r) const;

  const PointS2<FT>&     start() const;
  const PointS2<FT>&     source() const;

  PointS2<FT>     point(int i) const;
  PointS2<FT>     second_point() const;

  DirectionS2<FT> direction() const;
  LineS2<FT>      supporting_line() const;
  RayS2<FT>       opposite() const;

  RayS2<FT>       transform(const Aff_transformationS2<FT>& t) const;
 
  bool             is_horizontal() const;
  bool             is_vertical() const;
  bool             is_degenerate() const;
  bool             has_on(const PointS2<FT>& p) const;
  bool             collinear_has_on(const PointS2<FT>& p) const;

// private:
  PointS2<FT>     s;
  PointS2<FT>     t;
};


template < class FT >
CGAL_KERNEL_CTOR_INLINE
RayS2<FT>::RayS2(const PointS2<FT>& sp, const DirectionS2<FT> &d)
{ s = sp; t = sp + d.vector(); }

template < class FT >
CGAL_KERNEL_INLINE
bool 
RayS2<FT>::operator==(const RayS2<FT>& r) const
{ return ((source() == r.source()) && (direction() == r.direction()) ); }

template < class FT >
bool 
RayS2<FT>::operator!=(const RayS2<FT>& r) const
{ return !(*this == r); }

template < class FT >
inline
const PointS2<FT>&   
RayS2<FT>::start() const
{ return s; }

template < class FT >
inline
const PointS2<FT>&   
RayS2<FT>::source() const
{ return s; }

template < class FT >
inline
PointS2<FT>  
RayS2<FT>::second_point() const
{ return t; }

template < class FT >
CGAL_KERNEL_INLINE
PointS2<FT>  
RayS2<FT>::point(int i) const
{
  CGAL_kernel_precondition( i >= 0 );
  if (i == 0)
    return s;

  if (i == 1)
    return t;

  return source() + FT(i) * (second_point() - source());
}

template < class FT >
inline
DirectionS2<FT> 
RayS2<FT>::direction() const
{ return DirectionS2<FT>(  second_point() - source() ); }

template < class FT >
inline
LineS2<FT> RayS2<FT>::supporting_line() const
{ return LineS2<FT>(*this); }

template < class FT >
inline
RayS2<FT> RayS2<FT>::opposite() const
{ return RayS2<FT>( source(), - direction() ); }


template < class FT >
CGAL_KERNEL_INLINE
RayS2<FT> 
RayS2<FT>::transform(const Aff_transformationS2<FT>& t) const
{ return RayS2<FT>(t.transform(source()), t.transform(second_point())); }


template < class FT >
CGAL_KERNEL_INLINE
bool 
RayS2<FT>::is_horizontal() const
{ return (source().y() ==  second_point().y()); }

template < class FT >
CGAL_KERNEL_INLINE
bool 
RayS2<FT>::is_vertical() const
{ return  (source().x() == second_point().x()); }

template < class FT >
CGAL_KERNEL_INLINE
bool 
RayS2<FT>::is_degenerate() const
{ return (source() == second_point()); }

template < class FT >
CGAL_KERNEL_INLINE
bool 
RayS2<FT>::has_on(const PointS2<FT>& p) const
{
  return ( p == source()
           || ( collinear(source(), p, second_point())
          && ( DirectionS2<FT>(p - source()) == direction() )));
}

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
bool 
RayS2<FT>::collinear_has_on(const PointS2<FT>& p) const
{
    switch(compare_x(source(), second_point())){
    case SMALLER:
        return compare_x(source(), p) != LARGER;
    case LARGER:
        return compare_x(p, source()) != LARGER;
    default:
        switch(compare_y(source(), second_point())){
        case SMALLER:
            return compare_y(source(), p) != LARGER;
        case LARGER:
            return compare_y(p, source()) != LARGER;
        default:
            return true; // p == source()
        }
    }
}


#ifndef CGAL_NO_OSTREAM_INSERT_RAYS2
template < class FT >
std::ostream& operator<<(std::ostream &os, const RayS2<FT> &r)
{
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << r.source() << ' ' << r.direction();
    case IO::BINARY :
        return os << r.source() << r.direction();
    default:
        return os << "RayS2(" << r.source() <<  ", " << r.direction() << ")";
    }
}
#endif // CGAL_NO_OSTREAM_INSERT_RAYS2

#ifndef CGAL_NO_ISTREAM_EXTRACT_RAYS2
template < class FT >
std::istream& operator>>(std::istream &is, RayS2<FT> &r)
{
    PointS2<FT> p;
    DirectionS2<FT> d;

    is >> p >> d;

    r = RayS2<FT>(p, d);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_RAYS2



CGAL_END_NAMESPACE

#endif
