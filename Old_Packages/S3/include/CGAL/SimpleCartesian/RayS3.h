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
// file          : include/CGAL/SimpleCartesian/RayS3.h
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

#ifndef CGAL_RAYS3_H
#define CGAL_RAYS3_H

#include <CGAL/SimpleCartesian/LineS3.h>

CGAL_BEGIN_NAMESPACE

template < class FT >
class RayS3
{
public:
                  RayS3() {}
                  RayS3(const PointS3<FT>& sp, const PointS3<FT>& secondp);
                  RayS3(const PointS3<FT>& sp, const DirectionS3<FT>& d);

  bool            operator==(const RayS3<FT>& r) const;
  bool            operator!=(const RayS3<FT>& r) const;

  PointS3<FT>     start() const;
  PointS3<FT>     source() const;
  PointS3<FT>     second_point() const;
  PointS3<FT>     point(int i) const;

  DirectionS3<FT> direction() const;
  LineS3<FT>      supporting_line() const;
  RayS3           opposite() const;

  RayS3           transform(const Aff_transformationS3<FT>& t) const;

  bool            is_degenerate() const;
  bool            has_on(const PointS3<FT>& p) const;
  bool            collinear_has_on(const PointS3<FT>& p) const;

// private:
  PointS3<FT>     e0;
  PointS3<FT>     e1;
};


template < class FT >
RayS3<FT>::RayS3(const PointS3<FT>& sp, const PointS3<FT>& secondp)
{
  e0 = sp;
  e1 = secondp;
}


template < class FT >
RayS3<FT>::RayS3(const PointS3<FT>& sp, const DirectionS3<FT>& d)
{
  e0 = sp;
  e1 = sp + d.vector();
}


template < class FT >
inline 
bool 
RayS3<FT>::operator==(const RayS3<FT>& r) const
{ return (source() == r.source()) && (direction() == r.direction()); }


template < class FT >
inline 
bool 
RayS3<FT>::operator!=(const RayS3<FT>& r) const
{ return !(*this == r); }


template < class FT >
PointS3<FT>  
RayS3<FT>::start() const
{ return e0; }


template < class FT >
PointS3<FT>  
RayS3<FT>::source() const
{ return e0; }

template < class FT >
PointS3<FT>  
RayS3<FT>::second_point() const
{ return e1; }


template < class FT >
PointS3<FT>  
RayS3<FT>::point(int i) const
{
  CGAL_kernel_precondition( i >= 0 );
  if (i == 0)
    return e0;

  if (i == 1)
    return e1;

  return source() + FT(i) * (second_point() - source());
}


template < class FT >
inline
DirectionS3<FT>
RayS3<FT>::direction() const
{ return DirectionS3<FT>( second_point() - source() ); }


template < class FT >
inline
LineS3<FT>
RayS3<FT>::supporting_line() const
{ return LineS3<FT>(*this); }


template < class FT >
inline
RayS3<FT>
RayS3<FT>::opposite() const
{ return RayS3<FT>( source(), - direction() ); }


template < class FT >
inline
RayS3<FT>
RayS3<FT>::transform(const Aff_transformationS3<FT>& t) const
{ return RayS3<FT>(t.transform(source()), t.transform(second_point())); }


template < class FT >
bool
RayS3<FT>::has_on(const PointS3<FT>& p) const
{
  return (p == source()) ||
         ( collinear(source(), p, second_point())
           && ( DirectionS3<FT>(p - source()) == direction() ));
}


template < class FT >
inline
bool
RayS3<FT>::is_degenerate() const
{ return source() == second_point(); }


template < class FT >
inline
bool
RayS3<FT>::collinear_has_on(const PointS3<FT>& p) const
{
  CGAL_kernel_exactness_precondition( collinear(source(), p, second_point()) );

  Comparison_result cx = compare_x(source(), second_point());
  if (cx != EQUAL)
    return cx != compare_x(p, source());

  Comparison_result cy = compare_y(source(), second_point());
  if (cy != EQUAL)
    return cy != compare_y(p, source());

  Comparison_result cz = compare_z(source(), second_point());
  if (cz != EQUAL)
    return cz != compare_z(p, source());

  return true; // p == source()
}


#ifndef CGAL_NO_OSTREAM_INSERT_RAYS3
template < class FT >
std::ostream& 
operator<<(std::ostream& os, const RayS3<FT>& r)
{
    switch(os.iword(IO::mode)) 
    {
      case IO::ASCII :
        return os << r.start() << ' ' << r.direction();
      case IO::BINARY :
        return os<< r.start() << r.direction();
      default:
        return os << "RayS3(" << r.start() <<  ", " << r.direction() << ")";
    }
}
#endif // CGAL_NO_OSTREAM_INSERT_RAYS3

#ifndef CGAL_NO_ISTREAM_EXTRACT_RAYS3
template < class FT >
std::istream& 
operator>>(std::istream& is, RayS3<FT>& r)
{
    PointS3<FT> p;
    DirectionS3<FT> d;

    is >> p >> d;

    r = RayS3<FT>(p, d);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_RAYS3


CGAL_END_NAMESPACE

#endif
