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
// file          : include/CGAL/SimpleCartesian/LineS2.h
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


#ifndef CGAL_LINES2_H
#define CGAL_LINES2_H

#include <CGAL/SimpleCartesian/SegmentS2.h>
#include <CGAL/SimpleCartesian/RayS2.h>
#include <CGAL/SimpleCartesian/predicates_on_pointsS2.h>

CGAL_BEGIN_NAMESPACE

template < class FT >
class LineS2
{
public:
                  LineS2();
                  LineS2(const PointS2<FT>& p,
                         const PointS2<FT>& q);
                  LineS2(const FT& a, const FT &b, const FT &c);
                  LineS2(const SegmentS2<FT>& s);
                  LineS2(const RayS2<FT>& r);
                  LineS2(const PointS2<FT>& p,
                         const DirectionS2<FT>& d);

  bool            operator==(const LineS2<FT>& l) const;
  bool            operator!=(const LineS2<FT>& l) const;

  FT              a() const;
  FT              b() const;
  FT              c() const;

  FT              x_at_y(const FT& y) const;
  FT              y_at_x(const FT& x) const;

  LineS2<FT>     perpendicular(const PointS2<FT>& p) const;
  LineS2<FT>     opposite() const;
  PointS2<FT>    point(int i) const;

  PointS2<FT>    point() const;
  PointS2<FT>    projection(const PointS2<FT>& p) const;

  DirectionS2<FT> direction() const;

  Oriented_side   oriented_side(const PointS2<FT>& p) const;
  bool            has_on_boundary(const PointS2<FT>& p) const;
  bool            has_on_positive_side(const PointS2<FT>& p) const;
  bool            has_on_negative_side(const PointS2<FT>& p) const;

  bool            is_horizontal() const;
  bool            is_vertical() const;
  bool            is_degenerate() const;

  LineS2<FT>     transform(const Aff_transformationS2<FT>& t) const;

// private:
  void            new_rep(const PointS2<FT>& p, const PointS2<FT> &q);
  void            new_rep(const FT& a, const FT &b, const FT &c);

  FT e0;
  FT e1;
  FT e2;

};


template < class FT >
inline
LineS2<FT>::LineS2()
{}

template < class FT >
CGAL_KERNEL_INLINE
void LineS2<FT>::new_rep(const PointS2<FT>& p, const PointS2<FT> &q)
{
  e0 = p.y() - q.y();
  e1 = q.x() - p.x();
  e2 = p.x()*q.y() - p.y()*q.x();
}

template < class FT >
CGAL_KERNEL_INLINE
void LineS2<FT>::new_rep(const FT& a, const FT &b, const FT &c)
{
  e0 = a;
  e1 = b;
  e2 = c;
}

template < class FT >
inline
LineS2<FT>::LineS2(const PointS2<FT>& p, const PointS2<FT> &q)
{ new_rep(p,q); }

template < class FT >
inline
LineS2<FT>::LineS2(const FT& a, const FT &b, const FT &c)
{ new_rep(a,b,c); }

template < class FT >
inline
LineS2<FT>::LineS2(const SegmentS2<FT>& s)
{ new_rep( s.start(), s.end()); }

template < class FT >
inline
LineS2<FT>::LineS2(const RayS2<FT>& r)
{ new_rep(r.start(), r.second_point()); }

template < class FT >
CGAL_KERNEL_INLINE
LineS2<FT>::LineS2(const PointS2<FT>& p, const DirectionS2<FT> &d)
{ new_rep(-d.dy(), d.dx(), -d.dx()* p.y()  + d.dy() * p.x()); }


template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
bool LineS2<FT>::operator==(const LineS2<FT>& l) const
{
  if (  (a() * l.c() != l.a() * c())
      ||(b() * l.c() != l.b() * c()) )
      return false;
  int sc  = CGAL_NTS sign(c());
  int slc = CGAL_NTS sign(l.c());
  if ( sc == slc )
      return (sc == 0) ?  ( a()*l.b() ==  b()*l.a() )
                         && (CGAL_NTS sign(a() ) == CGAL_NTS sign( l.a() ))
                         && (CGAL_NTS sign(b() ) == CGAL_NTS sign( l.b() ))
                       : true;
  return false;
}

template < class FT >
inline
bool LineS2<FT>::operator!=(const LineS2<FT>& l) const
{ return !(*this == l); }

template < class FT >
inline
FT LineS2<FT>::a() const
{ return e0; }

template < class FT >
inline
FT LineS2<FT>::b() const
{ return e1; }

template < class FT >
inline
FT LineS2<FT>::c() const
{ return e2; }

template < class FT >
CGAL_KERNEL_INLINE
FT LineS2<FT>::x_at_y(const FT& y) const
{
  CGAL_kernel_precondition_msg( (a() != FT(0)),
               "Line::x_at_y(const FT& y) is undefined for horizontal line" );
  return ( -b()*y - c() ) / a();
}

template < class FT >
CGAL_KERNEL_INLINE
FT LineS2<FT>::y_at_x(const FT& x) const
{
  CGAL_kernel_precondition_msg( (b() != FT(0)),
              "Line::x_at_y(const FT& y) is undefined for vertical line");
  return ( -a()*x - c() ) / b();
}

template < class FT >
inline
LineS2<FT> LineS2<FT>::perpendicular(const PointS2<FT>& p) const
{ return LineS2<FT>( -b() , a() , b() * p.x() - a() * p.y()  ); }

template < class FT >
inline
LineS2<FT> LineS2<FT>::opposite() const
{ return LineS2<FT>( -a(), -b(), -c() ); }

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
PointS2<FT> LineS2<FT>::point(int i) const
{
  if (i == 0)
    return is_vertical() ? PointS2<FT>( (-b()-c())/a(), FT(1) )
                         : PointS2<FT>( FT(1), -(a()+c())/b());
  if (i == 1)
    return is_vertical() ? PointS2<FT>( (-b()-c())/a() + b(), FT(1) - a() )
                         : PointS2<FT>( FT(1) + b(), -(a()+c())/b() - a() );
  // we add i times the direction
  if (is_vertical())
    return PointS2<FT>( (-b()-c())/a() + FT(i)*b(), FT(1) - FT(i)*a() );
  return PointS2<FT>( FT(1) + FT(i)*b(), -(a()+c())/b() - FT(i)*a() );
}

template < class FT >
CGAL_KERNEL_INLINE
PointS2<FT> LineS2<FT>::point() const
{
  return is_vertical() ? PointS2<FT>( (-b()-c())/a(), FT(1) )
                       : PointS2<FT>( FT(1), -(a()+c())/b());
}

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
PointS2<FT> LineS2<FT>::projection(const PointS2<FT>& p) const
{
  if (is_horizontal())
      return PointS2<FT>(p.x(), -c()/b());

  if (is_vertical())
      return PointS2<FT>( -c()/a(), p.y());

  FT ab = a()/b(), ba = b()/a(), ca = c()/a();
  FT y = ( -p.x() + ab*p.y() - ca ) / ( ba + ab );
  return PointS2<FT>(-ba * y - ca, y);
}

template < class FT >
inline
DirectionS2<FT> LineS2<FT>::direction() const
{ return DirectionS2<FT>( b(), -a() ); }

template < class FT >
CGAL_KERNEL_INLINE
Oriented_side LineS2<FT>::oriented_side(const PointS2<FT>& p) const
{ return Oriented_side(CGAL_NTS sign(a()*p.x() + b()*p.y() + c())); }

template < class FT >
inline
bool LineS2<FT>::has_on_boundary(const PointS2<FT>& p) const
{ return (a()*p.x() + b()*p.y() + c()) == FT(0); }

template < class FT >
inline
bool LineS2<FT>::has_on_positive_side(const PointS2<FT>& p) const
{ return (a()*p.x() + b()*p.y() + c()) >  FT(0); }

template < class FT >
CGAL_KERNEL_INLINE
bool LineS2<FT>::has_on_negative_side(const PointS2<FT>& p) const
{ return (a()*p.x() + b()*p.y() + c()) <  FT(0); }

template < class FT >
inline
bool LineS2<FT>::is_horizontal() const
{ return a() == FT(0) ; }

template < class FT >
inline
bool LineS2<FT>::is_vertical() const
{ return b() == FT(0) ; }

template < class FT >
inline
bool LineS2<FT>::is_degenerate() const
{ return (a() == FT(0)) && (b() == FT(0)) ; }

template < class FT >
inline
LineS2<FT> LineS2<FT>::transform(const Aff_transformationS2<FT>& t) const
{ return LineS2<FT>( t.transform(point(0) ), t.transform(direction() )); }



#ifndef CGAL_NO_OSTREAM_INSERT_LINES2
template < class FT >
std::ostream& operator<<(std::ostream &os, const LineS2<FT> &l)
{

    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << l.a() << ' ' << l.b() << ' ' << l.c();
    case IO::BINARY :
        write(os, l.a());
        write(os, l.b());
        write(os, l.c());
        return os;
    default:
        return os << "LineS2(" << l.a() << ", " << l.b() << ", " << l.c() <<')';
    }
}
#endif // CGAL_NO_OSTREAM_INSERT_LINES2

#ifndef CGAL_NO_ISTREAM_EXTRACT_LINES2
template < class FT >
std::istream& operator>>(std::istream &is, LineS2<FT> &p)
{
    FT a, b, c;
    switch(is.iword(IO::mode)) {
    case IO::ASCII :
        is >> a >> b >> c;
        break;
    case IO::BINARY :
        read(is, a);
        read(is, b);
        read(is, c);
        break;
    default:
        std::cerr << "" << std::endl;
        std::cerr << "Stream must be in ascii or binary mode" << std::endl;
        break;
    }
    p = LineS2<FT>(a, b, c);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_LINES2



CGAL_END_NAMESPACE

#endif // CGAL_LINES2_H
