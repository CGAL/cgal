// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
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
// file          : LineH2.h
// package       : H2
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_LINEH2_H
#define CGAL_LINEH2_H

#include <CGAL/PointH2.h>

CGAL_BEGIN_NAMESPACE

template < class R_ >
class LineH2
  : public R_::Line_handle_2
{
public:
    typedef R_                                    R;
    typedef typename R::FT                        FT;
    typedef typename R::RT                        RT;

    typedef typename R::Line_handle_2             Line_handle_2_;
    typedef typename Line_handle_2_::element_type Line_ref_2;

    LineH2();
    LineH2(const PointH2<R>& p, const PointH2<R>& q);
    LineH2(const RT& a, const RT& b, const RT& c);
    LineH2(const SegmentH2<R>& s);
    LineH2(const RayH2<R>& r);
    LineH2(const PointH2<R>& p, const DirectionH2<R>& d);

    bool           operator==(const LineH2<R>& l) const ;
    bool           operator!=(const LineH2<R>& l) const ;

    RT             a() const;
    RT             b() const;
    RT             c() const;
    const RT&      a_ref() const { return Ptr()->e0; }
    const RT&      b_ref() const { return Ptr()->e1; }
    const RT&      c_ref() const { return Ptr()->e2; }

    FT             x_at_y(FT y) const;
    FT             y_at_x(FT x) const;

    LineH2<R>  perpendicular(const PointH2<R>& p ) const;
    LineH2<R>  opposite() const;
    PointH2<R> point() const;
    PointH2<R> point(int i) const;
    PointH2<R> projection(const PointH2<R>& p) const;
    DirectionH2<R>
                   direction() const;
    Oriented_side  oriented_side( const PointH2<R>& p ) const;
    bool           has_on( const PointH2<R>& p ) const;
    bool           has_on_boundary( const PointH2<R>& p ) const;
    bool           has_on_positive_side( const PointH2<R>& p ) const;
    bool           has_on_negative_side( const PointH2<R>& p ) const;
    bool           is_horizontal() const;
    bool           is_vertical()   const;
    bool           is_degenerate() const;

    LineH2<R>  transform(const Aff_transformationH2<R>&) const;
};

#ifdef CGAL_CFG_TYPENAME_BUG
#define typename
#endif

template < class R >
CGAL_KERNEL_CTOR_INLINE
LineH2<R>::LineH2()
 : Line_handle_2_ ( Line_ref_2() )
{}

template < class R >
CGAL_KERNEL_CTOR_MEDIUM_INLINE
LineH2<R>::LineH2(const PointH2<R>& p, const PointH2<R>& q)
 : Line_handle_2_ ( Line_ref_2(
  //  a() * X + b() * Y + c() * W() == 0
  //      |    X        Y       W     |
  //      |  p.hx()   p.hy()  p.hw()  |
  //      |  q.hx()   q.hy()  q.hw()  |

            p.hy()*q.hw() - p.hw()*q.hy(),
            p.hw()*q.hx() - p.hx()*q.hw(),
            p.hx()*q.hy() - p.hy()*q.hx() ))
{}

template < class R >
CGAL_KERNEL_CTOR_INLINE
LineH2<R>::LineH2(const RT& a, const RT& b, const RT& c)
 : Line_handle_2_ ( Line_ref_2(a,b,c) )
{}

template < class R >
CGAL_KERNEL_CTOR_INLINE
LineH2<R>::LineH2(const SegmentH2<R>& s)
{
  PointH2<R> p = s.start();
  PointH2<R> q = s.end();
  initialize_with( Line_ref_2 (
            p.hy()*q.hw() - p.hw()*q.hy(),
            p.hw()*q.hx() - p.hx()*q.hw(),
            p.hx()*q.hy() - p.hy()*q.hx() ) );
}

template < class R >
CGAL_KERNEL_CTOR_INLINE
LineH2<R>::LineH2(const RayH2<R>& r)
{
  PointH2<R> p = r.start();
  PointH2<R> q = r.second_point();
  initialize_with( Line_ref_2 (
            p.hy()*q.hw() - p.hw()*q.hy(),
            p.hw()*q.hx() - p.hx()*q.hw(),
            p.hx()*q.hy() - p.hy()*q.hx() ) );
}

template < class R >
CGAL_KERNEL_CTOR_INLINE
LineH2<R>::LineH2(const PointH2<R>& p,
                      const DirectionH2<R>& d)
{
  PointH2<R> q = p + VectorH2<R>(d);
  initialize_with( Line_ref_2 (
            p.hy()*q.hw() - p.hw()*q.hy(),
            p.hw()*q.hx() - p.hx()*q.hw(),
            p.hx()*q.hy() - p.hy()*q.hx() ) );
}

template < class R >
inline
typename LineH2<R>::RT
LineH2<R>::a() const
{ return Ptr()->e0; }

template < class R >
inline
typename LineH2<R>::RT
LineH2<R>::b() const
{ return Ptr()->e1; }

template < class R >
inline
typename LineH2<R>::RT
LineH2<R>::c() const
{ return Ptr()->e2; }


template < class R >
CGAL_KERNEL_INLINE
typename LineH2<R>::FT
LineH2<R>::x_at_y(FT y) const
{
  CGAL_kernel_precondition( !is_degenerate() );
  return (FT(-b())*y - FT(c()) )/FT(a());
}

template < class R >
CGAL_KERNEL_INLINE
typename LineH2<R>::FT
LineH2<R>::y_at_x(FT x) const
{
  CGAL_kernel_precondition( !is_degenerate() );
  return (FT(-a())*x - FT(c()) )/FT(b());
}

template < class R >
CGAL_KERNEL_INLINE
LineH2<R>
LineH2<R>::perpendicular(const PointH2<R>& p ) const
{
  CGAL_kernel_precondition( !is_degenerate() );
  return LineH2<R>( -b()*p.hw(), a()*p.hw(), b()*p.hx() - a()*p.hy() );
}

template < class R >
inline
LineH2<R>
LineH2<R>::opposite() const
{ return LineH2<R>( -a(), -b(), -c() ); }

template < class R >
CGAL_KERNEL_INLINE
PointH2<R>
LineH2<R>::point() const
{
  CGAL_kernel_precondition( !is_degenerate() );
  if (is_vertical() )
  {
      return PointH2<R>(-c(), RT(0)  , a() );
  }
  else
  {
      return PointH2<R>(RT(0)  , -c(), b() );
  }
}

template < class R >
CGAL_KERNEL_INLINE
PointH2<R>
LineH2<R>::point(int i) const
{ return point() + RT(i) * (direction().to_vector()); }

template < class R >
CGAL_KERNEL_INLINE
PointH2<R>
LineH2<R>::projection(const PointH2<R>& p) const
{
  CGAL_kernel_precondition( !is_degenerate() );
  LineH2<R>  l( p, DirectionH2<R>( a(), b() ));
  return PointH2<R>( b()*l.c() - l.b()*c(),
                              l.a()*c() - a()*l.c(),
                              a()*l.b() - l.a()*b() );
}

template < class R >
CGAL_KERNEL_INLINE
DirectionH2<R>
LineH2<R>::direction() const
{
  CGAL_kernel_precondition( !is_degenerate() );
  return DirectionH2<R>( b(), -a() );
}

template < class R >
CGAL_KERNEL_INLINE
LineH2<R>
LineH2<R>::transform(const Aff_transformationH2<R>& t) const
{
  CGAL_kernel_precondition( !is_degenerate() );
  PointH2<R> p = point() + direction().to_vector();
  return LineH2<R>( t.transform(point() ), t.transform(p) );
}

#ifndef CGAL_NO_OSTREAM_INSERT_LINEH2
template < class R >
std::ostream &
operator<<(std::ostream &os, const LineH2<R> &l)
{
  switch(os.iword(IO::mode))
  {
    case IO::ASCII :
        return os << l.a() << ' ' << l.b() << ' ' << l.c();
    case IO::BINARY :
        write(os, l.a());
        write(os, l.b());
        write(os, l.c());
        return os;
    default:
       return os << "LineH2(" << l.a() << ", " << l.b() << ", " << l.c() <<')';
  }
}
#endif // CGAL_NO_OSTREAM_INSERT_LINEH2

#ifndef CGAL_NO_ISTREAM_EXTRACT_LINEH2
template < class R >
std::istream &
operator>>(std::istream &is, LineH2<R> &p)
{
  typename R::RT a, b, c;
  switch(is.iword(IO::mode))
  {
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
  p = LineH2<R>(a, b, c);
  return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_LINEH2

template < class R >
CGAL_KERNEL_INLINE
bool
LineH2<R>::has_on( const PointH2<R>& p ) const
{
  CGAL_kernel_precondition( !is_degenerate() );
  return ( ( a()*p.hx() + b()*p.hy() + c()*p.hw() ) == RT(0)   );
}

template < class R >
CGAL_KERNEL_INLINE
bool
LineH2<R>::has_on_boundary( const PointH2<R>& p ) const
{
  CGAL_kernel_precondition( !is_degenerate() );
  return ( ( a()*p.hx() + b()*p.hy() + c()*p.hw() ) == RT(0)   );
}

template < class R >
CGAL_KERNEL_INLINE
bool
LineH2<R>::has_on_positive_side( const PointH2<R>& p ) const
{
  CGAL_kernel_precondition( !is_degenerate() );
  return ( ( a()*p.hx() + b()*p.hy() + c()*p.hw() ) > RT(0)   );
}

template < class R >
CGAL_KERNEL_INLINE
bool
LineH2<R>::has_on_negative_side( const PointH2<R>& p ) const
{
  CGAL_kernel_precondition( !is_degenerate() );
  return ( ( a()*p.hx() + b()*p.hy() + c()*p.hw() ) < RT(0)   );
}

template < class R >
CGAL_KERNEL_INLINE
Oriented_side
LineH2<R>::oriented_side( const PointH2<R>& p ) const
{
  CGAL_kernel_precondition( !is_degenerate() );
  RT v = a()*p.hx() + b()*p.hy() + c()*p.hw();
  if (v > RT(0)   )
  {
      return ON_POSITIVE_SIDE;
  }
  else
  {
      return (v < RT(0)   ) ? ON_NEGATIVE_SIDE : ON_ORIENTED_BOUNDARY;
  }
}

template < class R >
inline
bool
LineH2<R>::is_horizontal() const
{ return ( a() == RT(0)   ); }

template < class R >
inline
bool
LineH2<R>::is_vertical() const
{ return ( b() == RT(0)   ); }

template < class R >
inline
bool
LineH2<R>::is_degenerate() const
{ return (a() == RT(0)  )&&(b() == RT(0)  ) ; }

template < class R >
CGAL_KERNEL_MEDIUM_INLINE
bool
LineH2<R>::operator==(const LineH2<R>& l) const
{
  if (  (a() * l.c() != l.a() * c() )
      ||(b() * l.c() != l.b() * c() ) )
  {
      return false;
  }
  int sc  = static_cast<int>(CGAL_NTS sign(c()));
  int slc = static_cast<int>(CGAL_NTS sign(l.c()));
  if ( sc == slc )
  {
      if (sc == 0)
      {
          return (  (a()*l.b() == b()*l.a() )
                  &&(CGAL_NTS sign(a() )== CGAL_NTS sign( l.a() ))
                  &&(CGAL_NTS sign(b() )== CGAL_NTS sign( l.b() )) );
      }
      else
      {
          return true;
      }
  }
  else
  {
      return false;
  }
}

template < class R >
inline
bool
LineH2<R>::operator!=(const LineH2<R>& l) const
{ return !(*this == l); }

#ifdef CGAL_CFG_TYPENAME_BUG
#undef typename
#endif

CGAL_END_NAMESPACE

#endif // CGAL_LINEH2_H
