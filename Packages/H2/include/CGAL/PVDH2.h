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
// file          : PVDH2.h
// package       : H2
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_PVDH2_H
#define CGAL_PVDH2_H

#define CGAL_POINTH2_H
#define CGAL_VECTORH2_H
#define CGAL_DIRECTIONH2_H

#include <CGAL/homogeneous_classes.h>
#include <CGAL/Origin.h>
#include <CGAL/Bbox_2.h>

CGAL_BEGIN_NAMESPACE

template < class R_ >
class PointH2
  : public R_::Point_handle_2
{
public:
  typedef R_                                    R;
  typedef typename R::FT                        FT;
  typedef typename R::RT                        RT;

  typedef typename R::Point_handle_2            Point_handle_2_;
  typedef typename Point_handle_2_::element_type        Point_ref_2;

            PointH2();
            PointH2(const Origin & o);
            PointH2(const PointH2<R> & p);
            PointH2(const VectorH2<R>& v);
            PointH2(const RT& hx, const RT& hy );
            PointH2(const RT& hx, const RT& hy, const RT& hw );

    bool    operator==( const PointH2<R>& p) const;
    bool    operator!=( const PointH2<R>& p) const;

    RT      hx() const { return Ptr()->e0; };
    RT      hy() const { return Ptr()->e1; };
    RT      hw() const { return Ptr()->e2; };

    FT      x()  const { return FT(hx()) / FT(hw()); };
    FT      y()  const { return FT(hy()) / FT(hw()); };

    FT      cartesian(int i)   const;
    FT      operator[](int i)  const;
    RT      homogeneous(int i) const;

    // and for efficiency in the predicates:
    const RT&     hx_ref() const { return Ptr()->e0; };
    const RT&     hy_ref() const { return Ptr()->e1; };
    const RT&     hw_ref() const { return Ptr()->e2; };

    int     dimension() const;
    Bbox_2  bbox() const;

    PointH2<R> transform( const Aff_transformationH2<R> & t) const;
    DirectionH2<R> direction() const;
};

template < class R_ >
class VectorH2
  : public R_::Vector_handle_2
{
public:
  typedef R_                                    R;
  typedef typename R::FT                        FT;
  typedef typename R::RT                        RT;

  typedef typename R::Vector_handle_2           Vector_handle_2_;
  typedef typename Vector_handle_2_::element_type  Vector_ref_2;

            VectorH2();
            VectorH2(const VectorH2<R>& v);
            VectorH2(const PointH2<R>& a, const PointH2<R>& b);
            VectorH2(const Null_vector &);
            VectorH2(const RT& x, const RT& y);
            VectorH2(const RT& x, const RT& y, const RT& w );

    bool    operator==( const VectorH2<R>& v) const;
    bool    operator!=( const VectorH2<R>& v) const;
    bool    operator==( const Null_vector&) const;
    bool    operator!=( const Null_vector& v) const;

    RT      hx() const { return Ptr()->e0; };
    RT      hy() const { return Ptr()->e1; };
    RT      hw() const { return Ptr()->e2; };

    FT      x()  const { return FT(hx()) / FT(hw()); };
    FT      y()  const { return FT(hy()) / FT(hw()); };

    FT      cartesian(int i)   const;
    RT      homogeneous(int i) const;
    FT      operator[](int i)  const;

    int     dimension() const;
    DirectionH2<R> direction() const;
    VectorH2<R> transform(const Aff_transformationH2<R>& t ) const;
    VectorH2<R> perpendicular(const Orientation& o ) const;

    FT      operator*( const VectorH2<R>& v) const;
    VectorH2<R> operator-() const;
    VectorH2<R> opposite() const;

// undocumented:
            VectorH2(const DirectionH2<R> & dir);
protected:
            VectorH2(const PointH2<R> & p);
};

template < class R_ >
class DirectionH2
  : public R_::Direction_handle_2
{
public:
  typedef R_                                    R;
  typedef typename R::FT                        FT;
  typedef typename R::RT                        RT;

  typedef typename R::Direction_handle_2        Direction_handle_2_;
  typedef typename Direction_handle_2_::element_type  Direction_ref_2;

            DirectionH2();
            DirectionH2(const DirectionH2<R>& d );
            DirectionH2(const PointH2<R> & p );
            DirectionH2(const VectorH2<R> & v );
            DirectionH2(const LineH2<R> & l );
            DirectionH2(const RayH2<R> & r );
            DirectionH2(const SegmentH2<R> & s );
            DirectionH2(const RT& x, const RT& y);
            DirectionH2(const RT& x, const RT& y, const RT& w );

    bool    operator==( const DirectionH2<R>& d) const;
    bool    operator!=( const DirectionH2<R>& d) const;
    bool    operator< ( const DirectionH2<R>& d) const;
    bool    operator<=( const DirectionH2<R>& d) const;
    bool    operator> ( const DirectionH2<R>& d) const;
    bool    operator>=( const DirectionH2<R>& d) const;
    bool    counterclockwise_in_between( const DirectionH2<R>& d1,
                                         const DirectionH2<R>& d2 ) const;

    DirectionH2<R> operator-() const;

    VectorH2<R>    to_vector() const;

    RT      x() const { return Ptr()->e0; };
    RT      y() const { return Ptr()->e1; };

    RT      delta(int i) const;
    RT      dx() const { return Ptr()->e0; };
    RT      dy() const { return Ptr()->e1; };

    DirectionH2<R> perpendicular(const Orientation &o) const;
    DirectionH2<R> transform(const Aff_transformationH2<R> &) const;
};

#ifdef CGAL_CFG_TYPENAME_BUG
#define typename
#endif

template < class R >
CGAL_KERNEL_CTOR_INLINE
PointH2<R>::PointH2()
 : Point_handle_2_ ( Point_ref_2()) {}

template < class R >
CGAL_KERNEL_CTOR_INLINE
PointH2<R>::PointH2(const Origin&)
 : Point_handle_2_ ( Point_ref_2( RT(0), RT(0), RT(1))) {}

template < class R >
CGAL_KERNEL_CTOR_INLINE
PointH2<R>::PointH2(const PointH2<R>& p)
  : Point_handle_2_ (p)
{}

template < class R >
CGAL_KERNEL_CTOR_INLINE
PointH2<R>::PointH2(const RT& hx, const RT& hy)
 : Point_handle_2_ ( Point_ref_2( hx, hy, RT(1) )) {}

template < class R >
CGAL_KERNEL_CTOR_INLINE
PointH2<R>::PointH2(const RT& hx, const RT& hy, const RT& hw)
{
   if ( hw >= RT(0)   )
   { initialize_with( Point_ref_2( hx, hy, hw)); }
   else
   { initialize_with( Point_ref_2(-hx,-hy,-hw)); }
}

template < class R >
CGAL_KERNEL_CTOR_INLINE
PointH2<R>::PointH2(const VectorH2<R>& v)
  : Point_handle_2_ (v)
{}

template < class R >
CGAL_KERNEL_INLINE
bool
PointH2<R>::operator==( const PointH2<R>& p) const
{
  return (  (hx() * p.hw() == p.hx() * hw() )
          &&(hy() * p.hw() == p.hy() * hw() ) );
}

template < class R >
inline
bool
PointH2<R>::operator!=( const PointH2<R>& p) const
{ return !(*this == p); }   /* XXX */

template < class R >
CGAL_KERNEL_INLINE
typename PointH2<R>::FT
PointH2<R>::cartesian(int i) const
{
  CGAL_kernel_precondition( (i==0 || i==1) );
  if (i==0)
      return x();
  return y();
}

template < class R >
CGAL_KERNEL_INLINE
typename PointH2<R>::RT
PointH2<R>::homogeneous(int i) const
{
  CGAL_kernel_precondition( (i>=0) && (i<=2) );
  if (i==0)
      return hx();
  if (i==1)
      return hy();
  return hw();
}

template < class R >
inline
typename PointH2<R>::FT
PointH2<R>::operator[](int i) const
{ return cartesian(i); }


template < class R >
inline
int
PointH2<R>::dimension() const
{ return 2; }

template < class R >
CGAL_KERNEL_INLINE
DirectionH2<R>
PointH2<R>::direction() const
{ return DirectionH2<R>(*this); }


template < class R >
CGAL_KERNEL_CTOR_INLINE
VectorH2<R>::VectorH2()
 : Vector_handle_2_ ( Vector_ref_2()) {}

template < class R >
CGAL_KERNEL_CTOR_INLINE
VectorH2<R>::VectorH2(const Null_vector &)
 : Vector_handle_2_ ( Vector_ref_2(RT(0), RT(0), RT(1) )) {}

template < class R >
CGAL_KERNEL_CTOR_INLINE
VectorH2<R>::VectorH2(const VectorH2<R>& v)
  : Vector_handle_2_ (v) {}

template < class R >
CGAL_KERNEL_CTOR_INLINE
VectorH2<R>::VectorH2(const PointH2<R>& a, const PointH2<R>& b)
  : Vector_handle_2_ (b-a) {}

template < class R >
CGAL_KERNEL_CTOR_INLINE
VectorH2<R>::VectorH2(const RT& x, const RT& y)
 : Vector_handle_2_ ( Vector_ref_2( x,  y,  RT(1) )) {}

template < class R >
CGAL_KERNEL_CTOR_INLINE
VectorH2<R>::VectorH2(const RT& x, const RT& y, const RT& w)
{
  if ( w >= RT(0)   )
  { initialize_with( Vector_ref_2( x,  y,  w)); }
  else
  { initialize_with( Vector_ref_2(-x, -y, -w)); }
}

template < class R >
CGAL_KERNEL_CTOR_INLINE
VectorH2<R>::VectorH2(const PointH2<R> & p)
  : Vector_handle_2_ ( p) {}

template < class R >
CGAL_KERNEL_CTOR_INLINE
VectorH2<R>::VectorH2(const DirectionH2<R> & dir)
  : Vector_handle_2_ ( dir) {}

template < class R >
inline
bool
VectorH2<R>::operator==( const Null_vector&) const
{ return (hx() == RT(0)) && (hy() == RT(0)); }

template < class R >
inline
bool
VectorH2<R>::operator!=( const Null_vector& v) const
{ return !(*this == v); }

template < class R >
CGAL_KERNEL_INLINE
bool
VectorH2<R>::operator==( const VectorH2<R>& v) const
{
  return (  (hx() * v.hw() == v.hx() * hw() )
          &&(hy() * v.hw() == v.hy() * hw() ) );
}

template < class R >
inline
bool
VectorH2<R>::operator!=( const VectorH2<R>& v) const
{ return !(*this == v); }  /* XXX */

template < class R >
CGAL_KERNEL_INLINE
typename VectorH2<R>::FT
VectorH2<R>::cartesian(int i) const
{
  CGAL_kernel_precondition( (i==0 || i==1) );
  if (i==0)
      return x();
  return y();
}

template < class R >
CGAL_KERNEL_INLINE
typename VectorH2<R>::RT
VectorH2<R>::homogeneous(int i) const
{
  CGAL_kernel_precondition( (i>=0) && (i<=2) );
  if (i==0)
      return hx();
  if (i==1)
      return hy();
  return hw();
}

template < class R >
inline
typename VectorH2<R>::FT
VectorH2<R>::operator[](int i) const
{ return cartesian(i); }

template < class R >
inline
int
VectorH2<R>::dimension() const
{ return 2; }

template < class R >
CGAL_KERNEL_INLINE
DirectionH2<R>
VectorH2<R>::direction() const
{ return DirectionH2<R>(*this); }
template < class R >
inline
VectorH2<R>
VectorH2<R>::operator-() const
{ return VectorH2<R>(- hx(), - hy(), hw() ); }

template < class R >
inline
VectorH2<R>
VectorH2<R>::opposite() const
{ return VectorH2<R>(- hx(), - hy(), hw() ); }

template <class R >
CGAL_KERNEL_CTOR_INLINE
DirectionH2<R>::DirectionH2()
 : Direction_handle_2_ ( Direction_ref_2()) {}

template <class R >
CGAL_KERNEL_CTOR_INLINE
DirectionH2<R>::DirectionH2(const DirectionH2<R>& d )
  : Direction_handle_2_ ( d ) {}

template <class R >
CGAL_KERNEL_CTOR_INLINE
DirectionH2<R>::DirectionH2(const PointH2<R> & p )
  : Direction_handle_2_ ( p) {}

template <class R >
CGAL_KERNEL_CTOR_INLINE
DirectionH2<R>::DirectionH2(const VectorH2<R> & v )
  : Direction_handle_2_ ( v) {}

template <class R >
CGAL_KERNEL_CTOR_INLINE
DirectionH2<R>::DirectionH2(const LineH2<R> & l )
  : Direction_handle_2_ ( l.direction()) {}

template <class R >
CGAL_KERNEL_CTOR_INLINE
DirectionH2<R>::DirectionH2(const RayH2<R> & r )
  : Direction_handle_2_ ( r.direction()) {}

template <class R >
CGAL_KERNEL_CTOR_INLINE
DirectionH2<R>::DirectionH2(const SegmentH2<R> & s )
  : Direction_handle_2_ ( s.direction()) {}

template <class R >
CGAL_KERNEL_CTOR_INLINE
DirectionH2<R>::DirectionH2(const RT& x, const RT& y)
 : Direction_handle_2_ ( Direction_ref_2( x, y, RT(1) )) {}

template <class R >
CGAL_KERNEL_CTOR_INLINE
DirectionH2<R>::DirectionH2(const RT& x, const RT& y, const RT& w )
{
  if (w > RT(0)   )
  { initialize_with( Direction_ref_2( x, y, w)); }
  else
  { initialize_with( Direction_ref_2(-x,-y,-w)); }
}


template <class R >
CGAL_KERNEL_INLINE
bool
DirectionH2<R>::operator==( const DirectionH2<R>& d) const
{
  return (  ( x() * d.y() == y() * d.x() )
          &&( CGAL_NTS sign( x() ) == CGAL_NTS sign( d.x() ) )
          &&( CGAL_NTS sign( y() ) == CGAL_NTS sign( d.y() ) ) );
}

template <class R >
inline
bool
DirectionH2<R>::operator!=( const DirectionH2<R>& d) const
{ return !(*this == d); }

template <class R >
inline
DirectionH2<R>
DirectionH2<R>::operator-() const
{ return DirectionH2<R>( - x(), - y() ); }

template <class R >
CGAL_KERNEL_INLINE
typename DirectionH2<R>::RT
DirectionH2<R>::delta(int i) const
{
  CGAL_kernel_precondition( ( i == 0 ) || ( i == 1 ) );
  if (i == 0)
      return dx();
  return dy();
}


template <class R>
CGAL_KERNEL_INLINE
VectorH2<R>
operator+(const VectorH2<R>& u, const VectorH2<R>& v)
{
  return VectorH2<R>( u.hx()*v.hw() + v.hx()*u.hw(),
                          u.hy()*v.hw() + v.hy()*u.hw(),
                          u.hw()*v.hw() );
}

template <class R>
CGAL_KERNEL_INLINE
VectorH2<R>
operator-(const VectorH2<R>& u, const VectorH2<R>& v)
{
  return VectorH2<R>( u.hx()*v.hw() - v.hx()*u.hw(),
                          u.hy()*v.hw() - v.hy()*u.hw(),
                          u.hw()*v.hw() );
}

template <class R>
CGAL_KERNEL_INLINE
typename VectorH2<R>::FT
VectorH2<R>::operator*(const VectorH2<R>& v) const
{
  typedef typename R::RT RT;
  typedef typename R::FT FT;
  return FT( RT(hx()*v.hx() + hy()*v.hy()) ) / FT( RT(hw()*v.hw() ) );
}

template <class R>
CGAL_KERNEL_INLINE
VectorH2<R>
operator/(const VectorH2<R>& v, const typename R::RT& f)
{ return VectorH2<R>( v.hx(), v.hy(), v.hw()*f ); }

template <class R>
CGAL_KERNEL_INLINE
VectorH2<R>
operator*(const VectorH2<R>& v, const typename R::RT& f)
{ return VectorH2<R>( v.hx()*f, v.hy()*f, v.hw() ); }

template <class R>
CGAL_KERNEL_INLINE
VectorH2<R>
operator*(const typename R::RT& f, const VectorH2<R>& v)
{ return VectorH2<R>( v.hx()*f, v.hy()*f, v.hw() ); }

template <class R>
inline
PointH2<R>
origin_plus_vector(const VectorH2<R>& v)
{ return PointH2<R>( v ); }

template <class R>
inline
PointH2<R>
operator+(const Origin&, const VectorH2<R>& v)
{ return origin_plus_vector( v ); }

template <class R>
inline
PointH2<R>
origin_minus_vector(const VectorH2<R>& v)
{ return PointH2<R>( v.opposite() ); }

template <class R>
inline
PointH2<R>
operator-(const Origin&, const VectorH2<R>& v)
{ return origin_minus_vector( v ); }

template <class R>
inline
VectorH2<R>
point_minus_origin(const PointH2<R>& p)
{ return VectorH2<R>( p ); }

template <class R>
inline
VectorH2<R>
operator-(const PointH2<R>& p, const Origin&)
{ return point_minus_origin( p ); }

template <class R>
inline
VectorH2<R>
origin_minus_point(const PointH2<R>& p)
{ return  VectorH2<R>( p ).opposite(); }

template <class R>
inline
VectorH2<R>
operator-(const Origin&, const PointH2<R>& p)
{ return  origin_minus_point( p ); }


template <class R>
CGAL_KERNEL_INLINE
PointH2<R>
operator+(const PointH2<R>& p, const VectorH2<R>& v)
{
  return PointH2<R>( p.hx()*v.hw() + v.hx()*p.hw(),
                         p.hy()*v.hw() + v.hy()*p.hw(),
                         p.hw()*v.hw() );
}

template <class R>
CGAL_KERNEL_INLINE
PointH2<R>
operator-(const PointH2<R>& p, const VectorH2<R>& v)
{
  return PointH2<R>( p.hx()*v.hw() - v.hx()*p.hw(),
                         p.hy()*v.hw() - v.hy()*p.hw(),
                         p.hw()*v.hw() );
}

template <class R>
CGAL_KERNEL_INLINE
VectorH2<R>
operator-(const PointH2<R>& p, const PointH2<R>& q)
{
  return VectorH2<R>( p.hx()*q.hw() - q.hx()*p.hw(),
                          p.hy()*q.hw() - q.hy()*p.hw(),
                          p.hw()*q.hw() );
}

#ifdef CGAL_CFG_TYPENAME_BUG
#undef typename
#endif

CGAL_END_NAMESPACE

#include <CGAL/predicates_on_directionsH2.h>

CGAL_BEGIN_NAMESPACE

template <class R >
inline
bool
DirectionH2<R>::operator< (const DirectionH2<R>& d) const
{ return (compare_angle_with_x_axis(*this,d) == SMALLER); }

template <class R >
inline
bool
DirectionH2<R>::operator> (const DirectionH2<R>& d) const
{ return (compare_angle_with_x_axis(*this,d) == LARGER); }

template <class R >
inline
bool
DirectionH2<R>::operator>= (const DirectionH2<R>& d) const
{ return !(compare_angle_with_x_axis(*this,d) == SMALLER); }

template <class R >
inline
bool
DirectionH2<R>::operator<= (const DirectionH2<R>& d) const
{ return !(compare_angle_with_x_axis(*this,d) == LARGER); }

template <class R >
CGAL_KERNEL_INLINE
bool
DirectionH2<R>::
counterclockwise_in_between( const DirectionH2<R>& d1,
                             const DirectionH2<R>& d2) const
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

CGAL_END_NAMESPACE

#include <CGAL/Aff_transformationH2.h>

CGAL_BEGIN_NAMESPACE

template < class R >
CGAL_KERNEL_MEDIUM_INLINE
Bbox_2
PointH2<R>::bbox() const
{
#ifndef CGAL_CFG_NO_NAMESPACE
  using std::swap;
#endif // CGAL_CFG_NO_NAMESPACE

  // double eps  = exp2(-52);
  // the following is faster as it can be evaluated at compile time
  // and it is machine independent
  double eps  = double(1.0) /(double(1<<26) * double(1<<26));
  double hxd  = CGAL::to_double( hx() );
  double hyd  = CGAL::to_double( hy() );
  double hwd  = CGAL::to_double( hw() );
  double xmin = ( hxd - eps*hxd ) / ( hwd + eps*hwd );
  double xmax = ( hxd + eps*hxd ) / ( hwd - eps*hwd );
  double ymin = ( hyd - eps*hyd ) / ( hwd + eps*hwd );
  double ymax = ( hyd + eps*hyd ) / ( hwd - eps*hwd );
  if ( hx() < RT(0)   ) { swap(xmin, xmax); }
  if ( hy() < RT(0)   ) { swap(ymin, ymax); }
  return Bbox_2(xmin, ymin, xmax, ymax);
}

template < class R >
inline
PointH2<R>
PointH2<R>::transform(const Aff_transformationH2<R>& t) const
{ return t.transform(*this); }

#ifndef CGAL_NO_OSTREAM_INSERT_POINTH2
template < class R >
std::ostream &
operator<<(std::ostream &os, const PointH2<R> &p)
{
  switch(os.iword(IO::mode))
  {
    case IO::ASCII :
        return os << p.hx() << ' ' << p.hy() << ' ' << p.hw();
    case IO::BINARY :
        write(os, p.hx());
        write(os, p.hy());
        write(os, p.hw());
        return os;
    default:
        return os << "PointH2(" << p.hx() << ", "
                                << p.hy() << ", "
                                << p.hw() << ')';
  }
}
#endif // CGAL_NO_OSTREAM_INSERT_POINTH2

#ifndef CGAL_NO_ISTREAM_EXTRACT_POINTH2
template < class R >
std::istream &
operator>>(std::istream &is, PointH2<R> &p)
{
  typename R::RT hx, hy, hw;
  switch(is.iword(IO::mode))
  {
    case IO::ASCII :
        is >> hx >> hy >> hw;
        break;
    case IO::BINARY :
        read(is, hx);
        read(is, hy);
        read(is, hw);
        break;
    default:
        std::cerr << "" << std::endl;
        std::cerr << "Stream must be in ascii or binary mode" << std::endl;
        break;
  }
  p = PointH2<R>(hx, hy, hw);
  return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_POINTH2

template < class R >
CGAL_KERNEL_INLINE
VectorH2<R>
VectorH2<R>::perpendicular(const Orientation& o) const
{
  CGAL_kernel_precondition(o != COLLINEAR);
  if (o == COUNTERCLOCKWISE)
  {
      return VectorH2<R>(-hy(), hx(), hw());
  }
  else
  {
      return VectorH2<R>(hy(), -hx(), hw());
  }
}

template < class R >
inline
VectorH2<R>
VectorH2<R>::transform(const Aff_transformationH2<R>& t) const
{ return t.transform(*this); }


#ifndef CGAL_NO_OSTREAM_INSERT_VECTORH2
template < class R >
std::ostream &
operator<<(std::ostream &os, const VectorH2<R> &p)
{
  switch(os.iword(IO::mode))
  {
    case IO::ASCII :
        return os << p.hx() << ' ' << p.hy() << ' ' << p.hw();
    case IO::BINARY :
        write(os, p.hx());
        write(os, p.hy());
        write(os, p.hw());
        return os;
    default:
        return os << "VectorH2(" << p.hx() << ", "
                                 << p.hy() << ", "
                                 << p.hw() << ')';
  }
}
#endif // CGAL_NO_OSTREAM_INSERT_VECTORH2

#ifndef CGAL_NO_ISTREAM_EXTRACT_VECTORH2
template < class R >
std::istream &
operator>>(std::istream &is, VectorH2<R> &p)
{
  typename R::RT hx, hy, hw;
  switch(is.iword(IO::mode))
  {
    case IO::ASCII :
        is >> hx >> hy >> hw;
        break;
    case IO::BINARY :
        read(is, hx);
        read(is, hy);
        read(is, hw);
        break;
    default:
        std::cerr << "" << std::endl;
        std::cerr << "Stream must be in ascii or binary mode" << std::endl;
        break;
  }
  p = VectorH2<R>(hx, hy, hw);
  return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_VECTORH2

template <class R >
CGAL_KERNEL_INLINE
DirectionH2<R>
DirectionH2<R>::perpendicular(const Orientation& o) const
{
  CGAL_kernel_precondition(o != COLLINEAR);
  if (o == COUNTERCLOCKWISE)
  {
      return DirectionH2<R>(-dy(), dx());
  }
  else
  {
      return DirectionH2<R>(dy(), -dx());
  }
}

template <class R >
inline
DirectionH2<R>
DirectionH2<R>::
transform(const Aff_transformationH2<R>& t) const
{ return t.transform(*this); }

template <class R >
CGAL_KERNEL_INLINE
VectorH2<R>
DirectionH2<R>::to_vector() const
{ return VectorH2<R>( x(), y() ); }


#ifndef CGAL_NO_OSTREAM_INSERT_DIRECTIONH2
template < class R >
std::ostream &
operator<<(std::ostream &os, const DirectionH2<R> &p)
{
  switch(os.iword(IO::mode))
  {
    case IO::ASCII :
        return os << p.dx() << ' ' << p.dy();
    case IO::BINARY :
        write(os, p.dx());
        write(os, p.dy());
        return os;
    default:
        return os << "DirectionH2(" << p.dx() << ", "
                                    << p.dy() << ')';
  }
}
#endif // CGAL_NO_OSTREAM_INSERT_DIRECTIONH2

#ifndef CGAL_NO_ISTREAM_EXTRACT_DIRECTIONH2
template < class R >
std::istream &
operator>>(std::istream &is, DirectionH2<R> &p)
{
  typename R::RT x, y;
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
        std::cerr << "" << std::endl;
        std::cerr << "Stream must be in ascii or binary mode" << std::endl;
        break;
  }
  p = DirectionH2<R>(x, y);
  return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_DIRECTIONH2

CGAL_END_NAMESPACE

#endif // CGAL_PVDH2_H
