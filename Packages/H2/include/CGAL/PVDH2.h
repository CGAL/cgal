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
// release_date  : 2000, October 11
// 
// source        : PointVectorDirectionH2.fw
// file          : PVDH2.h
// package       : H2 (2.13)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// revision      : 2.13
// revision_date : 11 Oct 2000 
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_PVDH2
#define CGAL_PVDH2

#define CGAL_POINTH2_H
#define CGAL_VECTORH2_H
#define CGAL_DIRECTIONH2_H

#include <CGAL/Threetuple.h>
#include <CGAL/homogeneous_classes.h>
#include <CGAL/Origin.h>

#include <CGAL/point_vector_declarationsH2.h>

CGAL_BEGIN_NAMESPACE

template < class FT, class RT >
class PointH2 : public Handle_for< Threetuple<RT> >
{

public:
            PointH2();
            PointH2(const Origin & o);
            PointH2(const PointH2<FT,RT> & p);
            PointH2(const RT& hx, const RT& hy );
            PointH2(const RT& hx, const RT& hy, const RT& hw );
            ~PointH2();

    bool    operator==( const PointH2<FT,RT>& p) const;
    bool    operator!=( const PointH2<FT,RT>& p) const;

    FT      x()  const;
    FT      y()  const;

    RT      hx() const;
    RT      hy() const;
    RT      hw() const;

    FT      cartesian(int i)   const;
    FT      operator[](int i)  const;
    RT      homogeneous(int i) const;

    // and for efficiency in the predicates:
    const RT&     hx_ref() const;
    const RT&     hy_ref() const;
    const RT&     hw_ref() const;

    int     dimension() const;
    Bbox_2  bbox() const;

    PointH2<FT,RT> transform( const Aff_transformationH2<FT,RT> & t) const;
    DirectionH2<FT,RT>
            direction() const;

            PointH2(const VectorH2<FT,RT>& v);
 friend CGAL_FRIEND_INLINE
        PointH2<FT,RT>
        CGAL_SCOPE origin_plus_vector CGAL_NULL_TMPL_ARGS ( const VectorH2<FT,RT> & v);
 friend CGAL_FRIEND_INLINE
        PointH2<FT,RT>
        CGAL_SCOPE origin_minus_vector CGAL_NULL_TMPL_ARGS (const VectorH2<FT,RT> & v);
};

template < class FT, class RT >
class VectorH2 : public Handle_for<Threetuple<RT> >
{
public:
            VectorH2();
            VectorH2(const VectorH2<FT,RT>& v);
            VectorH2(const Null_vector &);
            VectorH2(const RT& x, const RT& y);
            VectorH2(const RT& x, const RT& y, const RT& w );
            ~VectorH2();

    bool    operator==( const VectorH2<FT,RT>& v) const;
    bool    operator!=( const VectorH2<FT,RT>& v) const;
    bool    operator==( const Null_vector&) const;
    bool    operator!=( const Null_vector& v) const;

    FT      x()  const;
    FT      y()  const;
    RT      hx() const;
    RT      hy() const;
    RT      hw() const;
    FT      cartesian(int i)   const;
    RT      homogeneous(int i) const;
    FT      operator[](int i)  const;

    int     dimension() const;
    DirectionH2<FT,RT>
            direction() const;
    VectorH2<FT,RT>
            transform(const Aff_transformationH2<FT,RT>& t ) const;
    VectorH2<FT,RT>
            perpendicular(const Orientation& o ) const;

    VectorH2<FT,RT>
            operator-() const;

    VectorH2<FT,RT>
            opposite() const;

// undocumented:
            VectorH2(const DirectionH2<FT,RT> & dir);
protected:
            VectorH2(const PointH2<FT,RT> & p);
// friends:

friend CGAL_FRIEND_INLINE
       VectorH2<FT,RT>
       CGAL_SCOPE point_minus_origin CGAL_NULL_TMPL_ARGS ( const PointH2<FT,RT>& p);
friend CGAL_FRIEND_INLINE
       VectorH2<FT,RT>
       CGAL_SCOPE origin_minus_point CGAL_NULL_TMPL_ARGS ( const PointH2<FT,RT>& p);
/*
friend CGAL_FRIEND_INLINE
       PointH2<FT,RT>
       origin_minus_vector CGAL_NULL_TMPL_ARGS ( const VectorH2<FT,RT>& v);

friend CGAL_FRIEND_INLINE
       PointH2<FT,RT>
       operator+ CGAL_NULL_TMPL_ARGS ( const Origin &,
                                       const VectorH2<FT,RT>& v );
friend CGAL_KERNEL_FRIEND_INLINE
       VectorH2<FT,RT>
       operator- CGAL_NULL_TMPL_ARGS ( const VectorH2<FT,RT>& v1,
                                       const VectorH2<FT,RT>& v2 );
*/

friend CGAL_KERNEL_FRIEND_INLINE
       VectorH2<FT,RT>
       operator+ CGAL_NULL_TMPL_ARGS ( const VectorH2<FT,RT>&,
                                       const VectorH2<FT,RT>& );
friend CGAL_KERNEL_FRIEND_INLINE
       FT
       operator* CGAL_NULL_TMPL_ARGS ( const VectorH2<FT,RT>&,
                                       const VectorH2<FT,RT>& );
friend CGAL_KERNEL_FRIEND_INLINE
       VectorH2<FT,RT>
       operator* CGAL_NULL_TMPL_ARGS ( const VectorH2<FT,RT>&,
                                       const RT& );
friend CGAL_KERNEL_FRIEND_INLINE
       VectorH2<FT,RT>
       operator* CGAL_NULL_TMPL_ARGS ( const RT&,
                                       const VectorH2<FT,RT>& );
friend CGAL_KERNEL_FRIEND_INLINE
       VectorH2<FT,RT>
       operator/ CGAL_NULL_TMPL_ARGS ( const VectorH2<FT,RT>&,
                                       const RT& );
};

template < class FT, class RT >
class DirectionH2 : public Handle_for<Threetuple<RT> >
{
public:
            DirectionH2();
            DirectionH2(const DirectionH2<FT,RT>& d );
            DirectionH2(const PointH2<FT,RT> & p );
            DirectionH2(const VectorH2<FT,RT> & v );
            DirectionH2(const RT& x, const RT& y);
            DirectionH2(const RT& x, const RT& y, const RT& w );
            ~DirectionH2();

    bool    operator==( const DirectionH2<FT,RT>& d) const;
    bool    operator!=( const DirectionH2<FT,RT>& d) const;
    bool    operator< ( const DirectionH2<FT,RT>& d) const;
    bool    operator<=( const DirectionH2<FT,RT>& d) const;
    bool    operator> ( const DirectionH2<FT,RT>& d) const;
    bool    operator>=( const DirectionH2<FT,RT>& d) const;
    bool    counterclockwise_in_between( const DirectionH2<FT,RT>& d1,
                                         const DirectionH2<FT,RT>& d2 ) const;

    DirectionH2<FT,RT>
            operator-() const;


    VectorH2<FT,RT>    to_vector() const;
    // VectorH2<FT,RT>    vector() const;

    RT      x() const;
    RT      y() const;

    RT      delta(int i) const;
    RT      dx() const;
    RT      dy() const;

    DirectionH2<FT,RT>
          perpendicular(const Orientation &o) const;
    DirectionH2<FT,RT>
          transform(const Aff_transformationH2<FT,RT> &) const ;

};


template < class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
PointH2<FT,RT>::PointH2()
 : Handle_for< Threetuple<RT> >( Threetuple<RT>()) {}

template < class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
PointH2<FT,RT>::PointH2(const Origin&)
 : Handle_for< Threetuple<RT> >( Threetuple<RT>( RT(0), RT(0), RT(1))) {}

template < class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
PointH2<FT,RT>::PointH2(const PointH2<FT,RT>& p)
  : Handle_for<Threetuple<RT> >(p)
{}

template < class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
PointH2<FT,RT>::PointH2(const RT& hx, const RT& hy)
 : Handle_for< Threetuple<RT> >( Threetuple<RT>( hx, hy, RT(1) )) {}

template < class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
PointH2<FT,RT>::PointH2(const RT& hx, const RT& hy, const RT& hw)
{
   if ( hw >= RT(0)   )
   { initialize_with( Threetuple<RT>( hx, hy, hw)); }
   else
   { initialize_with( Threetuple<RT>(-hx,-hy,-hw)); }
}

template < class FT, class RT >
inline
PointH2<FT,RT>::~PointH2()
{}
template < class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
PointH2<FT,RT>::PointH2(const VectorH2<FT,RT>& v)
  : Handle_for<Threetuple<RT> >(v)
{}

template < class FT, class RT >
CGAL_KERNEL_INLINE
bool
PointH2<FT,RT>::operator==( const PointH2<FT,RT>& p) const
{
  return (  (hx() * p.hw() == p.hx() * hw() )
          &&(hy() * p.hw() == p.hy() * hw() ) );
}

template < class FT, class RT >
inline
bool
PointH2<FT,RT>::operator!=( const PointH2<FT,RT>& p) const
{ return !(*this == p); }   /* XXX */

template < class FT, class RT >
inline
FT
PointH2<FT,RT>::x()  const
{ return ( FT( hx() ) / FT( hw() )); }

template < class FT, class RT >
inline
FT
PointH2<FT,RT>::y()  const
{ return ( FT( hy() ) / FT( hw() )); }

template < class FT, class RT >
inline
RT
PointH2<FT,RT>::hx() const
{ return  ptr->e0 ; }

template < class FT, class RT >
inline
RT
PointH2<FT,RT>::hy() const
{ return  ptr->e1 ; }

template < class FT, class RT >
inline
RT
PointH2<FT,RT>::hw() const
{ return  ptr->e2; }

template < class FT, class RT >
CGAL_KERNEL_INLINE
FT
PointH2<FT,RT>::cartesian(int i) const
{
  CGAL_kernel_precondition( (i==0 || i==1) );
  switch (i)
  {
    case 0:  return x();
    case 1:  return y();
  }
  return FT(RT());
}

template < class FT, class RT >
CGAL_KERNEL_INLINE
RT
PointH2<FT,RT>::homogeneous(int i) const
{
  CGAL_kernel_precondition( (i>=0) && (i<=2) );
  switch (i)
  {
    case 0:   return hx();
    case 1:   return hy();
    case 2:   return hw();
  }
  return RT();
}

template < class FT, class RT >
inline
FT
PointH2<FT,RT>::operator[](int i) const
{ return cartesian(i); }

template < class FT, class RT >
inline
const RT&
PointH2<FT,RT>::hx_ref() const
{ return  ptr->e0 ; }

template < class FT, class RT >
inline
const RT&
PointH2<FT,RT>::hy_ref() const
{ return  ptr->e1 ; }

template < class FT, class RT >
inline
const RT&
PointH2<FT,RT>::hw_ref() const
{ return  ptr->e2; }


template < class FT, class RT >
inline
int
PointH2<FT,RT>::dimension() const
{ return 2; }

template < class FT, class RT >
CGAL_KERNEL_INLINE
DirectionH2<FT,RT>
PointH2<FT,RT>::direction() const
{ return DirectionH2<FT,RT>(*this); }



template < class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
VectorH2<FT,RT>::VectorH2()
 : Handle_for< Threetuple<RT> >( Threetuple<RT>()) {}

template < class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
VectorH2<FT,RT>::VectorH2(const Null_vector &)
 : Handle_for< Threetuple<RT> >( Threetuple<RT>(RT(0), RT(0), RT(1) )) {}

template < class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
VectorH2<FT,RT>::VectorH2(const VectorH2<FT,RT>& v)
  : Handle_for<Threetuple<RT> >(v)
{}

template < class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
VectorH2<FT,RT>::VectorH2(const RT& x, const RT& y)
 : Handle_for< Threetuple<RT> >( Threetuple<RT>( x,  y,  RT(1) )) {}

template < class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
VectorH2<FT,RT>::VectorH2(const RT& x, const RT& y, const RT& w)
{
  if ( w >= RT(0)   )
  { initialize_with( Threetuple<RT>( x,  y,  w)); }
  else
  { initialize_with( Threetuple<RT>(-x, -y, -w)); }
}

template < class FT, class RT >
inline
VectorH2<FT,RT>::~VectorH2()
{}

template < class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
VectorH2<FT,RT>::VectorH2(const PointH2<FT,RT> & p)
  : Handle_for<Threetuple<RT> >( p)
{}

template < class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
VectorH2<FT,RT>::VectorH2(const DirectionH2<FT,RT> & dir)
  : Handle_for<Threetuple<RT> >( dir)
{}

template < class FT, class RT >
inline
bool
VectorH2<FT,RT>::operator==( const Null_vector&) const
{ return (hx() == RT(0)) && (hy() == RT(0)); }

template < class FT, class RT >
inline
bool
VectorH2<FT,RT>::operator!=( const Null_vector& v) const
{ return !(*this == v); }

template < class FT, class RT >
CGAL_KERNEL_INLINE
bool
VectorH2<FT,RT>::operator==( const VectorH2<FT,RT>& v) const
{
  return (  (hx() * v.hw() == v.hx() * hw() )
          &&(hy() * v.hw() == v.hy() * hw() ) );
}

template < class FT, class RT >
inline
bool
VectorH2<FT,RT>::operator!=( const VectorH2<FT,RT>& v) const
{ return !(*this == v); }  /* XXX */

template < class FT, class RT >
inline
FT
VectorH2<FT,RT>::x()  const
{ return FT(ptr->e0 )/FT(ptr->e2 ) ; }


template < class FT, class RT >
inline
FT
VectorH2<FT,RT>::y()  const
{ return FT(ptr->e1 )/FT(ptr->e2 ) ; }

template < class FT, class RT >
inline
RT
VectorH2<FT,RT>::hx() const
{ return  ptr->e0 ; }

template < class FT, class RT >
inline
RT
VectorH2<FT,RT>::hy() const
{ return  ptr->e1 ; }

template < class FT, class RT >
inline
RT
VectorH2<FT,RT>::hw() const
{ return  ptr->e2 ; }

template < class FT, class RT >
CGAL_KERNEL_INLINE
FT
VectorH2<FT,RT>::cartesian(int i) const
{
  CGAL_kernel_precondition( (i==0 || i==1) );
  switch (i)
  {
    case 0:  return x();
    case 1:  return y();
  }
  return FT(RT());
}

template < class FT, class RT >
CGAL_KERNEL_INLINE
RT
VectorH2<FT,RT>::homogeneous(int i) const
{
  CGAL_kernel_precondition( (i>=0) && (i<=2) );
  switch (i)
  {
    case 0:  return hx();
    case 1:  return hy();
    case 2:  return hw();
  }
  return RT();
}

template < class FT, class RT >
inline
FT
VectorH2<FT,RT>::operator[](int i) const
{ return cartesian(i); }

template < class FT, class RT >
inline
int
VectorH2<FT,RT>::dimension() const
{ return 2; }

template < class FT, class RT >
CGAL_KERNEL_INLINE
DirectionH2<FT,RT>
VectorH2<FT,RT>::direction() const
{ return DirectionH2<FT,RT>(*this); }
template < class FT, class RT >
inline
VectorH2<FT,RT>
VectorH2<FT,RT>::operator-() const
{ return VectorH2<FT,RT>(- hx(), - hy(), hw() ); }

template < class FT, class RT >
inline
VectorH2<FT,RT>
VectorH2<FT,RT>::opposite() const
{ return VectorH2<FT,RT>(- hx(), - hy(), hw() ); }

template <class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
DirectionH2<FT,RT>::DirectionH2()
 : Handle_for< Threetuple<RT> >( Threetuple<RT>()) {}

template <class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
DirectionH2<FT,RT>::DirectionH2(const DirectionH2<FT,RT>& d )
  : Handle_for<Threetuple<RT> >( d )
{}

template <class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
DirectionH2<FT,RT>::DirectionH2(const PointH2<FT,RT> & p )
  : Handle_for<Threetuple<RT> >( p)
{}

template <class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
DirectionH2<FT,RT>::DirectionH2(const VectorH2<FT,RT> & v )
  : Handle_for<Threetuple<RT> >( v)
{}

template <class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
DirectionH2<FT,RT>::DirectionH2(const RT& x, const RT& y)
 : Handle_for< Threetuple<RT> >( Threetuple<RT>( x, y, RT(1) )) {}

template <class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
DirectionH2<FT,RT>::DirectionH2(const RT& x, const RT& y, const RT& w )
{
  if (w > RT(0)   )
  { initialize_with( Threetuple<RT>( x, y, w)); }
  else
  { initialize_with( Threetuple<RT>(-x,-y,-w)); }
}

template <class FT, class RT >
inline
DirectionH2<FT,RT>::~DirectionH2()
{}


template <class FT, class RT >
CGAL_KERNEL_INLINE
bool
DirectionH2<FT,RT>::operator==( const DirectionH2<FT,RT>& d) const
{
  return (  ( x() * d.y() == y() * d.x() )
          &&( CGAL_NTS sign( x() ) == CGAL_NTS sign( d.x() ) )
          &&( CGAL_NTS sign( y() ) == CGAL_NTS sign( d.y() ) ) );
}

template <class FT, class RT >
inline
bool
DirectionH2<FT,RT>::operator!=( const DirectionH2<FT,RT>& d) const
{ return !(*this == d); }

template <class FT, class RT >
inline
DirectionH2<FT,RT>
DirectionH2<FT,RT>::operator-() const
{ return DirectionH2<FT,RT>( - x(), - y() ); }
template <class FT, class RT >
inline
RT
DirectionH2<FT,RT>::dx() const
{ return ptr->e0; }

template <class FT, class RT >
inline
RT
DirectionH2<FT,RT>::dy() const
{ return ptr->e1; }

template <class FT, class RT >
CGAL_KERNEL_INLINE
RT
DirectionH2<FT,RT>::delta(int i) const
{
  CGAL_kernel_precondition( ( i == 0 ) || ( i == 1 ) );
  if (i == 0)
  {
      return dx();
  }
  return dy();
}

template <class FT, class RT >
inline
RT
DirectionH2<FT,RT>::x() const
{ return ptr->e0; }

template <class FT, class RT >
inline
RT
DirectionH2<FT,RT>::y() const
{ return ptr->e1; }


template <class FT, class RT>
CGAL_KERNEL_INLINE
VectorH2<FT,RT>
operator+(const VectorH2<FT,RT>& u, const VectorH2<FT,RT>& v)
{
  return VectorH2<FT,RT>( u.hx()*v.hw() + v.hx()*u.hw(),
                          u.hy()*v.hw() + v.hy()*u.hw(),
                          u.hw()*v.hw() );
}

template <class FT, class RT>
CGAL_KERNEL_INLINE
VectorH2<FT,RT>
operator-(const VectorH2<FT,RT>& u, const VectorH2<FT,RT>& v)
{
  return VectorH2<FT,RT>( u.hx()*v.hw() - v.hx()*u.hw(),
                          u.hy()*v.hw() - v.hy()*u.hw(),
                          u.hw()*v.hw() );
}

template <class FT, class RT>
CGAL_KERNEL_INLINE
FT
operator*(const VectorH2<FT,RT>& u, const VectorH2<FT,RT>& v)
{
  return (   FT( RT(u.hx()*v.hx() + u.hy()*v.hy()) )  /
             FT( RT(u.hw()*v.hw() ) )    );
}

template <class FT, class RT>
CGAL_KERNEL_INLINE
VectorH2<FT,RT>
operator/(const VectorH2<FT,RT>& v, const RT& f)
{ return VectorH2<FT,RT>( v.hx(), v.hy(), v.hw()*f ); }

template <class FT, class RT>
CGAL_KERNEL_INLINE
VectorH2<FT,RT>
operator*(const VectorH2<FT,RT>& v, const RT& f)
{ return VectorH2<FT,RT>( v.hx()*f, v.hy()*f, v.hw() ); }

template <class FT, class RT>
CGAL_KERNEL_INLINE
VectorH2<FT,RT>
operator*(const RT& f, const VectorH2<FT,RT>& v)
{ return VectorH2<FT,RT>( v.hx()*f, v.hy()*f, v.hw() ); }

template <class FT, class RT>
inline
PointH2<FT,RT>
origin_plus_vector(const VectorH2<FT,RT>& v)
{ return PointH2<FT,RT>( v ); }

template <class FT, class RT>
inline
PointH2<FT,RT>
operator+(const Origin&, const VectorH2<FT,RT>& v)
{ return origin_plus_vector( v ); }

template <class FT, class RT>
inline
PointH2<FT,RT>
origin_minus_vector(const VectorH2<FT,RT>& v)
{ return PointH2<FT,RT>( v.opposite() ); }

template <class FT, class RT>
inline
PointH2<FT,RT>
operator-(const Origin&, const VectorH2<FT,RT>& v)
{ return origin_minus_vector( v ); }

template <class FT, class RT>
inline
VectorH2<FT,RT>
point_minus_origin(const PointH2<FT,RT>& p)
{ return VectorH2<FT,RT>( p ); }

template <class FT, class RT>
inline
VectorH2<FT,RT>
operator-(const PointH2<FT,RT>& p, const Origin&)
{ return point_minus_origin( p ); }

template <class FT, class RT>
inline
VectorH2<FT,RT>
origin_minus_point(const PointH2<FT,RT>& p)
{ return  VectorH2<FT,RT>( p ).opposite(); }

template <class FT, class RT>
inline
VectorH2<FT,RT>
operator-(const Origin&, const PointH2<FT,RT>& p)
{ return  origin_minus_point( p ); }


template <class FT, class RT>
CGAL_KERNEL_INLINE
PointH2<FT,RT>
operator+(const PointH2<FT,RT>& p, const VectorH2<FT,RT>& v)
{
  return PointH2<FT,RT>( p.hx()*v.hw() + v.hx()*p.hw(),
                         p.hy()*v.hw() + v.hy()*p.hw(),
                         p.hw()*v.hw() );
}

template <class FT, class RT>
CGAL_KERNEL_INLINE
PointH2<FT,RT>
operator-(const PointH2<FT,RT>& p, const VectorH2<FT,RT>& v)
{
  return PointH2<FT,RT>( p.hx()*v.hw() - v.hx()*p.hw(),
                         p.hy()*v.hw() - v.hy()*p.hw(),
                         p.hw()*v.hw() );
}

template <class FT, class RT>
CGAL_KERNEL_INLINE
VectorH2<FT,RT>
operator-(const PointH2<FT,RT>& p, const PointH2<FT,RT>& q)
{
  return VectorH2<FT,RT>( p.hx()*q.hw() - q.hx()*p.hw(),
                          p.hy()*q.hw() - q.hy()*p.hw(),
                          p.hw()*q.hw() );
}

CGAL_END_NAMESPACE


#ifndef CGAL_PREDICATES_ON_DIRECTIONSH2_H
#include <CGAL/predicates_on_directionsH2.h>
#endif // CGAL_PREDICATES_ON_DIRECTIONSH2_H

CGAL_BEGIN_NAMESPACE

template <class FT, class RT >
inline
bool
DirectionH2<FT,RT>::operator< (const DirectionH2<FT,RT>& d) const
{ return (compare_angles_with_x_axis(*this,d) == SMALLER); }

template <class FT, class RT >
inline
bool
DirectionH2<FT,RT>::operator> (const DirectionH2<FT,RT>& d) const
{ return (compare_angles_with_x_axis(*this,d) == LARGER); }

template <class FT, class RT >
inline
bool
DirectionH2<FT,RT>::operator>= (const DirectionH2<FT,RT>& d) const
{ return !(compare_angles_with_x_axis(*this,d) == SMALLER); }

template <class FT, class RT >
inline
bool
DirectionH2<FT,RT>::operator<= (const DirectionH2<FT,RT>& d) const
{ return !(compare_angles_with_x_axis(*this,d) == LARGER); }

template <class FT, class RT >
CGAL_KERNEL_INLINE
bool
DirectionH2<FT,RT>::
counterclockwise_in_between( const DirectionH2<FT,RT>& d1,
                             const DirectionH2<FT,RT>& d2) const
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


#ifndef CGAL_AFF_TRANSFORMATIONH2_H
#include <CGAL/Aff_transformationH2.h>
#endif // CGAL_AFF_TRANSFORMATIONH2_H
#ifndef CGAL_BBOX_2_H
#include <CGAL/Bbox_2.h>
#endif // CGAL_BBOX_2_H
#ifndef CGAL_MISC_H
#include <CGAL/misc.h>
#endif // CGAL_MISC_H

CGAL_BEGIN_NAMESPACE

template < class FT, class RT >
CGAL_KERNEL_MEDIUM_INLINE
Bbox_2
PointH2<FT,RT>::bbox() const
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

template < class FT, class RT >
inline
PointH2<FT,RT>
PointH2<FT,RT>::transform(const Aff_transformationH2<FT,RT>& t) const
{ return t.transform(*this); }


#ifndef NO_OSTREAM_INSERT_POINTH2
template < class FT, class RT >
std::ostream &
operator<<(std::ostream &os, const PointH2<FT,RT> &p)
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
#endif // NO_OSTREAM_INSERT_POINTH2

#ifndef NO_ISTREAM_EXTRACT_POINTH2
template < class FT, class RT >
std::istream &
operator>>(std::istream &is, PointH2<FT,RT> &p)
{
  RT hx, hy, hw;
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
  p = PointH2<FT,RT>(hx, hy, hw);
  return is;
}
#endif // NO_ISTREAM_EXTRACT_POINTH2

template < class FT, class RT >
CGAL_KERNEL_INLINE
VectorH2<FT,RT>
VectorH2<FT,RT>::perpendicular(const Orientation& o) const
{
  CGAL_kernel_precondition(o != COLLINEAR);
  if (o == COUNTERCLOCKWISE)
  {
      return VectorH2<FT,RT>(-hy(), hx(), hw());
  }
  else
  {
      return VectorH2<FT,RT>(hy(), -hx(), hw());
  }
}

template < class FT, class RT >
inline
VectorH2<FT,RT>
VectorH2<FT,RT>::transform(const Aff_transformationH2<FT,RT>& t) const
{ return t.transform(*this); }


#ifndef NO_OSTREAM_INSERT_VECTORH2
template < class FT, class RT >
std::ostream &
operator<<(std::ostream &os, const VectorH2<FT,RT> &p)
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
#endif // NO_OSTREAM_INSERT_VECTORH2

#ifndef NO_ISTREAM_EXTRACT_VECTORH2
template < class FT, class RT >
std::istream &
operator>>(std::istream &is, VectorH2<FT,RT> &p)
{
  RT hx, hy, hw;
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
  p = VectorH2<FT,RT>(hx, hy, hw);
  return is;
}
#endif // NO_ISTREAM_EXTRACT_VECTORH2

template <class FT, class RT >
CGAL_KERNEL_INLINE
DirectionH2<FT,RT>
DirectionH2<FT,RT>::perpendicular(const Orientation& o) const
{
  CGAL_kernel_precondition(o != COLLINEAR);
  if (o == COUNTERCLOCKWISE)
  {
      return DirectionH2<FT,RT>(-dy(), dx());
  }
  else
  {
      return DirectionH2<FT,RT>(dy(), -dx());
  }
}

template <class FT, class RT >
inline
DirectionH2<FT,RT>
DirectionH2<FT,RT>::
transform(const Aff_transformationH2<FT,RT>& t) const
{ return t.transform(*this); }

template <class FT, class RT >
CGAL_KERNEL_INLINE
VectorH2<FT,RT>
DirectionH2<FT,RT>::to_vector() const
{ return VectorH2<FT,RT>( x(), y() ); }


#ifndef NO_OSTREAM_INSERT_DIRECTIONH2
template < class FT, class RT >
std::ostream &
operator<<(std::ostream &os, const DirectionH2<FT,RT> &p)
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
#endif // NO_OSTREAM_INSERT_DIRECTIONH2

#ifndef NO_ISTREAM_EXTRACT_DIRECTIONH2
template < class FT, class RT >
std::istream &
operator>>(std::istream &is, DirectionH2<FT,RT> &p)
{
  RT x, y;
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
  p = DirectionH2<FT,RT>(x, y);
  return is;
}
#endif // NO_ISTREAM_EXTRACT_DIRECTIONH2

CGAL_END_NAMESPACE



#endif // CGAL_PVDH2
