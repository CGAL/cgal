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
// file          : PVDH3.h
// package       : H3
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_PVDH3_H
#define CGAL_PVDH3_H
#define CGAL_POINTH3_H
#define CGAL_VECTORH3_H
#define CGAL_DIRECTIONH3_H

#include <CGAL/homogeneous_classes.h>
#include <CGAL/Origin.h>
#include <CGAL/RepH3.h>
#include <CGAL/Bbox_3.h>

CGAL_BEGIN_NAMESPACE

template <class R>
inline
PointH3<R>
operator+  (const Origin &, const VectorH3<R> & v);

template <class R>
inline
PointH3<R>
operator-  (const Origin &, const VectorH3<R> & v);

template <class R>
inline
VectorH3<R>
operator-   ( const PointH3<R> &, const Origin & );

template <class R>
inline
VectorH3<R>
operator-   ( const Origin &, const PointH3<R> & );

template <class R>
inline
PointH3<R>
operator+   ( const Origin &, const VectorH3<R> & );

template <class R>
inline
PointH3<R>
operator-   ( const Origin &, const VectorH3<R> & );

template <class R>
CGAL_KERNEL_INLINE
VectorH3<R>
operator+   ( const VectorH3<R> &, const VectorH3<R> & );

template <class R>
CGAL_KERNEL_INLINE
VectorH3<R>
operator-   ( const VectorH3<R> &, const VectorH3<R> & );

template <class R>
CGAL_KERNEL_INLINE
typename R::FT
operator*   ( const VectorH3<R> &, const VectorH3<R> & );

template <class R>
CGAL_KERNEL_INLINE
VectorH3<R>
operator*   ( const VectorH3<R> &, const typename R::RT & );

template <class R>
CGAL_KERNEL_INLINE
VectorH3<R>
operator*   ( const typename R::RT &, const VectorH3<R> & );

template <class R>
CGAL_KERNEL_INLINE
VectorH3<R>
operator/   ( const VectorH3<R> &, const typename R::RT & );

template <class R>
CGAL_KERNEL_INLINE
VectorH3<R>
cross_product( const VectorH3<R>& a, const VectorH3<R>& b);

template <class R>
CGAL_KERNEL_INLINE
DirectionH3<R>
cross_product( const DirectionH3<R>& a, const DirectionH3<R>& b);



template < class R_ >
class PointH3
  : public R_::Point_handle_3
{
  public:
    typedef R_                 R;
    typedef typename R::RT     RT;
    typedef typename R::FT     FT;

  PointH3();
  PointH3(const PointH3<R> & tbc);
  PointH3(const Origin &);
  PointH3(const VectorH3<R>& v);
  PointH3(const RT& x, const RT& y, const RT& z);
  PointH3(const RT& x, const RT& y, const RT& z, const RT& w);

  FT    x()  const;
  FT    y()  const;
  FT    z()  const;
  RT    hx() const;
  RT    hy() const;
  RT    hz() const;
  RT    hw() const;
  FT    cartesian(int i) const;
  RT    homogeneous(int i) const;
  FT    operator[](int i) const;

  int   dimension() const;

  DirectionH3<R>
        direction() const;
  PointH3<R>
        transform( const Aff_transformationH3<R> & t) const;
  Bbox_3
        bbox() const;

  bool  operator==( const PointH3<R>& p) const;
  bool  operator!=( const PointH3<R>& p) const;

  const RT&     hx_ref() const;
  const RT&     hy_ref() const;
  const RT&     hz_ref() const;
  const RT&     hw_ref() const;

friend CGAL_FRIEND_INLINE
       PointH3<R>
       operator+  CGAL_NULL_TMPL_ARGS (const Origin &,
                                       const VectorH3<R> & v);
friend CGAL_FRIEND_INLINE
       PointH3<R>
       operator-  CGAL_NULL_TMPL_ARGS (const Origin &,
                                       const VectorH3<R> & v);
};


template < class R_ >
class VectorH3
  : public R_::Vector_handle_3
{
  public:
    typedef R_                 R;
    typedef typename R::RT     RT;
    typedef typename R::FT     FT;

  VectorH3();
  VectorH3(const VectorH3<R> & tbc);
  VectorH3(const Null_vector&);
  VectorH3(const RT& x, const RT& y, const RT& z);
  VectorH3(const RT& w, const RT& x, const RT& y, const RT& z);

  FT    x()  const;
  FT    y()  const;
  FT    z()  const;
  RT    hx() const;
  RT    hy() const;
  RT    hz() const;
  RT    hw() const;
  FT    cartesian(int i) const;
  RT    homogeneous(int i) const;
  FT    operator[](int i) const;

  int   dimension() const;

  DirectionH3<R>
        direction() const;
  VectorH3<R>
        transform(const Aff_transformationH3<R>& t ) const;

  VectorH3<R>
        operator-() const;

  bool  operator==( const VectorH3<R>& v) const;
  bool  operator!=( const VectorH3<R>& v) const;

// undocumented:

  VectorH3(const PointH3<R> & p);
  VectorH3(const DirectionH3<R> & dir);   /* XXX */

// friends:

/*
friend CGAL_FRIEND_INLINE
       VectorH3<R>
       operator-  CGAL_NULL_TMPL_ARGS ( const PointH3<R> &,
                                        const Origin & );
friend CGAL_FRIEND_INLINE
       VectorH3<R>
       operator-  CGAL_NULL_TMPL_ARGS ( const Origin &,
                                        const PointH3<R> & );
friend CGAL_FRIEND_INLINE
       PointH3<R>
       operator-  CGAL_NULL_TMPL_ARGS ( const Origin &,
                                        const VectorH3<R> & );
friend CGAL_KERNEL_FRIEND_INLINE
       VectorH3<R>
       operator-  CGAL_NULL_TMPL_ARGS ( const VectorH3<R> &,
                                        const VectorH3<R> & );
*/
friend CGAL_FRIEND_INLINE
       PointH3<R>
       operator+  CGAL_NULL_TMPL_ARGS ( const Origin &,
                                        const VectorH3<R> & );
friend CGAL_KERNEL_FRIEND_INLINE
       VectorH3<R>
       operator+  CGAL_NULL_TMPL_ARGS ( const VectorH3<R> &,
                                        const VectorH3<R> & );
friend CGAL_KERNEL_FRIEND_INLINE
       FT
       operator*  CGAL_NULL_TMPL_ARGS ( const VectorH3<R> &,
                                        const VectorH3<R> & );
friend CGAL_KERNEL_FRIEND_INLINE
       VectorH3<R>
       operator*  CGAL_NULL_TMPL_ARGS ( const VectorH3<R> &,
                                        const RT & );
friend CGAL_KERNEL_FRIEND_INLINE
       VectorH3<R>
       operator*  CGAL_NULL_TMPL_ARGS ( const RT &,
                                        const VectorH3<R> & );
friend CGAL_KERNEL_FRIEND_INLINE
       VectorH3<R>
       operator/  CGAL_NULL_TMPL_ARGS ( const VectorH3<R> &,
                                        const RT & );
};

template < class R_ >
class DirectionH3
  : public R_::Direction_handle_3
{
  public:
    typedef R_                 R;
    typedef typename R::RT     RT;
    typedef typename R::FT     FT;

  DirectionH3();
  DirectionH3(const DirectionH3<R>& tbc );
  DirectionH3(const PointH3<R> & p );
  DirectionH3(const VectorH3<R> & v );
  DirectionH3(const RT& x, const RT& y,
                   const RT& z, const RT& w = RT(1) );

  DirectionH3<R>
        transform(const Aff_transformationH3<R> &) const ;
  DirectionH3<R>
        operator-() const;

  bool  is_degenerate() const;

  bool  operator==( const DirectionH3<R>& d) const;
  bool  operator!=( const DirectionH3<R>& d) const;

  VectorH3<R>    to_vector() const;
  // VectorH3<R>    vector() const { return to_vector(); }

  RT    dx() const;
  RT    dy() const;
  RT    dz() const;
  RT    x()  const;
  RT    y()  const;
  RT    z()  const;
  RT    hx() const;
  RT    hy() const;
  RT    hz() const;

  RT    delta(int i) const;

/*
friend CGAL_KERNEL_FRIEND_INLINE
DirectionH3<R>
cross_product CGAL_NULL_TMPL_ARGS (const DirectionH3<R>& d1,
                                   const DirectionH3<R>& d2);
*/
};



template < class R >
CGAL_KERNEL_CTOR_INLINE
PointH3<R>::PointH3()
 : Handle_for< RepH3<RT> >( RepH3<RT>())
{}

template < class R >
CGAL_KERNEL_CTOR_INLINE
PointH3<R>::PointH3(const Origin&)
{
 const RT RT0(0);
 const RT RT1(1);
 initialize_with( RepH3<RT>( RT0, RT0, RT0, RT1 ));
}

template < class R >
CGAL_KERNEL_CTOR_INLINE
PointH3<R>::PointH3(const PointH3<R>& tbc)
 : Handle_for< RepH3<RT> >(tbc)
{}

template < class R >
CGAL_KERNEL_CTOR_INLINE
PointH3<R>::PointH3(const VectorH3<R>& v)
 : Handle_for< RepH3<RT> >(v)
{}

template < class R >
CGAL_KERNEL_CTOR_INLINE
PointH3<R>::PointH3(const RT& x, const RT& y, const RT& z)
 : Handle_for< RepH3<RT> >( RepH3<RT>(x,y,z, RT(1)) )
{}

template < class R >
CGAL_KERNEL_CTOR_INLINE
PointH3<R>::PointH3(const RT& x, const RT& y, const RT& z,
                        const RT& w)
{
  if ( w < RT(0) )
  { initialize_with( RepH3<RT>(-x,-y,-z,-w)); }
  else
  { initialize_with( RepH3<RT>(x,y,z,w)); }
}

template < class R >
CGAL_KERNEL_INLINE
typename R::FT
PointH3<R>::x()  const
{ return ( FT(ptr->hx() ) / FT(ptr->hw() )); }


template < class R >
CGAL_KERNEL_INLINE
typename R::FT
PointH3<R>::y()  const
{ return ( FT(ptr->hy() ) / FT(ptr->hw() )); }

template < class R >
CGAL_KERNEL_INLINE
typename R::FT
PointH3<R>::z()  const
{ return ( FT(ptr->hz() ) / FT(ptr->hw() )); }

template < class R >
inline
typename R::RT
PointH3<R>::hx() const
{ return  ptr->hx() ; }

template < class R >
inline
typename R::RT
PointH3<R>::hy() const
{ return  ptr->hy() ; }

template < class R >
inline
typename R::RT
PointH3<R>::hz() const
{ return  ptr->hz() ; }

template < class R >
inline
typename R::RT
PointH3<R>::hw() const
{ return  ptr->hw(); }

template < class R >
inline
const typename R::RT&
PointH3<R>::hx_ref() const
{ return  ptr->e0 ; }

template < class R >
inline
const typename R::RT&
PointH3<R>::hy_ref() const
{ return  ptr->e1 ; }

template < class R >
inline
const typename R::RT&
PointH3<R>::hz_ref() const
{ return  ptr->e2 ; }

template < class R >
inline
const typename R::RT&
PointH3<R>::hw_ref() const
{ return  ptr->e3; }

template < class R >
inline
int
PointH3<R>::dimension() const
{ return 3; }

template < class R >
CGAL_KERNEL_INLINE
typename R::FT
PointH3<R>::cartesian(int i) const
{
  switch (i)
  {
      case 0:  return x();
      case 1:  return y();
      case 2:  return z();
      default: return cartesian( i%3 );
  }
}

template < class R >
CGAL_KERNEL_INLINE
typename R::RT
PointH3<R>::homogeneous(int i) const
{
  CGAL_kernel_precondition(i == 0 || i == 1 || i == 2 || i == 3);
  switch (i)
  {
     case 0:   return hx();
     case 1:   return hy();
     case 2:   return hz();
  }
  return hw();
}

template < class R >
inline
typename R::FT
PointH3<R>::operator[](int i) const
{ return cartesian(i); }


template < class R >
inline
DirectionH3<R>
PointH3<R>::direction() const
{ return DirectionH3<R>(*this); }
template < class R >
CGAL_KERNEL_INLINE
bool
PointH3<R>::operator==( const PointH3<R> & p) const
{
  return ( (hx() * p.hw() == p.hx() * hw() )
         &&(hy() * p.hw() == p.hy() * hw() )
         &&(hz() * p.hw() == p.hz() * hw() ) );
}

template < class R >
inline
bool
PointH3<R>::operator!=( const PointH3<R> & p) const
{ return !(*this == p); }


#ifndef CGAL_NO_OSTREAM_INSERT_POINTH3
template < class R >
std::ostream &operator<<(std::ostream &os, const PointH3<R> &p)
{
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << p.hx() << ' ' << p.hy() << ' ' << p.hz() << ' ' << p.hw();
    case IO::BINARY :
        write(os, p.hx());
        write(os, p.hy());
        write(os, p.hz());
        write(os, p.hw());
        return os;
    default:
        return os << "PointH3(" << p.hx() << ", "
                                << p.hy() << ", "
                                << p.hz() << ", "
                                << p.hw() << ')';
    }
}
#endif // CGAL_NO_OSTREAM_INSERT_POINTH3

#ifndef CGAL_NO_ISTREAM_EXTRACT_POINTH3
template < class R >
std::istream &operator>>(std::istream &is, PointH3<R> &p)
{
  typename R::RT hx, hy, hz, hw;
  switch(is.iword(IO::mode)) {
  case IO::ASCII :
        is >> hx >> hy >> hz >> hw;
        break;
  case IO::BINARY :
        read(is, hx);
        read(is, hy);
        read(is, hz);
        read(is, hw);
        break;
  default:
        std::cerr << "" << std::endl;
        std::cerr << "Stream must be in ascii or binary mode" << std::endl;
        break;
  }
  p = PointH3<R>(hx, hy, hz, hw);
  return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_POINTH3

template < class R >
CGAL_KERNEL_CTOR_INLINE
VectorH3<R>::VectorH3()
 : Handle_for< RepH3<RT> >( RepH3<RT>() )
{}

template < class R >
CGAL_KERNEL_CTOR_INLINE
VectorH3<R>::VectorH3(const VectorH3<R>& tbc)
 : Handle_for< RepH3<RT> >(tbc)
{}


template < class R >
CGAL_KERNEL_CTOR_INLINE
VectorH3<R>::VectorH3(const Null_vector&)
 : Handle_for< RepH3<RT> >( RepH3<RT>(RT(0), RT(0), RT(0), RT(1)) )
{}

template < class R >
CGAL_KERNEL_CTOR_INLINE
VectorH3<R>::VectorH3(const PointH3<R> & p)
 : Handle_for< RepH3<RT> >(p)
{}

template < class R >
CGAL_KERNEL_CTOR_INLINE
VectorH3<R>::VectorH3(const DirectionH3<R> & d)
 : Handle_for< RepH3<RT> >(d)
{}

template < class R >
CGAL_KERNEL_CTOR_INLINE
VectorH3<R>::VectorH3(const RT& x, const RT& y, const RT& z,
                          const RT& w)
{
  if ( w >= RT(0) )
  { initialize_with( RepH3<RT>(x, y, z, w)); }
  else
  { initialize_with( RepH3<RT>(-x,-y,-z,-w)); }
}

template < class R >
CGAL_KERNEL_CTOR_INLINE
VectorH3<R>::VectorH3(const RT& x, const RT& y, const RT& z)
 : Handle_for< RepH3<RT> >( RepH3<RT>(x, y, z, RT(1)) )
{}
template < class R >
CGAL_KERNEL_INLINE
typename R::FT
VectorH3<R>::x()  const
{ return FT(ptr->hx() )/FT(ptr->hw() ) ; }

template < class R >
CGAL_KERNEL_INLINE
typename R::FT
VectorH3<R>::y()  const
{ return FT(ptr->hy() )/FT(ptr->hw() ) ; }

template < class R >
CGAL_KERNEL_INLINE
typename R::FT
VectorH3<R>::z()  const
{ return FT(ptr->hz() )/FT(ptr->hw() ) ; }

template < class R >
inline
typename R::RT
VectorH3<R>::hx() const
{ return  ptr->hx() ; }

template < class R >
inline
typename R::RT
VectorH3<R>::hy() const
{ return  ptr->hy() ; }

template < class R >
inline
typename R::RT
VectorH3<R>::hz() const
{ return  ptr->hz() ; }

template < class R >
inline
typename R::RT
VectorH3<R>::hw() const
{ return  ptr->hw() ; }

template < class R >
inline
int
VectorH3<R>::dimension() const
{ return 3; }

template < class R >
CGAL_KERNEL_INLINE
typename R::FT
VectorH3<R>::cartesian(int i) const
{
  switch (i)
  {
      case 0:   return x();
      case 1:   return y();
      case 2:   return z();
      default:  return cartesian( i%3 );
  }
}

template < class R >
CGAL_KERNEL_INLINE
typename R::RT
VectorH3<R>::homogeneous(int i) const
{
  CGAL_kernel_precondition(i == 0 || i == 1 || i == 2 || i == 3);
  switch (i)
  {
      case 0:   return hx();
      case 1:   return hy();
      case 2:   return hz();
  }
  return hw() ;
}
template < class R >
inline
DirectionH3<R>
VectorH3<R>::direction() const
{ return DirectionH3<R>(*this); }

template < class R >
CGAL_KERNEL_INLINE
bool
VectorH3<R>::operator==( const VectorH3<R>& v) const
{
 return ( (hx() * v.hw() == v.hx() * hw() )
        &&(hy() * v.hw() == v.hy() * hw() )
        &&(hz() * v.hw() == v.hz() * hw() ) );
}

template < class R >
inline
bool
VectorH3<R>::operator!=( const VectorH3<R>& v) const
{ return !(*this == v); }

template < class R >
inline
typename R::FT
VectorH3<R>::operator[](int i) const
{ return cartesian(i); }

template < class R >
CGAL_KERNEL_INLINE
VectorH3<R>
VectorH3<R>::operator-() const
{ return VectorH3<R>( - hx(), - hy(), -hz(), hw() ); }


template <class R >
CGAL_KERNEL_CTOR_INLINE
DirectionH3<R>::DirectionH3()
 : Handle_for< RepH3<RT> >( RepH3<RT>() )
{}

template <class R >
CGAL_KERNEL_CTOR_INLINE
DirectionH3<R>::DirectionH3(const DirectionH3<R>& tbc )
 : Handle_for< RepH3<RT> >(tbc)
{}

template <class R >
CGAL_KERNEL_CTOR_INLINE
DirectionH3<R>::DirectionH3(const PointH3<R> & p )
 : Handle_for< RepH3<RT> >( p )
{}

template <class R >
CGAL_KERNEL_CTOR_INLINE
DirectionH3<R>::DirectionH3(const VectorH3<R> & v )
 : Handle_for< RepH3<RT> >( v )
{}

template <class R >
CGAL_KERNEL_CTOR_INLINE
DirectionH3<R>::DirectionH3(const RT& x, const RT& y, const RT& z,
                                const RT& w)
{
  if ( w >= RT(0) )
  { initialize_with( RepH3<RT>(x,y,z,w)); }
  else
  { initialize_with( RepH3<RT>(-x,-y,-z,-w)); }
}

template <class R >
CGAL_KERNEL_INLINE
typename R::RT
DirectionH3<R>::delta(int i) const
{
  switch (i)
  {
      case 0:  return x();
      case 1:  return y();
      case 2:  return z();
      default: return delta( i%3 );
  }
}

template <class R >
inline
typename R::RT
DirectionH3<R>::dx() const
{ return ptr->e0; }

template <class R >
inline
typename R::RT
DirectionH3<R>::x() const
{ return ptr->e0; }

template <class R >
inline
typename R::RT
DirectionH3<R>::hx() const
{ return ptr->e0; }

template <class R >
inline
typename R::RT
DirectionH3<R>::dy() const
{ return ptr->e1; }

template <class R >
inline
typename R::RT
DirectionH3<R>::y() const
{ return ptr->e1; }

template <class R >
inline
typename R::RT
DirectionH3<R>::hy() const
{ return ptr->e1; }

template <class R >
inline
typename R::RT
DirectionH3<R>::dz() const
{ return ptr->e2; }

template <class R >
inline
typename R::RT
DirectionH3<R>::z() const
{ return ptr->e2; }

template <class R >
inline
typename R::RT
DirectionH3<R>::hz() const
{ return ptr->e2; }
template <class R >
CGAL_KERNEL_INLINE
bool
DirectionH3<R>::operator==( const DirectionH3<R>& d) const
{
  return ( ( ptr->hx()*d.ptr->hy() == ptr->hy()*d.ptr->hx() )
         &&( ptr->hx()*d.ptr->hz() == ptr->hz()*d.ptr->hx() )
         &&( ptr->hy()*d.ptr->hz() == ptr->hz()*d.ptr->hy() )
         &&( CGAL_NTS sign( ptr->hx() ) == CGAL_NTS sign( d.ptr->hx() ) )
         &&( CGAL_NTS sign( ptr->hy() ) == CGAL_NTS sign( d.ptr->hy() ) )
         &&( CGAL_NTS sign( ptr->hz() ) == CGAL_NTS sign( d.ptr->hz() ) ) );
}

template <class R >
inline
bool
DirectionH3<R>::operator!=( const DirectionH3<R>& d) const
{ return !operator==(d); }

template <class R >
CGAL_KERNEL_INLINE
bool
DirectionH3<R>::is_degenerate() const
{ return ((hx() == RT(0)) && (hy() == RT(0)) && (hz() == RT(0))); }

template <class R >
inline
DirectionH3<R>
DirectionH3<R>::operator-() const
{ return DirectionH3<R>(- ptr->hx(),- ptr->hy(),- ptr->hz() ); }


template <class R >
inline
VectorH3<R>
DirectionH3<R>::to_vector() const
{ return VectorH3<R>(*this); }


template <class R>
CGAL_KERNEL_INLINE
VectorH3<R>
operator+(const VectorH3<R>& u, const VectorH3<R>& v)
{
  return VectorH3<R>(u.hx()*v.hw() + v.hx()*u.hw(),
                         u.hy()*v.hw() + v.hy()*u.hw(),
                         u.hz()*v.hw() + v.hz()*u.hw(),
                         u.hw()*v.hw() );
}

template <class R>
CGAL_KERNEL_INLINE
VectorH3<R>
operator-(const VectorH3<R>& u, const VectorH3<R>& v)
{
  return VectorH3<R>(u.hx()*v.hw() - v.hx()*u.hw(),
                         u.hy()*v.hw() - v.hy()*u.hw(),
                         u.hz()*v.hw() - v.hz()*u.hw(),
                         u.hw()*v.hw() );
}

template <class R>
CGAL_KERNEL_INLINE
typename R::FT
operator*(const VectorH3<R>& u, const VectorH3<R>& v)
{
  typedef typename R::RT RT;
  typedef typename R::FT FT;
  CGAL_kernel_assertion( u.hw() != RT(0) );
  CGAL_kernel_assertion( v.hw() != RT(0) );
  return ( FT( u.hx()*v.hx() + u.hy()*v.hy() + u.hz()*v.hz() ) /
           FT( u.hw()*v.hw() ) );
}

template <class R>
CGAL_KERNEL_INLINE
VectorH3<R>
operator/(const VectorH3<R>& v, const typename R::RT& f)
{ return VectorH3<R>( v.hx(), v.hy(), v.hz(), v.hw()*f ); }

template <class R>
CGAL_KERNEL_INLINE
VectorH3<R>
operator*(const VectorH3<R>& v, const typename R::RT& f)
{ return VectorH3<R>( v.hx()*f, v.hy()*f, v.hz()*f, v.hw() ); }

template <class R>
CGAL_KERNEL_INLINE
VectorH3<R>
operator*(const typename R::RT& f, const VectorH3<R>& v)
{ return VectorH3<R>( v.hx()*f, v.hy()*f, v.hz()*f, v.hw() ); }

template <class R>
CGAL_KERNEL_INLINE
VectorH3<R>
cross_product(const VectorH3<R>& a, const VectorH3<R>& b)
{
 return VectorH3<R>(a.hy()*b.hz() - a.hz()*b.hy(),
                        a.hz()*b.hx() - a.hx()*b.hz(),
                        a.hx()*b.hy() - a.hy()*b.hx(),
                        a.hw()*b.hw() );
}

template <class R>
inline
PointH3<R>
operator+(const Origin& , const VectorH3<R>& v)
{ return PointH3<R>( v ); }

template <class R>
inline
PointH3<R>
operator-(const Origin& , const VectorH3<R>& v)
{ return  PointH3<R>(-v ); }

template <class R>
inline
VectorH3<R>
operator-(const PointH3<R>& p, const Origin& )
{ return VectorH3<R>( p ); }

template <class R>
inline
VectorH3<R>
operator-(const Origin& , const PointH3<R>& p)
{ return  - VectorH3<R>( p ); }


template <class R>
CGAL_KERNEL_INLINE
PointH3<R>
operator+(const PointH3<R>& p, const VectorH3<R>& v)
{
  return PointH3<R>(p.hx()*v.hw() + v.hx()*p.hw(),
                        p.hy()*v.hw() + v.hy()*p.hw(),
                        p.hz()*v.hw() + v.hz()*p.hw(),
                        p.hw()*v.hw() );
}

template <class R>
CGAL_KERNEL_INLINE
PointH3<R>
operator-(const PointH3<R>& p, const VectorH3<R>& v)
{
  return PointH3<R>( p.hx()*v.hw() - v.hx()*p.hw(),
                              p.hy()*v.hw() - v.hy()*p.hw(),
                              p.hz()*v.hw() - v.hz()*p.hw(),
                              p.hw()*v.hw() );
}

template <class R>
CGAL_KERNEL_INLINE
VectorH3<R>
operator-(const PointH3<R>& p, const PointH3<R>& q)
{
  return PointH3<R>( p.hx()*q.hw() - q.hx()*p.hw(),
                              p.hy()*q.hw() - q.hy()*p.hw(),
                              p.hz()*q.hw() - q.hz()*p.hw(),
                              p.hw()*q.hw() );
}



template <class R>
CGAL_KERNEL_INLINE
DirectionH3<R>
cross_product( const DirectionH3<R>& d1,
               const DirectionH3<R>& d2)
{ return cross_product(d1.to_vector(),d2.to_vector()).direction(); }

CGAL_END_NAMESPACE

#include <CGAL/Aff_transformationH3.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/misc.h>

CGAL_BEGIN_NAMESPACE

template < class R >
inline
PointH3<R>
PointH3<R>::transform(const Aff_transformationH3<R>& t) const
{ return t.transform(*this); }

template < class R >
CGAL_KERNEL_LARGE_INLINE
Bbox_3
PointH3<R>::bbox() const
{
#ifndef CGAL_CFG_NO_NAMESPACE
  using std::swap;
#endif // CGAL_CFG_NO_NAMESPACE

  // double bx = to_double(x());
  // double by = to_double(y());
  // double bz = to_double(z());
  // return Bbox_3(bx, by, bz, bx, by, bz);

  double eps  = double(1.0) /(double(1<<26) * double(1<<26));
  double hxd  = CGAL::to_double( hx() );
  double hyd  = CGAL::to_double( hy() );
  double hzd  = CGAL::to_double( hz() );
  double hwd  = CGAL::to_double( hw() );
  double xmin = ( hxd - eps*hxd ) / ( hwd + eps*hwd );
  double xmax = ( hxd + eps*hxd ) / ( hwd - eps*hwd );
  double ymin = ( hyd - eps*hyd ) / ( hwd + eps*hwd );
  double ymax = ( hyd + eps*hyd ) / ( hwd - eps*hwd );
  double zmin = ( hzd - eps*hzd ) / ( hwd + eps*hwd );
  double zmax = ( hzd + eps*hzd ) / ( hwd - eps*hwd );
  if ( hx() < RT(0)   )
  {
      swap(xmin, xmax);
  }
  if ( hy() < RT(0)   )
  {
      swap(ymin, ymax);
  }
  if ( hz() < RT(0)   )
  {
      swap(zmin, zmax);
  }
  return Bbox_3(xmin, ymin, zmin, xmax, ymax, zmax);
}

template < class R >
inline
VectorH3<R>
VectorH3<R>::transform(const Aff_transformationH3<R>&t ) const
{ return t.transform(*this); }


#ifndef CGAL_NO_OSTREAM_INSERT_VECTORH3
template < class R >
std::ostream& operator<<(std::ostream& os, const VectorH3<R>& v)
{
  switch(os.iword(IO::mode))
  {
    case IO::ASCII :
        return os << v.hx() << ' ' << v.hy() << ' ' << v.hz() << ' ' << v.hw();
    case IO::BINARY :
        write(os, v.hx());
        write(os, v.hy());
        write(os, v.hz());
        write(os, v.hw());
        return os;
    default:
        return os << "VectorH3(" << v.hx() << ", "
                                 << v.hy() << ", "
                                 << v.hz() << ", "
                                 << v.hw() << ')';
  }
}
#endif // CGAL_NO_OSTREAM_INSERT_VECTORH3

#ifndef CGAL_NO_ISTREAM_EXTRACT_VECTORH3
template < class R >
std::istream& operator>>(std::istream& is, VectorH3<R>& v)
{
  typename R::RT hx, hy, hz, hw;
  switch(is.iword(IO::mode))
  {
    case IO::ASCII :
        is >> hx >> hy >> hz >> hw;
        break;
    case IO::BINARY :
        read(is, hx);
        read(is, hy);
        read(is, hz);
        read(is, hw);
        break;
    default:
        std::cerr << "" << std::endl;
        std::cerr << "Stream must be in ascii or binary mode" << std::endl;
        break;
  }
  v = VectorH3<R>(hx, hy, hz, hw);
  return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_VECTORH3

template <class R >
inline
DirectionH3<R>
DirectionH3<R>::
transform(const Aff_transformationH3<R>& t) const
{ return t.transform(*this); }


#ifndef CGAL_NO_OSTREAM_INSERT_DIRECTIONH3
template < class R >
std::ostream &operator<<(std::ostream &os, const DirectionH3<R> &p)
{
  switch(os.iword(IO::mode))
  {
    case IO::ASCII :
        return os << p.dx() << ' ' << p.dy() << ' ' << p.dz();
    case IO::BINARY :
        write(os, p.dx());
        write(os, p.dy());
        write(os, p.dz());
        return os;
    default:
        return os << "DirectionH3(" << p.dx() << ", "
                                    << p.dy() << ", "
                                    << p.dz() << ')';
  }
}
#endif // CGAL_NO_OSTREAM_INSERT_DIRECTIONH3

#ifndef CGAL_NO_ISTREAM_EXTRACT_DIRECTIONH3
template < class R >
std::istream &operator>>(std::istream &is, DirectionH3<R> &p)
{
  typename R::RT x, y, z;
  switch(is.iword(IO::mode))
  {
    case IO::ASCII :
        is >> x >> y >> z;
        break;
    case IO::BINARY :
        read(is, x);
        read(is, y);
        read(is, z);
        break;
    default:
        std::cerr << "" << std::endl;
        std::cerr << "Stream must be in ascii or binary mode" << std::endl;
        break;
  }
  p = DirectionH3<R>(x, y, z);
  return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_DIRECTIONH3

/*
template <class R>
VectorH3<R>
operator/(const VectorH3<R>& v, const RT& f);
*/

CGAL_END_NAMESPACE

#endif // CGAL_PVDH3_H
