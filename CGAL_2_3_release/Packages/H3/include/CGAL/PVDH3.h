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
#include <CGAL/Bbox_3.h>

CGAL_BEGIN_NAMESPACE

template < class R_ >
class PointH3
  : public R_::Point_handle_3
{
  public:
    typedef R_                 R;
    typedef typename R::RT     RT;
    typedef typename R::FT     FT;

    typedef typename R::Point_handle_3       Point_handle_3_;
    typedef typename Point_handle_3_::element_type Point_ref_3;

  PointH3()
    : Point_handle_3_(Point_ref_3()) {}

  PointH3(const Origin &);

  PointH3(const VectorH3<R>& v)
    : Point_handle_3_(v) {}

  PointH3(const RT& x, const RT& y, const RT& z)
    : Point_handle_3_(Point_ref_3(x, y, z, RT(1))) {}

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
};


template < class R_ >
class VectorH3
  : public R_::Vector_handle_3
{
  public:
    typedef R_                 R;
    typedef typename R::RT     RT;
    typedef typename R::FT     FT;

  typedef typename R::Vector_handle_3           Vector_handle_3_;
  typedef typename Vector_handle_3_::element_type       Vector_ref_3;

  VectorH3()
    : Vector_handle_3_(Vector_ref_3()) {}

  VectorH3(const PointH3<R>& a, const PointH3<R>& b)
    : Vector_handle_3_(b-a) {}

  VectorH3(const Null_vector&)
    : Vector_handle_3_(Vector_ref_3(RT(0), RT(0), RT(0), RT(1))) {}

  VectorH3(const RT& x, const RT& y, const RT& z)
    : Vector_handle_3_(Vector_ref_3(x, y, z, RT(1))) {}

  VectorH3(const RT& w, const RT& x, const RT& y, const RT& z);

// undocumented:

  VectorH3(const PointH3<R> & p)
    : Vector_handle_3_(p) {}

  VectorH3(const DirectionH3<R> & d)   /* XXX */
    : Vector_handle_3_(d) {}

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

  FT    operator*( const VectorH3<R>& v) const;
};

template < class R_ >
class DirectionH3
  : public R_::Direction_handle_3
{
  public:
    typedef R_                 R;
    typedef typename R::RT     RT;
    typedef typename R::FT     FT;

    typedef typename R::Direction_handle_3   Direction_handle_3_;
    typedef typename Direction_handle_3_::element_type Direction_ref_3;

  DirectionH3()
    : Direction_handle_3_(Direction_ref_3()) {}

  DirectionH3(const PointH3<R> & p )
    : Direction_handle_3_(p) {}

  DirectionH3(const VectorH3<R> & v )
    : Direction_handle_3_(v) {}

  DirectionH3(const LineH3<R> & l )
    : Direction_handle_3_(l.direction()) {}

  DirectionH3(const RayH3<R> & r )
    : Direction_handle_3_(r.direction()) {}

  DirectionH3(const SegmentH3<R> & s )
    : Direction_handle_3_(s.direction()) {}

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
};

#ifdef CGAL_CFG_TYPENAME_BUG
#define typename
#endif

template < class R >
CGAL_KERNEL_CTOR_INLINE
PointH3<R>::PointH3(const Origin&)
{
 const RT RT0(0);
 const RT RT1(1);
 initialize_with( Point_ref_3( RT0, RT0, RT0, RT1 ));
}

template < class R >
CGAL_KERNEL_CTOR_INLINE
PointH3<R>::PointH3(const RT& x, const RT& y, const RT& z,
                        const RT& w)
{
  if ( w < RT(0) )
  { initialize_with( Point_ref_3(-x,-y,-z,-w)); }
  else
  { initialize_with( Point_ref_3(x,y,z,w)); }
}

template < class R >
CGAL_KERNEL_INLINE
typename PointH3<R>::FT
PointH3<R>::x()  const
{ return ( FT(Ptr()->hx() ) / FT(Ptr()->hw() )); }

template < class R >
CGAL_KERNEL_INLINE
typename PointH3<R>::FT
PointH3<R>::y()  const
{ return ( FT(Ptr()->hy() ) / FT(Ptr()->hw() )); }

template < class R >
CGAL_KERNEL_INLINE
typename PointH3<R>::FT
PointH3<R>::z()  const
{ return ( FT(Ptr()->hz() ) / FT(Ptr()->hw() )); }

template < class R >
inline
typename PointH3<R>::RT
PointH3<R>::hx() const
{ return  Ptr()->hx() ; }

template < class R >
inline
typename PointH3<R>::RT
PointH3<R>::hy() const
{ return  Ptr()->hy() ; }

template < class R >
inline
typename PointH3<R>::RT
PointH3<R>::hz() const
{ return  Ptr()->hz() ; }

template < class R >
inline
typename PointH3<R>::RT
PointH3<R>::hw() const
{ return  Ptr()->hw(); }

template < class R >
inline
const typename PointH3<R>::RT&
PointH3<R>::hx_ref() const
{ return  Ptr()->e0 ; }

template < class R >
inline
const typename PointH3<R>::RT&
PointH3<R>::hy_ref() const
{ return  Ptr()->e1 ; }

template < class R >
inline
const typename PointH3<R>::RT&
PointH3<R>::hz_ref() const
{ return  Ptr()->e2 ; }

template < class R >
inline
const typename PointH3<R>::RT&
PointH3<R>::hw_ref() const
{ return  Ptr()->e3; }

template < class R >
inline
int
PointH3<R>::dimension() const
{ return 3; }

template < class R >
CGAL_KERNEL_INLINE
typename PointH3<R>::FT
PointH3<R>::cartesian(int i) const
{
  CGAL_kernel_precondition(i == 0 || i == 1 || i == 2);
  switch (i)
  {
      case 0:  return x();
      case 1:  return y();
  }
  return z();
}

template < class R >
CGAL_KERNEL_INLINE
typename PointH3<R>::RT
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
typename PointH3<R>::FT
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
VectorH3<R>::VectorH3(const RT& x, const RT& y, const RT& z, const RT& w)
{
  if ( w >= RT(0) )
  { initialize_with( Vector_ref_3(x, y, z, w)); }
  else
  { initialize_with( Vector_ref_3(-x,-y,-z,-w)); }
}

template < class R >
CGAL_KERNEL_INLINE
typename VectorH3<R>::FT
VectorH3<R>::x()  const
{ return FT(Ptr()->hx() )/FT(Ptr()->hw() ) ; }

template < class R >
CGAL_KERNEL_INLINE
typename VectorH3<R>::FT
VectorH3<R>::y()  const
{ return FT(Ptr()->hy() )/FT(Ptr()->hw() ) ; }

template < class R >
CGAL_KERNEL_INLINE
typename VectorH3<R>::FT
VectorH3<R>::z()  const
{ return FT(Ptr()->hz() )/FT(Ptr()->hw() ) ; }

template < class R >
inline
typename VectorH3<R>::RT
VectorH3<R>::hx() const
{ return  Ptr()->hx() ; }

template < class R >
inline
typename VectorH3<R>::RT
VectorH3<R>::hy() const
{ return  Ptr()->hy() ; }

template < class R >
inline
typename VectorH3<R>::RT
VectorH3<R>::hz() const
{ return  Ptr()->hz() ; }

template < class R >
inline
typename VectorH3<R>::RT
VectorH3<R>::hw() const
{ return  Ptr()->hw() ; }

template < class R >
inline
int
VectorH3<R>::dimension() const
{ return 3; }

template < class R >
CGAL_KERNEL_INLINE
typename VectorH3<R>::FT
VectorH3<R>::cartesian(int i) const
{
  CGAL_kernel_precondition(i == 0 || i == 1 || i == 2);
  switch (i)
  {
      case 0:   return x();
      case 1:   return y();
  }
  return z();
}

template < class R >
CGAL_KERNEL_INLINE
typename VectorH3<R>::RT
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
typename VectorH3<R>::FT
VectorH3<R>::operator[](int i) const
{ return cartesian(i); }

template < class R >
CGAL_KERNEL_INLINE
VectorH3<R>
VectorH3<R>::operator-() const
{ return VectorH3<R>( - hx(), - hy(), -hz(), hw() ); }

template <class R >
CGAL_KERNEL_CTOR_INLINE
DirectionH3<R>::DirectionH3(const RT& x, const RT& y, const RT& z,
                                const RT& w)
{
  if ( w >= RT(0) )
  { initialize_with( Direction_ref_3(x,y,z,w)); }
  else
  { initialize_with( Direction_ref_3(-x,-y,-z,-w)); }
}

template <class R >
CGAL_KERNEL_INLINE
typename DirectionH3<R>::RT
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
typename DirectionH3<R>::RT
DirectionH3<R>::dx() const
{ return Ptr()->e0; }

template <class R >
inline
typename DirectionH3<R>::RT
DirectionH3<R>::x() const
{ return Ptr()->e0; }

template <class R >
inline
typename DirectionH3<R>::RT
DirectionH3<R>::hx() const
{ return Ptr()->e0; }

template <class R >
inline
typename DirectionH3<R>::RT
DirectionH3<R>::dy() const
{ return Ptr()->e1; }

template <class R >
inline
typename DirectionH3<R>::RT
DirectionH3<R>::y() const
{ return Ptr()->e1; }

template <class R >
inline
typename DirectionH3<R>::RT
DirectionH3<R>::hy() const
{ return Ptr()->e1; }

template <class R >
inline
typename DirectionH3<R>::RT
DirectionH3<R>::dz() const
{ return Ptr()->e2; }

template <class R >
inline
typename DirectionH3<R>::RT
DirectionH3<R>::z() const
{ return Ptr()->e2; }

template <class R >
inline
typename DirectionH3<R>::RT
DirectionH3<R>::hz() const
{ return Ptr()->e2; }

template <class R >
CGAL_KERNEL_INLINE
bool
DirectionH3<R>::operator==( const DirectionH3<R>& d) const
{
  return ( ( Ptr()->hx()*d.Ptr()->hy() == Ptr()->hy()*d.Ptr()->hx() )
        &&( Ptr()->hx()*d.Ptr()->hz() == Ptr()->hz()*d.Ptr()->hx() )
        &&( Ptr()->hy()*d.Ptr()->hz() == Ptr()->hz()*d.Ptr()->hy() )
        &&( CGAL_NTS sign( Ptr()->hx() ) == CGAL_NTS sign( d.Ptr()->hx() ) )
        &&( CGAL_NTS sign( Ptr()->hy() ) == CGAL_NTS sign( d.Ptr()->hy() ) )
        &&( CGAL_NTS sign( Ptr()->hz() ) == CGAL_NTS sign( d.Ptr()->hz() ) ) );
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
{ return DirectionH3<R>(- Ptr()->hx(),- Ptr()->hy(),- Ptr()->hz() ); }

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
typename VectorH3<R>::FT
VectorH3<R>::operator*(const VectorH3<R>& v) const
{
  typedef typename R::RT RT;
  typedef typename R::FT FT;
  CGAL_kernel_assertion( hw() != RT(0) );
  CGAL_kernel_assertion( hw() != RT(0) );
  return ( FT( hx()*v.hx() + hy()*v.hy() + hz()*v.hz() ) /
           FT( hw()*v.hw() ) );
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

#ifdef CGAL_CFG_TYPENAME_BUG
#undef typename
#endif

CGAL_END_NAMESPACE

#include <CGAL/Aff_transformationH3.h>

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

CGAL_END_NAMESPACE

#endif // CGAL_PVDH3_H
