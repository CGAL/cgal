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
// file          : include/CGAL/Homogeneous/PointH2.h
// package       : H2
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
// coordinator   : MPI, Saarbruecken 
// ======================================================================

#ifndef CGAL_HOMOGENEOUS_POINT_2_H
#define CGAL_HOMOGENEOUS_POINT_2_H

#include <CGAL/Origin.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Threetuple.h>
#include <CGAL/Kernel/Cartesian_coordinate_iterator_2.h>

CGAL_BEGIN_NAMESPACE

template < class R_ >
class PointH2
  : public R_::template Handle<Threetuple<typename R_::RT> >::type
{
CGAL_VC7_BUG_PROTECTED
  typedef typename R_::FT                   FT;
  typedef typename R_::RT                   RT;
  typedef typename R_::Vector_2             Vector_2;
  typedef typename R_::Point_2              Point_2;
  typedef typename R_::Direction_2          Direction_2;
  typedef typename R_::Aff_transformation_2 Aff_transformation_2;

  typedef Threetuple<RT>                           rep;
  typedef typename R_::template Handle<rep>::type  base;

public:
  typedef Cartesian_coordinate_iterator_2<R_> Cartesian_const_iterator;
  typedef R_                                    R;

    PointH2() 
       : base ( rep()) { }

    PointH2(const Origin &)  
       : base ( rep( RT(0), RT(0), RT(1))) { }

    PointH2(const PointH2<R> & p) 
       : base (p) { }

    PointH2(const Vector_2& v) 
       : base (v) { }

    PointH2(const RT& hx, const RT& hy )
      : base ( rep( hx, hy, RT(1) )) { }

    PointH2(const RT& hx, const RT& hy, const RT& hw)
    {
      if ( hw >= RT(0)   )
        initialize_with( rep( hx, hy, hw)); 
      else
        initialize_with( rep(-hx,-hy,-hw)); 
    }

    bool    operator==( const PointH2<R>& p) const;
    bool    operator!=( const PointH2<R>& p) const;

    const RT & hx() const { return Ptr()->e0; };
    const RT & hy() const { return Ptr()->e1; };
    const RT & hw() const { return Ptr()->e2; };

    FT      x()  const { return FT(hx()) / FT(hw()); };
    FT      y()  const { return FT(hy()) / FT(hw()); };

    FT      cartesian(int i)   const;
    FT      operator[](int i)  const;
    const RT & homogeneous(int i) const;

  Cartesian_const_iterator cartesian_begin() const 
  {
    return Cartesian_const_iterator(static_cast<const Point_2* const>(this), 0);
  }

  Cartesian_const_iterator cartesian_end() const 
  {
    return Cartesian_const_iterator(static_cast<const Point_2* const>(this), 2);
  }
    int     dimension() const;
    Bbox_2  bbox() const;

    PointH2<R> transform( const Aff_transformation_2 & t) const;
    Direction_2 direction() const;
};

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
const typename PointH2<R>::RT &
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
typename PointH2<R>::Direction_2
PointH2<R>::direction() const
{ return typename PointH2<R>::Direction_2(*this); }


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
PointH2<R>::transform(const typename PointH2<R>::Aff_transformation_2& t) const
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

CGAL_END_NAMESPACE

#endif // CGAL_HOMOGENEOUS_POINT_2_H
