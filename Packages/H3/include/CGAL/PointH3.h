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
// file          : PointH3.h
// package       : H3
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
// coordinator   : MPI, Saarbruecken  
// ======================================================================
 

#ifndef CGAL_HOMOGENEOUS_POINT_3_H
#define CGAL_HOMOGENEOUS_POINT_3_H

#include <CGAL/Origin.h>
#include <CGAL/Bbox_3.h>

CGAL_BEGIN_NAMESPACE

template < class R_ >
class PointH3
  : public R_::Point_handle_3
{
    typedef typename R_::RT     RT;
    typedef typename R_::FT     FT;
    typedef typename R_::Kernel_base::Vector_3             Vector_3;
    typedef typename R_::Kernel_base::Direction_3          Direction_3;
    typedef typename R_::Kernel_base::Aff_transformation_3 Aff_transformation_3;

    typedef typename R_::Point_handle_3            Point_handle_3_;
    typedef typename Point_handle_3_::element_type Point_ref_3;

public:
    typedef R_                 R;

  PointH3()
    : Point_handle_3_(Point_ref_3()) {}

  PointH3(const Origin &)
    : Point_handle_3_ (Point_ref_3( RT(0), RT(0), RT(0), RT(1) )) { }

  PointH3(const Vector_3& v)
    : Point_handle_3_(v) {}

  PointH3(const RT& x, const RT& y, const RT& z)
    : Point_handle_3_(Point_ref_3(x, y, z, RT(1))) {}

  PointH3(const RT& x, const RT& y, const RT& z, const RT& w)
  {
    if ( w < RT(0) )
      initialize_with( Point_ref_3(-x,-y,-z,-w)); 
    else
      initialize_with( Point_ref_3(x,y,z,w)); 
  }

  FT    x()  const;
  FT    y()  const;
  FT    z()  const;
  const RT & hx() const;
  const RT & hy() const;
  const RT & hz() const;
  const RT & hw() const;
  const RT & homogeneous(int i) const;
  FT    cartesian(int i) const;
  FT    operator[](int i) const;

  int   dimension() const;

  Direction_3
        direction() const;
  PointH3<R>
        transform( const Aff_transformation_3 & t) const;
  Bbox_3
        bbox() const;

  bool  operator==( const PointH3<R>& p) const;
  bool  operator!=( const PointH3<R>& p) const;
};


#ifdef CGAL_CFG_TYPENAME_BUG
#define typename
#endif

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
const typename PointH3<R>::RT &
PointH3<R>::hx() const
{ return  Ptr()->hx() ; }

template < class R >
inline
const typename PointH3<R>::RT &
PointH3<R>::hy() const
{ return  Ptr()->hy() ; }

template < class R >
inline
const typename PointH3<R>::RT &
PointH3<R>::hz() const
{ return  Ptr()->hz() ; }

template < class R >
inline
const typename PointH3<R>::RT &
PointH3<R>::hw() const
{ return  Ptr()->hw(); }

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
const typename PointH3<R>::RT &
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
typename PointH3<R>::Direction_3
PointH3<R>::direction() const
{ return Direction_3(*this); }

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
inline
PointH3<R>
PointH3<R>::transform(const typename PointH3<R>::Aff_transformation_3& t) const
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

CGAL_END_NAMESPACE

#endif // CGAL_HOMOGENEOUS_POINT_3_H
