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
// file          : include/CGAL/Homogeneous/VectorH3.h
// package       : H3
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
// coordinator   : MPI, Saarbruecken  
// ======================================================================
 
#ifndef CGAL_HOMOGENEOUS_VECTOR_3_H
#define CGAL_HOMOGENEOUS_VECTOR_3_H

#include <CGAL/Origin.h>
#include <CGAL/Fourtuple.h>

CGAL_BEGIN_NAMESPACE

template < class R_ >
class VectorH3
  : public R_::template Handle<Fourtuple<typename R_::RT> >::type
{
CGAL_VC7_BUG_PROTECTED
  typedef typename R_::RT                   RT;
  typedef typename R_::FT                   FT;
  typedef typename R_::Point_3              Point_3;
  typedef typename R_::Direction_3          Direction_3;
  typedef typename R_::Vector_3             Vector_3;
  typedef typename R_::Aff_transformation_3 Aff_transformation_3;

  typedef Fourtuple<RT>                            rep;
  typedef typename R_::template Handle<rep>::type  base;

public:
  typedef R_                 R;

  VectorH3()
    : base(rep()) {}

  VectorH3(const Point_3& a, const Point_3& b)
    : base(b-a) {}

  VectorH3(const Null_vector&)
    : base(rep(RT(0), RT(0), RT(0), RT(1))) {}

  VectorH3(const RT& x, const RT& y, const RT& z)
    : base(rep(x, y, z, RT(1))) {}

  VectorH3(const RT& w, const RT& x, const RT& y, const RT& z);

// undocumented:

  VectorH3(const Point_3 & p)
    : base(p) {}

  VectorH3(const Direction_3 & d)   /* XXX */
    : base(d) {}

  const RT & hx() const { return  Ptr()->e0 ; }
  const RT & hy() const { return  Ptr()->e1 ; }
  const RT & hz() const { return  Ptr()->e2 ; }
  const RT & hw() const { return  Ptr()->e3 ; }
  FT    x()  const { return FT(hx())/FT(hw()) ; }
  FT    y()  const { return FT(hy())/FT(hw()) ; }
  FT    z()  const { return FT(hz())/FT(hw()) ; }
  const RT & homogeneous(int i) const;
  FT    cartesian(int i) const;
  FT    operator[](int i) const;

  int   dimension() const { return 3; };

  Direction_3 direction() const;
  Vector_3 transform(const Aff_transformation_3& t ) const;

  Vector_3 operator-() const;

  bool  operator==( const VectorH3<R>& v) const;
  bool  operator!=( const VectorH3<R>& v) const;

  Vector_3 operator+( const VectorH3 &v) const;
  Vector_3 operator-( const VectorH3 &v) const;
  FT       operator*( const VectorH3 &v) const;
  Vector_3 operator*( const RT &f) const;
  Vector_3 operator*( const FT &f) const;
  Vector_3 operator/( const RT &f) const;
  Vector_3 operator/( const FT &f) const;
};

template < class R >
CGAL_KERNEL_INLINE
VectorH3<R>::VectorH3(const RT& x, const RT& y, const RT& z, const RT& w)
{
  if ( w >= RT(0) )
  { initialize_with( rep(x, y, z, w)); }
  else
  { initialize_with( rep(-x,-y,-z,-w)); }
}


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
const typename VectorH3<R>::RT &
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
typename VectorH3<R>::Direction_3
VectorH3<R>::direction() const
{ return Direction_3(*this); }

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
typename VectorH3<R>::Vector_3
VectorH3<R>::operator-() const
{ return Vector_3( - hx(), - hy(), -hz(), hw() ); }

template <class R>
CGAL_KERNEL_INLINE
typename VectorH3<R>::Vector_3
VectorH3<R>::operator+(const VectorH3<R>& v) const
{
  return VectorH3<R>(hx()*v.hw() + v.hx()*hw(),
                     hy()*v.hw() + v.hy()*hw(),
                     hz()*v.hw() + v.hz()*hw(),
                     hw()*v.hw() );
}

template <class R>
CGAL_KERNEL_INLINE
typename VectorH3<R>::Vector_3
VectorH3<R>::operator-(const VectorH3<R>& v) const
{
  return VectorH3<R>(hx()*v.hw() - v.hx()*hw(),
                     hy()*v.hw() - v.hy()*hw(),
                     hz()*v.hw() - v.hz()*hw(),
                     hw()*v.hw() );
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
typename VectorH3<R>::Vector_3
VectorH3<R>::operator/(const typename VectorH3<R>::RT& f) const
{ return VectorH3<R>( hx(), hy(), hz(), hw()*f ); }

template <class R>
CGAL_KERNEL_INLINE
typename VectorH3<R>::Vector_3
VectorH3<R>::operator/(const typename VectorH3<R>::FT& f) const
{ return VectorH3<R>( hx()*f.denominator(), hy()*f.denominator(),
		      hz()*f.denominator(), hw()*f.numerator() ); }

template <class R>
CGAL_KERNEL_INLINE
typename VectorH3<R>::Vector_3
VectorH3<R>::operator*(const typename VectorH3<R>::RT& f) const
{ return VectorH3<R>( hx()*f, hy()*f, hz()*f, hw() ); }

template <class R>
CGAL_KERNEL_INLINE
typename VectorH3<R>::Vector_3
VectorH3<R>::operator*(const typename VectorH3<R>::FT& f) const
{ return VectorH3<R>( hx()*f.numerator(), hy()*f.numerator(),
		      hz()*f.numerator(), hw()*f.denominator() ); }

template <class R>
CGAL_KERNEL_INLINE
typename R::Vector_3
operator*(const typename R::RT& f, const VectorH3<R>& v)
{ return VectorH3<R>( v.hx()*f, v.hy()*f, v.hz()*f, v.hw() ); }

template <class R>
CGAL_KERNEL_INLINE
typename R::Vector_3
cross_product(const VectorH3<R>& a, const VectorH3<R>& b)
{
 return VectorH3<R>(a.hy()*b.hz() - a.hz()*b.hy(),
                    a.hz()*b.hx() - a.hx()*b.hz(),
                    a.hx()*b.hy() - a.hy()*b.hx(),
                    a.hw()*b.hw() );
}

template <class R>
inline
typename R::Point_3
operator+(const Origin& , const VectorH3<R>& v)
{ return PointH3<R>( v ); }

template <class R>
inline
typename R::Point_3
operator-(const Origin& , const VectorH3<R>& v)
{ return  PointH3<R>(-v ); }

template <class R>
inline
typename R::Vector_3
operator-(const PointH3<R>& p, const Origin& )
{ return VectorH3<R>( p ); }

template <class R>
inline
typename R::Vector_3
operator-(const Origin& , const PointH3<R>& p)
{ return  - VectorH3<R>( p ); }

template <class R>
CGAL_KERNEL_INLINE
typename R::Point_3
operator+(const PointH3<R>& p, const VectorH3<R>& v)
{
  return PointH3<R>(p.hx()*v.hw() + v.hx()*p.hw(),
                    p.hy()*v.hw() + v.hy()*p.hw(),
                    p.hz()*v.hw() + v.hz()*p.hw(),
                    p.hw()*v.hw() );
}

template <class R>
CGAL_KERNEL_INLINE
typename R::Point_3
operator-(const PointH3<R>& p, const VectorH3<R>& v)
{
  return PointH3<R>( p.hx()*v.hw() - v.hx()*p.hw(),
                     p.hy()*v.hw() - v.hy()*p.hw(),
                     p.hz()*v.hw() - v.hz()*p.hw(),
                     p.hw()*v.hw() );
}

template <class R>
CGAL_KERNEL_INLINE
typename R::Vector_3
operator-(const PointH3<R>& p, const PointH3<R>& q)
{
  return VectorH3<R>(p.hx()*q.hw() - q.hx()*p.hw(),
                     p.hy()*q.hw() - q.hy()*p.hw(),
                     p.hz()*q.hw() - q.hz()*p.hw(),
                     p.hw()*q.hw() );
}

template < class R >
inline
typename VectorH3<R>::Vector_3
VectorH3<R>::
transform(const typename VectorH3<R>::Aff_transformation_3&t ) const
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

CGAL_END_NAMESPACE

#endif // CGAL_HOMOGENEOUS_VECTOR_3_H
