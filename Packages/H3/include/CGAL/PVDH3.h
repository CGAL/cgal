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

#ifndef CGAL_HOMOGENEOUS_CLASSES_H
#include <CGAL/homogeneous_classes.h>
#endif // CGAL_HOMOGENEOUS_CLASSES_H
#ifndef CGAL_ORIGIN_H
#include <CGAL/Origin.h>
#endif // CGAL_ORIGIN_H

#ifndef CGAL_REPH3_H
#include <CGAL/RepH3.h>
#endif // CGAL_REPH3_H

CGAL_BEGIN_NAMESPACE

template <class FT, class RT>
inline
PointH3<FT,RT>
operator+  (const Origin &, const VectorH3<FT,RT> & v);

template <class FT, class RT>
inline
PointH3<FT,RT>
operator-  (const Origin &, const VectorH3<FT,RT> & v);

template <class FT, class RT>
inline
VectorH3<FT,RT>
operator-   ( const PointH3<FT,RT> &, const Origin & );

template <class FT, class RT>
inline
VectorH3<FT,RT>
operator-   ( const Origin &, const PointH3<FT,RT> & );

template <class FT, class RT>
inline
PointH3<FT,RT>
operator+   ( const Origin &, const VectorH3<FT,RT> & );

template <class FT, class RT>
inline
PointH3<FT,RT>
operator-   ( const Origin &, const VectorH3<FT,RT> & );

template <class FT, class RT>
CGAL_KERNEL_INLINE
VectorH3<FT,RT>
operator+   ( const VectorH3<FT,RT> &, const VectorH3<FT,RT> & );

template <class FT, class RT>
CGAL_KERNEL_INLINE
VectorH3<FT,RT>
operator-   ( const VectorH3<FT,RT> &, const VectorH3<FT,RT> & );

template <class FT, class RT>
CGAL_KERNEL_INLINE
FT
operator*   ( const VectorH3<FT,RT> &, const VectorH3<FT,RT> & );

template <class FT, class RT>
CGAL_KERNEL_INLINE
VectorH3<FT,RT>
operator*   ( const VectorH3<FT,RT> &, const RT & );

template <class FT, class RT>
CGAL_KERNEL_INLINE
VectorH3<FT,RT>
operator*   ( const RT &, const VectorH3<FT,RT> & );

template <class FT, class RT>
CGAL_KERNEL_INLINE
VectorH3<FT,RT>
operator/   ( const VectorH3<FT,RT> &, const RT & );

template <class FT, class RT>
CGAL_KERNEL_INLINE
VectorH3<FT,RT>
cross_product( const VectorH3<FT,RT>& a, const VectorH3<FT,RT>& b);

template <class FT, class RT>
CGAL_KERNEL_INLINE
DirectionH3<FT,RT>
cross_product( const DirectionH3<FT,RT>& a, const DirectionH3<FT,RT>& b);



template < class FT, class RT >
class PointH3 : public Handle_for< RepH3<RT> >
{

public:
  PointH3();
  PointH3(const PointH3<FT,RT> & tbc);
  PointH3(const Origin &);
  PointH3(const VectorH3<FT,RT>& v);
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

  DirectionH3<FT,RT>
        direction() const;
  PointH3<FT,RT>
        transform( const Aff_transformationH3<FT,RT> & t) const;
  Bbox_3
        bbox() const;

  bool  operator==( const PointH3<FT,RT>& p) const;
  bool  operator!=( const PointH3<FT,RT>& p) const;

  const RT&     hx_ref() const;
  const RT&     hy_ref() const;
  const RT&     hz_ref() const;
  const RT&     hw_ref() const;

friend CGAL_FRIEND_INLINE
       PointH3<FT,RT>
       operator+  CGAL_NULL_TMPL_ARGS (const Origin &,
                                       const VectorH3<FT,RT> & v);
friend CGAL_FRIEND_INLINE
       PointH3<FT,RT>
       operator-  CGAL_NULL_TMPL_ARGS (const Origin &,
                                       const VectorH3<FT,RT> & v);
};


template < class FT, class RT >
class VectorH3 : public Handle_for< RepH3<RT> >
{
public:
  VectorH3();
  VectorH3(const VectorH3<FT,RT> & tbc);
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

  DirectionH3<FT,RT>
        direction() const;
  VectorH3<FT,RT>
        transform(const Aff_transformationH3<FT,RT>& t ) const;

  VectorH3<FT,RT>
        operator-() const;

  bool  operator==( const VectorH3<FT,RT>& v) const;
  bool  operator!=( const VectorH3<FT,RT>& v) const;

// undocumented:

  VectorH3(const PointH3<FT,RT> & p);
  VectorH3(const DirectionH3<FT,RT> & dir);   /* XXX */

// friends:

/*
friend CGAL_FRIEND_INLINE
       VectorH3<FT,RT>
       operator-  CGAL_NULL_TMPL_ARGS ( const PointH3<FT,RT> &,
                                        const Origin & );
friend CGAL_FRIEND_INLINE
       VectorH3<FT,RT>
       operator-  CGAL_NULL_TMPL_ARGS ( const Origin &,
                                        const PointH3<FT,RT> & );
friend CGAL_FRIEND_INLINE
       PointH3<FT,RT>
       operator-  CGAL_NULL_TMPL_ARGS ( const Origin &,
                                        const VectorH3<FT,RT> & );
friend CGAL_KERNEL_FRIEND_INLINE
       VectorH3<FT,RT>
       operator-  CGAL_NULL_TMPL_ARGS ( const VectorH3<FT,RT> &,
                                        const VectorH3<FT,RT> & );
*/
friend CGAL_FRIEND_INLINE
       PointH3<FT,RT>
       operator+  CGAL_NULL_TMPL_ARGS ( const Origin &,
                                        const VectorH3<FT,RT> & );
friend CGAL_KERNEL_FRIEND_INLINE
       VectorH3<FT,RT>
       operator+  CGAL_NULL_TMPL_ARGS ( const VectorH3<FT,RT> &,
                                        const VectorH3<FT,RT> & );
friend CGAL_KERNEL_FRIEND_INLINE
       FT
       operator*  CGAL_NULL_TMPL_ARGS ( const VectorH3<FT,RT> &,
                                        const VectorH3<FT,RT> & );
friend CGAL_KERNEL_FRIEND_INLINE
       VectorH3<FT,RT>
       operator*  CGAL_NULL_TMPL_ARGS ( const VectorH3<FT,RT> &,
                                        const RT & );
friend CGAL_KERNEL_FRIEND_INLINE
       VectorH3<FT,RT>
       operator*  CGAL_NULL_TMPL_ARGS ( const RT &,
                                        const VectorH3<FT,RT> & );
friend CGAL_KERNEL_FRIEND_INLINE
       VectorH3<FT,RT>
       operator/  CGAL_NULL_TMPL_ARGS ( const VectorH3<FT,RT> &,
                                        const RT & );
};


template < class FT, class RT >
class DirectionH3 : public Handle_for< RepH3<RT> >
{
public:
  DirectionH3();
  DirectionH3(const DirectionH3<FT,RT>& tbc );
  DirectionH3(const PointH3<FT,RT> & p );
  DirectionH3(const VectorH3<FT,RT> & v );
  DirectionH3(const RT& x, const RT& y,
                   const RT& z, const RT& w = RT(1) );

  DirectionH3<FT,RT>
        transform(const Aff_transformationH3<FT,RT> &) const ;
  DirectionH3<FT,RT>
        operator-() const;

  bool  is_degenerate() const;

  bool  operator==( const DirectionH3<FT,RT>& d) const;
  bool  operator!=( const DirectionH3<FT,RT>& d) const;

  VectorH3<FT,RT>    to_vector() const;
  // VectorH3<FT,RT>    vector() const { return to_vector(); }

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
DirectionH3<FT,RT>
cross_product CGAL_NULL_TMPL_ARGS (const DirectionH3<FT,RT>& d1,
                                   const DirectionH3<FT,RT>& d2);
*/
};



template < class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
PointH3<FT,RT>::PointH3()
 : Handle_for< RepH3<RT> >( RepH3<RT>())
{}

template < class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
PointH3<FT,RT>::PointH3(const Origin&)
{
 const RT RT0(0);
 const RT RT1(1);
 initialize_with( RepH3<RT>( RT0, RT0, RT0, RT1 ));
}

template < class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
PointH3<FT,RT>::PointH3(const PointH3<FT,RT>& tbc)
 : Handle_for< RepH3<RT> >(tbc)
{}

template < class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
PointH3<FT,RT>::PointH3(const VectorH3<FT,RT>& v)
 : Handle_for< RepH3<RT> >(v)
{}

template < class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
PointH3<FT,RT>::PointH3(const RT& x, const RT& y, const RT& z)
 : Handle_for< RepH3<RT> >( RepH3<RT>(x,y,z, RT(1)) )
{}

template < class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
PointH3<FT,RT>::PointH3(const RT& x, const RT& y, const RT& z,
                        const RT& w)
{
  if ( w < RT(0) )
  { initialize_with( RepH3<RT>(-x,-y,-z,-w)); }
  else
  { initialize_with( RepH3<RT>(x,y,z,w)); }
}

template < class FT, class RT >
CGAL_KERNEL_INLINE
FT
PointH3<FT,RT>::x()  const
{ return ( FT(ptr->hx() ) / FT(ptr->hw() )); }


template < class FT, class RT >
CGAL_KERNEL_INLINE
FT
PointH3<FT,RT>::y()  const
{ return ( FT(ptr->hy() ) / FT(ptr->hw() )); }

template < class FT, class RT >
CGAL_KERNEL_INLINE
FT
PointH3<FT,RT>::z()  const
{ return ( FT(ptr->hz() ) / FT(ptr->hw() )); }

template < class FT, class RT >
inline
RT
PointH3<FT,RT>::hx() const
{ return  ptr->hx() ; }

template < class FT, class RT >
inline
RT
PointH3<FT,RT>::hy() const
{ return  ptr->hy() ; }

template < class FT, class RT >
inline
RT
PointH3<FT,RT>::hz() const
{ return  ptr->hz() ; }

template < class FT, class RT >
inline
RT
PointH3<FT,RT>::hw() const
{ return  ptr->hw(); }

template < class FT, class RT >
inline
const RT&
PointH3<FT,RT>::hx_ref() const
{ return  ptr->e0 ; }

template < class FT, class RT >
inline
const RT&
PointH3<FT,RT>::hy_ref() const
{ return  ptr->e1 ; }

template < class FT, class RT >
inline
const RT&
PointH3<FT,RT>::hz_ref() const
{ return  ptr->e2 ; }

template < class FT, class RT >
inline
const RT&
PointH3<FT,RT>::hw_ref() const
{ return  ptr->e3; }

template < class FT, class RT >
inline
int
PointH3<FT,RT>::dimension() const
{ return 3; }

template < class FT, class RT >
CGAL_KERNEL_INLINE
FT
PointH3<FT,RT>::cartesian(int i) const
{
  switch (i)
  {
      case 0:  return x();
      case 1:  return y();
      case 2:  return z();
      default: return cartesian( i%3 );
  }
}

template < class FT, class RT >
CGAL_KERNEL_INLINE
RT
PointH3<FT,RT>::homogeneous(int i) const
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

template < class FT, class RT >
inline
FT
PointH3<FT,RT>::operator[](int i) const
{ return cartesian(i); }


template < class FT, class RT >
inline
DirectionH3<FT,RT>
PointH3<FT,RT>::direction() const
{ return DirectionH3<FT,RT>(*this); }
template < class FT, class RT >
CGAL_KERNEL_INLINE
bool
PointH3<FT,RT>::operator==( const PointH3<FT,RT> & p) const
{
  return ( (hx() * p.hw() == p.hx() * hw() )
         &&(hy() * p.hw() == p.hy() * hw() )
         &&(hz() * p.hw() == p.hz() * hw() ) );
}

template < class FT, class RT >
inline
bool
PointH3<FT,RT>::operator!=( const PointH3<FT,RT> & p) const
{ return !(*this == p); }


#ifndef NO_OSTREAM_INSERT_POINTH3
template < class FT, class RT >
std::ostream &operator<<(std::ostream &os, const PointH3<FT,RT> &p)
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
#endif // NO_OSTREAM_INSERT_POINTH3

#ifndef NO_ISTREAM_EXTRACT_POINTH3
template < class FT, class RT >
std::istream &operator>>(std::istream &is, PointH3<FT,RT> &p)
{
  RT hx, hy, hz, hw;
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
  p = PointH3<FT,RT>(hx, hy, hz, hw);
  return is;
}
#endif // NO_ISTREAM_EXTRACT_POINTH3

template < class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
VectorH3<FT,RT>::VectorH3()
 : Handle_for< RepH3<RT> >( RepH3<RT>() )
{}

template < class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
VectorH3<FT,RT>::VectorH3(const VectorH3<FT,RT>& tbc)
 : Handle_for< RepH3<RT> >(tbc)
{}


template < class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
VectorH3<FT,RT>::VectorH3(const Null_vector&)
 : Handle_for< RepH3<RT> >( RepH3<RT>(RT(0), RT(0), RT(0), RT(1)) )
{}

template < class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
VectorH3<FT,RT>::VectorH3(const PointH3<FT,RT> & p)
 : Handle_for< RepH3<RT> >(p)
{}

template < class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
VectorH3<FT,RT>::VectorH3(const DirectionH3<FT,RT> & d)
 : Handle_for< RepH3<RT> >(d)
{}

template < class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
VectorH3<FT,RT>::VectorH3(const RT& x, const RT& y, const RT& z,
                          const RT& w)
{
  if ( w >= RT(0) )
  { initialize_with( RepH3<RT>(x, y, z, w)); }
  else
  { initialize_with( RepH3<RT>(-x,-y,-z,-w)); }
}

template < class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
VectorH3<FT,RT>::VectorH3(const RT& x, const RT& y, const RT& z)
 : Handle_for< RepH3<RT> >( RepH3<RT>(x, y, z, RT(1)) )
{}
template < class FT, class RT >
CGAL_KERNEL_INLINE
FT
VectorH3<FT,RT>::x()  const
{ return FT(ptr->hx() )/FT(ptr->hw() ) ; }

template < class FT, class RT >
CGAL_KERNEL_INLINE
FT
VectorH3<FT,RT>::y()  const
{ return FT(ptr->hy() )/FT(ptr->hw() ) ; }

template < class FT, class RT >
CGAL_KERNEL_INLINE
FT
VectorH3<FT,RT>::z()  const
{ return FT(ptr->hz() )/FT(ptr->hw() ) ; }

template < class FT, class RT >
inline
RT
VectorH3<FT,RT>::hx() const
{ return  ptr->hx() ; }

template < class FT, class RT >
inline
RT
VectorH3<FT,RT>::hy() const
{ return  ptr->hy() ; }

template < class FT, class RT >
inline
RT
VectorH3<FT,RT>::hz() const
{ return  ptr->hz() ; }

template < class FT, class RT >
inline
RT
VectorH3<FT,RT>::hw() const
{ return  ptr->hw() ; }

template < class FT, class RT >
inline
int
VectorH3<FT,RT>::dimension() const
{ return 3; }

template < class FT, class RT >
CGAL_KERNEL_INLINE
FT
VectorH3<FT,RT>::cartesian(int i) const
{
  switch (i)
  {
      case 0:   return x();
      case 1:   return y();
      case 2:   return z();
      default:  return cartesian( i%3 );
  }
}

template < class FT, class RT >
CGAL_KERNEL_INLINE
RT
VectorH3<FT,RT>::homogeneous(int i) const
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
template < class FT, class RT >
inline
DirectionH3<FT,RT>
VectorH3<FT,RT>::direction() const
{ return DirectionH3<FT,RT>(*this); }

template < class FT, class RT >
CGAL_KERNEL_INLINE
bool
VectorH3<FT,RT>::operator==( const VectorH3<FT,RT>& v) const
{
 return ( (hx() * v.hw() == v.hx() * hw() )
        &&(hy() * v.hw() == v.hy() * hw() )
        &&(hz() * v.hw() == v.hz() * hw() ) );
}

template < class FT, class RT >
inline
bool
VectorH3<FT,RT>::operator!=( const VectorH3<FT,RT>& v) const
{ return !(*this == v); }

template < class FT, class RT >
inline
FT
VectorH3<FT,RT>::operator[](int i) const
{ return cartesian(i); }

template < class FT, class RT >
CGAL_KERNEL_INLINE
VectorH3<FT,RT>
VectorH3<FT,RT>::operator-() const
{ return VectorH3<FT,RT>( - hx(), - hy(), -hz(), hw() ); }


template <class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
DirectionH3<FT,RT>::DirectionH3()
 : Handle_for< RepH3<RT> >( RepH3<RT>() )
{}

template <class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
DirectionH3<FT,RT>::DirectionH3(const DirectionH3<FT,RT>& tbc )
 : Handle_for< RepH3<RT> >(tbc)
{}

template <class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
DirectionH3<FT,RT>::DirectionH3(const PointH3<FT,RT> & p )
 : Handle_for< RepH3<RT> >( p )
{}

template <class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
DirectionH3<FT,RT>::DirectionH3(const VectorH3<FT,RT> & v )
 : Handle_for< RepH3<RT> >( v )
{}

template <class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
DirectionH3<FT,RT>::DirectionH3(const RT& x, const RT& y, const RT& z,
                                const RT& w)
{
  if ( w >= RT(0) )
  { initialize_with( RepH3<RT>(x,y,z,w)); }
  else
  { initialize_with( RepH3<RT>(-x,-y,-z,-w)); }
}

template <class FT, class RT >
CGAL_KERNEL_INLINE
RT
DirectionH3<FT,RT>::delta(int i) const
{
  switch (i)
  {
      case 0:  return x();
      case 1:  return y();
      case 2:  return z();
      default: return delta( i%3 );
  }
}

template <class FT, class RT >
inline
RT
DirectionH3<FT,RT>::dx() const
{ return ptr->e0; }

template <class FT, class RT >
inline
RT
DirectionH3<FT,RT>::x() const
{ return ptr->e0; }

template <class FT, class RT >
inline
RT
DirectionH3<FT,RT>::hx() const
{ return ptr->e0; }

template <class FT, class RT >
inline
RT
DirectionH3<FT,RT>::dy() const
{ return ptr->e1; }

template <class FT, class RT >
inline
RT
DirectionH3<FT,RT>::y() const
{ return ptr->e1; }

template <class FT, class RT >
inline
RT
DirectionH3<FT,RT>::hy() const
{ return ptr->e1; }

template <class FT, class RT >
inline
RT
DirectionH3<FT,RT>::dz() const
{ return ptr->e2; }

template <class FT, class RT >
inline
RT
DirectionH3<FT,RT>::z() const
{ return ptr->e2; }

template <class FT, class RT >
inline
RT
DirectionH3<FT,RT>::hz() const
{ return ptr->e2; }
template <class FT, class RT >
CGAL_KERNEL_INLINE
bool
DirectionH3<FT,RT>::operator==( const DirectionH3<FT,RT>& d) const
{
  return ( ( ptr->hx()*d.ptr->hy() == ptr->hy()*d.ptr->hx() )
         &&( ptr->hx()*d.ptr->hz() == ptr->hz()*d.ptr->hx() )
         &&( ptr->hy()*d.ptr->hz() == ptr->hz()*d.ptr->hy() )
         &&( CGAL_NTS sign( ptr->hx() ) == CGAL_NTS sign( d.ptr->hx() ) )
         &&( CGAL_NTS sign( ptr->hy() ) == CGAL_NTS sign( d.ptr->hy() ) )
         &&( CGAL_NTS sign( ptr->hz() ) == CGAL_NTS sign( d.ptr->hz() ) ) );
}

template <class FT, class RT >
inline
bool
DirectionH3<FT,RT>::operator!=( const DirectionH3<FT,RT>& d) const
{ return !operator==(d); }

template <class FT, class RT >
CGAL_KERNEL_INLINE
bool
DirectionH3<FT,RT>::is_degenerate() const
{ return ((hx() == RT(0)) && (hy() == RT(0)) && (hz() == RT(0))); }

template <class FT, class RT >
inline
DirectionH3<FT,RT>
DirectionH3<FT,RT>::operator-() const
{ return DirectionH3<FT,RT>(- ptr->hx(),- ptr->hy(),- ptr->hz() ); }


template <class FT, class RT >
inline
VectorH3<FT,RT>
DirectionH3<FT,RT>::to_vector() const
{ return VectorH3<FT,RT>(*this); }


template <class FT, class RT>
CGAL_KERNEL_INLINE
VectorH3<FT,RT>
operator+(const VectorH3<FT,RT>& u, const VectorH3<FT,RT>& v)
{
  return VectorH3<FT,RT>(u.hx()*v.hw() + v.hx()*u.hw(),
                         u.hy()*v.hw() + v.hy()*u.hw(),
                         u.hz()*v.hw() + v.hz()*u.hw(),
                         u.hw()*v.hw() );
}

template <class FT, class RT>
CGAL_KERNEL_INLINE
VectorH3<FT,RT>
operator-(const VectorH3<FT,RT>& u, const VectorH3<FT,RT>& v)
{
  return VectorH3<FT,RT>(u.hx()*v.hw() - v.hx()*u.hw(),
                         u.hy()*v.hw() - v.hy()*u.hw(),
                         u.hz()*v.hw() - v.hz()*u.hw(),
                         u.hw()*v.hw() );
}

template <class FT, class RT>
CGAL_KERNEL_INLINE
FT
operator*(const VectorH3<FT,RT>& u, const VectorH3<FT,RT>& v)
{
  CGAL_kernel_assertion( u.hw() != RT(0) );
  CGAL_kernel_assertion( v.hw() != RT(0) );
  return ( FT( u.hx()*v.hx() + u.hy()*v.hy() + u.hz()*v.hz() ) /
           FT( u.hw()*v.hw() ) );
}

template <class FT, class RT>
CGAL_KERNEL_INLINE
VectorH3<FT,RT>
operator/(const VectorH3<FT,RT>& v, const RT& f)
{ return VectorH3<FT,RT>( v.hx(), v.hy(), v.hz(), v.hw()*f ); }

template <class FT, class RT>
CGAL_KERNEL_INLINE
VectorH3<FT,RT>
operator*(const VectorH3<FT,RT>& v, const RT& f)
{ return VectorH3<FT,RT>( v.hx()*f, v.hy()*f, v.hz()*f, v.hw() ); }

template <class FT, class RT>
CGAL_KERNEL_INLINE
VectorH3<FT,RT>
operator*(const RT& f, const VectorH3<FT,RT>& v)
{ return VectorH3<FT,RT>( v.hx()*f, v.hy()*f, v.hz()*f, v.hw() ); }

template <class FT, class RT>
CGAL_KERNEL_INLINE
VectorH3<FT,RT>
cross_product(const VectorH3<FT,RT>& a, const VectorH3<FT,RT>& b)
{
 return VectorH3<FT,RT>(a.hy()*b.hz() - a.hz()*b.hy(),
                        a.hz()*b.hx() - a.hx()*b.hz(),
                        a.hx()*b.hy() - a.hy()*b.hx(),
                        a.hw()*b.hw() );
}

template <class FT, class RT>
inline
PointH3<FT,RT>
operator+(const Origin& , const VectorH3<FT,RT>& v)
{ return PointH3<FT,RT>( v ); }

template <class FT, class RT>
inline
PointH3<FT,RT>
operator-(const Origin& , const VectorH3<FT,RT>& v)
{ return  PointH3<FT,RT>(-v ); }

template <class FT, class RT>
inline
VectorH3<FT,RT>
operator-(const PointH3<FT,RT>& p, const Origin& )
{ return VectorH3<FT,RT>( p ); }

template <class FT, class RT>
inline
VectorH3<FT,RT>
operator-(const Origin& , const PointH3<FT,RT>& p)
{ return  - VectorH3<FT,RT>( p ); }


template <class FT, class RT>
CGAL_KERNEL_INLINE
PointH3<FT,RT>
operator+(const PointH3<FT,RT>& p, const VectorH3<FT,RT>& v)
{
  return PointH3<FT,RT>(p.hx()*v.hw() + v.hx()*p.hw(),
                        p.hy()*v.hw() + v.hy()*p.hw(),
                        p.hz()*v.hw() + v.hz()*p.hw(),
                        p.hw()*v.hw() );
}

template <class FT, class RT>
CGAL_KERNEL_INLINE
PointH3<FT,RT>
operator-(const PointH3<FT,RT>& p, const VectorH3<FT,RT>& v)
{
  return PointH3<FT,RT>( p.hx()*v.hw() - v.hx()*p.hw(),
                              p.hy()*v.hw() - v.hy()*p.hw(),
                              p.hz()*v.hw() - v.hz()*p.hw(),
                              p.hw()*v.hw() );
}

template <class FT, class RT>
CGAL_KERNEL_INLINE
VectorH3<FT,RT>
operator-(const PointH3<FT,RT>& p, const PointH3<FT,RT>& q)
{
  return PointH3<FT,RT>( p.hx()*q.hw() - q.hx()*p.hw(),
                              p.hy()*q.hw() - q.hy()*p.hw(),
                              p.hz()*q.hw() - q.hz()*p.hw(),
                              p.hw()*q.hw() );
}



template <class FT, class RT>
CGAL_KERNEL_INLINE
DirectionH3<FT,RT>
cross_product( const DirectionH3<FT,RT>& d1,
               const DirectionH3<FT,RT>& d2)
{ return cross_product(d1.to_vector(),d2.to_vector()).direction(); }

CGAL_END_NAMESPACE


#ifndef CGAL_AFF_TRANSFORMATIONH3_H
#include <CGAL/Aff_transformationH3.h>
#endif // CGAL_AFF_TRANSFORMATIONH3_H

#ifndef CGAL_BBOX_3_H
#include <CGAL/Bbox_3.h>
#endif // CGAL_BBOX_3_H
#ifndef CGAL_MISC_H
#include <CGAL/misc.h>
#endif // CGAL_MISC_H

CGAL_BEGIN_NAMESPACE


template < class FT, class RT >
inline
PointH3<FT,RT>
PointH3<FT,RT>::transform(const Aff_transformationH3<FT,RT>& t) const
{ return t.transform(*this); }

template < class FT, class RT >
CGAL_KERNEL_LARGE_INLINE
Bbox_3
PointH3<FT,RT>::bbox() const
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

template < class FT, class RT >
inline
VectorH3<FT,RT>
VectorH3<FT,RT>::transform(const Aff_transformationH3<FT,RT>&t ) const
{ return t.transform(*this); }


#ifndef NO_OSTREAM_INSERT_VECTORH3
template < class FT, class RT >
std::ostream& operator<<(std::ostream& os, const VectorH3<FT,RT>& v)
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
#endif // NO_OSTREAM_INSERT_VECTORH3

#ifndef NO_ISTREAM_EXTRACT_VECTORH3
template < class FT, class RT >
std::istream& operator>>(std::istream& is, VectorH3<FT,RT>& v)
{
  RT hx, hy, hz, hw;
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
  v = VectorH3<FT,RT>(hx, hy, hz, hw);
  return is;
}
#endif // NO_ISTREAM_EXTRACT_VECTORH3

template <class FT, class RT >
inline
DirectionH3<FT,RT>
DirectionH3<FT,RT>::
transform(const Aff_transformationH3<FT,RT>& t) const
{ return t.transform(*this); }


#ifndef NO_OSTREAM_INSERT_DIRECTIONH3
template < class FT, class RT >
std::ostream &operator<<(std::ostream &os, const DirectionH3<FT,RT> &p)
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
#endif // NO_OSTREAM_INSERT_DIRECTIONH3

#ifndef NO_ISTREAM_EXTRACT_DIRECTIONH3
template < class FT, class RT >
std::istream &operator>>(std::istream &is, DirectionH3<FT,RT> &p)
{
  RT x, y, z;
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
  p = DirectionH3<FT,RT>(x, y, z);
  return is;
}
#endif // NO_ISTREAM_EXTRACT_DIRECTIONH3


/*
template <class FT, class RT>
VectorH3<FT,RT>
operator/(const VectorH3<FT,RT>& v, const RT& f);
*/
CGAL_END_NAMESPACE


#endif // CGAL_PVDH3_H
