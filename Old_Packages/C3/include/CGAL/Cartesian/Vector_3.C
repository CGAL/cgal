// ======================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.1-I-28 $
// release_date  : $CGAL_Date: 1999/10/12 $
//
// file          : include/CGAL/Cartesian/Vector_3.C
// package       : C3 (3.6.2)
// source        : include/CGAL/Cartesian/Vector_3.C
// revision      : $Revision$
// revision_date : $Date$
// author        : Andreas.Fabri@sophia.inria.fr
//
// coordinator   : INRIA Sophia-Antipolis (Herve.Bronnimann@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_3_H
#define CGAL_CTAG
#endif

#ifdef CGAL_CFG_TYPENAME_BUG
#define typename
#endif

#ifndef CGAL_CARTESIAN_VECTOR_3_C
#define CGAL_CARTESIAN_VECTOR_3_C

#ifndef CGAL_CARTESIAN_DIRECTION_3_H
#include <CGAL/Cartesian/Direction_3.h>
#endif // CGAL_CARTESIAN_DIRECTION_3_H

CGAL_BEGIN_NAMESPACE

template < class R >
inline _Threetuple<typename R::FT>*
VectorC3<R CGAL_CTAG>::ptr() const
{
  return (_Threetuple<FT>*)PTR;
}

template < class R >
VectorC3<R CGAL_CTAG>::VectorC3()
{
  PTR = new _Threetuple<FT>(FT(0), FT(0), FT(0));
}

template < class R >
VectorC3<R CGAL_CTAG>::
VectorC3(const VectorC3<R CGAL_CTAG>  &v)
  : Handle(v)
{}

template < class R >
VectorC3<R CGAL_CTAG>::
VectorC3(const Null_vector  &)
{
  PTR = new _Threetuple<FT>(FT(0), FT(0), FT(0));
}

template < class R >
VectorC3<R CGAL_CTAG>::
VectorC3(const typename VectorC3<R CGAL_CTAG>::FT &x,
         const typename VectorC3<R CGAL_CTAG>::FT &y,
         const typename VectorC3<R CGAL_CTAG>::FT &z)
{
  PTR = new _Threetuple<FT>(x, y, z);
}

template < class R >
VectorC3<R CGAL_CTAG>::
VectorC3(const typename VectorC3<R CGAL_CTAG>::FT &x,
         const typename VectorC3<R CGAL_CTAG>::FT &y,
         const typename VectorC3<R CGAL_CTAG>::FT &z,
         const typename VectorC3<R CGAL_CTAG>::FT &w)
{
  if (w != FT(1))
    PTR = new _Threetuple<FT>(x/w, y/w, z/w);
  else
    PTR = new _Threetuple<FT>(x, y, z);
}

template < class R >
VectorC3<R CGAL_CTAG>::~VectorC3()
{}

template < class R >
VectorC3<R CGAL_CTAG> &
VectorC3<R CGAL_CTAG>::operator=(const VectorC3<R CGAL_CTAG> &v)
{
  Handle::operator=(v);
  return *this;
}

template < class R >
inline
VectorC3<R CGAL_CTAG>::
VectorC3(const typename VectorC3<R CGAL_CTAG>::Point_3 &p)
  : Handle((Handle&)p)
{
}

template < class R >
inline VectorC3<R CGAL_CTAG>::
VectorC3(const typename VectorC3<R CGAL_CTAG>::Direction_3 &d) :
  Handle((Handle&)d)
{
}

template < class R >
bool VectorC3<R CGAL_CTAG>::operator==(const VectorC3<R CGAL_CTAG> &v) const
{
  return (x() == v.x()) && (y() == v.y()) && (z() == v.z()) ;
}

template < class R >
inline bool VectorC3<R CGAL_CTAG>::operator!=(const VectorC3<R CGAL_CTAG> &v) const
{
  return !(*this == v);
}


template < class R >
bool VectorC3<R CGAL_CTAG>::operator==(const Null_vector &) const
{
  return (x() == FT(0)) && (y() == FT(0)) && (z() == FT(0)) ;
}

template < class R >
inline bool VectorC3<R CGAL_CTAG>::operator!=(const Null_vector &v) const
{
  return !(*this == v);
}


template < class R >
inline
long VectorC3<R CGAL_CTAG>::id() const
{
  return (long) PTR;
}
template < class R >
inline
typename VectorC3<R CGAL_CTAG>::FT
VectorC3<R CGAL_CTAG>::x()  const
{
  return ptr()->e0;
}

template < class R >
inline
typename VectorC3<R CGAL_CTAG>::FT
VectorC3<R CGAL_CTAG>::y()  const
{
  return  ptr()->e1;
}

template < class R >
inline
typename VectorC3<R CGAL_CTAG>::FT
VectorC3<R CGAL_CTAG>::z()  const
{
  return  ptr()->e2;
}

template < class R >
inline
typename VectorC3<R CGAL_CTAG>::FT
VectorC3<R CGAL_CTAG>::cartesian(int i) const
{
  CGAL_kernel_precondition( (i>=0) && (i<3) );
  return (i==0) ? x() :
         (i==1) ? y() : z();
}

template < class R >
inline
typename VectorC3<R CGAL_CTAG>::FT
VectorC3<R CGAL_CTAG>::operator[](int i) const
{
  return cartesian(i);
}

template < class R >
inline
int
VectorC3<R CGAL_CTAG>::dimension() const
{
  return 3;
}

template < class R >
inline
typename VectorC3<R CGAL_CTAG>::FT
VectorC3<R CGAL_CTAG>::hx()  const
{
  return ptr()->e0;
}

template < class R >
inline
typename VectorC3<R CGAL_CTAG>::FT
VectorC3<R CGAL_CTAG>::hy()  const
{
  return ptr()->e1;
}

template < class R >
inline
typename VectorC3<R CGAL_CTAG>::FT
VectorC3<R CGAL_CTAG>::hz()  const
{
  return ptr()->e2;
}

template < class R >
typename VectorC3<R CGAL_CTAG>::FT
VectorC3<R CGAL_CTAG>::hw()  const
{
  return FT(1);
}

template < class R >
typename VectorC3<R CGAL_CTAG>::FT
VectorC3<R CGAL_CTAG>::homogeneous(int i) const
{
  return (i==3) ? FT(1) : cartesian(i);
}

template < class R >
inline
VectorC3<R CGAL_CTAG>
VectorC3<R CGAL_CTAG>::operator+(const VectorC3<R CGAL_CTAG> &w) const
{
  return VectorC3<R CGAL_CTAG>(x() + w.x(), y() + w.y(), z() + w.z()) ;
}

template < class R >
inline
VectorC3<R CGAL_CTAG>
VectorC3<R CGAL_CTAG>::operator-(const VectorC3<R CGAL_CTAG> &w) const
{
  return VectorC3<R CGAL_CTAG>(x() - w.x(), y() - w.y(), z() - w.z()) ;
}

template < class R >
inline
VectorC3<R CGAL_CTAG> VectorC3<R CGAL_CTAG>::operator-() const
{

  return VectorC3<R CGAL_CTAG>(-x(), -y(), -z()) ;
}

template < class R >
inline
typename VectorC3<R CGAL_CTAG>::FT
VectorC3<R CGAL_CTAG>::operator*(const VectorC3<R CGAL_CTAG> &w) const
{
  return x() * w.x() + y() * w.y() + z() * w.z() ;
}

template < class R >
inline
VectorC3<R CGAL_CTAG>
VectorC3<R CGAL_CTAG>::operator/(const typename VectorC3<R CGAL_CTAG>::FT &c) const
{
  return VectorC3<R CGAL_CTAG>( x()/c, y()/c, z()/c) ;
}

template < class R >
inline
typename VectorC3<R CGAL_CTAG>::Direction_3
VectorC3<R CGAL_CTAG>::direction() const
{
  return Direction_3(*this);
}

template < class R >
VectorC3<R CGAL_CTAG>
VectorC3<R CGAL_CTAG>::
transform(const typename VectorC3<R CGAL_CTAG>::Aff_transformation_3 &t) const
{
  return t.transform(*this);
}


#ifndef CGAL_CARTESIAN_NO_OSTREAM_INSERT_VECTORC3
template < class R >
std::ostream &operator<<(std::ostream &os, const VectorC3<R CGAL_CTAG> &v)
{
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << v.x() << ' ' << v.y()  << ' ' << v.z();
    case IO::BINARY :
        write(os, v.x());
        write(os, v.y());
        write(os, v.z());
        return os;
    default:
        os << "VectorC3(" << v.x() << ", " << v.y() <<  ", " << v.z() << ")";
        return os;
    }
}
#endif // CGAL_CARTESIAN_NO_OSTREAM_INSERT_VECTORC3

#ifndef CGAL_CARTESIAN_NO_ISTREAM_EXTRACT_VECTORC3
template < class R >
std::istream &operator>>(std::istream &is, VectorC3<R CGAL_CTAG> &p)
{
    typename R::FT x, y, z;
    switch(is.iword(IO::mode)) {
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
    p = VectorC3<R CGAL_CTAG>(x, y, z);
    return is;
}
#endif // CGAL_CARTESIAN_NO_ISTREAM_EXTRACT_VECTORC3

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_TYPENAME_BUG
#undef typename
#endif

#endif // CGAL_CARTESIAN_VECTOR_3_C
