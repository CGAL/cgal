// ======================================================================
//
// Copyright (c) 1998, 1999 The CGAL Consortium
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
// file          : include/CGAL/Cartesian/Point_3.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
//
// coordinator   : INRIA Sophia-Antipolis (Herve.Bronnimann@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_POINT_3_C
#define CGAL_CARTESIAN_POINT_3_C

#include <CGAL/Origin.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/number_utils.h>

#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_3_H
#define CGAL_CTAG
#endif

#ifdef CGAL_CFG_TYPENAME_BUG
#define typename
#endif

CGAL_BEGIN_NAMESPACE

template < class R >
inline
_Threetuple<typename R::FT>*
PointC3<R CGAL_CTAG>::ptr() const
{
  return (_Threetuple<FT>*)PTR;
}

template < class R >
PointC3<R CGAL_CTAG>::PointC3()
{
  PTR = new _Threetuple<typename R::FT>(FT(0), FT(0), FT(0));
}

template < class R >
PointC3<R CGAL_CTAG>::PointC3(const Origin &)
{
  PTR = new _Threetuple<FT>(FT(0), FT(0), FT(0));
}

template < class R >
PointC3<R CGAL_CTAG>::PointC3(const PointC3<R CGAL_CTAG> &p)
  : Handle((Handle&)p)
{}

template < class R >
PointC3<R CGAL_CTAG>::PointC3(const typename PointC3<R CGAL_CTAG>::FT &x,
                              const typename PointC3<R CGAL_CTAG>::FT &y,
			      const typename PointC3<R CGAL_CTAG>::FT &z)
{
  PTR = new _Threetuple<FT>(x, y, z);
}

template < class R >
PointC3<R CGAL_CTAG>::PointC3(const typename PointC3<R CGAL_CTAG>::FT &x,
                              const typename PointC3<R CGAL_CTAG>::FT &y,
			      const typename PointC3<R CGAL_CTAG>::FT &z,
			      const FT &w)
{
  if (w != FT(1))
    PTR = new _Threetuple<FT>(x/w, y/w, z/w);
  else
    PTR = new _Threetuple<FT>(x, y, z);
}

template < class R >
PointC3<R CGAL_CTAG>::~PointC3()
{}

template < class R >
PointC3<R CGAL_CTAG> &
PointC3<R CGAL_CTAG>::operator=(const PointC3<R CGAL_CTAG> &p)
{
  Handle::operator=(p);
  return *this;
}

template < class R >
PointC3<R CGAL_CTAG>::
PointC3(const typename PointC3<R CGAL_CTAG>::Vector_3 &v)
  : Handle((Handle&)v)
{}

template < class R >
inline
bool
PointC3<R CGAL_CTAG>::operator==(const PointC3<R CGAL_CTAG>& p) const
{
  return (x() == p.x()) && (y() == p.y()) && (z() == p.z());
}

template < class R >
inline
bool
PointC3<R CGAL_CTAG>::operator!=(const PointC3<R CGAL_CTAG>& p) const
{
  return !(*this == p);
}

template < class R >
inline
long
PointC3<R CGAL_CTAG>::id() const
{
  return (long) PTR;
}

template < class R >
inline
typename PointC3<R CGAL_CTAG>::FT
PointC3<R CGAL_CTAG>::x()  const
{
  return ptr()->e0;
}

template < class R >
inline
typename PointC3<R CGAL_CTAG>::FT
PointC3<R CGAL_CTAG>::y()  const
{
  return  ptr()->e1;
}

template < class R >
inline
typename PointC3<R CGAL_CTAG>::FT
PointC3<R CGAL_CTAG>::z()  const
{
  return  ptr()->e2;
}

template < class R >
inline
int
PointC3<R CGAL_CTAG>::dimension() const
{
  return 3;
}

template < class R >
inline
typename PointC3<R CGAL_CTAG>::FT
PointC3<R CGAL_CTAG>::cartesian(int i) const
{
  CGAL_kernel_precondition( (i>=0) && (i<=2) );
  return (i==0) ? x() :
         (i==1) ? y() : z();
}

template < class R >
inline
typename PointC3<R CGAL_CTAG>::FT
PointC3<R CGAL_CTAG>::operator[](int i) const
{
  return cartesian(i);
}

template < class R >
inline
typename PointC3<R CGAL_CTAG>::FT
PointC3<R CGAL_CTAG>::hx()  const
{
  return ptr()->e0;
}

template < class R >
inline
typename PointC3<R CGAL_CTAG>::FT
PointC3<R CGAL_CTAG>::hy()  const
{
  return ptr()->e1;
}

template < class R >
inline
typename PointC3<R CGAL_CTAG>::FT
PointC3<R CGAL_CTAG>::hz()  const
{
  return ptr()->e2;
}

template < class R >
inline
typename PointC3<R CGAL_CTAG>::FT
PointC3<R CGAL_CTAG>::hw()  const
{
  return FT(1);
}

template < class R >
typename PointC3<R CGAL_CTAG>::FT
PointC3<R CGAL_CTAG>::homogeneous(int i) const
{
  CGAL_kernel_precondition((i>=0) && (i<=3));
  return (i<3) ? cartesian(i) : FT(1);
}

template < class R >
inline
PointC3<R CGAL_CTAG>
PointC3<R CGAL_CTAG>::
transform(const typename PointC3<R CGAL_CTAG>::Aff_transformation_3 &t) const
{
  return t.transform(*this);
}

template < class R >
Bbox_3
PointC3<R CGAL_CTAG>::bbox() const
{
  // Not robust...
  double bx = CGAL::to_double(x());
  double by = CGAL::to_double(y());
  double bz = CGAL::to_double(z());
  return Bbox_3(bx, by, bz, bx, by, bz);
}

#ifndef CGAL_CARTESIAN_NO_OSTREAM_INSERT_POINTC3
template < class R >
std::ostream &
operator<<(std::ostream &os, const PointC3<R CGAL_CTAG> &p)
{
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << p.x() << ' ' << p.y()  << ' ' << p.z();
    case IO::BINARY :
        write(os, p.x());
        write(os, p.y());
        write(os, p.z());
        return os;
    default:
        os << "PointC3(" << p.x() << ", " << p.y()  << ", " << p.z() << ")";
        return os;
    }
}
#endif // CGAL_CARTESIAN_NO_OSTREAM_INSERT_POINTC3

#ifndef CGAL_CARTESIAN_NO_ISTREAM_EXTRACT_POINTC3
template < class R >
std::istream &
operator>>(std::istream &is, PointC3<R CGAL_CTAG> &p)
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
    p = PointC3<R CGAL_CTAG>(x, y, z);
    return is;
}
#endif // CGAL_CARTESIAN_NO_ISTREAM_EXTRACT_POINTC3

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_TYPENAME_BUG
#undef typename
#endif

#endif // CGAL_CARTESIAN_POINT_3_C
