// ======================================================================
//
// Copyright (c) 2000 The CGAL Consortium
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
// file          : include/CGAL/Cartesian/Direction_3.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_DIRECTION_3_H
#define CGAL_CARTESIAN_DIRECTION_3_H

#include <CGAL/Cartesian/redefine_names_3.h>
#include <CGAL/Threetuple.h>

CGAL_BEGIN_NAMESPACE

template < class R_ >
class DirectionC3 CGAL_ADVANCED_KERNEL_PARTIAL_SPEC
  : public Handle_for< Threetuple<typename R_::FT> >
{
public:
  typedef R_                               R;
  typedef typename R::FT                   FT;
  typedef typename R::RT                   RT;
#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
  typedef DirectionC3<R CGAL_CTAG>         Self;
  typedef typename R::Vector_3             Vector_3;
  typedef typename R::Aff_transformation_3 Aff_transformation_3;
#else
  typedef DirectionC3<R>                   Self;
  typedef typename R::Vector_3_base        Vector_3;
  typedef typename R::Aff_transformation_3_base Aff_transformation_3;
#endif

  DirectionC3();
  DirectionC3(const Self &d);
  DirectionC3(const Vector_3 &v);
  DirectionC3(const FT &x, const FT &y, const FT &z);
  ~DirectionC3() {}

  bool           operator==(const Self &d) const;
  bool           operator!=(const Self &d) const;

  Vector_3       to_vector() const;
  Self           transform(const Aff_transformation_3 &t) const;
  Self           operator-() const;

  FT delta(int i) const;
  FT dx() const
  {
      return ptr->e0;
  }
  FT dy() const
  {
      return ptr->e1;
  }
  FT dz() const
  {
      return ptr->e2;
  }

  FT hdx() const
  {
      return dx();
  }
  FT hdy() const
  {
      return dy();
  }
  FT hdz() const
  {
      return dz();
  }
  FT hw() const
  {
      return FT(1);
  }
};

template < class R >
CGAL_KERNEL_CTOR_INLINE
DirectionC3<R CGAL_CTAG>::
DirectionC3()
{
   new ( static_cast< void*>(ptr)) Threetuple<FT>();
}

template < class R >
CGAL_KERNEL_CTOR_INLINE
DirectionC3<R CGAL_CTAG>::
DirectionC3(const DirectionC3<R CGAL_CTAG> &d)
  : Handle_for<Threetuple<typename R::FT> >(d)
{}

template < class R >
CGAL_KERNEL_CTOR_INLINE
DirectionC3<R CGAL_CTAG>::
DirectionC3(const typename DirectionC3<R CGAL_CTAG>::Vector_3 &v)
  : Handle_for<Threetuple<typename R::FT> >(v)
{}

template < class R >
CGAL_KERNEL_CTOR_INLINE
DirectionC3<R CGAL_CTAG>::
DirectionC3(const typename DirectionC3<R CGAL_CTAG>::FT &x,
            const typename DirectionC3<R CGAL_CTAG>::FT &y,
            const typename DirectionC3<R CGAL_CTAG>::FT &z)
{
  new ( static_cast< void*>(ptr)) Threetuple<FT>(x, y, z);
}

template < class R >
inline
bool
DirectionC3<R CGAL_CTAG>::operator==(const DirectionC3<R CGAL_CTAG> &d) const
{
  if (ptr == d.ptr) return true;
  return equal_directionC3(dx(), dy(), dz(), d.dx(), d.dy(), d.dz());
}

template < class R >
inline
bool
DirectionC3<R CGAL_CTAG>::operator!=(const DirectionC3<R CGAL_CTAG> &d) const
{
  return !(*this == d);
}

template < class R >
inline
typename DirectionC3<R CGAL_CTAG>::Vector_3
DirectionC3<R CGAL_CTAG>::to_vector() const
{
  return Vector_3(*this);
}

template < class R >
inline
DirectionC3<R CGAL_CTAG>
DirectionC3<R CGAL_CTAG>::
transform
  (const typename DirectionC3<R CGAL_CTAG>::Aff_transformation_3 &t) const
{
  return t.transform(*this);
}

template < class R >
inline
DirectionC3<R CGAL_CTAG> 
DirectionC3<R CGAL_CTAG>::operator-() const
{
  return DirectionC3<R>(-dx(), -dy(), -dz());
}

template < class R >
typename DirectionC3<R CGAL_CTAG>::FT 
DirectionC3<R CGAL_CTAG>::delta(int i) const
{
  CGAL_kernel_precondition( i >= 0 && i <= 2 );
  if (i==0) return dx();
  if (i==1) return dy();
  return dz();
}

#ifndef CGAL_NO_OSTREAM_INSERT_DIRECTIONC3
template < class R >
std::ostream &
operator<<(std::ostream &os, const DirectionC3<R CGAL_CTAG> &d)
{
  typename DirectionC3<R CGAL_CTAG>::Vector_3 v = d.to_vector();
  switch(os.iword(IO::mode)) {
    case IO::ASCII :
      return os << v.x() << ' ' << v.y()  << ' ' << v.z();
    case IO::BINARY :
      write(os, v.x());
      write(os, v.y());
      write(os, v.z());
      return os;
    default:
      os << "DirectionC3(" << v.x() << ", " << v.y() << ", " << v.z() << ")";
      return os;
  }
}
#endif // CGAL_NO_OSTREAM_INSERT_DIRECTIONC3

#ifndef CGAL_NO_ISTREAM_EXTRACT_DIRECTIONC3
template < class R >
std::istream &
operator>>(std::istream &is, DirectionC3<R CGAL_CTAG> &p)
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
  p = DirectionC3<R CGAL_CTAG>(x, y, z);
  return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_DIRECTIONC3

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_DIRECTION_3_H
