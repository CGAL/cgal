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
// file          : include/CGAL/Cartesian/Point_3.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri and Hervé Brönnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_POINT_3_H
#define CGAL_CARTESIAN_POINT_3_H

#include <CGAL/Threetuple.h>
#include <CGAL/Origin.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Kernel/Cartesian_coordinate_iterator_3.h>

CGAL_BEGIN_NAMESPACE

template < class R_ >
class PointC3
  : public R_::template Handle<Threetuple<typename R_::FT> >::type
{
CGAL_VC7_BUG_PROTECTED
  typedef typename R_::FT                   FT;
  typedef typename R_::Vector_3             Vector_3;
  typedef typename R_::Point_3              Point_3;
  typedef typename R_::Aff_transformation_3 Aff_transformation_3;

  typedef Threetuple<FT>                           rep;
  typedef typename R_::template Handle<rep>::type  base;

public:
  typedef Cartesian_coordinate_iterator_3<R_> Cartesian_const_iterator;
  typedef R_                                R;

  PointC3()
    : base(rep()) {}

  PointC3(const Origin &)
    : base(rep(FT(0), FT(0), FT(0))) {}

  PointC3(const Vector_3 &v)
    : base(v) {}

  PointC3(const FT &x, const FT &y, const FT &z)
    : base(rep(x, y, z)) {}

  PointC3(const FT &x, const FT &y, const FT &z, const FT &w)
  {
    if (w != FT(1))
      initialize_with(rep(x/w, y/w, z/w));
    else
      initialize_with(rep(x, y, z));
  }

  bool operator==(const PointC3 &p) const
  {
      if (identical(p))
	  return true;
      return x_equal(*this, p) && y_equal(*this, p) && z_equal(*this, p);
  }
  bool operator!=(const PointC3 &p) const
  {
      return !(*this == p);
  }

  const FT & x() const
  {
      return Ptr()->e0;
  }
  const FT & y() const
  {
      return Ptr()->e1;
  }
  const FT & z() const
  {
      return Ptr()->e2;
  }

  const FT & hx() const
  {
      return x();
  }
  const FT & hy() const
  {
      return y();
  }
  const FT & hz() const
  {
      return z();
  }
  FT hw() const
  {
      return FT(1);
  }

  const FT & cartesian(int i) const;
  const FT & operator[](int i) const;
  FT homogeneous(int i) const;

  Cartesian_const_iterator cartesian_begin() const 
  {
    return Cartesian_const_iterator(static_cast<const Point_3* >(this),0);
    //return Cartesian_const_iterator(this,0);
  }

  Cartesian_const_iterator cartesian_end() const 
  {
    return Cartesian_const_iterator(static_cast<const Point_3* >(this), 3);
    //return Cartesian_const_iterator(this, 3);
  }

  int dimension() const
  {
      return 3;
  }
  Bbox_3 bbox() const;

  PointC3 transform(const Aff_transformation_3 &t) const
  {
    return t.transform(*this);
  }
};

template < class R >
inline
const typename PointC3<R>::FT &
PointC3<R>::cartesian(int i) const
{
  CGAL_kernel_precondition( (i>=0) && (i<=2) );
  // return (i==0) ? x() :
//          (i==1) ? y() : z();
  if (i==0) return x();
  if (i==1) return y();
  return z();
}

template < class R >
inline
const typename PointC3<R>::FT &
PointC3<R>::operator[](int i) const
{
  return cartesian(i);
}

template < class R >
inline
typename PointC3<R>::FT
PointC3<R>::homogeneous(int i) const
{
  CGAL_kernel_precondition(i>=0 && i<=3);
  if (i<3) return cartesian(i);
  return FT(1);
}

template < class R >
Bbox_3
PointC3<R>::bbox() const
{
  // FIXME: to_interval
  double bx = CGAL::to_double(x());
  double by = CGAL::to_double(y());
  double bz = CGAL::to_double(z());
  return Bbox_3(bx, by, bz, bx, by, bz);
}

#ifndef CGAL_CARTESIAN_NO_OSTREAM_INSERT_POINTC3
template < class R >
std::ostream &
operator<<(std::ostream &os, const PointC3<R> &p)
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
operator>>(std::istream &is, PointC3<R> &p)
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
    if (is)
	p = PointC3<R>(x, y, z);
    return is;
}
#endif // CGAL_CARTESIAN_NO_ISTREAM_EXTRACT_POINTC3

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_POINT_3_H
