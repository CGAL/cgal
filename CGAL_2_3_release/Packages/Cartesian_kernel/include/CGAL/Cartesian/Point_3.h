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

#include <CGAL/Origin.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Cartesian/redefine_names_3.h>

CGAL_BEGIN_NAMESPACE

template < class R_ >
class PointC3 CGAL_ADVANCED_KERNEL_PARTIAL_SPEC
  : public R_::Point_handle_3
{
public:
  typedef R_                               R;
  typedef typename R::FT                   FT;
  typedef typename R::RT                   RT;

  typedef typename R::Point_handle_3       Point_handle_3_;
  typedef typename Point_handle_3_::element_type Point_ref_3;

#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
  typedef PointC3<R CGAL_CTAG>             Self;
  typedef typename R::Vector_3             Vector_3;
  typedef typename R::Aff_transformation_3 Aff_transformation_3;
#else
  typedef PointC3<R>                       Self;
  typedef typename R::Vector_3_base        Vector_3;
  typedef typename R::Aff_transformation_3_base Aff_transformation_3;
#endif

  PointC3()
    : Point_handle_3_(Point_ref_3()) {}

  PointC3(const Origin &)
    : Point_handle_3_(Point_ref_3(FT(0), FT(0), FT(0))) {}

  PointC3(const Vector_3 &v)
    : Point_handle_3_(v) {}

  PointC3(const FT &x, const FT &y, const FT &z)
    : Point_handle_3_(Point_ref_3(x, y, z)) {}

  PointC3(const FT &x, const FT &y, const FT &z, const FT &w)
  {
    if (w != FT(1))
      initialize_with(Point_ref_3(x/w, y/w, z/w));
    else
      initialize_with(Point_ref_3(x, y, z));
  }

  bool operator==(const Self &p) const
  {
      if (identical(p))
	  return true;
      return x_equal(*this, p) && y_equal(*this, p) && z_equal(*this, p);
  }
  bool operator!=(const Self &p) const
  {
      return !(*this == p);
  }

  FT x() const
  {
      return Ptr()->e0;
  }
  FT y() const
  {
      return Ptr()->e1;
  }
  FT z() const
  {
      return Ptr()->e2;
  }

  FT hx() const
  {
      return x();
  }
  FT hy() const
  {
      return y();
  }
  FT hz() const
  {
      return z();
  }
  FT hw() const
  {
      return FT(1);
  }

  FT cartesian(int i) const;
  FT operator[](int i) const;
  FT homogeneous(int i) const;

  int dimension() const
  {
      return 3;
  }
  Bbox_3 bbox() const;

  Self transform(const Aff_transformation_3 &t) const
  {
    return t.transform(*this);
  }
};

#ifdef CGAL_CFG_TYPENAME_BUG
#define typename
#endif

template < class R >
inline
typename PointC3<R CGAL_CTAG>::FT
PointC3<R CGAL_CTAG>::cartesian(int i) const
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
typename PointC3<R CGAL_CTAG>::FT
PointC3<R CGAL_CTAG>::operator[](int i) const
{
  return cartesian(i);
}

template < class R >
inline
typename PointC3<R CGAL_CTAG>::FT
PointC3<R CGAL_CTAG>::homogeneous(int i) const
{
  CGAL_kernel_precondition(i>=0 && i<=3);
  if (i<3) return cartesian(i);
  return FT(1);
}

template < class R >
Bbox_3
PointC3<R CGAL_CTAG>::bbox() const
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
    if (is)
	p = PointC3<R CGAL_CTAG>(x, y, z);
    return is;
}
#endif // CGAL_CARTESIAN_NO_ISTREAM_EXTRACT_POINTC3

#ifdef CGAL_CFG_TYPENAME_BUG
#undef typename
#endif

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_POINT_3_H
