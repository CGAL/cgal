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
// file          : include/CGAL/Cartesian/Point_2.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri, Herve Bronnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_POINT_2_H
#define CGAL_CARTESIAN_POINT_2_H

#include <CGAL/Cartesian/redefine_names_2.h>
#include <CGAL/Origin.h>
#include <CGAL/Bbox_2.h>

CGAL_BEGIN_NAMESPACE

template < class R_ >
class PointC2 CGAL_ADVANCED_KERNEL_PARTIAL_SPEC
  : public R_::Point_handle_2
{
public:
  typedef R_                                    R;
  typedef typename R::FT                        FT;
  typedef typename R::RT                        RT;

  // Guess why we have the trailing underscore ?  Yes, that's it : VC++ !
  typedef typename R::Point_handle_2		Point_handle_2_;
  typedef typename Point_handle_2_::element_type	Point_ref_2;

#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
  typedef PointC2<R,Cartesian_tag>              Self;
  typedef typename R::Vector_2                  Vector_2;
  typedef typename R::Direction_2               Direction_2;
  typedef typename R::Line_2                    Line_2;
  typedef typename R::Ray_2                     Ray_2;
  typedef typename R::Triangle_2                Triangle_2;
  typedef typename R::Segment_2                 Segment_2;
  typedef typename R::Iso_rectangle_2           Iso_rectangle_2;
  typedef typename R::Aff_transformation_2      Aff_transformation_2;
  typedef typename R::Circle_2                  Circle_2;
#else
  typedef PointC2<R>                            Self;
  typedef typename R::Vector_2_base             Vector_2;
  typedef typename R::Direction_2_base          Direction_2;
  typedef typename R::Line_2_base               Line_2;
  typedef typename R::Ray_2_base                Ray_2;
  typedef typename R::Triangle_2_base           Triangle_2;
  typedef typename R::Segment_2_base            Segment_2;
  typedef typename R::Iso_rectangle_2_base      Iso_rectangle_2;
  typedef typename R::Aff_transformation_2_base Aff_transformation_2;
  typedef typename R::Circle_2_base             Circle_2;
#endif

  PointC2()
    : Point_handle_2_(Point_ref_2()) {}

  PointC2(const Origin &)
    : Point_handle_2_(Point_ref_2(FT(0), FT(0))) {}

  PointC2(const FT &x, const FT &y)
    : Point_handle_2_(Point_ref_2(x, y)) {}

  PointC2(const FT &hx, const FT &hy, const FT &hw)
  {
    if (hw != FT(1))
      initialize_with( Point_ref_2(hx/hw, hy/hw) );
    else
      initialize_with( Point_ref_2(hx, hy) );
  }

  PointC2(const Vector_2 &v)
    : Point_handle_2_(v) {}

  FT x() const
  {
      return Ptr()->e0;
  }
  FT y() const
  {
      return Ptr()->e1;
  }

  FT hx() const
  {
      return x();
  }
  FT hy() const
  {
      return y();
  }
  FT hw() const
  {
      return FT(1);
  }

  FT cartesian(int i) const;
  FT homogeneous(int i) const;
  FT operator[](int i) const
  {
      return cartesian(i);
  }

  int dimension() const
  {
      return 2;
  }

  bool operator==(const Self &p) const
  {
      if (identical(p))
	  return true;
      return equal_xy(*this, p);
  }
  bool operator!=(const Self &p) const
  {
      return !(*this == p);
  }

  Bbox_2 bbox() const;

  Self transform(const Aff_transformation_2 &t) const
  {
    return t.transform(*this);
  }
};

#ifdef CGAL_CFG_TYPENAME_BUG
#define typename
#endif

template < class R >
CGAL_KERNEL_INLINE
typename PointC2<R CGAL_CTAG>::FT
PointC2<R CGAL_CTAG>::cartesian(int i) const
{
  CGAL_kernel_precondition( (i == 0) || (i == 1) );
  return (i == 0) ? x() : y();
}

template < class R >
CGAL_KERNEL_INLINE
typename PointC2<R CGAL_CTAG>::FT
PointC2<R CGAL_CTAG>::homogeneous(int i) const
{
  CGAL_kernel_precondition( (i>=0) && (i<=2) );
  if (i<2)
    return cartesian(i);
  return FT(1);
}

template < class R >
CGAL_KERNEL_INLINE
Bbox_2
PointC2<R CGAL_CTAG>::bbox() const
{
  double bx = CGAL::to_double(x());
  double by = CGAL::to_double(y());
  return Bbox_2(bx,by, bx,by);
}

#ifndef CGAL_NO_OSTREAM_INSERT_POINTC2
template < class R >
std::ostream &
operator<<(std::ostream &os, const PointC2<R CGAL_CTAG> &p)
{
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << p.x() << ' ' << p.y();
    case IO::BINARY :
        write(os, p.x());
        write(os, p.y());
        return os;
    default:
        return os << "PointC2(" << p.x() << ", " << p.y() << ')';
    }
}
#endif // CGAL_NO_OSTREAM_INSERT_POINTC2

#ifndef CGAL_NO_ISTREAM_EXTRACT_POINTC2
template < class R >
std::istream &
operator>>(std::istream &is, PointC2<R CGAL_CTAG> &p)
{
    typename PointC2<R CGAL_CTAG>::FT x, y;
    switch(is.iword(IO::mode)) {
    case IO::ASCII :
        is >> x >> y;
        break;
    case IO::BINARY :
        read(is, x);
        read(is, y);
        break;
    default:
        std::cerr << "" << std::endl;
        std::cerr << "Stream must be in ascii or binary mode" << std::endl;
        break;
    }
    if (is)
	p = PointC2<R CGAL_CTAG>(x, y);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_POINTC2

#ifdef CGAL_CFG_TYPENAME_BUG
#undef typename
#endif

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_POINT_2_H
