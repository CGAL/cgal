// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Cartesian/Point_2.h
// source        : include/CGAL/Cartesian/Point_2.h
// revision      : $2.14$
// revision_date : $Date$
// author(s)     : Andreas.Fabri@sophia.inria.fr
//                 Herve.Bronnimann@sophia.inria.fr
//
// coordinator   : INRIA Sophia-Antipolis (Herve.Bronnimann@sophia.inria.fr)
//
// ============================================================================


#ifndef CGAL_CARTESIAN_POINT_2_H
#define CGAL_CARTESIAN_POINT_2_H

#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_2_H
#include <CGAL/Cartesian/redefine_names_2.h>
#endif

#ifndef CGAL_TWOTUPLE_H
#include <CGAL/Twotuple.h>
#endif // CGAL_TWOTUPLE_H

CGAL_BEGIN_NAMESPACE

// THIS FILE IS AN EXAMPLE THAT RUNS WITH VC++
// ALL THE FUNCTION DEFINITIONS ARE PUT INSIDE THE CLASS

template < class _R >
class PointC2
#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
// This is a partial specialization
<_R,Cartesian_tag>
#endif
 : public Handle // Later we will use the handles of Lutz and Michael
{
public:
  typedef _R                                    R;
  typedef typename R::FT                        FT;
  typedef typename R::RT                        RT;
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

  PointC2() { PTR = new _Twotuple<FT>; }
  PointC2(const Origin &) { PTR = new _Twotuple<FT>(FT(0), FT(0)); }
  PointC2(const Self &p) : Handle((Handle&)p) {}
  PointC2(const FT &x, const FT &y) { PTR = new _Twotuple<FT>(x, y); }
  PointC2(const FT &hx, const FT &hy, const FT &hw) { PTR = new _Twotuple<FT>(hx/hw, hy/hw); }
  PointC2(const Vector_2 &v) : Handle((Handle&)v) {}
  ~PointC2() {}

  Self    &operator=(const Self &p)
{
  Handle::operator=(p);
  return *this;
}

  bool    operator==(const Self &p) const { return ((x() == p.x())
                                                 && (y() == p.y())) ; }
  bool    operator!=(const Self &p) const { return !(*this == p); }
  int     id() const { return (int)PTR; }

  FT      x() const { return ptr()->e0; }
  FT      y() const { return  ptr()->e1 ; }
  FT      cartesian(int i) const
{
  CGAL_kernel_precondition( (i == 0) || (i == 1) );
  return (i == 0) ? x() : y();
}
  FT      operator[](int i) const
{
  return cartesian(i);
}

  FT      hx() const { return ptr()->e0; }
  FT      hy() const { return ptr()->e1; }
  FT      hw() const { return FT(1); }
  FT      homogeneous(int i) const
{
  CGAL_kernel_precondition( (i>=0) && (i<=2) );
  if (i<2) return cartesian(i);
  return FT(1);
}

  int     dimension() const { return 2; }
  Bbox_2  bbox() const
{
  double bx = CGAL::to_double(x());
  double by = CGAL::to_double(y());
  return Bbox_2(bx,by, bx,by);
}


  Self transform(const Aff_transformation_2 &) const
{ return t.transform(*this); }


private:
  _Twotuple<FT>*  ptr() const
{
  return (_Twotuple<FT>*)PTR;
}

};

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_POINT_2_H
