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
#include <CGAL/Twotuple.h>

CGAL_BEGIN_NAMESPACE

template < class R_ >
class PointC2 CGAL_ADVANCED_KERNEL_PARTIAL_SPEC
 : public Handle_for< Twotuple<typename R_::FT> >
{
public:
  typedef R_                                    R;
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

  PointC2();
  PointC2(const Origin &);
  PointC2(const Self &p);
  PointC2(const FT &x, const FT &y);
  PointC2(const FT &hx, const FT &hy, const FT &hw);
  PointC2(const Vector_2 &v);
  ~PointC2();

  bool    operator==(const Self &p) const;
  bool    operator!=(const Self &p) const;
 
  FT      x() const;
  FT      y() const;
  FT      cartesian(int i) const;
  FT      operator[](int i) const;

  FT      hx() const;
  FT      hy() const;
  FT      hw() const;
  FT      homogeneous(int i) const;

  int     dimension() const;
  Bbox_2  bbox() const;


  Self transform(const Aff_transformation_2 &) const;

};

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_POINT_2_H
