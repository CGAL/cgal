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

  PointC2();
  PointC2(const Origin &);
  PointC2(const Self &p);
  PointC2(const FT &x, const FT &y);
  PointC2(const FT &hx, const FT &hy, const FT &hw);
  PointC2(const Vector_2 &v);
  ~PointC2();

  Self    &operator=(const Self &p);

  bool    operator==(const Self &p) const;
  bool    operator!=(const Self &p) const;
  int     id() const;

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

private:
  _Twotuple<FT>*  ptr() const;
};

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_POINT_2_H
