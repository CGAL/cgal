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
// file          : include/CGAL/Cartesian/Iso_rectangle_2.h
// source        : include/CGAL/Cartesian/Cartesian_2/Iso_rectangle_2.h
// revision      : include/CGAL/Cartesian/Cartesian_2/Iso_rectangle_2.h
// revision_date : $Date$
// author(s)     : Andreas.Fabri@sophia.inria.fr
//                 Herve Bronnimann@sophia.inria.fr
//
// coordinator   : INRIA Sophia-Antipolis (Herve.Bronnimann@sophia.inria.fr)
//
// ============================================================================


#ifndef CGAL_CARTESIAN_ISO_RECTANGLE_2_H
#define CGAL_CARTESIAN_ISO_RECTANGLE_2_H

#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_2_H
#include <CGAL/Cartesian/redefine_names_2.h>
#endif

#ifndef CGAL_TWOTUPLE_H
#include <CGAL/Twotuple.h>
#endif // CGAL_TWOTUPLE_H
#ifndef CGAL_CARTESIAN_POINT_2_H
#include <CGAL/Point_2.h>
#endif // CGAL_CARTESIAN_POINT_2_H

CGAL_BEGIN_NAMESPACE

template <class _R>
class Iso_rectangleC2
#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
// This is a partial specialization
<_R,Cartesian_tag>
#endif
  : public Handle
{
public:
  typedef _R                                    R;
  typedef typename R::FT                        FT;
  typedef typename R::RT                        RT;
#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
  typedef Iso_rectangleC2<R,Cartesian_tag>      Self;
  typedef typename R::Point_2                   Point_2;
  typedef typename R::Vector_2                  Vector_2;
  typedef typename R::Direction_2               Direction_2;
  typedef typename R::Line_2                    Line_2;
  typedef typename R::Ray_2                     Ray_2;
  typedef typename R::Triangle_2                Triangle_2;
  typedef typename R::Segment_2                 Segment_2;
  typedef typename R::Aff_transformation_2      Aff_transformation_2;
  typedef typename R::Circle_2                  Circle_2;
#else
  typedef Iso_rectangleC2<R>                    Self;
  typedef typename R::Point_2_base              Point_2;
  typedef typename R::Vector_2_base             Vector_2;
  typedef typename R::Direction_2_base          Direction_2;
  typedef typename R::Line_2_base               Line_2;
  typedef typename R::Ray_2_base                Ray_2;
  typedef typename R::Triangle_2_base           Triangle_2;
  typedef typename R::Segment_2_base            Segment_2;
  typedef typename R::Aff_transformation_2_base Aff_transformation_2;
  typedef typename R::Circle_2_base             Circle_2;
#endif

  Iso_rectangleC2();
  Iso_rectangleC2(const Self &);
  Iso_rectangleC2(const Point_2 &p, const Point_2 &q);
  ~Iso_rectangleC2();

  Self           &operator=(const Self &r);

  bool            operator==(const Self &s) const;
  bool            operator!=(const Self &s) const;
  int             id() const;

  Point_2         min() const;
  Point_2         max() const;
  Point_2         vertex(int i) const;
  Point_2         operator[](int i) const;

  Self            transform(const Aff_transformation_2 &t) const;

  Bounded_side    bounded_side(const Point_2 &p) const;
  bool            has_on_boundary(const Point_2 &p) const;
  bool            has_on_bounded_side(const Point_2 &p) const;
  bool            has_on_unbounded_side(const Point_2 &p) const;

  bool            is_degenerate() const;

  Bbox_2          bbox() const;

  FT              xmin() const;
  FT              ymin() const;
  FT              xmax() const;
  FT              ymax() const;

private:
  _Twotuple< Point_2 >*   ptr() const;
};

CGAL_END_NAMESPACE


#ifndef CGAL_CARTESIAN_CLASS_DEFINED
#ifndef CGAL_CARTESIAN_ISO_RECTANGLE_2_C
#include <CGAL/Cartesian/Iso_rectangle_2.C>
#endif CGAL_CARTESIAN_ISO_RECTANGLE_2_C
#endif 

#endif // CGAL_CARTESIAN_ISO_RECTANGLE_2_H
