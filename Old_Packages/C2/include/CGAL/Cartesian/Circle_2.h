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
// file          : include/CGAL/Cartesian/Circle_2.h
// source        : include/CGAL/Cartesian/Circle_2.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas.Fabri@sophia.inria.fr
//                 Herve.Bronnimann@sophia.inria.fr
//
// coordinator   : INRIA Sophia-Antipolis (Herve.Bronnimann@sophia.inria.fr)
//
// ============================================================================


#ifndef CGAL_CARTESIAN_CIRCLE_2_H
#define CGAL_CARTESIAN_CIRCLE_2_H

#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_2_H
#include <CGAL/Cartesian/redefine_names_2.h>
#endif

#ifndef CGAL_CARTESIAN_CIRCLE_REP_2_H
#include <CGAL/Cartesian/Circle_rep_2.h>
#endif // CGAL_CARTESIAN_CIRCLE_REP_2_H

CGAL_BEGIN_NAMESPACE

template <class _R >
class CircleC2
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
  typedef CircleC2<R,Cartesian_tag>             Self;
  typedef typename R::Point_2                   Point_2;
  typedef typename R::Vector_2                  Vector_2;
  typedef typename R::Direction_2               Direction_2;
  typedef typename R::Line_2                    Line_2;
  typedef typename R::Ray_2                     Ray_2;
  typedef typename R::Triangle_2                Triangle_2;
  typedef typename R::Segment_2                 Segment_2;
  typedef typename R::Iso_rectangle_2           Iso_rectangle_2;
  typedef typename R::Aff_transformation_2      Aff_transformation_2;
#else
  typedef CircleC2<R>                           Self;
  typedef typename R::Point_2_base              Point_2;
  typedef typename R::Vector_2_base             Vector_2;
  typedef typename R::Direction_2_base          Direction_2;
  typedef typename R::Line_2_base               Line_2;
  typedef typename R::Ray_2_base                Ray_2;
  typedef typename R::Triangle_2_base           Triangle_2;
  typedef typename R::Segment_2_base            Segment_2;
  typedef typename R::Iso_rectangle_2_base      Iso_rectangle_2;
  typedef typename R::Aff_transformation_2_base Aff_transformation_2;
#endif

  CircleC2();
  CircleC2(const Self &);
  CircleC2(const Point_2 &center, const FT &squared_radius = FT(0),
           const Orientation &orient = COUNTERCLOCKWISE); // Is this new?
  CircleC2(const Point_2 &center, const Orientation &orient); // Is this new?
  CircleC2(const Point_2 &p, const Point_2 &q,
           const Orientation &orient = COUNTERCLOCKWISE); // And this too?
  CircleC2(const Point_2 &p, const Point_2 &q, const Point_2 &r);
  ~CircleC2();

  Self   &operator=(const Self &t);

  bool           operator==(const Self &s) const;
  bool           operator!=(const Self &s) const;
  int            id() const;

  Point_2    center() const;
  FT             squared_radius() const;

  Self   opposite() const;

//  EllipseC2<FT> transform(const Aff_transformation_2 &t) const;

  Self   orthogonal_transform(const Aff_transformation_2 &t) const;

  Orientation    orientation() const;

  Oriented_side  oriented_side(const Point_2 &p) const;
  Bounded_side   bounded_side(const Point_2 &p) const;

  bool           has_on_boundary(const Point_2 &p) const;
  bool           has_on_negative_side(const Point_2 &p) const;
  bool           has_on_positive_side(const Point_2 &p) const;

  bool           has_on_bounded_side(const Point_2 &p) const;
  bool           has_on_unbounded_side(const Point_2 &p) const;

  bool           is_degenerate() const;

  Bbox_2         bbox() const;

private:
  Circle_repC2<R> *ptr() const;
};

CGAL_END_NAMESPACE

#ifndef CGAL_CARTESIAN_CLASS_DEFINED
#ifndef CGAL_CARTESIAN_CIRCLE_2_C
#include <CGAL/Cartesian/Circle_2.C>
#endif // CGAL_CARTESIAN_CIRCLE_2_C
#endif 

#endif // CGAL_CARTESIAN_CIRCLE_2_H
