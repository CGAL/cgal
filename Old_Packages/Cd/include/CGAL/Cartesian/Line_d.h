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
// file          : include/CGAL/Cartesian/Line_d.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve.Bronnimann@sophia.inria.fr
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_LINE_D_H
#define CGAL_CARTESIAN_LINE_D_H

#include <CGAL/Cartesian/redefine_names_d.h>
#include <CGAL/Twotuple.h>

CGAL_BEGIN_NAMESPACE

template < class R_ >
class LineCd CGAL_ADVANCED_KERNEL_PARTIAL_SPEC
  : public Handle
{
public:
  typedef R_                               R;
  typedef typename R::FT                   FT;
  typedef typename R::RT                   RT;
#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
  typedef LineCd<R CGAL_CTAG>              Self;
  typedef typename R::Point_d              Point_d;
  typedef typename R::Vector_d             Vector_d;
  typedef typename R::Direction_d          Direction_d;
  typedef typename R::Plane_d              Plane_d;
  typedef typename R::Ray_d                Ray_d;
  typedef typename R::Segment_d            Segment_d;
  typedef typename R::Aff_transformation_d Aff_transformation_d;
#else
  typedef LineCd<R>                             Self;
  typedef typename R::Point_d_base              Point_d;
  typedef typename R::Vector_d_base             Vector_d;
  typedef typename R::Direction_d_base          Direction_d;
  typedef typename R::Plane_d_base              Plane_d;
  typedef typename R::Ray_d_base                Ray_d;
  typedef typename R::Segment_d_base            Segment_d;
  typedef typename R::Aff_transformation_d_base Aff_transformation_d;
#endif

  LineCd();
  LineCd(const Self  &l);
  LineCd(const Point_d &p, const Point_d &q);
  LineCd(const Segment_d &s);
  LineCd(const Ray_d &r);
  LineCd(const Point_d &p, const Direction_d &d);
  ~LineCd();

  Self        &operator=(const Self &l);

  bool        operator==(const Self &l) const;
  bool        operator!=(const Self &l) const;
  long        id() const;
  int         dimension() const;

  Plane_d     perpendicular_plane(const Point_d &p) const;
  Self        opposite() const;

  Point_d     point() const;
  Point_d     point(int i) const;

  Point_d     projection(const Point_d &p) const;

  Direction_d direction() const;

  bool        has_on(const Point_d &p) const;
  bool        is_degenerate() const;

  Self        transform(const Aff_transformation_d &t) const;

private:
  _Twotuple< Point_d >* ptr() const;
  void         new_rep(const Point_d &p,
                       const Vector_d &v);
};

CGAL_END_NAMESPACE

#ifndef CGAL_CARTESIAN_CLASS_DEFINED
#include <CGAL/Cartesian/Line_d.C>
#endif 

#endif // CGAL_CARTESIAN_LINE_D_H
