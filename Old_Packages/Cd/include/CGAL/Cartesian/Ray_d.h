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
// file          : include/CGAL/Cartesian/Ray_d.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve.Bronnimann@sophia.inria.fr
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_RAY_D_H
#define CGAL_CARTESIAN_RAY_D_H

#include <CGAL/Cartesian/redefine_names_d.h>
#include <CGAL/Twotuple.h>

CGAL_BEGIN_NAMESPACE

template < class R_ >
class RayCd CGAL_ADVANCED_KERNEL_PARTIAL_SPEC
  : public Handle
{
public:
  typedef R_                                    R;
  typedef typename R::FT                        FT;
  typedef typename R::RT                        RT;
#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
  typedef RayCd<R CGAL_CTAG>                    Self;
  typedef typename R::Point_d                   Point_d;
  typedef typename R::Direction_d               Direction_d;
  typedef typename R::Line_d                    Line_d;
  typedef typename R::Aff_transformation_d      Aff_transformation_d;
#else
  typedef RayCd<R>                              Self;
  typedef typename R::Point_d_base              Point_d;
  typedef typename R::Direction_d_base          Direction_d;
  typedef typename R::Line_d_base               Line_d;
  typedef typename R::Aff_transformation_d_base Aff_transformation_d;
#endif

  RayCd();
  RayCd(const Self &r);
  RayCd(const Point_d &sp, const Point_d &secondp);
  RayCd(const Point_d &sp, const Direction_d &d);
  ~RayCd();

  Self        &operator=(const Self &r);

  bool        operator==(const Self &r) const;
  bool        operator!=(const Self &r) const;
  long        id() const;

  Point_d     start() const;
  Point_d     source() const;
  Point_d     second_point() const;
  Point_d     point(int i) const;

  Direction_d direction() const;
  Line_d      supporting_line() const;
  Self        opposite() const;

  Self        transform(const Aff_transformation_d &t) const;

  bool        is_degenerate() const;
  bool        has_on(const Point_d &p) const;
  bool        collinear_has_on(const Point_d &p) const;

private:
  _Twotuple< Point_d > *ptr() const;
};

CGAL_END_NAMESPACE

#ifndef CGAL_CARTESIAN_CLASS_DEFINED
#include <CGAL/Cartesian/Ray_d.C>
#endif 

#endif // CGAL_CARTESIAN_RAY_D_H
