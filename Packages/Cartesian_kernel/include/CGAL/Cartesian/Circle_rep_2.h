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
// file          : include/CGAL/Cartesian/Circle_rep_2.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri, Herve Bronnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_CIRCLE_REP_2_H
#define CGAL_CARTESIAN_CIRCLE_REP_2_H

CGAL_BEGIN_NAMESPACE

template < class R >
class Simple_Circle_repC2
{
public:
  typedef typename R::FT                        FT;
#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
  typedef typename R::Point_2                   Point_2;
#else
  typedef typename R::Point_2_base              Point_2;
#endif

  Simple_Circle_repC2() {}

  Simple_Circle_repC2(const Point_2 & c, const FT & r, const Orientation &o)
    : center(c), squared_radius(r), orient(o) {}

  Point_2      center;
  FT           squared_radius;
  Orientation  orient;
};

template < class R >
class Circle_repC2 : public Ref_counted
{
public:
  typedef typename R::FT                        FT;
#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
  typedef typename R::Point_2                   Point_2;
#else
  typedef typename R::Point_2_base              Point_2;
#endif

  Circle_repC2() {}

  Circle_repC2(const Point_2 & c, const FT & r, const Orientation &o)
    : center(c), squared_radius(r), orient(o) {}

  Point_2      center;
  FT           squared_radius;
  Orientation  orient;
};

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_CIRCLE_REP_2_H
