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
// file          : include/CGAL/Cartesian/Line_rep_3.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Sylvain Pion
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_LINE_REP_3_H
#define CGAL_CARTESIAN_LINE_REP_3_H

CGAL_BEGIN_NAMESPACE

template < class R >
class Simple_Line_repC3
{
public:
  typedef typename R::FT                        FT;
#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
  typedef typename R::Point_3                   Point_3;
  typedef typename R::Direction_3               Direction_3;
#else
  typedef typename R::Point_3_base              Point_3;
  typedef typename R::Direction_3_base          Direction_3;
#endif

  Simple_Line_repC3() {}

  Simple_Line_repC3(const Point_3 &p, const Direction_3 &d)
    : basepoint(p), direction(d) {}

// private:
  Point_3       basepoint;
  Direction_3   direction;
};

template < class R >
class Line_repC3 : public Ref_counted
{
public:
  typedef typename R::FT                        FT;
#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
  typedef typename R::Point_3                   Point_3;
  typedef typename R::Direction_3               Direction_3;
#else
  typedef typename R::Point_3_base              Point_3;
  typedef typename R::Direction_3_base          Direction_3;
#endif

  Line_repC3() {}

  Line_repC3(const Point_3 &p, const Direction_3 &d)
    : basepoint(p), direction(d) {}

// private:
  Point_3       basepoint;
  Direction_3   direction;
};

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_LINE_REP_3_H
