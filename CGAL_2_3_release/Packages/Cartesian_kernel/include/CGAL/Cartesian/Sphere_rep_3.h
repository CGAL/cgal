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
// file          : include/CGAL/Cartesian/Sphere_rep_3.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_SPHERE_REP_3_H
#define CGAL_CARTESIAN_SPHERE_REP_3_H

CGAL_BEGIN_NAMESPACE

template < class R >
class Simple_Sphere_repC3
{
public:
  typedef typename R::FT                        FT;
#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
  typedef typename R::Point_3                   Point_3;
#else
  typedef typename R::Point_3_base              Point_3;
#endif

  Simple_Sphere_repC3() {}

  Simple_Sphere_repC3(const Point_3 & c, const FT & r, const Orientation &o)
    : center(c), squared_radius(r), orient(o) {}

  Point_3      center;
  FT           squared_radius;
  Orientation  orient;
};

template < class R >
class Sphere_repC3 : public Ref_counted
{
public:
  typedef typename R::FT                        FT;
#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
  typedef typename R::Point_3                   Point_3;
#else
  typedef typename R::Point_3_base              Point_3;
#endif

  Sphere_repC3() {}

  Sphere_repC3(const Point_3 & c, const FT & r, const Orientation &o)
    : center(c), squared_radius(r), orient(o) {}

  Point_3      center;
  FT           squared_radius;
  Orientation  orient;
};

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_SPHERE_REP_3_H
