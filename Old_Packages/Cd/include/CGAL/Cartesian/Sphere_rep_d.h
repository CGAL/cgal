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
// file          : include/CGAL/Cartesian/Sphere_rep_d.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve.Bronnimann@sophia.inria.fr
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_SPHERE_REP_D_H
#define CGAL_CARTESIAN_SPHERE_REP_D_H

CGAL_BEGIN_NAMESPACE

template < class R >
class Sphere_repC3 : public Rep
{
public:
  typedef typename R::FT                        FT;
#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
  typedef typename R::Point_d                   Point_d;
#else
  typedef typename R::Point_D_base              Point_d;
#endif

  Point_d      center;
  FT           squared_radius;
  Orientation  orient;


  Sphere_repC3() {}

  Sphere_repC3(const Point_d & c, const FT & r, const Orientation &o)
    : center(c), squared_radius(r), orient(o) {}

  ~Sphere_repC3() {}
};

CGAL_END_NAMESPACE

#endif  // CGAL_CARTESIAN_SPHERE_REP_D_H
