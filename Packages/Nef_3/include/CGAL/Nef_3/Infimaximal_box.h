// ============================================================================
//
// Copyright (c) 1997-2002 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: $
// release_date  : $CGAL_Date: $
//
// file          : include/CGAL/Nef_3/Infimaximal_box.h
// package       : Nef_3
// chapter       : 3D-Nef Polyhedra
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Peter Hachenberger    <hachenberger@mpi-sb.mpg.de>
// maintainer    : Peter Hachenberger    <hachenberger@mpi-sb.mpg.de>
// coordinator   : MPI Saarbruecken
//
// Infimaximal_box.h    answers queries related to the infimaximal box no matter
//                      if it exists (extended kernel) or not (otherwise)
// ============================================================================
#ifndef CGAL_INFIMAXIMAL_BOX_H
#define CGAL_INFIMAXIMAL_BOX_H

#undef _DEBUG
#define _DEBUG 191
#include <CGAL/Nef_3/debug.h>

#include <CGAL/Simple_homogeneous.h>
#include <CGAL/Extended_homogeneous_3.h>

CGAL_BEGIN_NAMESPACE

template <class T, class Kernel>
class Infimaximal_box {

  typedef typename Kernel::Point_3   Point_3;

 public:
  static bool is_standard(Point_3& p) {
    return true;
  }

  static Point_3 simplify(Point_3& p) {
    return p;
  }

  static int degree(const typename Kernel::RT& n) {
    return 0;
  }

  template <typename SNC_structure_>
    static Point_3 target_for_ray_shot(SNC_structure_* sncp, Point_3& p) {
    return sncp->vertices_last();
  }

};

template <class Kernel>
class Infimaximal_box<Tag_true, Kernel > {

  typedef typename Kernel::Point_3   Point_3;
  typedef typename Kernel::RT::NT    NT;
  
 public:
  static bool is_standard(Point_3& p) {
    return Kernel::is_standard(p);
  }

  static Point_3 simplify(Point_3& p) {
    CGAL_assertion(p.hw().degree() == 0);
    int deg = p.hx().degree() > p.hy().degree() 
      ? p.hx().degree() 
      : p.hy().degree();
    deg = p.hz().degree() > deg 
      ? p.hz().degree() 
      : deg;
    return Point_3(p.hx()(deg),p.hy()(deg),p.hz()(deg),p.hw()[0]);
  }

  static int degree(const typename Kernel::RT& n) {
    return n.degree();
  }

  template <typename SNC_structure_>
  static Point_3 target_for_ray_shot(SNC_structure_* sncp, Point_3 p) {
    return Kernel::epoint(0, p.hx()[0], 0, p.hy()[0], p.hw()[0], 0, p.hw()[0]);
  }

  static Point_3 create_extended_point(NT x, NT y, NT z) {
    return Kernel::epoint(x,0,y,0,z,0,1);
  }

};

CGAL_END_NAMESPACE
#endif // CGAL_INFIMAXIMAL_BOX_H
