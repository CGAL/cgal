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

#include <CGAL/Nef_S2/Sphere_point.h>

CGAL_BEGIN_NAMESPACE

template<class Kernel>
struct Is_extended_kernel {
       typedef Tag_false value_type;
};

template<class NT>
struct Is_extended_kernel<Extended_homogeneous_3<NT> > {
       typedef Tag_true value_type;
};

template <class T, class Kernel>
class Infimaximal_box {

  typedef typename Kernel::RT               NT;
  typedef Kernel                            Standard_kernel;
  typedef typename Kernel::Point_3          Point_3;
  typedef typename Kernel::Plane_3          Plane_3;
  typedef Point_3                           Standard_point;

 public:
  static bool is_standard(Point_3& p) {
    return true;
  }

  static Point_3 simplify(Point_3& p) {
    return p;
  }

  static Point_3 box_point(Point_3& p, NT d=10000) {
    return p;
  }

  static Standard_point standard_point(Point_3 p, NT d=1) {
    return p;
  }

  static int degree(const typename Kernel::RT& n) {
    return 0;
  }


  static Point_3 create_extended_point(NT x, NT y, NT z) {
    cerr << "function should not be called" << std::endl;
    return Point_3(0,0,0);
  }

  template <typename SNC_decorator_>
    static Point_3 target_for_ray_shot(SNC_decorator_* deco, Point_3& p) {
    return deco->vertices_begin()->point();
  }

  template <typename SNC_constructor_>
  static void create_vertices_of_box_with_plane(SNC_constructor_& C, const Plane_3& h, bool b) {
    cerr << "Constructor not available for this Kernel" << std::endl;
  }

};

template <class Kernel>
class Infimaximal_box<Tag_true, Kernel > {

  typedef typename Kernel::RT               RT;
  typedef typename Kernel::RT::NT           NT;
  typedef typename Kernel::Standard_kernel  Standard_kernel;
  typedef typename Standard_kernel::Point_3 Standard_point;
  typedef typename Kernel::Point_3          Point_3;
  typedef typename Kernel::Plane_3          Plane_3;
  typedef typename Kernel::Vector_3         Vector_3;

  enum Boundary { EXCLUDED=0, INCLUDED=1 };
  
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

  static int degree(const RT& n) {
    return n.degree();
  }

  template <typename SNC_decorator_>
  static Point_3 target_for_ray_shot(SNC_decorator_& deco, Point_3 p) {
    return Kernel::epoint(0, p.hx()[0], 0, p.hy()[0], p.hw()[0], 0, p.hw()[0]);
  }

  static Standard_point standard_point(Point_3 p, NT d=1) {
    return Standard_point(p.hx().eval_at(d),
			  p.hy().eval_at(d),
			  p.hz().eval_at(d),
			  p.hw().eval_at(1));
  }


  static Point_3 box_point(Point_3 p, NT d=10000) {
    return Point_3(p.hx().eval_at(d),
		   p.hy().eval_at(d),
		   p.hz().eval_at(d),
		   p.hw().eval_at(1));
  }

  static Point_3 create_extended_point(NT x, NT y, NT z) {
    return Kernel::epoint(x,0,y,0,z,0,1);
  }

  template <typename SNC_constructor_>
  static void create_vertices_of_box_with_plane(SNC_constructor_& C, const Plane_3& h, bool b) {
    C.create_vertices_of_box_with_plane(h, b);
  }

};

CGAL_END_NAMESPACE
#endif // CGAL_INFIMAXIMAL_BOX_H
