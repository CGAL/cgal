// Copyright (c) 1997-2002  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Peter Hachenberger    <hachenberger@mpi-sb.mpg.de>
#ifndef CGAL_INFIMAXIMAL_BOX_H
#define CGAL_INFIMAXIMAL_BOX_H

#undef _DEBUG
#define _DEBUG 191
#include <CGAL/Nef_3/debug.h>

#include <CGAL/Simple_homogeneous.h>
#include <CGAL/Extended_homogeneous_3.h>

// #include <CGAL/Nef_S2/Sphere_circle.h>
// #include <CGAL/Nef_3/SNC_SM_decorator.h>

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
  typedef typename Kernel::RT               RT;
  typedef Kernel                            Standard_kernel;
  typedef typename Kernel::Point_3          Point_3;
  typedef typename Kernel::Plane_3          Plane_3;
  typedef typename Kernel::Ray_3            Ray_3;
  typedef typename Kernel::Direction_3      Direction_3;
  typedef Point_3                           Standard_point;

 public:

  static bool standard_Kernel() {
    return true;
  }

  static bool extended_Kernel() {
    return false;
  }

  static bool is_standard(const Point_3& p) {
    return true;
  }

  static bool is_standard(Plane_3& p) {
    return true;
  }

  static Point_3 simplify(Point_3& p) {
    return p;
  }

  static Point_3 box_point(const Point_3& p, NT d=10000) {
    return p;
  }

  static Standard_point standard_point(Point_3 p, NT d=1) {
    return p;
  }

  static int degree(const typename Kernel::RT& n) {
    return 0;
  }

  static NT get_coeff(const RT& p, unsigned int i) {
    return p;
  }

  template <typename SNC_structure>
  static typename SNC_structure::Volume_handle getNirvana(SNC_structure& snc) {
    return snc.volumes_begin();
  }

  template <typename SNC_structure>
  static bool is_beyond_Infibox(typename SNC_structure::SFace_handle sf, 
				SNC_structure& snc) {
    return false;
  }

  static Point_3 create_extended_point(NT x, NT y, NT z) {
    cerr << "function should not be called" << std::endl;
    while(1);
    return Point_3(0,0,0);
  }

  template <typename SNC_decorator, typename Point>
    static Ray_3 get_ray(SNC_decorator& D, Point& p) {
    //    return D.point(D.vertex(D.shells_begin(D.volumes_begin())));
    return Ray_3(p, Direction_3(-1,0,0));
  }

  template <typename SNC_constructor_>
  static void create_vertices_of_box_with_plane(SNC_constructor_& C, const Plane_3& h, bool b) {
    cerr << "Constructor not available for this Kernel" << std::endl;
  }

  template <typename SNC_constructor>
  static void initialize_infibox_vertices(SNC_constructor& C, bool space) {
  }

  template <typename SHalfedge_handle, typename SNC_decorator>
  static bool is_sedge_on_infibox(SHalfedge_handle sh, SNC_decorator D) {
    return false;
  }

  template <typename Halfedge_handle, typename SNC_decorator>
  static bool is_edge_on_infibox(Halfedge_handle e, SNC_decorator D) {
    return false; 
  }

};

template <class Kernel>
class Infimaximal_box<Tag_true, Kernel> {

  typedef typename Kernel::RT                    RT;
  typedef typename Kernel::RT::NT                NT;
  typedef typename Kernel::Standard_kernel       Standard_kernel;
  typedef typename Standard_kernel::Point_3      Standard_point;
  typedef typename Kernel::Point_3               Point_3;
  typedef typename Kernel::Plane_3               Plane_3;
  typedef typename Kernel::Vector_3              Vector_3;
  //  typedef typename SNC_structure::Sphere_circle  Sphere_circle;
  //  typedef SNC_SM_decorator<SNC_structure>        SM_decorator;

  enum Boundary { EXCLUDED=0, INCLUDED=1 };
  
 public:


  static bool standard_Kernel() {
    return false;
  }

  static bool extended_Kernel() {
    return true;
  }

  static bool is_standard(const Point_3& p) {
    return Kernel::is_standard(p);
  }

  static bool is_standard(Plane_3& p) {
    return (p.d().degree() == 0);
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

  static int degree(const RT& p) {
    return p.degree();
  }

  static NT get_coeff(const RT& p, unsigned int i) {
    return p[i];
  }

  static Point_3 target_for_ray_shoot_on_minus_x_direction(Point_3 p) {
    CGAL_warning( Kernel::is_standard(p));
    return Kernel::epoint( -1, 0, 0, p.hy()[0], 0, p.hz()[0], p.hw()[0]);
  }

  static Standard_point standard_point(Point_3 p, NT d=1) {
    return Standard_point(p.hx().eval_at(d),
			  p.hy().eval_at(d),
			  p.hz().eval_at(d),
			  p.hw().eval_at(1));
  }


  static Point_3 box_point(const Point_3& p, NT d=10000) {
    return Point_3(p.hx().eval_at(d),
		   p.hy().eval_at(d),
		   p.hz().eval_at(d),
		   p.hw().eval_at(1));
  }

  static Point_3 create_extended_point(NT x, NT y, NT z) {
    return Kernel::epoint(x,0,y,0,z,0,1);
  }

  template <typename SNC_structure>
  static typename SNC_structure::Volume_handle getNirvana(SNC_structure& snc) {
    return  ++(snc.volumes_begin());
  }

  template <typename SNC_structure>
  static bool is_beyond_Infibox(typename SNC_structure::SFace_handle sf, 
				SNC_structure& snc) {
    typename SNC_structure::SNC_decorator D(snc);
    return (D.volume(sf) == snc.volumes_begin());
  }

  template <typename SNC_constructor>
  static void create_vertices_of_box_with_plane(SNC_constructor& C, const Plane_3& h, bool b) {
    C.create_vertices_of_box_with_plane(h, b);
  }

  template <typename SNC_constructor>
  static void initialize_infibox_vertices(SNC_constructor& C, bool space) {
    C.create_extended_box_corner( 1, 1, 1, space );
    C.create_extended_box_corner(-1, 1, 1, space );
    C.create_extended_box_corner( 1,-1, 1, space );
    C.create_extended_box_corner(-1,-1, 1, space );
    C.create_extended_box_corner( 1, 1,-1, space );
    C.create_extended_box_corner(-1, 1,-1, space );
    C.create_extended_box_corner( 1,-1,-1, space );
    C.create_extended_box_corner(-1,-1,-1, space ); 
  }

  template <typename SNC_decorator>
  static bool is_edge_on_infibox(typename SNC_decorator::Halfedge_handle e, SNC_decorator& D) {

    Point_3 p  = D.point(D.vertex(e));
    if(Kernel::is_standard(p)) return false;

    Vector_3 v(D.tmp_point(e));
    CGAL_assertion(p.hw().degree() == 0);
    RT Outer(0,CGAL_NTS abs(p.hw()[0]));

    if(CGAL_NTS abs(p.hx()) == Outer && 
       ((p.hx() > 0 && v.hx() > 0)||(p.hx() < 0 && v.hx() < 0))) return false;
    if(CGAL_NTS abs(p.hy()) == Outer && 
       ((p.hy() > 0 && v.hy() > 0)||(p.hy() < 0 && v.hy() < 0))) return false;
    if(CGAL_NTS abs(p.hz()) == Outer && 
       ((p.hz() > 0 && v.hz() > 0)||(p.hz() < 0 && v.hz() < 0))) return false;

    if(CGAL_NTS abs(p.hx()) == Outer && v.hx() == 0) return true;
    if(CGAL_NTS abs(p.hy()) == Outer && v.hy() == 0) return true;
    if(CGAL_NTS abs(p.hz()) == Outer && v.hz() == 0) return true;

    return false; 
  }

  template <typename SHalfedge_handle, typename SNC_decorator>
  static bool is_sedge_on_infibox(SHalfedge_handle sh, SNC_decorator D) {
    typedef typename SNC_decorator::SNC_structure  SNC_structure;
    typedef typename SNC_structure::Sphere_circle  Sphere_circle;
    typedef typename SNC_decorator::SM_decorator   SM_decorator;

    Point_3 p = D.point(D.vertex(sh));
    TRACEN("Point " << p);
    if(Kernel::is_standard(p)) return false;

    Sphere_circle c(sh->tmp_circle());
    TRACEN("Circle " << c << " has signum " << sign_of(c));
    CGAL_assertion(p.hw().degree() == 0);
    RT R(0,CGAL_NTS abs(p.hw()[0]));

    SM_decorator SD(D.vertex(sh));
    if((c.a() == 0 && c.b() == 0 && CGAL_NTS abs(p.hz())== R) || 
       (c.a() == 0 && c.c() == 0 && CGAL_NTS abs(p.hy())== R) ||
       (c.b() == 0 && c.c() == 0 && CGAL_NTS abs(p.hx())== R))
      if(is_edge_on_infibox(SD.source(sh),D) && 
	 is_edge_on_infibox(SD.target(sh),D))
	return true;

    return false;
  }

};

CGAL_END_NAMESPACE
#endif // CGAL_INFIMAXIMAL_BOX_H
