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

#include <CGAL/Extended_cartesian.h>
#include <CGAL/Extended_homogeneous.h>

#include <CGAL/Nef_3/SNC_intersection.h>

CGAL_BEGIN_NAMESPACE

template<class Kernel>
struct Is_extended_kernel {
       typedef Tag_false value_type;
};

template<class NT>
struct Is_extended_kernel<Extended_homogeneous<NT> > {
       typedef Tag_true value_type;
};

template<class NT>
struct Is_extended_kernel<Extended_cartesian<NT> > {
       typedef Tag_true value_type;
};

template <class T, class Kernel>
class Infimaximal_box {

 public:
  typedef typename Kernel::RT               NT;
  typedef typename Kernel::RT               RT;
  typedef Kernel                            Standard_kernel;
  typedef typename Kernel::Point_3          Point_3;
  typedef typename Kernel::Vector_3         Vector_3;
  typedef typename Kernel::Plane_3          Plane_3;
  typedef typename Kernel::Ray_3            Ray_3;
  typedef typename Kernel::Aff_transformation_3   Aff_transformation_3;
  typedef Plane_3                           Standard_plane;
  typedef Vector_3                          Standard_vector;
  typedef Point_3                           Standard_point;

  static bool standard_kernel() {
    return true;
  }

  static bool extended_kernel() {
    return false;
  }

  static bool is_standard(const Point_3& p) {
    return true;
  }

  static bool is_standard(const Plane_3& p) {
    return true;
  }

  static bool check_point_on_plane(Point_3 p, Plane_3 h) {
    return (p.hx()*h.a()+p.hy()*h.b()+p.hz()*h.c()+p.hw()*h.d() == 0);
  }

  static Point_3 simplify(Point_3& p) {
    return p;
  }

  static Standard_point standard_point(Point_3 p, NT d=1) {
    return p;
  }

  static Standard_plane standard_plane(Plane_3 p, NT d=1) {
    return p;
  }

  static Standard_vector standard_vector(Vector_3 p) {
    return p;
  }

  static void set_size_of_infimaximal_box(const NT& size) {}

  static int degree(const RT& n) {
    return 0;
  }

  static NT get_coeff(const RT& p, unsigned int i) {
    return p;
  }

  static NT eval_at(const RT& p) {
    return p;
  }

  static bool is_infibox_corner(const Point_3& p) {
    return false;
  }

  template <typename Sphere_map>
  static bool is_complex_facet_infibox_intersection(const Sphere_map& sm) {
    return false;
  }

  static int type_of_infibox_point(const Point_3& p) {
    return 0;
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

  static bool x_on_box(const Point_3& p) {
    return false;
  }

  static bool y_on_box(const Point_3& p) {
    return false;
  }

  static bool z_on_box(const Point_3& p) {
    return false;
  }

  template<typename SNC_structure>
  static NT compute_evaluation_constant_for_halfedge_pairup(const SNC_structure& snc) {
    return NT(1);
  }

  static void compute_min_max(const Plane_3& h, NT orth_coords[3], int& min, int& max) { }

  static Point_3 scale_infibox_vertex(const Point_3 pin) {
    return Point_3();
  }

  static Point_3 normalize_transformed_vertex(const Point_3& p) {
    return Point_3();
  }

  template <typename SNC_constructor, typename Mark>
  static std::list<typename SNC_constructor::Vertex_handle> 
    create_vertices_on_infibox(SNC_constructor& C, 
			       const Plane_3& h, const std::list<Point_3>& points, 
			       const Mark& bnd, const Mark& inside, const Mark& outside) {
    return std::list<typename SNC_constructor::Vertex_handle>();
  }

  template <typename SNC_constructor>
  static std::list<Point_3> find_points_of_box_with_plane(SNC_constructor& C, const Plane_3& h) {
    return std::list<Point_3>();
  }

  static typename std::list<Point_3>::const_iterator segment_on_side(int side_of_point, 
							      const std::list<Point_3>& segs) {  
    return segs.begin();
  }

  static Point_3 create_extended_point(NT x, NT y, NT z) {
    std::cerr << "function should not be called" << std::endl;
    CGAL_assertion_msg(0,"function should not be called");
    return Point_3(0,0,0);
  }

  static Plane_3 create_extended_plane(NT a, NT b, NT c, NT d) {
    std::cerr << "function should not be called" << std::endl;
    return Plane_3(1,0,0,0);
  }

  template <typename SNC_decorator, typename Point>
    static Ray_3 get_ray(SNC_decorator& D, Point& p) {
    //    return D.point(D.vertex(D.shells_begin(D.volumes_begin())));
    return Ray_3(p, Vector_3(-1,0,0));
  }

  template <typename SNC_constructor_>
  static void create_vertices_of_box_with_plane(SNC_constructor_& C, const Plane_3& h, bool b) {
    std::cerr << "Constructor not available for this Kernel" << std::endl;
  }

  template <typename SNC_constructor>
  static void initialize_infibox_vertices(SNC_constructor& C, bool space) {
  }

  template <typename SHalfedge_handle>
  static bool is_sedge_on_infibox(SHalfedge_handle sh) {
    return false;
  }

  template <typename Halfedge_handle>
  static bool is_edge_on_infibox(Halfedge_handle e) {
    return false; 
  }

  template <typename Halfedge_handle>
  static bool is_type4(Halfedge_handle e) {return false;}

  template <typename Halfedge_handle>
  static bool is_type3(Halfedge_handle e) {return false;}
};

template <class Kernel>
class Infimaximal_box<Tag_true, Kernel> {

 public:
  typedef typename Kernel::RT                    RT;
  typedef typename Kernel::RT::NT                NT;
  typedef typename Kernel::Standard_kernel       Standard_kernel;
  typedef typename Standard_kernel::Point_3      Standard_point;
  typedef typename Standard_kernel::Plane_3      Standard_plane;
  typedef typename Standard_kernel::Vector_3     Standard_vector;
  typedef typename Kernel::Point_3               Point_3;
  typedef typename Kernel::Plane_3               Plane_3;
  typedef typename Kernel::Vector_3              Vector_3;
  typedef typename Kernel::Segment_3             Segment_3;
  typedef typename Kernel::Direction_3           Direction_3;
  typedef typename Kernel::Aff_transformation_3  Aff_transformation_3;

  //  typedef typename SNC_structure::Sphere_point   Sphere_point;
  //  typedef typename SNC_structure::Sphere_circle  Sphere_circle;

  //  typedef SM_decorator<SNC_structure>        SM_decorator;

  enum Boundary { EXCLUDED=0, INCLUDED=1 };
  
  static const int RADIUS = 10000000;

  static bool standard_kernel() {
    return false;
  }

  static bool extended_kernel() {
    return true;
  }

  static bool is_standard(const Point_3& p) {
    CGAL_assertion(p.hw().degree()==0);
    return (p.hx().degree()==0 && p.hy().degree()==0 && p.hz().degree()==0);
  }

  static bool is_standard(const Plane_3& p) {
    return (p.d().degree() == 0);
  }

  static int type_of_infibox_point(const Point_3& p) {
    int res = 0;
    RT W(NT(0),p.hw()[0]);
    if(CGAL_NTS abs(p.hx()) == W) ++res;
    if(CGAL_NTS abs(p.hy()) == W) ++res;
    if(CGAL_NTS abs(p.hz()) == W) ++res;
    return res;
  }

  static bool is_infibox_corner(const Point_3& p) {
    return type_of_infibox_point(p) == 3;
  }

  static Point_3 simplify(Point_3& p) {
    CGAL_assertion(p.hw().degree() == 0);
    int deg = p.hx().degree() > p.hy().degree() 
      ? p.hx().degree() 
      : p.hy().degree();
    deg = p.hz().degree() > deg 
      ? p.hz().degree() 
      : deg;
    return Point_3(p.hx().degree() == deg ? p.hx()[deg] : 0,
		   p.hy().degree() == deg ? p.hy()[deg] : 0,
		   p.hz().degree() == deg ? p.hz()[deg] : 0,
		   p.hw()[0]);
  }

  static int degree(const RT& p) {
    return p.degree();
  }

  static NT get_coeff(const RT& p, unsigned int i) {
    return p[i];
  }

  static NT eval_at(const RT& p, NT d=RADIUS) {
    return p.eval_at(d);
  }

  static Point_3 target_for_ray_shoot_on_minus_x_direction(Point_3 p) {
    CGAL_warning(is_standard(p));
    return Point(RT(0,-1), RT(p.hy()[0]), RT(p.hz()[0]), RT(p.hw()[0]));
  }

  static bool check_point_on_plane(Point_3 p, Plane_3 h) {
    NT x(p.hx().eval_at(100));
    NT y(p.hy().eval_at(100));
    NT z(p.hz().eval_at(100));
    NT w(p.hw().eval_at(100));
    NT d(h.d().eval_at(100));
    return (x*h.a()+y*h.b()+z*h.c()+w*d == 0);
  }

  static Standard_point standard_point(Point_3 p, NT d=1) {
    return Standard_point(p.hx().eval_at(d),
			  p.hy().eval_at(d),
			  p.hz().eval_at(d),
			  p.hw().eval_at(1));
  }

  static Standard_plane standard_plane(Plane_3 p, NT d=1) {
    return Standard_plane(p.a().eval_at(1),
			  p.b().eval_at(1),
			  p.c().eval_at(1),
			  p.d().eval_at(d));
  }

  static Standard_vector standard_vector(Vector_3 p) {
    return Standard_vector(p.hx().eval_at(1),
			   p.hy().eval_at(1),
			   p.hz().eval_at(1));
  }

  static void set_size_of_infimaximal_box(NT size) {
    RT::set_R(size);
  }

  static bool x_on_box(const Point_3& p) {
    return CGAL_NTS abs(p.hx()) == RT(0,p.hw()[0]);
  }

  static bool y_on_box(const Point_3& p) {
    return CGAL_NTS abs(p.hy()) == RT(0,p.hw()[0]);
  }

  static bool z_on_box(const Point_3& p) {
    return CGAL_NTS abs(p.hz()) == RT(0,p.hw()[0]);
  }

  static Point_3 create_extended_point(NT x, NT y, NT z) {
    return Point_3(RT(0,x), RT(0,y), RT(0,z), RT(1));
  }

  static Plane_3 create_extended_plane(NT a, NT b, NT c, NT d) {
    return Plane_3(a,b,c,RT(0,d));
  }

  static Point_3 create_extended_point(NT a0,NT a1,NT b0,NT b1,NT c0,NT c1,NT d) {
    return Point_3(RT(a0,a1),RT(b0,b1),RT(c0,c1),RT(d));
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

  template <typename SNC_constructor, typename Mark>
  static std::list<typename SNC_constructor::Vertex_handle> 
    create_vertices_on_infibox(SNC_constructor& C, 
			       const Plane_3& h, const std::list<Point_3>& points, 
			       const Mark& bnd, const Mark& inside, const Mark& outside) {
    return C.create_vertices_on_infibox(h,points,bnd,inside,outside);
  }

  template <typename SNC_constructor>
  static void create_vertices_of_box_with_plane(SNC_constructor& C, const Plane_3& h, bool b) {
    C.create_vertices_of_box_with_plane(h, b);
  }

  static void compute_min_max(const Plane_3& h, NT orth_coords[3], int& min, int& max) {
    Vector_3 orth = h.orthogonal_vector();
    
    orth_coords[0] = CGAL_NTS abs(orth.hx()[0]);
    orth_coords[1] = CGAL_NTS abs(orth.hy()[0]);
    orth_coords[2] = CGAL_NTS abs(orth.hz()[0]);
    
    max = 0;
    if(orth_coords[1] > orth_coords[0])
      max = 1;
    if(orth_coords[2] > orth_coords[max])
      max = 2;   
    
    min = 0;
    if(orth_coords[1] < orth_coords[0])
      min = 1;
    if(orth_coords[2] < orth_coords[min])
      min = 2;
  }

  template<typename SNC_structure>
  static NT compute_evaluation_constant_for_halfedge_pairup(const SNC_structure& snc) {
    NT eval = 0;
    typename SNC_structure::Vertex_const_iterator v;
    CGAL_forall_vertices(v, snc) {
      Point_3 p(v->point());
      if(p.hx()[0] > eval) eval = p.hx()[0];
      if(p.hy()[0] > eval) eval = p.hy()[0];
      if(p.hz()[0] > eval) eval = p.hz()[0];
    }
    eval *= 4;
    if(eval == 0) return 1;
    return eval;
  }

  static Point_3 scale_infibox_vertex(const Point_3 pin, const Aff_transformation_3& aff) {
    RT lx(pin.hx()[0]);
    RT ly(pin.hy()[0]);
    RT lz(pin.hz()[0]);
    RT hx(pin.hx()-lx);
    RT hy(pin.hy()-ly);
    RT hz(pin.hz()-lz);
    RT hw(pin.hw());
    Point_3 p(Point_3(lx,ly,lz,hw).transform(aff));
    return Point_3(hx+p.hx(),hy+p.hy(),hz+p.hz(),hw);    
  }

  static Point_3 normalize_transformed_vertex(const Point_3& p) {
    int i=0;
    if(CGAL_NTS abs(p.hx()) < CGAL_NTS abs(p.hy()))
      if(CGAL_NTS abs(p.hy()) < CGAL_NTS abs(p.hz()))
	i = 2;
      else
	i = 1;
    else if(CGAL_NTS abs(p.hx()) < CGAL_NTS abs(p.hz()))
      i = 2;

    switch(i) {
    case 0:
      CGAL_assertion(p.hx().degree() == 1);
      if(p.hx()[1] > 0)
	return Point_3(RT(0,p.hx()[1]*p.hw()[0]),
		       RT(p.hy()[0]*p.hx()[1]-p.hx()[0]*p.hy()(1),p.hy()(1)*p.hw()[0]),
		       RT(p.hz()[0]*p.hx()[1]-p.hx()[0]*p.hz()(1),p.hz()(1)*p.hw()[0]),
		       RT(p.hw()[0]*p.hx()[1]));
      else
	return Point_3(RT(0,-p.hx()[1]*p.hw()[0]),
		       RT(p.hy()[0]*p.hx()[1]-p.hx()[0]*p.hy()(1),-p.hy()(1)*p.hw()[0]),
		       RT(p.hz()[0]*p.hx()[1]-p.hx()[0]*p.hz()(1),-p.hz()(1)*p.hw()[0]),
		       RT(p.hw()[0]*p.hx()[1]));
    case 1:
      CGAL_assertion(p.hy().degree() == 1);
      if(p.hy()[1] > 0)
	return Point_3(RT(p.hx()[0]*p.hy()[1]-p.hy()[0]*p.hx()(1),p.hx()(1)*p.hw()[0]),
		       RT(0,p.hy()[1]*p.hw()[0]),
		       RT(p.hz()[0]*p.hy()[1]-p.hy()[0]*p.hz()(1),p.hz()(1)*p.hw()[0]),
		       RT(p.hw()[0]*p.hy()[1]));
      else
	return Point_3(RT(p.hx()[0]*p.hy()[1]-p.hy()[0]*p.hx()(1),-p.hx()(1)*p.hw()[0]),
		       RT(0,-p.hy()[1]*p.hw()[0]),
		       RT(p.hz()[0]*p.hy()[1]-p.hy()[0]*p.hz()(1),-p.hz()(1)*p.hw()[0]),
		       RT(p.hw()[0]*p.hy()[1]));
    case 2:
      CGAL_assertion(p.hz().degree() == 1);
      if(p.hz()[1] > 0)
	return Point_3(RT(p.hx()[0]*p.hz()[1]-p.hz()[0]*p.hx()(1),p.hx()(1)*p.hw()[0]),
		       RT(p.hy()[0]*p.hz()[1]-p.hz()[0]*p.hy()(1),p.hy()(1)*p.hw()[0]),
		       RT(0,p.hz()[1]*p.hw()[0]),
		       RT(p.hw()[0]*p.hz()[1]));
      else
	return Point_3(RT(p.hx()[0]*p.hz()[1]-p.hz()[0]*p.hx()(1),-p.hx()(1)*p.hw()[0]),
		       RT(p.hy()[0]*p.hz()[1]-p.hz()[0]*p.hy()(1),-p.hy()(1)*p.hw()[0]),
		       RT(0,-p.hz()[1]*p.hw()[0]),
		       RT(p.hw()[0]*p.hz()[1]));
    default: CGAL_assertion_msg(false, "wrong value");
    }
    return Point_3();
  }

  static typename std::list<Point_3>::const_iterator segment_on_side(int side_of_point, 
								     const std::list<Point_3>& segs) {  

    typename std::list<Point_3>::const_iterator s1,t1;
    for(s1 = segs.begin(); s1 != segs.end(); ++s1) {
      t1 = s1;
      ++t1;
      if(t1 == segs.end()) t1 = segs.begin();
      switch(side_of_point) {
      case  1: 
	if( s1->hx()(1) != s1->hw()) continue; 
	if( t1->hx()(1) != t1->hw()) continue;
	return s1;
      case -1: 
	if(-s1->hx()(1) != s1->hw()) continue;
	if(-t1->hx()(1) != t1->hw()) continue; 
	return s1;
      case  2: 
	if( s1->hy()(1) != s1->hw()) continue;
	if( t1->hy()(1) != t1->hw()) continue;
	return s1;	break;
      case -2: 
	if(-s1->hy()(1) != s1->hw()) continue; 
	if(-t1->hy()(1) != t1->hw()) continue; 
	return s1;
      case  3: 
	if( s1->hz()(1) != s1->hw()) continue; 
	if( t1->hz()(1) != t1->hw()) continue; 
	return s1;
      case -3: 
	if(-s1->hz()(1) != s1->hw()) continue; 
	if(-t1->hz()(1) != t1->hw()) continue; 
	return s1;
      default: CGAL_assertion_msg(false, "wrong value");
      }
    }
    CGAL_assertion_msg(false, "this line of shall not be reached");
    return s1;
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

  template <typename SHalfedge_handle>
  static bool is_edge_on_infibox(SHalfedge_handle e) {

    Point_3 p  = e->center_vertex()->point();
    if(is_standard(p)) return false;

    Vector_3 v(e->vector());
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

  template <typename SHalfedge_handle>
  static bool is_sedge_on_infibox(SHalfedge_handle sh) {

    Point_3 p = sh->source()->center_vertex()->point();
    TRACEN("Point " << p);
    if(is_standard(p)) return false;

    TRACEN("Circle " << sh->circle() << 
	   " has signum " << sign_of(sh->circle()));
    CGAL_assertion(p.hw().degree() == 0);
    RT R(0,CGAL_NTS abs(p.hw()[0]));

    if((sh->circle().a() == 0 && sh->circle().b() == 0 
	&& CGAL_NTS abs(p.hz())== R) || 
       (sh->circle().a() == 0 && sh->circle().c() == 0 
	&& CGAL_NTS abs(p.hy())== R) ||
       (sh->circle().b() == 0 && sh->circle().c() == 0 
	&& CGAL_NTS abs(p.hx())== R))
      if(is_edge_on_infibox(sh->source()) && 
	 is_edge_on_infibox(sh->twin()->source()))
	return true;

    return false;
  }
  
  template <typename Sphere_map>
  static bool is_complex_facet_infibox_intersection(const Sphere_map& sm) {
    
    typename Sphere_map::SHalfedge_const_iterator sei;
    bool found = false;
    CGAL_forall_sedges(sei, sm) {
      if(!is_sedge_on_infibox(sei))
	if(found)
	  return true;
	else
	  found = true;
    }
    return false;
  }

  template <typename SNC_constructor>
  static std::list<Point_3> find_points_of_box_with_plane(SNC_constructor& C, const Plane_3& h) {
    return C.find_points_of_box_with_plane(h);
  }

  template <typename Halfedge_handle>
  static bool is_type4(Halfedge_handle e) {

    Point_3 p(e->center_vertex()->point());
    Direction_3 d(e->vector());

    if((CGAL_NTS abs(p.hx()) == CGAL_NTS abs(p.hw()) || 
	d == Direction_3(1,0,0) ||
	d == Direction_3(-1,0,0)) &&
       (CGAL_NTS abs(p.hy()) == CGAL_NTS abs(p.hw()) || 
	d == Direction_3(0,1,0) ||
	d == Direction_3(0,-1,0)) &&
       (CGAL_NTS abs(p.hz()) == CGAL_NTS abs(p.hw()) || 
	d == Direction_3(0,0,1) ||
	d == Direction_3(0,0,-1)))
      return true;
    return false;
  }

  template <typename Halfedge_handle>
  static bool is_type3(Halfedge_handle e) {

    Point_3 p(e->center_vertex()->point());
    Direction_3 d(e->vector());
    
    if(d == Direction_3(1,0,0) || d == Direction_3(-1,0,0)) {
      if(CGAL_NTS abs(p.hy()) == CGAL_NTS abs(p.hw()))
	return true;
      return false;
    }
    if(d == Direction_3(0,1,0) || d == Direction_3(0,-1,0)) {
      if(CGAL_NTS abs(p.hx()) == CGAL_NTS abs(p.hw()))
	return true;
      return false;
    }
    if(d == Direction_3(0,0,1) || d == Direction_3(0,0,-1)) {
      if(CGAL_NTS abs(p.hy()) == CGAL_NTS abs(p.hw()))
	return true;
      return false;
    }

    return false;
  }

};

CGAL_END_NAMESPACE
#endif // CGAL_INFIMAXIMAL_BOX_H
