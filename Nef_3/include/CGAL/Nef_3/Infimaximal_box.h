// Copyright (c) 1997-2002  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Peter Hachenberger    <hachenberger@mpi-sb.mpg.de>

#ifndef CGAL_INFIMAXIMAL_BOX_H
#define CGAL_INFIMAXIMAL_BOX_H

#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 191
#include <CGAL/Nef_2/debug.h>

#include <CGAL/Extended_cartesian.h>
#include <CGAL/Extended_homogeneous.h>

#include <CGAL/Nef_3/SNC_intersection.h>

namespace CGAL {

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

  static bool is_standard(const Point_3& ) {
    return true;
  }

  static bool is_standard(const Plane_3& ) {
    return true;
  }

  static bool check_point_on_plane(Point_3 p, Plane_3 h) {
    return (p.hx()*h.a()+p.hy()*h.b()+p.hz()*h.c()+p.hw()*h.d() == 0);
  }

  static Point_3 simplify(Point_3& p) {
    return p;
  }

  static Standard_point standard_point(Point_3 p, NT =1) {
    return p;
  }

  static Standard_plane standard_plane(Plane_3 p, NT =1) {
    return p;
  }

  static Standard_vector standard_vector(Vector_3 p) {
    return p;
  }

  static void set_size_of_infimaximal_box(const NT& ) {}

  static int degree(const RT& ) {
    return 0;
  }

  static NT get_coeff(const RT& p, unsigned int ) {
    return p;
  }

  static NT eval_at(const RT& p) {
    return p;
  }

  static bool is_infibox_corner(const Point_3& ) {
    return false;
  }

  template <typename Sphere_map>
  static bool is_complex_facet_infibox_intersection(const Sphere_map& ) {
    return false;
  }

  static int type_of_infibox_point(const Point_3& ) {
    return 0;
  }

  template <typename SNC_structure>
  static typename SNC_structure::Volume_const_handle getNirvana(SNC_structure& snc) {
    return snc.volumes_begin();
  }

  template <typename SNC_structure>
  static bool is_beyond_Infibox(typename SNC_structure::SFace_handle , 
				SNC_structure& ) {
    return false;
  }

  static bool x_on_box(const Point_3& ) {
    return false;
  }

  static bool y_on_box(const Point_3& ) {
    return false;
  }

  static bool z_on_box(const Point_3& ) {
    return false;
  }

  template<typename SNC_structure>
  static NT compute_evaluation_constant_for_halfedge_pairup(const SNC_structure& ) {
    return NT(1);
  }

  static void compute_min_max(const Plane_3& , NT orth_coords[3], int& /* min */, int& /* max */) { 
    (void)orth_coords;
  }

  static Point_3 scale_infibox_vertex(const Point_3& ) {
    return Point_3();
  }

  static Point_3 normalize_transformed_vertex(const Point_3& ) {
    return Point_3();
  }

  template <typename SNC_constructor, typename Mark>
  static std::list<typename SNC_constructor::Vertex_handle> 
    create_vertices_on_infibox(SNC_constructor&, 
			       const Plane_3&, const std::list<Point_3>&, 
			       const Mark&, const Mark&, const Mark&) {
    // TODO: warning oder assertion einbauen
    return std::list<typename SNC_constructor::Vertex_handle>();
  }

  template <typename SNC_constructor>
  static std::list<Point_3> find_points_of_box_with_plane(SNC_constructor&, const Plane_3&) {
    return std::list<Point_3>();
  }

  static typename std::list<Point_3>::const_iterator segment_on_side(int /*side_of_point*/, 
							      const std::list<Point_3>& segs) {  
    return segs.begin();
  }

  static Point_3 create_extended_point(NT, NT, NT) {
    std::cerr << "function should not be called" << std::endl;
    CGAL_error_msg("function should not be called");
    return Point_3(0,0,0);
  }

  static Plane_3 create_extended_plane(NT, NT, NT, NT) {
    std::cerr << "function should not be called" << std::endl;
    return Plane_3(1,0,0,0);
  }

  template <typename SNC_decorator, typename Point>
    static Ray_3 get_ray(SNC_decorator& , Point& p) {
    //    return D.point(D.vertex(D.shells_begin(D.volumes_begin())));
    return Ray_3(p, Vector_3(-1,0,0));
  }

  template <typename SNC_constructor_>
  static void create_vertices_of_box_with_plane(SNC_constructor_&, const Plane_3&, bool) {
    std::cerr << "Constructor not available for this Kernel" << std::endl;
  }

  template <typename SNC_constructor>
  static void initialize_infibox_vertices(SNC_constructor& , bool ) {
  }

  template <typename SHalfedge_handle>
  static bool is_sedge_on_infibox(SHalfedge_handle ) {
    return false;
  }

  template <typename Halfedge_handle>
  static bool is_edge_on_infibox(Halfedge_handle ) {
    return false; 
  }

  template <typename Vertex_handle>
  static bool is_redundant_box_vertex(Vertex_handle ) {
    return false;
  }

  template <typename Halfedge_handle>
  static bool is_type4(Halfedge_handle ) {return false;}

  template <typename Halfedge_handle>
  static bool is_type3(Halfedge_handle ) {return false;}
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

  static bool check_point_on_plane(Point_3 p, Plane_3 h) {
    NT x(p.hx().eval_at(100));
    NT y(p.hy().eval_at(100));
    NT z(p.hz().eval_at(100));
    NT w(p.hw().eval_at(100));
    NT a(h.a().eval_at(100));
    NT b(h.b().eval_at(100));
    NT c(h.c().eval_at(100));
    NT d(h.d().eval_at(100));
    return (x*a+y*b+z*c+w*d == 0);
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
    RT::infi_maximal_value() = size;
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
    default: CGAL_error_msg( "wrong value");
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
      default: CGAL_error_msg( "wrong value");
      }
    }
    CGAL_error_msg( "this line of shall not be reached");
    return s1;
  }

  template <typename SNC_constructor>
  static void initialize_infibox_vertices(SNC_constructor& C, bool space) {

#ifdef CGAL_NEF_INDEXED_ITEMS
    int base = Index_generator::get_unique_index();
    for(int i=0; i<11; ++i)
      Index_generator::get_unique_index();
#endif

    C.create_extended_box_corner( 1, 1, 1, space, true
#ifdef CGAL_NEF_INDEXED_ITEMS
				  , base
#endif
				  );
    C.create_extended_box_corner(-1, 1, 1, space, true
#ifdef CGAL_NEF_INDEXED_ITEMS
				  , base
#endif
				  );
    C.create_extended_box_corner( 1,-1, 1, space, true
#ifdef CGAL_NEF_INDEXED_ITEMS
				  , base
#endif
				  );
    C.create_extended_box_corner(-1,-1, 1, space, true
#ifdef CGAL_NEF_INDEXED_ITEMS
				  , base
#endif
				  );
    C.create_extended_box_corner( 1, 1,-1, space, true
#ifdef CGAL_NEF_INDEXED_ITEMS
				  , base
#endif
				  );
    C.create_extended_box_corner(-1, 1,-1, space, true
#ifdef CGAL_NEF_INDEXED_ITEMS
				  , base
#endif
				  );
    C.create_extended_box_corner( 1,-1,-1, space, true
#ifdef CGAL_NEF_INDEXED_ITEMS
				  , base
#endif
				  );
    C.create_extended_box_corner(-1,-1,-1, space, true
#ifdef CGAL_NEF_INDEXED_ITEMS
				  , base
#endif
				  ); 
  }

  template <typename Halfedge_handle>
  static bool is_edge_on_infibox(Halfedge_handle e) {

    Point_3 p  = e->center_vertex()->point();
    if(is_standard(p)) return false;

    Vector_3 v(e->vector());
    CGAL_assertion(p.hw().degree() == 0);
    RT Outer(0,CGAL_NTS abs(p.hw()[0]));

    // TODO: are these lines really redundant??
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
    CGAL_NEF_TRACEN("Point " << p);
    if(is_standard(p)) return false;

    CGAL_NEF_TRACEN("Circle " << sh->circle() << 
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
  static bool is_redundant_box_vertex(const Sphere_map& sm) {
    if(is_standard(sm.point())) return false;
    if(is_infibox_corner(sm.point())) return false;
    typename Sphere_map::SVertex_const_iterator svi;
    for(svi = sm.svertices_begin();
	svi != sm.svertices_end(); ++svi)
      if(!is_edge_on_infibox(svi)) {
	return false;
      }
    return true;
  }
  
  template <typename Sphere_map>
  static bool is_complex_facet_infibox_intersection(const Sphere_map& sm) {
    
    typename Sphere_map::SHalfedge_const_iterator sei;
    bool found = false;
    CGAL_forall_sedges(sei, sm) {
      if(!is_sedge_on_infibox(sei)) {
	if(found)
	  return true;
	else
	  found = true;
      }
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

/*
template <typename RT_>
class Infimaximal_box<Tag_true, CGAL::Pseudo_extended_homogeneous<RT_> > {

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

  static bool check_point_on_plane(Point_3 p, Plane_3 h) {
    std::cerr << p << ", " << h << std::endl;

    NT x(p.hx().eval_at(100));
    NT y(p.hy().eval_at(100));
    NT z(p.hz().eval_at(100));
    NT w(p.hw().eval_at(100));
    NT a(h.a().eval_at(100));
    NT b(h.b().eval_at(100));
    NT c(h.c().eval_at(100));
    NT d(h.d().eval_at(100));
    return (x*a+y*b+z*c+w*d == 0);
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
    RT::infi_maximal_value() = size;
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
    default: CGAL_error_msg( "wrong value");
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
      default: CGAL_error_msg( "wrong value");
      }
    }
    CGAL_error_msg( "this line of shall not be reached");
    return s1;
  }

  template <typename SNC_constructor>
  static void initialize_infibox_vertices(SNC_constructor& C, bool space) {

#ifdef CGAL_NEF_INDEXED_ITEMS
    int base = Index_generator::get_unique_index();
    for(int i=0; i<11; ++i)
      Index_generator::get_unique_index();
#endif

    C.create_extended_box_corner( 1, 1, 1, space, true
#ifdef CGAL_NEF_INDEXED_ITEMS
				  , base
#endif
				  );
    C.create_extended_box_corner(-1, 1, 1, space, true
#ifdef CGAL_NEF_INDEXED_ITEMS
				  , base
#endif
				  );
    C.create_extended_box_corner( 1,-1, 1, space, true
#ifdef CGAL_NEF_INDEXED_ITEMS
				  , base
#endif
				  );
    C.create_extended_box_corner(-1,-1, 1, space, true
#ifdef CGAL_NEF_INDEXED_ITEMS
				  , base
#endif
				  );
    C.create_extended_box_corner( 1, 1,-1, space, true
#ifdef CGAL_NEF_INDEXED_ITEMS
				  , base
#endif
				  );
    C.create_extended_box_corner(-1, 1,-1, space, true
#ifdef CGAL_NEF_INDEXED_ITEMS
				  , base
#endif
				  );
    C.create_extended_box_corner( 1,-1,-1, space, true
#ifdef CGAL_NEF_INDEXED_ITEMS
				  , base
#endif
				  );
    C.create_extended_box_corner(-1,-1,-1, space, true
#ifdef CGAL_NEF_INDEXED_ITEMS
				  , base
#endif
				  ); 
  }

  template <typename Halfedge_handle>
  static bool is_edge_on_infibox(Halfedge_handle e) {

    Point_3 p  = e->center_vertex()->point();
    if(is_standard(p)) return false;

    Vector_3 v(e->vector());
    CGAL_assertion(p.hw().degree() == 0);
    RT Outer(0,CGAL_NTS abs(p.hw()[0]));

    // TODO: are these lines really redundant??
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
    CGAL_NEF_TRACEN("Point " << p);
    if(is_standard(p)) return false;

    CGAL_NEF_TRACEN("Circle " << sh->circle() << 
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
  static bool is_redundant_box_vertex(const Sphere_map& sm) {
    if(is_standard(sm.point())) return false;
    if(is_infibox_corner(sm.point())) return false;
    typename Sphere_map::SVertex_const_iterator svi;
    for(svi = sm.svertices_begin();
	svi != sm.svertices_end(); ++svi)
      if(!is_edge_on_infibox(svi)) {
	return false;
      }
    return true;
  }
  
  template <typename Sphere_map>
  static bool is_complex_facet_infibox_intersection(const Sphere_map& sm) {
    
    typename Sphere_map::SHalfedge_const_iterator sei;
    bool found = false;
    CGAL_forall_sedges(sei, sm) {
      if(!is_sedge_on_infibox(sei)) {
	if(found)
	  return true;
	else
	  found = true;
      }
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
*/

} //namespace CGAL
#endif // CGAL_INFIMAXIMAL_BOX_H
