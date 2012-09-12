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
// Author(s)     : Michael Seel       <seel@mpi-sb.mpg.de> 
//                 Miguel Granados    <granados@mpi-sb.mpg.de>
//                 Susan Hert         <hert@mpi-sb.mpg.de>
//                 Lutz Kettner       <kettner@mpi-sb.mpg.de>
//                 Peter Hachenberger <hachenberger@mpi-sb.mpg.de>
#ifndef CGAL_SNC_CONSTRUCTOR_H
#define CGAL_SNC_CONSTRUCTOR_H

#include <CGAL/function_objects.h> 
#include <CGAL/Circulator_project.h>
#include <CGAL/Nef_S2/Normalizing.h>
#include <CGAL/Nef_3/bounded_side_3.h>
#include <CGAL/Nef_3/Pluecker_line_3.h>
#include <CGAL/Nef_3/SNC_decorator.h>
#include <CGAL/Nef_3/SNC_SM_overlayer.h>
#include <CGAL/Nef_S2/SM_point_locator.h>
#include <CGAL/Nef_3/SNC_sphere_map.h>
#include <CGAL/Nef_3/SNC_intersection.h>
#include <CGAL/Nef_3/SNC_external_structure.h>
#ifdef SM_VISUALIZOR
#include <CGAL/Nef_3/SNC_SM_visualizor.h>
#endif // SM_VISUALIZOR
#include <map>
#include <list>
#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 43
#include <CGAL/Nef_2/debug.h>

namespace CGAL {

template <typename Infi_box, typename Vertex_handle>
struct Frame_point_lt {

  Frame_point_lt() {}
  bool operator()(Vertex_handle v1, Vertex_handle v2) const {
    if(v1->point() != v2->point())
      return CGAL::lexicographically_xyz_smaller(v1->point(),v2->point());
    if(!Infi_box::is_complex_facet_infibox_intersection(*v1) && 
       Infi_box::is_complex_facet_infibox_intersection(*v2))
      return true;
    return false;
  }
};

template <typename T>
struct circle_lt {
  
  int m_max;
  typedef typename T::Point_3                   Point_3;
  typedef typename T::RT                        RT;

  circle_lt(int m) :m_max(m) {};
  bool operator()(const Point_3& p1, const Point_3& p2) const { 
        
    const Quotient<RT> zero(RT(0));
    Quotient<RT> x[2];
    Quotient<RT> y[2];

    switch(m_max) {
    case 0:
      x[0] = p1.y(); 
      y[0] = p1.z();
      x[1] = p2.y(); 
      y[1] = p2.z();  
      break;
    case 1:
      x[0] = p1.x(); 
      y[0] = p1.z();
      x[1] = p2.x(); 
      y[1] = p2.z();  
      break;
    case 2:
      x[0] = p1.x(); 
      y[0] = p1.y();
      x[1] = p2.x(); 
      y[1] = p2.y();  
      break;
    }
    
    if(y[0] >= zero) {
      if(y[1] < zero) return false;
      if(x[0] != x[1]) return (x[0]<x[1]);
      if(x[0] > zero) return (y[0]>y[1]);
      else return (y[0]<y[1]);
    }
    else {
      if(y[1] >= zero) return true;
      if(x[0]!=x[1]) return(x[0]>x[1]);
      if(x[0] > zero) return (y[0]>y[1]);
      else return  (y[0]<y[1]);
    }
    CGAL_error_msg( "control should not reach this line");
    return false;
  }
};

// ----------------------------------------------------------------------------
// SNC_constructor 
// ----------------------------------------------------------------------------

/*{\Manpage{SNC_constructor}{SNC}{overlay functionality}{O}}*/

template <typename Items, typename SNC_structure_>
class SNC_constructor_base : public SNC_decorator<SNC_structure_>
{ 
public:
  typedef SNC_structure_ SNC_structure;
  typedef typename SNC_structure::Infi_box                   Infi_box;
  typedef typename SNC_structure_::Sphere_kernel             Sphere_kernel;
  typedef typename Infi_box::Standard_kernel                 Standard_kernel;
  typedef typename Standard_kernel::Point_3                  Standard_point_3;
  typedef typename SNC_structure_::Kernel                    Kernel;
  typedef typename Kernel::RT                                RT;
  typedef typename Infi_box::NT                              NT;
  typedef CGAL::SNC_constructor_base<Items, SNC_structure>        Self;
  typedef CGAL::SNC_decorator<SNC_structure>                 SNC_decorator;
  typedef typename CGAL::SNC_const_decorator<SNC_structure>  SNC_const_decorator;
  typedef CGAL::SNC_intersection<SNC_structure>              SNC_intersection;

  typedef typename SNC_structure::Sphere_map             Sphere_map;
  typedef CGAL::SM_decorator<Sphere_map>                 SM_decorator;  
  typedef CGAL::SNC_SM_overlayer<Items, SM_decorator>    SM_overlayer;
  typedef CGAL::SM_const_decorator<Sphere_map>           SM_const_decorator;
  typedef CGAL::SM_point_locator<SM_decorator>           SM_point_locator;

  typedef typename SNC_structure::Vertex_iterator Vertex_iterator;
  typedef typename SNC_structure::Halfedge_iterator Halfedge_iterator;
  typedef typename SNC_structure::Halffacet_iterator Halffacet_iterator;

  typedef typename SNC_structure::Vertex_handle Vertex_handle;
  typedef typename SNC_structure::Halfedge_handle Halfedge_handle;
  typedef typename SNC_structure::Halffacet_handle Halffacet_handle;
  typedef typename SNC_structure::Volume_handle Volume_handle;

  typedef typename SNC_structure::Vertex_const_handle Vertex_const_handle;
  typedef typename SNC_structure::Halfedge_const_handle Halfedge_const_handle;
  typedef typename SNC_structure::Halffacet_const_handle Halffacet_const_handle;
  typedef typename SNC_structure::Halffacet_const_iterator Halffacet_const_iterator;

  typedef typename SNC_structure::SVertex_iterator SVertex_iterator;
  typedef typename SNC_structure::SHalfedge_iterator SHalfedge_iterator;
  typedef typename SNC_structure::SFace_iterator SFace_iterator;
  typedef typename SNC_structure::SHalfloop_iterator SHalfloop_iterator;

  typedef typename SNC_structure::SVertex_handle SVertex_handle;
  typedef typename SNC_structure::SHalfedge_handle SHalfedge_handle;
  typedef typename SNC_structure::SFace_handle SFace_handle;
  typedef typename SNC_structure::SHalfloop_handle SHalfloop_handle;

  typedef typename SNC_structure::SVertex_const_handle SVertex_const_handle; 
  typedef typename SNC_structure::SHalfedge_const_handle SHalfedge_const_handle; 
  typedef typename SNC_structure::SHalfloop_const_handle SHalfloop_const_handle; 
  typedef typename SNC_structure::SFace_const_handle SFace_const_handle; 

  typedef typename SNC_structure::SHalfedge_around_facet_circulator 
    SHalfedge_around_facet_circulator;
  typedef typename SNC_structure::SHalfedge_around_facet_const_circulator 
    SHalfedge_around_facet_const_circulator;
  typedef typename SNC_structure::SFace_cycle_iterator SFace_cycle_iterator;
  typedef typename SNC_structure::SFace_cycle_const_iterator SFace_cycle_const_iterator;
  typedef typename SNC_structure::Halffacet_cycle_iterator Halffacet_cycle_iterator;
  typedef typename SNC_structure::Halffacet_cycle_const_iterator Halffacet_cycle_const_iterator;

  typedef typename SNC_structure::Point_3 Point_3;
  typedef typename SNC_structure::Vector_3 Vector_3;
  typedef typename SNC_structure::Segment_3 Segment_3;
  typedef typename SNC_structure::Plane_3 Plane_3;
  typedef typename SNC_structure::Aff_transformation_3 Aff_transformation_3;

  typedef typename SNC_structure::Sphere_point Sphere_point;
  typedef typename SNC_structure::Sphere_segment Sphere_segment;
  typedef typename SNC_structure::Sphere_circle Sphere_circle;

  typedef typename SNC_structure::Mark Mark;

  typedef typename SM_decorator::SHalfedge_around_svertex_circulator 
                                 SHalfedge_around_svertex_circulator;
  typedef typename SM_decorator::SHalfedge_around_sface_circulator 
                                 SHalfedge_around_sface_circulator;
  typedef typename SM_const_decorator::SHalfedge_around_svertex_const_circulator 
                                       SHalfedge_around_svertex_const_circulator; 

  using SNC_decorator::is_standard;
  
  enum{NORMAL, CORNER, DEGENERATE};

  SNC_constructor_base( SNC_structure& W) : SNC_decorator(W) {}
  /*{\Mcreate makes |\Mvar| a decorator of |W|.}*/

 private:
  int compute_index(NT a, NT b) const {
    int i = 0;
    if(a > 0) ++i;
    if(b > 0) i+=2;
    return i;
  }

 public:

  Vertex_handle 
  create_extended_box_corner(NT x, NT y, NT z,
			     bool space, 
			     bool boundary
#ifdef CGAL_NEF_INDEXED_ITEMS
			     , int base
#endif
			     ) const { 

    CGAL_assertion(CGAL_NTS abs(x) == CGAL_NTS abs(y) &&
		   CGAL_NTS abs(y) == CGAL_NTS abs(z));
    
    CGAL_assertion(boundary == true);
    
    CGAL_NEF_TRACEN("  constructing box corner on "<<Point_3(x,y,z)<<"...");
    Point_3 p = Infi_box::create_extended_point(x,y,z); 
    Vertex_handle v = this->sncp()->new_vertex(p , boundary);
    CGAL_NEF_TRACEN( v->point());
    SM_decorator SD(&*v);
    Sphere_point sp[] = { Sphere_point(-x, 0, 0), 
			  Sphere_point(0, -y, 0), 
			  Sphere_point(0, 0, -z) };
    
  /* create box vertices */
    SVertex_handle sv[3];
    for(int vi=0; vi<3; ++vi) {
      sv[vi] = SD.new_svertex(sp[vi]);
      sv[vi]->mark() = boundary;
    }

#ifdef CGAL_NEF_INDEXED_ITEMS
    sv[0]->set_index(base+compute_index(y,z));
    sv[1]->set_index(base+4+compute_index(x,z));
    sv[2]->set_index(base+8+compute_index(x,y));
#endif

    /* create facet's edge uses */
    //  Sphere_segment ss[3];
    SHalfedge_handle she[3];
    for(int si=0; si<3; ++si)
      she[si] = SD.new_shalfedge_pair(sv[si], sv[(si+1)%3]);
    
    for(int i=0; i<3;++i) {
      she[i]->circle() = 
      Sphere_circle(Plane_3(sp[i],sp[(i+1)%3],Point_3(0,0,0)));
      she[i]->twin()->circle() =  she[i]->circle().opposite();
      she[i]->mark() = she[i]->twin()->mark() = boundary;
    }

#ifdef CGAL_NEF_INDEXED_ITEMS
    int bit = she[0]->circle().c()>0?0:1;
    int idx = base + (z>0?0:2);
    she[0]->set_index(idx+bit);
    she[0]->twin()->set_index(idx+1-bit);
    bit = she[1]->circle().a()>0?0:1;
    idx = base + 4 + (x>0?0:2);
    she[1]->set_index(idx+bit);
    she[1]->twin()->set_index(idx+1-bit);
    bit = she[2]->circle().b()>0?0:1;
    idx = base + 8 + (y>0?0:2);
    she[2]->set_index(idx+bit);
    she[2]->twin()->set_index(idx+1-bit);
#endif
    
    /* create facets */
    SFace_handle fi = SD.new_sface();
    SFace_handle fe = SD.new_sface();
    SD.link_as_face_cycle(she[0], fi);
    SD.link_as_face_cycle(she[0]->twin(), fe);
    
    Sphere_point p1 = she[0]->source()->point();
    Sphere_point p2 = she[0]->twin()->source()->point();
    Sphere_point p3 = she[0]->snext()->twin()->source()->point();
    if ( spherical_orientation(p1,p2,p3) > 0 ) {
      fi->mark() = space;
      fe->mark() = 0;
    }
    else {
      fi->mark() = 0;
      fe->mark() = space;
    }
    
    return v;
  }

  /*{\Mop produces the sphere map representing thp,e box corner in
          direction $(x,y,z)$.}*/

  template<typename Forward_iterator>  
  SHalfedge_handle add_outer_sedge_cycle(Vertex_handle v, 
					 Forward_iterator start,
					 Forward_iterator end, 
					 bool orient) {
    
    CGAL_assertion(start!=end);
    
    v->mark() = true;
    SM_decorator SD(&*v);
    
    Forward_iterator si;
    for(si=start; si!=end; ++si)
      SD.new_svertex(*si);
    
    SHalfedge_handle se,se_prev;
    std::list<SHalfedge_handle> se_list;
    SVertex_iterator sv(SD.svertices_begin()), sv_next(sv);
    for(;sv!=SD.svertices_end();++sv) {
      ++sv_next;
      sv->mark()=true;
      if(sv_next==SD.svertices_end()) sv_next=SD.svertices_begin();
      se=SD.new_shalfedge_pair(sv,sv_next);
      se_list.push_back(se);
      se->mark() = se->twin()->mark() = true;
     if(sv!=SD.svertices_begin()) {
	se->sprev() = se_prev;
	se_prev->snext() = se;
	se->twin()->snext() = se_prev->twin();
	se_prev->twin()->sprev() = se->twin();
      }
     se_prev = se;
    }
    
    se=*se_list.begin();
    se->sprev() = se_prev;
    se_prev->snext() = se;
    se->twin()->snext() = se_prev->twin();
    se_prev->twin()->sprev() = se->twin();

    typename std::list<SHalfedge_handle>::iterator seli;
    for(seli=se_list.begin();seli!=se_list.end();++seli) {
      se = *seli;
      se->circle() = Sphere_circle(se->source()->point(), se->snext()->source()->point());
      if(orient && seli==se_list.begin())
	 se->circle() = Sphere_circle(se->snext()->source()->point(), se->source()->point());
      se->circle() = normalized(se->circle());
      se->twin()->circle() = se->circle().opposite();
    }
    
    SFace_handle sfa = SD.new_sface();
    SFace_handle sfi = SD.new_sface();

    sfi->mark()=true;
    sfa->mark()=false;

    se=*se_list.begin();
    SD.link_as_face_cycle(se,sfi);
    SD.link_as_face_cycle(se->twin(),sfa);

    return se;
  }

  template<typename Forward_iterator>  
  SHalfedge_handle add_inner_sedge_cycle(Vertex_handle v, 
					 Forward_iterator start,
					 Forward_iterator end, 
					 bool orient, bool camera) {
    CGAL_assertion(start!=end);

    v->mark() = true;
    SM_decorator SD(&*v);

    Forward_iterator si;
    std::list<SVertex_handle> sv_list;
    for(si=start; si!=end; ++si)
      sv_list.push_back(SD.new_svertex(*si));
    
    SHalfedge_handle se, se_prev;
    typename std::list<SVertex_handle>::iterator 
      sv(sv_list.begin()), sv_next(sv);
    std::list<SHalfedge_handle> se_list;
    
    for(;sv!=sv_list.end();++sv) {
      ++sv_next;
      (*sv)->mark()=true;
      if(sv_next==sv_list.end()) sv_next=sv_list.begin();
      se=SD.new_shalfedge_pair(*sv,*sv_next);
      se_list.push_back(se);
      se->mark() = se->twin()->mark() = true;
      if(sv!=sv_list.begin()) {
	se->prev() = se_prev;
	se_prev->next() = se;
	se->twin()->next() = se_prev->twin();
	se_prev->twin()->prev() = se->twin();
      }
      se_prev = se;
    }

    se=*se_list.begin();
    se->prev() = se_prev;
    se_prev->next() = se;
    se->twin()->next() = se_prev->twin();
    se_prev->twin()->prev() = se->twin();

    typename std::list<SHalfedge_handle>::iterator seli;
    for(seli=se_list.begin();seli!=se_list.end();++seli) {
      se = *seli;
      se->circle() = Sphere_circle(se->source()->point(), se->snext()->source()->point());
      if(orient && seli==se_list.begin())
	 se->circle() = Sphere_circle(se->snext()->source()->point(), se->source()->point());
      se->circle() = normalized(se->circle());
      se->twin()->circle() = se->circle().opposite();
    }  

    SFace_handle sfa,sfi;
    if(camera) {
      sfa = ++SD.sfaces_begin();
      sfi = SD.new_sface();
      sfi->mark()=false;
    } else {
      sfa = SD.new_sface();
      sfi = SD.new_sface();
      sfa->mark()=true;
      sfi->mark()=false;
    }

    se=*se_list.begin();
    SD.link_as_face_cycle(se,sfi);
    SD.link_as_face_cycle(se->twin(),sfa);

    return se;
  }

  Vertex_handle create_for_infibox_overlay(Vertex_const_handle vin) const {

    Unique_hash_map<SHalfedge_handle, Mark> mark_of_right_sface;
    
    Vertex_handle v = this->sncp()->new_vertex(vin->point(), vin->mark());
    SM_decorator D(&*v);
    SM_const_decorator E(&*vin);
    
    SVertex_const_handle e;
    CGAL_forall_svertices(e, E)
      if(!Infi_box::is_edge_on_infibox(e))
	break;
    
    Sphere_point ps = e->point();
    ps = normalized(ps);
    SVertex_handle v1 = D.new_svertex(ps);
    SVertex_handle v2 = D.new_svertex(ps.antipode());
    CGAL_NEF_TRACEN("new svertex 1 " << ps);
    CGAL_NEF_TRACEN("new svertex 2 " << ps.antipode());
    v1->mark() = v2->mark() = e->mark();
    CGAL_NEF_TRACEN("svertex 1 " << ps << " has mark " << v1->mark());
    CGAL_NEF_TRACEN("svertex 2 " << ps.antipode() << " has mark " << v2->mark());
    
    if(E.is_isolated(e)) {
      CGAL_NEF_TRACEN("edge is isolated");
      SFace_handle f = D.new_sface();
      f->mark() = e->incident_sface()->mark();    
      D.link_as_isolated_vertex(v1, f);
      D.link_as_isolated_vertex(v2, f);
    }
    else {
      CGAL_NEF_TRACEN("edge is not isolated");
      SHalfedge_handle se;
      SFace_handle sf;
      
      SHalfedge_around_svertex_const_circulator ec(E.out_edges(e)), ee(ec);
      CGAL_For_all(ec,ee) {
	Sphere_segment seg(ec->source()->point(), 
			   ec->source()->point().antipode(), 
			   ec->circle());
	se = D.new_shalfedge_pair(v1, v2);
	se->mark() = se->twin()->mark() = ec->mark();
	mark_of_right_sface[se] = ec->incident_sface()->mark();
	se->circle() = ec->circle();
	se->twin()->circle() = se->circle().opposite();
      }
      
      SHalfedge_around_svertex_circulator ec2(D.out_edges(v1)), ee2(ec2);
      CGAL_For_all(ec2,ee2) {
	sf = D.new_sface();
	sf->mark() = mark_of_right_sface[ec2];
	D.link_as_face_cycle(SHalfedge_handle(ec2),sf);
      }   
    }
    
    return v;
  }

  Vertex_handle create_from_plane(const Plane_3& pl, const Point_3& p, 
				  const Mark& bnd, 
				  const Mark& in, const Mark& out) const {
    Vertex_handle v = this->sncp()->new_vertex( p, bnd);
    v->point() = p;
    Sphere_circle c(pl); // circle through origin parallel to h
    SM_decorator D(&*v);
    SHalfloop_handle l = D.new_shalfloop_pair();
    SFace_handle f1 = D.new_sface(), f2 = D.new_sface();
    D.link_as_loop(l,f1);
    D.link_as_loop(l->twin(),f2);
    
    l->circle() = c; 
    l->twin()->circle() = c.opposite();
    f1->mark() = out;
    f2->mark() = in;
    l->mark() = l->twin()->mark() = bnd;
    return v;
  }

  Vertex_handle create_from_point_on_infibox_facet(const Point_3& p) const {
    
    CGAL_assertion(!Infi_box::is_standard(p));
    if(Infi_box::x_on_box(p)) {
      if(p.hx() > 0) 
	return create_from_plane(Plane_3(1,0,0,1), p, 1, 1, 0);
      else
	return create_from_plane(Plane_3(-1,0,0,1), p, 1, 1, 0);
    }
    if(Infi_box::y_on_box(p)) {
      if(p.hy() > 0) 
	return create_from_plane(Plane_3(0,1,0,1), p, 1, 1, 0);
      else
	return create_from_plane(Plane_3(0,-1,0,1), p, 1, 1, 0);
    }
    if(Infi_box::z_on_box(p)) {
      if(p.hz() > 0) 
	return create_from_plane(Plane_3(0,0,1,1), p, 1, 1, 0);
      else
	return create_from_plane(Plane_3(0,0,-1,1), p, 1, 1, 0);
    }
    CGAL_error_msg( "something is wrong");
    return Vertex_handle();
  }
 
  Vertex_handle create_from_point_on_infibox_edge(const Point_3& p) const {

    Sphere_point sp;
    Sphere_circle c1,c2;
    
    if(!Infi_box::x_on_box(p)) {
      sp = Sphere_point(1,0,0);
      c1 = (p.hy() < RT(0)) ? Sphere_circle(0,-1,0) : Sphere_circle(0,1,0);
      c2 = (p.hz() < RT(0)) ? Sphere_circle(0,0,-1) : Sphere_circle(0,0,1);
    } else if(!Infi_box::y_on_box(p)) {
      sp = Sphere_point(0,1,0);
      c1 = (p.hx() < RT(0)) ? Sphere_circle(-1,0,0) : Sphere_circle(1,0,0);
      c2 = (p.hz() < RT(0)) ? Sphere_circle(0,0,-1) : Sphere_circle(0,0,1);
    } else if(!Infi_box::z_on_box(p)) {
      sp = Sphere_point(0,0,1);
      c1 = (p.hx() < RT(0)) ? Sphere_circle(-1,0,0) : Sphere_circle(1,0,0);
      c2 = (p.hy() < RT(0)) ? Sphere_circle(0,-1,0) : Sphere_circle(0,1,0);
    } else
      CGAL_error_msg( "line of code shall not be reached");
    
    Vertex_handle v = this->sncp()->new_vertex( p, true);
    SM_decorator D(&*v);
    SVertex_handle v1 = D.new_svertex(sp);
    SVertex_handle v2 = D.new_svertex(sp.antipode());
    SHalfedge_handle e1 = D.new_shalfedge_pair(v1,v2);
    SHalfedge_handle e2 = D.new_shalfedge_pair(v1,v2);
    SFace_handle f1 = D.new_sface();
    SFace_handle f2 = D.new_sface();
    v1->mark() = v2->mark() = 
      e1->mark() = e1->twin()->mark() =
      e2->mark() = e2->twin()->mark() = true;
    f1->mark() = false;
    f2->mark() = true;
    D.link_as_face_cycle(e1,f1);  
    D.link_as_face_cycle(e2,f2);
    e1->circle() = c1;
    e1->twin()->circle() = c1.opposite();
    e1->snext()->circle() = c2;
    e1->snext()->twin()->circle() = c2.opposite();
    
    return v;
  }
 
  Vertex_handle create_from_point_on_infibox_vertex(const Point_3& p) const {

    typedef typename Infi_box::Standard_point Standard_point;
    Standard_point ps(Infi_box::standard_point(p,1));
    return create_extended_box_corner(ps.hx(),
				      ps.hy(),
				      ps.hz(),true,true);
  }

  Vertex_handle create_from_facet(Halffacet_const_handle f,
				  const Point_3& p) {
    return create_from_plane(f->plane(), p, 
			     f->mark(), 
			     f->twin()->incident_volume()->mark(), 
			     f->incident_volume()->mark());
  }

  /*{\Mop produces the sphere map at point $p$ representing the local
     view of $f$. \precond $p$ is part of $f$.}*/

  Vertex_handle create_from_edge(Halfedge_const_handle e,
				 const Point_3& p) {
    //    typedef typename CGAL::SNC_const_decorator<SNC_structure> SNC_const_decorator;
    
    //    CGAL_assertion(SNC_const_decorator::segment(e).has_on(p));
    Vertex_handle v = this->sncp()->new_vertex( p, e->mark());
    SM_decorator D(&*v);
    SM_const_decorator E(&*e->source());
    Sphere_point ps = e->point();
    SVertex_handle v1 = D.new_svertex(ps);
    SVertex_handle v2 = D.new_svertex(ps.antipode());
    v1->mark() = v2->mark() = e->mark();

    bool first = true;
    
    // CGAL_NEF_SETDTHREAD(19*43*131);
    
    SHalfedge_const_handle ceee;
    CGAL_NEF_TRACEN("---------------------" << e->center_vertex()->point());
    CGAL_forall_shalfedges(ceee,E)
      CGAL_NEF_TRACEN("|" << ceee->circle() <<
	     "|" << ceee->mark() << 
	     " " << ceee->incident_sface()->mark());
    CGAL_NEF_TRACEN(" ");
    
    if(E.is_isolated(e)) {
      SFace_handle f = D.new_sface();
      D.link_as_isolated_vertex(v1,f);
      D.link_as_isolated_vertex(v2,f);
      f->mark() = e->incident_sface()->mark();
    }
    
    SHalfedge_around_svertex_const_circulator 
      ec1(e->out_sedge()), ee(ec1);
    SHalfedge_handle e1,e2;
    CGAL_For_all(ec1,ee) {
      if (first) e1 = D.new_shalfedge_pair(v1,v2);
      else       e1 = D.new_shalfedge_pair(e1, e2, SM_decorator::AFTER, 
					   SM_decorator::BEFORE);
      e2 = e1->twin(); 
      first = false;
    }
    
    SHalfedge_handle eee;
    CGAL_forall_shalfedges(eee,D)
      CGAL_NEF_TRACEN("|" << eee->circle());
    CGAL_NEF_TRACEN(" ");
       
    ec1 = e->out_sedge();
    SHalfedge_around_svertex_circulator ec2(v1->out_sedge());
    CGAL_For_all(ec1,ee) {
      CGAL_NEF_TRACEN("|" << ec1->circle() <<
	     "|" << ec1->mark() << 
	     " " << ec1->incident_sface()->mark());
      ec2->mark() = ec2->twin()->mark() = ec1->mark();
      ec2->circle() = ec1->circle();
      ec2->twin()->circle() = ec1->twin()->circle();
      SFace_handle f = D.new_sface();
      D.link_as_face_cycle(ec2,f);
      f->mark() = ec1->incident_sface()->mark();
      ++ec2;
    }
    
    CGAL_NEF_TRACEN(" ");
    CGAL_NEF_TRACEN("new vertex ");
    
    CGAL_forall_svertices(v1, D) {
      CGAL_NEF_TRACE("|" << v1->point() << "|" << v1->mark());
      CGAL_NEF_TRACEN(" ");
    }
    CGAL_NEF_TRACEN(" ");

    CGAL_forall_shalfedges(eee,D)
      CGAL_NEF_TRACEN("|" << eee->circle() <<
	     "|" << eee->mark() << 
	     " " << eee->incident_sface()->mark());
    CGAL_NEF_TRACEN("---------------------");
    
    return v;
  }

  /*{\Mop produces the sphere map at point $p$ representing the local
     view of $e$. \precond $p$ is part of $e$.}*/

  Vertex_handle clone_SM( Vertex_const_handle vin) {
    
#ifdef CGAL_NEF3_TIMER_SPHERE_SWEEPS
    ++number_of_clones;
#endif

    CGAL::Unique_hash_map<SVertex_const_handle, SVertex_handle>         VM;
    CGAL::Unique_hash_map<SHalfedge_const_handle, SHalfedge_handle>     EM;
    CGAL::Unique_hash_map<SHalfloop_const_handle, SHalfloop_handle>     LM;
    CGAL::Unique_hash_map<SFace_const_handle, SFace_handle>             FM;
    
    SM_const_decorator E(&*vin);
    Vertex_handle vout = this->sncp()->new_vertex(vin->point(), vin->mark());
    SM_decorator D(&*vout);
    
    SVertex_const_handle sv;
    CGAL_forall_svertices(sv, E) {
      VM[sv] = D.new_svertex(sv->point());
    }
    
    SHalfedge_const_handle se;
    CGAL_forall_sedges(se, E) {
      EM[se] = D.new_shalfedge_pair();
      EM[se->twin()] = EM[se]->twin();
    }
    
    SFace_const_handle sf;
    CGAL_forall_sfaces(sf, E)
      FM[sf] = D.new_sface();
    
    SHalfloop_handle sl;
    if(E.has_shalfloop()) {
      sl = LM[E.shalfloop()] = D.new_shalfloop_pair();
      LM[E.shalfloop()->twin()] = sl->twin();
      D.set_face(sl, FM[E.shalfloop()->incident_sface()]);
      D.set_face(sl->twin(), FM[E.shalfloop()->twin()->incident_sface()]);
      sl->circle() = E.shalfloop()->circle();
      sl->twin()->circle() = E.shalfloop()->twin()->circle();
      sl->mark() = sl->twin()->mark() = E.shalfloop()->mark();
    }
    
    CGAL_forall_svertices(sv, E) {
      D.set_first_out_edge(VM[sv], EM[E.first_out_edge(sv)]);
      D.set_face(VM[sv], FM[sv->incident_sface()]);
      VM[sv]->mark() = sv->mark();
    }
    
    CGAL_forall_shalfedges(se, E) {
      EM[se]->mark() = EM[se]->twin()->mark() = se->mark();
      D.set_source(EM[se], VM[se->source()]);
      D.set_prev(EM[se], EM[se->sprev()]); 
      D.set_next(EM[se], EM[se->snext()]);
      D.set_face(EM[se], FM[se->incident_sface()]);
      EM[se]->circle() = se->circle();
    }
    
    CGAL_forall_sfaces(sf, E) {
      FM[sf]->mark() = sf->mark();
      SFace_cycle_const_iterator sfc;
      for(sfc = sf->sface_cycles_begin(); sfc != sf->sface_cycles_end(); ++sfc) {
	if(sfc.is_svertex())
	  D.store_sm_boundary_object(VM[SVertex_const_handle(sfc)], FM[sf]);
	else if(sfc.is_shalfedge())
	  D.store_sm_boundary_object(EM[SHalfedge_const_handle(sfc)], FM[sf]);
	else if(sfc.is_shalfloop())
	  D.store_sm_boundary_object(LM[SHalfloop_const_handle(sfc)], FM[sf]);
	else CGAL_error_msg("damn wrong handle.");
      }
    }
    
    return vout;
  }

  template<typename Selection, typename Association>
  Sphere_map* create_edge_facet_overlay( Halfedge_const_handle e, 
					 Halffacet_const_handle f,
					 const Point_3& p,
					 const Selection& BOP, bool inv,
					 Association& ) {

#ifdef CGAL_NEF3_TIMER_SPHERE_SWEEPS
    ++number_of_edge_facet_overlays;
#endif

    CGAL_NEF_TRACEN("edge facet overlay " << p);

    Unique_hash_map<SHalfedge_handle, Mark> mark_of_right_sface;
    
    SM_decorator D(&*this->sncp()->new_vertex(p, BOP(e->mark(), f->mark())));
    SM_const_decorator E(&*e->source());

    Sphere_point ps = e->point();
    ps = normalized(ps);
    SVertex_handle v1 = D.new_svertex(ps);
    SVertex_handle v2 = D.new_svertex(ps.antipode());

    CGAL_NEF_TRACEN("new svertex 1 " << ps);
    CGAL_NEF_TRACEN("new svertex 2 " << ps.antipode());
    Halffacet_const_handle faces_p(f);
    Vector_3 vec(ps-CGAL::ORIGIN);
    if(faces_p->plane().oriented_side(p+vec) == ON_NEGATIVE_SIDE)
      faces_p = faces_p->twin();
    v1->mark() = BOP(e->mark(), faces_p->incident_volume()->mark(), inv);
    v2->mark() = BOP(e->mark(), faces_p->twin()->incident_volume()->mark(), inv);
    CGAL_NEF_TRACEN("svertex 1 " << ps << " has mark " << v1->mark());
    CGAL_NEF_TRACEN("svertex 2 " << ps.antipode() << " has mark " << v2->mark());

    if(E.is_isolated(e)) {
      CGAL_NEF_TRACEN("edge is isolated");
      Mark mf1 = BOP(e->incident_sface()->mark(), faces_p->incident_volume()->mark(), inv);
      Mark mf2 = BOP(e->incident_sface()->mark(), faces_p->twin()->incident_volume()->mark(), inv);
      Mark ml = BOP(e->incident_sface()->mark(), faces_p->mark(), inv);
      
      SFace_handle f1 = D.new_sface();
      D.link_as_isolated_vertex(v1, f1);
      f1->mark() = mf1;
      
      if(mf1 == mf2 && mf1 == ml) {
	D.link_as_isolated_vertex(v2, f1);
      }
      else {
	SHalfloop_handle l = D.new_shalfloop_pair();
	SFace_handle f2 = D.new_sface();    
	D.link_as_isolated_vertex(v2, f2);
	D.link_as_loop(l,f1);
	D.link_as_loop(l->twin(),f2);
	l->circle() = Sphere_circle(faces_p->plane()); 
	l->twin()->circle() = l->circle().opposite();
	f2->mark() = mf2;
	l->mark() = l->twin()->mark() = ml;
      }
    }
    else {
      CGAL_NEF_TRACEN("edge is not isolated");
      SVertex_handle sv;
      SHalfedge_handle se1;
      SHalfedge_handle se2;
      SFace_handle sf;
      Sphere_circle c(f->plane());
      
      SHalfedge_handle next_edge;
      SHalfedge_around_svertex_const_circulator ec(E.out_edges(e)), ee(ec);
      CGAL_For_all(ec,ee) {
	Sphere_segment seg(ec->source()->point(), 
			   ec->source()->point().antipode(), 
			   ec->circle());
	Sphere_point sp(intersection(c, seg.sphere_circle()));
	CGAL_NEF_TRACEN(seg <<" has_on " << sp);
	if(!seg.has_on(sp))
	  sp = sp.antipode();
	sv = D.new_svertex(sp);
	CGAL_NEF_TRACEN("new svertex 3 " << normalized(sp));
	sv->mark() = BOP(ec->mark(), f->mark(), inv);
	se1 = D.new_shalfedge_pair(v1, sv);
	if(next_edge == SHalfedge_handle())
	  se2 = D.new_shalfedge_pair(sv, v2); 
	else
	  se2 = D.new_shalfedge_pair(sv, next_edge, -1);
	next_edge = se2->twin();
	se1->mark() = se1->twin()->mark() = BOP(ec->mark(), faces_p->incident_volume()->mark(), inv);
	se2->mark() = se2->twin()->mark() = BOP(ec->mark(), faces_p->twin()->incident_volume()->mark(), inv);
	mark_of_right_sface[se1] = ec->incident_sface()->mark();
	se1->circle() = se2->circle() = ec->circle();
	se1->twin()->circle() = se2->twin()->circle() = se1->circle().opposite();
      }
      
      SHalfedge_around_svertex_circulator ec2(D.out_edges(v1)), ee2(ec2);
      CGAL_For_all(ec2,ee2) {
	SHalfedge_around_svertex_circulator en(ec2);
	++en;
	se1 = D.new_shalfedge_pair(ec2->twin(), en->twin(), -1, 1);
	CGAL_NEF_TRACEN("new edge pair " << ec2->twin()->source()->vector() << 
			" -> " << en->twin()->source()->vector());
	se1->circle() = Sphere_circle(faces_p->plane());
	se1->twin()->circle() = se1->circle().opposite();
        se1->mark() = se1->twin()->mark() = BOP(mark_of_right_sface[ec2], faces_p->mark(), inv);

	sf = D.new_sface();
	sf->mark() = BOP(mark_of_right_sface[ec2], faces_p->incident_volume()->mark(), inv);
	D.link_as_face_cycle(se1,sf);
	sf = D.new_sface();
	sf->mark() = BOP(mark_of_right_sface[ec2], faces_p->twin()->incident_volume()->mark(), inv);
	D.link_as_face_cycle(se1->twin(),sf);
      }   
    }
    
    return D.sphere_map();
  }

 public:
  Point_3 get_transformed_coords_of_vertex(const Point_3& p,
					   const std::list<Point_3>& segs1,
  					   const std::list<Point_3>& segs2) {
    CGAL_assertion(!segs1.empty() && !segs2.empty());
    
    int side_of_point=1;
    RT max = p.hx();
    if(CGAL_NTS abs(p.hy()) > CGAL_NTS abs(max)) {
      max = p.hy();
      side_of_point=2;
    }
    if(CGAL_NTS abs(p.hz()) > CGAL_NTS abs(max)) {
      max = p.hz();
      side_of_point=3;
    }
    if(max < RT(0)) 
      side_of_point = -side_of_point;
    
    typename std::list<Point_3>::const_iterator s1,s2,t1,t2;
    s1 = Infi_box::segment_on_side(side_of_point, segs1);
    t1 = s1;
    ++t1;
    if(t1 == segs1.end()) t1=segs1.begin();
    s2 = Infi_box::segment_on_side(side_of_point, segs2);
    t2 = s2;
    ++t2;
    if(t2 == segs2.end()) t2=segs2.begin();
    
    SNC_intersection is;
    Point_3 ip;
    bool flag=is.does_intersect_internally(Segment_3(*s1,*t1),Segment_3(*s2,*t2),ip);
    if(!flag) {
      if(*s1 == *s2) return normalized(*s1);
      else if(*s1 == *t2) return normalized(*s1);
      else if(*t1 == *s2) return normalized(*t1);
      else if(*t1 == *t2) return normalized(*t1);
    }
    return normalized(ip);
  }

  Point_3 transform_point_on_infibox(const Point_3& p, 
				     const Aff_transformation_3& aff) {

    Point_3 res(p.transform(aff));
    res = normalized(p);
    RT hw(CGAL_NTS abs(res.hx()));
    if(CGAL_NTS abs(res.hy()) > hw) hw = CGAL_NTS abs(res.hy());
    if(CGAL_NTS abs(res.hz()) > hw) hw = CGAL_NTS abs(res.hz());
    CGAL_assertion(hw.degree() == 1);
    CGAL_assertion(hw[0] == 0);
    return Point_3(res.hx(), res.hy(), res.hz(), hw(1));
  }

  bool erase_redundant_vertices() {
    
    std::list<Vertex_handle> redundant_points;
    Vertex_iterator v;
    CGAL_forall_vertices(v, *this->sncp()) 
      if(!is_standard(v))
	redundant_points.push_back(v);
    
    redundant_points.sort(Frame_point_lt<Infi_box,Vertex_handle>());
    
    typename std::list<Vertex_handle>::iterator vi, vinext;
    for(vi = redundant_points.begin(); vi != redundant_points.end(); ++vi) {
      vinext = vi;
      ++vinext;
      if(vinext == redundant_points.end()) break;
      if((*vi)->point() == (*vinext)->point()) {
	CGAL_assertion(!Infi_box::is_complex_facet_infibox_intersection(**vi));
	this->sncp()->delete_vertex(*vi);
      }
    }
    
    NT eval(Infi_box::compute_evaluation_constant_for_halfedge_pairup(*this->sncp()));
    CGAL_NEF_TRACEN("eval at " << eval);

    bool res;
    do {
    
      res = false;
      typedef Halfedge_key< Point_3, Halfedge_handle>
	Halfedge_key;
      typedef Halfedge_key_lt< Point_3, Halfedge_handle, SNC_decorator> 
	Halfedge_key_lt;
      typedef std::list<Halfedge_key>  Halfedge_list;
      
      typedef typename Standard_kernel::Kernel_tag Kernel_tag;
      typedef CGAL::Pluecker_line_3<Kernel_tag,Standard_kernel> Pluecker_line_3;
      typedef CGAL::Pluecker_line_lt        Pluecker_line_lt;
      typedef std::map< Pluecker_line_3, Halfedge_list, Pluecker_line_lt> 
	Pluecker_line_map;
      
      Unique_hash_map<Vertex_handle, bool> erase_vertex(false);
      std::list<Point_3> recreate;
      
      SNC_decorator D(*this);
      Pluecker_line_map M2;
      Pluecker_line_map M3;
      Pluecker_line_map M4;
      
      Halfedge_iterator e;
      CGAL_forall_halfedges(e,*this->sncp()) {
	Point_3 p = e->source()->point();
	Point_3 q = p + e->vector();
	Standard_point_3 sp = Infi_box::standard_point(p,eval);
	Standard_point_3 sq = Infi_box::standard_point(q,eval);
	Pluecker_line_3 l( sp, sq);
	
	int inverted;
	l = categorize( l, inverted);
	
	CGAL_NEF_TRACEN(" segment("<<p<<", "<<q<<")"<<
	       " direction("<<e->vector()<<")"<<
	       " line("<<l<<")"<<" inverted="<<inverted);
	
	if(Infi_box::is_edge_on_infibox(e)) {
	  if(Infi_box::is_type4(e))
	    M4[l].push_back(Halfedge_key(p,inverted,e));
	  else
	    if(Infi_box::is_type3(e))
	      M3[l].push_back(Halfedge_key(p,inverted,e));
	    else
	      M2[l].push_back(Halfedge_key(p,inverted,e));
	}
      }
      
      typename Pluecker_line_map::iterator it;
      
      CGAL_forall_iterators(it,M4) {
	CGAL_NEF_TRACEN("search opposite4  "<<it->first);
	it->second.sort(Halfedge_key_lt());
	typename Halfedge_list::iterator itl;
	CGAL_forall_iterators(itl,it->second) {
	  Halfedge_handle e1 = itl->e;
	  CGAL_NEF_TRACE("    " << e1->source()->point() << " -> ");
	  ++itl; 
	  if(itl == it->second.end()) {
	    erase_vertex[e1->source()] = true;
	    break;
	  }
	  Halfedge_handle e2 = itl->e;
	  CGAL_NEF_TRACE(e2->source()->point());
	  if(normalized(e1->point())!=normalized(e2->point().antipode())) {
	    erase_vertex[e1->source()] = true;
	    --itl;
	    CGAL_NEF_TRACE("   failed ");
	  }
	  CGAL_NEF_TRACEN("");
	}
	CGAL_NEF_TRACEN("");
	CGAL_NEF_TRACEN("");
      }
      
      CGAL_forall_iterators(it,M3) {
	CGAL_NEF_TRACEN("search opposite3  "<<it->first); 
	it->second.sort(Halfedge_key_lt());
	typename Halfedge_list::iterator itl;
	CGAL_forall_iterators(itl,it->second) {
	  Halfedge_handle e1 = itl->e;
	  CGAL_NEF_TRACE("    " << e1->source()->point() << " -> ");
	  ++itl; 
	  if(itl == it->second.end()) {
	    erase_vertex[e1->source()] = true;
	    break;
	  }
	  Halfedge_handle e2 = itl->e;
	  CGAL_NEF_TRACE(e2->source()->point());
	  if(normalized(e1->point())!=normalized(e2->point().antipode())) {
	    erase_vertex[e1->source()] = true;
	    --itl;
	    CGAL_NEF_TRACE("   failed ");
	  }
	  CGAL_NEF_TRACEN("");
	}
	CGAL_NEF_TRACEN("");
	CGAL_NEF_TRACEN("");
      }
      
      CGAL_forall_iterators(it,M2) {
	CGAL_NEF_TRACEN("search opposite2  "<<it->first); 
	it->second.sort(Halfedge_key_lt());
	typename Halfedge_list::iterator itl;
	CGAL_forall_iterators(itl,it->second) {
	  Halfedge_handle e1 = itl->e;
	  CGAL_NEF_TRACE("    " << e1->source()->point() << " -> ");
	  ++itl; 
	  if(itl == it->second.end()) {
	    erase_vertex[e1->source()] = true;
	    break;
	  }
	  Halfedge_handle e2 = itl->e;
	  CGAL_NEF_TRACE(e2->source()->point());
	  if(normalized(e1->point())!=normalized(e2->point().antipode())) {
	    erase_vertex[e1->source()] = true;
	    --itl;	    
	    CGAL_NEF_TRACE("   failed ");
	  }
	  CGAL_NEF_TRACEN("");
	}
	CGAL_NEF_TRACEN("");
	CGAL_NEF_TRACEN("");
      }
      
      std::vector<Vertex_iterator> vertices_to_delete ;
      Vertex_iterator v;
      CGAL_forall_vertices(v, *this->sncp()) {
	if(erase_vertex[v]) {
	  CGAL_NEF_TRACEN("erase " << v->point());
	  if(Infi_box::is_infibox_corner(v->point()))
	    recreate.push_back(v->point());
	  vertices_to_delete.push_back(v);
	  res = true;
	}
      }

      for( typename std::vector<Vertex_iterator>::iterator vit = vertices_to_delete.begin()
         ; vit != vertices_to_delete.end()
         ; ++ vit 
         )
        this->sncp()->delete_vertex(*vit);
         
      
      typename std::list<Point_3>::const_iterator pi;
      for(pi = recreate.begin(); pi != recreate.end(); ++pi)
	create_from_point_on_infibox_vertex(*pi);
      
    } while(res);
    
    return res;
  }

  std::list<Point_3> find_points_of_box_with_plane(const Plane_3& h) {
    Vector_3 orth = h.orthogonal_vector();
    
    NT orth_coords[3];
    orth_coords[0] = orth.hx()[0];
    orth_coords[1] = orth.hy()[0];
    orth_coords[2] = orth.hz()[0];

    int add_corners = 0;
    while(orth_coords[add_corners] == 0){
      CGAL_assertion(add_corners < 2);
      ++add_corners;
    }
    
    std::list<Point_3> points;
    for(int dir=0; dir<3;++dir) {

      NT cnst[3];
      for(int i=0; i<3;++i)
	cnst[i] = (i==dir? -h.d()[0] : 0);
      
      NT cross[4][4];
      cross[0][dir] = -orth_coords[(dir+1)%3]-orth_coords[(dir+2)%3];
      cross[1][dir] =  orth_coords[(dir+1)%3]-orth_coords[(dir+2)%3];
      cross[2][dir] =  orth_coords[(dir+1)%3]+orth_coords[(dir+2)%3];  
      cross[3][dir] = -orth_coords[(dir+1)%3]+orth_coords[(dir+2)%3];
  
      for(int i=0;i<4;++i)
	cross[i][3] = orth_coords[dir];

      cross[0][(dir+1)%3] = cross[3][(dir+1)%3] =  orth_coords[dir];
      cross[1][(dir+1)%3] = cross[2][(dir+1)%3] = -orth_coords[dir];
      
      cross[0][(dir+2)%3] = cross[1][(dir+2)%3] =  orth_coords[dir];
      cross[2][(dir+2)%3] = cross[3][(dir+2)%3] = -orth_coords[dir];

      for(int i=0; i<4; ++i)
	if(CGAL_NTS abs(RT(cnst[dir],cross[i][dir])) < 
	   CGAL_NTS abs(RT(0,orth_coords[dir])) ||
	   (CGAL_NTS abs(RT(cnst[dir],cross[i][dir])) == 
	    CGAL_NTS abs(RT(0,orth_coords[dir])) && 
	    dir == add_corners))

	  points.push_back(Point_3(RT(cnst[0], cross[i][0]),
				   RT(cnst[1], cross[i][1]),
				   RT(cnst[2], cross[i][2]),
				   RT(cross[i][3])));
    }

    for(int i=0;i<3;++i)
      orth_coords[i] = CGAL_NTS abs(orth_coords[i]);

    int max = 0;
    if(orth_coords[1] > orth_coords[0])
      max = 1;
    if(orth_coords[2] > orth_coords[max])
      max = 2;   

    int min = 0;
    if(orth_coords[1] < orth_coords[0])
      min = 1;
    if(orth_coords[2] < orth_coords[min])
      min = 2;   

    points.sort(circle_lt<Kernel>(max));

    typename std::list<Point_3>::const_iterator p;
    for(p=points.begin();p!=points.end();++p)
      CGAL_NEF_TRACEN(*p);
    
    return points;
  }  

  std::list<Point_3> find_facet_infibox_intersections(Halffacet_handle fi,
						      std::list<Point_3> points) {

    // points is a list of the required points, but with the inverse transformation applied

    std::list<Point_3> res;
    Halffacet_cycle_iterator fc = fi->facet_cycles_begin();
    CGAL_assertion(fc.is_shalfedge());
    SHalfedge_handle sh(fc);
    SHalfedge_around_facet_circulator fcc(sh), fend(fcc);
    CGAL_For_all(fcc,fend) {
      Point_3 src(fcc->source()->source()->point());
      Point_3 trg(fcc->source()->twin()->source()->point());
      typename std::list<Point_3>::iterator pi,pprev;
      pi = points.begin();
      while(pi != points.end()) {
	*pi = normalized(*pi);
	CGAL_NEF_TRACEN(src << "->" << trg << " has on " << *pi << src.x() << " : " );
	CGAL_NEF_TRACEN((src.x()-pi->x() <= 0) << "|" << (src.y()-pi->y() <= 0) << "|" << (src.z()-pi->z() <= 0) );
	CGAL_NEF_TRACEN((pi->x()-trg.x() <= 0) << "|" << (pi->y()-trg.y() <= 0) << "|" << (pi->z()-trg.z() <= 0) );
	if(( (src.x()-pi->x() <= 0 && pi->x()-trg.x() <= 0) || 
             (src.x()-pi->x() >= 0 && pi->x()-trg.x() >= 0)) &&
	   ((src.y()-pi->y() <= 0 && pi->y()-trg.y() <= 0) || 
	    (src.y()-pi->y() >= 0 && pi->y()-trg.y() >= 0)) &&
	   ((src.z()-pi->z() <= 0 && pi->z()-trg.z() <= 0) ||
	    (src.z()-pi->z() >= 0 && pi->z()-trg.z() >= 0))) {
	  pprev = pi;
	  ++pi;
	  res.push_back(*pprev);
	  points.erase(pprev);
	}
	else {
	  ++pi;
	}
      }
    }
    return res;
  }

  std::list<Vertex_handle> 
  create_vertices_on_infibox(const Plane_3& h, const std::list<Point_3> points, 
			     const Mark& bnd, const Mark& inside, const Mark& outside) {
    std::list<Vertex_handle> res;
    NT orth_coords[3];
    int min,max;
    Infi_box::compute_min_max(h,orth_coords,min,max);

    typename std::list<Point_3>::const_iterator p,prev,next;
    for(p=points.begin();p!=points.end();++p){
      
      if(p==points.begin()) prev = --points.end();
      else { prev = p; prev--;}
      if(p==--points.end()) next=points.begin();
      else {next = p; ++next;}
      CGAL_NEF_TRACEN("points " << *prev << "           " << *p << "      " << *next);

      Vector_3 v= *prev - *p;
      Sphere_point sp1(v);
      sp1 = normalized(sp1);
      CGAL_assertion(Infi_box::degree(sp1.hx()) == 0);
      CGAL_assertion(Infi_box::degree(sp1.hy()) == 0);
      CGAL_assertion(Infi_box::degree(sp1.hz()) == 0);
      CGAL_assertion(Infi_box::degree(sp1.hw()) == 0);
      
      v= *next - *p;
      Sphere_point sp2(v);
      sp2 = normalized(sp2);
      CGAL_assertion(Infi_box::degree(sp2.hx()) == 0);
      CGAL_assertion(Infi_box::degree(sp2.hy()) == 0);
      CGAL_assertion(Infi_box::degree(sp2.hz()) == 0);
      CGAL_assertion(Infi_box::degree(sp2.hw()) == 0);
      
      CGAL_NEF_TRACEN("sps " << sp1 << "     " << sp2);
      CGAL_NEF_TRACEN(orth_coords[min] << "|" << 
	     orth_coords[(min+1)%3] << "|" << 
	     orth_coords[(min+2)%3]);
      
      if(orth_coords[min]==0 && orth_coords[(min+1)%3] == 
	 orth_coords[(min+2)%3] && h.d() == 0) 
	res.push_back(create_degenerate_corner_frame_point(*p,sp1,sp2,min,max,h,bnd,inside,outside));
      else if(CGAL_NTS abs(p->hx()) == CGAL_NTS abs(p->hy()) && 
	      CGAL_NTS abs(p->hz()) == CGAL_NTS abs(p->hy()))
	res.push_back(create_corner_frame_point(*p,sp1,sp2,max,h,bnd,inside,outside));
      else
	res.push_back(create_frame_point(*p,sp1,sp2,h,bnd,inside,outside));
    }
    return res;
  }
    
  void create_vertices_for_halfspace
    (const Plane_3& h, const Mark& boundary) 
  {
    std::list<Point_3> points(find_points_of_box_with_plane(h));

    NT orth_coords[3];
    int min,max;
    Infi_box::compute_min_max(h,orth_coords,min,max);
    Mark inside = orth_coords[max] < 0;

#ifdef CGAL_NEF_INDEXED_ITEMS
    int base = Index_generator::get_unique_index();
    for(int i=0; i<3; ++i)
      Index_generator::get_unique_index();
    int plus = 0;
#endif

    typename std::list<Point_3>::const_iterator p,prev,next;
    for(p=points.begin();p!=points.end();++p){
      
      if(p==points.begin()) prev = --points.end();
      else { prev = p; prev--;}
      if(p==--points.end()) next=points.begin();
      else {next = p; ++next;}
      CGAL_NEF_TRACEN("points " << *prev << "           " << *p << "      " << *next);

      Vector_3 v= *prev - *p;
      Sphere_point sp1(v);
      sp1 = normalized(sp1);
      CGAL_assertion(Infi_box::degree(sp1.hx()) == 0);
      CGAL_assertion(Infi_box::degree(sp1.hy()) == 0);
      CGAL_assertion(Infi_box::degree(sp1.hz()) == 0);
      CGAL_assertion(Infi_box::degree(sp1.hw()) == 0);
      
      v= *next - *p;
      Sphere_point sp2(v);
      sp2 = normalized(sp2);
      CGAL_assertion(Infi_box::degree(sp2.hx()) == 0);
      CGAL_assertion(Infi_box::degree(sp2.hy()) == 0);
      CGAL_assertion(Infi_box::degree(sp2.hz()) == 0);
      CGAL_assertion(Infi_box::degree(sp2.hw()) == 0);
      
#ifdef CGAL_NEF_INDEXED_ITEMS
      create_halfspace_vertex(*p,sp1,sp2,h,boundary,
			      inside, !inside,
			      base, plus);
#else
      create_halfspace_vertex(*p,sp1,sp2,h,boundary,
			      inside, !inside);
#endif

    }
  }

  void create_vertices_of_box_with_plane(const Plane_3& h, bool b) {
    // CGAL_NEF_SETDTHREAD(19*43*11);
  
    std::list<Point_3> points(find_points_of_box_with_plane(h));
    create_vertices_on_infibox(h,points,b,true,false);

#ifdef CGAL_NEF_INDEXED_ITEMS
    int base = Index_generator::get_unique_index();
    for(int i=0; i<11; ++i)
      Index_generator::get_unique_index();
#endif

    RT sum= h.a()+h.b()+h.c(); 
    if(h.d()!=0 || sum!= 0) { 
      CGAL_NEF_TRACEN(sum); 
      create_extended_box_corner
	( 1, 1, 1, (sum<0 || (sum == 0 && h.d()<0)), true
#ifdef CGAL_NEF_INDEXED_ITEMS
	  , base
#endif
	  );
    }
    sum=-h.a()+h.b()+h.c(); 
    if(h.d()!=0 || sum!= 0) { 
      CGAL_NEF_TRACEN(sum); 
      create_extended_box_corner
	(-1, 1, 1, (sum<0 || (sum == 0 && h.d()<0)), true
#ifdef CGAL_NEF_INDEXED_ITEMS
	  , base
#endif
	  );
    }
    sum= h.a()-h.b()+h.c(); 
    if(h.d()!=0 || sum!= 0) { 
      CGAL_NEF_TRACEN(sum); 
      create_extended_box_corner
	( 1,-1, 1, (sum<0 || (sum == 0 && h.d()<0)), true
#ifdef CGAL_NEF_INDEXED_ITEMS
	  , base
#endif
	  );
    }
    sum=-h.a()-h.b()+h.c(); 
    if(h.d()!=0 || sum!= 0) { 
      CGAL_NEF_TRACEN(sum); 
      create_extended_box_corner
	(-1,-1, 1, (sum<0 || (sum == 0 && h.d()<0)), true
#ifdef CGAL_NEF_INDEXED_ITEMS
	  , base
#endif
	  );
    }
    sum= h.a()+h.b()-h.c(); 
    if(h.d()!=0 || sum!= 0) { 
      CGAL_NEF_TRACEN(sum); 
      create_extended_box_corner
	( 1, 1,-1, (sum<0 || (sum == 0 && h.d()<0)), true
#ifdef CGAL_NEF_INDEXED_ITEMS
	  , base
#endif
	  );
    }
    sum=-h.a()+h.b()-h.c(); 
    if(h.d()!=0 || sum!= 0) { 
      CGAL_NEF_TRACEN(sum); 
      create_extended_box_corner
	(-1, 1,-1, (sum<0 || (sum == 0 && h.d()<0)), true
#ifdef CGAL_NEF_INDEXED_ITEMS
	  , base
#endif
	  );
    }
    sum= h.a()-h.b()-h.c(); 
    if(h.d()!=0 || sum!= 0) { 
      CGAL_NEF_TRACEN(sum); 
      create_extended_box_corner
	( 1,-1,-1, (sum<0 || (sum == 0 && h.d()<0)), true
#ifdef CGAL_NEF_INDEXED_ITEMS
	  , base
#endif
	  );
    }
    sum=-h.a()-h.b()-h.c(); 
    if(h.d()!=0 || sum!= 0) { 
      CGAL_NEF_TRACEN(sum); 
      create_extended_box_corner
	(-1,-1,-1, (sum<0 || (sum == 0 && h.d()<0)), true
#ifdef CGAL_NEF_INDEXED_ITEMS
	  , base
#endif
	  );
    }
  }

  Vertex_handle create_halfspace_vertex
    (const Point_3& p, const Point_3& sp0, const Point_3& sp1,
     const Plane_3& h, const Mark& boundary, 
     const Mark& inside, const Mark& outside
#ifdef CGAL_NEF_INDEXED_ITEMS
     , int base, int& plus
#endif
     ) const
  {
    // Question: why can this be a const fucnction??

    if(h.d() == 0) {
      CGAL_assertion(CGAL_NTS abs(p.hy()) != 
		     CGAL_NTS abs(p.hx()) ||
		     CGAL_NTS abs(p.hz()) != 
		     CGAL_NTS abs(p.hx()));
    }
    
    Vertex_handle v=
      this->sncp()->new_vertex(normalized(p), boundary);
    SM_decorator SD(&*v); 

    std::vector<SVertex_handle> sv(2);
    sv[0] = SD.new_svertex(sp0);
    sv[1] = SD.new_svertex(sp1);
    sv[0]->mark() = sv[1]->mark() = boundary;
    
    SHalfedge_handle se =
      SD.new_shalfedge_pair(sv[0], sv[1]);
    se->circle() = 
      Sphere_circle(Plane_3(sp0, sp1, Point_3(0,0,0)));
    se->twin()->circle() = se->circle().opposite();
    se->mark() = se->twin()->mark() = boundary;
    
    SHalfedge_handle se2 =
      SD.new_shalfedge_pair(sv[0], sv[1]);
    se2->circle() = se->twin()->circle();
    se2->twin()->circle() = se->circle();
    se2->mark() = se2->twin()->mark() = boundary;

#ifdef CGAL_NEF_INDEXED_ITEMS
    sv[0]->set_index(base+plus);
    plus = (plus+1) % 4;
    sv[1]->set_index(base+plus);
    se->set_index(base);
    se->twin()->set_index(base+1);    
#endif

    SFace_handle sf0 = SD.new_sface();
    SFace_handle sf1 = SD.new_sface();
    sf0->mark()= inside;
    sf1->mark()= outside;
    
    SD.link_as_face_cycle(se,sf0);
    SD.link_as_face_cycle(se->twin(),sf1);
    
    return v;
  }

  Vertex_handle 
    create_frame_point(Point_3 p, Point_3 sp1, Point_3 sp2, Plane_3 h, 
		       const Mark& boundary, const Mark& inside, const Mark& outside) const 
  {
    if(h.d() == 0) {
      CGAL_assertion(CGAL_NTS abs(p.hy()) != CGAL_NTS abs(p.hx()) ||
		     CGAL_NTS abs(p.hz()) != CGAL_NTS abs(p.hx()));
    }
    
    int max = 0;
    if(CGAL_NTS abs(p.hx()) > CGAL_NTS abs(p.hy()))
      max = 1;
    if(CGAL_NTS abs(p.hx()) > CGAL_NTS abs(p.hz()))
      max = 2;
    
    CGAL_NEF_TRACEN("create frame point ");
    
    CGAL_NEF_TRACEN("create spoints");
    Sphere_point SP[4];
    switch(max) {
    case 0: SP[2] = Sphere_point(1,0,0); break;
    case 1: SP[2] = Sphere_point(0,1,0); break;
    case 2: SP[2] = Sphere_point(0,0,1); break;
    default: CGAL_error_msg("wrong value");
    }
    
    SP[1]=sp1;
    SP[0]=sp2;
    
    if (spherical_orientation(SP[0],SP[1],SP[2]) < 0) {
      SP[3] = SP[2];
      SP[2] = CGAL::ORIGIN-SP[3];
    }
    else
      SP[3] = CGAL::ORIGIN-SP[2];
    
    RT delta = h.a()*SP[2].hx()+h.b()*SP[2].hy()+h.c()*SP[2].hz();
    CGAL_assertion(delta !=0);
    Mark swtch = (delta <  0);
    Mark fmark0 = (swtch && inside)  || (!swtch && !inside);
    Mark fmark1 = (swtch && outside) || (!swtch && !outside);
    
    return create_SM_on_infibox(p, SP, 4, boundary, fmark0, fmark1);
  }    
  
  Vertex_handle 
    create_corner_frame_point(Point_3 p, Point_3 sp1, Point_3 sp2, 
			      int max, Plane_3 h, 
			      const Mark& boundary, const Mark& inside, const Mark& outside) const {
    CGAL_assertion(h.d() == 0);
    
    Vector_3 vec = h.orthogonal_vector();
    
    CGAL_assertion(CGAL_NTS abs(vec.hx()) != CGAL_NTS abs(vec.hy()) &&
		   CGAL_NTS abs(vec.hy()) != CGAL_NTS abs(vec.hz()) && 
		   CGAL_NTS abs(vec.hx()) != CGAL_NTS abs(vec.hz()));
    
    CGAL_assertion(vec.hx() + vec.hy() == vec.hz() ||
		   vec.hx() + vec.hz() == vec.hy() ||
		   vec.hy() + vec.hz() == vec.hx());
    
    CGAL_NEF_TRACEN("create corner frame point ");
    
    RT vp[3];
    vp[0] = -p.hx()[1];
    vp[1] = -p.hy()[1];
    vp[2] = -p.hz()[1];
    
    CGAL_NEF_TRACEN("create spoints");
    Sphere_point SP[5];
    switch(max) {
    case 0: 
      SP[3] = Sphere_point(0,vp[1],0); 
      SP[2]= Sphere_point(0,0,vp[2]); 
      SP[4] = Sphere_point(vp[0],0,0); 
      break;
    case 1: 
      SP[3] = Sphere_point(vp[0],0,0); 
      SP[2]= Sphere_point(0,0,vp[2]); 
      SP[4] = Sphere_point(0,vp[1],0); 
      break;
    case 2: 
      SP[3] = Sphere_point(vp[0],0,0); 
      SP[2]= Sphere_point(0,vp[1],0); 
      SP[4] = Sphere_point(0,0,vp[2]); 
      break;
    default: CGAL_error_msg("wrong value");
    }
    
    if (spherical_orientation(SP[3],Sphere_point(sp1),Sphere_point(sp2)) > 0) {
      SP[0] = sp1;
      SP[1] = sp2;
    }
    else {
      SP[0] = sp2;
      SP[1] = sp1;
    }
    
    if (spherical_orientation(SP[2],SP[3],SP[0]) < 0) {
      Sphere_point sx = SP[2];
      SP[2] = SP[3];
      SP[3] = sx;
    }
    
    RT delta = h.a()*SP[4].hx()+h.b()*SP[4].hy()+h.c()*SP[4].hz();
    CGAL_assertion(delta !=0);
    Mark swtch = (delta >  0);
    Mark fmark0 = (swtch && inside)  || (!swtch && !inside);
    Mark fmark1 = (swtch && outside) || (!swtch && !outside);
    
    return create_SM_on_infibox(p, SP, 5, boundary, fmark0, fmark1);
  }


  Vertex_handle 
    create_degenerate_corner_frame_point(Point_3 p, Point_3 sp1,Point_3 sp2, 
					 int min,int max, Plane_3 h,
					 const Mark& boundary, 
					 const Mark& inside, const Mark& outside)const {
    CGAL_assertion(h.d() == 0);
    
    Vector_3 vec = h.orthogonal_vector();
    
    CGAL_assertion(
       (CGAL_NTS abs(vec.hx()) == CGAL_NTS abs(vec.hy()) && vec.hz() == 0) ||
       (CGAL_NTS abs(vec.hy()) == CGAL_NTS abs(vec.hz()) && vec.hx() == 0) || 
       (CGAL_NTS abs(vec.hx()) == CGAL_NTS abs(vec.hz()) && vec.hy() == 0));
    
    RT vp[3];
    vp[0] = -p.hx()[1];
    vp[1] = -p.hy()[1];
    vp[2] = -p.hz()[1];
    
    CGAL_NEF_TRACEN("create degenerate corner frame point ");
    
    CGAL_NEF_TRACEN("create spoints");
    Sphere_point SP[4];
    
    switch(max) { 
    case 0: SP[2] = Sphere_point(vp[0],0,0); break; // plane(x,x,0), plane(x,0,x)
    case 1: SP[2] = Sphere_point(0,vp[1],0); break; // plane(0,x,x)
    default: CGAL_error_msg("wrong value \"max\"");
    }
    
    switch(min+max) {
    case 1: SP[3] = Sphere_point(0,0,vp[2]); break; // plane(0,x,x), plane(x,0,x)
    case 2: SP[3] = Sphere_point(0,vp[1],0); break; // plane(x,x,0)
    default: CGAL_error_msg("wrong value \"min+max\"");
    }
    
    if (spherical_orientation(SP[2],Sphere_point(sp1),Sphere_point(sp2)) > 0) {
      SP[0] = sp1;
      SP[1] = sp2;
    }
    else {
      SP[0] = sp2;
      SP[1] = sp1;
    }
    
    RT delta = h.a()*SP[2].hx()+h.b()*SP[2].hy()+h.c()*SP[2].hz();
    CGAL_assertion(delta !=0);
    Mark swtch = (delta <  0);
    Mark fmark0 = (swtch && inside)  || (!swtch && !inside);
    Mark fmark1 = (swtch && outside) || (!swtch && !outside);

    return create_SM_on_infibox(p, SP, 4, boundary, fmark0, fmark1);
  }   



  Vertex_handle 
    create_SM_on_infibox(const Point_3& center, Sphere_point* SP, 
			 int size, const Mark& boundary, 
			 const Mark& fmark0, const Mark& fmark1) const {
    Vertex_handle v=this->sncp()->new_vertex(normalized(center), boundary);
    SM_decorator SD(&*v); 
    
    CGAL_NEF_TRACEN("create_SM_on_infibox: center = " << center);
    
    CGAL_NEF_TRACEN("create svertices");
    std::vector<SVertex_handle> sv(size);
    CGAL_NEF_TRACEN("create_SM_on_infibox: size = " << size);
    for(int i=0; i<size; ++i) {
      CGAL_NEF_TRACEN("                      SP["<< i << "]=" << SP[i]);
      sv[i] = SD.new_svertex(SP[i]);
      sv[i]->mark() = (i < 2 ? boundary : 1);
    }
    
    CGAL_NEF_TRACEN("create sedges");
    std::vector<SHalfedge_handle> she(size+1);
    for(int si=0; si<size-1;++si) {
      she[si]=SD.new_shalfedge_pair(sv[si], sv[(si+1)%(size-1)]);
      she[si]->circle()= 
	Sphere_circle(Plane_3(SP[si],SP[(si+1)%(size-1)],Point_3(0,0,0)));
      she[si]->twin()->circle() = she[si]->circle().opposite();
      she[si]->mark() = she[si]->twin()->mark() = (si == 0 ? boundary : 1);
    }
    
    she[size-1] = SD.new_shalfedge_pair(she[0], sv[size-1], -1);
    CGAL_assertion(she[size-1]->snext() == she[size-1]->twin());
    CGAL_assertion(she[size-1]->twin()->sprev() == she[size-1]);
    CGAL_assertion(she[size-1]->twin()->sprev()->twin() == she[size-1]->twin());
    
    she[size]   = SD.new_shalfedge_pair(she[size-1]->twin(), she[0]->twin(), 1, 1);
    
    CGAL_assertion(she[0]->snext() == she[1]);
    CGAL_assertion(she[1]->snext() == she[2]);
    if(size == 4)
      CGAL_assertion(she[2]->snext() == she[0]);  
    else {
      CGAL_assertion(she[2]->snext() == she[3]);
      CGAL_assertion(she[3]->snext() == she[0]);
    }
    
    CGAL_assertion(she[0]->twin()->snext() == she[size-1]);
    CGAL_assertion(she[size-1]->snext() == she[size]);
    CGAL_assertion(she[size]->snext() == she[0]->twin());
    
    CGAL_assertion(she[4]->twin()->snext() == she[3]->twin());
    CGAL_assertion(she[3]->twin()->snext() == she[2]->twin());
    CGAL_assertion(she[2]->twin()->snext() == she[1]->twin());
    if(size == 4)
      CGAL_assertion(she[1]->twin()->snext() == she[4]->twin());
    else {
      CGAL_assertion(she[5]->twin()->snext() == she[4]->twin());
      CGAL_assertion(she[1]->twin()->snext() == she[5]->twin());
    }
    
    she[size-1]->circle()= Sphere_circle(Plane_3(SP[0],SP[size-1],Point_3(0,0,0)));
    she[size-1]->twin()->circle() =  she[size-1]->circle().opposite();	  
    she[size-1]->mark() = she[size-1]->twin()->mark() = 1;
    
    she[size]->circle()= Sphere_circle(Plane_3(SP[size-1],SP[1],Point_3(0,0,0)));
    she[size]->twin()->circle() =  she[size]->circle().opposite();
    she[size]->mark() = she[size]->twin()->mark() = 1;
    
    CGAL_NEF_TRACEN("create sfaces");
    SFace_handle sf[3];
    for(int i=0; i<3; ++i)
      sf[i] = SD.new_sface();
    sf[0]->mark()= fmark0;
    sf[1]->mark()= fmark1;
    sf[2]->mark()=0;
    
    SD.link_as_face_cycle(she[0],sf[0]);
    SD.link_as_face_cycle(she[0]->twin(),sf[1]);
    SD.link_as_face_cycle(she[1]->twin(),sf[2]);
    
    return v;
  }

  void correct_infibox_sface_marks() {
    Vertex_iterator v;
    CGAL_forall_vertices(v,*this->sncp()) {
      if(is_standard(v) || Infi_box::is_infibox_corner(v->point())) continue;
      SM_decorator SD(&*v);
      SFace_iterator sf;
      CGAL_forall_sfaces(sf, SD)
	sf->volume()->mark() = sf->mark();
    }
    CGAL_forall_vertices(v,*this->sncp()) {
      if(is_standard(v) || !Infi_box::is_infibox_corner(v->point())) continue;
      SM_decorator SD(&*v);
      SFace_iterator sf;
      CGAL_forall_sfaces(sf, SD)
	sf->mark() = sf->volume()->mark();
    }    
  }

  void create_facets_and_volumes_of_halfspace(const Plane_3& h) {
    SHalfedge_handle se =
      this->sncp()->vertices_begin()->shalfedges_begin();
    Halffacet_handle f =
      this->sncp()->new_halffacet_pair(h);
    Volume_handle c0 = 
      this->sncp()->new_volume();
    Volume_handle c1 = 
      this->sncp()->new_volume();

    link_as_facet_cycle(se, f->twin());
    link_as_facet_cycle(se->twin(), f);
    f->mark() = f->twin()->mark() = se->mark();
    store_boundary_object(se->incident_sface(), c0);
    store_boundary_object(se->twin()->incident_sface(), c1);
    c0->mark() = se->twin()->incident_sface()->mark();
    c1->mark() = se->incident_sface()->mark();
    f->incident_volume() = c1;
    f->twin()->incident_volume() = c0;

    SHalfedge_around_facet_circulator
      sfc(se), send(sfc);
    CGAL_For_all(sfc, send) {
      sfc->incident_sface()->volume() = c1;
      sfc->twin()->incident_sface()->volume() = c0;
    }

    SNC_io_parser<SNC_structure>::dump(*this->sncp(), std::cerr);
  }

  void correct_infibox_sedge_marks() {
    SHalfedge_iterator se;
    Vertex_iterator v;
    CGAL_forall_vertices(v,*this->sncp()) {
      if(is_standard(v)) continue;
      SM_decorator SD(&*v);
      CGAL_forall_shalfedges(se,SD)
	if(Infi_box::is_sedge_on_infibox(se)) {
	  se->mark() = true;
	  se->source()->mark() = true;
	}
    }
  }

  void assign_indices() {}

}; // SNC_constructor_base<SNC>

template <typename Items, typename SNC_structure_>
class SNC_constructor : public SNC_constructor_base<Items, SNC_structure_>
{
public:
  SNC_constructor( SNC_structure_& W)
    : SNC_constructor_base<Items, SNC_structure_>(W)
  {}

};

  

template<typename SNC_structure_>
class SNC_constructor<SNC_indexed_items, SNC_structure_>
  : public SNC_constructor_base<int, SNC_structure_> {

    typedef SNC_structure_                               SNC_structure;
  typedef SNC_constructor_base<int, SNC_structure>       Base;
  typedef typename SNC_structure::Sphere_map             Sphere_map;
  typedef CGAL::SM_const_decorator<Sphere_map>           SM_const_decorator;
  typedef CGAL::SM_decorator<Sphere_map>                 SM_decorator;

  public:  
  typedef typename SNC_structure::SVertex_iterator SVertex_iterator;
  typedef typename SNC_structure::Halfedge_iterator Halfedge_iterator;
  typedef typename SNC_structure::Halffacet_iterator Halffacet_iterator;

  typedef typename SNC_structure::Vertex_handle Vertex_handle;
  typedef typename SNC_structure::Halfedge_handle Halfedge_handle;
  typedef typename SNC_structure::Halffacet_handle Halffacet_handle;

  typedef typename SNC_structure::Vertex_const_handle Vertex_const_handle;
  typedef typename SNC_structure::Halfedge_const_handle Halfedge_const_handle;
  typedef typename SNC_structure::Halffacet_const_handle Halffacet_const_handle;
  typedef typename SNC_structure::SVertex_const_handle SVertex_const_handle; 
  typedef typename SNC_structure::SHalfedge_const_handle SHalfedge_const_handle; 
  typedef typename SNC_structure::SHalfloop_const_handle SHalfloop_const_handle; 
  typedef typename SNC_structure::SFace_const_handle SFace_const_handle; 

  typedef typename SNC_structure::SVertex_handle SVertex_handle;
  typedef typename SNC_structure::SHalfedge_handle SHalfedge_handle;
  typedef typename SNC_structure::SFace_handle SFace_handle;
  typedef typename SNC_structure::SHalfloop_handle SHalfloop_handle;

  typedef typename SNC_structure::SHalfedge_around_facet_circulator 
    SHalfedge_around_facet_circulator;
  typedef typename SNC_structure::SHalfedge_around_facet_const_circulator 
    SHalfedge_around_facet_const_circulator;
  typedef typename SNC_structure::SHalfedge_around_svertex_circulator 
    SHalfedge_around_svertex_circulator;
  typedef typename SNC_structure::SHalfedge_around_svertex_const_circulator 
    SHalfedge_around_svertex_const_circulator;
  typedef typename SNC_structure::SFace_cycle_const_iterator SFace_cycle_const_iterator;
  typedef typename SNC_structure::Halffacet_cycle_iterator Halffacet_cycle_iterator;
  typedef typename SNC_structure::Halffacet_cycle_const_iterator Halffacet_cycle_const_iterator;

  typedef typename SNC_structure::Point_3   Point_3;
  typedef typename SNC_structure::Vector_3  Vector_3;

  typedef typename SNC_structure::Sphere_point Sphere_point;
  typedef typename SNC_structure::Sphere_segment Sphere_segment;
  typedef typename SNC_structure::Sphere_circle Sphere_circle;

  typedef typename SNC_structure::Mark      Mark;

  using Base::create_from_plane;

  public:
  SNC_constructor( SNC_structure& W) : Base(W) {}

  Vertex_handle create_from_facet(Halffacet_const_handle f,
				  const Point_3& p) {

    Vertex_handle v = create_from_plane(f->plane(), p, 
					f->mark(), 
					f->twin()->incident_volume()->mark(), 
					f->incident_volume()->mark());
    
    v->shalfloop()->set_index_facet(f->twin());
    v->shalfloop()->twin()->set_index_facet(f);
    
    SHalfedge_const_handle se = f->twin()->facet_cycles_begin();
    CGAL_assertion(v->shalfloop()->circle() == se->circle());
    v->shalfloop()->set_index(se->get_index());
    v->shalfloop()->twin()->set_index(se->twin()->get_index());
    return v;
  }


  Vertex_handle create_from_edge(Halfedge_const_handle e,
				 const Point_3& p) {

    Vertex_handle v = Base::create_from_edge(e,p);
    SVertex_iterator sv = v->svertices_begin();

    SHalfedge_around_svertex_const_circulator
      ec1(e->out_sedge()), ee(ec1);
    SHalfedge_around_svertex_circulator
      ec2(sv->out_sedge());
    CGAL_For_all(ec1,ee) {
      ec2->set_index_facet(ec1->facet());
      ec2->twin()->set_index_facet(ec1->twin()->facet());
      ec2->set_index(ec1->get_index());
      ec2->twin()->set_index(ec1->twin()->get_index());
      ++ec2;
    }

    sv->set_index(e->get_index());
    ++sv;
    sv->set_index(e->get_index());
    
    return v;
  }

  Vertex_handle clone_SM( Vertex_const_handle vin) {
    
#ifdef CGAL_NEF3_TIMER_SPHERE_SWEEPS
    ++number_of_clones;
#endif

    CGAL::Unique_hash_map<SVertex_const_handle, SVertex_handle>         VM;
    CGAL::Unique_hash_map<SHalfedge_const_handle, SHalfedge_handle>     EM;
    CGAL::Unique_hash_map<SHalfloop_const_handle, SHalfloop_handle>     LM;
    CGAL::Unique_hash_map<SFace_const_handle, SFace_handle>             FM;
    
    SM_const_decorator E(&*vin);
    Vertex_handle vout = this->sncp()->new_vertex(vin->point(), vin->mark());
    SM_decorator D(&*vout);
    
    SVertex_const_handle sv;
    CGAL_forall_svertices(sv, E) {
      VM[sv] = D.new_svertex(sv->point());
      VM[sv]->set_index(sv->get_index());
    }
    
    SHalfedge_const_handle se;
    CGAL_forall_sedges(se, E) {
      EM[se] = D.new_shalfedge_pair();
      EM[se->twin()] = EM[se]->twin();
      EM[se]->set_index(se->get_index());
      EM[se->twin()]->set_index(se->twin()->get_index());
    }
    
    SFace_const_handle sf;
    CGAL_forall_sfaces(sf, E)
      FM[sf] = D.new_sface();
    
    SHalfloop_handle sl;
    if(E.has_shalfloop()) {
      sl = LM[E.shalfloop()] = D.new_shalfloop_pair();
      LM[E.shalfloop()->twin()] = sl->twin();
      D.set_face(sl, FM[E.shalfloop()->incident_sface()]);
      D.set_face(sl->twin(), FM[E.shalfloop()->twin()->incident_sface()]);
      sl->circle() = E.shalfloop()->circle();
      sl->twin()->circle() = E.shalfloop()->twin()->circle();
      sl->mark() = sl->twin()->mark() = E.shalfloop()->mark();
      sl->set_index(E.shalfloop()->get_index());
      sl->twin()->set_index(E.shalfloop()->twin()->get_index());
    }
    
    CGAL_forall_svertices(sv, E) {
      D.set_first_out_edge(VM[sv], EM[E.first_out_edge(sv)]);
      D.set_face(VM[sv], FM[sv->incident_sface()]);
      VM[sv]->mark() = sv->mark();
    }
    
    CGAL_forall_shalfedges(se, E) {
      EM[se]->mark() = EM[se]->twin()->mark() = se->mark();
      D.set_source(EM[se], VM[se->source()]);
      D.set_prev(EM[se], EM[se->sprev()]); 
      D.set_next(EM[se], EM[se->snext()]);
      D.set_face(EM[se], FM[se->incident_sface()]);
      EM[se]->circle() = se->circle();
    }
    
    CGAL_forall_sfaces(sf, E) {
      FM[sf]->mark() = sf->mark();
      SFace_cycle_const_iterator sfc;
      for(sfc = sf->sface_cycles_begin(); sfc != sf->sface_cycles_end(); ++sfc) {
	if(sfc.is_svertex())
	  D.store_sm_boundary_object(VM[SVertex_const_handle(sfc)], FM[sf]);
	else if(sfc.is_shalfedge())
	  D.store_sm_boundary_object(EM[SHalfedge_const_handle(sfc)], FM[sf]);
	else if(sfc.is_shalfloop())
	  D.store_sm_boundary_object(LM[SHalfloop_const_handle(sfc)], FM[sf]);
	else CGAL_error_msg("damn wrong handle.");
      }
    }
    
    return vout;
  }

  template<typename Selection, typename Association>
  Sphere_map* create_edge_facet_overlay( Halfedge_const_handle e, 
					 Halffacet_const_handle f,
					 const Point_3& p,
					 const Selection& BOP, bool inv,
					 Association& A) {

#ifdef CGAL_NEF3_TIMER_SPHERE_SWEEPS
    ++number_of_edge_facet_overlays;
#endif

    CGAL_NEF_TRACEN("edge facet overlay " << p);

    Unique_hash_map<SHalfedge_handle, Mark> mark_of_right_sface;
    
    SM_decorator D(&*this->sncp()->new_vertex(p, BOP(e->mark(), f->mark())));
    SM_const_decorator E(&*e->source());

    Sphere_point ps = e->point();
    ps = normalized(ps);
    SVertex_handle v1 = D.new_svertex(ps);
    SVertex_handle v2 = D.new_svertex(ps.antipode());
    v1->set_index(e->get_index());
    v2->set_index(e->get_index());

    CGAL_NEF_TRACEN("new svertex 1 " << ps);
    CGAL_NEF_TRACEN("new svertex 2 " << ps.antipode());
    Halffacet_const_handle faces_p(f);
    Vector_3 vec(ps-CGAL::ORIGIN);
    if(faces_p->plane().oriented_side(p+vec) == ON_NEGATIVE_SIDE)
      faces_p = faces_p->twin();
    v1->mark() = BOP(e->mark(), faces_p->incident_volume()->mark(), inv);
    v2->mark() = BOP(e->mark(), faces_p->twin()->incident_volume()->mark(), inv);
    CGAL_NEF_TRACEN("svertex 1 " << ps << " has mark " << v1->mark());
    CGAL_NEF_TRACEN("svertex 2 " << ps.antipode() << " has mark " << v2->mark());

    if(E.is_isolated(e)) {
      CGAL_NEF_TRACEN("edge is isolated");
      Mark mf1 = BOP(e->incident_sface()->mark(), faces_p->incident_volume()->mark(), inv);
      Mark mf2 = BOP(e->incident_sface()->mark(), faces_p->twin()->incident_volume()->mark(), inv);
      Mark ml = BOP(e->incident_sface()->mark(), faces_p->mark(), inv);
      
      SFace_handle f1 = D.new_sface();
      D.link_as_isolated_vertex(v1, f1);
      f1->mark() = mf1;
      
      if(mf1 == mf2 && mf1 == ml) {
	D.link_as_isolated_vertex(v2, f1);
      }
      else {
	SHalfloop_handle l = D.new_shalfloop_pair();
	SFace_handle f2 = D.new_sface();    
	D.link_as_isolated_vertex(v2, f2);
	D.link_as_loop(l,f1);
	D.link_as_loop(l->twin(),f2);
	l->circle() = Sphere_circle(faces_p->plane()); 
	l->twin()->circle() = l->circle().opposite();
	f2->mark() = mf2;
	l->mark() = l->twin()->mark() = ml;
	Halffacet_cycle_const_iterator fci = faces_p->facet_cycles_begin();
	SHalfedge_around_facet_const_circulator sea(fci);
	l->set_index(sea->twin()->get_index());
	l->twin()->set_index(sea->get_index());
      }
    }
    else {
      CGAL_NEF_TRACEN("edge is not isolated");
      SVertex_handle sv;
      SHalfedge_handle se1;
      SHalfedge_handle se2;
      SFace_handle sf;
      Sphere_circle c(f->plane());
      
      SHalfedge_handle next_edge;
      SHalfedge_around_svertex_const_circulator ec(E.out_edges(e)), ee(ec);
      CGAL_For_all(ec,ee) {
	Sphere_segment seg(ec->source()->point(), 
			   ec->source()->point().antipode(), 
			   ec->circle());
	Sphere_point sp(intersection(c, seg.sphere_circle()));
	CGAL_NEF_TRACEN(seg <<" has_on " << sp);
	if(!seg.has_on(sp))
	  sp = sp.antipode();
	sv = D.new_svertex(sp);
	CGAL_NEF_TRACEN("new svertex 3 " << normalized(sp));
	sv->mark() = BOP(ec->mark(), f->mark(), inv);
	Halffacet_const_handle f1, f2;
	if(inv) {
	  f1 = f;
	  f2 = ec->facet();
	} else {
	  f1 = ec->facet();
	  f2 = f;
	}
	if(f1->is_twin()) f1 = f1->twin();
	if(f2->is_twin()) f2 = f2->twin();
	A.hash_facet_pair(sv, f1, f2);

	se1 = D.new_shalfedge_pair(v1, sv);
	if(next_edge == SHalfedge_handle())
	  se2 = D.new_shalfedge_pair(sv, v2); 
	else
	  se2 = D.new_shalfedge_pair(sv, next_edge, -1);
	next_edge = se2->twin();
	se1->mark() = se1->twin()->mark() = BOP(ec->mark(), faces_p->incident_volume()->mark(), inv);
	se2->mark() = se2->twin()->mark() = BOP(ec->mark(), faces_p->twin()->incident_volume()->mark(), inv);
	mark_of_right_sface[se1] = ec->incident_sface()->mark();
	se1->circle() = se2->circle() = ec->circle();
	se1->twin()->circle() = se2->twin()->circle() = se1->circle().opposite();
	se1->set_index(ec->get_index());
	se2->set_index(ec->get_index());
	se1->twin()->set_index(ec->twin()->get_index());
	se2->twin()->set_index(ec->twin()->get_index());
      }
      
      SHalfedge_around_svertex_circulator ec2(D.out_edges(v1)), ee2(ec2);
      CGAL_For_all(ec2,ee2) {
	SHalfedge_around_svertex_circulator en(ec2);
	++en;
	se1 = D.new_shalfedge_pair(ec2->twin(), en->twin(), -1, 1);
	CGAL_NEF_TRACEN("new edge pair " << ec2->twin()->source()->vector() << 
			" -> " << en->twin()->source()->vector());
	se1->circle() = Sphere_circle(faces_p->plane());
	se1->twin()->circle() = se1->circle().opposite();
        se1->mark() = se1->twin()->mark() = BOP(mark_of_right_sface[ec2], faces_p->mark(), inv);

	Halffacet_cycle_const_iterator fci = faces_p->facet_cycles_begin();
	SHalfedge_around_facet_const_circulator sea(fci);
	se1->set_index(sea->twin()->get_index());
	se1->twin()->set_index(sea->get_index());
	
	sf = D.new_sface();
	sf->mark() = BOP(mark_of_right_sface[ec2], faces_p->incident_volume()->mark(), inv);
	D.link_as_face_cycle(se1,sf);
	sf = D.new_sface();
	sf->mark() = BOP(mark_of_right_sface[ec2], faces_p->twin()->incident_volume()->mark(), inv);
	D.link_as_face_cycle(se1->twin(),sf);
      }   
    }
    
    return D.sphere_map();
  }

  void assign_indices() {

    Halfedge_iterator e;
    CGAL_forall_edges(e, *this->sncp()) {
      e->set_index();
      e->twin()->set_index(e->get_index());
    }

    Halffacet_iterator f;
    CGAL_forall_halffacets(f, *this->sncp()) {
      Halffacet_cycle_iterator fci(f->facet_cycles_begin());
      SHalfedge_handle se(fci);
      se->set_index();
      int index(se->get_index());
      for(; fci != f->facet_cycles_end(); ++fci) {
	if(fci.is_shalfedge()) {
	  SHalfedge_around_facet_circulator c1(fci), c2(c1);
	  CGAL_For_all(c1,c2) {
	    c1->set_index(index);
	  }
	} else if(fci.is_shalfloop()) {
	  SHalfloop_handle sl(fci);
	  sl->set_index(index);
	}
      }
    }
  }

};

} //namespace CGAL
#endif //CGAL_SNC_CONSTRUCTOR_H
