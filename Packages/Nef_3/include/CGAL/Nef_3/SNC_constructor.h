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
//file:///usr/share/doc/mozilla-browser/localstart.html
// Author(s)     : Michael Seel       <seel@mpi-sb.mpg.de> 
//                 Miguel Granados    <granados@mpi-sb.mpg.de>
//                 Susan Hert         <hert@mpi-sb.mpg.de>
//                 Lutz Kettner       <kettner@mpi-sb.mpg.de>
//                 Peter Hachenberger <hachenberger@mpi-sb.mpg.de>
#ifndef CGAL_SNC_CONSTRUCTOR_H
#define CGAL_SNC_CONSTRUCTOR_H

#include <CGAL/basic.h>
#include <CGAL/functional.h> 
#include <CGAL/function_objects.h> 
#include <CGAL/Circulator_project.h>
#include <CGAL/Nef_S2/Normalizing.h>
#include <CGAL/Nef_3/bounded_side_3.h>
#include <CGAL/Nef_3/Pluecker_line_3.h>
#include <CGAL/Nef_3/SNC_decorator.h>
#include <CGAL/Nef_3/SNC_SM_overlayer.h>
#include <CGAL/Nef_S2/SM_point_locator.h>
#include <CGAL/Nef_3/SNC_FM_decorator.h>
#include <CGAL/Nef_3/SNC_point_locator.h>
#include <CGAL/Nef_3/SNC_sphere_map.h>
#include <CGAL/Nef_3/SNC_io_parser.h>
#include <CGAL/Nef_3/SNC_intersection.h>
#ifdef SM_VISUALIZOR
#include <CGAL/Nef_3/SNC_SM_visualizor.h>
#endif // SM_VISUALIZOR
#include <map>
#include <list>
#undef _DEBUG
#define _DEBUG 43
#include <CGAL/Nef_3/debug.h>
#include <CGAL/Nef_3/Progress_indicator.h>

CGAL_BEGIN_NAMESPACE

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
  
  int max;
  typedef typename T::Point_3                   Point_3;
  typedef typename T::RT                        RT;

  circle_lt(int m) :max(m) {};
  bool operator()(const Point_3& p1, const Point_3& p2) const { 
        
    const Quotient<RT> zero(RT(0));
    Quotient<RT> x[2];
    Quotient<RT> y[2];

    switch(max) {
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
      if(x[0] < zero) return (y[0]<y[1]);
      else return (y[0]>y[1]);
    }
    else {
      if(y[1] >= zero) return true;
      if(x[0]!=x[1]) return(x[0]>x[1]);
      if(x[0] > zero) return (y[0]<y[1]);
      else return  (y[0]>y[1]);
    }
    CGAL_assertion_msg(false, "control should not reach this line");
    return false;
  }
};

template <typename Point, typename Edge, class Decorator>
struct Halfedge_key {
  typedef Halfedge_key<Point,Edge,Decorator> Self;
  Point p; int i; Edge e;
  Decorator& D;
  Halfedge_key(Point pi, int ii, Edge ei, Decorator& Di ) : 
    p(pi), i(ii), e(ei), D(Di) {}
  Halfedge_key(const Self& k) : p(k.p), i(k.i), e(k.e), D(k.D) {}
  Self& operator=(const Self& k) { p=k.p; i=k.i; e=k.e; return *this; }
  bool operator==(const Self& k) const { return p==k.p && i==k.i; }
  bool operator!=(const Self& k) const { return !operator==(k); }
};

template <typename Point, typename Edge, class Decorator>
struct Halfedge_key_lt {
  typedef Halfedge_key<Point,Edge,Decorator> Key;
  typedef typename Point::R R;
  typedef typename R::Vector_3 Vector;
  typedef typename R::Direction_3 Direction;
  bool operator()( const Key& k1, const Key& k2) const { 
    if ( k1.p == k2.p ) 
      return (k1.i < k2.i);
    /* previous code: 
       else return CGAL::lexicographically_xyz_smaller(k1.p,k2.p); */
    Direction l(k1.e->vector());
    if( k1.i < 0) l = -l;
    return (Direction( k2.p - k1.p) == l); 
  }
};

template <typename Point, typename Edge, class Decorator>
std::ostream& operator<<(std::ostream& os, 
                         const Halfedge_key<Point,Edge,Decorator>& k )
{ os << k.p << " " << k.i; return os; }

template <typename R>
int sign_of(const CGAL::Plane_3<R>& h)
{ if ( h.c() != 0 ) return CGAL_NTS sign(h.c());
  if ( h.b() != 0 ) return CGAL_NTS sign(-h.b());
  return CGAL_NTS sign(h.a());
}

struct Plane_lt {
  template <typename R>
  bool operator()(const CGAL::Plane_3<R>& h1,
                  const CGAL::Plane_3<R>& h2) const
  { 
    typedef typename R::RT     RT;
    RT diff = h1.a()-h2.a();
    if ( (diff) != 0 ) return CGAL_NTS sign(diff) < 0;
    diff = h1.b()-h2.b();
    if ( (diff) != 0 ) return CGAL_NTS sign(diff) < 0;
    diff = h1.c()-h2.c();
    if ( (diff) != 0 ) return CGAL_NTS sign(diff) < 0;
    diff = h1.d()-h2.d(); return CGAL_NTS sign(diff) < 0;
  }
};

// ----------------------------------------------------------------------------
// SNC_constructor 
// ----------------------------------------------------------------------------

/*{\Manpage{SNC_constructor}{SNC}{overlay functionality}{O}}*/

template <typename SNC_structure_>
class SNC_constructor : public SNC_decorator<SNC_structure_>
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
  typedef CGAL::SNC_constructor<SNC_structure>               Self;
  typedef CGAL::SNC_decorator<SNC_structure>                 Base;
  typedef CGAL::SNC_decorator<SNC_structure>                 SNC_decorator;
  typedef typename CGAL::SNC_const_decorator<SNC_structure>  SNC_const_decorator;
  typedef CGAL::SNC_point_locator<SNC_decorator>       SNC_point_locator;
  typedef CGAL::SNC_FM_decorator<SNC_structure>              FM_decorator;
  typedef CGAL::SNC_intersection<SNC_structure>              SNC_intersection;

  typedef typename SNC_structure::Items                   Items;
  typedef typename SNC_structure::Sphere_map              Sphere_map;
  typedef CGAL::SM_decorator<Sphere_map>                 SM_decorator;  
  typedef CGAL::SNC_SM_overlayer<SM_decorator>           SM_overlayer;
  typedef CGAL::SM_const_decorator<Sphere_map>           SM_const_decorator;
  typedef CGAL::SM_point_locator<SM_const_decorator>     SM_point_locator;

  typedef typename SNC_structure::Vertex Vertex;
  typedef typename SNC_structure::Halfedge Halfedge;
  typedef typename SNC_structure::Halffacet Halffacet;
  typedef typename SNC_structure::Volume Volume;
  
  typedef typename SNC_structure::Vertex_iterator Vertex_iterator;
  typedef typename SNC_structure::Halfedge_iterator Halfedge_iterator;
  typedef typename SNC_structure::Halffacet_iterator Halffacet_iterator;
  typedef typename SNC_structure::Volume_iterator Volume_iterator;

  typedef typename SNC_structure::Vertex_handle Vertex_handle;
  typedef typename SNC_structure::Halfedge_handle Halfedge_handle;
  typedef typename SNC_structure::Halffacet_handle Halffacet_handle;
  typedef typename SNC_structure::Volume_handle Volume_handle;

  typedef typename SNC_structure::Vertex_const_handle Vertex_const_handle;
  typedef typename SNC_structure::Halfedge_const_handle Halfedge_const_handle;
  typedef typename SNC_structure::Halffacet_const_handle Halffacet_const_handle;
  typedef typename SNC_structure::Halffacet_const_iterator Halffacet_const_iterator;
  typedef typename SNC_structure::Volume_const_handle Volume_const_handle;

  typedef typename SNC_structure::SVertex_iterator SVertex_iterator;
  typedef typename SNC_structure::SHalfedge_iterator SHalfedge_iterator;
  typedef typename SNC_structure::SFace_iterator SFace_iterator;
  typedef typename SNC_structure::SHalfloop_iterator SHalfloop_iterator;

  typedef typename SNC_structure::SVertex SVertex;
  typedef typename SNC_structure::SHalfedge SHalfedge;
  typedef typename SNC_structure::SFace SFace;
  typedef typename SNC_structure::SHalfloop SHalfloop;

  typedef typename SNC_structure::SVertex_handle SVertex_handle;
  typedef typename SNC_structure::SHalfedge_handle SHalfedge_handle;
  typedef typename SNC_structure::SFace_handle SFace_handle;
  typedef typename SNC_structure::SHalfloop_handle SHalfloop_handle;

  typedef typename SNC_structure::SVertex_const_handle SVertex_const_handle; 
  typedef typename SNC_structure::SHalfedge_const_handle SHalfedge_const_handle; 
  typedef typename SNC_structure::SHalfloop_const_handle SHalfloop_const_handle; 
  typedef typename SNC_structure::SFace_const_handle SFace_const_handle; 

  typedef typename SNC_structure::Object_handle Object_handle;

  typedef typename SNC_structure::SHalfedge_around_facet_circulator SHalfedge_around_facet_circulator;
  typedef typename SNC_structure::SFace_cycle_iterator SFace_cycle_iterator;
  typedef typename SNC_structure::SFace_cycle_const_iterator SFace_cycle_const_iterator;
  typedef typename SNC_structure::Halffacet_cycle_iterator Halffacet_cycle_iterator;
  typedef typename SNC_structure::Halffacet_cycle_const_iterator Halffacet_cycle_const_iterator;
  typedef typename SNC_structure::Shell_entry_iterator Shell_entry_iterator;
  typedef typename SNC_structure::Shell_entry_const_iterator Shell_entry_const_iterator;

  typedef typename SNC_structure::Point_3 Point_3;
  typedef typename SNC_structure::Vector_3 Vector_3;
  typedef typename SNC_structure::Direction_3 Direction_3;
  typedef typename SNC_structure::Segment_3 Segment_3;
  typedef typename SNC_structure::Line_3 Line_3;
  typedef typename SNC_structure::Plane_3 Plane_3;
  typedef typename SNC_structure::Ray_3 Ray_3;
  typedef typename SNC_structure::Aff_transformation_3 Aff_transformation_3;

  typedef typename SNC_structure::Sphere_point Sphere_point;
  typedef typename SNC_structure::Sphere_segment Sphere_segment;
  typedef typename SNC_structure::Sphere_circle Sphere_circle;
  typedef typename SNC_structure::Sphere_direction Sphere_direction;

  typedef typename SNC_structure::Mark Mark;

  typedef typename SM_decorator::SHalfedge_around_svertex_circulator 
                                 SHalfedge_around_svertex_circulator;
  typedef typename SM_const_decorator::SHalfedge_around_svertex_const_circulator 
                                       SHalfedge_around_svertex_const_circulator;

  enum{NORMAL, CORNER, DEGENERATE};

  SNC_point_locator* pl;

  Vertex_handle minVertex;

  typedef void* GenPtr;

  typedef CGAL::Unique_hash_map<SFace_const_handle,unsigned int>  
                                                         Sface_shell_hash;
  typedef CGAL::Unique_hash_map<Halffacet_const_handle,unsigned int>  
                                                         Face_shell_hash;
  typedef CGAL::Unique_hash_map<SFace_const_handle,bool> SFace_visited_hash;
  typedef CGAL::Unique_hash_map<SFace_const_handle,bool> Shell_closed_hash;

  struct Shell_explorer {
    const SNC_decorator& D;
    Sface_shell_hash&  ShellSf;
    Face_shell_hash&   ShellF;
    //    Shell_closed_hash& Closed;
    SFace_visited_hash& Done;
    SFace_handle sf_min;
    int n;

    Shell_explorer(const SNC_decorator& Di, Sface_shell_hash& SSf, 
                   Face_shell_hash& SF, SFace_visited_hash& Vi) 
      : D(Di), ShellSf(SSf), ShellF(SF), Done(Vi), n(0) {}

    void visit(SFace_handle h) { 
      TRACEN("visit sf "<<D.point(D.vertex(h)));
      ShellSf[h]=n;
      Done[h]=true;
      if ( CGAL::lexicographically_xyz_smaller(
		 D.point(D.vertex(h)),D.point(D.vertex(sf_min)) )) 
	sf_min = h; 
    }

    void visit(Vertex_handle h) { 
      TRACEN("visit v  "<<D.point(h));
    }

    void visit(Halfedge_handle h) { 
      TRACEN("visit he "<<D.point(D.source(h)));
      //      SM_decorator SD(D.vertex(h));
      //      SFace_handle sf = D.source(h)->sfaces_begin();
      //      if( Closed[sf] ) {
      //	if(SD.is_isolated(h)){
      //	  if(!SD.has_shalfloop()) Closed[sf] = false;
      //	}
      //	else {
      //	  if(SD.first_out_edge(h) == SD.last_out_edge(h)) Closed[sf] = false;
      //	}
      //      }
    }

    void visit(Halffacet_handle h) { 
      TRACEN(h->plane()); 
      ShellF[h]=n; 
    }

    void visit(SHalfedge_handle se) {}
    void visit(SHalfloop_handle sl) {}

    SFace_handle& minimal_sface() { return sf_min; }

    void increment_shell_number() { 
      TRACEN("leaving shell "<<n);
      ++n; 
    }
  };

  SNC_constructor( SNC_structure& W, SNC_point_locator* spl = NULL) 
    : Base(W), pl(spl) {}
  /*{\Mcreate makes |\Mvar| a decorator of |W|.}*/

  void set_min_Vertex(Vertex_handle v_min) {minVertex=v_min;}
  Vertex_handle get_min_Vertex() {return minVertex;}

  Vertex_handle create_extended_box_corner(NT x, NT y, NT z,
                                  bool space=true, bool boundary=true) const { 

    CGAL_assertion(CGAL_NTS abs(x) == CGAL_NTS abs(y) &&
		   CGAL_NTS abs(y) == CGAL_NTS abs(z));
    
    CGAL_assertion(boundary == true);
    
    TRACEN("  constructing box corner on "<<Point_3(x,y,z)<<"...");
    Point_3 p = Infi_box::create_extended_point(x,y,z); 
    Vertex_handle v = this->sncp()->new_vertex(p , boundary);
    TRACEN( point(v));
    SM_decorator SD(&*v);
    Sphere_point sp[] = { Sphere_point(-x, 0, 0), 
			  Sphere_point(0, -y, 0), 
			  Sphere_point(0, 0, -z) };
    
  /* create box vertices */
    SVertex_handle sv[3];
    for(int vi=0; vi<3; ++vi) {
      sv[vi] = SD.new_svertex(sp[vi]);
      mark(sv[vi]) = boundary;
    }
    /* create facet's edge uses */
    //  Sphere_segment ss[3];
    SHalfedge_handle she[3];
    for(int si=0; si<3; ++si)
      she[si] = SD.new_shalfedge_pair(sv[si], sv[(si+1)%3]);
    
    for(int i=0; i<3;++i) {
    SD.circle(she[i]) = 
      Sphere_circle(Plane_3(sp[i],sp[(i+1)%3],Point_3(0,0,0)));
    SD.circle(SD.twin(she[i])) =  SD.circle(she[i]).opposite();
    SD.mark(she[i]) = SD.mark(SD.twin(she[i])) = boundary;
    }
    
    /* create facets */
    SFace_handle fi = SD.new_sface();
    SFace_handle fe = SD.new_sface();
    SD.link_as_face_cycle(she[0], fi);
    SD.link_as_face_cycle(SD.twin(she[0]), fe);
    
    Sphere_point p1 = SD.point(SD.source(she[0]));
    Sphere_point p2 = SD.point(SD.target(she[0]));
    Sphere_point p3 = SD.point(SD.target(SD.next(she[0])));
    if ( spherical_orientation(p1,p2,p3) > 0 ) {
      SD.mark(fi) = space;
      SD.mark(fe) = 0;
    }
    else {
      SD.mark(fi) = 0;
      SD.mark(fe) = space;
    }
    
    return v;
  }

  /*{\Mop produces the sphere map representing thp,e box corner in
          direction $(x,y,z)$.}*/

  template<typename Forward_iterator>  
  void add_outer_sedge_cycle(Vertex_handle v, 
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
  }

  template<typename Forward_iterator>  
  void add_inner_sedge_cycle(Vertex_handle v, 
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
    
    Sphere_point ps = E.point(e);
    ps = normalized(ps);
    SVertex_handle v1 = D.new_svertex(ps);
    SVertex_handle v2 = D.new_svertex(ps.antipode());
    TRACEN("new svertex 1 " << ps);
    TRACEN("new svertex 2 " << ps.antipode());
    D.mark(v1) = D.mark(v2) = E.mark(e);
    TRACEN("svertex 1 " << ps << " has mark " << D.mark(v1));
    TRACEN("svertex 2 " << ps.antipode() << " has mark " << D.mark(v2));
    
    if(E.is_isolated(e)) {
      TRACEN("edge is isolated");
      SFace_handle f = D.new_sface();
      D.mark(f) = E.mark(E.face(e));    
      D.link_as_isolated_vertex(v1, f);
      D.link_as_isolated_vertex(v2, f);
    }
    else {
      TRACEN("edge is not isolated");
      SHalfedge_handle se;
      SFace_handle sf;
      
      SHalfedge_around_svertex_const_circulator ec(E.out_edges(e)), ee(ec);
      CGAL_For_all(ec,ee) {
	Sphere_segment seg(E.point(E.source(ec)), 
			   E.point(E.source(ec)).antipode(), 
			   E.circle(ec));
	se = D.new_shalfedge_pair(v1, v2);
	D.mark(se) = D.mark(D.twin(se)) = E.mark(ec);
	mark_of_right_sface[se] = E.mark(E.face(ec));
	D.circle(se) = E.circle(ec);
	D.circle(D.twin(se)) = D.circle(se).opposite();
      }
      
      SHalfedge_around_svertex_circulator ec2(D.out_edges(v1)), ee2(ec2);
      CGAL_For_all(ec2,ee2) {
	sf = D.new_sface();
	D.mark(sf) = mark_of_right_sface[ec2];
	D.link_as_face_cycle(SHalfedge_handle(ec2),sf);
      }   
    }
    
    return v;
  }

  Vertex_handle create_from_plane(const Plane_3& pl, const Point_3& p, 
				  const Mark& bnd, 
				  const Mark& in, const Mark& out) const {
    typedef typename CGAL::SNC_const_decorator<SNC_structure> SNC_const_decorator;
    Vertex_handle v = this->sncp()->new_vertex( p, bnd);
    point(v) = p;
    Sphere_circle c(pl); // circle through origin parallel to h
    SM_decorator D(&*v);
    SHalfloop_handle l = D.new_shalfloop_pair();
    SFace_handle f1 = D.new_sface(), f2 = D.new_sface();
    D.link_as_loop(l,f1);
    D.link_as_loop(twin(l),f2);
    
    D.circle(l) = c; 
    D.circle(twin(l)) = c.opposite();
    D.mark(f1) = out;
    D.mark(f2) = in;
    D.mark(l) = D.mark(D.twin(l)) = bnd;
    return v;
  }

  Vertex_handle create_from_point_on_infibox_facet(const Point_3& p) const {
    
    CGAL_assertion(!Infi_box::is_standard(p));
    if(Infi_box::x_on_box(p))
      if(p.hx() > 0) 
	return create_from_plane(Plane_3(1,0,0,1), p, 1, 1, 0);
      else
	return create_from_plane(Plane_3(-1,0,0,1), p, 1, 1, 0);
    if(Infi_box::y_on_box(p))
      if(p.hy() > 0) 
	return create_from_plane(Plane_3(0,1,0,1), p, 1, 1, 0);
      else
	return create_from_plane(Plane_3(0,-1,0,1), p, 1, 1, 0);
    if(Infi_box::z_on_box(p))
      if(p.hz() > 0) 
	return create_from_plane(Plane_3(0,0,1,1), p, 1, 1, 0);
      else
	return create_from_plane(Plane_3(0,0,-1,1), p, 1, 1, 0);
    CGAL_assertion_msg(0, "something is wrong");
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
      CGAL_assertion_msg(false, "line of code shall not be reached");
    
    Vertex_handle v = this->sncp()->new_vertex( p, true);
    SM_decorator D(&*v);
    SVertex_handle v1 = D.new_svertex(sp);
    SVertex_handle v2 = D.new_svertex(sp.antipode());
    SHalfedge_handle e1 = D.new_shalfedge_pair(v1,v2);
    SHalfedge_handle e2 = D.new_shalfedge_pair(v1,v2);
    SFace_handle f1 = D.new_sface();
    SFace_handle f2 = D.new_sface();
    D.mark(v1) = D.mark(v2) = 
      D.mark(e1) = D.mark(e1->twin()) = 
      D.mark(e2) = D.mark(e2->twin()) = true;
    D.mark(f1) = false;
    D.mark(f2) = true;
    D.link_as_face_cycle(e1,f1);  
    D.link_as_face_cycle(e2,f2);
    D.circle(e1) = c1;
    D.circle(e1->twin()) = c1.opposite();
    D.circle(e1->snext()) = c2;
    D.circle(e1->snext()->twin()) = c2.opposite();
    
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
				  const Point_3& p) const {
    return create_from_plane(f->plane(), p, 
			     f->mark(), 
			     f->twin()->incident_volume()->mark(), 
			     f->incident_volume()->mark());
  }

  /*{\Mop produces the sphere map at point $p$ representing the local
     view of $f$. \precond $p$ is part of $f$.}*/

  Vertex_handle create_from_edge(Halfedge_const_handle e,
				 const Point_3& p) const {
    typedef typename CGAL::SNC_const_decorator<SNC_structure> SNC_const_decorator;
    
    CGAL_assertion(SNC_const_decorator::segment(e).has_on(p));
    Vertex_handle v = this->sncp()->new_vertex( p, e->mark());
    SM_decorator D(&*v);
    SM_const_decorator E(&*e->source());
    Sphere_point ps = e->point();
    SVertex_handle v1 = D.new_svertex(ps);
    SVertex_handle v2 = D.new_svertex(ps.antipode());
    v1->mark() = v2->mark() = e->mark();
    bool first = true;
    
    // SETDTHREAD(19*43*131);
    
    SHalfedge_const_handle ceee;
    TRACEN("---------------------" << e->center_vertex()->point());
    CGAL_forall_shalfedges(ceee,E)
      TRACEN("|" << ceee->circle() <<
	     "|" << ceee->mark() << 
	     " " << ceee->incident_sface()->mark());
    TRACEN(" ");
    
    if(E.is_isolated(e)) {
      SFace_handle f = D.new_sface();
      D.link_as_isolated_vertex(v1,f);
      D.link_as_isolated_vertex(v2,f);
      f->mark() = e->incident_sface()->mark();
    }
    
    SHalfedge_around_svertex_const_circulator ec1(E.out_edges(e)), ee(ec1);
    SHalfedge_handle e1,e2;
    CGAL_For_all(ec1,ee) {
      if (first) e1 = D.new_shalfedge_pair(v1,v2);
      else       e1 = D.new_shalfedge_pair(e1, e2, SM_decorator::AFTER, 
					   SM_decorator::BEFORE);
      e2 = D.twin(e1); 
      first = false;
    }
    
    SHalfedge_handle eee;
    CGAL_forall_shalfedges(eee,D)
      TRACEN("|" << eee->circle());
    TRACEN(" ");
       
    ec1 = E.out_edges(e);
    SHalfedge_around_svertex_circulator ec2(D.out_edges(v1));
    CGAL_For_all(ec1,ee) {
      TRACEN("|" << ec1->circle() <<
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
    
    TRACEN(" ");
    TRACEN("new vertex ");
    
    CGAL_forall_svertices(v1, D)
      TRACEN("|" << v1->point() << "|" << v1->mark());
    TRACEN(" ");

    CGAL_forall_shalfedges(eee,D)
      TRACEN("|" << eee->circle() <<
	     "|" << eee->mark() << 
	     " " << eee->incident_sface()->mark());
    TRACEN("---------------------");
    
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
    CGAL_forall_svertices(sv, E)
      VM[sv] = D.new_svertex(E.point(sv));
    
    SHalfedge_const_handle se;
    CGAL_forall_sedges(se, E) {
      EM[se] = D.new_shalfedge_pair();
      EM[E.twin(se)] = D.twin(EM[se]);
    }
    
    SFace_const_handle sf;
    CGAL_forall_sfaces(sf, E)
      FM[sf] = D.new_sface();
    
    SHalfloop_handle sl;
    if(E.has_shalfloop()) {
      sl = LM[E.shalfloop()] = D.new_shalfloop_pair();
      LM[E.twin(E.shalfloop())] = D.twin(sl);
      D.set_face(sl, FM[E.face(E.shalfloop())]);
      D.set_face(D.twin(sl), FM[E.face(E.twin(E.shalfloop()))]);
      D.circle(sl) = E.circle(E.shalfloop());
      D.circle(D.twin(sl)) = E.circle(E.twin(E.shalfloop()));
    }
    
    CGAL_forall_svertices(sv, E) {
      D.set_first_out_edge(VM[sv], EM[E.first_out_edge(sv)]);
      D.set_face(VM[sv], FM[E.face(sv)]);
      D.mark(VM[sv]) = E.mark(sv);
    }
    
    CGAL_forall_shalfedges(se, E) {
      D.mark(EM[se]) = D.mark(D.twin(EM[se])) = E.mark(se);
      D.set_source(EM[se], VM[E.source(se)]);
      D.set_prev(EM[se], EM[E.previous(se)]); 
      D.set_next(EM[se], EM[E.next(se)]);
      D.set_face(EM[se], FM[E.face(se)]);
      D.circle(EM[se]) = E.circle(se);
    }
    
    CGAL_forall_sfaces(sf, E) {
      D.mark(FM[sf]) = E.mark(sf);
      SFace_cycle_const_iterator sfc;
      for(sfc = sf->sface_cycles_begin(); sfc != sf->sface_cycles_end(); ++sfc) {
	if(sfc.is_svertex())
	  D.store_sm_boundary_object(VM[SVertex_const_handle(sfc)], FM[sf]);
	else if(sfc.is_shalfedge())
	  D.store_sm_boundary_object(EM[SHalfedge_const_handle(sfc)], FM[sf]);
	else if(sfc.is_shalfloop())
	  D.store_sm_boundary_object(LM[SHalfloop_const_handle(sfc)], FM[sf]);
	else CGAL_assertion_msg(0,"damn wrong handle.");
      }
    }
    
    return vout;
  }

  template<typename Selection>
  Sphere_map* create_edge_facet_overlay( Halfedge_const_handle e, 
		                	   Halffacet_const_handle f,
					   const Point_3& p,
					 const Selection& BOP, bool inv) {

#ifdef CGAL_NEF3_TIMER_SPHERE_SWEEPS
    ++number_of_edge_facet_overlays;
#endif

    TRACEN("edge facet overlay " << p);
    
    Unique_hash_map<SHalfedge_handle, Mark> mark_of_right_sface;
    
    SM_decorator D(&*this->sncp()->new_vertex(p, BOP(e->mark(), f->mark())));
    SM_const_decorator E(&*e->source());
    
    Sphere_point ps = e->point();
    ps = normalized(ps);
    SVertex_handle v1 = D.new_svertex(ps);
    SVertex_handle v2 = D.new_svertex(ps.antipode());
    TRACEN("new svertex 1 " << ps);
    TRACEN("new svertex 2 " << ps.antipode());
    Halffacet_const_handle faces_p(f);
    Vector_3 vec(ps-CGAL::ORIGIN);
    if(faces_p->plane().oriented_side(p+vec) == ON_NEGATIVE_SIDE)
      faces_p = faces_p->twin();
    v1->mark() = BOP(e->mark(), faces_p->incident_volume()->mark(), inv);
    v2->mark() = BOP(e->mark(), faces_p->twin()->incident_volume()->mark(), inv);
    TRACEN("svertex 1 " << ps << " has mark " << v1->mark());
    TRACEN("svertex 2 " << ps.antipode() << " has mark " << v2->mark());
    
    if(E.is_isolated(e)) {
      TRACEN("edge is isolated");
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
      TRACEN("edge is not isolated");
      SVertex_handle sv;
      SHalfedge_handle se1;
      SHalfedge_handle se2;
      SFace_handle sf;
      Sphere_circle c(f->plane());
      
      SHalfedge_handle next_edge;
      SHalfedge_around_svertex_const_circulator ec(E.out_edges(e)), ee(ec);
      CGAL_For_all(ec,ee) {
	Sphere_segment seg(E.point(E.source(ec)), 
			   E.point(E.source(ec)).antipode(), 
			   E.circle(ec));
	Sphere_point sp(intersection(c, seg.sphere_circle()));
	TRACEN(seg <<" has_on " << sp);
	if(!seg.has_on(sp))
	  sp = sp.antipode();
	sv = D.new_svertex(sp);
	TRACEN("new svertex 3 " << normalized(sp));
	sv->mark() = BOP(ec->mark(), f->mark(), inv);
	se1 = D.new_shalfedge_pair(v1, sv);
	if(next_edge == SHalfedge_handle())
	  se2 = D.new_shalfedge_pair(sv, v2); 
	else
	  se2 = D.new_shalfedge_pair(sv, next_edge, -1);
	next_edge = twin(se2);
	se1->mark() = se1->twin()->mark() = BOP(E.mark(ec), faces_p->incident_volume()->mark(), inv);
	se2->mark() = se2->twin()->mark() = BOP(E.mark(ec), faces_p->twin()->incident_volume()->mark(), inv);
	mark_of_right_sface[se1] = E.mark(E.face(ec));
	D.circle(se1) = D.circle(se2) = E.circle(ec);
	D.circle(D.twin(se1)) = D.circle(D.twin(se2)) = D.circle(se1).opposite();
      }
      
      SHalfedge_around_svertex_circulator ec2(D.out_edges(v1)), ee2(ec2);
      CGAL_For_all(ec2,ee2) {
	SHalfedge_around_svertex_circulator en(ec2);
	++en;
	se1 = D.new_shalfedge_pair(twin(ec2), twin(en), -1, 1);
	TRACEN("new edge pair " << ssource(twin(ec2))->vector() << 
	       " -> " << ssource(twin(en))->vector());
	D.circle(se1) = Sphere_circle(faces_p->plane());
	D.circle(D.twin(se1)) = D.circle(se1).opposite();
	D.mark(se1) = D.mark(D.twin(se1)) = BOP(mark_of_right_sface[ec2], faces_p->mark(), inv);
	
	sf = D.new_sface();
	D.mark(sf) = BOP(mark_of_right_sface[ec2], faces_p->incident_volume()->mark(), inv);
	D.link_as_face_cycle(se1,sf);
	sf = D.new_sface();
	D.mark(sf) = BOP(mark_of_right_sface[ec2], faces_p->twin()->incident_volume()->mark(), inv);
	D.link_as_face_cycle(D.twin(se1),sf);
      }   
    }
    
    return D.sphere_map();
  }

 private:
  void pair_up_halfedges() const {
  /*{\Mop pairs all halfedge stubs to create the edges in 3-space.}*/
    
  //  SETDTHREAD(43*61);
    TRACEN(">>>>>pair_up_halfedges");
    typedef Halfedge_key< Point_3, Halfedge_handle, SNC_decorator>
      Halfedge_key;
    typedef Halfedge_key_lt< Point_3, Halfedge_handle, SNC_decorator> 
      Halfedge_key_lt;
    typedef std::list<Halfedge_key>  Halfedge_list;
    
    typedef typename Standard_kernel::Kernel_tag   Kernel_tag;
    typedef CGAL::Pluecker_line_3<Kernel_tag,Standard_kernel> Pluecker_line_3;
    typedef CGAL::Pluecker_line_lt        Pluecker_line_lt;
    typedef std::map< Pluecker_line_3, Halfedge_list, Pluecker_line_lt> 
      Pluecker_line_map;
    
    SNC_decorator D(*this);
    Pluecker_line_map M;
    Pluecker_line_map M2;
    Pluecker_line_map M3;
    Pluecker_line_map M4;
    
    //  Progress_indicator_clog progress( 2*this->sncp()->number_of_halfedges(), 
    //				    "SNC_constructor: pairing up edges...");
    
    NT eval(Infi_box::compute_evaluation_constant_for_halfedge_pairup(*this->sncp()));;
    
    Halfedge_iterator e;
    CGAL_forall_halfedges(e,*this->sncp()) {
      //    progress++;
      Point_3 p = point(vertex(e));
      Point_3 q = p + e->vector();
      TRACE(" segment("<<p<<", "<<q<<")"<<
	    " direction("<<e->vector()<<")");
      Standard_point_3 sp = Infi_box::standard_point(p,eval);
      Standard_point_3 sq = Infi_box::standard_point(q,eval);
      Pluecker_line_3  l( sp, sq);
      
      int inverted;
      l = categorize( l, inverted);
      
      if(Infi_box::is_edge_on_infibox(e))
	if(Infi_box::is_type4(e))
	  M4[l].push_back(Halfedge_key(p,inverted,e, D));
	else
	  if(Infi_box::is_type3(e))
	    M3[l].push_back(Halfedge_key(p,inverted,e, D));
	  else
	    M2[l].push_back(Halfedge_key(p,inverted,e, D));
      else
	M[l].push_back(Halfedge_key(p,inverted,e,D));
      
      // the following trace crashes when compiling with optimizations (-O{n})
      //TRACEN(Infi_box::standard_point(point(vertex(e)))+
      //   Vector_3(point(e))); 
      
      TRACEN(" line("<<l<<")"<<" inverted="<<inverted);
    }
    
    typename Pluecker_line_map::iterator it;
    
    CGAL_forall_iterators(it,M4) {
      //    progress++;
      it->second.sort(Halfedge_key_lt());
      TRACEN("search opposite  "<<it->first); 
      typename Halfedge_list::iterator itl;
      CGAL_forall_iterators(itl,it->second) {
	Halfedge_handle e1 = itl->e;
	++itl; 
	CGAL_assertion(itl != it->second.end());
	Halfedge_handle e2 = itl->e;
	TRACEN("    " << point(vertex(e1)) << " -> " << point(vertex(e2)));
	TRACEN(e1->vector()<<" -> "<<-e2->vector());
	CGAL_assertion(normalized(e1->vector())==normalized(-e2->vector()));
	make_twins(e1,e2);
	CGAL_assertion(mark(e1)==mark(e2));
	
	// discard temporary sphere_point ?
      }
    }
    
    CGAL_forall_iterators(it,M3) {
      //    progress++;
      it->second.sort(Halfedge_key_lt());
      TRACEN("search opposite  "<<it->first); 
      typename Halfedge_list::iterator itl;
      CGAL_forall_iterators(itl,it->second) {
	Halfedge_handle e1 = itl->e;
	++itl; 
	CGAL_assertion(itl != it->second.end());
	Halfedge_handle e2 = itl->e;
	TRACEN("    " << point(vertex(e1)) << " -> " << point(vertex(e2)));
	TRACEN(e1->vector()<<" -> "<<-e2->vector());
	CGAL_assertion(normalized(e1->vector())==normalized(-e2->vector()));
	make_twins(e1,e2);
	CGAL_assertion(mark(e1)==mark(e2));
	
	// discard temporary sphere_point ?
      }
    }
    
    CGAL_forall_iterators(it,M2) {
      //    progress++;
      it->second.sort(Halfedge_key_lt());
      TRACEN("search opposite  "<<it->first); 
      typename Halfedge_list::iterator itl;
      CGAL_forall_iterators(itl,it->second) {
	Halfedge_handle e1 = itl->e;
	++itl; 
	CGAL_assertion(itl != it->second.end());
	Halfedge_handle e2 = itl->e;
	TRACEN("    " << point(vertex(e1)) << " -> " << point(vertex(e2)));
	TRACEN(e1->vector()<<" -> "<<-e2->vector());
	CGAL_assertion(normalized(e1->vector())==normalized(-e2->vector()));
	make_twins(e1,e2);
	CGAL_assertion(mark(e1)==mark(e2));
	
	// discard temporary sphere_point ?
      }
    }
    
    CGAL_forall_iterators(it,M) {
      //    progress++;
      it->second.sort(Halfedge_key_lt());
      TRACEN("search opposite  "<<it->first); 
      typename Halfedge_list::iterator itl;
      CGAL_forall_iterators(itl,it->second) {
	Halfedge_handle e1 = itl->e;
	++itl; 
	CGAL_assertion(itl != it->second.end());
	Halfedge_handle e2 = itl->e;
	TRACEN("    " << point(vertex(e1)) << " -> " << point(vertex(e2)));
	TRACEN(e1->vector()<<" -> "<< -e2->vector());
	//	CGAL_assertion(normalized(e1->vector())==normalized(-e2->vector()));
	CGAL_assertion(point(vertex(e1)) != point(vertex(e2)));
	CGAL_assertion(mark(e1)==mark(e2));
	make_twins(e1,e2);
	
	// discard temporary sphere_point ?
      }
    }
  }

  void link_shalfedges_to_facet_cycles() const {
  /*{\Mop creates all non-trivial facet cycles from sedges. 
  \precond |pair_up_halfedges()| was called before.}*/
    
  //  SETDTHREAD(43*31);
    TRACEN(">>>>>link_shalfedges_to_facet_cycles");
    
    //  Progress_indicator_clog progress
    //    (this->sncp()->number_of_halfedges(), 
    //     "SNC_constructor: building facet cycles...");
    
    Halfedge_iterator e;
    CGAL_forall_edges(e,*this->sncp()) {
      //    progress++;
      TRACEN("");
      TRACEN(PH(e));
      Halfedge_iterator et = twin(e);
      SM_decorator D(&*vertex(e)), Dt(&*vertex(et));
      TRACEN(point(vertex(e)));
      if ( D.is_isolated(e) ) continue;
      SHalfedge_around_svertex_circulator ce(D.first_out_edge(e)),cee(ce);
      SHalfedge_around_svertex_circulator cet(Dt.first_out_edge(et)),cete(cet);
      
      CGAL_For_all(cet,cete)
	if ( Dt.circle(cet) == D.circle(ce).opposite() && 
	     twin(Dt.source(cet)) == D.source(ce) ) 
          break;
      
      /*    DEBUG 
	    
      if( Dt.circle(cet) != D.circle(ce).opposite() )
	TRACEN("assertion failed!");
      
      TRACEN("vertices " << point(vertex(e)) << 
      "    "      << point(vertex(et)));
      
      
      SHalfedge_around_svertex_circulator sc(D.first_out_edge(e));
      SHalfedge_around_svertex_circulator sct(Dt.first_out_edge(et));
      CGAL_For_all(sc,cee)
      TRACEN("sseg@E addr="<<&*sc<<
      " src="<<D.point(D.source(sc))<<
      " tgt="<<D.point(D.target(sc))<<std::endl<<
      " circle=" << D.circle(sc));
      CGAL_For_all(sct,cete)
      TRACEN("sseg@ET addr="<<&*sct<<
      " src="<<Dt.point(Dt.source(sct))<<
      " tgt="<<Dt.point(Dt.target(sct))<<std::endl<<
      " circle=" << Dt.circle(sct));
      
      #ifdef SM_VISUALIZOR
      typedef CGAL::SNC_SM_visualizor<SNC_structure> SMV;
      CGAL::OGL::add_sphere();
      SMV V(vertex(e), CGAL::OGL::spheres_.back());
      V.draw_map();
      SMV Vt(vertex(et), CGAL::OGL::spheres_.back());
      Vt.draw_map();
      CGAL::OGL::start_viewer();
      char c;
      cin >> c;
      #endif
      */

      CGAL_assertion( Dt.circle(cet) == D.circle(ce).opposite() ); 
      CGAL_assertion( twin(Dt.source(cet)) == D.source(ce)); 
      CGAL_For_all(ce,cee) { 
	TRACEN("circles " << Dt.circle(cet) << "   " << D.circle(ce) << 
	       " sources " << Dt.point(Dt.target(cet)) << 
	       "   " << D.point(D.target(ce)));
	CGAL_assertion( Dt.circle(cet) == D.circle(ce).opposite() ); 
	CGAL_assertion( twin(Dt.source(cet)) == D.source(ce)); 
	CGAL_assertion(ce->mark()==cet->mark());
	link_as_prev_next_pair(Dt.twin(cet),ce);
	link_as_prev_next_pair(D.twin(ce),cet);
	--cet; // ce moves ccw, cet moves cw
      }
    }
  }

  void categorize_facet_cycles_and_create_facets() const {
  /*{\Mop collects all facet cycles incident to a facet and creates
  the facets. \precond |link_shalfedges_to_facet_cycles()| was called
  before.}*/

    //    SETDTHREAD(43*31*23);
    TRACEN(">>>>>categorize_facet_cycles_and_create_facets");
    
    typedef std::list<Object_handle> Object_list;
    typedef std::map<Plane_3, Object_list, Plane_lt> 
      Map_planes;
    
    //  Progress_indicator_clog progress
    //    (this->sncp()->number_of_shalfedges()+this->sncp()->number_of_shalfloops(), 
    //     "SNC_constructor: categorizing facet cycles...");
    
    Map_planes M;
    SHalfedge_iterator e;
    CGAL_forall_shalfedges(e,*this->sncp()) {
      //    progress++;
      Sphere_circle c(circle(e));
      Plane_3 h = c.plane_through(point(vertex(e)));
      SM_decorator SD(&*vertex(e));
      TRACEN("\n" << point(target(SD.source(e))) <<" - "<< point(vertex(e)) <<" - "<< 
	     point(target(SD.target(e))) << 
	     " has plane " << h << " has circle " << circle(e) << 
	     " has signum " << sign_of(h));
      if ( sign_of(h)<0 ) continue;
      M[normalized(h)].push_back(Object_handle(twin(e)));
      TRACEN(" normalized as " << normalized(h)); 
    }
    SHalfloop_iterator l;
    CGAL_forall_shalfloops(l,*this->sncp()) {
      Sphere_circle c(circle(l));
      Plane_3 h = c.plane_through(point(vertex(l))); 
      if ( sign_of(h)<0 ) continue;
      // CGAL_assertion( h == normalized(h));
      M[normalized(h)].push_back(Object_handle(twin(l)));
    }
    
    //  Progress_indicator_clog progress2 
    //    ( M.size(), "SNC_constructor: creating facets...");
    
#ifdef CGAL_NEF3_TIMER_PLANE_SWEEPS
  number_of_plane_sweeps=0;
  timer_plane_sweeps.reset();
#endif

    typename Map_planes::iterator it;
    CGAL_forall_iterators(it,M) { 
      //    progress2++;
      TRACEN("  plane "<<it->first<<"             "<<(it->first).point());
      FM_decorator D(*this->sncp());
      D.create_facet_objects(it->first,it->second.begin(),it->second.end());
    }
    //    SETDTHREAD(1);
  }

  void create_volumes() {
  /*{\Mop collects all shells incident to a volume and creates the
  volumes.  \precond |categorize_facet_cycles_and_creating_facets()| was
  called before.}*/

#ifdef CGAL_NEF3_TIMER_POINT_LOCATION
    timer_ray_shooting.reset();
#endif 

    //    SETDTHREAD(37*43*19*47);
    TRACEN(">>>>>create_volumes");
    Sface_shell_hash     ShellSf(-1);
    Face_shell_hash      ShellF(-1);
    SFace_visited_hash Done(false);
    Shell_explorer V(*this,ShellSf,ShellF,Done);
    std::vector<SFace_handle> MinimalSFace;
    std::vector<SFace_handle> EntrySFace;
    std::vector<bool> Closed;
    
  //  Progress_indicator_clog progress
  //    (this->sncp()->number_of_sfaces(), "SNC_construction: building shells...");
    
    SFace_iterator f;
    // First, we classify all the Shere Faces per Shell.  For each Shell we
    //     determine its minimum lexicographyly vertex and we check wheter the
    //     Shell encloses a region (closed surface) or not. 
    CGAL_forall_sfaces(f,*this->sncp()) {
      //    progress++;
      TRACEN("sface in " << ShellSf[f]);
      if ( Done[f] ) 
	continue;
      V.minimal_sface() = f;
      visit_shell_objects(f,V);
      
      MinimalSFace.push_back(V.minimal_sface());
      EntrySFace.push_back(f);
      V.increment_shell_number();
      TRACEN("sface out " << ShellSf[f]);
    }
    
    for(unsigned int i=0; i<EntrySFace.size(); ++i)
      Closed.push_back(false);
    
    Halffacet_iterator hf;
    CGAL_forall_facets(hf,*this) 
      if(ShellF[hf] != ShellF[twin(hf)]) {
	Closed[ShellF[hf]] = true;
	Closed[ShellF[twin(hf)]] = true;
      }
    
    CGAL_assertion( pl != NULL);

#ifdef CGAL_NEF3_TIMER_INITIALIZE_KDTREE
    CGAL::Timer timer_initialize_kdtree;
    timer_initialize_kdtree.start();
#endif
    pl->initialize(this->sncp()); // construct the point locator 
#ifdef CGAL_NEF3_TIMER_INITIALIZE_KDTREE
    timer_initialize_kdtree.stop();
    if(cgal_nef3_timer_on)
      std::cout << "Runtime_initialize_kdtree: " 
		<< timer_initialize_kdtree.time() << std::endl;
#endif

    //  Progress_indicator_clog progress2
    //    ( MinimalSFace.size(), "SNC_constructor: creating outer volumes...");
    
    //   then, we determine the Shells which correspond to Volumes via a ray
    //   shootting in the direction (-1,0,0) over the Sphere_map of the minimal 
    //   vertex.  The Shell corresponds to a Volume if the object hit belongs 
    //   to another Shell. 
    
    this->sncp()->new_volume(); // outermost volume (nirvana)
    if(MinimalSFace.size() == 0) return;
    Vertex_handle v_min = vertex(MinimalSFace[0]);
    for( unsigned int i = 0; i < MinimalSFace.size(); ++i) {
      //    progress2++;
      Vertex_handle v = vertex(MinimalSFace[i]);
      if(CGAL::lexicographically_xyz_smaller(point(v),point(v_min)))
	v_min=v;
      TRACEN( "Shell #" << i << " minimal vertex: " << point(v));
      SM_point_locator D((Sphere_map*) &*v);
      Object_handle o = D.locate(Sphere_point(-1,0,0));
      SFace_const_handle sfc;
      if( !CGAL::assign(sfc, o) || ShellSf[sfc] != i) {
	TRACEN("the shell encloses a volume");
	TRACEN("sface hit? "<<CGAL::assign(sfc,o));
	if( CGAL::assign(sfc, o)) TRACEN("sface on another shell? "<<ShellSf[sfc]);
	SFace_handle f = MinimalSFace[i];
	CGAL_assertion( ShellSf[f] == i );
	if( Closed[i] ) {
	  TRACEN("Shell #" << i << " is closed");
	  SM_decorator SD(&*v);
	  Volume_handle c = this->sncp()->new_volume();
	  mark(c) = SD.mark(f); // TODO: test if line is redundant
	  link_as_inner_shell(f, c );
	  TRACE( "Shell #" << i <<" linked as inner shell");
	  TRACEN( "(sface" << (CGAL::assign(sfc,o)?"":" not") << " hit case)");
	}
      }
    }
    minVertex = v_min;
    
    // finaly, we go through all the Shells which do not correspond to a Volume 
    //     and we assign them to its enclosing Volume determined via a facet below
    //     check. 
    
    //  Progress_indicator_clog progress3
    //    ( this->sncp()->number_of_sfaces(), 
    //      "SNC_constructor: finding nesting structure...");
    
    CGAL_forall_sfaces(f,*this->sncp()) {
      //    progress3++;
      if ( volume(f) != Volume_handle() ) 
	continue;
      TRACEN( "Outer shell #" << ShellSf[f] << " volume?");
      Volume_handle c = determine_volume( MinimalSFace[ShellSf[f]], 
					  MinimalSFace, ShellSf );
      c->mark() = f->mark();
      link_as_outer_shell( f, c );
    } 
  }

  Halffacet_handle get_facet_below( Vertex_handle vi, 
				    const std::vector< SFace_handle>& MinimalSFace, 
				    const Sface_shell_hash&  Shell) const {
    // {\Mop determines the facet below a vertex |vi| via ray shooting. }
    
    Halffacet_handle f_below;
    Point_3 p = point(vi);
    if(!Infi_box::is_standard(p))
      return Halffacet_handle();
    
    Ray_3 ray = Ray_3(p, Direction_3(-1,0,0));
#ifdef CGAL_NEF3_TIMER_POINT_LOCATION
    number_of_ray_shooting_queries++;
    timer_ray_shooting.start();
#endif 
    Object_handle o = pl->shoot(ray);
#ifdef CGAL_NEF3_TIMER_POINT_LOCATION
    timer_ray_shooting.stop();
#endif 
    // The ray here has an special property since it is shooted from the lowest
    // vertex in a shell, so it would be expected that the ray goes along the
    // interior of a volume before it hits a 2-skeleton element.
    // Unfortunatelly, it seems to be possible that several shells are incident
    // to this lowest vertex, and in consequence, the ray could also go along
    // an edge or a facet belonging to a different shell.
    // This fact invalidates the precondition of the get_visible_facet method,
    // (the ray must pierce the local view of the hit object in a sface).
    // This situation has not been analyzed and has to be verified. Granados.
    Vertex_handle v;
    Halfedge_handle e;
    Halffacet_handle f;
    TRACEN("get_facet_below");
    if( CGAL::assign(v, o)) {
      TRACEN("facet below from from vertex...");
      f_below = get_visible_facet(v, ray);
      if( f_below == Halffacet_handle()) {
	CGAL_assertion(v->sfaces_begin() == v->sfaces_last());
	f_below = get_facet_below(vertex(MinimalSFace[Shell[v->sfaces_begin()]]), 
				  MinimalSFace,Shell);
      }
    }
    else if( CGAL::assign(e, o)) {
      TRACEN("facet below from from edge...");
      f_below = get_visible_facet(e, ray);
      if( f_below == Halffacet_handle()) {
	CGAL_assertion(vertex(e)->sfaces_begin() == vertex(e)->sfaces_last());
	f_below = get_facet_below(vertex(MinimalSFace[Shell[vertex(e)->sfaces_begin()]]), 
				  MinimalSFace, Shell);
      }
    }
    else if( CGAL::assign(f, o)) {
      TRACEN("facet below from from facet...");
      f_below = get_visible_facet(f, ray);
      CGAL_assertion( f_below != Halffacet_handle());
    }
    else { TRACEN("no facet below found..."); }
    return f_below;
  }  

  Volume_handle determine_volume( SFace_handle sf, 
                const std::vector< SFace_handle>& MinimalSFace, 
				  const Sface_shell_hash&  Shell ) const {
    //{\Mop determines the volume |C| that a shell |S| pointed by |sf| 
    //  belongs to.  \precondition |S| separates the volume |C| from an enclosed
    //  volume.}
    
    TRACEN("determine volume");
    Vertex_handle v_min = vertex(MinimalSFace[Shell[sf]]);
    
    Halffacet_handle f_below = get_facet_below(v_min, MinimalSFace, Shell);
    if ( f_below == Halffacet_handle())
      return Base(*this).volumes_begin();
    Volume_handle c = volume(f_below);
    if( c != Volume_handle()) {
      TRACE( "Volume " << &*c << " hit ");
      TRACEN("(Shell #" << Shell[adjacent_sface(f_below)] << ")");
      return c;
    }
    SFace_handle sf_below = adjacent_sface(f_below);
    TRACE( "Shell not assigned to a volume hit ");
    TRACEN( "(Inner shell #" << Shell[sf_below] << ")");
    c = determine_volume( sf_below, MinimalSFace, Shell);
    link_as_inner_shell( sf_below, c);
    return c;
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
    if(!flag)
      if(*s1 == *s2) return normalized(*s1);
      else if(*s1 == *t2) return normalized(*s1);
      else if(*t1 == *s2) return normalized(*t1);
      else if(*t1 == *t2) return normalized(*t1);
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
      if(point(*vi) == point(*vinext)) {
	CGAL_assertion(!Infi_box::is_complex_facet_infibox_intersection(**vi));
	this->sncp()->delete_vertex(*vi);
      }
    }
    
    NT eval(Infi_box::compute_evaluation_constant_for_halfedge_pairup(*this->sncp()));
    
    bool res;
    do {
    
      res = false;
      typedef Halfedge_key< Point_3, Halfedge_handle, SNC_decorator>
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
	Point_3 p = point(vertex(e));
	Point_3 q = p + e->vector();
	Standard_point_3 sp = Infi_box::standard_point(p,eval);
	Standard_point_3 sq = Infi_box::standard_point(q,eval);
	Pluecker_line_3 l( sp, sq);
	
	int inverted;
	l = categorize( l, inverted);
	
	TRACEN(" segment("<<p<<", "<<q<<")"<<
	       " direction("<<e->vector()<<")"<<
	       " line("<<l<<")"<<" inverted="<<inverted);
	
	if(Infi_box::is_edge_on_infibox(e))
	  if(Infi_box::is_type4(e))
	    M4[l].push_back(Halfedge_key(p,inverted,e, D));
	  else
	    if(Infi_box::is_type3(e))
	      M3[l].push_back(Halfedge_key(p,inverted,e, D));
	    else
	      M2[l].push_back(Halfedge_key(p,inverted,e, D));
      }
      
      typename Pluecker_line_map::iterator it;
      
      CGAL_forall_iterators(it,M4) {
	TRACEN("search opposite  "<<it->first);
	it->second.sort(Halfedge_key_lt());
	typename Halfedge_list::iterator itl;
	CGAL_forall_iterators(itl,it->second) {
	  Halfedge_handle e1 = itl->e;
	  TRACE("    " << point(vertex(e1)) << " -> ");
	  ++itl; 
	  if(itl == it->second.end()) {
	    erase_vertex[e1->source()] = true;
	    break;
	  }
	  Halfedge_handle e2 = itl->e;
	  TRACE(point(vertex(e2)));
	  if(normalized(e1->vector())!=normalized(-e2->vector())) {
	    erase_vertex[e1->source()] = true;
	    --itl;
	    TRACE("   failed ");
	  }
	  TRACEN("");
	}
	TRACEN("");
	TRACEN("");
      }
      
      CGAL_forall_iterators(it,M3) {
	TRACEN("search opposite  "<<it->first); 
	it->second.sort(Halfedge_key_lt());
	typename Halfedge_list::iterator itl;
	CGAL_forall_iterators(itl,it->second) {
	  Halfedge_handle e1 = itl->e;
	  TRACE("    " << point(vertex(e1)) << " -> ");
	  ++itl; 
	  if(itl == it->second.end()) {
	    erase_vertex[e1->source()] = true;
	    break;
	  }
	  Halfedge_handle e2 = itl->e;
	  TRACE(point(vertex(e2)));
	  if(normalized(e1->vector())!=normalized(-e2->vector())) {
	    erase_vertex[e1->source()] = true;
	    --itl;
	    TRACE("   failed ");
	  }
	  TRACEN("");
	}
	TRACEN("");
	TRACEN("");
      }
      
      CGAL_forall_iterators(it,M2) {
	TRACEN("search opposite  "<<it->first); 
	it->second.sort(Halfedge_key_lt());
	typename Halfedge_list::iterator itl;
	CGAL_forall_iterators(itl,it->second) {
	  Halfedge_handle e1 = itl->e;
	  TRACE("    " << point(vertex(e1)) << " -> ");
	  ++itl; 
	  if(itl == it->second.end()) {
	    erase_vertex[e1->source()] = true;
	    break;
	  }
	  Halfedge_handle e2 = itl->e;
	  TRACE(point(vertex(e2)));
	  if(normalized(e1->vector())!=normalized(-e2->vector())) {
	    erase_vertex[e1->source()] = true;
	    --itl;
	    TRACE("   failed ");
	  }
	  TRACEN("");
	}
	TRACEN("");
	TRACEN("");
      }
      
      Vertex_iterator v;
      CGAL_forall_vertices(v, *this->sncp())
	if(erase_vertex[v]) {
	  TRACEN("erase " << v->point());
	  if(Infi_box::is_infibox_corner(v->point()))
	    recreate.push_back(v->point());
	  this->sncp()->delete_vertex(v);
	  res = true;
	}
      
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
      TRACEN(*p);
    
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
	//      std::cerr << src << "->" << trg << " has on " << *pi << src.x() << " : " << std::endl;
	//      std::cerr << (src.x()-pi->x() <= 0) << "|" << (src.y()-pi->y() <= 0) << "|" << (src.z()-pi->z() <= 0) << std::endl;
	//      std::cerr << (pi->x()-trg.x() <= 0) << "|" << (pi->y()-trg.y() <= 0) << "|" << (pi->z()-trg.z() <= 0) << std::endl;
	if((src.x()-pi->x() <= 0 && pi->x()-trg.x() <= 0 || 
	    src.x()-pi->x() >= 0 && pi->x()-trg.x() >= 0) &&
	   (src.y()-pi->y() <= 0 && pi->y()-trg.y() <= 0 || 
	    src.y()-pi->y() >= 0 && pi->y()-trg.y() >= 0) &&
	   (src.z()-pi->z() <= 0 && pi->z()-trg.z() <= 0 ||
	    src.z()-pi->z() >= 0 && pi->z()-trg.z() >= 0)) {
	  //	std::cerr << "true" << std::endl;
	  pprev = pi;
	  ++pi;
	  res.push_back(*pprev);
	  points.erase(pprev);
	}
	else {
	  //	std::cerr << "false" << std::endl;
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
      TRACEN("points " << *prev << "           " << *p << "      " << *next);
      
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
      
      TRACEN("sps " << sp1 << "     " << sp2);
      TRACEN(orth_coords[min] << "|" << 
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
    
  void create_vertices_of_box_with_plane(const Plane_3& h, bool b) {
    // SETDTHREAD(19*43*11);
  
    std::list<Point_3> points(find_points_of_box_with_plane(h));
    create_vertices_on_infibox(h,points,b,true,false);

    RT sum= h.a()+h.b()+h.c(); 
    if(h.d()!=0 || sum!= 0) { 
      TRACEN(sum); 
      create_extended_box_corner( 1, 1, 1, (sum<0 || (sum == 0 && h.d()<0)));
    }
    sum=-h.a()+h.b()+h.c(); 
    if(h.d()!=0 || sum!= 0) { 
      TRACEN(sum); 
      create_extended_box_corner(-1, 1, 1, (sum<0 || (sum == 0 && h.d()<0)));
    }
    sum= h.a()-h.b()+h.c(); 
    if(h.d()!=0 || sum!= 0) { 
      TRACEN(sum); 
      create_extended_box_corner( 1,-1, 1, (sum<0 || (sum == 0 && h.d()<0)));
    }
    sum=-h.a()-h.b()+h.c(); 
    if(h.d()!=0 || sum!= 0) { 
      TRACEN(sum); 
      create_extended_box_corner(-1,-1, 1, (sum<0 || (sum == 0 && h.d()<0)));
    }
    sum= h.a()+h.b()-h.c(); 
    if(h.d()!=0 || sum!= 0) { 
      TRACEN(sum); 
      create_extended_box_corner( 1, 1,-1, (sum<0 || (sum == 0 && h.d()<0)));
    }
    sum=-h.a()+h.b()-h.c(); 
    if(h.d()!=0 || sum!= 0) { 
      TRACEN(sum); 
      create_extended_box_corner(-1, 1,-1, (sum<0 || (sum == 0 && h.d()<0)));
    }
    sum= h.a()-h.b()-h.c(); 
    if(h.d()!=0 || sum!= 0) { 
      TRACEN(sum); 
      create_extended_box_corner( 1,-1,-1, (sum<0 || (sum == 0 && h.d()<0)));
    }
    sum=-h.a()-h.b()-h.c(); 
    if(h.d()!=0 || sum!= 0) { 
      TRACEN(sum); 
      create_extended_box_corner(-1,-1,-1, (sum<0 || (sum == 0 && h.d()<0)));
    }
  }

  Vertex_handle 
    create_frame_point(Point_3 p, Point_3 sp1, Point_3 sp2, Plane_3 h, 
		       const Mark& boundary, const Mark& inside, const Mark& outside) const {
    if(h.d() == 0) {
      CGAL_assertion(CGAL_NTS abs(p.hy()) != CGAL_NTS abs(p.hx()) ||
		     CGAL_NTS abs(p.hz()) != CGAL_NTS abs(p.hx()));
    }
    
    int max = 0;
    if(CGAL_NTS abs(p.hx()) > CGAL_NTS abs(p.hy()))
      max = 1;
    if(CGAL_NTS abs(p.hx()) > CGAL_NTS abs(p.hz()))
      max = 2;
    
    TRACEN("create frame point ");
    
    TRACEN("create spoints");
    Sphere_point SP[4];
    switch(max) {
    case 0: SP[2] = Sphere_point(1,0,0); break;
    case 1: SP[2] = Sphere_point(0,1,0); break;
    case 2: SP[2] = Sphere_point(0,0,1); break;
    default: CGAL_assertion_msg(0,"wrong value");
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
    
    TRACEN("create corner frame point ");
    
    RT vp[3];
    vp[0] = -p.hx()[1];
    vp[1] = -p.hy()[1];
    vp[2] = -p.hz()[1];
    
    TRACEN("create spoints");
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
    default: CGAL_assertion_msg(0,"wrong value");
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
    
    TRACEN("create degenerate corner frame point ");
    
    TRACEN("create spoints");
    Sphere_point SP[4];
    
    switch(max) { 
    case 0: SP[2] = Sphere_point(vp[0],0,0); break; // plane(x,x,0), plane(x,0,x)
    case 1: SP[2] = Sphere_point(0,vp[1],0); break; // plane(0,x,x)
    default: CGAL_assertion_msg(0,"wrong value \"max\"");
    }
    
    switch(min+max) {
    case 1: SP[3] = Sphere_point(0,0,vp[2]); break; // plane(0,x,x), plane(x,0,x)
    case 2: SP[3] = Sphere_point(0,vp[1],0); break; // plane(x,x,0)
    default: CGAL_assertion_msg(0,"wrong value \"min+max\"");
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
    
    TRACEN("create_SM_on_infibox: center = " << center);
    
    TRACEN("create svertices");
    std::vector<SVertex_handle> sv(size);
    TRACEN("create_SM_on_infibox: size = " << size);
    for(int i=0; i<size; ++i) {
      TRACEN("                      SP["<< i << "]=" << SP[i]);
      sv[i] = SD.new_svertex(SP[i]);
      mark(sv[i]) = (i < 2 ? boundary : 1);
    }
    
    TRACEN("create sedges");
    std::vector<SHalfedge_handle> she(size+1);
    for(int si=0; si<size-1;++si) {
      she[si]=SD.new_shalfedge_pair(sv[si], sv[(si+1)%(size-1)]);
      SD.circle(she[si])= 
	Sphere_circle(Plane_3(SP[si],SP[(si+1)%(size-1)],Point_3(0,0,0)));
      SD.circle(SD.twin(she[si])) =  SD.circle(she[si]).opposite();
      SD.mark(she[si]) = SD.mark(SD.twin(she[si])) = (si == 0 ? boundary : 1);
    }
    
    she[size-1] = SD.new_shalfedge_pair(she[0], sv[size-1], -1);
    CGAL_assertion(SD.next(she[size-1]) == SD.twin(she[size-1]));
    CGAL_assertion(SD.previous(SD.twin(she[size-1])) == she[size-1]);
    CGAL_assertion(SD.twin(SD.previous(SD.twin(she[size-1]))) == SD.twin(she[size-1]));
    
    she[size]   = SD.new_shalfedge_pair(SD.twin(she[size-1]), SD.twin(she[0]), 1, 1);
    
    CGAL_assertion(SD.next(she[0]) == she[1]);
    CGAL_assertion(SD.next(she[1]) == she[2]);
    if(size == 4)
      CGAL_assertion(SD.next(she[2]) == she[0]);  
    else {
      CGAL_assertion(SD.next(she[2]) == she[3]);
      CGAL_assertion(SD.next(she[3]) == she[0]);
    }
    
    CGAL_assertion(SD.next(SD.twin(she[0])) == she[size-1]);
    CGAL_assertion(SD.next(she[size-1]) == she[size]);
    CGAL_assertion(SD.next(she[size]) == SD.twin(she[0]));
    
    CGAL_assertion(SD.next(SD.twin(she[4])) == SD.twin(she[3]));
    CGAL_assertion(SD.next(SD.twin(she[3])) == SD.twin(she[2]));
    CGAL_assertion(SD.next(SD.twin(she[2])) == SD.twin(she[1]));
    if(size == 4)
      CGAL_assertion(SD.next(SD.twin(she[1])) == SD.twin(she[4]));
    else {
      CGAL_assertion(SD.next(SD.twin(she[5])) == SD.twin(she[4]));
      CGAL_assertion(SD.next(SD.twin(she[1])) == SD.twin(she[5]));
    }
    
    SD.circle(she[size-1])= Sphere_circle(Plane_3(SP[0],SP[size-1],Point_3(0,0,0)));
    SD.circle(SD.twin(she[size-1])) =  SD.circle(she[size-1]).opposite();	  
    SD.mark(she[size-1]) = SD.mark(SD.twin(she[size-1])) = 1;
    
    SD.circle(she[size])= Sphere_circle(Plane_3(SP[size-1],SP[1],Point_3(0,0,0)));
    SD.circle(SD.twin(she[size])) =  SD.circle(she[size]).opposite();
    SD.mark(she[size]) = SD.mark(SD.twin(she[size])) = 1;
    
    TRACEN("create sfaces");
    SFace_handle sf[3];
    for(int i=0; i<3; ++i)
      sf[i] = SD.new_sface();
    SD.mark(sf[0])= fmark0;
    SD.mark(sf[1])= fmark1;
    SD.mark(sf[2])=0;
    
    SD.link_as_face_cycle(she[0],sf[0]);
    SD.link_as_face_cycle(SD.twin(she[0]),sf[1]);
    SD.link_as_face_cycle(SD.twin(she[1]),sf[2]);
    
    return v;
  }

  void build_external_structure() {
    //    SETDTHREAD(19*43*131);

#ifdef CGAL_NEF3_TIMER_EXTERNAL_STRUCTURE
    CGAL::Timer timer_external_structure;
    timer_external_structure.start();
#endif

#ifdef CGAL_NEF3_TIMER_PLUECKER
    CGAL::Timer timer_pluecker;
    timer_pluecker.start();
#endif
    pair_up_halfedges();
#ifdef CGAL_NEF3_TIMER_PLUECKER
    timer_pluecker.stop();
     if(cgal_nef3_timer_on)
      std::cout << "Runtime_pluecker: " 
		<< timer_pluecker.time() << std::endl;
#endif   
     //    SNC_io_parser<SNC_structure> O0(std::cout,*this->sncp());
     //    O0.print();
    link_shalfedges_to_facet_cycles();
    categorize_facet_cycles_and_create_facets();
    create_volumes();

#ifdef CGAL_NEF3_TIMER_EXTERNAL_STRUCTURE
    timer_external_structure.stop();
    if(cgal_nef3_timer_on)
      std::cout << "Runtime_external_structure: " 
		<< timer_external_structure.time() << std::endl;
#endif
  }

  void clear_external_structure() {
    this->sncp()->clear_boundary();

    while(this->sncp()->volumes_begin()!= this->sncp()->volumes_end())
      this->sncp()->delete_volume(this->sncp()->volumes_begin());

    while(this->sncp()->halffacets_begin()!= this->sncp()->halffacets_end())
      this->sncp()->delete_halffacet_pair(this->sncp()->halffacets_begin());
  }

  void correct_infibox_sface_marks() {
    Vertex_iterator v;
    CGAL_forall_vertices(v,*this->sncp()) {
      if(is_standard(v) || Infi_box::is_infibox_corner(point(v))) continue;
      SM_decorator SD(&*v);
      SFace_iterator sf;
      CGAL_forall_sfaces(sf, SD)
	sf->volume()->mark() = sf->mark();
    }
    CGAL_forall_vertices(v,*this->sncp()) {
      if(is_standard(v) || !Infi_box::is_infibox_corner(point(v))) continue;
      SM_decorator SD(&*v);
      SFace_iterator sf;
      CGAL_forall_sfaces(sf, SD)
	sf->mark() = sf->volume()->mark();
    }    
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

}; // SNC_constructor<SNC>

CGAL_END_NAMESPACE
#endif //CGAL_SNC_CONSTRUCTOR_H

