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
#ifndef CGAL_SNC_EXTERNAL_STRUCTURE_H
#define CGAL_SNC_EXTERNAL_STRUCTURE_H

#include <CGAL/license/Nef_3.h>


#include <CGAL/Nef_S2/Normalizing.h>
#include <CGAL/Nef_3/Pluecker_line_3.h>
#include <CGAL/Nef_3/SNC_decorator.h>
#include <CGAL/Nef_3/SNC_point_locator.h>
#include <CGAL/Nef_S2/SM_point_locator.h>
#include <CGAL/Nef_3/SNC_FM_decorator.h>
#include <CGAL/Nef_3/SNC_io_parser.h>
#include <CGAL/Nef_3/SNC_indexed_items.h>
#include <CGAL/Nef_3/SNC_simplify.h>
#include <map>
#include <list>

#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 43
#include <CGAL/Nef_2/debug.h>

#include <CGAL/use.h>

namespace CGAL {

struct int_lt {
  bool operator()(const int& i1, const int& i2) const { return i1<i2; }
};
template <typename Edge_handle>
struct Halfedge_key_lt4 {

  bool operator()(const Edge_handle& e1, const Edge_handle& e2) const {
    if(CGAL::sign(e1->point().x()) != 0) {
      if(e1->source() != e2->source())
	return CGAL::compare_x(e1->source()->point(), e2->source()->point()) < 0; 
      else 
	return e1->point().x() < 0;
    }
    if(CGAL::sign(e1->point().y()) != 0) {
      if(e1->source() != e2->source())
	return CGAL::compare_y(e1->source()->point(), e2->source()->point()) < 0; 
      else 
	return e1->point().y() < 0;
    }
    if(e1->source() != e2->source())
      return CGAL::compare_z(e1->source()->point(), e2->source()->point()) < 0; 
    return e1->point().z() < 0;
  }
};

template <typename Edge_handle>
struct Halfedge_key_lt3 {

  bool operator()(const Edge_handle& e1, const Edge_handle& e2) const {
    if(e1->source() != e2->source())
      return CGAL::lexicographically_xyz_smaller(e1->source()->point(), e2->source()->point()); 
    if(CGAL::sign(e1->point().x()) != 0)
      return e1->point().x() < 0;
    if(CGAL::sign(e1->point().y()) != 0)
      return e1->point().y() < 0;
    return e1->point().z() < 0;
  }
};

template <typename Point, typename Edge>
struct Halfedge_key {
  typedef Halfedge_key<Point,Edge> Self;
  Point p; int i; Edge e;
  Halfedge_key(Point pi, int ii, Edge ei) : 
    p(pi), i(ii), e(ei) {}
  Halfedge_key(const Self& k) : p(k.p), i(k.i), e(k.e) {}
  Self& operator=(const Self& k) { p=k.p; i=k.i; e=k.e; return *this; }
  bool operator==(const Self& k) const { return p==k.p && i==k.i; }
  bool operator!=(const Self& k) const { return !operator==(k); }
};

template <typename Point, typename Edge, class Decorator>
struct Halfedge_key_lt {
  typedef Halfedge_key<Point,Edge> Key;
  typedef typename Point::R R;
  typedef typename R::Vector_3 Vector;
  typedef typename R::Direction_3 Direction;
  bool operator()( const Key& k1, const Key& k2) const { 
    if( k1.e->source() == k2.e->source())
      return (k1.i < k2.i);
    Direction l(k1.e->vector());
    if( k1.i < 0) l = -l;
    return (Direction( k2.p - k1.p) == l); 
  }
};

template <typename Point, typename Edge>
std::ostream& operator<<(std::ostream& os, 
                         const Halfedge_key<Point,Edge>& k )
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
        typedef typename R::FT     FT;

        bool flag=false;
        FT ratioa,ratiob,ratioc,ratiod;
        FT a1,b1,c1,d1,a2,b2,c2,d2;
        a1=h1.a();
        a2=h2.a();
        b1=h1.b();
        b2=h2.b();
        c1=h1.c();
        c2=h2.c();
        d1=h1.d();
        d2=h2.d();
        if(a2==0 || a1==0)
        {
            if(a2==a1) ratioa=1;
            else flag=true;
        }
        else{ratioa=a1/a2;}
        if(b2==0 || b1==0)
        {
            if(b2==b1) ratiob=ratioa;
            else flag=true;
        }
        else{ ratiob=b1/b2;}
        if(c2==0 || c1==0)
        {
            if(c2==c1) ratioc=ratioa;
            else flag=true;
        }
        else{ ratioc=c1/c2;}
        if(d2==0 || d1==0)
        {
            if(d2==d1) ratiod=ratioc;
            else flag=true;
        }
        else{ ratiod=d1/d2;}
        if(flag || !(ratioa==ratiob && ratioc==ratiod && ratioa==ratioc))
        {
            RT diff = h1.a()-h2.a();
            if ( (diff) != 0 ) return CGAL_NTS sign(diff) < 0;
            diff = h1.b()-h2.b();
            if ( (diff) != 0 ) return CGAL_NTS sign(diff) < 0;
            diff = h1.c()-h2.c();
            if ( (diff) != 0 ) return CGAL_NTS sign(diff) < 0;
            diff = h1.d()-h2.d();
            if ( (diff) != 0 ) return CGAL_NTS sign(diff) < 0;
        }
        return 0;
    }
};


struct Plane_RT_lt {
  template <typename R>
  bool operator()(const CGAL::Plane_3<R>& h1,
                  const CGAL::Plane_3<R>& h2) const
  { return (&(h1.a())) < (&(h2.a())); }
};

// ----------------------------------------------------------------------------
// SNC_external_structure 
// ----------------------------------------------------------------------------

template <typename Items_, typename SNC_structure_>
class SNC_external_structure_base : public SNC_decorator<SNC_structure_>
{ 
public:
  typedef SNC_structure_ SNC_structure;

  typedef typename SNC_structure::Infi_box                   Infi_box;
  typedef typename Infi_box::Standard_kernel                 Standard_kernel;
  typedef typename Standard_kernel::Point_3                  Standard_point_3;
  typedef typename Infi_box::NT                              NT;

  typedef CGAL::SNC_decorator<SNC_structure>                 SNC_decorator;
  typedef CGAL::SNC_point_locator<SNC_decorator>             SNC_point_locator;
  typedef CGAL::SNC_FM_decorator<SNC_structure>              FM_decorator;
  typedef CGAL::SNC_simplify<Items_, SNC_structure>          SNC_simplify;

  typedef typename SNC_structure::Sphere_map                 Sphere_map;
  typedef CGAL::SM_decorator<Sphere_map>                     SM_decorator;  
  typedef CGAL::SM_const_decorator<Sphere_map>               SM_const_decorator;
  typedef CGAL::SM_point_locator<SM_const_decorator>         SM_point_locator;

  typedef typename SNC_structure::Halfedge_iterator Halfedge_iterator;
  typedef typename SNC_structure::Halffacet_iterator Halffacet_iterator;
  typedef typename SNC_structure::Volume_iterator Volume_iterator;

  typedef typename SNC_structure::Vertex_handle Vertex_handle;
  typedef typename SNC_structure::Halfedge_handle Halfedge_handle;
  typedef typename SNC_structure::Halffacet_handle Halffacet_handle;
  typedef typename SNC_structure::Volume_handle Volume_handle;

  typedef typename SNC_structure::Halffacet_const_handle Halffacet_const_handle;
  typedef typename SNC_structure::SFace_const_handle SFace_const_handle; 

  typedef typename SNC_structure::SHalfedge_iterator SHalfedge_iterator;
  typedef typename SNC_structure::SFace_iterator SFace_iterator;
  typedef typename SNC_structure::SHalfloop_iterator SHalfloop_iterator;

  typedef typename SNC_structure::SVertex_handle SVertex_handle;
  typedef typename SNC_structure::SHalfedge_handle SHalfedge_handle;
  typedef typename SNC_structure::SFace_handle SFace_handle;
  typedef typename SNC_structure::SHalfloop_handle SHalfloop_handle;

  typedef typename SNC_structure::Object_handle Object_handle;

  typedef typename SNC_structure::SHalfedge_around_facet_circulator 
    SHalfedge_around_facet_circulator;

  typedef typename SNC_structure::Point_3 Point_3;
  typedef typename SNC_structure::Direction_3 Direction_3;
  typedef typename SNC_structure::Plane_3 Plane_3;
  typedef typename SNC_structure::Ray_3 Ray_3;

  typedef typename SNC_structure::Sphere_point Sphere_point;
  typedef typename SNC_structure::Sphere_circle Sphere_circle;

  typedef typename SM_decorator::SHalfedge_around_svertex_circulator 
    SHalfedge_around_svertex_circulator;

  SNC_point_locator* pl;

  typedef CGAL::Unique_hash_map<SFace_const_handle,unsigned int>  
                                                         Sface_shell_hash;
  typedef CGAL::Unique_hash_map<Halffacet_const_handle,unsigned int>  
                                                         Face_shell_hash;
  typedef CGAL::Unique_hash_map<SFace_const_handle,bool> SFace_visited_hash;
  typedef CGAL::Unique_hash_map<SFace_const_handle,bool> Shell_closed_hash;

  using SNC_decorator::visit_shell_objects;
  using SNC_decorator::link_as_inner_shell;
  using SNC_decorator::link_as_outer_shell;
  using SNC_decorator::link_as_prev_next_pair;
  using SNC_decorator::get_visible_facet;
  using SNC_decorator::adjacent_sface;
  using SNC_decorator::make_twins;

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
      CGAL_NEF_TRACEN("visit sf "<<h->center_vertex()->point());
      ShellSf[h]=n;
      Done[h]=true;
      if ( CGAL::lexicographically_xyz_smaller(h->center_vertex()->point(),
					       sf_min->center_vertex()->point())) 
	sf_min = h; 
    }

    void visit(Vertex_handle h) { 
      CGAL_USE(h);
      CGAL_NEF_TRACEN("visit v  "<<h->point());
    }

    void visit(Halfedge_handle h) { 
      CGAL_USE(h);
      CGAL_NEF_TRACEN("visit he "<< h->source()->point());
    }

    void visit(Halffacet_handle h) { 
      CGAL_NEF_TRACEN(h->plane()); 
      ShellF[h]=n; 
    }

    void visit(SHalfedge_handle ) {}
    void visit(SHalfloop_handle ) {}

    SFace_handle& minimal_sface() { return sf_min; }

    void increment_shell_number() { 
      CGAL_NEF_TRACEN("leaving shell "<<n);
      ++n; 
    }
  };

  SNC_external_structure_base( SNC_structure& W, SNC_point_locator* spl = NULL) 
    : SNC_decorator(W), pl(spl) {}
  /*{\Mcreate makes |\Mvar| a decorator of |W|.}*/

 public:
  //#define CGAL_NEF_NO_HALFEDGE_KEYS
#ifdef CGAL_NEF_NO_HALFEDGE_KEYS
  void pair_up_halfedges() const {
  /*{\Mop pairs all halfedge stubs to create the edges in 3-space.}*/

  //  CGAL_NEF_SETDTHREAD(43*61);
    CGAL_NEF_TRACEN(">>>>>pair_up_halfedges");
    typedef Halfedge_key_lt3<Halfedge_handle> 
      Halfedge_key_lt;
    typedef std::list<Halfedge_handle>  Halfedge_list;
    
    typedef typename Standard_kernel::Kernel_tag   Kernel_tag;
    typedef CGAL::Pluecker_line_3<Kernel_tag,Standard_kernel> Pluecker_line_3;
    typedef CGAL::Pluecker_line_lt        Pluecker_line_lt;
    typedef std::map< Pluecker_line_3, Halfedge_list, Pluecker_line_lt> 
      Pluecker_line_map;
    
    Pluecker_line_map M;
    Pluecker_line_map M2;
    Pluecker_line_map M3;
    Pluecker_line_map M4;
    
    NT eval(Infi_box::compute_evaluation_constant_for_halfedge_pairup(*this->sncp()));
    
    Halfedge_iterator e;
    CGAL_forall_halfedges(e,*this->sncp()) {
      //    progress++;
      Point_3 p = e->source()->point();
      Point_3 q = p + e->vector();
      CGAL_NEF_TRACE(" segment("<<p<<", "<<q<<")"<<
	    " direction("<<e->vector()<<")");
      Standard_point_3 sp = Infi_box::standard_point(p,eval);
      Standard_point_3 sq = Infi_box::standard_point(q,eval);
      Pluecker_line_3  l( sp, sq);
      
      int inverted;
      l = categorize( l, inverted);
      
      if(Infi_box::is_edge_on_infibox(e))
	if(Infi_box::is_type4(e))
	  M4[l].push_back(e);
	else
	  if(Infi_box::is_type3(e))
	    M3[l].push_back(e);
	  else
	    M2[l].push_back(e);
      else
	M[l].push_back(e);
      
      // the following trace crashes when compiling with optimizations (-O{n})
      //CGAL_NEF_TRACEN(Infi_box::standard_point(point(vertex(e)))+
      
      CGAL_NEF_TRACEN(" line("<<l<<")"<<" inverted="<<inverted);
    }
    
    typename Pluecker_line_map::iterator it;
    
    CGAL_forall_iterators(it,M4) {
      //    progress++;
      it->second.sort(Halfedge_key_lt());
      CGAL_NEF_TRACEN("search opposite  "<<it->first); 
      typename Halfedge_list::iterator itl;
      CGAL_forall_iterators(itl,it->second) {
	Halfedge_handle e1 = *itl;
	++itl; 
	CGAL_assertion(itl != it->second.end());
	Halfedge_handle e2 = *itl;
	while(normalized(e1->vector()) != normalized(-e2->vector())) {
	  ++itl;
	  make_twins(e1,*itl);
	  e1 = e2;
	  ++itl;
	  e2 = *itl;
	}
	CGAL_NEF_TRACEN("    " << e1->source()->point() 
			<< " -> " << e2->source()->point());
	CGAL_NEF_TRACEN(e1->vector()<<" -> "<<-e2->vector());
	make_twins(e1,e2);
	CGAL_assertion(e1->mark()==e2->mark());
	
	// discard temporary sphere_point ?
      }
    }

    CGAL_forall_iterators(it,M3) {
      //    progress++;
      it->second.sort(Halfedge_key_lt());
      CGAL_NEF_TRACEN("search opposite  "<<it->first); 
      typename Halfedge_list::iterator itl;
      CGAL_forall_iterators(itl,it->second) {
	Halfedge_handle e1 = *itl;
	++itl; 
	CGAL_assertion(itl != it->second.end());
	Halfedge_handle e2 = *itl;
	while(normalized(e1->vector()) != normalized(-e2->vector())) {
	  ++itl;
	  make_twins(e1,*itl);
	  e1 = e2;
	  ++itl;
	  e2 = *itl;
	}
	CGAL_NEF_TRACEN("    " << e1->source()->point() 
			<< " -> " << e2->source()->point());
	CGAL_NEF_TRACEN(e1->vector()<<" -> "<<-e2->vector());
	make_twins(e1,e2);
	CGAL_assertion(e1->mark()==e2->mark());
	
	// discard temporary sphere_point ?
      }
    }    

    CGAL_forall_iterators(it,M2) {
      //    progress++;
      it->second.sort(Halfedge_key_lt());
      CGAL_NEF_TRACEN("search opposite  "<<it->first); 
      typename Halfedge_list::iterator itl;
      CGAL_forall_iterators(itl,it->second) {
	Halfedge_handle e1 = *itl;
	++itl; 
	CGAL_assertion(itl != it->second.end());
	Halfedge_handle e2 = *itl;
	while(normalized(e1->vector()) != normalized(-e2->vector())) {
	  ++itl;
	  this->make_twins(e1,*itl);
	  e1 = e2;
	  ++itl;
	  e2 = *itl;
	}
	CGAL_NEF_TRACEN("    " << e1->source()->point() 
			<< " -> " << e2->source()->point());
	CGAL_NEF_TRACEN(e1->vector()<<" -> "<<-e2->vector());
	this->make_twins(e1,e2);
	CGAL_assertion(e1->mark()==e2->mark());
	
	// discard temporary sphere_point ?
      }
    }
    
    CGAL_forall_iterators(it,M) {
      //    progress++;
      it->second.sort(Halfedge_key_lt());
      CGAL_NEF_TRACEN("search opposite  "<<it->first); 
      typename Halfedge_list::iterator itl;
      CGAL_forall_iterators(itl,it->second) {
	Halfedge_handle e1 = *itl;
	++itl; 
	CGAL_assertion(itl != it->second.end());
	Halfedge_handle e2 = *itl;
	while(normalized(e1->vector()) != normalized(-e2->vector())) {
	  ++itl;
	  this->make_twins(e1,*itl);
	  e1 = e2;
	  ++itl;
	  e2 = *itl;
	}
	CGAL_NEF_TRACEN("    " << e1->source()->point() 
			<< " -> " << e2->source()->point());
	CGAL_NEF_TRACEN(e1->vector()<<" -> "<<-e2->vector());
	this->make_twins(e1,e2);
	CGAL_assertion(e1->mark()==e2->mark());
	
	// discard temporary sphere_point ?
      }
    }

  }
#else
  void pair_up_halfedges() const {
  /*{\Mop pairs all halfedge stubs to create the edges in 3-space.}*/

//    CGAL_NEF_SETDTHREAD(43*61);
    CGAL_NEF_TRACEN(">>>>>pair_up_halfedges");
    typedef Halfedge_key< Point_3, Halfedge_handle>
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
    
    NT eval(Infi_box::compute_evaluation_constant_for_halfedge_pairup(*this->sncp()));;
    
    Halfedge_iterator e;
    CGAL_forall_halfedges(e,*this->sncp()) {
      //    progress++;
      Point_3 p = e->source()->point();
      Point_3 q = p + e->vector();
      CGAL_NEF_TRACE(" segment("<<p<<", "<<q<<")"<<
	    " direction("<<e->vector()<<")");
      Standard_point_3 sp = Infi_box::standard_point(p,eval);
      Standard_point_3 sq = Infi_box::standard_point(q,eval);
      Pluecker_line_3  l( sp, sq);
      
      int inverted;
      l = categorize( l, inverted);
      
      if(Infi_box::is_edge_on_infibox(e))
	if(Infi_box::is_type4(e))
	  M4[l].push_back(Halfedge_key(p,inverted,e));
	else
	  if(Infi_box::is_type3(e))
	    M3[l].push_back(Halfedge_key(p,inverted,e));
	  else
	    M2[l].push_back(Halfedge_key(p,inverted,e));
      else
	M[l].push_back(Halfedge_key(p,inverted,e));
      
      // the following trace crashes when compiling with optimizations (-O{n})
      //CGAL_NEF_TRACEN(Infi_box::standard_point(point(vertex(e)))+
      
      CGAL_NEF_TRACEN(" line("<<l<<")"<<" inverted="<<inverted);
    }
    
    typename Pluecker_line_map::iterator it;
    
    CGAL_forall_iterators(it,M4) {
      //    progress++;
      it->second.sort(Halfedge_key_lt());
      CGAL_NEF_TRACEN("search opposite  "<<it->first); 
      typename Halfedge_list::iterator itl;
      CGAL_forall_iterators(itl,it->second) {
	Halfedge_handle e1 = itl->e;
	++itl; 
	CGAL_assertion(itl != it->second.end());
	Halfedge_handle e2 = itl->e;
	CGAL_NEF_TRACEN("    " << e1->source()->point() 
			<< " -> " << e2->source()->point());
	CGAL_NEF_TRACEN(e1->vector()<<" -> "<<-e2->vector());
	make_twins(e1,e2);
	CGAL_assertion(e1->mark()==e2->mark());
	
	// discard temporary sphere_point ?
      }
    }
    
    CGAL_forall_iterators(it,M3) {
      //    progress++;
      it->second.sort(Halfedge_key_lt());
      CGAL_NEF_TRACEN("search opposite  "<<it->first); 
      typename Halfedge_list::iterator itl;
      CGAL_forall_iterators(itl,it->second) {
	Halfedge_handle e1 = itl->e;
	++itl; 
	CGAL_assertion(itl != it->second.end());
	Halfedge_handle e2 = itl->e;
	CGAL_NEF_TRACEN("    " << e1->source()->point() 
			<< " -> " << e2->source()->point());
	CGAL_NEF_TRACEN(e1->vector()<<" -> "<<-e2->vector());
	make_twins(e1,e2);
	CGAL_assertion(e1->mark()==e2->mark());
	
	// discard temporary sphere_point ?
      }
    }
    
    CGAL_forall_iterators(it,M2) {
      //    progress++;
      it->second.sort(Halfedge_key_lt());
      CGAL_NEF_TRACEN("search opposite  "<<it->first); 
      typename Halfedge_list::iterator itl;
      CGAL_forall_iterators(itl,it->second) {
	Halfedge_handle e1 = itl->e;
	++itl; 
	CGAL_assertion(itl != it->second.end());
	Halfedge_handle e2 = itl->e;
	CGAL_NEF_TRACEN("    " << e1->source()->point() 
			<< " -> " << e2->source()->point());
	CGAL_NEF_TRACEN(e1->vector()<<" -> "<<-e2->vector());
	make_twins(e1,e2);
	CGAL_assertion(e1->mark()==e2->mark());
	
	// discard temporary sphere_point ?
      }
    }
    
    CGAL_forall_iterators(it,M) {
      //    progress++;
      it->second.sort(Halfedge_key_lt());
      CGAL_NEF_TRACEN("search opposite  "<<it->first); 
      typename Halfedge_list::iterator itl;
      CGAL_forall_iterators(itl,it->second) {
	Halfedge_handle e1 = itl->e;
	++itl; 
	CGAL_assertion(itl != it->second.end());
	Halfedge_handle e2 = itl->e;
	CGAL_NEF_TRACEN("    " << e1->source()->point() 
			<< " -> " << e2->source()->point());
	CGAL_NEF_TRACEN(e1->vector()<<" -> "<< -e2->vector());
	CGAL_assertion(e1->source()->point() != e2->source()->point());
	CGAL_assertion(e1->mark()==e2->mark());
	make_twins(e1,e2);
	
	// discard temporary sphere_point ?
      }
    }
  }
#endif

  void link_shalfedges_to_facet_cycles() const {
  /*{\Mop creates all non-trivial facet cycles from sedges. 
  \precond |pair_up_halfedges()| was called before.}*/
    
  //  CGAL_NEF_SETDTHREAD(43*31);
    CGAL_NEF_TRACEN(">>>>>link_shalfedges_to_facet_cycles");

#ifdef CGAL_NEF_EXPLOIT_REFERENCE_COUNTING
    Point_3 p1(1,2,7), p2(p1);     
    bool reference_counted = (&(p1.hx()) == &(p2.hx()));
#endif

    Halfedge_iterator e;
    CGAL_forall_edges(e,*this->sncp()) {
      //    progress++;
      CGAL_NEF_TRACEN("");
      CGAL_NEF_TRACEN(PH(e));
      Halfedge_iterator et = e->twin();
      SM_decorator D(&*e->source()), Dt(&*et->source());
      CGAL_NEF_TRACEN(e->source()->point());
      if ( D.is_isolated(e) ) continue;
      SHalfedge_around_svertex_circulator ce(D.first_out_edge(e)),cee(ce);
      SHalfedge_around_svertex_circulator cet(Dt.first_out_edge(et)),cete(cet);
      
#ifdef CGAL_NEF_EXPLOIT_REFERENCE_COUNTING
      if(reference_counted) {
	CGAL_For_all(cet,cete)
	  if ( &(cet->circle().a()) == &(ce->circle().opposite().a()) && 
	       cet->source()->twin() == ce->source() )
	    break;
      } else
#endif
        CGAL_For_all(cet,cete)
	  if ( cet->circle() == ce->circle().opposite() && 
	       cet->source()->twin() == ce->source() ) 
            break;

#ifndef NDEBUG
      if( cet->circle() != ce->circle().opposite() )
	CGAL_NEF_TRACEN("assertion failed!");
      
      CGAL_NEF_TRACEN("vertices " << e->source()->point() << 
      "    "      << et->source()->point());
      
      
      SHalfedge_around_svertex_circulator sc(D.first_out_edge(e));
      SHalfedge_around_svertex_circulator sct(Dt.first_out_edge(et));

      CGAL_NEF_TRACEN("");
      CGAL_For_all(sc,cee)
	CGAL_NEF_TRACEN("sseg@E addr="<<&*sc<<
			" src="<< sc->source()->point()<<
			" tgt="<< sc->target()->point()<<std::endl<<
			" circle=" << normalized(sc->circle()));
      CGAL_NEF_TRACEN("");

      CGAL_For_all(sct,cete)
      CGAL_NEF_TRACEN("sseg@ET addr="<<&*sct<<
		      " src="<< sct->source()->point()<<
		      " tgt="<<sct->target()->point()<<std::endl<<
		      " circle=" << normalized(sct->circle()));
      CGAL_NEF_TRACEN("");
#endif

      CGAL_assertion( normalized(cet->circle()) == normalized(ce->circle().opposite()) ); 
      CGAL_assertion( cet->source()->twin() == ce->source()); 
      CGAL_For_all(ce,cee) { 
	CGAL_NEF_TRACEN("circles " << normalized(cet->circle()) << "   " << normalized(ce->circle()) << 
			" sources " << cet->target()->point() << 
			"   " << ce->target()->point());
	CGAL_assertion( normalized(cet->circle()) == normalized(ce->circle().opposite())); 
	CGAL_assertion( cet->source()->twin() == ce->source()); 
	CGAL_assertion(ce->mark()==cet->mark());
	link_as_prev_next_pair(cet->twin(),ce);
	link_as_prev_next_pair(ce->twin(),cet);
	--cet; // ce moves ccw, cet moves cw
      }
    }
  }

  void categorize_facet_cycles_and_create_facets() const {
  /*{\Mop collects all facet cycles incident to a facet and creates
  the facets. \precond |link_shalfedges_to_facet_cycles()| was called
  before.}*/

    //    CGAL_NEF_SETDTHREAD(43*31);
    CGAL_NEF_TRACEN(">>>>>categorize_facet_cycles_and_create_facets");
    
    typedef std::list<Object_handle> Object_list;
#ifdef CGAL_NEF_EXPLOIT_REFERENCE_COUNTING
    typedef std::map<Plane_3, Object_list, Plane_RT_lt> 
      Map_planes;
#else
    typedef std::map<Plane_3, Object_list, Plane_lt> 
      Map_planes;
#endif 
   
    Map_planes M;
    SHalfedge_iterator e;
    CGAL_forall_shalfedges(e,*this->sncp()) {
      Sphere_circle c(e->circle());
      Plane_3 h = c.plane_through(e->source()->source()->point());
      CGAL_NEF_TRACEN("\n" << e->source()->twin()->source()->point() <<" - "
		      << e->source()->source()->point() <<" - "<< 
		      e->twin()->source()->twin()->source()->point() << 
		      " has plane " << h << " has circle " << e->circle() << 
		      " has signum " << sign_of(h));
      if ( sign_of(h)<0 ) continue;
      M[normalized(h)].push_back(make_object(e->twin())); 
      CGAL_NEF_TRACEN(" normalized as " << normalized(h));
      /*
	Unique_hash_map<SHalfedge_handle, bool> Done(false);
	SHalfedge_iterator ei;
	CGAL_forall_sedges(ei,*this->sncp()) {
	if(Done[ei]) continue;
	Sphere_circle c(ei->circle());
	Plane_3 h = c.plane_through(ei->source()->source()->point());
	CGAL_NEF_TRACEN("\n" << ei->source()->twin()->source()->point() <<" - "
	<< ei->source()->source()->point() <<" - "<< 
		      ei->twin()->source()->twin()->source()->point() << 
		      " has plane " << h << " has circle " << ei->circle() << 
		      " has signum " << sign_of(h));
      SHalfedge_handle e(ei);
      if ( sign_of(h)<0 ) {
	h = h.opposite();
	e = e->twin();
      }
      SHalfedge_around_facet_circulator sfc(e), send(sfc);
      CGAL_For_all(sfc, send) {
	M[normalized(h)].push_back(make_object(e->twin())); 
	Done[sfc] = true;
	Done[sfc->twin()] = true;
	CGAL_NEF_TRACEN(" normalized as " << normalized(h)); 
      }
      */
    }
    SHalfloop_iterator l;
    CGAL_forall_shalfloops(l,*this->sncp()) {
      Sphere_circle c(l->circle());
      Plane_3 h = c.plane_through(l->incident_sface()->center_vertex()->point()); 
      if ( sign_of(h)<0 ) continue;
      // CGAL_assertion( h == normalized(h));
      M[normalized(h)].push_back(make_object(l->twin()));
    }
    
#ifdef CGAL_NEF3_TIMER_PLANE_SWEEPS
  number_of_plane_sweeps=0;
  timer_plane_sweeps.reset();
#endif

    typename Map_planes::iterator it;
    CGAL_forall_iterators(it,M) { 
      //    progress2++;
      //      CGAL_NEF_TRACEN("  plane "<<it->first<<"             "<<(it->first).point());
      FM_decorator D(*this->sncp());
      D.create_facet_objects(it->first,it->second.begin(),it->second.end());
    }
    //       CGAL_NEF_SETDTHREAD(1);
  }

  void create_volumes() {
  /*{\Mop collects all shells incident to a volume and creates the
  volumes.  \precond |categorize_facet_cycles_and_creating_facets()| was
  called before.}*/

#ifdef CGAL_NEF3_TIMER_POINT_LOCATION
    number_of_ray_shooting_queries=0;
    timer_ray_shooting.reset();
#endif 

    //    CGAL_NEF_SETDTHREAD(37*43*503*509);

    CGAL_NEF_TRACEN(">>>>>create_volumes");
    Sface_shell_hash     ShellSf(0);
    Face_shell_hash      ShellF(0);
    SFace_visited_hash Done(false);
    Shell_explorer V(*this,ShellSf,ShellF,Done);
    std::vector<SFace_handle> MinimalSFace;
    std::vector<SFace_handle> EntrySFace;
    std::vector<bool> Closed;
    
    SFace_iterator f;
    // First, we classify all the Shere Faces per Shell.  For each Shell we
    //     determine its minimum lexicographyly vertex and we check wheter the
    //     Shell encloses a region (closed surface) or not. 
    CGAL_forall_sfaces(f,*this->sncp()) {
      //    progress++;
      CGAL_NEF_TRACEN("sface in " << ShellSf[f]);
      if ( Done[f] ) 
	continue;
      V.minimal_sface() = f;
      visit_shell_objects(f,V);
      
      CGAL_NEF_TRACEN("minimal vertex " << V.minimal_sface()->center_vertex()->point());

      MinimalSFace.push_back(V.minimal_sface());
      EntrySFace.push_back(f);
      V.increment_shell_number();
      CGAL_NEF_TRACEN("sface out " << ShellSf[f]);
    }
    
    for(unsigned int i=0; i<EntrySFace.size(); ++i)
      Closed.push_back(false);
    
    Halffacet_iterator hf;
    CGAL_forall_facets(hf,*this) 
      if(ShellF[hf] != ShellF[hf->twin()]) {
	Closed[ShellF[hf]] = true;
	Closed[ShellF[hf->twin()]] = true;
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

    //   then, we determine the Shells which correspond to Volumes via a ray
    //   shootting in the direction (-1,0,0) over the Sphere_map of the minimal 
    //   vertex.  The Shell corresponds to a Volume if the object hit belongs 
    //   to another Shell. 
    
    this->sncp()->new_volume(); // outermost volume (nirvana)
    if(MinimalSFace.size() == 0) return;
    Vertex_handle v_min = MinimalSFace[0]->center_vertex();
    for( unsigned int i = 0; i < MinimalSFace.size(); ++i) {
      //    progress2++;
      Vertex_handle v = MinimalSFace[i]->center_vertex();
      if(CGAL::lexicographically_xyz_smaller(v->point(),v_min->point()))
	v_min=v;
      CGAL_NEF_TRACEN( "Shell #" << i << " minimal vertex: " << v->point());
      SM_point_locator D((Sphere_map*) &*v);
      Object_handle o = D.locate(Sphere_point(-1,0,0));
      SFace_const_handle sfc;
      if( !CGAL::assign(sfc, o) || ShellSf[sfc] != i) {
	CGAL_NEF_TRACEN("the shell encloses a volume");
	CGAL_NEF_TRACEN("sface hit? "<<CGAL::assign(sfc,o));
	if( CGAL::assign(sfc, o)) { CGAL_NEF_TRACEN("sface on another shell? "<<ShellSf[sfc]); }
	SFace_handle f = MinimalSFace[i];
	CGAL_assertion( ShellSf[f] == i );
	if( Closed[i] ) {
	  CGAL_NEF_TRACEN("Shell #" << i << " is closed");
	  SM_decorator SD(&*v);
	  Volume_handle c = this->sncp()->new_volume();
	  c->mark() = f->mark(); // TODO test if line is redundant
	  link_as_inner_shell(f, c );
	  CGAL_NEF_TRACE( "Shell #" << i <<" linked as inner shell");
	  CGAL_NEF_TRACEN( "(sface" << (CGAL::assign(sfc,o)?"":" not") << " hit case)");
	}
      }
    }
    
    // finaly, we go through all the Shells which do not correspond to a Volume 
    //     and we assign them to its enclosing Volume determined via a facet below
    //     check. 
    
    CGAL_forall_sfaces(f,*this->sncp()) {
      //    progress3++;
      if ( f->volume() != Volume_handle() ) 
	continue;
      CGAL_NEF_TRACEN( "Outer shell #" << ShellSf[f] << " volume?");
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
    Point_3 p = vi->point();
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
    CGAL_NEF_TRACEN("get_facet_below");
    if( CGAL::assign(v, o)) {
      CGAL_NEF_TRACEN("facet below from from vertex...");
      f_below = get_visible_facet(v, ray);
      if( f_below == Halffacet_handle()) {
	CGAL_assertion(v->sfaces_begin() == v->sfaces_last());
	f_below = get_facet_below(MinimalSFace[Shell[v->sfaces_begin()]]->center_vertex(), 
				  MinimalSFace,Shell);
      }
    }
    else if( CGAL::assign(e, o)) {
      CGAL_NEF_TRACEN("facet below from from edge...");
      f_below = get_visible_facet(e, ray);
      if( f_below == Halffacet_handle()) {
	CGAL_assertion(e->source()->sfaces_begin() == e->source()->sfaces_last());
	f_below = get_facet_below(MinimalSFace[Shell[e->source()->sfaces_begin()]]->center_vertex(), 
				  MinimalSFace, Shell);
      }
    }
    else if( CGAL::assign(f, o)) {
      CGAL_NEF_TRACEN("facet below from from facet...");
      f_below = get_visible_facet(f, ray);
      CGAL_assertion( f_below != Halffacet_handle());
    }
    else { CGAL_NEF_TRACEN("no facet below found..."); }
    return f_below;
  }  

  Volume_handle determine_volume( SFace_handle sf, 
                const std::vector< SFace_handle>& MinimalSFace, 
				  const Sface_shell_hash&  Shell ) const {
    //{\Mop determines the volume |C| that a shell |S| pointed by |sf| 
    //  belongs to.  \precondition |S| separates the volume |C| from an enclosed
    //  volume.}
    
    CGAL_NEF_TRACEN("determine volume");
    Vertex_handle v_min = MinimalSFace[Shell[sf]]->center_vertex();
    
    Halffacet_handle f_below = get_facet_below(v_min, MinimalSFace, Shell);
    if ( f_below == Halffacet_handle())
      return SNC_decorator(*this).volumes_begin();
    Volume_handle c = f_below->incident_volume();
    if( c != Volume_handle()) {
      CGAL_NEF_TRACE( "Volume " << &*c << " hit ");
      CGAL_NEF_TRACEN("(Shell #" << Shell[adjacent_sface(f_below)] << ")");
      return c;
    }
    SFace_handle sf_below = adjacent_sface(f_below);
    CGAL_NEF_TRACE( "Shell not assigned to a volume hit ");
    CGAL_NEF_TRACEN( "(Inner shell #" << Shell[sf_below] << ")");
    c = determine_volume( sf_below, MinimalSFace, Shell);
    link_as_inner_shell( sf_below, c);
    return c;
  }
   
  void build_external_structure() {

#ifdef CGAL_NEF3_TIMER_EXTERNAL_STRUCTURE
    CGAL::Timer timer_external_structure;
    timer_external_structure.start();
#endif

#ifdef CGAL_NEF3_TIMER_PLUECKER
    CGAL::Timer timer_pluecker;
    timer_pluecker.start();
#endif
    //    SNC_io_parser<SNC_structure> O0(std::cerr,*this->sncp());
    //    O0.print();
    pair_up_halfedges();
#ifdef CGAL_NEF3_TIMER_PLUECKER
    timer_pluecker.stop();
     if(cgal_nef3_timer_on)
      std::cout << "Runtime_pluecker: " 
		<< timer_pluecker.time() << std::endl;
#endif
    link_shalfedges_to_facet_cycles();
    //    SNC_io_parser<SNC_structure> O0(std::cerr,*this->sncp());
    //    O0.print();
    categorize_facet_cycles_and_create_facets();
    create_volumes();

#ifdef CGAL_NEF3_TIMER_EXTERNAL_STRUCTURE
    timer_external_structure.stop();
    if(cgal_nef3_timer_on)
      std::cout << "Runtime_external_structure: " 
		<< timer_external_structure.time() << std::endl;
#endif
  }

  template<typename Association>
  void build_after_binary_operation(Association) {

#ifdef CGAL_NEF3_TIMER_SIMPLIFICATION
    CGAL::Timer timer_simplification;
    timer_simplification.start();
#endif
    SNC_simplify simp(*this->sncp());
    simp.vertex_simplification(SNC_simplify::NO_SNC);
#ifdef CGAL_NEF3_TIMER_SIMPLIFICATION
    timer_simplification.stop();
    if(cgal_nef3_timer_on)
      std::cout << "Runtime_simplification: " 
		<< timer_simplification.time() << std::endl;
#endif
    
    CGAL_NEF_TRACEN("\nnumber of vertices (so far...) = "
		    << this->sncp()->number_of_vertices());
    
    CGAL_NEF_TRACEN("=> resultant vertices (after simplification): ");
    
    build_external_structure();
  }

  void clear_external_structure() {
    this->sncp()->clear_snc_boundary();

    while(this->sncp()->volumes_begin()!= this->sncp()->volumes_end())
      this->sncp()->delete_volume(this->sncp()->volumes_begin());

    while(this->sncp()->halffacets_begin()!= this->sncp()->halffacets_end())
      this->sncp()->delete_halffacet_pair(this->sncp()->halffacets_begin());

    SHalfedge_iterator se;
    CGAL_forall_shalfedges(se,*this)
      se->facet() = Halffacet_handle();

    SFace_iterator sf;
    CGAL_forall_sfaces(sf,*this)
      sf->volume() = Volume_handle();
  }
}; 



template <typename Items_, typename SNC_structure_>
class SNC_external_structure : public SNC_external_structure_base<Items_, SNC_structure_>
{ 
  typedef CGAL::SNC_decorator<SNC_structure_>                 SNC_decorator;
  typedef CGAL::SNC_point_locator<SNC_decorator>             SNC_point_locator;
public:

  SNC_external_structure( SNC_structure_& W, SNC_point_locator* spl = NULL) 
    : SNC_external_structure_base<Items_, SNC_structure_>(W, spl) 
  {}
};




template <typename SNC_structure_>
class SNC_external_structure<SNC_indexed_items, SNC_structure_> 
  : public SNC_external_structure_base<int, SNC_structure_> {
 
public:
  typedef SNC_structure_                                SNC_structure;
  typedef SNC_external_structure_base<int, SNC_structure>    Base;

  typedef CGAL::SNC_decorator<SNC_structure>                 SNC_decorator;
  typedef CGAL::SNC_point_locator<SNC_decorator>             SNC_point_locator;
  typedef CGAL::SNC_FM_decorator<SNC_structure>              FM_decorator;
  typedef typename SNC_structure::Sphere_map                 Sphere_map;
  typedef CGAL::SM_decorator<Sphere_map>                     SM_decorator;  
  typedef CGAL::SNC_simplify<SNC_indexed_items, SNC_structure> SNC_simplify;

  typedef typename SNC_structure::Halfedge_iterator Halfedge_iterator;
  typedef typename SNC_structure::SHalfedge_iterator SHalfedge_iterator;
  typedef typename SNC_structure::SHalfloop_iterator SHalfloop_iterator;

  typedef typename SNC_structure::Halfedge_handle Halfedge_handle;
  typedef typename SNC_structure::SHalfedge_handle SHalfedge_handle;

  typedef typename SNC_structure::Object_handle Object_handle;

  typedef typename SNC_structure::SHalfedge_around_facet_circulator 
    SHalfedge_around_facet_circulator;
  typedef typename SM_decorator::SHalfedge_around_svertex_circulator 
    SHalfedge_around_svertex_circulator;
  
  typedef typename SNC_structure::Plane_3 Plane_3;

  using Base::make_twins;
  using Base::link_as_prev_next_pair;
  using Base::link_as_inner_shell;


  SNC_external_structure( SNC_structure& W, SNC_point_locator* spl = NULL) 
    : Base(W, spl) {}
  /*{\Mcreate makes |\Mvar| a decorator of |W|.}*/

 public:
  void pair_up_halfedges() const {
    typedef Halfedge_key_lt4<Halfedge_handle>  Halfedge_key_lt;
    typedef std::list<Halfedge_handle>  Halfedge_list;
    typedef std::map<int, Halfedge_list, int_lt> index_map;

    CGAL_NEF_TRACEN("pair up by indexes");

    index_map i2he;
    Halfedge_iterator ei;
    CGAL_forall_halfedges(ei, *this->sncp()) 
      i2he[ei->get_index()].push_back(ei);
    
    typename index_map::iterator it;
    CGAL_forall_iterators(it,i2he) {
      CGAL_NEF_TRACEN("pair up " << it->first);
      it->second.sort(Halfedge_key_lt());
      typename Halfedge_list::iterator itl;
      CGAL_forall_iterators(itl,it->second) {
	Halfedge_handle e1 = *itl;
	CGAL_NEF_TRACEN(e1->source()->point() << ", " << e1->vector());
	++itl; 
	CGAL_assertion(itl != it->second.end());
	Halfedge_handle e2 = *itl;
	CGAL_NEF_TRACEN(" + " << e2->source()->point() << ", " << e2->vector());
	make_twins(e1,e2);
	//	SNC_io_parser<SNC_structure> O0(std::cerr,*this->sncp());
	//	O0.print();
	CGAL_assertion(e1->mark()==e2->mark());
      }
    }
  }

  void link_shalfedges_to_facet_cycles() const {
  /*{\Mop creates all non-trivial facet cycles from sedges. 
  \precond |pair_up_halfedges()| was called before.}*/
    
  //  CGAL_NEF_SETDTHREAD(43*31);
    CGAL_NEF_TRACEN(">>>>>link_shalfedges_to_facet_cycles");

    Halfedge_iterator e;
    CGAL_forall_edges(e,*this->sncp()) {
      //    progress++;
      CGAL_NEF_TRACEN("");
      CGAL_NEF_TRACEN(PH(e));
      Halfedge_iterator et = e->twin();
      SM_decorator D(&*e->source()), Dt(&*et->source());
      CGAL_NEF_TRACEN(e->source()->point());
      if ( D.is_isolated(e) ) continue;
      SHalfedge_around_svertex_circulator ce(D.first_out_edge(e)),cee(ce);
      SHalfedge_around_svertex_circulator cet(Dt.first_out_edge(et)),cete(cet);
      
      CGAL_For_all(cet,cete) {
	//	std::cerr << cet->get_index() << ", " << ce->twin()->get_index() << std::endl;
	if (cet->get_forward_index() == ce->twin()->get_backward_index())
	    //	    cet->source()->twin() == ce->source())
	  break;
      }

      CGAL_NEF_TRACEN("vertices " << e->source()->point() << 
      "    "      << et->source()->point());
      
      SHalfedge_around_svertex_circulator sc(D.first_out_edge(e));
      SHalfedge_around_svertex_circulator sct(Dt.first_out_edge(et));

      CGAL_NEF_TRACEN("");
      CGAL_For_all(sc,cee)
	CGAL_NEF_TRACEN("sseg@E addr="<<&*sc<<
			" src="<< sc->source()->point()<< 
			" tgt="<< sc->target()->point()<< std::endl << 
			" circle=" << sc->circle()<< std::endl <<
			" indexes=" << sc->get_forward_index() << 
			"," << sc->get_backward_index() << std::endl <<
			"         " << sc->twin()->get_forward_index() <<
			"," << sc->twin()->get_backward_index());
 
      CGAL_NEF_TRACEN("");

      CGAL_For_all(sct,cete)
      CGAL_NEF_TRACEN("sseg@ET addr="<<&*sct<<
		      " src="<< sct->source()->point()<<
		      " tgt="<<sct->target()->point() <<std::endl<<
		      " circle=" << sct->circle() << std::endl<<
		      " indexes=" << sct->get_forward_index() << 
		      "," << sct->get_backward_index() << std::endl <<
		      "         " << sct->twin()->get_forward_index() <<
	              "," << sct->twin()->get_backward_index());
      CGAL_NEF_TRACEN("");

      CGAL_assertion( cet->get_index() == ce->twin()->get_index());
      CGAL_assertion( cet->twin()->get_index() == ce->get_index());
      CGAL_assertion( cet->source()->twin() == ce->source()); 
      CGAL_For_all(ce,cee) { 
	CGAL_NEF_TRACEN("circles " << cet->circle() << "   " << ce->circle() << 
			" sources " << cet->target()->point() << 
			"   " << ce->target()->point());
	CGAL_assertion( cet->source()->twin() == ce->source()); 
	CGAL_assertion(ce->mark()==cet->mark());
	link_as_prev_next_pair(cet->twin(),ce);
	link_as_prev_next_pair(ce->twin(),cet);
	--cet; // ce moves ccw, cet moves cw
      }
    }
  }

  void categorize_facet_cycles_and_create_facets() const {
  /*{\Mop collects all facet cycles incident to a facet and creates
  the facets. \precond |link_shalfedges_to_facet_cycles()| was called
  before.}*/

    //    CGAL_NEF_SETDTHREAD(43*31);
    CGAL_NEF_TRACEN(">>>>>categorize_facet_cycles_and_create_facets");
    
    typedef std::list<Object_handle> Object_list;
    typedef std::map<int, Object_list> 
      Map_planes;
   
    Map_planes M;
    SHalfedge_iterator e;
    CGAL_forall_shalfedges(e,*this->sncp()) {
      if(e->get_index() > e->twin()->get_index())
	continue;
      M[e->get_index()].push_back(make_object(e));
    }
    SHalfloop_iterator l;
    CGAL_forall_shalfloops(l,*this->sncp()) {
      if(l->get_index() > l->twin()->get_index())
	continue;
      M[l->get_index()].push_back(make_object(l));
    }
    
#ifdef CGAL_NEF3_TIMER_PLANE_SWEEPS
  number_of_plane_sweeps=0;
  timer_plane_sweeps.reset();
#endif

    typename Map_planes::iterator it;
    CGAL_forall_iterators(it,M) { 
      CGAL_NEF_TRACEN("  plane "<< it->first);
      CGAL_NEF_TRACEN("  size "<< it->second.size());
      FM_decorator D(*this->sncp());
      Plane_3 h;
      Object_handle o(*it->second.begin());
      if(CGAL::assign(e, o))
	h = e->circle().opposite().plane_through(e->source()->source()->point());
      else if(CGAL::assign(l, o))
	h = l->circle().opposite().plane_through(l->incident_sface()->center_vertex()->point());
      else
	CGAL_error_msg( "wrong handle");

      D.create_facet_objects(h,it->second.begin(),it->second.end());
    }
    //       CGAL_NEF_SETDTHREAD(1);
  }

  void create_volumes() {
    Base::create_volumes();
  }
   
  void build_external_structure() {
//    CGAL_NEF_SETDTHREAD(503*509);

#ifdef CGAL_NEF3_TIMER_EXTERNAL_STRUCTURE
    CGAL::Timer timer_external_structure;
    timer_external_structure.start();
#endif

#ifdef CGAL_NEF3_TIMER_PLUECKER
    CGAL::Timer timer_pluecker;
    timer_pluecker.start();
#endif
    //    SNC_io_parser<SNC_structure> O0(std::cerr,*this->sncp());
    //    O0.print();
    pair_up_halfedges();
#ifdef CGAL_NEF3_TIMER_PLUECKER
    timer_pluecker.stop();
     if(cgal_nef3_timer_on)
      std::cout << "Runtime_pluecker: " 
		<< timer_pluecker.time() << std::endl;
#endif
     //     SNC_io_parser<SNC_structure> O0(std::cerr,*this->sncp());
     //     O0.print();
    link_shalfedges_to_facet_cycles();

    std::map<int, int> hash;
    CGAL::Unique_hash_map<SHalfedge_handle, bool> done(false);

    SHalfedge_iterator sei;
    CGAL_forall_shalfedges(sei, *this->sncp()) {
      hash[sei->get_index()] = sei->get_index();
    }
    
    CGAL_forall_shalfedges(sei, *this->sncp()) {
      if(done[sei])
	continue;
      SHalfedge_around_facet_circulator circ(sei), end(circ);
      int index = circ->get_index();
      ++circ;
      CGAL_For_all(circ, end)
	if(circ->get_index() < index)
	  index = circ->get_index();      
      index = hash[index];
      CGAL_For_all(circ, end) {
	hash[circ->get_index()] = index;
	circ->set_index(index);
	done[circ] = true;
      }
    }
    
    SHalfloop_iterator sli;
    CGAL_forall_shalfloops(sli, *this->sncp())
      sli->set_index(hash[sli->get_index()]);

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
    Base::clear_external_structure();    
  }

  template<typename Association>
  void build_after_binary_operation(Association& A) {
    
    SHalfedge_iterator sei;
    CGAL_forall_shalfedges(sei, *this->sncp()) {
      CGAL_NEF_TRACEN("hash sedge " << sei->get_index() 
		      << "->" << A.get_hash(sei->get_index()));
      sei->set_index(A.get_hash(sei->get_index()));
    }

    SHalfloop_iterator sli;
    CGAL_forall_shalfloops(sli, *this->sncp()) {
      CGAL_NEF_TRACEN("hash sloop " << sli->get_index() 
		      << "->" << A.get_hash(sli->get_index()));
      sli->set_index(A.get_hash(sli->get_index()));
    }

    //    CGAL_NEF_SETDTHREAD(43);
    /*
    {    CGAL::SNC_io_parser<SNC_structure> O
      (std::cerr, *this->sncp(), false, true);
      O.print();}
    */
    pair_up_halfedges();
    /*
    {      CGAL::SNC_io_parser<SNC_structure> O
	(std::cerr, *this->sncp(), false, true);
      O.print();}
    */
    link_shalfedges_to_facet_cycles();

    SNC_simplify simp(*this->sncp());
    simp.vertex_simplificationI();

    //    std::map<int, int> hash;
    CGAL::Unique_hash_map<SHalfedge_handle, bool> 
      done(false);

    /*
    SHalfedge_iterator sei;
    CGAL_forall_shalfedges(sei, *this->sncp()) {
      hash[sei->get_forward_index()] = sei->get_forward_index();
      hash[sei->get_backward_index()] = sei->get_backward_index();
    }
    */

    categorize_facet_cycles_and_create_facets();
    create_volumes();
  }
}; 

} //namespace CGAL
#endif //CGAL_SNC_EXTERNAL_STRUCTURE_H
