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
// Author(s)     : Michael Seel    <seel@mpi-sb.mpg.de>
//                 Miguel Granados <granados@mpi-sb.mpg.de>
//                 Susan Hert      <hert@mpi-sb.mpg.de>
//                 Lutz Kettner    <kettner@mpi-sb.mpg.de>

#ifndef CGAL_SNC_IO_PARSER_H
#define CGAL_SNC_IO_PARSER_H

#include <CGAL/basic.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/Nef_S2/SM_decorator.h>
#include <CGAL/Nef_3/SNC_structure.h>
#include <CGAL/Nef_3/SNC_decorator.h>
#include <CGAL/Nef_2/Object_index.h>
#include <CGAL/Nef_3/Normalizing.h>
#include <vector>

#undef _DEBUG
#define _DEBUG 293
#include <CGAL/Nef_3/debug.h>

CGAL_BEGIN_NAMESPACE

template<typename T>
class moreLeft : public T {
  
  typedef typename T::SM_decorator      SM_decorator;
  typedef typename T::SHalfedge_handle  SHalfedge_handle;
  typedef typename T::Vector_3          Vector_3;
  typedef typename T::FT                FT;
  typedef typename T::RT                RT;
 
 public:
  moreLeft(T D) : T(D) {}
  
  bool operator()(SHalfedge_handle se1, SHalfedge_handle se2) {

    CGAL_assertion(se1 != SHalfedge_handle());
    if(se2 == SHalfedge_handle())
      return true;

    SM_decorator SM(&*vertex(se1));
    Vector_3 vec1 = SM.circle(se1).orthogonal_vector();
    Vector_3 vec2 = SM.circle(se2).orthogonal_vector();

    if(vec1 != vec2) {
      if(vec1.x() == RT(0) && vec2.x() == RT(0)) {
	if(vec1.y() != vec2.y())
	  return vec1.y() < vec2.y();
	return vec1.z() < vec2.z();
      }
      
      Vector_3 minus(-1,0,0);
      FT sk1(minus*vec1),  sk2(minus*vec2);
      if((sk1 >= FT(0) && sk2 <= FT(0)) ||
	 (sk1 <= FT(0) && sk2 >= FT(0)))
	return (sk1 > FT(0) || sk2 < FT(0));
      
      FT len1 = vec1.x()*vec1.x()+vec1.y()*vec1.y()+vec1.z()*vec1.z();
      FT len2 = vec2.x()*vec2.x()+vec2.y()*vec2.y()+vec2.z()*vec2.z();
      FT diff = len1*sk2*sk2 - len2*sk1*sk1;
      
      return ((sk1>FT(0) && diff<FT(0)) || (sk1<FT(0) && diff>FT(0)));
    }
    
    return lexicographically_xyz_smaller(SM.point(ssource(se1)), 
                                         SM.point(ssource(se2)));
  }
};

template <typename T>
class sort_vertices : public SNC_decorator<T> {
  
  typedef T SNC_structure;
  typedef SNC_decorator<T>          Base;
  typedef typename T::Vertex_handle Vertex_handle;
  typedef typename T::Point_3       Point_3;
  
 public:
  sort_vertices(T& D) : Base(D) {}
  
  bool operator() (Vertex_handle v1, Vertex_handle v2) const {
    return lexicographically_xyz_smaller(point(v1), point(v2));
  }
};
  
template <typename T>
class sort_edges : public SNC_decorator<T> {
  
  typedef T SNC_structure;  
  typedef SNC_decorator<T>            Base;
  typedef typename T::Halfedge_handle Halfedge_handle;
  
 public:
  sort_edges(T& D) : Base(D) {}
  
  bool operator() (Halfedge_handle e1, Halfedge_handle e2) const {
    sort_vertices<T> SORT(*sncp());
    if(source(e1) != source(e2))
      return SORT(source(e1),source(e2));
    return SORT(source(twin(e1)), source(twin(e2)));
  }
};

template <typename T>
class sort_facets : public SNC_decorator<T> {
  
  typedef T SNC_structure;  
  typedef SNC_decorator<T>             Base;
  typedef typename T::Halffacet_handle Halffacet_handle;
  typedef typename T::SHalfedge_handle SHalfedge_handle;
  typedef typename T::Vector_3         Vector_3;
  typedef typename T::Plane_3          Plane_3;

 public:
  sort_facets(T& D) : Base(D) {}
  
  bool operator() (Halffacet_handle f1, Halffacet_handle f2) const {
    
    Plane_3 p1(plane(f1));
    Plane_3 p2(plane(f2));
    
    if(p1.d() != p2.d())
      return p1.d() < p2.d();
    else if(p1.a() != p2.a())
      return p1.a() < p2.a();
    else if(p1.b() != p2.b())
      return p1.b() < p2.b();    
    else if(p1.c() != p2.c())
      return p1.c() < p2.c();    

    SHalfedge_handle se1 = SHalfedge_handle(f1->facet_cycles_begin());
    SHalfedge_handle se2 = SHalfedge_handle(f2->facet_cycles_begin());
    
    sort_vertices<T> SORT(*sncp());
    if(vertex(se1) != vertex(se2))
      return SORT(vertex(se1), vertex(se2));

    se1 = next(se1);
    se2 = next(se2);
    CGAL_assertion(vertex(se1) != vertex(se2));
    return SORT(vertex(se1), vertex(se2));

 
  }
};

template <typename T>
class sort_sedges : public SNC_decorator<T> {
  
  typedef T SNC_structure;  
  typedef SNC_decorator<T>             Base;
  typedef SM_decorator<T>          SM_decorator;
  typedef typename T::Vertex_handle    Vertex_handle;
  typedef typename T::SHalfedge_handle SHalfedge_handle;
  typedef typename T::Sphere_circle    Sphere_circle;

 public:
  sort_sedges(T& D) : Base(D) {}
  
  bool operator() (SHalfedge_handle se1, SHalfedge_handle se2) const {
    sort_vertices<T> SORT(*sncp());
    if(source(se1) != source(se2))
      return SORT(source(se1),source(se2));
    SM_decorator SD;
    if(se1 == twin(se2))
      if(ssource(se1) == ssource(se2)) {
	Sphere_circle vec1 = SD.circle(se1);
	Sphere_circle vec2 = SD.circle(se2);
	if(vec1.a() != vec2.a())
	  return vec1.a() < vec2.a();
	else if(vec1.b() != vec2.b())
	  return vec1.b() < vec2.b();
	return vec1.c() < vec2.c();     	
      }
      else
	return SORT(vertex(twin(ssource(se1))), vertex(twin(ssource(se2))));
    if(SORT(vertex(twin(ssource(twin(se1)))), 
	    vertex(twin(ssource(se1)))))
      se1 = twin(se1);
    if(SORT(vertex(twin(ssource(twin(se2)))), 
	    vertex(twin(ssource(se2)))))
      se2 = twin(se2);      
    if(ssource(se1) != ssource(se2))
      return SORT(vertex(twin(ssource(se1))), vertex(twin(ssource(se2))));
    if(target(se1) != target(se2))
      return SORT(target(se1), target(se2));

    CGAL_assertion(SD.circle(se1) != SD.circle(se2));
    Sphere_circle vec1 = SD.circle(se1);
    Sphere_circle vec2 = SD.circle(se2);
 
    if(vec1.a() != vec2.a())
      return vec1.a() < vec2.a();
    else if(vec1.b() != vec2.b())
      return vec1.b() < vec2.b();
    return vec1.c() < vec2.c();     
  }
};


template <typename T>
class sort_sloops : public SNC_decorator<T> {
  
  typedef T SNC_structure;  
  typedef SNC_decorator<T>             Base;
  typedef typename T::SHalfloop_handle SHalfloop_handle;
  
 public:
  sort_sloops(T& D) : Base(D) {}
  
  bool operator() (SHalfloop_handle sl1, SHalfloop_handle sl2) const {
    sort_vertices<T> SORTV(*sncp());
    sort_facets<T> SORTF(*sncp());
    if(vertex(sl1) != vertex(sl2))
      return SORTV(vertex(sl1),vertex(sl2));
    return SORTF(facet(sl1), facet(sl2));
  }
};

template <typename T>
class sort_sfaces : public SNC_decorator<T> {
  
  typedef T SNC_structure;
  typedef SNC_decorator<T>                  Base;
  typedef typename T::SM_decorator          SM_decorator;
  typedef typename T::Point_3               Point_3;
  typedef typename T::Vector_3              Vector_3;
  typedef typename T::SVertex_handle        SVertex_handle;
  typedef typename T::SHalfedge_handle      SHalfedge_handle;
  typedef typename T::SHalfloop_handle      SHalfloop_handle;
  typedef typename T::SFace_handle          SFace_handle;
  typedef typename T::SFace_cycle_iterator  SFace_cycle_iterator;
  typedef typename T::SHalfedge_around_sface_circulator
                      SHalfedge_around_sface_circulator;
  
 public:
  sort_sfaces(T& D) : Base(D) {}

  bool operator() (SFace_handle sf1, SFace_handle sf2) const {

    TRACEN("");
    TRACEN("sort sfaces");

    sort_vertices<T> SORT(*sncp());
    
    if(vertex(sf1) != vertex(sf2))
      return SORT(vertex(sf1), vertex(sf2));
    
    if(sf1 == sf2)
      return false;

    SM_decorator SD(&*vertex(sf1));
    moreLeft<Base> ml((Base) *this);
    Vector_3 plus(1,0,0);

    CGAL_assertion_code(SVertex_handle sv);
    SHalfedge_handle se;
    SFace_cycle_iterator fc;
  
    TRACEN("vertex " << point(vertex(sf1)));
    TRACEN("sface 1");

    SHalfedge_handle se1;
    SHalfloop_handle sl1;
    CGAL_forall_sface_cycles_of(fc,sf1) {
      if(assign(se,fc)) {
	SHalfedge_around_sface_circulator ec(se),ee(se);
	CGAL_For_all(ec,ee) {
	  TRACEN("   " << SD.point(ssource(ec)) << 
		 " | " << SD.circle(ec).orthogonal_vector());
	  if(ml(ec, se1))
	    se1 = ec;
	}
      }
      else if(!assign(sl1,fc))
	CGAL_assertion(assign(sv,fc));
    }

    TRACEN("sface 2");
  
    SHalfedge_handle se2;
    SHalfloop_handle sl2;
    CGAL_forall_sface_cycles_of(fc,sf2) {
      if(assign(se,fc)) {
	SHalfedge_around_sface_circulator ec(se),ee(se);
	CGAL_For_all(ec,ee) { 
	  TRACEN("   " << SD.point(ssource(ec)) << 
		 " | " << SD.circle(ec).orthogonal_vector());
	  if(ml(ec, se2))
	    se2 = ec;
	}
      }
      else if(!assign(sl2,fc))
	CGAL_assertion(assign(sv,fc));
    }
  
    if(se1 != SHalfedge_handle() && se2 == SHalfedge_handle())
      return true;
    if(se1 == SHalfedge_handle() && se2 != SHalfedge_handle())
      return false;  
  
    if(se1 == SHalfedge_handle() && se2 == SHalfedge_handle()) {
      Vector_3 vec1 = SD.circle(sl1).orthogonal_vector();
      Vector_3 vec2 = SD.circle(sl2).orthogonal_vector();
      if(vec1.x() != vec2.x())
	return vec1.x() < vec2.x();
      else if(vec1.y() != vec2.y())
	return vec1.y() < vec2.y();
      else if(vec1.z() != vec2.z())
	return vec1.z() < vec2.z();     
    }
    
    CGAL_assertion(se1 != SHalfedge_handle() && se2 != SHalfedge_handle());

    TRACEN("minimal sedge in sface 1:" << SD.point(ssource(se1)) << 
	   " | " << SD.circle(se1).orthogonal_vector());
    TRACEN("minimal sedge in sface 2:" << SD.point(ssource(se2)) << 
	   " | " << SD.circle(se2).orthogonal_vector());
    TRACE("result " << ml(se1,se2));
    return ml(se1, se2);
  }
 
};

template <typename T>
class sort_volumes : public SNC_decorator<T> {
  
  typedef T SNC_structure;
  typedef SNC_decorator<T>          Base;
  typedef typename T::Volume_handle Volume_handle;
  typedef typename T::SFace_handle  SFace_handle;
  
 public:
  sort_volumes(T& D) : Base(D) {}
  
  bool operator() (Volume_handle c1, Volume_handle c2) const {
    SFace_handle sf1 = SFace_handle(c1->shells_begin()); 
    SFace_handle sf2 = SFace_handle(c2->shells_begin()); 

    sort_sfaces<T> SORT(*sncp());
    return SORT(sf1, sf2);
  }
};

template <typename T>
class sort_facet_cycle_entries : public T {
  
  typedef typename T::SNC_structure     SNC_structure;
  typedef typename T::SM_decorator      SM_decorator;
  typedef typename T::Object_handle     Object_handle;
  typedef typename T::SHalfedge_handle  SHalfedge_handle;
  typedef typename T::SHalfloop_handle  SHalfloop_handle;
  typedef typename T::SFace_handle      SFace_handle;
  typedef typename T::Point_3           Point_3;
  typedef typename T::Vector_3          Vector_3;
  
 public:
  sort_facet_cycle_entries(T D) : T(D) {}
  
  bool operator() (Object_handle o1, Object_handle o2) const {
    
    SHalfedge_handle se1, se2;
    SHalfloop_handle sl1, sl2;

    if(!assign(se1,o1) && !assign(sl1,o1))
      CGAL_assertion_msg(0,"wrong handle");
    
    if(!assign(se2,o2) && !assign(sl2,o2))
      CGAL_assertion_msg(0,"wrong handle");    
    
    if(se1 != SHalfedge_handle() && se2 != SHalfedge_handle()) {
      sort_vertices<SNC_structure> SORT(*sncp());
      return SORT(source(se1), source(se2));
    }

    if(se1 != SHalfedge_handle())
      return true;
    if(se2 != SHalfedge_handle())
      return false;

    CGAL_assertion(sl1 != SHalfloop_handle() && 
			sl2 != SHalfloop_handle());

    SM_decorator SD(&*vertex(sl1));
    Vector_3 vec1(SD.circle(sl1).orthogonal_vector());
    Vector_3 vec2(SD.circle(sl2).orthogonal_vector());
    //    CGAL_assertion(vec1 == vec2.antipode());
    if(vec1.x() != vec2.x())
      return vec1.x() < vec2.x();
    else if(vec1.y() != vec2.y())
      return vec1.y() < vec2.y();
    else
      return vec1.z() < vec2.z();          
  }
};

template <typename T>
class sort_sface_cycle_entries : public T {
  
  typedef typename T::SNC_structure     SNC_structure;
  typedef typename T::SM_decorator      SM_decorator;
  typedef typename T::Object_handle Object_handle;
  typedef typename T::SVertex_handle  SVertex_handle;
  typedef typename T::SHalfedge_handle  SHalfedge_handle;
  typedef typename T::SHalfloop_handle  SHalfloop_handle;
  typedef typename T::SFace_handle  SFace_handle;
  typedef typename T::Point_3       Point_3;
  typedef typename T::Vector_3       Vector_3;
  
 public:
  sort_sface_cycle_entries(T D) : T(D) {}
  
  bool operator() (Object_handle o1, Object_handle o2) const {
    
    SVertex_handle sv1, sv2;
    SHalfedge_handle se1, se2;
    SHalfloop_handle sl1, sl2;

    if(!assign(se1,o1) && !assign(sl1,o1) && !assign(sv1,o1))
      CGAL_assertion_msg(0,"wrong handle");
    
    if(!assign(se2,o1) && !assign(sl2,o2) && !assign(sv2,o2))
      CGAL_assertion_msg(0,"wrong handle");    
    
    if(se1 != SHalfedge_handle() && se2 == SHalfedge_handle())
      return true;

    if(se1 == SHalfedge_handle() && se2 != SHalfedge_handle())
      return false;

    if(sl1 != SHalfloop_handle() && sv2 != SVertex_handle())
      return true;

    if(sl2 != SHalfloop_handle() && sv1 != SVertex_handle())
      return false;

    if(se1 != SHalfedge_handle() && se2 != SHalfedge_handle()) {
      sort_vertices<SNC_structure> SORT(*sncp());
      if(ssource(se1) != ssource(se2))
	return SORT(vertex(twin(ssource(se1))), vertex(twin(ssource(se2))));
      else
	return SORT(target(se1), target(se2));
    }

    if(sl1 != SHalfloop_handle() && sl2 != SHalfloop_handle()) {
      SM_decorator SD(&*vertex(sl1));
      Vector_3 vec1(SD.circle(sl1).orthogonal_vector());
      Vector_3 vec2(SD.circle(sl2).orthogonal_vector());
      //      CGAL_assertion(vec1 == vec2.antipode());
      if(vec1.x() != vec2.x())
	return vec1.x() < vec2.x();
      else if(vec1.y() != vec2.y())
	return vec1.y() < vec2.y();
      else if(vec1.z() != vec2.z())
	return vec1.z() < vec2.z();          
    }

    CGAL_assertion(sv1 != SVertex_handle() && sv2 != SVertex_handle());
    sort_vertices<SNC_structure> SORT(*sncp());
    return SORT(target(sv1), target(sv2));
  }
};

template <typename T>
class sort_shell_entries : public T {
  
  typedef typename T::Object_handle Object_handle;
  typedef typename T::Shell_entry_iterator  Shell_entry_iterator;
  typedef typename T::SFace_handle  SFace_handle;
  typedef typename T::Point_3       Point_3;
  
 public:
  sort_shell_entries(T D) : T(D) {}
  
  bool operator() (Object_handle o1, Object_handle o2) const {
    SFace_handle sf1, sf2;
    assign(sf1, o1);
    assign(sf2, o2);
    Point_3 p1(point(vertex(sf1))), p2(point(vertex(sf2)));
    if(p1.x() != p2.x())
      return p1.x() < p2.x();
    else if(p1.y() != p2.y())
      return p1.y() < p2.y();
    return p1.z() < p2.z();
  }
};

template<typename T>
struct find_minimal_sface_of_shell : public SNC_decorator<T> {
  
  typedef T                        SNC_structure;
  typedef SNC_decorator<T>          Base;
  typedef typename T::Vertex_handle  Vertex_handle;
  typedef typename T::Halfedge_handle  Halfedge_handle;
  typedef typename T::Halffacet_handle  Halffacet_handle;
  typedef typename T::SFace_handle  SFace_handle;
  typedef CGAL::Unique_hash_map<SFace_handle,bool> SFace_visited_hash;

  SFace_visited_hash& Done;
  SFace_handle sf_min;
  sort_sfaces<T> SORT;
  
  find_minimal_sface_of_shell(T& D, SFace_visited_hash& Vi) 
    : Base(D), Done(Vi), SORT(D) {}
  
  void visit(SFace_handle h) { 
    Done[h]=true;
    if(sf_min == SFace_handle())
      sf_min = h;
    else {
      if(SORT(h,sf_min))
	sf_min = h;
    }
  }
  
  void visit(Vertex_handle h) {}
  void visit(Halfedge_handle h) {}
  void visit(Halffacet_handle h) {}
  
  SFace_handle& minimal_sface() { return sf_min; }
};

template <typename SNC_structure_>
class SNC_io_parser : public SNC_decorator<SNC_structure_>
{ typedef SNC_structure_ SNC_structure;
  typedef CGAL::SNC_io_parser<SNC_structure_> Self;
  typedef CGAL::SNC_decorator<SNC_structure_> Base;
  typedef typename SNC_structure::Sphere_map  Sphere_map;
  typedef CGAL::SM_decorator<Sphere_map>      SM_decorator;
  typedef typename SNC_structure::Infi_box    Infi_box;
  typedef typename Infi_box::Standard_kernel  Standard_kernel;
public:
  typedef typename SNC_structure::Vertex_iterator Vertex_iterator; 
  typedef typename SNC_structure::Vertex_handle Vertex_handle;
  typedef typename SNC_structure::Halfedge_iterator Halfedge_iterator; 
  typedef typename SNC_structure::Halfedge_handle Halfedge_handle;
  typedef typename SNC_structure::Halffacet_iterator Halffacet_iterator; 
  typedef typename SNC_structure::Halffacet_handle Halffacet_handle;
  typedef typename SNC_structure::Volume_iterator Volume_iterator; 
  typedef typename SNC_structure::Volume_handle Volume_handle;
  typedef typename SNC_structure::SVertex_iterator SVertex_iterator; 
  typedef typename SNC_structure::SVertex_handle SVertex_handle;
  typedef typename SNC_structure::SHalfedge_iterator SHalfedge_iterator; 
  typedef typename SNC_structure::SHalfedge_handle SHalfedge_handle;
  typedef typename SNC_structure::SFace_iterator SFace_iterator; 
  typedef typename SNC_structure::SFace_handle SFace_handle;
  typedef typename SNC_structure::SHalfloop_iterator SHalfloop_iterator; 
  typedef typename SNC_structure::SHalfloop_handle SHalfloop_handle;
  typedef typename SNC_structure::Object_iterator Object_iterator; 
  typedef typename SNC_structure::Object_handle Object_handle;
  typedef typename SNC_structure::SFace_cycle_iterator SFace_cycle_iterator;
  typedef typename SNC_structure::Halffacet_cycle_iterator Halffacet_cycle_iterator;
  typedef typename SNC_structure::Shell_entry_iterator Shell_entry_iterator;
  typedef typename SNC_structure::SHalfedge_around_svertex_circulator 
                                  SHalfedge_around_svertex_circulator;
  typedef typename SNC_structure::SHalfedge_around_sface_circulator 
                                  SHalfedge_around_sface_circulator;
  typedef typename SNC_structure::Point_3 Point_3;
  typedef typename SNC_structure::Plane_3 Plane_3;
  typedef typename SNC_structure::Vector_3 Vector_3;
  typedef typename SNC_structure::Sphere_point Sphere_point;
  typedef typename SNC_structure::Sphere_segment Sphere_segment;
  typedef typename SNC_structure::Sphere_circle Sphere_circle;
  typedef typename SNC_structure::Mark Mark;
  typedef typename SNC_structure::Kernel Kernel;
  typedef typename Kernel::RT RT;
  typedef void* GenPtr;

 private:
  std::istream& in; std::ostream& out;
  bool verbose;
  bool reduce;
  bool sorted;
  bool addInfiBox;

  CGAL::Object_index<Vertex_iterator> VI;  
  CGAL::Object_index<Halfedge_iterator> EI;
  CGAL::Object_index<Halffacet_iterator>    FI;
  CGAL::Object_index<Volume_iterator>   CI;
  CGAL::Object_index<SHalfedge_iterator> SEI;
  CGAL::Object_index<SHalfloop_iterator>   SLI;
  CGAL::Object_index<SFace_iterator>     SFI;
  std::list<Vertex_iterator> VL;
  std::list<Halfedge_iterator> EL;
  std::list<Halffacet_iterator> FL;
  std::list<Volume_iterator> CL;
  std::list<SHalfedge_iterator> SEL;
  std::list<SHalfloop_iterator> SLL;
  std::list<SFace_iterator> SFL;
  std::vector<Vertex_iterator>   Vertex_of;
  std::vector<Halfedge_iterator> Edge_of;
  std::vector<Halffacet_iterator>    Halffacet_of;
  std::vector<Volume_iterator>   Volume_of;
  std::vector<SVertex_iterator>  SVertex_of; 
  std::vector<SHalfedge_iterator> SEdge_of;
  std::vector<SHalfloop_iterator> SLoop_of;
  std::vector<SFace_iterator>     SFace_of;
  long i,vn,en,fn,cn,sen,sln,sfn;

public:
  SNC_io_parser(std::istream& is, SNC_structure& W);
  SNC_io_parser(std::ostream& os, SNC_structure& W, bool sort);

  std::string index(Vertex_iterator v) const
  { return VI(v,verbose); } 
  std::string index(Halfedge_iterator e) const 
  { return EI(e,verbose); }
  std::string index(Halffacet_iterator f) const 
  { return FI(f,verbose); }
  std::string index(Volume_iterator c) const 
  { return CI(c,verbose); }
  std::string index(SHalfedge_iterator e) const 
  { return SEI(e,verbose); }
  std::string index(SHalfloop_iterator l) const 
  { return SLI(l,verbose); }
  std::string index(SFace_iterator f) const 
  { return SFI(f,verbose); }
  std::string index(Object_iterator o) const
  { if( o == 0 )
      return string("undef");
    Vertex_iterator v;
    Halfedge_iterator e;
    Halffacet_iterator f;
    Volume_iterator c;
    SHalfedge_iterator se;
    SHalfloop_iterator sl;
    SFace_iterator sf;
    if( assign( v, *o))
      return index(v);
    else if( assign( e, *o))
      return index(e);
    else if( assign( f, *o))
      return index(f);
    else if( assign( c, *o))
      return index(c);
    else if( assign( se, *o))
      return index(se);
    else if( assign( sl, *o))
      return index(sl);
    else if( assign( sf, *o))
      return index(sf);
    return string("unknown object");
  }

  bool check_sep(char* sep) const;
  bool test_string(std::string s) const;
  void print_vertex(Vertex_handle) const;
  void print_edge(Halfedge_handle) const;
  void print_facet(Halffacet_handle) const;
  void print_volume(Volume_handle) const;
  void print_sedge(SHalfedge_handle) const;
  void print_sloop(SHalfloop_handle) const;
  void print_sface(SFace_handle) const;
  void print() const;
  void print_local_graph(Vertex_handle) const;

  template <typename NT> bool read_vertex(Vertex_handle);
  template <typename NT> bool read_edge(Halfedge_handle);
  template <typename NT> bool read_facet(Halffacet_handle);
  template <typename NT> bool read_volume(Volume_handle);
  template <typename NT> bool read_svertex(SVertex_handle);
  template <typename NT> bool read_sedge(SHalfedge_handle);
  template <typename NT> bool read_sloop(SHalfloop_handle);
  template <typename NT> bool read_sface(SFace_handle);
  void add_infi_box();
  void read();

  static void dump(SNC_structure& W, std::ostream& os = std::cerr, bool sort = false)
  { Self O(os,W, sort); O.print(); }

  template <typename Iter, typename Index>
    void output_sorted_indexes(Iter begin, Iter end, Index i) const {
    int low = i[begin];
    int high = low;
    for(Iter it=begin; it != end; it++) {
      if(i[it] < low) low = i[it];
      if(i[it] > high) high = i[it];
    }
    out << low << " " << high << ", ";
  }

};

template <typename EW>
SNC_io_parser<EW>::SNC_io_parser(std::istream& is, SNC_structure& W) : 
  Base(W), in(is), out(std::cout) { 
  W.clear();
  CGAL_assertion(W.is_empty());
  verbose = false; 
}


template <typename EW>
SNC_io_parser<EW>::SNC_io_parser(std::ostream& os, SNC_structure& W, 
                                                   bool sort = false) : 
  Base(W), in(std::cin), out(os),
  FI(W.halffacets_begin(),W.halffacets_end(),'F'),
  CI(W.volumes_begin(),W.volumes_end(),'C'),
  SEI(W.shalfedges_begin(),W.shalfedges_end(),'e'),
  SLI(W.shalfloops_begin(),W.shalfloops_end(),'l'),
  SFI(W.sfaces_begin(),W.sfaces_end(),'f'),
  vn(W.number_of_vertices()), 
  en(W.number_of_halfedges()), 
  fn(W.number_of_halffacets()),
  cn(W.number_of_volumes()),
  sen(W.number_of_shalfedges()),
  sln(W.number_of_shalfloops()),
  sfn(W.number_of_sfaces())
{ 

  verbose = (out.iword(CGAL::IO::mode) != CGAL::IO::ASCII &&
             out.iword(CGAL::IO::mode) != CGAL::IO::BINARY);  
  reduce = false; 
  sorted = sort;

  if(sorted && Infi_box::extended_kernel()) {
    int i = 0;
    Halffacet_handle hf;
    CGAL_forall_facets(hf, *sncp()) {
      if(!Infi_box::is_standard(plane(hf)))
	i ++;
    }
    //    TRACEN("There are " << i  << " extended facets");
    if(i < 7)
      reduce = true;
  }
  
  Vertex_iterator vi; 
  CGAL_forall_vertices(vi, *sncp()) {
    VL.push_back(vi);
    if(sorted)
      point(vi) = normalized(point(vi));
  }
  if(sorted) {
    VL.sort(sort_vertices<SNC_structure>(*sncp()));
  }  
  if(reduce)
    for(int k=0; k<4; k++){
      VL.pop_front(); VL.pop_back();
    }
  int i = 0;
  typename std::list<Vertex_iterator>::iterator vl;
  for(vl = VL.begin(); vl != VL.end(); vl++)
    VI[*vl] = i++;
  
  SM_decorator SD;
  Halfedge_iterator ei; 
  CGAL_forall_halfedges(ei, *sncp()) {
    EL.push_back(ei);
    if(sorted) {
      //      std::cerr << point(ei) << " | " << normalized(point(ei)) << " |";
      SD.point(ei) = normalized(SD.point(ei));
      //      std::cerr << point(ei) << std::endl;
      SHalfedge_handle new_outedge = ei->out_sedge();
      SHalfedge_around_svertex_circulator cb(new_outedge), ce(cb);
      CGAL_For_all(cb,ce) {
	if(source(cb) == target(cb) && source(cb) == target(new_outedge)) {
	  Plane_3 c(SD.circle(cb));
	  if(c.a() != 0) {
	    if(c.a() < 0)
	      new_outedge = cb;
	  }
	  else if(c.b() != 0) {
	    if(c.b() < 0)
	      new_outedge = cb;
	  }
	  else if(c.c() != 0) {
	    if(c.c() < 0)
	      new_outedge = cb;
	  }
	} else
	  if(lexicographically_xyz_smaller(point(target(cb)), point(target(new_outedge))))													
	    new_outedge = cb;
      }
      ei->out_sedge() = new_outedge;
    }
  }
  if(sorted) EL.sort(sort_edges<SNC_structure>(*sncp()));
  if(reduce)
    for(int k=0; k<12; k++){
      EL.pop_front(); EL.pop_back();
    }
  i = 0;
  typename std::list<Halfedge_iterator>::iterator el;
  for(el = EL.begin(); el != EL.end(); el++)
    EI[*el] = i++;

  Halffacet_iterator fi; 
  CGAL_forall_halffacets(fi, *sncp()){
    if(sorted) {
      plane(fi) = normalized(plane(fi));
      fi->boundary_entry_objects_.sort(sort_facet_cycle_entries<Base>((Base) *this));
    }
    FL.push_back(fi);
  }
  if(sorted) FL.sort(sort_facets<SNC_structure>(*sncp()));
  if(reduce) {
    for(int k=0; k<6; k++){
      FL.pop_front();
      FL.pop_back();
    }
  }
  i = 0;
  typename std::list<Halffacet_iterator>::iterator fl;
  for(fl = FL.begin(); fl != FL.end(); fl++)
    FI[*fl] = i++;

  Volume_iterator ci; 
  CGAL::Unique_hash_map<SFace_handle,bool> Done(false);
  find_minimal_sface_of_shell<SNC_structure> findMinSF(*sncp(),Done);
  CGAL_forall_volumes(ci, *sncp()) {
    if(sorted) {
      Shell_entry_iterator it;
      CGAL_forall_shells_of(it,ci) {
	findMinSF.minimal_sface() = SFace_handle(it);
	visit_shell_objects(SFace_handle(it),findMinSF);
	*it = findMinSF.minimal_sface();
      }
      ci->shell_entry_objects().sort(sort_shell_entries<Base>((Base)*this));
    } 
    CL.push_back(ci);
  }
 
  if(sorted) CL.sort(sort_volumes<SNC_structure>(*sncp()));
  if(reduce)
    CL.pop_front();
  i = 0;
  typename std::list<Volume_iterator>::iterator cl;
  for(cl = CL.begin(); cl != CL.end(); cl++)
    CI[*cl] = i++;

  SHalfedge_iterator sei; 
  CGAL_forall_shalfedges(sei, *sncp()) {
    SEL.push_back(sei);
    if(sorted)
      SD.circle(sei) = normalized(SD.circle(sei));
  }
  if(sorted) SEL.sort(sort_sedges<SNC_structure>(*sncp()));
  if(reduce)
    for(int k=0; k<24; k++){
      SEL.pop_front(); SEL.pop_back();
    }
  i = 0;
  typename std::list<SHalfedge_iterator>::iterator sel;
  for(sel = SEL.begin(); sel != SEL.end(); sel++)
    SEI[*sel] = i++;

  SHalfloop_iterator sli; 
  CGAL_forall_shalfloops(sli, *sncp()) {
    SLL.push_back(sli);
    if(sorted)
      SD.circle(sli) = normalized(SD.circle(sli));
  }
  if(sorted) SLL.sort(sort_sloops<SNC_structure>(*sncp()));
  i = 0;
  typename std::list<SHalfloop_iterator>::iterator sll;
  for(sll = SLL.begin(); sll != SLL.end(); sll++)
    SLI[*sll] = i++;

  SFace_iterator sfi; 
  CGAL_forall_sfaces(sfi, *sncp()) {
    if(sorted) {
      SFace_cycle_iterator fc;
      CGAL_forall_sface_cycles_of(fc, sfi) {
	SHalfedge_handle se;
	if(assign(se,fc)) {
	  SHalfedge_around_sface_circulator cb(se), ce(cb);
	  CGAL_For_all(cb,ce) {
	    if(ssource(cb) != ssource(se)) {
	      if(lexicographically_xyz_smaller(point(vertex(twin(ssource(cb)))), 
					       point(vertex(twin(ssource(se))))))
		se = cb;
	    }
	    else 
	      if(lexicographically_xyz_smaller(point(target(cb)), point(target(se))))
		se = cb;
	  }
	  *fc = se;
	}
      }
      sfi->boundary_entry_objects_.sort(sort_sface_cycle_entries<Base>((Base) *this));
    }
    SFL.push_back(sfi);
  }
  if(sorted) SFL.sort(sort_sfaces<SNC_structure>(*sncp()));
  if(reduce)
    for(int k=0; k<8; k++){
      SFL.pop_front(); SFL.pop_back();
    }
  i = 0;
  typename std::list<SFace_iterator>::iterator sfl;
  for(sfl = SFL.begin(); sfl != SFL.end(); sfl++)
    SFI[*sfl] = i++;

  VI[W.vertices_end()]=-2;
  EI[W.halfedges_end()]=-2;
  FI[W.halffacets_end()]=-2;
  CI[W.volumes_end()]=-2;
  SEI[W.shalfedges_end()]=-2;
  SLI[W.shalfloops_end()]=-2;
  SFI[W.sfaces_end()]=-2;
}

template <typename EW>
bool SNC_io_parser<EW>::check_sep(char* sep) const
{
  char c; 
  do in.get(c); while (isspace(c));
  while (*sep != '\0') { 
    if (*sep != c) {
      in.putback(c);
      return false;
    }
    ++sep; in.get(c);
  }
  in.putback(c);
  return true;  
}

template <typename EW>
bool SNC_io_parser<EW>::test_string(std::string s) const {
  std::string s2;
  in >> s2;
  return (s==s2);
}

template <typename EW>
void SNC_io_parser<EW>::print() const
{ 
  out << "Selective Nef Complex" << std::endl;
  out << "vertices   " << VL.size() << std::endl;
  out << "halfedges  " << EL.size() << std::endl;
  out << "facets     " << FL.size() << std::endl;
  out << "volumes    " << CL.size() << std::endl;
  out << "shalfedges " << SEL.size() << std::endl;
  out << "shalfloops " << SLL.size() << std::endl;
  out << "sfaces     " << SFL.size() << std::endl;

  if (verbose) 
    out << "/* Vertex: index { svs sve ses see sfs sfe sl,"
        << " mark, point } */\n";
  typename std::list<Vertex_iterator>::const_iterator v;
  for(v=VL.begin();v!=VL.end();v++)
    print_vertex(*v);

  if (verbose) 
  out << "/* Edge: index { twin, source, isolated incident_object,"
      << " mark } */\n";
  typename std::list<Halfedge_iterator>::const_iterator e;
  for(e=EL.begin();e!=EL.end();e++)
    print_edge(*e);

  if (verbose) 
  out << "/* Facet: index { twin, fclist, ivlist, volume | plane } mark */\n";
  typename std::list<Halffacet_iterator>::const_iterator f;
  for(f=FL.begin();f!=FL.end();f++)
    print_facet(*f);

  if (verbose) 
  out << "/* Volume: index { shlist } mark  */\n";
  typename std::list<Volume_iterator>::const_iterator c;
  for(c=CL.begin();c!=CL.end();c++)
    print_volume(*c);

  if (verbose) 
  out << "/* SEdge: index { twin, sprev, snext, source, sface,"
      << " prev, next, facet } */\n";
  typename std::list<SHalfedge_iterator>::const_iterator se;
  for(se=SEL.begin();se!=SEL.end();se++)
    print_sedge(*se);

  if (verbose) 
  out << "/* SLoop: index { twin, sface, facet } */" << std::endl;
  typename std::list<SHalfloop_iterator>::const_iterator sl;
  for(sl=SLL.begin();sl!=SLL.end();sl++)
    print_sloop(*sl);

  if (verbose) 
  out << "/* SFace: index { fclist, ivlist, sloop, volume } */" << std::endl;
  typename std::list<SFace_iterator>::const_iterator sf;
  for(sf=SFL.begin();sf!=SFL.end();sf++)
    print_sface(*sf);

  out << "/* end Selective Nef complex */" << std::endl;
}

template <typename EW>
void SNC_io_parser<EW>::read()
{ 
  if ( !check_sep("Selective Nef Complex") )  
    CGAL_assertion_msg(0,"SNC_io_parser::read: no SNC header.");
  if ( !(check_sep("vertices") && (in >> vn)) ) 
    CGAL_assertion_msg(0,"SNC_io_parser::read: wrong vertex line.");
  if ( !(check_sep("halfedges") && (in >> en) && (en%2==0)) )
    CGAL_assertion_msg(0,"SNC_io_parser::read: wrong edge line.");
  if ( !(check_sep("facets") && (in >> fn) && (fn%2==0)) )
    CGAL_assertion_msg(0,"SNC_io_parser::read: wrong facet line.");
  if ( !(check_sep("volumes") && (in >> cn)) )
    CGAL_assertion_msg(0,"SNC_io_parser::read: wrong volume line.");
  if ( !(check_sep("shalfedges") && (in >> sen)) )
    CGAL_assertion_msg(0,"SNC_io_parser::read: wrong sedge line.");
  if ( !(check_sep("shalfloops") && (in >> sln)) )
    CGAL_assertion_msg(0,"SNC_io_parser::read: wrong sloop line.");
  if ( !(check_sep("sfaces") && (in >> sfn)) )
    CGAL_assertion_msg(0,"SNC_io_parser::read: wrong sface line.");

  char c;
  int spaces = 0;
  streampos pos(in.tellg());
  in.get(c);
  while(in.get(c) && c != '\n')
    if(c == ' ')
      ++spaces;
  in.seekg(pos,ios::beg); 

  CGAL_assertion_msg(c=='\n' && (spaces==15 || spaces==26),"wrong file format");
  CGAL_assertion_msg(!(Infi_box::standard_kernel() && spaces==26), 
		     "you cannot load an inifinite polyhedron without an extended kernel.");
  addInfiBox = (spaces==15 && Infi_box::extended_kernel());

  if(addInfiBox) {
    Vertex_of.reserve(vn+8);
    Edge_of.reserve(en+24);
    Halffacet_of.reserve(fn+12);
    Volume_of.reserve(cn+1);
    SEdge_of.reserve(sen+48);
    SLoop_of.reserve(sln);
    SFace_of.reserve(sfn+16);
    for(i=0; i<vn+8; ++i)   Vertex_of[i] = sncp()->new_vertex_only();
    for(i=0; i<en+24; ++i)  Edge_of[i] = sncp()->new_halfedge_only();
    for(i=0; i<fn+12; ++i)  Halffacet_of[i] = sncp()->new_halffacet_only();
    for(i=0; i<cn+1; ++i)   Volume_of[i] = sncp()->new_volume_only();
    for(i=0; i<sen+48; ++i) SEdge_of[i] = sncp()->new_shalfedge_only();
    for(i=0; i<sln; ++i)    SLoop_of[i] = sncp()->new_shalfloop_only();
    for(i=0; i<sfn+16; ++i) SFace_of[i] = sncp()->new_sface_only();
  }
  else {
    Vertex_of.reserve(vn);
    Edge_of.reserve(en);
    Halffacet_of.reserve(fn);
    Volume_of.reserve(cn);
    SEdge_of.reserve(sen);
    SLoop_of.reserve(sln);
    SFace_of.reserve(sfn);
    for(i=0; i<vn; ++i)  Vertex_of[i] = sncp()->new_vertex_only();
    for(i=0; i<en; ++i)  Edge_of[i] = sncp()->new_halfedge_only();
    for(i=0; i<fn; ++i)  Halffacet_of[i] = sncp()->new_halffacet_only();
    for(i=0; i<cn; ++i)  Volume_of[i] = sncp()->new_volume_only();
    for(i=0; i<sen; ++i) SEdge_of[i] = sncp()->new_shalfedge_only();
    for(i=0; i<sln; ++i) SLoop_of[i] = sncp()->new_shalfloop_only();
    for(i=0; i<sfn; ++i) SFace_of[i] = sncp()->new_sface_only();
  }



  if(addInfiBox) {
    for(i=0; i<vn; ++i) {
      if (!read_vertex<typename Standard_kernel::RT>(Vertex_of[i]))
	CGAL_assertion_msg(0,"SNC_io_parser::read: error in node line");
    }
    for(i=0; i<en; ++i) {
      if (!read_edge<typename Standard_kernel::RT>(Edge_of[i]))
	CGAL_assertion_msg(0,"SNC_io_parser::read: error in edge line");
    }
    for(i=0; i<fn; ++i) {
      if (!read_facet<typename Standard_kernel::RT>(Halffacet_of[i]))
	CGAL_assertion_msg(0,"SNC_io_parser::read: error in facet line");
    }
    for(i=0; i<cn; ++i) {
      if (!read_volume<typename Standard_kernel::RT>(Volume_of[i+1]))
	CGAL_assertion_msg(0,"SNC_io_parser::read: error in volume line");
    }
    for(i=0; i<sen; ++i) {
      if (!read_sedge<typename Standard_kernel::RT>(SEdge_of[i]))
	CGAL_assertion_msg(0,"SNC_io_parser::read: error in sedge line");
    }
    for(i=0; i<sln; ++i) {
      if (!read_sloop<typename Standard_kernel::RT>(SLoop_of[i]))
	CGAL_assertion_msg(0,"SNC_io_parser::read: error in sloop line");
    }
    for(i=0; i<sfn; ++i) {
      if (!read_sface<typename Standard_kernel::RT>(SFace_of[i]))
	CGAL_assertion_msg(0,"SNC_io_parser::read: error in sface line");
    }
    add_infi_box();
  }
  else {
    for(i=0; i<vn; ++i) {
      if (!read_vertex<RT>(Vertex_of[i]))
	CGAL_assertion_msg(0,"SNC_io_parser::read: error in node line");
    }
    for(i=0; i<en; ++i) {
      if (!read_edge<RT>(Edge_of[i]))
	CGAL_assertion_msg(0,"SNC_io_parser::read: error in edge line");
    }
    for(i=0; i<fn; ++i) {
      if (!read_facet<RT>(Halffacet_of[i]))
	CGAL_assertion_msg(0,"SNC_io_parser::read: error in facet line");
    }
    for(i=0; i<cn; ++i) {
      if (!read_volume<RT>(Volume_of[i]))
	CGAL_assertion_msg(0,"SNC_io_parser::read: error in volume line");
    }
    for(i=0; i<sen; ++i) {
      if (!read_sedge<RT>(SEdge_of[i]))
	CGAL_assertion_msg(0,"SNC_io_parser::read: error in sedge line");
    }
    for(i=0; i<sln; ++i) {
      if (!read_sloop<RT>(SLoop_of[i]))
	CGAL_assertion_msg(0,"SNC_io_parser::read: error in sloop line");
    }
    for(i=0; i<sfn; ++i) {
      if (!read_sface<RT>(SFace_of[i]))
	CGAL_assertion_msg(0,"SNC_io_parser::read: error in sface line");
    }
  }  
}


template <typename EW>
void SNC_io_parser<EW>::print_vertex(Vertex_handle v) const
{ // syntax: index { svs sve, ses see, sfs sfe, sl | point } mark
  SM_decorator SD(&*v);
  out << index(v) << " { ";
  if(sorted) {
    
    output_sorted_indexes(v->svertices_begin(), 
			  v->svertices_end(), EI);
    output_sorted_indexes(v->shalfedges_begin(), 
			  v->shalfedges_end(), SEI);
    output_sorted_indexes(v->sfaces_begin(), 
			  v->sfaces_end(), SFI);
    out << index(SD.shalfloop()) << " | ";
  }
  else {
    out 
    << index(v->svertices_begin()) << " "
    << index(v->svertices_last()) << ", "
    << index(v->shalfedges_begin()) << " "
    << index(v->shalfedges_last()) << ", "
    << index(v->sfaces_begin()) << " "
    << index(v->sfaces_last()) << ", "
    << index(SD.shalfloop()) << " | ";
  }
  if(reduce)
    out << Infi_box::standard_point(point(v));
  else
    out << point(v);
  out << " } "  << mark(v) << std::endl;
}

/*
template <typename EW>
bool SNC_io_parser<EW>::
read_vertex(Vertex_handle v) const
{ // syntax: index { svs sve ses see sfs sfe sl, mark, point }
  int n, svs, sve, ses, see, sfs, sfe, sl;  Mark m;  Point_3 p;
  if ( !(in >> n) ||
       !check_sep("{") ||
       !(in >> svs) || !(in >> sve) || 
       !(in >> ses) || !(in >> see) || 
       !(in >> sfs) || !(in >> sfe) || 
       !(in >> sl) ||
       !check_sep(",") ||
       !(in >> m) || 
       !check_sep(",") ||
       !(in >> p) || 
       !check_sep("}") ) return false;
  CGAL_assertion_msg(
      Vertex_of[n] == v &&
      svs >= 0 && svs < en && sve >= 0 && sve < en && 
      ses >= 0 && ses < sen && see >= 0 && see < sen && 
      sfs >= 0 && sfs < sfn && sfe >= 0 && sfe < sfn && 
      sl >= 0 && sl < sln , "wrong index in read_vertex");

  v->svertices_begin_ = Edge_of[svs]; v->svertices_last_ = Edge_of[sve];
  v->shalfedges_begin_ = SEdge_of[ses]; v->shalfedges_last_ = SEdge_of[see];
  v->sfaces_begin_ = SFace_of[sfs]; v->sfaces_last_ = SFace_of[sfe];
  v->shalfloop_ = SLoop_of[sl];
  mark(v) = m; point(v) = p;
  return true; 
}
*/

template <typename EW>
template <typename NT>
bool SNC_io_parser<EW>::
read_vertex(Vertex_handle vh) {

  bool OK = true;
  int index;
  NT hx, hy, hz, hw;

  in >> index;
  OK = OK && test_string("{");
  vh->sncp() = sncp();
  
  in >> index;
  vh->svertices_begin() = (index >= 0 ? Edge_of[index] : svertices_end());
  in >> index;
  vh->svertices_last()  = index >= 0 ? Edge_of[index] : svertices_end();
  OK = OK && test_string(",");
  in >> index;
  vh->shalfedges_begin() = index >= 0 ? SEdge_of[index] : shalfedges_end();
  in >> index;
  vh->shalfedges_last()  = index >= 0 ? SEdge_of[index] : shalfedges_end();
  OK = OK && test_string(",");
  in >> index;
  vh->sfaces_begin() = index >= 0 ? SFace_of[index] : sfaces_end();
  in >> index;
  vh->sfaces_last()  = index >= 0 ? SFace_of[index] : sfaces_end();
  OK = OK && test_string(",");
  in >> index;
  vh->shalfloop() = index >= 0 ? SLoop_of[index] : shalfloops_end();
  OK = OK && test_string("|");
  in >> hx >> hy >> hz >> hw;
  vh->point() = Point_3(hx, hy, hz, hw);
  OK = OK && test_string("}");
  in >> vh->mark();
  
  return OK;
}

template <typename EW>
void SNC_io_parser<EW>::print_edge(Halfedge_handle e) const
{ // syntax: index { twin, source, isolated incident_object | spoint } mark
  SM_decorator D(&*vertex(e));
  out << index(e) << " { " << index(twin(e)) << ", "
      << index(source(e)) << ", ";
  if ( D.is_isolated(e) ) out << "1 " << index(D.face(e));
  else out << "0 " << index(D.first_out_edge(e));
  out << " | ";
  if (reduce)
    out << Infi_box::standard_vector(vector(e));
  else
    out << vector(e);
  out << " } "<< mark(e) << std::endl;
}

/*
template <typename EW>
bool SNC_io_parser<EW>::
read_edge(Halfedge_handle e) const
{ // syntax: index { twin, source, isolated incident_object, mark }
  int n, et, vs, ef, efm; bool iso; Mark m;
  if ( !(in >> n) ||
       !check_sep("{") ||
       !(in >> et) || !check_sep(",") ||
       !(in >> vs) || !check_sep(",") ||
       !(in >> iso) || !(in >> ef) || !check_sep(",") ||
       !(in >> m) || !check_sep("}") )
    return false;
  
  if (iso) efm=sfn; else efm=sen;
  CGAL_assertion_msg (
      Edge_of[n] == e &&
      et >= 0 && et < en && vs >= 0 && vs < vn && 
      ef >= 0 && et < efm , "wrong index in read_edge");
  
  e->twin_ = Edge_of[et];
  e->center_vertex() = Vertex_of[vs];
  if ( iso ) e->incident_sface_ = SFace_of[ef];
  else       e->out_sedge() = SEdge_of[ef];
  mark(e) = m;
  return true;
}
*/
template <typename EW>
template <typename NT>
bool SNC_io_parser<EW>::
read_edge(Halfedge_handle eh) {

  bool OK = true;
  int index;
  NT hx, hy, hz, hw;
  
  in >> index;
  OK = OK && test_string("{");
  
  in >> index;
  eh->twin() = Edge_of[index];
  OK = OK && test_string(",");
  in >> index;
  eh->center_vertex() = Vertex_of[index];
  OK = OK && test_string(",");
  in >> index;
  if(index == 0) {
    in >> index;
    eh->out_sedge() = SEdge_of[index];
  } else { 
    in >> index;
    eh->incident_sface() = SFace_of[index];
  }
  OK = OK && test_string("|");
  in >> hx >> hy >> hz >> hw;
  eh->vector() = Sphere_point(hx,hy,hz);
  OK = OK && test_string("}");
  in >> eh->mark();

  return OK;
}

template <typename EW>
void SNC_io_parser<EW>::print_facet(Halffacet_handle f) const
{ // syntax: index { twin, fclist, ivlist, volume | plane } mark
  out << index(f) << " { "; 
  out << index(twin(f)) << ", ";
  Halffacet_cycle_iterator it; 
  CGAL_forall_facet_cycles_of(it,f)
    if ( it.is_shalfedge() ) out << index(SHalfedge_handle(it)) << ' ';
  out << ", ";
  CGAL_forall_facet_cycles_of(it,f)
    if ( it.is_shalfloop() ) out << index(SHalfloop_handle(it)) << ' ';
  out << ", " << index(volume(f)) << " | ";
  if(reduce)
    out << Infi_box::standard_plane(plane(f));
  else 
    out << plane(f);
  out << " } " << mark(f) << std::endl;
}

/*
template <typename EW>
bool SNC_io_parser<EW>::
read_facet(Halffacet_handle f) const
{ // syntax: index { fclist, ivlist, loop, mark }
  int n, ei, li; Mark m;
  if ( !(in >> n) || !check_sep("{") ) return false;
  while (in >> ei) { 
    CGAL_assertion_msg(ei >= 0 && ei < sen, 
      "wrong index in facet cycle list.");
    store_boundary_object(SEdge_of[ei],f);
  } in.clear();
  if (!check_sep(",")) { return false; }
  while (in >> li) { 
    CGAL_assertion_msg(li >= 0 && li < sln, 
      "wrong index in facet cycle list.");
    store_boundary_object(SLoop_of[li],f);
  } in.clear();
  if (!check_sep(",") || !(in >> m) || !check_sep("}") ) 
    return false;
  mark(f) = m;
  return true;
}
*/

template <typename EW>
template <typename NT>
bool SNC_io_parser<EW>::
read_facet(Halffacet_handle fh) {

  bool OK = true;
  int index;
  char cc;
  NT hx, hy, hz, hw;

  in >> index;
  OK = OK && test_string("{");
  
  in >> index;
  fh->twin() = Halffacet_of[index];
  OK = OK && test_string(",");
  
  in >> cc;
  while(isdigit(cc)) {
    in.putback(cc);
    in >> index;
    fh->boundary_entry_objects().push_back(SEdge_of[index]);
    in >> cc;
  }
  
  in >> cc;
  while(isdigit(cc)) {
    in.putback(cc);
    in >> index;
    fh->boundary_entry_objects().push_back(SLoop_of[index]);
    in >> cc;
  }
  
  in >> index;
  fh->volume() = Volume_of[index+addInfiBox];
  OK = OK && test_string("|");
  in >> hx >> hy >> hz >> hw;
  fh->supporting_plane_ = Plane_3(hx,hy,hz,hw);
  OK = OK && test_string("}");
  in >> fh->mark();

  return OK;
}

template <typename EW>
void SNC_io_parser<EW>::print_volume(Volume_handle c) const
{ // syntax: index { shlist } mark
  out << index(c) << " { "; 
  Shell_entry_iterator it;
  CGAL_forall_shells_of(it,c) 
    if(!reduce || Infi_box::is_standard(point(vertex(SFace_handle(it)))))
      out << index(SFace_handle(it)) << ' ';
  out << "} " << mark(c) << std::endl;
}

/*
template <typename EW>
bool SNC_io_parser<EW>::
read_volume(Volume_handle c) const
{ // syntax: index { shlist, mark }
  int n, fi; Mark m;
  if ( !(in >> n) || !check_sep("{") ) return false;
  while (in >> fi) { 
    CGAL_assertion_msg(fi >= 0 && fi < sfn, 
      "wrong index in shell list.");
    store_boundary_object(SFace_of[fi],c);
  } in.clear();
  if (!check_sep(",") || !(in >> m) || !check_sep("}") ) 
    return false;
  mark(c) = m;
  return true;
}
*/

template <typename EW>
template <typename NT>
bool SNC_io_parser<EW>::
read_volume(Volume_handle ch) {

  bool OK = true;
  int index;
  char cc;
  NT hx, hy, hz, hw;
  
  in >> index;
  OK = OK && test_string("{");
  
  in >> cc;
  while(isdigit(cc)) {
    in.putback(cc);
    in >> index;
    ch->shell_entry_objects().push_back(SFace_of[index]);
    in >> cc;
  }
  in >> ch->mark();

  return OK;
}

template <typename EW>
void SNC_io_parser<EW>::
print_sedge(SHalfedge_handle e) const { 
//index { twin, sprev, snext, source, sface, prev, next, facet | circle } mark
  SM_decorator D(&*vertex(e));
  out << index(e) << " { "
      << index(D.twin(e)) << ", " 
      << index(D.previous(e)) << ", " << index(D.next(e)) << ", "
      << index(D.source(e)) << ", " << index(D.face(e)) << ", "
      << index(previous(e)) << ", " << index(next(e)) << ", "
      << index(facet(e)) 
      << " | ";
  if(reduce)
    out << Infi_box::standard_plane((Plane_3) D.circle(e));
  else
    out << D.circle(e);
  out << " } " << D.mark(e) << "\n";
}

/*
template <typename EW>
bool SNC_io_parser<EW>::
read_sedge(SHalfedge_handle e) const
{ // syntax: index { twin, sprev, snext, source, sface, prev, next, facet }
  int n, et, sp, sn, vs, sf, ep, en, ft; 
  if ( !(in >> n) ||
       !check_sep("{") ||
       !(in >> et) || !check_sep(",") ||
       !(in >> sp) || !check_sep(",") ||
       !(in >> sn) || !check_sep(",") ||
       !(in >> vs) || !check_sep(",") ||
       !(in >> sf) || !check_sep(",") ||
       !(in >> ep) || !check_sep(",") ||
       !(in >> en) || !check_sep(",") ||
       !(in >> ft) || !check_sep("}") )
    return false;
  CGAL_assertion_msg 
     (et >= 0 && et < sen && sp >= 0 && sp < sen && 
      sn >= 0 && sn < sen && vs >= 0 && vs < en && 
      sf >= 0 && sf < sfn && ep >= 0 && ep < sen &&
      en >= 0 && en < sen && ft >= 0 && ft < fn ,
      "wrong index in read_sedge");
  
  // precond: features exist!
  CGAL_assertion(SEdge_of[n]==e);
  e->sprev() = SEdge_of[sp];
  e->snext() = SEdge_of[sn];
  e->source() = Edge_of[vs];
  e->incident_sface() = SFace_of[sf];
  e->prev() = SEdge_of[ep];
  e->next() = SEdge_of[en];
  e->incident_facet() = Halffacet_of[ft];
  return true;
}
*/

template <typename EW>
template <typename NT>
bool SNC_io_parser<EW>::
read_sedge(SHalfedge_handle seh) {

  bool OK = true;
  int index;
  NT hx, hy, hz, hw;
  
  in >> index;
  OK = OK && test_string("{");
  
  in >> index;
  seh->twin() = SEdge_of[index];
  OK = OK && test_string(",");
  in >> index;
  seh->sprev() = SEdge_of[index];
  OK = OK && test_string(",");
  in >> index;
  seh->snext() = SEdge_of[index];
  OK = OK && test_string(",");
  in >> index;
  seh->source() = Edge_of[index];
  OK = OK && test_string(",");
  in >> index;
  seh->incident_sface() = SFace_of[index];
  OK = OK && test_string(",");
  in >> index;
  seh->prev() = SEdge_of[index];
  OK = OK && test_string(",");
  in >> index;
  seh->next() = SEdge_of[index];
  OK = OK && test_string(",");
  in >> index;
  seh->incident_facet() = Halffacet_of[index];
  OK = OK && test_string("|");
  in >> hx >> hy >> hz >> hw;
  seh->circle() = Sphere_circle(Plane_3(hx,hy,hz,hw));
  OK = OK && test_string("}");
  in >> seh->mark();

  return OK;
}

template <typename EW>
void SNC_io_parser<EW>::
print_sloop(SHalfloop_handle l) const
{ // syntax: index { twin, sface, facet | circle } mark
  SM_decorator D(&*vertex(l));
  out << index(l) << " { "
      << index(D.twin(l)) << ", " << index(D.face(l)) << ", "
      << index(facet(l)) 
      << " | ";  
  if(reduce)
    out << Infi_box::standard_plane((Plane_3) D.circle(l));
  else
    out << D.circle(l);
  out << " } " << D.mark(l) << "\n";
}

/*
template <typename EW>
bool SNC_io_parser<EW>::
read_sloop(SHalfloop_handle l) const
{ // syntax: index { twin, sface, facet }
  int n, lt, sf, ft; 
  if ( !(in >> n) ||
       !check_sep("{") ||
       !(in >> lt) || !check_sep(",") ||
       !(in >> sf) || !check_sep(",") ||
       !(in >> ft) || !check_sep("}") )
    return false;
  CGAL_assertion_msg 
     (lt >= 0 && lt < sen && sf >= 0 && sf < sfn && 
      ft >= 0 && ft < fn ,
      "wrong index in read_sedge");
  
  CGAL_assertion(SLoop_of[n]==l);
  l->twin_ = SLoop_of[lt];
  l->incident_sface_ = SFace_of[sf];
  l->incident_facet_ = Halffacet_of[ft];
  return true;
}
*/

template <typename EW>
template <typename NT>
bool SNC_io_parser<EW>::
read_sloop(SHalfloop_handle slh) {

  bool OK = true;
  int index;
  NT hx, hy, hz, hw;

  in >> index;
  OK = OK && test_string("{");
  
  in >> index;
  slh->twin() = SLoop_of[index];
  OK = OK && test_string(",");
  in >> index;
  slh->incident_sface() = SFace_of[index];
  OK = OK && test_string(",");
  in >> index;
  slh->incident_facet() = Halffacet_of[index];
  OK = OK && test_string("|");	
  in >> hx >> hy >> hz >> hw;
  slh->circle() = Sphere_circle(Plane_3(hx,hy,hz,hw));	
  OK = OK && test_string("}");	
  in >> slh->mark();
  
  return OK;
}

template <typename EW>
void SNC_io_parser<EW>::
print_sface(SFace_handle f) const
{ // syntax: index { vertex, fclist, ivlist, sloop, volume }
  SM_decorator D(&*f->center_vertex());
  out << index(f) << " { " << index(f->center_vertex()) << ", "; 
  SFace_cycle_iterator it;
  CGAL_forall_sface_cycles_of(it,f)
    if ( it.is_shalfedge() ) out << index(SHalfedge_handle(it)) << ' ';
  out << ", ";
  CGAL_forall_sface_cycles_of(it,f)
    if ( it.is_svertex() ) out << index(SVertex_handle(it)) << ' ';
  out << ", ";
  CGAL_forall_sface_cycles_of(it,f)
    if ( it.is_shalfloop() ) out << index(SHalfloop_handle(it));
  out << ", " << index(volume(f)) << " } " << D.mark(f) <<"\n";
}

/*
template <typename EW>
bool SNC_io_parser<EW>::
read_sface(SFace_handle f) const
{ // syntax: index { vertex, fclist, ivlist, sloop, volume }
  int n, vc, ei, vi, li, c;
  if ( !(in >> n) || !check_sep("{") ||
       !(in >> vc) || !check_sep(",") ) 
    return false;
  CGAL_assertion(vc >= 0 && vc < vn);
  f->center_vertex() = Vertex_of[vc];
  SM_decorator D(&*Vertex_of[vc]);
  while (in >> ei) { 
    CGAL_assertion_msg(ei >= 0 && ei < sen, 
      "wrong index in sface cycle list.");
    D.store_boundary_object(SEdge_of[ei],f);
  } in.clear();
  while (in >> vi) { 
    CGAL_assertion_msg(vi >= 0 && vi < en, 
      "wrong index in sface cycle list.");
    D.store_boundary_object(Edge_of[vi],f);
  } in.clear();
  if (!check_sep(",")) { return false; }
  while (in >> li) { 
    CGAL_assertion_msg(li >= 0 && li < sln, 
      "wrong index in sface cycle list.");
    D.store_boundary_object(SLoop_of[li],f);
  } in.clear();
  if (!check_sep(",") || !(in >> c) || !check_sep("}") ) 
    return false;
  f->incident_volume_ = Volume_of[c];
  return true;
}
*/

template <typename EW>
template <typename NT>
bool SNC_io_parser<EW>::
read_sface(SFace_handle sfh) {

  bool OK = true;
  int index;
  char cc;
  NT hx, hy, hz, hw;

  in >> index;
  OK = OK && test_string("{");
  
  in >> index;
  sfh->center_vertex() = Vertex_of[index];
  OK = OK && test_string(",");
  
  in >> cc;
  while(isdigit(cc)) {
    in.putback(cc);
    in >> index;
    sfh->boundary_entry_objects().push_back(SEdge_of[index]);
    in >> cc;
  }
  
  in >> cc;
  while(isdigit(cc)) {
    in.putback(cc);
    in >> index;
    sfh->boundary_entry_objects().push_back(Edge_of[index]);
    in >> cc;
  }

  in >> cc;
  while(isdigit(cc)) {
    in.putback(cc);
    in >> index;
    sfh->boundary_entry_objects().push_back(SLoop_of[index]);
    in >> cc;
  }
  
  in >> index;
  sfh->incident_volume() = Volume_of[index+addInfiBox];
  OK = OK && test_string("}");	
  in >> sfh->mark();
  
  return OK;
}

template <typename EW>
void SNC_io_parser<EW>::print_local_graph(Vertex_handle v) const
{ SM_decorator D(&*v);
  out << "Local Graph " 
      << D.number_of_vertices() << " " << D.number_of_edges() << " "
      << D.number_of_loops() << " " << D.number_of_faces() << " "
      << std::endl;
  if (verbose) 
    out << "/* index { twin, source, isolated incident_object, mark } */\n";
  SVertex_iterator vit;
  CGAL_forall_svertices_of(vit,v) print_edge(vit);
  if (verbose) 
    out << "/* index { twin, sprev, snext, source, sface,"
        << " prev, next, facet } */\n";
  SHalfedge_iterator eit;
  CGAL_forall_shalfedges_of(eit,v) print_sedge(eit);
  if (verbose) 
    out << "/* index { twin, sface, facet } */" << std::endl;
  if ( D.has_sloop() ) 
  { print_sloop(D.loop()); print_sloop(D.twin(D.loop())); }
  if (verbose) 
    out << "/* index { fclist, ivlist, sloop, volume } */" << std::endl;
  SFace_iterator fit;
  CGAL_forall_sfaces_of(fit,v) print_sface(fit);
  out.flush();
}

template <typename EW>
void SNC_io_parser<EW>::add_infi_box() {

  typename Standard_kernel::RT hx,hy,hz,hw;
  for(int i=0; i<8; i++) {
    Vertex_handle vh = Vertex_of[vn+i];
    vh->svertices_begin() = Edge_of[en+3*i];
    vh->svertices_last()  = Edge_of[en+3*i+2];
    vh->shalfedges_begin() = SEdge_of[sen+6*i];
    vh->shalfedges_last()  = SEdge_of[sen+6*i+5];
    vh->sfaces_begin() = SFace_of[sfn+2*i];
    vh->sfaces_last()  = SFace_of[sfn+2*i+1];
    vh->shalfloop() = shalfloops_end();
    hx = i % 2 ? -1 : 1;
    hy = i % 4 > 1 ? -1 : 1;
    hz = i > 3 ? -1 : 1;
    vh->point() = Infi_box::create_extended_point(hx, hy, hz);
    vh->mark() = 1;		
    vh->sncp() = sncp();
  }
  
  int seOff[3] = {0, 1, 3};
  int twinIdx[24] = { 3, 7,14,
		      0,10,17,
		      9, 1,20,
		      6, 4,23,
		      15,19, 2,
		      12,22, 5,
		      21,13, 8,
		      18,16,11};
  
  for(int i = 0; i < 24; i++) {
    Halfedge_handle eh = Edge_of[en+i];
    eh->twin() = Edge_of[en+twinIdx[i]];
    eh->center_vertex() = Vertex_of[vn+(i/3)];
    eh->out_sedge() = SEdge_of[sen+(i/3*6)+seOff[i%3]];
    switch(i%3) {
    case 0 : 
      hx = i % 6 ? 1 : -1;
      hy = hz = 0;
      break;
    case 1:
      hy = i % 12 >= 6 ? 1 : -1;
      hx = hz = 0;
      break;
    case 2:
      hz = i >= 12 ? 1 : -1;
      hx = hy = 0;
      break;
    }
    eh->vector() = Sphere_point(hx,hy,hz);
    eh->mark() = 1;
  }
  
  int bnd[12] = {19, 18, 43, 42, 35, 34,
		 47, 46, 39, 38, 45, 44};
  for(int i = 0; i < 12; i++) {
    Halffacet_handle fh = Halffacet_of[fn+i];
    fh->twin() = Halffacet_of[fn+(i/2*2)+((i+1)%2)];
    fh->boundary_entry_objects().push_back(SEdge_of[sen+bnd[i]]);
    fh->volume() = Volume_of[((i%4) == 1 || (i%4 == 2)) ? 1 : 0];
    if(i<4) {
      hz = i % 2 ? -1 : 1;
      hx = hy = 0;
    }
    else if(i<8) {
      hy = i % 2 ? -1 : 1;
      hx = hz = 0;
    }
    else {
      hx = i % 2 ? -1 : 1;
      hz = hy = 0;
    }
    hw = ((i%4) == 1 || (i%4) == 2) ? 1 : -1;
    fh->supporting_plane_ = Infi_box::create_extended_plane(hx,hy,hz,hw);
    fh->mark() = 1;
  }
  
  Volume_of[0]->shell_entry_objects().push_back(SFace_of[sfn]);
  Volume_of[0]->mark() = 0;
  Volume_of[1]->shell_entry_objects().push_front(SFace_of[sfn+1]);
  
  int sprevOff[6] = {4,3,0,5,2,1};
  int snextOff[6] = {2,5,4,1,0,3};
  int prevIdx[48] = {7,12,15,26,29,10,
		     1,18,21,32,35,4,
		     19,0,3,38,41,22,
		     13,6,9,44,47,16,
		     31,36,39,2,5,34,
		     25,42,45,8,11,28,
		     43,24,27,14,17,46,
		     37,30,33,20,23,40};
  int nextIdx[48] = {13,6,27,14,11,28,
		     19,0,33,20,5,34,
		     1,18,39,2,23,40,
		     7,12,45,8,17,46,
		     37,30,3,38,35,4,
		     43,24,9,44,29,10,
		     25,42,15,26,47,16,
		     31,36,21,32,41,22};
  int factIdx[48] = {1,0,9,8,5,4,
		     0,1,11,10,4,5,
		     0,1,8,9,7,6,
		     1,0,10,11,6,7,
		     3,2,8,9,4,5,
		     2,3,10,11,5,4,
		     2,3,9,8,6,7,
		     3,2,11,10,7,6};
  int sgn[24] = {1,1,1,-1,1,-1,
		 -1,-1,1,1,-1,-1,
		 1,-1,-1,-1,-1,1,
		 -1,1,-1,1,1,1};
  
  for(int i = 0; i < 48; i++) {
    SHalfedge_handle seh = SEdge_of[sen+i];
    
    seh->twin() = SEdge_of[sen+(i/2*2)+((i+1)%2)];
    seh->sprev() = SEdge_of[sen+sprevOff[i%6]+(i/6*6)];
    seh->snext() = SEdge_of[sen+snextOff[i%6]+(i/6*6)];
    seh->source() = Edge_of[en+((i+1)%6)/2+(i/6)*3];
    seh->incident_sface() = SFace_of[sfn+(i%2)+(i/6)*2];
    seh->prev() = SEdge_of[sen+prevIdx[i]];
    seh->next() = SEdge_of[sen+nextIdx[i]];
    seh->incident_facet() = Halffacet_of[fn+factIdx[i]];
    if(i%6 < 2) {
      hz = (i%2) ? sgn[i/2] * (-1) : sgn[i/2];
      hx = hy = 0;
    }
    else if(i%6 < 4) {
      hx = (i%2) ? sgn[i/2] * (-1) : sgn[i/2];
      hz = hy = 0;
    }
    else {
      hy = (i%2) ? sgn[i/2] * (-1) : sgn[i/2];
      hx = hz = 0;
    }
    seh->circle() = Sphere_circle(Plane_3(RT(hx),RT(hy),RT(hz),RT(0)));
    seh->mark() = 1;
  }
  
  int volIdx[8] = {0,1,1,0,1,0,0,1};
  
  for(int i = 0; i < 16; i++) {
    SFace_handle sfh = SFace_of[sfn+i];
    sfh->center_vertex() = Vertex_of[vn+(i/2)];
    sfh->boundary_entry_objects().push_back(SEdge_of[sen+(i/2*6)+(i%2)]);
    int cIdx = i%2 ? 1-volIdx[i/2] : volIdx[i/2];
    sfh->incident_volume() = Volume_of[cIdx];
    sfh->mark() = cIdx ? Volume_of[1]->mark() : 0;
  }
}

CGAL_END_NAMESPACE
#endif //CGAL_SNC_IO_PARSER_H

