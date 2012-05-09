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
// Author(s)     : Michael Seel    <seel@mpi-sb.mpg.de>
//                 Miguel Granados <granados@mpi-sb.mpg.de>
//                 Susan Hert      <hert@mpi-sb.mpg.de>
//                 Lutz Kettner    <kettner@mpi-sb.mpg.de>
#ifndef CGAL_SNC_STRUCTURE_H
#define CGAL_SNC_STRUCTURE_H

#include <CGAL/basic.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/In_place_list.h>
#include <CGAL/Nef_3/SNC_list.h>
#include <CGAL/Nef_S2/Generic_handle_map.h>
#include <CGAL/Nef_2/Object_handle.h>
#include <CGAL/Nef_3/SNC_iteration.h>
#include <CGAL/Nef_2/iterator_tools.h>
#include <CGAL/Nef_S2/Sphere_geometry.h>
#ifdef CGAL_NEF3_FACET_WITH_BOX
#include <CGAL/Box_intersection_d/Box_d.h>
#endif
#include <list>

#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 41
#include <CGAL/Nef_2/debug.h>
#include <CGAL/Nef_2/Object_index.h>

#include <boost/optional.hpp>
#include <boost/none.hpp>

namespace CGAL {

template <typename HE>
struct move_shalfedge_around_facet {
  void forward(HE& e) const { e = (e->next()); }
  void backward(HE& e) const { e = (e->prev()); }
};

template <class Object, class Hash_map, class Union_find>
void merge_sets( Object o1, Object o2, Hash_map& hash, Union_find& uf) {
  //CGAL_assertion( hash[o1] != 0 && hash[o2] != 0);
  CGAL_assertion( hash.is_defined(o1) && hash.is_defined(o2));
  if( !uf.same_set( hash[o1], hash[o2]))
    uf.unify_sets( hash[o1], hash[o2]);
}

template <typename K, typename I, typename M> class SNC_sphere_map;
template <typename S> class SM_decorator;
template <typename S> class SNC_decorator;

/*{\Manpage {SNC_structure}{Items}{Selective Nef Complex}{C}}*/

template <typename Kernel_, typename Items_, typename Mark_>
class SNC_structure {
/*{\Mdefinition The extended Wuerzburg structure is the topological
structure of Nef polyhedra. It is programmed around the local
graphs of the vertices of a Nef polyhedron, which describe the
point set completely. All other concepts are either derived from
the local graph or added for the comfort of the user.}*/

public:
  /*{\Mtypes 7}*/

  typedef Items_                             Items;
  typedef Kernel_                            Kernel;
  typedef Mark_                              Mark;
  typedef SNC_structure<Kernel,Items,Mark>   Self;

  //  typedef bool                          Mark;
  typedef CGAL::SNC_decorator<Self>           SNC_decorator;

  typedef typename Kernel::FT           FT;
  typedef typename Kernel::RT           RT;
  typedef CGAL::Sphere_geometry<Kernel> Sphere_kernel;
  
  typedef SNC_sphere_map<Kernel, Items, Mark> Sphere_map;
  typedef CGAL::SM_decorator<Sphere_map>  SM_decorator;

  typedef typename Kernel::Point_3      Point_3;
  /*{\Mtypemember embedding vertices.}*/
  typedef typename Kernel::Plane_3      Plane_3;
  /*{\Mtypemember supporting facets.}*/
  typedef typename Kernel::Vector_3     Vector_3;
  /*{\Mtypemember normal vectors.}*/
  typedef typename Kernel::Direction_3  Direction_3;
  /*{\Mtypemember normal directions.}*/
  typedef typename Kernel::Segment_3    Segment_3;
  /*{\Mtypemember segments in space.}*/
  typedef typename Kernel::Line_3       Line_3;
  /*{\Mtypemember lines in space.}*/
  typedef typename Kernel::Ray_3       Ray_3;
  /*{\Mtypemember rays in space.}*/
  typedef typename Kernel::Triangle_3       Triangle_3;
  /*{\Mtypemember triangles in space.}*/

  typedef typename Kernel::Aff_transformation_3 Aff_transformation_3;

  typedef typename Sphere_kernel::Sphere_point     Sphere_point;
  /*{\Mtypemember points on the unit sphere.}*/
  typedef typename Sphere_kernel::Sphere_segment   Sphere_segment;
  /*{\Mtypemember segments on the unit sphere.}*/
  typedef typename Sphere_kernel::Sphere_circle    Sphere_circle;
  /*{\Mtypemember segments on the unit sphere.}*/
  typedef typename Sphere_kernel::Sphere_direction Sphere_direction;
  /*{\Mtypemember directions on the unit sphere.}*/
  typedef size_t Size_type;
  /*{\Mtypemember size type.}*/

 private:
  friend class CGAL::SM_decorator<Sphere_map>;
  friend class CGAL::SNC_decorator<Self>;

  /*{\Mtext For all objects |Vertex|, |Halfedge|, |Halffacet|, |Volume|
  there are handle and iterator types |xxx_handle|, |xxx_iterator|.
  Additionally all objects of the local graph of a vertex 
  |SVertex|, |SHalfedge|, |SHalfloop|, |SFace| are accessed via handles
  and iterators. There's no type |SHalfloop_iterator|, as there is
  at most one |SLoop| pair per vertex.}*/

#ifdef CGAL_NEF3_FACET_WITH_BOX
  template<class Refs>
  class Facet_with_box : public Items::template Halffacet<Refs> {
  public:
    typedef typename Items::template Halffacet<Refs> Halffacet;
    typedef typename Refs::FT FT;
    typedef typename Box_intersection_d::Box_d<FT,3> Box;

    Box b;
    
    Facet_with_box() : Halffacet(), b() {}     
    Facet_with_box(const Plane_3& h, Mark m) : Halffacet(h,m) {}
    Facet_with_box(const Facet_with_box<Refs>& f) : Halffacet(f) {
      b = f.b;
    }
    
    Facet_with_box<Refs>& operator=(const Facet_with_box<Refs>& f) {
      (Halffacet) *this = (Halffacet) f;
      b = f.b;
      return *this;
    }
  };
#endif

 public:
  typedef Sphere_map                                        Vertex_base;
  typedef SNC_in_place_list_sm<Vertex_base>                 Vertex; 
  typedef CGAL::In_place_list<Vertex,false>                 Vertex_list;
  typedef CGAL_ALLOCATOR(Vertex)                            Vertex_alloc;
  typedef typename Vertex_list::iterator                    Vertex_handle;
  typedef typename Vertex_list::const_iterator              Vertex_const_handle;
  typedef typename Vertex_list::iterator                    Vertex_iterator;
  typedef typename Vertex_list::const_iterator              Vertex_const_iterator;

#ifdef CGAL_NEF3_FACET_WITH_BOX
  typedef Facet_with_box<SNC_structure>                     Halffacet_base;
#else
  typedef typename Items::template Halffacet<SNC_structure> Halffacet_base;
#endif
  typedef SNC_in_place_list_halffacet<Halffacet_base>       Halffacet;
  typedef CGAL::In_place_list<Halffacet,false>              Halffacet_list;
  typedef CGAL_ALLOCATOR(Halffacet)                         Halffacet_alloc;
  typedef typename Halffacet_list::iterator                 Halffacet_handle;
  typedef typename Halffacet_list::const_iterator           Halffacet_const_handle;
  typedef typename Halffacet_list::iterator                 Halffacet_iterator;
  typedef typename Halffacet_list::const_iterator           Halffacet_const_iterator;

  typedef typename Items::template Volume<SNC_structure>    Volume_base;
  typedef SNC_in_place_list_volume<Volume_base>             Volume;
  typedef CGAL::In_place_list<Volume,false>                 Volume_list;
  typedef CGAL_ALLOCATOR(Volume)                            Volume_alloc;
  typedef typename Volume_list::iterator                    Volume_handle;
  typedef typename Volume_list::const_iterator              Volume_const_handle;
  typedef typename Volume_list::iterator                    Volume_iterator;
  typedef typename Volume_list::const_iterator              Volume_const_iterator;

  typedef typename Items::template SVertex<SNC_structure>   SVertex_base;
  typedef SNC_in_place_list_svertex<SVertex_base>           SVertex;
  typedef CGAL::In_place_list<SVertex,false>                SVertex_list;
  typedef CGAL_ALLOCATOR(SVertex)                           SVertex_alloc;
  typedef typename SVertex_list::iterator                   SVertex_handle;
  typedef typename SVertex_list::const_iterator             SVertex_const_handle;
  typedef typename SVertex_list::iterator                   SVertex_iterator;
  typedef typename SVertex_list::const_iterator             SVertex_const_iterator;

  typedef typename Items::template SVertex<SNC_structure>   Halfedge_base;
  typedef SNC_in_place_list_svertex<SVertex_base>           Halfedge;
  typedef CGAL::In_place_list<SVertex,false>                Halfedge_list;
  typedef CGAL_ALLOCATOR(SVertex)                           Halfedge_alloc;
  typedef typename SVertex_list::iterator                   Halfedge_handle;
  typedef typename SVertex_list::const_iterator             Halfedge_const_handle;
  typedef typename SVertex_list::iterator                   Halfedge_iterator;
  typedef typename SVertex_list::const_iterator             Halfedge_const_iterator;

  typedef typename Items::template SHalfedge<SNC_structure> SHalfedge_base;
  typedef SNC_in_place_list_shalfedge<SHalfedge_base>       SHalfedge;
  typedef CGAL::In_place_list<SHalfedge,false>              SHalfedge_list;
  typedef CGAL_ALLOCATOR(SHalfedge)                         SHalfedge_alloc;
  typedef typename SHalfedge_list::iterator                 SHalfedge_handle;
  typedef typename SHalfedge_list::const_iterator           SHalfedge_const_handle;
  typedef typename SHalfedge_list::iterator                 SHalfedge_iterator;
  typedef typename SHalfedge_list::const_iterator           SHalfedge_const_iterator;

  typedef typename Items::template SHalfloop<SNC_structure> SHalfloop_base;
  typedef SNC_in_place_list_shalfloop<SHalfloop_base>       SHalfloop;
  typedef CGAL::In_place_list<SHalfloop,false>              SHalfloop_list;
  typedef CGAL_ALLOCATOR(SHalfloop)                         SHalfloop_alloc;
  typedef typename SHalfloop_list::iterator                 SHalfloop_handle;
  typedef typename SHalfloop_list::const_iterator           SHalfloop_const_handle;
  typedef typename SHalfloop_list::iterator                 SHalfloop_iterator;
  typedef typename SHalfloop_list::const_iterator           SHalfloop_const_iterator;

  typedef typename Items::template SFace<SNC_structure>     SFace_base;
  typedef SNC_in_place_list_sface<SFace_base>               SFace;
  typedef CGAL::In_place_list<SFace,false>                  SFace_list;
  typedef CGAL_ALLOCATOR(SFace)                             SFace_alloc;
  typedef typename SFace_list::iterator                     SFace_handle;
  typedef typename SFace_list::const_iterator               SFace_const_handle;
  typedef typename SFace_list::iterator                     SFace_iterator;
  typedef typename SFace_list::const_iterator               SFace_const_iterator;

  typedef CGAL::Object_handle         Object_handle;
  typedef std::list<Object_handle>    Object_list;
  typedef Object_list::iterator       Object_iterator;
  typedef Object_list::const_iterator Object_const_iterator;
  typedef Object_list::const_iterator Object_const_handle;

  typedef typename Sphere_map::SHalfedge_around_svertex_circulator
                               SHalfedge_around_svertex_circulator;
  typedef typename Sphere_map::SHalfedge_around_sface_circulator
                               SHalfedge_around_sface_circulator;
  typedef typename Sphere_map::SHalfedge_around_svertex_const_circulator
                               SHalfedge_around_svertex_const_circulator;
  typedef typename Sphere_map::SHalfedge_around_sface_const_circulator
                               SHalfedge_around_sface_const_circulator;

  typedef typename Sphere_map::Infi_box   Infi_box;
  typedef typename Infi_box::Standard_kernel  Standard_kernel;

  typedef Vertex_handle Constructor_parameter;
  typedef Vertex_const_handle Constructor_const_parameter;

  // Halffacet triangle

#ifdef CGAL_NEF_LIST_OF_TRIANGLES
  class Halffacet_triangle_handle : public Halffacet_handle {
    typedef Halffacet_handle Base;
    Triangle_3* triangle;
  public:
    Halffacet_triangle_handle() : Base() {}
    Halffacet_triangle_handle( Halffacet_handle h, Triangle_3& t) :
      Base(h), triangle(&t) {}
    Triangle_3& get_triangle() { return *triangle; }
    void transform(const Aff_transformation_3& t) { 
      *triangle = Triangle_3((*triangle)[0].transform(t),
    	 		     (*triangle)[1].transform(t),
		 	     (*triangle)[2].transform(t)); 
    }
  };
#else
  class Halffacet_triangle_const_handle : public Halffacet_const_handle {
    typedef Halffacet_const_handle Base;
    Triangle_3 triangle;
  public:
    Halffacet_triangle_const_handle() : Base() {}
    Halffacet_triangle_const_handle( Halffacet_const_handle h, Triangle_3& t) :
      Base(h), triangle(t) {}
    Triangle_3 get_triangle() { return triangle; }
    void transform(const Aff_transformation_3& t) { 
      triangle = Triangle_3(triangle[0].transform(t),
    	 		     triangle[1].transform(t),
		 	     triangle[2].transform(t)); 
    }
  };

  class Halffacet_triangle_handle : public Halffacet_handle {
    typedef Halffacet_handle Base;
    Triangle_3 triangle;
  public:
    Halffacet_triangle_handle() : Base() {}
    Halffacet_triangle_handle( Halffacet_handle h, Triangle_3& t) :
      Base(h), triangle(t) {}
    Triangle_3 get_triangle() { return triangle; }
    void transform(const Aff_transformation_3& t) { 
      triangle = Triangle_3(triangle[0].transform(t),
    	 		     triangle[1].transform(t),
		 	     triangle[2].transform(t)); 
    }
  };
#endif
  
  class Halffacet_cycle_iterator : public Object_iterator 
  /*{\Mtypemember a generic handle to an object in the boundary
  of a facet. Convertible to |Object_handle|.}*/
  { typedef Object_iterator Ibase;
  public:
    Halffacet_cycle_iterator() : Ibase() {}
    Halffacet_cycle_iterator(const Ibase& b) : Ibase(b) {}
    Halffacet_cycle_iterator(const Halffacet_cycle_iterator& i) 
      : Ibase(i) {}
    bool is_shalfedge() const
    { SHalfedge_handle e; return CGAL::assign(e,Ibase::operator*()); }
    bool is_shalfloop() const
    { SHalfloop_handle l; return CGAL::assign(l,Ibase::operator*()); }
 
   operator SHalfedge_handle() const 
    { SHalfedge_handle e; CGAL::assign(e,Ibase::operator*()); return e; }
    operator SHalfloop_handle() const 
    { SHalfloop_handle l; CGAL::assign(l,Ibase::operator*()); return l; }

    operator Object_handle() const { return Ibase::operator*(); }
    Object_handle& operator*() const { return Ibase::operator*(); }
    Object_handle  operator->() const 
      { CGAL_error_msg("not impl."); return Object_handle();}
  };

  class Halffacet_cycle_const_iterator : public Object_const_iterator 
  /*{\Mtypemember a generic handle to an object in the boundary
  of a facet. Convertible to |Object_handle|.}*/
  { typedef Object_const_iterator Ibase;
  public:
    Halffacet_cycle_const_iterator() : Ibase() {}
    Halffacet_cycle_const_iterator(const Ibase& b) : Ibase(b) {}
    Halffacet_cycle_const_iterator(const Halffacet_cycle_const_iterator& i) 
      : Ibase(i) {}  
    bool is_shalfedge() const
    { SHalfedge_handle e; return CGAL::assign(e,Ibase::operator*()); }
    bool is_shalfloop() const
    { SHalfloop_handle l; return CGAL::assign(l,Ibase::operator*()); }
 
   operator SHalfedge_const_handle() const 
    { SHalfedge_handle e; CGAL::assign(e,Ibase::operator*()); 
      return SHalfedge_const_handle(e); }
    operator SHalfloop_const_handle() const 
    { SHalfloop_handle l; CGAL::assign(l,Ibase::operator*()); 
      return SHalfloop_const_handle(l); }

    operator Object_handle() const { return Ibase::operator*(); }
    const Object_handle& operator*() const { return Ibase::operator*(); }
    Object_handle  operator->() const 
    { CGAL_error_msg("not impl."); return Object_handle();}
  };

  class SFace_cycle_iterator : public Object_iterator 
  /*{\Mtypemember a generic iterator to an object in the boundary
  of a sface. Convertible to |Object_handle|.}*/
  { typedef Object_iterator Ibase;
  public:
    SFace_cycle_iterator() : Ibase() {}
    SFace_cycle_iterator(const Ibase& b) : Ibase(b) {}
    SFace_cycle_iterator(const SFace_cycle_iterator& i) : Ibase(i) {}  
    bool is_svertex() const 
    { SVertex_handle v; return CGAL::assign(v,Ibase::operator*()); }
    bool is_shalfedge() const
    { SHalfedge_handle e; return CGAL::assign(e,Ibase::operator*()); }
    bool is_shalfloop() const
    { SHalfloop_handle l; return CGAL::assign(l,Ibase::operator*()); }
    operator SVertex_handle() const 
    { SVertex_handle v; CGAL::assign(v,Ibase::operator*()); return v; }
    operator SHalfedge_handle() const 
    { SHalfedge_handle e; CGAL::assign(e,Ibase::operator*()); return e; }
    operator SHalfloop_handle() const 
    { SHalfloop_handle l; CGAL::assign(l,Ibase::operator*()); return l; }

    operator Object_handle() const { return Ibase::operator*(); }
    Object_handle& operator*() const { return Ibase::operator*(); }
    Object_handle  operator->() const 
    {  return Object_handle(); }
  };

  class SFace_cycle_const_iterator : public Object_const_iterator 
  /*{\Mtypemember a generic iterator to an object in the boundary
  of a facet. Convertible to |Object_handle|.}*/
  { typedef Object_const_iterator Ibase;
  public:
    SFace_cycle_const_iterator() : Ibase() {}
    SFace_cycle_const_iterator(const Ibase& b) : Ibase(b) {}
    SFace_cycle_const_iterator(const SFace_cycle_const_iterator& i) 
      : Ibase(i) {}  
    bool is_svertex() const 
    { SVertex_handle v; return CGAL::assign(v,Ibase::operator*()); }
    bool is_shalfedge() const
    { SHalfedge_handle e; return CGAL::assign(e,Ibase::operator*()); }
    bool is_shalfloop() const
    { SHalfloop_handle l; return CGAL::assign(l,Ibase::operator*()); }
    operator SVertex_const_handle() const 
    { SVertex_handle v; CGAL::assign(v,Ibase::operator*()); 
      return SVertex_const_handle(v); }
    operator SHalfedge_const_handle() const 
    { SHalfedge_handle e; CGAL::assign(e,Ibase::operator*()); 
      return SHalfedge_const_handle(e); }
    operator SHalfloop_const_handle() const 
    { SHalfloop_handle l; CGAL::assign(l,Ibase::operator*()); 
      return SHalfloop_const_handle(l); }

    operator Object_handle() const { return Ibase::operator*(); }
    const Object_handle& operator*() const { return Ibase::operator*(); }
    Object_handle  operator->() const 
    { return Object_handle(); }
  };


  class Shell_entry_iterator : public Object_iterator 
  /*{\Mtypemember a generic iterator to an object in the boundary
  of a shell. Convertible to |SFace_handle|.}*/
  { typedef Object_iterator Ibase; 
  public:
    Shell_entry_iterator() : Ibase() {}
    Shell_entry_iterator(const Ibase& b) : Ibase(b) {}
    Shell_entry_iterator(const Shell_entry_iterator& i) : Ibase(i) {}  

    operator SFace_handle() const 
    { SFace_handle f; 
      CGAL_assertion( CGAL::assign(f,Ibase::operator*()) );
      CGAL::assign(f,Ibase::operator*()); return f; }

    operator Object_handle() const { return Ibase::operator*(); }
    Object_handle& operator*() const { return Ibase::operator*(); }
    Object_handle  operator->() const 
    {  return Object_handle(); }
  };

  class Shell_entry_const_iterator : public Object_const_iterator 
  { typedef Object_const_iterator Ibase; 
  public:
    Shell_entry_const_iterator() : Ibase() {}
    Shell_entry_const_iterator(const Ibase& b) : Ibase(b) {}
    Shell_entry_const_iterator(const Shell_entry_const_iterator& i) :
      Ibase(i) {}  

    operator SFace_const_handle() const 
    { SFace_handle f; 
      CGAL_assertion( CGAL::assign(f,Ibase::operator*()) );
      CGAL::assign(f,Ibase::operator*()); 
      return SFace_const_handle(f); }

    operator Object_handle() const { return Ibase::operator*(); }
    const Object_handle& operator*() const { return Ibase::operator*(); }
    Object_handle  operator->() const 
    {  return Object_handle(); }
  };

  typedef CircFromIt<SHalfedge_const_iterator, 
          move_shalfedge_around_facet<SHalfedge_const_iterator> > 
          SHalfedge_around_facet_const_circulator;

  // Mutable Circulators: 

  typedef CircFromIt<SHalfedge_iterator, 
          move_shalfedge_around_facet<SHalfedge_iterator> > 
          SHalfedge_around_facet_circulator;

#ifdef CGAL_NEF3_FACET_WITH_BOX
  typedef std::pair<SHalfedge_around_facet_circulator,
                    SHalfedge_around_facet_circulator> Outer_cycle;
  typedef std::pair<SHalfedge_around_facet_circulator,
                    SHalfedge_around_facet_circulator> Inner_cycle;
  
  class Partial_facet {
  public:
    Halffacet_handle f;

    std::list<Outer_cycle> outer_cycles;
    std::list<Inner_cycle> inner_cycles;
    std::list<Point_3>     isolated_vertices;

    typedef typename std::list<Outer_cycle>::iterator Outer_cycle_iterator;
    typedef typename std::list<Inner_cycle>::iterator Inner_cycle_iterator;
    typedef typename std::list<Point_3>::iterator Isolated_vertex_iterator;

    Partial_facet() {}

    Partial_facet(const Partial_facet& pf) {
      f = pf.f;
      outer_cycles = pf.outer_cycles;
      inner_cycles = pf.inner_cycles;
      isolated_vertices = pf.isolated_vertices;
    }

    Partial_facet& operator=(const Partial_facet& pf) {
      f = pf.f;
      outer_cycles = pf.outer_cycles;
      inner_cycles = pf.inner_cycles;
      isolated_vertices = pf.isolated_vertices;
      return *this;
    }

    explicit Partial_facet(Halffacet_handle fin) : f(fin) {
      Halffacet_cycle_iterator fc = f->facet_cycles_begin();
      for(;fc != f->facet_cycles_end();++fc) {
	if(fc.is_shalfedge()) {
	  SHalfedge_around_facet_circulator se(fc), se_next(se);
	  ++se_next;
	  if(fc == f->facet_cycles_begin()) {
	    outer_cycles.push_back(Outer_cycle(se, se));
	    //	    outer_cycles.push_back(Outer_cycle(se_next, se));
	  } else {
	    inner_cycles.push_back(Inner_cycle(se, se));
	    //	    inner_cycles.push_back(Inner_cycle(se_next, se));
	  }
	} else if(fc.is_shalfloop()) {
	  SHalfloop_handle l(fc);
	  isolated_vertices.push_back(l->incident_sface()->center_vertex()->point());
	} else
	  CGAL_error_msg( "wrong value");
      }
    }

    Outer_cycle_iterator outer_cycles_begin() { return outer_cycles.begin(); }
    Inner_cycle_iterator inner_cycles_begin() { return inner_cycles.begin(); }
    Isolated_vertex_iterator isolated_vertices_begin() { return isolated_vertices.begin(); }

    Outer_cycle_iterator outer_cycles_end() { return outer_cycles.end(); }
    Inner_cycle_iterator inner_cycles_end() { return inner_cycles.end(); }
    Isolated_vertex_iterator isolated_vertices_end() { return isolated_vertices.end(); }

    bool divide(const Plane_3& p, Partial_facet& pf1, Partial_facet& pf2) {
      //      std::cerr << "divide " << std::endl;
      //      debug();
      pf1.f = pf2.f = f;
      Outer_cycle_iterator oc = outer_cycles.begin();
      for(;oc != outer_cycles.end(); ++oc) {
	bool next = false;
	//	CGAL_assertion(oc->first != oc->second);
	SHalfedge_around_facet_circulator se = oc->first, se_begin(se), se_new(se), se_end;
	Oriented_side ref = p.oriented_side(se->source()->source()->point()), cur;
	//	std::cerr << "start " << se->source()->source()->point() << ":" << ref << std::endl;
	++se;
	while(ref == ON_ORIENTED_BOUNDARY && se != oc->second) {
	  ref = p.oriented_side(se->source()->source()->point());
	  ++se;
	}
	if(se == oc->second) 
	  return false;

	for(;se != oc->second;++se) {
	  cur = p.oriented_side(se->source()->source()->point());
	  //	  std::cerr << "current " << se->source()->source()->point() << ":" << cur  << std::endl;
	  if(cur != ref) {
	    CGAL_assertion(ref != ON_ORIENTED_BOUNDARY);
	    if(cur == ON_ORIENTED_BOUNDARY) {
	      next = true;
	      continue;
	    }
	    se_end = se;
	    if(next)
	      --se_end;
	    if(cur == ON_NEGATIVE_SIDE) {
	      pf2.outer_cycles.push_back(Outer_cycle(se_begin, se_end));
	    } else if(cur == ON_POSITIVE_SIDE) {
	      pf1.outer_cycles.push_back(Outer_cycle(se_begin, se_end));
	    }
	    se_begin = se_new;
	    if(next)
	      ++se_begin;
	    ref = cur;
	  } else
	      se_new = se;
	  next = false;
	}
	//	std::cerr << "end of cycle " << ref << std::endl;
	if(next) {
	  if(ref == ON_POSITIVE_SIDE) {
	    pf1.outer_cycles.push_back(Outer_cycle(se_new, se));
	    pf2.outer_cycles.push_back(Outer_cycle(se_begin, --se));
	  } else {
	    pf2.outer_cycles.push_back(Outer_cycle(se_new, se));	 
	    pf1.outer_cycles.push_back(Outer_cycle(se_begin, --se));	 
	  }
	} else {
	  if(ref == ON_POSITIVE_SIDE)
	    pf2.outer_cycles.push_back(Outer_cycle(se_begin, se));
	  else if(ref == ON_NEGATIVE_SIDE)
	    pf1.outer_cycles.push_back(Outer_cycle(se_begin, se));
	}
      }

      Inner_cycle_iterator ic = inner_cycles.begin();
      for(;ic != inner_cycles.end(); ++ic) {
	bool next = false;
	SHalfedge_around_facet_circulator se = ic->first, se_begin(se), se_new(se), se_end;
	Oriented_side ref = p.oriented_side(se->source()->source()->point()), cur;
	++se;
	while(ref == ON_ORIENTED_BOUNDARY && se != ic->second) {
	  ref = p.oriented_side(se->source()->source()->point());
	  ++se;
	}
	if(se == ic->second) 
	  return false;

	for(;se != ic->second; ++se) {
	  cur = p.oriented_side(se->source()->source()->point());
	  if(cur != ref) {
	    CGAL_assertion(ref != ON_ORIENTED_BOUNDARY);
	    if(cur == ON_ORIENTED_BOUNDARY) {
	      next = true;
	      continue;
	    }
	    se_end = se;
	    if(next)
	      --se_end;
	    if(cur == ON_NEGATIVE_SIDE) {
	      pf2.inner_cycles.push_back(Inner_cycle(se_begin, se_end));
	    } else if(cur == ON_POSITIVE_SIDE) {
	      pf1.inner_cycles.push_back(Inner_cycle(se_begin, se_end));
	    }
	    se_begin = se_new;
	    if(next)
	      ++se_begin;
	    ref = cur;
	  } else
	      se_new = se;
	  next = false;
	}
	if(next) {
	  if(ref == ON_POSITIVE_SIDE) {
	    pf1.inner_cycles.push_back(Inner_cycle(se_new, se));
	    pf2.inner_cycles.push_back(Inner_cycle(se_begin, --se));
	  } else {
	    pf2.inner_cycles.push_back(Inner_cycle(se_new, se));	 
	    pf1.inner_cycles.push_back(Inner_cycle(se_begin, --se));	 
	  }
	} else {
	  if(ref == ON_POSITIVE_SIDE)
	    pf2.inner_cycles.push_back(Inner_cycle(se_begin, se));
	  else if(ref == ON_NEGATIVE_SIDE)
	    pf1.inner_cycles.push_back(Inner_cycle(se_begin, se));
	}
      }

      Isolated_vertex_iterator iv = isolated_vertices.begin();
      for(;iv != isolated_vertices.end();++iv) {
	Oriented_side side = p.oriented_side(*iv);
	if( side == ON_NEGATIVE_SIDE || side == ON_ORIENTED_BOUNDARY)
	  pf1.isolated_vertices.push_back(*iv);
	if( side == ON_POSITIVE_SIDE || side == ON_ORIENTED_BOUNDARY)
	  pf1.isolated_vertices.push_back(*iv);
      }
      //      std::cerr << "into " << std::endl;
      //      pf1.debug();      
      //      pf2.debug();
      return true;
    }

    void debug() {
      std::cerr << "Partial_facet " << std::endl;
      std::cerr << "Box " << std::endl;
      std::cerr << " " << f->b.min_coord(0) << std::endl;
      std::cerr << " " << f->b.min_coord(1) << std::endl;
      std::cerr << " " << f->b.min_coord(2) << std::endl;
      std::cerr << " " << f->b.max_coord(0) << std::endl;
      std::cerr << " " << f->b.max_coord(1) << std::endl;
      std::cerr << " " << f->b.max_coord(2) << std::endl;

      Outer_cycle_iterator oc = outer_cycles_begin();
      for(; oc != outer_cycles_end(); ++oc) {
	std::cerr << "Outer cycle " << std::endl;
	SHalfedge_around_facet_circulator sb(oc->first), se(oc->second);
	CGAL_For_all(sb,se) {
	  std::cerr << "  " << sb->source()->source()->point() << std::endl;
	}
      }

      Inner_cycle_iterator ic = inner_cycles_begin();
      for(; ic != inner_cycles_end(); ++ic) {
	std::cerr << "Inner cycle " << std::endl;
	SHalfedge_around_facet_circulator sb(ic->first), se(ic->second);
	CGAL_For_all(sb,se) {
	  std::cerr << "  " << sb->source()->source()->point() << std::endl;
	}
      }
    }  
  };

#endif

  /*{\Mcreation 3}*/
  /*{\Mtext |\Mname| is default and copy constructible. Note that copy
  construction means cloning an isomorphic structure and is thus an
  expensive operation.}*/

  SNC_structure() : 
    boundary_item_(boost::none), sm_boundary_item_(boost::none),
    vertices_(), halfedges_(), halffacets_(), volumes_(),
    shalfedges_(), shalfloops_(), sfaces_() {}
  ~SNC_structure() { CGAL_NEF_TRACEN("~SNC_structure: clearing "<<this); clear(); }

  SNC_structure(const Self& D) : 
    boundary_item_(boost::none), sm_boundary_item_(boost::none),
    vertices_(D.vertices_), halfedges_(D.halfedges_), 
    halffacets_(D.halffacets_), volumes_(D.volumes_),
    shalfedges_(D.shalfedges_), shalfloops_(D.shalfloops_), 
    sfaces_(D.sfaces_)
  { pointer_update(D); }

  Self& operator=(const Self& D) { 
    if ( this == &D ) 
      return *this;
    clear();
    boundary_item_.clear(boost::none);
    sm_boundary_item_.clear(boost::none);
    vertices_ = D.vertices_;
    halfedges_ = D.halfedges_;
    halffacets_ = D.halffacets_;
    volumes_ = D.volumes_;
    shalfedges_ = D.shalfedges_;
    shalfloops_ = D.shalfloops_;
    sfaces_ = D.sfaces_;
    pointer_update(D);
    return *this;
  }

  void clear_boundary() {
    boundary_item_.clear(boost::none);
    sm_boundary_item_.clear(boost::none);
  }

  void clear_snc_boundary() {
    boundary_item_.clear(boost::none);
  }

  void clear() { 
    boundary_item_.clear(boost::none);
    sm_boundary_item_.clear(boost::none);
    vertices_.destroy();
    halfedges_.destroy();
    halffacets_.destroy();
    volumes_.destroy();
    shalfedges_.destroy();
    shalfloops_.destroy();
    sfaces_.destroy();
  }

  template <typename H>
  bool is_boundary_object(H h) 
  { return boundary_item_[h]!=boost::none; }
  template <typename H>
  bool is_sm_boundary_object(H h) 
  { return sm_boundary_item_[h]!=boost::none; }

  template <typename H>
  Object_iterator& boundary_item(H h)
  { return *boundary_item_[h]; }
  template <typename H>
  Object_iterator& sm_boundary_item(H h)
  { return *sm_boundary_item_[h]; }

  template <typename H>
  void store_boundary_item(H h, Object_iterator o)
  { boundary_item_[h] = o; }
  template <typename H>
  void store_sm_boundary_item(H h, Object_iterator o)
  { sm_boundary_item_[h] = o; }

  template <typename H>
  void undef_boundary_item(H h)
  { CGAL_assertion(boundary_item_[h]!=boost::none);
    boundary_item_[h] = boost::none; }
  template <typename H>
  void undef_sm_boundary_item(H h)
  { CGAL_assertion(sm_boundary_item_[h]!=boost::none);
    sm_boundary_item_[h] = boost::none; }

  void reset_iterator_hash(Object_iterator it)
  { SVertex_handle sv;
    SHalfedge_handle se;
    SHalfloop_handle sl;
    if ( CGAL::assign(se,*it) ) { 
      if( is_boundary_object(se)) 
	undef_boundary_item(se); 
      return; 
    }
    if ( CGAL::assign(sl,*it) ) { 
      if( is_boundary_object(sl)) 
	undef_boundary_item(sl);
      return; 
    }
    if ( CGAL::assign(sv,*it) ) { 
      if( is_boundary_object(sv)) 
	undef_boundary_item(sv); 
      return; 
    }  
  }

  void reset_sm_object_list(Object_list& L)
  { Object_iterator oit;
    CGAL_forall_iterators(oit,L) reset_sm_iterator_hash(oit);
    L.clear();
  }

  void reset_sm_iterator_hash(Object_iterator it)
  { SVertex_handle sv;
    SHalfedge_handle se;
    SHalfloop_handle sl;
    if ( CGAL::assign(se,*it) ) { 
      if( is_sm_boundary_object(se)) 
	undef_sm_boundary_item(se); 
      return; 
    }
    if ( CGAL::assign(sl,*it) ) { 
      if( is_sm_boundary_object(sl)) 
	undef_sm_boundary_item(sl);
      return; 
    }
    if ( CGAL::assign(sv,*it) ) { 
      if( is_sm_boundary_object(sv)) 
	undef_sm_boundary_item(sv); 
      return; 
    }
  }

  void reset_object_list(Object_list& L)
  { Object_iterator oit;
    CGAL_forall_iterators(oit,L) reset_iterator_hash(oit);
    L.clear();
  }

  /*{\Moperations 2.5 3}*/

  // The constant iterators and circulators.
  Vertex_const_iterator vertices_begin() const 
    { return vertices_.begin();}
  Vertex_const_iterator vertices_end() const 
    { return vertices_.end();}
  Halfedge_const_iterator halfedges_begin() const 
    { return halfedges_.begin();}
  Halfedge_const_iterator halfedges_end() const 
    { return halfedges_.end();}
  Halffacet_const_iterator halffacets_begin() const 
    { return halffacets_.begin();}
  Halffacet_const_iterator halffacets_end() const 
    { return halffacets_.end();}
  Volume_const_iterator   volumes_begin() const 
    { return volumes_.begin();}
  Volume_const_iterator   volumes_end() const 
    { return volumes_.end();}

  SVertex_const_iterator svertices_begin() const 
    { return halfedges_.begin();}
  SVertex_const_iterator svertices_end() const 
    { return halfedges_.end();}
  SHalfedge_const_iterator shalfedges_begin() const 
    { return shalfedges_.begin();}
  SHalfedge_const_iterator shalfedges_end() const 
    { return shalfedges_.end();}
  SHalfloop_const_iterator shalfloops_begin() const 
    { return shalfloops_.begin();}
  SHalfloop_const_iterator shalfloops_end() const 
    { return shalfloops_.end();}
  SFace_const_iterator sfaces_begin() const 
    { return sfaces_.begin();}
  SFace_const_iterator sfaces_end() const 
    { return sfaces_.end();}

  Vertex_iterator    vertices_begin()   { return vertices_.begin();}
  Vertex_iterator    vertices_end()     { return vertices_.end();}
  Halfedge_iterator  halfedges_begin()  { return halfedges_.begin();}
  Halfedge_iterator  halfedges_end()    { return halfedges_.end();}
  Halffacet_iterator halffacets_begin() { return halffacets_.begin();}
  Halffacet_iterator halffacets_end()   { return halffacets_.end();}
  Volume_iterator    volumes_begin()    { return volumes_.begin();}
  Volume_iterator    volumes_end()      { return volumes_.end();}

  SVertex_iterator   svertices_begin()  { return halfedges_.begin();}
  SVertex_iterator   svertices_end()    { return halfedges_.end();}
  SHalfedge_iterator shalfedges_begin() { return shalfedges_.begin();}
  SHalfedge_iterator shalfedges_end()   { return shalfedges_.end();}
  SHalfloop_iterator shalfloops_begin() { return shalfloops_.begin();}
  SHalfloop_iterator shalfloops_end()   { return shalfloops_.end();}
  SFace_iterator     sfaces_begin()     { return sfaces_.begin();}
  SFace_iterator     sfaces_end()       { return sfaces_.end();}

  /*{\Mtext The list of all objects can be accessed via iterator ranges.
  For comfortable iteration we also provide iterations macros. 
  The iterator range access operations are of the following kind:\\
  |Vertex_iterator vertices_begin()/vertices_end()|\\
  |Halfedge_iterator halfedges_begin()/halfedges_end()|\\
  |Halffacet_iterator halffacets_begin()/halffacets_end()|\\
  |Volume_iterator volumes_begin()/volumes_end()|

  The macros are then |CGAL_forall_vertices(v,\Mvar)|, 
  |CGAL_forall_halfedges(e,\Mvar)|, |CGAL_forall_edges(e,\Mvar)|,
  |CGAL_forall_halffacets(f,\Mvar)|, |CGAL_forall_facets(f,\Mvar)|,
  |CGAL_forall_volumes(w,\Mvar)|.}*/

  Size_type number_of_vertices() const  { return vertices_.size(); }
  /*{\Mop returns the number of vertices.}*/
  Size_type number_of_halfedges() const { return halfedges_.size(); }
  /*{\Mop returns the number of (directed edges).}*/
  Size_type number_of_edges() const { return halfedges_.size()/2; }
  /*{\Mop returns the number of (directed edges).}*/
  Size_type number_of_halffacets() const { return halffacets_.size();}
  /*{\Mop returns the number of halffacets.}*/
  Size_type number_of_facets() const { return halffacets_.size()/2;}
  /*{\Mop returns the number of facets.}*/
  Size_type number_of_volumes() const   { return volumes_.size();}
  /*{\Mop returns the number of volumes.}*/
  Size_type number_of_svertices() const   { return halfedges_.size();}
  /*{\Mop returns the number of svertices.}*/
  Size_type number_of_shalfedges() const { return shalfedges_.size();}
  /*{\Mop returns the number of sedges.}*/
  Size_type number_of_shalfloops() const { return shalfloops_.size();}
  /*{\Mop returns the number of sloops.}*/
  Size_type number_of_sfaces() const { return sfaces_.size();}
  /*{\Mop returns the number of sfaces.}*/

  void print_statistics(std::ostream& os = std::cout) const
  /*{\Mop print the statistics of |P|: the number of vertices, edges, 
      faces, volumes and sobjects.}*/
  {
    os << "Selective Nef Complex - Statistics\n";
    os << "|V| = " << number_of_vertices() << std::endl;
    os << "|E| = " << number_of_halfedges() << std::endl;
    os << "|F| = " << number_of_halffacets() << std::endl;
    os << "|C| = " << number_of_volumes() << std::endl;
    os << "|VS| = " << number_of_svertices() << std::endl;
    os << "|ES| = " << number_of_shalfedges() << std::endl;
    os << "|LS| = " << number_of_shalfloops() << std::endl;
    os << "|FS| = " << number_of_sfaces() << std::endl;
    os << std::endl;
  }

  bool is_empty() const {
  /*{\Mop returns true if |\Mvar| is empty, false otherwise.}*/
    return number_of_vertices() == 0 &&
           number_of_halfedges() == 0 &&
           number_of_halffacets() == 0 &&
           number_of_volumes() == 0 &&
           number_of_shalfedges() == 0 &&
           number_of_shalfloops() == 0 &&
           number_of_sfaces() == 0;
  }

  bool has_bbox_only() const  {
  /*{\Mop returns true if |\Mvar| is only the infimaximal box, 
    false otherwise.}*/
    return (number_of_vertices() == 8 &&
	    number_of_edges() == 12 &&
	    number_of_facets() == 6 &&
	    number_of_volumes() == 2 &&
	    (++volumes_begin())->mark() == false);
  }

  Vertex_handle new_vertex(const Point_3& p = Point_3(), Mark m = Mark())
  /*{\Mop returns a new vertex at point |p| marked by |m|.}*/ { 
    Vertex_handle vh = new_vertex_only();
    vh->point() = p;
    vh->mark() = m;
    vh->sncp() = this;
    vh->svertices_begin() = vh->svertices_last() = svertices_end();
    vh->shalfedges_begin() = vh->shalfedges_last() = shalfedges_end();
    vh->sfaces_begin() = vh->sfaces_last() = sfaces_end();
    vh->shalfloop() = shalfloops_end();
    return vh;
  }

  Halfedge_handle new_halfedge_pair(Vertex_handle v1, Vertex_handle v2,
				    Mark m = Mark())
  /*{\Mop creates a new halfedge pair between the vertices $v_1$
  and $v_2$. The edge is marked by |m|.}*/ { 
    SM_decorator D1(&*v1);
    SM_decorator D2(&*v2);
    SVertex_handle e1 = D1.new_vertex();
    SVertex_handle e2 = D2.new_vertex();
    make_twins(e1,e2);
    e1->mark() = m;
    return e1;
  }

  Halffacet_handle new_halffacet_pair(const Plane_3& h = Plane_3(), 
				      Mark m = Mark())
  /*{\Mop creates a new facet supported by the plane |h| and
  marked with |m|.}*/ {
    Halffacet_handle f1 = new_halffacet_only();
    Halffacet_handle f2 = new_halffacet_only();
    f1->plane() = h; f2->plane() = h.opposite();
    make_twins(f1,f2);
    f1->mark() = f2->mark() = m;
    return f1;
  }

  Volume_handle new_volume(Mark m = Mark())
  /*{\Mop creates a new volume marked with |m|.}*/ { 
    Volume_handle vh = new_volume_only();
    vh->mark() = m;
    return vh;
  }

  template <typename H>
  void make_twins(H h1, H h2) { 
    h1->twin() = h2; h2->twin() = h1; 
  }

  void delete_vertex(Vertex_handle v)
  /*{\Mop deletes the vertex including the objects in its local graph.}*/  { 
    CGAL_NEF_TRACEN("~ deleting vertex "<<&*v<<" from "<<&*this);
    v->clear(true); 
    delete_vertex_only(v);
    CGAL_NEF_TRACEN("~~ vertex deleted"<<&*v);
  }

  void delete_halfedge_pair(Halfedge_handle e)
  /*{\Mop deletes the halfedge pair of |e,twin(e)|.  Does not care about
  incident objects in the local graph of |source(e)|.}*/ { 
    CGAL_NEF_TRACEN("~ deleting halfedges pair "<<&*e<<", "<<&*(e->twin())<<
	   " from "<<&*this);
    Halfedge_handle et = e->twin();
    SM_decorator D1(&*e->center_vertex()), D2(&*et->center_vertex());
    D1.delete_vertex(e);
    D2.delete_vertex(et);
  }

  void delete_halffacet_pair(Halffacet_handle f)
  /*{\Mop deletes the halffacet pair |f,twin(f)|. Does not care about
  boundary cycle objects.}*/ { 
    CGAL_NEF_TRACEN("~ deleting halffacets pair "<<&*f<<", "<<&*(f->twin())<<
	   " from "<<&*this);
    reset_object_list(f->boundary_entry_objects());
    reset_object_list(f->twin()->boundary_entry_objects());
    delete_halffacet_only(f->twin());
    delete_halffacet_only(f);
  }

  void delete_volume(Volume_handle c)
  /*{\Mop deletes the volume |c|. Does not care about shell objects.}*/ { 
    CGAL_NEF_TRACEN("~ deleting volume "<<&*c<<" from "<<&*this);
    reset_object_list(c->shell_entry_objects());
    delete_volume_only(c);
  }

  Vertex_alloc vertex_allocator;
  Vertex* get_vertex_node( const Vertex& ) {
    Vertex* p = vertex_allocator.allocate(1);
    vertex_allocator.construct( p, Vertex());
    return p;
  }
  void put_vertex_node( Vertex* p) {
    vertex_allocator.destroy(p);
    vertex_allocator.deallocate( p, 1);
  }

  Halfedge_alloc halfedge_allocator;
  Halfedge* get_halfedge_node( const Halfedge&) {
    Halfedge* p = halfedge_allocator.allocate(1);
    halfedge_allocator.construct( p, Halfedge());
    return p;
  }
  void put_halfedge_node( Halfedge* p) {
    halfedge_allocator.destroy(p);
    halfedge_allocator.deallocate( p, 1);
  }

  Halffacet_alloc halffacet_allocator;
  Halffacet* get_halffacet_node( const Halffacet& ) {
    Halffacet* p = halffacet_allocator.allocate(1);
    halffacet_allocator.construct( p, Halffacet());
    return p;
  }
  void put_halffacet_node( Halffacet* p) {
    halffacet_allocator.destroy(p);
    halffacet_allocator.deallocate( p, 1);
  }

  Volume_alloc volume_allocator;
  Volume* get_volume_node( const Volume& ) {
    Volume* p = volume_allocator.allocate(1);
    volume_allocator.construct( p, Volume());
    return p;
  }
  void put_volume_node( Volume* p) {
    volume_allocator.destroy(p);
    volume_allocator.deallocate( p, 1);
  }

  SHalfedge_alloc shalfedge_allocator;
  SHalfedge* get_shalfedge_node( const SHalfedge& ) {
    SHalfedge* p = shalfedge_allocator.allocate(1);
    shalfedge_allocator.construct( p, SHalfedge());
    return p;
  }
  void put_shalfedge_node( SHalfedge* p) {
    shalfedge_allocator.destroy(p);
    shalfedge_allocator.deallocate( p, 1);
  }

  SHalfloop_alloc shalfloop_allocator;
  SHalfloop* get_shalfloop_node( const SHalfloop& ) {
    SHalfloop* p = shalfloop_allocator.allocate(1);
    shalfloop_allocator.construct( p, SHalfloop());
    return p;
  }
  void put_shalfloop_node( SHalfloop* p) {
    shalfloop_allocator.destroy(p);
    shalfloop_allocator.deallocate( p, 1);
  }

  SFace_alloc sface_allocator;
  SFace* get_sface_node( const SFace& ) {
    SFace* p = sface_allocator.allocate(1);
    sface_allocator.construct( p, SFace());
    return p;
  }
  void put_sface_node( SFace* p) {
    sface_allocator.destroy(p);
    sface_allocator.deallocate( p, 1);
  }


  Vertex_handle new_vertex_only() { 
    vertices_.push_back(* get_vertex_node(Vertex()));
    CGAL_NEF_TRACEN("  new vertex only "<<&*(--vertices_end()));
    return --vertices_end(); 
  }
  Halfedge_handle new_halfedge_only(Halfedge_handle e)  { 
    Halfedge_handle ne = halfedges_.insert(e, * get_halfedge_node(Halfedge()));
    CGAL_NEF_TRACEN("  after "<<&*e<<" new halfedge only "<<&*ne);
    return ne;
  }
  Halfedge_handle new_halfedge_only()  { 
    CGAL_NEF_TRACEN("  new halfedge only "<<&*(--halfedges_end()));
    halfedges_.push_back( * get_halfedge_node(Halfedge()));
    return --halfedges_end();
  }
  Halffacet_handle new_halffacet_only()  { 
    halffacets_.push_back( * get_halffacet_node(Halffacet()));
    CGAL_NEF_TRACEN("  new halffacet only "<<&*(--halffacets_end()));
    return --halffacets_end(); 
  } 
  Volume_handle new_volume_only()  { 
    volumes_.push_back( * get_volume_node(Volume()));
    CGAL_NEF_TRACEN("  new volume only "<<&*(--volumes_end()));
    return --volumes_end(); 
  }
  SHalfedge_handle new_shalfedge_only()  {
    shalfedges_.push_back( * get_shalfedge_node(SHalfedge()));
    CGAL_NEF_TRACEN("  new shalfedge only "<<&*(--shalfedges_end()));
    return --shalfedges_end();
  }
  SHalfedge_handle new_shalfedge_only(SHalfedge_handle se) {
    SHalfedge_handle nse = shalfedges_.insert(se, * get_shalfedge_node(SHalfedge()));
    CGAL_NEF_TRACEN("  after " << &*se << " new shalfedge only " << &*nse);
    return nse;
  }
  SHalfloop_handle new_shalfloop_only()  {
    shalfloops_.push_back( * get_shalfloop_node(SHalfloop()));
    CGAL_NEF_TRACEN("  new shalfloop only "<<&*(--shalfloops_end()));
    return --shalfloops_end(); 
  }
  SFace_handle new_sface_only() {
    sfaces_.push_back( * get_sface_node(SFace()));
    CGAL_NEF_TRACEN("  new sface only "<<&*(--sfaces_end()));
    return --sfaces_end(); 
  }
  SFace_handle new_sface_only(SFace_handle sf) {
    SFace_handle nsf = sfaces_.insert(sf, * get_sface_node(SFace()));
    CGAL_NEF_TRACEN("  after " << &*sf << " new sface only " << &*nsf);
    return nsf;
  }

  void delete_vertex_only(Vertex_handle h) {
    CGAL_NEF_TRACEN("~ deleting vertex only "<<&*h<<" from "<<&*this);
    vertices_.erase(h);
    put_vertex_node(&*h);
  }
  void delete_halfedge_only(Halfedge_handle h) { 
    CGAL_NEF_TRACEN("~ deleting halfedge only "<<&*h<<" from "<<&*this);
    CGAL_assertion(!is_sm_boundary_object(h));
    halfedges_.erase(h);
    put_halfedge_node(&*h);
  }
  void delete_halffacet_only(Halffacet_handle h) { 
    CGAL_NEF_TRACEN("~ deleting halffacet only "<<&*h<<" from "<<&*this);
    halffacets_.erase(h);         
    put_halffacet_node(&*h);
  }
  void delete_volume_only(Volume_handle h) {
    CGAL_NEF_TRACEN("~ deleting volume only "<<&*h<<" from "<<&*this);
    volumes_.erase(h); 
    put_volume_node(&*h);
  }
  void delete_shalfedge_only(SHalfedge_handle h)  { 
    CGAL_NEF_TRACEN("~ deleting shalfedge only "<<&*h<<" from "<<&*this);
    CGAL_assertion(!is_sm_boundary_object(h));
    shalfedges_.erase(h);
    put_shalfedge_node(&*h);
  }
  void delete_shalfloop_only(SHalfloop_handle h)  { 
    CGAL_NEF_TRACEN("~ deleting shalfloop only "<<&*h<<" from "<<&*this);
    CGAL_assertion(!is_sm_boundary_object(h));
    shalfloops_.erase(h); 
    put_shalfloop_node(&*h);
  }
  void delete_sface_only(SFace_handle h)  { 
    CGAL_NEF_TRACEN("~ deleting sface only "<<&*h<<" from "<<&*this);
    CGAL_assertion(!is_boundary_object(h));
    sfaces_.erase(h);
    put_sface_node(&*h);
  }
    //SL: in the following function, I guess the sizeof(void*) is related to the void* info that was
    //used together with geninfo to store an arbitrary type. I replaced that with any and did not changed that
    std::size_t bytes() {
      // bytes used for the SNC_structure

      std::size_t result = sizeof(Self);
      result += number_of_vertices()   * (sizeof(Vertex) 
					  - sizeof(Point_3));
      result += number_of_halfedges()  * (sizeof(Halfedge) 
					  - sizeof(Sphere_point));
      result += number_of_halffacets() * (sizeof(Halffacet)
					  - sizeof(Plane_3));
      result += number_of_volumes()    * sizeof(Volume);
      result += number_of_shalfedges() * (sizeof(SHalfedge)
					  - sizeof(Sphere_circle));
      result += number_of_shalfloops() * sizeof(SHalfloop);
      result += number_of_sfaces()     * sizeof(SFace);

      Halffacet_iterator hf;
      CGAL_forall_halffacets(hf, *this) {
	Halffacet_cycle_iterator fc;
	CGAL_forall_facet_cycles_of(fc, hf)
	  result += sizeof(*fc) + 2 * sizeof(void*);
      }

      Volume_iterator c;
      CGAL_forall_volumes(c, *this) {
	Shell_entry_iterator sei;
	CGAL_forall_shells_of(sei,c) 
	  result += sizeof(*sei) + 2 * sizeof(void*);
      }

      SFace_iterator sf;
      CGAL_forall_sfaces(sf, *this) {
	SFace_cycle_iterator sfc;
	CGAL_forall_sface_cycles_of(sfc,sf)
	  result += sizeof(*sfc) + 2 * sizeof(void*);
      }

      return result;
    }

    //SL: in the following function, I guess the sizeof(void*) is related to the void* info that was
    //used together with geninfo to store an arbitrary type. I replaced that with any and did not changed that
    std::size_t bytes_reduced() {
      // bytes used for the SNC_structure

      std::size_t result = sizeof(Self);
      result += number_of_vertices()   * (sizeof(Vertex) 
					  - sizeof(Point_3)
					  - sizeof(SHalfloop_iterator)
					  - 2 * sizeof(Mark)
					  - sizeof(void*));
      result += number_of_halfedges()  * (sizeof(Halfedge) 
					  - sizeof(SHalfedge_handle)
					  - sizeof(SFace_handle)
					  - sizeof(void*)
					  - sizeof(Sphere_point));
      result += number_of_halffacets() * (sizeof(Halffacet)
					  - sizeof(Plane_3));
      result += number_of_volumes()    * sizeof(Volume);
      result += number_of_shalfedges() * (sizeof(SHalfedge)
					  - sizeof(SVertex_handle)
					  - 3 * sizeof(SHalfedge_handle)
					  - sizeof(SFace_handle)
					  - sizeof(void*)
					  - sizeof(Mark)
					  - sizeof(Sphere_circle));
      result += number_of_sfaces()     * (sizeof(SFace)
					  - sizeof(void*)
					  - sizeof(Mark)
					  - sizeof(Object_list));

      Halffacet_iterator hf;
      CGAL_forall_halffacets(hf, *this) {
	Halffacet_cycle_iterator fc;
	CGAL_forall_facet_cycles_of(fc, hf)
	  result += sizeof(*fc) + 2 * sizeof(void*);
      }

      Volume_iterator c;
      CGAL_forall_volumes(c, *this) {
	Shell_entry_iterator sei;
	CGAL_forall_shells_of(sei,c) 
	  result += sizeof(*sei) + 2 * sizeof(void*);
      }

      return result;
    }

    std::size_t bytes_reduced2() {
      // bytes used for the SNC_structure

      std::size_t result = sizeof(Self);
      result += number_of_vertices()   * (sizeof(Mark) 
					  + sizeof(SNC_structure*)
					  + sizeof(Object_list)
					  + 2 * sizeof(SFace_handle));
      result += number_of_halfedges()  * (sizeof(Vertex_handle)
					  + sizeof(SVertex_handle)
					  + sizeof(Mark)
					  + 2 * sizeof(Object_handle));
      result += number_of_halffacets() * (sizeof(Halffacet)
					  - sizeof(Plane_3));
      result += number_of_volumes()    * sizeof(Volume);
      result += number_of_shalfedges() * (2 * sizeof(SHalfedge_handle)
					  + sizeof(Halffacet_handle));
      result += number_of_shalfloops() * sizeof(SHalfloop);
      result += number_of_sfaces()     * (sizeof(Vertex_handle)
					  + sizeof(Volume_handle));

      Halffacet_iterator hf;
      CGAL_forall_halffacets(hf, *this) {
	Halffacet_cycle_iterator fc;
	CGAL_forall_facet_cycles_of(fc, hf)
	  result += sizeof(*fc) + 2 * sizeof(void*);
      }

      Volume_iterator c;
      CGAL_forall_volumes(c, *this) {
	Shell_entry_iterator sei;
	CGAL_forall_shells_of(sei,c) 
	  result += sizeof(*sei) + 2 * sizeof(void*);
      }

      return result;
    }

protected:
  void pointer_update(const Self& D);
  
  typedef boost::optional<Object_iterator> Optional_object_iterator ;
 private:
  Generic_handle_map<Optional_object_iterator> boundary_item_;
  Generic_handle_map<Optional_object_iterator> sm_boundary_item_;
 protected:
  Vertex_list    vertices_;
  Halfedge_list  halfedges_;
  Halffacet_list halffacets_;
  Volume_list    volumes_;
  SHalfedge_list shalfedges_;
  SHalfloop_list shalfloops_;
  SFace_list     sfaces_;

}; // SNC_structure


template <typename Kernel, typename Items, typename Mark>
void SNC_structure<Kernel,Items,Mark>::
pointer_update(const SNC_structure<Kernel,Items,Mark>& D)
{
  CGAL::Unique_hash_map<Vertex_const_handle,Vertex_handle>       VM;
  CGAL::Unique_hash_map<Halfedge_const_handle,Halfedge_handle>   EM;
  CGAL::Unique_hash_map<Halffacet_const_handle,Halffacet_handle> FM;
  CGAL::Unique_hash_map<Volume_const_handle,Volume_handle>       CM;
  CGAL::Unique_hash_map<SHalfedge_const_handle,SHalfedge_handle> SEM;
  CGAL::Unique_hash_map<SHalfloop_const_handle,SHalfloop_handle> SLM;
  CGAL::Unique_hash_map<SFace_const_handle,SFace_handle>         SFM;
  Vertex_const_iterator vc = D.vertices_begin();
  Vertex_iterator v = vertices_begin();
  for ( ; vc != D.vertices_end(); ++vc,++v) VM[vc] = v;
  VM[D.vertices_end()] = vertices_end();
  Halfedge_const_iterator ec = D.halfedges_begin();
  Halfedge_iterator e = halfedges_begin();
  for ( ; ec != D.halfedges_end(); ++ec,++e) EM[ec] = e;
  EM[D.halfedges_end()] = halfedges_end();
  Halffacet_const_iterator fc = D.halffacets_begin();
  Halffacet_iterator f = halffacets_begin();
  for ( ; fc != D.halffacets_end(); ++fc,++f) FM[fc] = f;
  FM[D.halffacets_end()] = halffacets_end();
  Volume_const_iterator cc = D.volumes_begin();
  Volume_iterator c = volumes_begin();
  for ( ; cc != D.volumes_end(); ++cc,++c) CM[cc] = c;
  CM[D.volumes_end()] = volumes_end();
  SHalfedge_const_iterator sec = D.shalfedges_begin();
  SHalfedge_iterator se = shalfedges_begin();
  for ( ; sec != D.shalfedges_end(); ++sec,++se) SEM[sec] = se;
  SEM[D.shalfedges_end()] = shalfedges_end();
  SHalfloop_const_iterator slc = D.shalfloops_begin();
  SHalfloop_iterator sl = shalfloops_begin();
  for ( ; slc != D.shalfloops_end(); ++slc,++sl) SLM[slc] = sl;
  SLM[D.shalfloops_end()] = shalfloops_end();
  SFace_const_iterator sfc = D.sfaces_begin();
  SFace_iterator sf = sfaces_begin();
  for ( ; sfc != D.sfaces_end(); ++sfc,++sf) SFM[sfc] = sf;
  SFM[D.sfaces_end()] = sfaces_end();

  CGAL_forall_vertices(v,*this) {
    // Local Graph update: (SVertices are postponed/updated as Edges)
    v->sncp() = this;
    v->svertices_begin() = EM[v->svertices_begin()];
    v->svertices_last() = EM[v->svertices_last()];
    v->shalfedges_begin() = SEM[v->shalfedges_begin()];
    v->shalfedges_last() = SEM[v->shalfedges_last()];
    v->sfaces_begin() = SFM[v->sfaces_begin()];
    v->sfaces_last() = SFM[v->sfaces_last()];
    v->shalfloop() = SLM[v->shalfloop()];
  }
  // Halfedge update:
  CGAL_forall_halfedges(e,*this) {
    e->center_vertex() = VM[e->center_vertex()];
    e->twin() = EM[e->twin()];
    e->out_sedge() = SEM[e->out_sedge()];
    e->incident_sface() = SFM[e->incident_sface()];
  }
  // Halffacet update
  CGAL_forall_halffacets(f,*this) {
    f->twin() = FM[f->twin()];
    f->incident_volume() = CM[f->incident_volume()];
    Halffacet_cycle_iterator ftc;
    for(ftc = f->facet_cycles_begin(); ftc !=  f->facet_cycles_end(); ++ftc) {
      if (ftc.is_shalfedge() ) {
	se = SHalfedge_handle(ftc);
	*ftc = make_object(SEM[se]); 
	store_boundary_item(se,ftc); 
      } else if (ftc.is_shalfloop() ) {
	sl = SHalfloop_handle(ftc);
	*ftc = make_object(SLM[sl]); 
	store_boundary_item(sl,ftc); 
      } else CGAL_error_msg("damn wrong boundary item in facet.");
    }
  }

  // Volume update
  CGAL_forall_volumes(c,*this) {
    Shell_entry_iterator sei;
    CGAL_forall_shells_of(sei,c) {
      sf = sei; // conversion from generic iterator to sface const handle
      *sei = make_object(SFM[sf]); 
      store_boundary_item(sf,sei); 
    }
  }

  CGAL_forall_shalfedges(se,*this) {
    se->source() = EM[se->source()];
    se->sprev() = SEM[se->sprev()]; se->snext() = SEM[se->snext()];
    se->incident_sface() = SFM[se->incident_sface()];
    se->twin() = SEM[se->twin()];
    se->prev() = SEM[se->prev()]; se->next() = SEM[se->next()];
    se->facet() = FM[se->facet()];
  }

  CGAL_forall_shalfloops(sl,*this) {
    sl->twin() = SLM[sl->twin()];
    sl->incident_sface() = SFM[sl->incident_sface()];
    sl->facet() = FM[sl->facet()];
  }
  for ( slc = D.shalfloops_begin(), sl = shalfloops_begin();
	slc != D.shalfloops_end(); ++slc, ++sl) {
    /* It is possible that the is_twin() property differs for equivalent 
       sloops on both SNC structures.  So, we need to store the correct
       selection mark in the correct (non-twin) facet of a shalfloop pair. */
    CGAL_assertion_code( if( slc->is_twin() == sl->is_twin())
			 CGAL_assertion( slc->mark() == sl->mark()));
    if( !sl->is_twin() && slc->is_twin()) sl->mark() = sl->twin()->mark();
  }

  CGAL_forall_sfaces(sf,*this) {
	  sf->center_vertex() = VM[sf->center_vertex()];
	  sf->volume() = CM[sf->volume()];
	  SFace_cycle_iterator sfc;
	  for(sfc = sf->sface_cycles_begin(); 
		  sfc != sf->sface_cycles_end(); ++sfc) 
	  {
		  if (sfc.is_svertex()) { 
			  SVertex_handle sv(sfc);
			  sv = EM[sv];
			  *sfc = make_object(sv);
			  store_sm_boundary_item(sv,sfc);
		  } else if (sfc.is_shalfedge()) { 
			  se = SHalfedge_handle(sfc);
			  se = SEM[se];
			  *sfc = make_object(se);
			  store_sm_boundary_item(se,sfc);
		  } else if (sfc.is_shalfloop()) { 
			  sl = SHalfloop_handle(sfc);
			  sl = SLM[sl];
			  *sfc = make_object(sl);
			  store_sm_boundary_item(sl,sfc);
		  } else CGAL_error_msg("damn wrong boundary item in sface.");
	  }
  }
}

} //namespace CGAL
#endif // CGAL_SNC_STRUCTURE_H
