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
#ifndef CGAL_SNC_STRUCTURE_H
#define CGAL_SNC_STRUCTURE_H

#include <CGAL/basic.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/Nef_S2/Generic_handle_map.h>
#include <CGAL/Nef_S2/Sphere_geometry.h>
#include <CGAL/Nef_3/SNC_iteration.h>
#include <CGAL/Nef_3/SNC_items.h>
#include <CGAL/Nef_3/nef3_assertions.h>
#include <CGAL/Nef_3/Infimaximal_box.h>
#include <CGAL/Nef_3/SNC_const_decorator.h>
#include <CGAL/Nef_3/SNC_SM_point_locator.h>
#include <CGAL/Nef_2/iterator_tools.h>
#include <CGAL/Union_find.h>
#include <list>

#include <CGAL/Extended_homogeneous_3.h>

#undef _DEBUG
#define _DEBUG 41
#include <CGAL/Nef_3/debug.h>
#include <CGAL/Nef_2/Object_index.h>

CGAL_BEGIN_NAMESPACE

// Const Circulators: 
template <typename HE>
struct move_shalfedge_around_svertex {
  void forward(HE& e) const { e = (e->sprev_->twin_); }
  void backward(HE& e) const { e = (e->twin_->snext_); }
};

template <typename HE>
struct move_shalfedge_around_sface {
  void forward(HE& e)  const { e = (e->snext_); }
  void backward(HE& e) const { e = (e->sprev_); }
};

template <typename HE>
struct move_shalfedge_around_facet {
  void forward(HE& e) const { e = (e->next_); }
  void backward(HE& e) const { e = (e->prev_); }
};

template <class Object, class Hash_map, class Union_find>
void merge_sets( Object o1, Object o2, Hash_map& hash, Union_find& uf) {
  //CGAL_nef3_assertion( hash[o1] != 0 && hash[o2] != 0);
  CGAL_nef3_assertion( hash.is_defined(o1) && hash.is_defined(o2));
  if( !uf.same_set( hash[o1], hash[o2]))
    uf.unify_sets( hash[o1], hash[o2]);
}

/*{\Manpage {SNC_structure}{Items}{Selective Nef Complex}{C}}*/

template <typename Items_>
class SNC_structure {
/*{\Mdefinition The extended Wuerzburg structure is the topological
structure of Nef polyhedra. It is programmed around the local
graphs of the vertices of a Nef polyhedron, which describe the
point set completely. All other concepts are either derived from
the local graph or added for the comfort of the user.}*/

  typedef Items_                 Items;
  typedef SNC_structure<Items>   Self;

  friend class SNC_SM_decorator<Self>;
  friend class SNC_decorator<Self>;
  friend class SNC_io_parser<Self>;

public:
  /*{\Mtypes 7}*/

  typedef SNC_SM_decorator<Self> SM_decorator;
  typedef SNC_decorator<Self>    SNC_decorator;
  typedef SNC_const_decorator<Self> SNC_const_decorator;
  typedef SNC_SM_point_locator<Self> SM_point_locator;

  typedef typename Items::Kernel        Kernel;
  typedef typename Kernel::FT           FT;
  typedef typename Kernel::RT           RT;
  typedef typename Items::Sphere_kernel Sphere_kernel;
  
  typedef Infimaximal_box<typename Is_extended_kernel<Kernel>::value_type, Kernel> Infi_box;
  typedef typename Infi_box::Standard_kernel  Standard_kernel;
  typedef Infimaximal_box<typename Is_extended_kernel<Standard_kernel>::value_type, Standard_kernel> No_box;

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
  typedef typename Items::Mark Mark;
  /*{\Mtypemember attributes of all objects.}*/
  typedef size_t Size_type;
  /*{\Mtypemember size type.}*/

  /*{\Mtext For all objects |Vertex|, |Halfedge|, |Halffacet|, |Volume|
  there are handle and iterator types |xxx_handle|, |xxx_iterator|.
  Additionally all objects of the local graph of a vertex 
  |SVertex|, |SHalfedge|, |SHalfloop|, |SFace| are accessed via handles
  and iterators. There's no type |SHalfloop_iterator|, as there is
  at most one |SLoop| pair per vertex.}*/

  typedef typename Items::template Vertex<Self>    Vertex;
  typedef CGAL::In_place_list<Vertex,false>        Vertex_list;
  typedef CGAL_ALLOCATOR(Vertex)                   Vertex_alloc;
  typedef typename Vertex_list::iterator           Vertex_handle;
  typedef typename Vertex_list::const_iterator     Vertex_const_handle;
  typedef typename Vertex_list::iterator           Vertex_iterator;
  typedef typename Vertex_list::const_iterator     Vertex_const_iterator;

  typedef typename Items::template Halfedge<Self>  Halfedge;
  typedef CGAL::In_place_list<Halfedge,false>      Halfedge_list;
  typedef CGAL_ALLOCATOR(Halfedge)                 Halfedge_alloc;
  typedef typename Halfedge_list::iterator         Halfedge_handle;
  typedef typename Halfedge_list::const_iterator   Halfedge_const_handle;
  typedef typename Halfedge_list::iterator         Halfedge_iterator;
  typedef typename Halfedge_list::const_iterator   Halfedge_const_iterator;

  typedef typename Items::template Halffacet<Self> Halffacet;
  typedef CGAL::In_place_list<Halffacet,false>     Halffacet_list;
  typedef CGAL_ALLOCATOR(Halffacet)                Halffacet_alloc;
  typedef typename Halffacet_list::iterator        Halffacet_handle;
  typedef typename Halffacet_list::const_iterator  Halffacet_const_handle;
  typedef typename Halffacet_list::iterator        Halffacet_iterator;
  typedef typename Halffacet_list::const_iterator  Halffacet_const_iterator;

  typedef typename Items::template Volume<Self>    Volume;
  typedef CGAL::In_place_list<Volume,false>        Volume_list;
  typedef CGAL_ALLOCATOR(Volume)                   Volume_alloc;
  typedef typename Volume_list::iterator           Volume_handle;
  typedef typename Volume_list::const_iterator     Volume_const_handle;
  typedef typename Volume_list::iterator           Volume_iterator;
  typedef typename Volume_list::const_iterator     Volume_const_iterator;

  typedef typename Items::template Halfedge<Self>  SVertex;
  // typedef CGAL_ALLOCATOR(SVertex)                  Svertex_alloc;
  typedef typename Halfedge_list::iterator         SVertex_handle;
  typedef typename Halfedge_list::const_iterator   SVertex_const_handle;
  typedef typename Halfedge_list::iterator         SVertex_iterator;
  typedef typename Halfedge_list::const_iterator   SVertex_const_iterator;

  typedef typename Items::template SHalfedge<Self> SHalfedge;
  typedef CGAL::In_place_list<SHalfedge,false>     SHalfedge_list;
  typedef CGAL_ALLOCATOR(SHalfedge)                SHalfedge_alloc;
  typedef typename SHalfedge_list::iterator        SHalfedge_handle;
  typedef typename SHalfedge_list::const_iterator  SHalfedge_const_handle;
  typedef typename SHalfedge_list::iterator        SHalfedge_iterator;
  typedef typename SHalfedge_list::const_iterator  SHalfedge_const_iterator;

  typedef typename Items::template SHalfloop<Self> SHalfloop;
  typedef CGAL::In_place_list<SHalfloop,false>     SHalfloop_list;
  typedef CGAL_ALLOCATOR(SHalfloop)                SHalfloop_alloc;
  typedef typename SHalfloop_list::iterator        SHalfloop_handle;
  typedef typename SHalfloop_list::const_iterator  SHalfloop_const_handle;
  typedef typename SHalfloop_list::iterator        SHalfloop_iterator;
  typedef typename SHalfloop_list::const_iterator  SHalfloop_const_iterator;

  typedef typename Items::template SFace<Self>     SFace;
  typedef CGAL::In_place_list<SFace,false>         SFace_list;
  typedef CGAL_ALLOCATOR(SFace)                    SFace_alloc;
  typedef typename SFace_list::iterator            SFace_handle;
  typedef typename SFace_list::const_iterator      SFace_const_handle;
  typedef typename SFace_list::iterator            SFace_iterator;
  typedef typename SFace_list::const_iterator      SFace_const_iterator;

  typedef CGAL::Object_handle Object_handle;
  typedef CGAL::Object_handle SObject_handle;

  // a generic handle 

  typedef std::list<Object_handle>    Object_list;
  typedef std::list<Object_handle>    SObject_list;
  typedef Object_list::iterator       Object_iterator;
  typedef Object_list::const_iterator Object_const_iterator;
  typedef Object_list::const_iterator Object_const_handle;
  typedef Object_list::const_iterator SObject_const_handle;

  // Halffacet triangle
  class Halffacet_triangle_handle : public Halffacet_handle {
    typedef Halffacet_handle Base;
    Triangle_3 triangle;
  public:
    Halffacet_triangle_handle() : Base() {}
    Halffacet_triangle_handle( Halffacet_handle h, Triangle_3 t = Triangle_3()) :
      Base(h), triangle(t) {}
    Triangle_3 get_triangle() { return triangle; }
  };
  
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
    { CGAL_nef3_assertion_msg(0,"not impl."); }
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
    Object_handle& operator*() const { return Ibase::operator*(); }
    Object_handle  operator->() const 
    { CGAL_nef3_assertion_msg(0,"not impl."); }
  };

  class SFace_cycle_iterator : public Object_iterator 
  /*{\Mtypemember a generic iterator to an object in the boundary
  of a sface. Convertible to |SObject_handle|.}*/
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

    operator SObject_handle() const { return Ibase::operator*(); }
    SObject_handle& operator*() const { return Ibase::operator*(); }
    SObject_handle  operator->() const 
    { CGAL_nef_assertion_msg(0,"not impl."); }
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

    operator SObject_handle() const { return Ibase::operator*(); }
    const SObject_handle& operator*() const { return Ibase::operator*(); }
    SObject_handle  operator->() const 
    { CGAL_nef_assertion_msg(0,"not impl."); }
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
      CGAL_nef3_assertion( CGAL::assign(f,Ibase::operator*()) );
      CGAL::assign(f,Ibase::operator*()); return f; }

    operator SObject_handle() const { return Ibase::operator*(); }
    SObject_handle& operator*() const { return Ibase::operator*(); }
    SObject_handle  operator->() const 
    { CGAL_nef_assertion_msg(0,"not impl."); }
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
      CGAL_nef3_assertion( CGAL::assign(f,Ibase::operator*()) );
      CGAL::assign(f,Ibase::operator*()); 
      return SFace_const_handle(f); }

    operator SObject_handle() const { return Ibase::operator*(); }
    SObject_handle& operator*() const { return Ibase::operator*(); }
    SObject_handle  operator->() const 
    { CGAL_nef_assertion_msg(0,"not impl."); }
  };

  typedef CircFromIt<
        SHalfedge_const_iterator, 
        move_shalfedge_around_svertex<SHalfedge_const_iterator> > 
        SHalfedge_around_svertex_const_circulator;

  typedef CircFromIt<
        SHalfedge_const_iterator, 
        move_shalfedge_around_sface<SHalfedge_const_iterator> > 
        SHalfedge_around_sface_const_circulator;

  typedef CircFromIt<SHalfedge_const_iterator, 
          move_shalfedge_around_facet<SHalfedge_const_iterator> > 
          SHalfedge_around_facet_const_circulator;

  // Mutable Circulators: 

  typedef CircFromIt<
        SHalfedge_iterator, 
        move_shalfedge_around_svertex<SHalfedge_iterator> > 
        SHalfedge_around_svertex_circulator;

  typedef CircFromIt<
        SHalfedge_iterator, 
        move_shalfedge_around_sface<SHalfedge_iterator> > 
        SHalfedge_around_sface_circulator;

  typedef CircFromIt<SHalfedge_iterator, 
          move_shalfedge_around_facet<SHalfedge_iterator> > 
          SHalfedge_around_facet_circulator;


  /*{\Mcreation 3}*/
  /*{\Mtext |\Mname| is default and copy constructible. Note that copy
  construction means cloning an isomorphic structure and is thus an
  expensive operation.}*/

  SNC_structure() : 
    boundary_item_(undef_), sm_boundary_item_(undef_),
    vertices_(), halfedges_(), halffacets_(), volumes_(),
    shalfedges_(), shalfloops_(), sfaces_(), IO(NULL) {}
  ~SNC_structure() { TRACEN("~SNC_structure: clearing "<<this); clear(); }

  SNC_structure(const Self& D) : 
    boundary_item_(undef_), sm_boundary_item_(undef_),
    vertices_(D.vertices_), halfedges_(D.halfedges_), 
    halffacets_(D.halffacets_), volumes_(D.volumes_),
    shalfedges_(D.shalfedges_), shalfloops_(D.shalfloops_), 
    sfaces_(D.sfaces_), IO(NULL)
  { pointer_update(D); }

  Self& operator=(const Self& D) { 
    if ( this == &D ) 
      return *this;
    clear();
    boundary_item_.clear(undef_);
    sm_boundary_item_.clear(undef_);
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

  void clear() { 
    boundary_item_.clear();
    sm_boundary_item_.clear();
    vertices_.destroy();
    halfedges_.destroy();
    halffacets_.destroy();
    volumes_.destroy();
    shalfedges_.destroy();
    shalfloops_.destroy();
    sfaces_.destroy();
  }

  template <typename H>
  bool is_boundary_object(H h) const
  { return boundary_item_[h]!=undef_; }
  template <typename H>
  bool is_sm_boundary_object(H h) 
  { return sm_boundary_item_[h]!=undef_; }

  template <typename H>
  Object_iterator& boundary_item(H h)
  { return boundary_item_[h]; }
  template <typename H>
  Object_iterator& sm_boundary_item(H h)
  { return sm_boundary_item_[h]; }

  template <typename H>
  void store_boundary_item(H h, Object_iterator o)
  { boundary_item_[h] = o; }
  template <typename H>
  void store_sm_boundary_item(H h, Object_iterator o)
  { sm_boundary_item_[h] = o; }

  template <typename H>
  void undef_boundary_item(H h)
  { CGAL_nef3_assertion(boundary_item_[h]!=undef_);
    boundary_item_[h] = undef_; }
  template <typename H>
  void undef_sm_boundary_item(H h)
  { CGAL_nef3_assertion(sm_boundary_item_[h]!=undef_);
    sm_boundary_item_[h] = undef_; }

  void reset_iterator_hash(Object_iterator it)
  { SVertex_handle sv;
    SHalfedge_handle se;
    SHalfloop_handle sl;

    if ( assign(se,*it) ) { 
      if( is_boundary_object(se)) 
	undef_boundary_item(se); 
      return; 
    }
    if ( assign(sl,*it) ) { 
      if( is_boundary_object(sl)) 
	undef_boundary_item(sl);
      return; 
    }
    if ( assign(sv,*it) ) { 
      if( is_boundary_object(sv)) 
	undef_boundary_item(sv); 
      return; 
    }
  }

  void reset_sm_object_list(Object_list& L)
  { Object_iterator oit;
    CGAL_nef3_forall_iterators(oit,L) reset_sm_iterator_hash(oit);
    L.clear();
  }

  void reset_sm_iterator_hash(Object_iterator it)
  { SVertex_handle sv;
    SHalfedge_handle se;
    SHalfloop_handle sl;
    if ( assign(se,*it) ) { 
      if( is_sm_boundary_object(se)) 
	undef_sm_boundary_item(se); 
      return; 
    }
    if ( assign(sl,*it) ) { 
      if( is_sm_boundary_object(sl)) 
	undef_sm_boundary_item(sl);
      return; 
    }
    if ( assign(sv,*it) ) { 
      if( is_sm_boundary_object(sv)) 
	undef_sm_boundary_item(sv); 
      return; 
    }
  }

  void reset_object_list(Object_list& L)
  { Object_iterator oit;
    CGAL_nef3_forall_iterators(oit,L) reset_iterator_hash(oit);
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

  The macros are then |CGAL_nef3_forall_vertices(v,\Mvar)|, 
  |CGAL_nef3_forall_halfedges(e,\Mvar)|, |CGAL_nef3_forall_edges(e,\Mvar)|,
  |CGAL_nef3_forall_halffacets(f,\Mvar)|, |CGAL_nef3_forall_facets(f,\Mvar)|,
  |CGAL_nef3_forall_volumes(w,\Mvar)|.}*/

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
	    (++volumes_begin())->mark_ == false);
  }

  Vertex_handle new_vertex(const Point_3& p = Point_3(), Mark m = Mark())
  /*{\Mop returns a new vertex at point |p| marked by |m|.}*/ { 
    Vertex_handle vh = new_vertex_only();
    vh->point_at_center_ = p;
    vh->mark_ = m;
    vh->sncp_ = this;
    vh->svertices_begin_ = vh->svertices_last_ = svertices_end();
    vh->shalfedges_begin_ = vh->shalfedges_last_ = shalfedges_end();
    vh->sfaces_begin_ = vh->sfaces_last_ = sfaces_end();
    vh->shalfloop_ = shalfloops_end();
    return vh;
  }

  Halfedge_handle new_halfedge_pair(Vertex_handle v1, Vertex_handle v2,
				    Mark m = Mark())
  /*{\Mop creates a new halfedge pair between the vertices $v_1$
  and $v_2$. The edge is marked by |m|.}*/ { 
    SM_decorator D1(v1);
    SM_decorator D2(v2);
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
    f1->supporting_plane_ = h; f2->supporting_plane_ = h.opposite();
    f1->mark_ = f2->mark_ = m;
    make_twins(f1,f2);
    return f1;
  }

  Volume_handle new_volume(Mark m = Mark())
  /*{\Mop creates a new volume marked with |m|.}*/ { 
    Volume_handle vh = new_volume_only();
    vh->mark_ = m;
    return vh;
  }

  template <typename H>
  void make_twins(H h1, H h2) { 
    h1->twin_ = h2; h2->twin_ = h1; 
  }

  void delete_vertex(Vertex_handle v)
  /*{\Mop deletes the vertex including the objects in its local graph.}*/  { 
    TRACEN("~ deleting vertex "<<&*v<<" from "<<&*this);
    v->clear_local_graph(); 
    delete_vertex_only(v);
    TRACEN("~~ vertex deleted"<<&*v);
  }

  void delete_halfedge_pair(Halfedge_handle e)
  /*{\Mop deletes the halfedge pair of |e,twin(e)|.  Does not care about
  incident objects in the local graph of |source(e)|.}*/ { 
    TRACEN("~ deleting halfedges pair "<<&*e<<", "<<&*(e->twin_)<<
	   " from "<<&*this);
    Halfedge_handle et = e->twin_;
    SM_decorator D1(e->center_vertex_), D2(et->center_vertex_);
    D1.delete_vertex(e);
    D2.delete_vertex(et);
  }

  void delete_halffacet_pair(Halffacet_handle f)
  /*{\Mop deletes the halffacet pair |f,twin(f)|. Does not care about
  boundary cycle objects.}*/ { 
    TRACEN("~ deleting halffacets pair "<<&*f<<", "<<&*(f->twin_)<<
	   " from "<<&*this);
    reset_object_list(f->boundary_entry_objects_);
    reset_object_list(f->twin_->boundary_entry_objects_);
    delete_halffacet_only(f->twin_);
    delete_halffacet_only(f);
  }

  void delete_volume(Volume_handle c)
  /*{\Mop deletes the volume |c|. Does not care about shell objects.}*/ { 
    TRACEN("~ deleting volume "<<&*c<<" from "<<&*this);
    reset_object_list(c->shell_entry_objects_);
    delete_volume_only(c);
  }

  Vertex_alloc vertex_allocator;
  Vertex* get_vertex_node( const Vertex& t) {
    Vertex* p = vertex_allocator.allocate(1);
    vertex_allocator.construct( p, Vertex());
    return p;
  }
  void put_vertex_node( Vertex* p) {
    vertex_allocator.destroy(p);
    vertex_allocator.deallocate( p, 1);
  }

  Halfedge_alloc halfedge_allocator;
  Halfedge* get_halfedge_node( const Halfedge& t) {
    Halfedge* p = halfedge_allocator.allocate(1);
    halfedge_allocator.construct( p, Halfedge());
    return p;
  }
  void put_halfedge_node( Halfedge* p) {
    halfedge_allocator.destroy(p);
    halfedge_allocator.deallocate( p, 1);
  }

  Halffacet_alloc halffacet_allocator;
  Halffacet* get_halffacet_node( const Halffacet& t) {
    Halffacet* p = halffacet_allocator.allocate(1);
    halffacet_allocator.construct( p, Halffacet());
    return p;
  }
  void put_halffacet_node( Halffacet* p) {
    halffacet_allocator.destroy(p);
    halffacet_allocator.deallocate( p, 1);
  }

  Volume_alloc volume_allocator;
  Volume* get_volume_node( const Volume& t) {
    Volume* p = volume_allocator.allocate(1);
    volume_allocator.construct( p, Volume());
    return p;
  }
  void put_volume_node( Volume* p) {
    volume_allocator.destroy(p);
    volume_allocator.deallocate( p, 1);
  }

  SHalfedge_alloc shalfedge_allocator;
  SHalfedge* get_shalfedge_node( const SHalfedge& t) {
    SHalfedge* p = shalfedge_allocator.allocate(1);
    shalfedge_allocator.construct( p, SHalfedge());
    return p;
  }
  void put_shalfedge_node( SHalfedge* p) {
    shalfedge_allocator.destroy(p);
    shalfedge_allocator.deallocate( p, 1);
  }

  SHalfloop_alloc shalfloop_allocator;
  SHalfloop* get_shalfloop_node( const SHalfloop& t) {
    SHalfloop* p = shalfloop_allocator.allocate(1);
    shalfloop_allocator.construct( p, SHalfloop());
    return p;
  }
  void put_shalfloop_node( SHalfloop* p) {
    shalfloop_allocator.destroy(p);
    shalfloop_allocator.deallocate( p, 1);
  }

  SFace_alloc sface_allocator;
  SFace* get_sface_node( const SFace& t) {
    SFace* p = sface_allocator.allocate(1);
    sface_allocator.construct( p, SFace());
    return p;
  }
  void put_sface_node( SFace* p) {
    sface_allocator.destroy(p);
    sface_allocator.deallocate( p, 1);
  }


  Vertex_handle new_vertex_only() { 
    vertices_.push_back( * get_vertex_node(Vertex()));
    TRACEN("  new vertex only "<<&*(--vertices_end()));
    return --vertices_end(); 
  }
  Halfedge_handle new_halfedge_only(Halfedge_handle e)  { 
    Halfedge_handle ne = halfedges_.insert(e, * get_halfedge_node(Halfedge()));
    TRACEN("  after "<<&*e<<" new halfedge only "<<&*ne);
    return ne;
  }
  Halfedge_handle new_halfedge_only()  { 
    TRACEN("  new halfedge only "<<&*(--halfedges_end()));
    halfedges_.push_back( * get_halfedge_node(Halfedge()));
    return --halfedges_end();
  }
  Halffacet_handle new_halffacet_only()  { 
    halffacets_.push_back( * get_halffacet_node(Halffacet()));
    TRACEN("  new halffacet only "<<&*(--halffacets_end()));
    return --halffacets_end(); 
  } 
  Volume_handle new_volume_only()  { 
    volumes_.push_back( * get_volume_node(Volume()));
    TRACEN("  new volume only "<<&*(--volumes_end()));
    return --volumes_end(); 
  }
  SHalfedge_handle new_shalfedge_only()  {
    shalfedges_.push_back( * get_shalfedge_node(SHalfedge()));
    TRACEN("  new shalfedge only "<<&*(--shalfedges_end()));
    return --shalfedges_end();
  }
  SHalfloop_handle new_shalfloop_only()  {
    shalfloops_.push_back( * get_shalfloop_node(SHalfloop()));
    TRACEN("  new shalfloop only "<<&*(--shalfloops_end()));
    return --shalfloops_end(); 
  }
  SFace_handle new_sface_only() {
    sfaces_.push_back( * get_sface_node(SFace()));
    TRACEN("  new sface only "<<&*(--sfaces_end()));
    return --sfaces_end(); 
  }

 
  void delete_vertex_only(Vertex_handle h) {
    TRACEN("~ deleting vertex only "<<&*h<<" from "<<&*this);
    vertices_.erase(h);
    put_vertex_node(&*h);
  }
  void delete_halfedge_only(Halfedge_handle h) { 
    TRACEN("~ deleting halfedge only "<<&*h<<" from "<<&*this);
    CGAL_nef3_assertion(!is_sm_boundary_object(h));
    halfedges_.erase(h);
    put_halfedge_node(&*h);
  }
  void delete_halffacet_only(Halffacet_handle h) { 
    TRACEN("~ deleting halffacet only "<<&*h<<" from "<<&*this);
    halffacets_.erase(h);         
    put_halffacet_node(&*h);
  }
  void delete_volume_only(Volume_handle h) {
    TRACEN("~ deleting volume only "<<&*h<<" from "<<&*this);
    volumes_.erase(h); 
    put_volume_node(&*h);
  }
  void delete_shalfedge_only(SHalfedge_handle h)  { 
    TRACEN("~ deleting shalfedge only "<<&*h<<" from "<<&*this);
    CGAL_nef3_assertion(!is_sm_boundary_object(h));
    shalfedges_.erase(h);
    put_shalfedge_node(&*h);
  }
  void delete_shalfloop_only(SHalfloop_handle h)  { 
    TRACEN("~ deleting shalfloop only "<<&*h<<" from "<<&*this);
    CGAL_nef3_assertion(!is_sm_boundary_object(h));
    shalfloops_.erase(h); 
    put_shalfloop_node(&*h);
  }
  void delete_sface_only(SFace_handle h)  { 
    TRACEN("~ deleting sface only "<<&*h<<" from "<<&*this);
    CGAL_nef3_assertion(!is_boundary_object(h));
    sfaces_.erase(h);
    put_sface_node(&*h);
  }
  
  char PSE(SHalfedge_handle h)
    /* prints a sphere segment */ {
    SNC_decorator D;
    SM_decorator SD;
    TRACE(IO->index(h)<<" @ "<<IO->index(D.vertex(h))<<
	  " "<<D.point(D.vertex(h))<<
	  ", prev " <<IO->index(D.previous(h))<<
	  ", next " <<IO->index(D.next(h))<<
	  ", sprev "<<IO->index(SD.previous(h))<<
	  ", snext "<<IO->index(SD.next(h))<<
	  ", twin " <<IO->index(D.twin(h)));
    return ' ';
  }

  char PFC(SHalfedge_handle e)
    /* prints a facet cycle */ {
    TRACEN("--> Facet cycle begin");
    SHalfedge_around_facet_circulator c(e), cend(c);
    CGAL_For_all(c, cend) {
      TRACEN(PSE(c));
    }
    TRACE("--> Facet cycle end"); 
    return ' ';
  }

  char PFB(Halffacet_handle f)
    /* prints facet boundary entry points */ {
    Halffacet_cycle_iterator fc;
    SHalfedge_handle e; SHalfloop_handle l;
    CGAL_nef3_forall_facet_cycles_of(fc, f) {
      if( assign(e, fc)) {
	TRACE(' '<<IO->index(e)); }
      else if( assign(l, fc)) {
	TRACE(' '<<IO->index(l)<<"(sl)");
      }
    }
    return '.';
  }

  char PSFB(SFace_handle f)
    /* prints sphere face boundary entry points */ {
    SFace_cycle_iterator it;
    CGAL_nef3_forall_sface_cycles_of(it,f)
      if ( it.is_shalfedge() ) TRACE(IO->index(SHalfedge_handle(it))<<' ');
    TRACE(", ");
    CGAL_nef3_forall_sface_cycles_of(it,f)
      if ( it.is_svertex() ) TRACE(IO->index(SVertex_handle(it))<<' ');
    TRACE(", ");
    CGAL_nef3_forall_sface_cycles_of(it,f)
      if ( it.is_shalfloop() ) TRACE(IO->index(SHalfloop_handle(it)));
    return '.';
  }

  typedef typename Union_find< Volume_handle>::handle UFH_volume;
  typedef typename Union_find< Halffacet_handle>::handle UFH_facet;
  typedef typename Union_find< SFace_handle>::handle UFH_sface;

  void remove_f_including_all_edge_uses_in_its_boundary_cycles
    ( Halffacet_handle f,
      Unique_hash_map< SFace_handle, UFH_sface>& hash,
      Union_find< SFace_handle>& uf )
    /* removes f and its boundary cycles, and merges up the sphere facets
       incident to them. */ {
    SNC_decorator D;
    Halffacet_cycle_iterator fc;
    CGAL_nef3_forall_facet_cycles_of(fc, f) {
      SHalfedge_handle e;
      SHalfloop_handle l;
      if( assign(e, fc) ) {
	SHalfedge_around_facet_circulator u(e), eend(e);
	CGAL_For_all(u, eend) {
	  SFace_handle fu = D.sface(u), ftu = D.sface(D.twin(u));
	  TRACEN("sfUNION of "<<IO->index(fu)<<" & "<<IO->index(ftu));
	  merge_sets( fu, ftu, hash, uf);
	  SM_decorator SD(D.vertex(u));
	  TRACEN("removing "<<IO->index(u)<<" & "<<IO->index(SD.twin(u)));
	  Halfedge_handle src(SD.source(u)), tgt(SD.target(u));
	  if ( SD.is_closed_at_source(u) ) 
	    SD.set_face( src, fu);
	  if ( SD.is_closed_at_source( SD.twin(u)) ) 
 	    SD.set_face( tgt, fu);
	  /* TO VERIFY: does is_closed_at_source(u) imply is_isolated(src)?
	     if it is true, the svertex face update is not necesary. */
	  SD.delete_edge_pair(u); 
	  if( SD.is_isolated(src))
	    //	    SD.delete_vertex_only(src);
	    SD.set_face(src,fu);
	  SM_decorator SD2(D.vertex(tgt));
	  if( SD2.is_isolated(tgt))
	    // SD.delete_vertex_only(tgt);
	    SD2.set_face(tgt,fu);
	  /* TO VERIFY: can both svertices be isolated at the same time? */
      }
      }
      else if( assign(l, fc)) {
	// this code is currenlty not used, but it is potentially need 
	// in the future, e.g for complex marks or a relative interior 
	// function 
	SFace_handle fu = D.sface(l), ftu = D.sface(D.twin(l));
	TRACEN("UNION of "<<IO->index(fu)<<" & "<<IO->index(ftu));
	merge_sets( fu, ftu, hash, uf);
	SM_decorator SD(D.vertex(l));
	TRACEN("removing "<<IO->index(l)<<" & "<<IO->index(SD.twin(l)));
	SD.delete_loop_only();
      }
    }
    TRACEN("removing "<<IO->index(f)<<" & "<<IO->index(D.twin(f)));
    delete_halffacet_pair(f);
    return;
  }

  bool is_part_of_volume(Vertex_handle v) 
    /* determines if a vertex v is part of a volume, cheking if its local
       graph is trivial (only one sface with no boundary). */  {
    SM_decorator SD(v);
    CGAL_nef3_assertion( !is_empty_range( SD.sfaces_begin(), SD.sfaces_end()));
    if( is_empty_range( SD.svertices_begin(), SD.svertices_end()) &&
	is_empty_range( SD.shalfedges_begin(), SD.shalfedges_end()) &&
	!SD.has_loop())
      return true;
    return false;
  }

  bool is_part_of_facet(Vertex_handle v) 
    /* determines if a vertex v is part of a the relative interior of a 
       facet, checking if its local graph consists just of a sloop and
       two incident sfaces. */ {
    SM_decorator SD(v);
    CGAL_nef3_assertion( !is_empty_range( SD.svertices_begin(),
					  SD.svertices_end()) ||
			 is_empty_range( SD.shalfedges_begin(),
					 SD.shalfedges_end()));
    return( SD.has_loop() &&
	    is_empty_range( SD.svertices_begin(), SD.svertices_end()));
  }

  bool is_part_of_edge(Vertex_handle v) {
    /* determines if a vertex v is part of a edge, checking at its local 
       graph for exactly two antipodal vertices  */

    SM_decorator SD(v);
    if(SD.has_loop())
      return false;
    if(SD.svertices_begin() == SD.svertices_end())
      return false;
    if(++(SD.svertices_begin()) == SD.svertices_end())
      return false;
    
    TRACE(SNC_decorator(*this).point(v)<<" is in edge interior? ");
    SVertex_iterator sv(SD.svertices_begin());
    SVertex_handle p1(sv++), p2(sv++);
    if( sv != SD.svertices_end())
      return false;
    
    TRACE("has two svertices ");
    Sphere_point sp1(SD.point(p1)), sp2(SD.point(p2));
    return (sp1 == sp2.antipode());
  }

  bool vertex_simplification(bool snc_computed = true) {
    bool simplified = false;

#ifdef _DEBUG
    delete IO;
    IO = new SNC_io_parser<SNC_structure>(std::cerr, *this);
#endif

    SNC_decorator D(*this);

    Vertex_iterator v = (*this).vertices_begin();
    while( v != (*this).vertices_end()) {
      SM_decorator SD(v);
      Vertex_iterator v_next(v);
      v_next++;
      if( is_part_of_volume(v)) {
	TRACEN("mark("<<IO->index(v)<<")="<<D.mark(v)<<", "<<
	       "mark("<<IO->index(D.volume(SD.sfaces_begin()))<<")="<<
	       SD.mark(SD.sfaces_begin()));
	if(D.mark(v) == SD.mark(SD.sfaces_begin())) {
	  TRACEN("removing isolated vertex "<<IO->index(v));
	  delete_vertex(v);
	  simplified = true;
	}
      }
      else if( is_part_of_facet(v)) {
	if( D.mark(v) == SD.mark(SD.shalfloop())) {
	  TRACEN("removing "<<IO->index(v)<<
		 " on facet ");
	  delete_vertex(v);
	  simplified = true;
	}
      }
      else if( is_part_of_edge(v)) {
	SVertex_iterator sv(SD.svertices_begin());
	Halfedge_handle e1(sv++), e2(sv++);
	CGAL_nef3_assertion( sv == SD.svertices_end());
	if( D.mark(e1) == D.mark(v) && D.mark(v) == D.mark(e2)) {
	  TRACEN("merging "<<IO->index(e1)<<" & "<<IO->index(e2)<<
		 " in "<<IO->index(v));
	  if(snc_computed)
	    merge_halfedge_pairs( e1, e2);
	  else
	    delete_vertex(v);
	  simplified = true;
	}
      }
      v = v_next;
    }
    return simplified;
  }
  
  bool simplify() {

    bool update_facets  =  false;
    bool update_sfaces  =  false;
    bool update_volumes =  false;

    TRACEN(">>> simplifying");
    SNC_decorator D(*this);

#ifdef _DEBUG
    delete IO;
    IO = new SNC_io_parser<SNC_structure>(std::cerr, *this);
#endif
    
    Unique_hash_map< Volume_handle, UFH_volume> hash_volume;
    Unique_hash_map< Halffacet_handle, UFH_facet> hash_facet;
    Unique_hash_map< SFace_handle, UFH_sface> hash_sface;
    Union_find< Volume_handle> uf_volume;
    Union_find< Halffacet_handle> uf_facet;
    Union_find< SFace_handle> uf_sface;

    /* We discard  the information about boundary entry points, first
       on volumes, facets on sfacets.  Since during the volumes simplification
       is required the remotion of facet cycles, the information about those
       cycles is keep until the this simplification step is performed. */

    boundary_item_.clear();
    sm_boundary_item_.clear();

    Volume_iterator c;
    CGAL_nef3_forall_volumes( c, *this) {
      hash_volume[c] = uf_volume.make_set(c);
      reset_object_list(c->shell_entry_objects_);
    }
    SFace_iterator sf;
    CGAL_nef3_forall_sfaces( sf, *this) {
      hash_sface[sf] = uf_sface.make_set(sf);
      reset_sm_object_list(sf->boundary_entry_objects_);
    }

    /* 
     * Volumes simplification 
     */

    Halffacet_handle f(D.halffacets_begin());
    while( f != D.halffacets_end() && f->is_twin())
      f++;
    while( f != D.halffacets_end()) {
      CGAL_nef3_assertion( !f->is_twin());
      Halffacet_iterator f_next(f);
      do
	f_next++;
      while( f_next != D.halffacets_end() && f_next->is_twin());
      CGAL_nef3_assertion( f != D.twin(f));
      Volume_handle c1 = D.volume(f), c2 = D.volume(D.twin(f));
      TRACEN(" mark("<<IO->index(c1)<<")="<<D.mark(c1)<<
      	     " mark("<<IO->index(f) <<")="<<D.mark(f) <<
	     " mark("<<IO->index(c2)<<")="<<D.mark(c2)<<
	     " is_twin(f)="<<f->is_twin());
      if( D.mark(c1) == D.mark(f) && D.mark(f) == D.mark(c2)
	  && !D.is_infbox_facet(f)) {
	merge_sets( c1, c2, hash_volume, uf_volume);
	remove_f_including_all_edge_uses_in_its_boundary_cycles
	  (f, hash_sface, uf_sface);
	update_sfaces = update_volumes = true;
	TRACEN("UNION of "<<IO->index(c1)<<" & "<<IO->index(c2));
      }
      f = f_next;
    }
    
    CGAL_nef3_forall_halffacets( f, *this) {
      hash_facet[f] = uf_facet.make_set(f);
      reset_object_list(f->boundary_entry_objects_);
    }
    
    /* 
     * Edges simplification
     */

    Halfedge_iterator e(D.halfedges_begin());
    while( e != D.halfedges_end() && e->is_twin())
      e++;
    while( e != (*this).halfedges_end()) {
      CGAL_nef3_assertion( !e->is_twin());
      Halfedge_iterator e_next(e);
      do 
	e_next++;
      while( e_next != D.halfedges_end() && e_next->is_twin());
      SM_decorator SD(D.source(e));
      if( SD.is_isolated(e)) {
	if(D.mark(e) == D.mark(D.volume(D.sface(e)))) {
	  TRACEN("removing pair "<<IO->index(e)<<' '<<IO->index(D.twin(e)));
	  delete_halfedge_pair(e);
	  update_facets = true;
	}
      } 
      else { 
	if( D.has_outdeg_two(e)) {
	  SHalfedge_handle e1(SD.first_out_edge(e)); 
	  SHalfedge_handle e2(SD.cyclic_adj_succ(e1));
	  if( SD.circle(e1)==SD.circle(SD.twin(e2)) &&
	      D.mark(e1)==D.mark(e) && D.mark(e)==D.mark(e2)) {
	    Halffacet_handle f1(D.facet(e1)); 
	    Halffacet_handle f2(D.facet(e2));
	    TRACEN("UNION of "<<IO->index(f1)<<" & "<<IO->index(D.twin(f2))<<
		   " ("<<IO->index(D.twin(f1))<<" & "<<IO->index(f2)<<")");
	    merge_sets( f1, D.twin(f2), hash_facet, uf_facet);
	    merge_sets( D.twin(f1), f2, hash_facet, uf_facet);
	    TRACEN("BEFORE "<<PFC(e1)<<std::endl<<PFC(D.twin(e2)));
	    TRACEN("removing "<<IO->index(e));
	    remove_edge_and_merge_facet_cycles(e);
	    update_facets = true;
	    // TRACEN("AFTER "<<PFC(e1)); // e1 not valid after sloop creation
	  }
	}
      }
      e = e_next;
    }

    update_facets = vertex_simplification() || update_facets;
    purge_no_find_objects(hash_volume, hash_facet, hash_sface, uf_volume, 
			  uf_facet, uf_sface);
    create_boundary_links_forall_sfaces( hash_sface, uf_sface);
    create_boundary_links_forall_facets( hash_facet, uf_facet);
    create_boundary_links_forall_volumes( hash_volume, uf_volume);
    
    TRACEN(">>> simplifying done ");

    return update_sfaces || update_facets || update_volumes;
  }
   
  void remove_edge_and_merge_facet_cycles( Halfedge_handle e) {
     SNC_decorator D(*this);
     CGAL_nef3_assertion( D.has_outdeg_two(e));
     Halfedge_handle et = D.twin(e);
     CGAL_nef3_assertion( D.has_outdeg_two(et));
     SM_decorator SD1(D.vertex(e));
     SM_decorator SD2(D.vertex(et));
     SHalfedge_handle e1 = SD1.first_out_edge(e);
     SHalfedge_handle e2 = SD2.next(D.previous(e1));
     merge_sedges_at_target_and_remove_svertex( D.twin(e1), e);
     merge_sedges_at_target_and_remove_svertex( D.twin(e2), et);
   }

   void merge_sedges_at_target_and_remove_svertex( SHalfedge_handle s1,
						   SVertex_handle v) {
     SNC_decorator D(*this);
     SM_decorator SD(D.vertex(v));
     CGAL_nef3_assertion( SD.target(s1) == v);
     SHalfedge_handle s2(SD.next(s1));
     CGAL_nef3_assertion( SD.source(s2) == v);
     if( s1 == s2) {
       TRACEN(IO->index(s1)<<'('<<IO->index(D.twin(s2))<<") to sloop");
       SD.convert_edge_to_loop(s1);
       CGAL_nef3_assertion(SD.shalfloop() != SHalfloop_handle());
       D.add_sloop_to_facet( SD.shalfloop(), D.facet(s1));
       TRACEN(IO->index(s2)<<" removed");
     }
     else {
       CGAL_nef3_assertion( D.has_outdeg_two(v));
       D.link_as_prev_next_pair( s1, D.next(s2));
       TRACEN(IO->index(s1)<<" "<<IO->index(D.next(s2))<<" linked.");
       D.link_as_prev_next_pair( D.twin(D.next(s2)), D.twin(s1));
       TRACEN(IO->index(D.twin(D.next(s2)))<<" "<<
	      IO->index(D.twin(s1))<<" linked.");
       SD.merge_edge_pairs_at_target( s1); // s2 is removed
       TRACEN(IO->index(s2)<<" removed");
     }
   }

   void merge_halfedge_pairs( SVertex_handle p, SVertex_handle q) {
     SNC_decorator D(*this);
     CGAL_nef3_assertion( D.vertex(p) == D.vertex(q));
     Vertex_handle v(D.vertex(p)); 
     CGAL_nef3_assertion( is_part_of_edge(v));
     SM_decorator SD(v);
     SHalfedge_around_svertex_circulator s(SD.first_out_edge(p)), se(s);
     CGAL_For_all( s, se) {
       D.link_as_prev_next_pair( D.previous(s), D.next(s));
       D.link_as_prev_next_pair( D.previous(SD.twin(s)), D.next(SD.twin(s)));
     }
     D.make_twins( D.twin(p), D.twin(q));
     SD.delete_vertex(p);
     SD.delete_vertex(q);
     delete_vertex(v);
   }

   void purge_no_find_objects( 
      Unique_hash_map< Volume_handle, UFH_volume>& hash_volume,
      Unique_hash_map< Halffacet_handle, UFH_facet>& hash_facet,
      Unique_hash_map< SFace_handle, UFH_sface>& hash_sface,
      Union_find< Volume_handle>& uf_volume,
      Union_find< Halffacet_handle>& uf_facet,
      Union_find< SFace_handle>& uf_sface ) {
     SNC_decorator D(*this);
     SFace_iterator sf;
     std::list<SFace_handle> sflist;
     CGAL_nef3_forall_sfaces( sf, *this) {
       if( uf_sface.find(hash_sface[sf]) != hash_sface[sf]) {
	 TRACEN("no find object "<<IO->index(sf));
	 sflist.push_back(sf);
       }
     }

     typename std::list<SFace_handle>::const_iterator sfli;
     for(sfli = sflist.begin(); sfli != sflist.end(); sfli++){     
       SM_decorator SD(D.vertex(*sfli));
       SD.delete_face_only(*sfli);
     }

     Halffacet_iterator f;
     std::list<Halffacet_handle> flist;
     CGAL_nef3_forall_facets( f, *this) {
       TRACEN("facet "<<IO->index(f));
       if( uf_facet.find(hash_facet[f]) != hash_facet[f]) {
	 TRACEN("no find object "<<IO->index(f));
	 flist.push_back(f);
       }
     }
     
     typename std::list<Halffacet_handle>::const_iterator fli;
     for(fli = flist.begin(); fli != flist.end(); fli++)
       delete_halffacet_pair(*fli);

     Volume_iterator c;
     std::list<Volume_handle> clist;
     CGAL_nef3_forall_volumes( c, *this) {
       if( uf_volume.find(hash_volume[c]) != hash_volume[c]) {
	 TRACEN("no find object "<<IO->index(c));
	 clist.push_back(c);
       }
     }

     typename std::list<Volume_handle>::const_iterator cli;
     for(cli = clist.begin(); cli != clist.end(); cli++){     
       delete_volume(*cli);
     }
   }

  void create_boundary_links_forall_sfaces(
      Unique_hash_map< SFace_handle, UFH_sface>& hash,
      Union_find< SFace_handle>& uf ) {
    Unique_hash_map< SHalfedge_handle, bool> linked(false);
    SNC_decorator D(*this);
    SHalfedge_iterator e;
    CGAL_nef3_forall_shalfedges(e, *this) {
      if( linked[e])
	continue;
      SM_decorator SD(D.vertex(e));
      SFace_handle sf = *(uf.find(hash[D.sface(e)]));
      CGAL_nef3_assertion( sf != SFace_handle());
      SHalfedge_around_sface_circulator c(e), cend(c);
      CGAL_For_all( c, cend) {
	SD.set_face(c, sf);
	linked[c] = true;
      }
      SD.store_boundary_object( e, sf);
    }

    SVertex_handle sv;
    CGAL_nef3_forall_svertices(sv, *this) {
      SM_decorator SD(D.vertex(sv));
      if( SD.is_isolated(sv)) {
	SFace_handle sf = *(uf.find(hash[SD.face(sv)])); 
	CGAL_nef3_assertion( sf != SFace_handle());
	SD.set_face( sv, sf);
	SD.store_boundary_object( sv, sf);
      }
    }

    SHalfloop_handle sl;
    CGAL_nef3_forall_shalfloops(sl, *this) {
      SM_decorator SD(D.vertex(sl));
      SD.store_boundary_object( sl, SD.face(sl));
    }
  }

  void create_boundary_links_forall_facets(
      Unique_hash_map< Halffacet_handle, UFH_facet>& hash,
      Union_find< Halffacet_handle>& uf) {
    Unique_hash_map< SHalfedge_handle, bool> linked(false);
    SNC_decorator D(*this);
    SHalfedge_iterator u;
    CGAL_nef3_forall_shalfedges(u, *this) {
      if( linked[u])
	continue;
      /* set find(f) as incident facet of every edge use on the cycle of u */
      SHalfedge_handle u_min = u;
      Halffacet_handle f = *(uf.find(hash[D.facet(u)]));
      SHalfedge_around_facet_circulator c(u), cend(c);
      CGAL_For_all( c, cend) {
	D.set_facet( c, f);
	if( lexicographically_xyz_smaller(D.point(D.vertex(c)), D.point(D.vertex(u_min))))
	  u_min = c;
	linked[c] = true;
      }
      /* store the edge use at the lexicographicaly minimum facet vertex, as
	 a cycle entry of f.  The outermost cycle is stored at first
	 on the facet's cycles list. */
      if( is_empty_range( f->boundary_entry_objects_.begin(),
			  f->boundary_entry_objects_.end())) {
	D.store_boundary_object( u_min, f);
	TRACEN("new outer cycle min. vertex: "<<D.point(D.vertex(u_min)));
      }
      else {
	SHalfedge_handle f_sedge;
	CGAL_nef3_assertion( assign( f_sedge, 
				     f->boundary_entry_objects_.front()));
	assign( f_sedge, f->boundary_entry_objects_.front());
	if( lexicographically_xyz_smaller(D.point(D.vertex(u_min)), D.point(D.vertex(f_sedge))))
	  D.store_as_first_boundary_object( u_min, f);
	else
	  D.store_boundary_object( u_min, f);
      }
    }
    SHalfloop_iterator l;
    CGAL_nef3_forall_shalfloops( l, *this) {
      Halffacet_handle f = *(uf.find(hash[D.facet(l)]));
      D.set_facet( l, f);
      D.store_boundary_object( l, f);
    }
  }

  void create_boundary_links_forall_volumes( 
      Unique_hash_map< Volume_handle, UFH_volume>& hash,
      Union_find< Volume_handle>& uf) {
    typedef typename SNC_decorator::Shell_volume_setter Volume_setter;
    //   typedef Unique_hash_map< SFace_handle, bool> SFace_map;
    //  SFace_map linked(false);

#ifdef _DEBUG
    delete IO;
    IO = new SNC_io_parser<SNC_structure>(std::cerr, *this);
#endif

    SNC_decorator D(*this);
    Volume_setter setter(D);

    SFace_iterator sf;
    Volume_handle c;
    CGAL_nef3_forall_sfaces(sf, *this) {
      TRACEN("SFace " << IO->index(sf));
      if( setter.is_linked(sf)) continue;
      c = *(uf.find(hash[D.volume(sf)]));
      TRACEN("Volume " << IO->index(c));
      setter.set_volume(c);
      D.visit_shell_objects( sf, setter );      
      D.store_boundary_object( sf, c);
    }
  }
  
    // Returns the bounding box of the finite vertices of the polyhedron.
    // Returns $[-1,+1]^3$ as bounding box if no finite vertex exists.
    Bbox_3  bounded_bbox() const {
        SNC_const_decorator deco(*this);
        Vertex_const_iterator vi = vertices_begin();
        bool first_vertex = true;
        Bbox_3 bbox( -1.0, -1.0, -1.0, 1.0, 1.0, 1.0);
        for ( ; vi != vertices_end(); ++vi) {
            if ( ! deco.is_infbox_vertex(vi)) {
                if ( first_vertex) {
                    bbox = vi->point().bbox();
                    first_vertex = false;
                } else {
                    bbox = bbox + vi->point().bbox();
                }
            }
        }
        return bbox;
    }

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
      CGAL_nef3_forall_halffacets(hf, *this) {
	Halffacet_cycle_iterator fc;
	CGAL_nef3_forall_facet_cycles_of(fc, hf)
	  result += sizeof(*fc) + 2 * sizeof(void*);
      }

      Volume_iterator c;
      CGAL_nef3_forall_volumes(c, *this) {
	Shell_entry_iterator sei;
	CGAL_nef3_forall_shells_of(sei,c) 
	  result += sizeof(*sei) + 2 * sizeof(void*);
      }

      SFace_iterator sf;
      CGAL_nef3_forall_sfaces(sf, *this) {
	SFace_cycle_iterator sfc;
	CGAL_nef3_forall_sface_cycles_of(sfc,sf)
	  result += sizeof(*sfc) + 2 * sizeof(void*);
      }

      return result;
    }


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
      CGAL_nef3_forall_halffacets(hf, *this) {
	Halffacet_cycle_iterator fc;
	CGAL_nef3_forall_facet_cycles_of(fc, hf)
	  result += sizeof(*fc) + 2 * sizeof(void*);
      }

      Volume_iterator c;
      CGAL_nef3_forall_volumes(c, *this) {
	Shell_entry_iterator sei;
	CGAL_nef3_forall_shells_of(sei,c) 
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
      CGAL_nef3_forall_halffacets(hf, *this) {
	Halffacet_cycle_iterator fc;
	CGAL_nef3_forall_facet_cycles_of(fc, hf)
	  result += sizeof(*fc) + 2 * sizeof(void*);
      }

      Volume_iterator c;
      CGAL_nef3_forall_volumes(c, *this) {
	Shell_entry_iterator sei;
	CGAL_nef3_forall_shells_of(sei,c) 
	  result += sizeof(*sei) + 2 * sizeof(void*);
      }

      return result;
    }

    bool test_string(std::string s, std::ifstream& in) {
      std::string s2;
      in >> s2;
      return (s==s2);
    }

    bool load(std::ifstream& in) {

      TRACEN("load");
      
      int i;
      char cc;
      bool OK = true;
      OK = OK && test_string("Selective", in);
      OK = OK && test_string("Nef", in);
      OK = OK && test_string("Complex", in);

      CGAL_nef3_assertion(OK);

      int v;
      OK = OK && test_string("vertices", in);
      in >> v;
      Vertex_handle* V = new Vertex_handle[v];
      for(i = 0; i < v; i++)
	V[i] = new_vertex_only();

      int e;
      OK = OK && test_string("halfedges", in);
      in >> e;
      Halfedge_handle* E = new Halfedge_handle[e];
      for(i = 0; i < e; i++)
	E[i] = new_halfedge_only();

      int f;
      OK = OK && test_string("facets", in);
      in >> f;
      Halffacet_handle* F = new Halffacet_handle[f];
      for(i = 0; i < f; i++)
	F[i] = new_halffacet_only();

      int c;
      OK = OK && test_string("volumes", in);
      in >> c;
      Volume_handle* C = new Volume_handle[c];
      for(i = 0; i < c; i++)
	C[i] = new_volume_only();

      int se;
      OK = OK && test_string("shalfedges", in);
      in >> se;
      SHalfedge_handle* SE = new SHalfedge_handle[se];
      for(i = 0; i < se; i++)
	SE[i] = new_shalfedge_only();
      
      int sl;
      OK = OK && test_string("shalfloops", in);
      in >> sl;
      SHalfloop_handle* SL = new SHalfloop_handle[sl];
      for(i = 0; i < sl; i++)
	SL[i] = new_shalfloop_only();

      int sf;
      OK = OK && test_string("sfaces", in);
      in >> sf;
      SFace_handle* SF = new SFace_handle[sf];
      for(i = 0; i < sf; i++)
	SF[i] = new_sface_only();

      CGAL_nef3_assertion(OK);

      int index;
      RT hx, hy, hz, hw;
      for(i = 0; i < v; i++) {
	in >> index;
	OK = OK && (index == i);
	OK = OK && test_string("{",in);
	Vertex_handle vh = V[i];
	vh->sncp_ = this;

	in >> index;
	vh->svertices_begin_ = (index >= 0 ? E[index] : svertices_end());
	in >> index;
	vh->svertices_last_  = index >= 0 ? E[index] : svertices_end();
	OK = OK && test_string(",",in);
	in >> index;
	vh->shalfedges_begin_ = index >= 0 ? SE[index] : shalfedges_end();
	in >> index;
	vh->shalfedges_last_  = index >= 0 ? SE[index] : shalfedges_end();
	OK = OK && test_string(",",in);
	in >> index;
	vh->sfaces_begin_ = index >= 0 ? SF[index] : sfaces_end();
	in >> index;
	vh->sfaces_last_  = index >= 0 ? SF[index] : sfaces_end();
	OK = OK && test_string(",",in);
	in >> index;
	vh->shalfloop_ = index >= 0 ? SL[index] : shalfloops_end();
	OK = OK && test_string("|",in);
	in >> hx >> hy >> hz >> hw;
	vh->point_at_center_ = Point_3(hx, hy, hz, hw);
	OK = OK && test_string("}",in);
	in >> vh->mark_;
      }

      CGAL_nef3_assertion(OK);

      for(i = 0; i < e; i++) {
	in >> index;
	OK = OK && (index == i);
	OK = OK && test_string("{",in);
	Halfedge_handle eh = E[i];
	
	in >> index;
	eh->twin_ = E[index];
	OK = OK && test_string(",",in);
        in >> index;
	eh->center_vertex_ = V[index];
	OK = OK && test_string(",",in);
	in >> index;
	if(index == 0) {
	  in >> index;
	  eh->out_sedge_ = SE[index];
	} else { 
	  in >> index;
	  eh->incident_sface_ = SF[index];
	}
	OK = OK && test_string("|",in);
	in >> hx >> hy >> hz >> hw;
	eh->point_on_surface_ = Sphere_point(hx,hy,hz);
	OK = OK && test_string("}",in);
	in >> eh->mark_;
      }

      CGAL_nef3_assertion(OK);
      
      for(i = 0; i < f; i++) {
	in >> index;
	OK = OK && (index == i);
	OK = OK && test_string("{",in);
	Halffacet_handle fh = F[i];
      
	in >> index;
	fh->twin_ = F[index];
	OK = OK && test_string(",",in);

	in >> cc;
	while(isdigit(cc)) {
	  in.putback(cc);
	  in >> index;
	  fh->boundary_entry_objects_.push_back(SE[index]);
	  in >> cc;
	}
	
	in >> cc;
	while(isdigit(cc)) {
	  in.putback(cc);
	  in >> index;
	  fh->boundary_entry_objects_.push_back(SL[index]);
	  in >> cc;
	}

	in >> index;
	fh->volume_ = C[index];
	OK = OK && test_string("|",in);
	in >> hx >> hy >> hz >> hw;
	fh->supporting_plane_ = Plane_3(hx,hy,hz,hw);
	OK = OK && test_string("}",in);
	in >> fh->mark_;
      }

      CGAL_nef3_assertion(OK);

      for(i = 0; i < c; i++) {
	in >> index;
	OK = OK && (index == i);
	OK = OK && test_string("{",in);
	Volume_handle ch = C[i];

	in >> cc;
	while(isdigit(cc)) {
	  in.putback(cc);
	  in >> index;
	  ch->shell_entry_objects_.push_back(SF[index]);
	  in >> cc;
	}
	in >> ch->mark_;
      }	

      CGAL_nef3_assertion(OK);

      for(i = 0; i < se; i++) {
	in >> index;
	OK = OK && (index == i);
	OK = OK && test_string("{",in);
	SHalfedge_handle seh = SE[i];
	
	in >> index;
	seh->twin_ = SE[index];
	OK = OK && test_string(",",in);
	in >> index;
	seh->sprev_ = SE[index];
	OK = OK && test_string(",",in);
	in >> index;
	seh->snext_ = SE[index];
	OK = OK && test_string(",",in);
	in >> index;
	seh->source_ = E[index];
	OK = OK && test_string(",",in);
	in >> index;
	seh->incident_sface_ = SF[index];
	OK = OK && test_string(",",in);
	in >> index;
	seh->prev_ = SE[index];
	OK = OK && test_string(",",in);
	in >> index;
	seh->next_ = SE[index];
	OK = OK && test_string(",",in);
	in >> index;
	seh->incident_facet_ = F[index];
	OK = OK && test_string("|",in);
	in >> hx >> hy >> hz >> hw;
	seh->circle_ = Sphere_circle(Plane_3(hx,hy,hz,hw));
	OK = OK && test_string("}",in);
	in >> seh->mark_;
      }

      CGAL_nef3_assertion(OK);

      for(i = 0; i < sl; i++) {
	in >> index;
	OK = OK && (index == i);
	OK = OK && test_string("{",in);
	SHalfloop_handle slh = SL[i];
	
	in >> index;
	slh->twin_ = SL[index];
	OK = OK && test_string(",",in);
	in >> index;
	slh->incident_sface_ = SF[index];
	OK = OK && test_string(",",in);
	in >> index;
	slh->incident_facet_ = F[index];
	OK = OK && test_string("|",in);	
	in >> hx >> hy >> hz >> hw;
	slh->circle_ = Sphere_circle(Plane_3(hx,hy,hz,hw));	
	OK = OK && test_string("}",in);	
	in >> slh->mark_;
      }

      CGAL_nef3_assertion(OK);

      for(i = 0; i < sf; i++) {
	in >> index;
	OK = OK && (index == i);
	OK = OK && test_string("{",in);
	SFace_handle sfh = SF[i];
	
	in >> index;
	sfh->center_vertex_ = V[index];
	OK = OK && test_string(",",in);

	in >> cc;
	while(isdigit(cc)) {
	  in.putback(cc);
	  in >> index;
	  sfh->boundary_entry_objects_.push_back(SE[index]);
	  in >> cc;
	}

	in >> cc;
	while(isdigit(cc)) {
	  in.putback(cc);
	  in >> index;
	  sfh->boundary_entry_objects_.push_back(E[index]);
	  in >> cc;
	}

	in >> cc;
	while(isdigit(cc)) {
	  in.putback(cc);
	  in >> index;
	  sfh->boundary_entry_objects_.push_back(SL[index]);
	  in >> cc;
	}

	in >> index;
	sfh->incident_volume_ = C[index];
	OK = OK && test_string("}",in);	
	in >> sfh->mark_;
      }

      CGAL_nef3_assertion(OK);

      for(i = 0; i<v; i++) {
	SM_point_locator PL(V[i]);
	PL.init_marks_of_halfspheres(); 
      }

      CGAL_nef3_assertion(OK);

      delete[] V;
      delete[] E;
      delete[] F;
      delete[] C;
      delete[] SE;
      delete[] SL;
      delete[] SF;

      return OK;
    }

    bool load_simple(std::ifstream& in) {

      bool addInfiBox = Infi_box::extended_Kernel();

      int i,z;
      char cc;
      bool OK = true;
      OK = OK && test_string("Selective", in);
      OK = OK && test_string("Nef", in);
      OK = OK && test_string("Complex", in);

      int v;
      OK = OK && test_string("vertices", in);
      in >> v;
      z = addInfiBox ? v+8 : v;
      Vertex_handle* V = new Vertex_handle[z];
      for(i = 0; i < z; i++) {
	V[i] = new_vertex_only();
	V[i]->sncp_ = this;
      }

      int e;
      OK = OK && test_string("halfedges", in);
      in >> e;
      z = addInfiBox ? e+24 : e;
      Halfedge_handle* E = new Halfedge_handle[z];
      for(i = 0; i < z; i++)
	E[i] = new_halfedge_only();

      int f;
      OK = OK && test_string("facets", in);
      in >> f;
      z = addInfiBox ? f+12 : f;
      Halffacet_handle* F = new Halffacet_handle[z];
      for(i = 0; i < z; i++)
	F[i] = new_halffacet_only();

      int c;
      OK = OK && test_string("volumes", in);
      in >> c;
      z = addInfiBox ? c+1 : c;
      Volume_handle* C = new Volume_handle[z];

      for(i = 0; i < z; i++)
	C[i] = new_volume_only();

      int se;
      OK = OK && test_string("shalfedges", in);
      in >> se;
      z = addInfiBox ? se+48 : se;
      SHalfedge_handle* SE = new SHalfedge_handle[z];
      for(i = 0; i < z; i++)
	SE[i] = new_shalfedge_only();
      
      int sl;
      OK = OK && test_string("shalfloops", in);
      in >> sl;
      SHalfloop_handle* SL = new SHalfloop_handle[sl];
      for(i = 0; i < sl; i++)
	SL[i] = new_shalfloop_only();

      int sf;
      OK = OK && test_string("sfaces", in);
      in >> sf;
      z = addInfiBox ? sf+16 : sf;
      SFace_handle* SF = new SFace_handle[z];
      for(i = 0; i < z; i++)
	SF[i] = new_sface_only();

      int index;
      typename Standard_kernel::RT hx, hy, hz, hw;
      for(i = 0; i < v; i++) {
	in >> index;
	OK = OK && (index == i);
	OK = OK && test_string("{",in);
	Vertex_handle vh = V[i];

	in >> index;
	vh->svertices_begin_ = (index >= 0 ? E[index] : svertices_end());
	in >> index;
	vh->svertices_last_  = index >= 0 ? E[index] : svertices_end();
	OK = OK && test_string(",",in);
	in >> index;
	vh->shalfedges_begin_ = index >= 0 ? SE[index] : shalfedges_end();
	in >> index;
	vh->shalfedges_last_  = index >= 0 ? SE[index] : shalfedges_end();
	OK = OK && test_string(",",in);
	in >> index;
	vh->sfaces_begin_ = index >= 0 ? SF[index] : sfaces_end();
	in >> index;
	vh->sfaces_last_  = index >= 0 ? SF[index] : sfaces_end();
	OK = OK && test_string(",",in);
	in >> index;
	vh->shalfloop_ = index >= 0 ? SL[index] : shalfloops_end();
	OK = OK && test_string("|",in);
	in >> hx >> hy >> hz >> hw;
	vh->point_at_center_ = Point_3(hx, hy, hz, hw);
	OK = OK && test_string("}",in);
	in >> vh->mark_;
      }

      for(i = 0; i < e; i++) {
	in >> index;
	OK = OK && (index == i);
	OK = OK && test_string("{",in);
	Halfedge_handle eh = E[i];
	
	in >> index;
	eh->twin_ = E[index];
	OK = OK && test_string(",",in);
        in >> index;
	eh->center_vertex_ = V[index];
	OK = OK && test_string(",",in);
	in >> index;
	if(index == 0) {
	  in >> index;
	  eh->out_sedge_ = SE[index];
	} else { 
	  in >> index;
	  eh->incident_sface_ = SF[index];
	}
	OK = OK && test_string("|",in);
	in >> hx >> hy >> hz >> hw;
	eh->point_on_surface_ = Sphere_point(hx,hy,hz);
	OK = OK && test_string("}",in);
	in >> eh->mark_;
      }
      
      for(i = 0; i < f; i++) {
	in >> index;
	OK = OK && (index == i);
	OK = OK && test_string("{",in);
	Halffacet_handle fh = F[i];
      
	in >> index;
	fh->twin_ = F[index];
	OK = OK && test_string(",",in);

	in >> cc;
	while(isdigit(cc)) {
	  in.putback(cc);
	  in >> index;
	  fh->boundary_entry_objects_.push_back(SE[index]);
	  in >> cc;
	}
	
	in >> cc;
	while(isdigit(cc)) {
	  in.putback(cc);
	  in >> index;
	  fh->boundary_entry_objects_.push_back(SL[index]);
	  in >> cc;
	}

	in >> index;
	fh->volume_ = C[index+addInfiBox];
	OK = OK && test_string("|",in);
	in >> hx >> hy >> hz >> hw;
	fh->supporting_plane_ = Plane_3(hx,hy,hz,hw);
	OK = OK && test_string("}",in);
	in >> fh->mark_;
      }

      for(i = 0; i < c; i++) {
	in >> index;
	OK = OK && (index == i);
	OK = OK && test_string("{",in);
	Volume_handle ch = C[i+addInfiBox];

	in >> cc;
	while(isdigit(cc)) {
	  in.putback(cc);
	  in >> index;
	  ch->shell_entry_objects_.push_back(SF[index]);
	  in >> cc;
	}
	in >> ch->mark_;
      }	

      for(i = 0; i < se; i++) {
	in >> index;
	OK = OK && (index == i);
	OK = OK && test_string("{",in);
	SHalfedge_handle seh = SE[i];
	
	in >> index;
	seh->twin_ = SE[index];
	OK = OK && test_string(",",in);
	in >> index;
	seh->sprev_ = SE[index];
	OK = OK && test_string(",",in);
	in >> index;
	seh->snext_ = SE[index];
	OK = OK && test_string(",",in);
	in >> index;
	seh->source_ = E[index];
	OK = OK && test_string(",",in);
	in >> index;
	seh->incident_sface_ = SF[index];
	OK = OK && test_string(",",in);
	in >> index;
	seh->prev_ = SE[index];
	OK = OK && test_string(",",in);
	in >> index;
	seh->next_ = SE[index];
	OK = OK && test_string(",",in);
	in >> index;
	seh->incident_facet_ = F[index];
	OK = OK && test_string("|",in);
	in >> hx >> hy >> hz >> hw;
	seh->circle_ = Sphere_circle(Plane_3(hx,hy,hz,hw));
	OK = OK && test_string("}",in);
	in >> seh->mark_;
      }

      for(i = 0; i < sl; i++) {
	in >> index;
	OK = OK && (index == i);
	OK = OK && test_string("{",in);
	SHalfloop_handle slh = SL[i];
	
	in >> index;
	slh->twin_ = SL[index];
	OK = OK && test_string(",",in);
	in >> index;
	slh->incident_sface_ = SF[index];
	OK = OK && test_string(",",in);
	in >> index;
	slh->incident_facet_ = F[index];
	OK = OK && test_string("|",in);	
	in >> hx >> hy >> hz >> hw;
	slh->circle_ = Sphere_circle(Plane_3(hx,hy,hz,hw));	
	OK = OK && test_string("}",in);	
	in >> slh->mark_;
      }

      for(i = 0; i < sf; i++) {
	in >> index;
	OK = OK && (index == i);
	OK = OK && test_string("{",in);
	SFace_handle sfh = SF[i];
	
	in >> index;
	sfh->center_vertex_ = V[index];
	OK = OK && test_string(",",in);

	in >> cc;
	while(isdigit(cc)) {
	  in.putback(cc);
	  in >> index;
	  sfh->boundary_entry_objects_.push_back(SE[index]);
	  in >> cc;
	}

	in >> cc;
	while(isdigit(cc)) {
	  in.putback(cc);
	  in >> index;
	  sfh->boundary_entry_objects_.push_back(E[index]);
	  in >> cc;
	}

	in >> cc;
	while(isdigit(cc)) {
	  in.putback(cc);
	  in >> index;
	  sfh->boundary_entry_objects_.push_back(SL[index]);
	  in >> cc;
	}

	in >> index;
	sfh->incident_volume_ = C[index+addInfiBox];
	OK = OK && test_string("}",in);	
	in >> sfh->mark_;
      }

      if(addInfiBox) {
	
	for(i=0; i<8; i++) {
	  Vertex_handle vh = V[v+i];
	  vh->svertices_begin_ = E[e+3*i];
	  vh->svertices_last_  = E[e+3*i+2];
	  vh->shalfedges_begin_ = SE[se+6*i];
	  vh->shalfedges_last_  = SE[se+6*i+5];
	  vh->sfaces_begin_ = SF[sf+2*i];
	  vh->sfaces_last_  = SF[sf+2*i+1];
	  vh->shalfloop_ = shalfloops_end();
	  hx = i % 2 ? -1 : 1;
	  hy = i % 4 > 1 ? -1 : 1;
	  hz = i > 3 ? -1 : 1;
	  vh->point_at_center_ = Infi_box::create_extended_point(hx, hy, hz);
	  vh->mark_ = 1;		
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

	for(i = 0; i < 24; i++) {
	  Halfedge_handle eh = E[e+i];
	  eh->twin_ = E[e+twinIdx[i]];
	  eh->center_vertex_ = V[v+(i/3)];
	  eh->out_sedge_ = SE[se+(i/3*6)+seOff[i%3]];
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
	  eh->point_on_surface_ = Sphere_point(hx,hy,hz);
	  eh->mark_ = 1;
	}

	int bnd[12] = {19, 18, 43, 42, 35, 34,
		       47, 46, 39, 38, 45, 44};
	for(i = 0; i < 12; i++) {
	  Halffacet_handle fh = F[f+i];
	  fh->twin_ = F[f+(i/2*2)+((i+1)%2)];
	  fh->boundary_entry_objects_.push_back(SE[se+bnd[i]]);
	  fh->volume_ = C[((i%4) == 1 || (i%4 == 2)) ? 1 : 0];
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
	  fh->mark_ = 1;
	}

	C[0]->shell_entry_objects_.push_back(SF[sf]);
	C[0]->mark_ = 0;
	C[1]->shell_entry_objects_.push_front(SF[sf+1]);

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

	for(i = 0; i < 48; i++) {
	  SHalfedge_handle seh = SE[se+i];
	  
	  seh->twin_ = SE[se+(i/2*2)+((i+1)%2)];
	  seh->sprev_ = SE[se+sprevOff[i%6]+(i/6*6)];
	  seh->snext_ = SE[se+snextOff[i%6]+(i/6*6)];
	  seh->source_ = E[e+((i+1)%6)/2+(i/6)*3];
	  seh->incident_sface_ = SF[sf+(i%2)+(i/6)*2];
	  seh->prev_ = SE[se+prevIdx[i]];
	  seh->next_ = SE[se+nextIdx[i]];
	  seh->incident_facet_ = F[f+factIdx[i]];
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
	  seh->circle_ = Sphere_circle(Plane_3(RT(hx),RT(hy),RT(hz),RT(0)));
	  seh->mark_ = 1;
	}
	
	int volIdx[8] = {0,1,1,0,1,0,0,1};
		
	for(i = 0; i < 16; i++) {
	  SFace_handle sfh = SF[sf+i];
	  sfh->center_vertex_ = V[v+(i/2)];
	  sfh->boundary_entry_objects_.push_back(SE[se+(i/2*6)+(i%2)]);
	  int cIdx = i%2 ? 1-volIdx[i/2] : volIdx[i/2];
	  sfh->incident_volume_ = C[cIdx];
	  sfh->mark_ = cIdx ? C[1]->mark_ : 0;
	}
	
      }

      for(i = 0; i<v+addInfiBox*8; i++) {
	SM_point_locator PL(V[i]);
	PL.init_marks_of_halfspheres(); 
      }

      delete[] V;
      delete[] E;
      delete[] F;
      delete[] C;
      delete[] SE;
      delete[] SL;
      delete[] SF;

      return OK;
    }

protected:
  void pointer_update(const Self& D);
  static Object_iterator              undef_;
  Generic_handle_map<Object_iterator> boundary_item_;
  Generic_handle_map<Object_iterator> sm_boundary_item_;

  Vertex_list    vertices_;
  Halfedge_list  halfedges_;
  Halffacet_list halffacets_;
  Volume_list    volumes_;
  SHalfedge_list shalfedges_;
  SHalfloop_list shalfloops_;
  SFace_list     sfaces_;

  SNC_io_parser<SNC_structure> *IO;
}; // SNC_structure


template <typename Items>
void SNC_structure<Items>::
pointer_update(const SNC_structure<Items>& D)
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

  CGAL_nef3_forall_vertices(v,*this) {
    // Local Graph update: (SVertices are postponed/updated as Edges)
    v->sncp_ = this;
    v->svertices_begin_ = EM[v->svertices_begin_];
    v->svertices_last_ = EM[v->svertices_last_];
    v->shalfedges_begin_ = SEM[v->shalfedges_begin_];
    v->shalfedges_last_ = SEM[v->shalfedges_last_];
    v->sfaces_begin_ = SFM[v->sfaces_begin_];
    v->sfaces_last_ = SFM[v->sfaces_last_];
    v->shalfloop_ = SLM[v->shalfloop_];
  }
  // Halfedge update:
  CGAL_nef3_forall_halfedges(e,*this) {
    e->center_vertex_ = VM[e->center_vertex_];
    e->twin_ = EM[e->twin_];
    e->out_sedge_ = SEM[e->out_sedge_];
    e->incident_sface_ = SFM[e->incident_sface_];
  }
  // Halffacet update
  CGAL_nef3_forall_halffacets(f,*this) {
    f->twin_ = FM[f->twin_];
    f->volume_ = CM[f->volume_];
    Halffacet_cycle_iterator ftc;
    for(ftc = f->boundary_entry_objects_.begin(); 
        ftc !=  f->boundary_entry_objects_.end(); ++ftc) {
      if ( assign( se, ftc) ) 
      { *ftc = SObject_handle(SEM[se]); store_boundary_item(se,ftc); }
      else if ( assign( sl, ftc) ) 
      { *ftc = SObject_handle(SLM[sl]); store_boundary_item(sl,ftc); }
      else CGAL_nef3_assertion_msg(0,"damn wrong boundary item in facet.");
    }
  }
  for ( fc = D.halffacets_begin(), f = halffacets_begin();
	fc != D.halffacets_end(); ++fc, ++f) {
    /* It is possible that the is_twin() property differs for equivalent 
       facets on both SNC structures.  So, we need to store the correct
       selection mark in the correct (non-twin) facet of a halffacet pair. */
    CGAL_nef3_assertion_code( if( fc->is_twin() == f->is_twin())
			 CGAL_nef3_assertion( fc->mark_ == f->mark_));
    if( !f->is_twin() && fc->is_twin()) f->mark_ = f->twin_->mark_;
  }
  // Volume update
  CGAL_nef3_forall_volumes(c,*this) {
    Shell_entry_iterator sei;
    CGAL_nef3_forall_shells_of(sei,c) {
      sf = sei; // conversion from generic iterator to sface const handle
      *sei = Object_handle(SFM[sf]); 
      store_boundary_item(sf,sei); 
    }
  }

  CGAL_nef3_forall_shalfedges(se,*this) {
    se->source_ = EM[se->source_];
    se->sprev_ = SEM[se->sprev_]; se->snext_ = SEM[se->snext_];
    se->incident_sface_ = SFM[se->incident_sface_];
    se->twin_ = SEM[se->twin_];
    se->prev_ = SEM[se->prev_]; se->next_ = SEM[se->next_];
    se->incident_facet_ = FM[se->incident_facet_];
  }
  for ( sec = D.shalfedges_begin(), se = shalfedges_begin();
	sec != D.shalfedges_end(); ++sec, ++se) {
    /* It is possible that the is_twin() property differs for equivalent 
       sedges on both SNC structures.  So, we need to store the correct
       selection mark in the correct (non-twin) facet of a shalfedge pair. */
    CGAL_nef3_assertion_code( if( sec->is_twin() == se->is_twin())
			 CGAL_nef3_assertion( sec->mark_ == se->mark_));
    if( !se->is_twin() && sec->is_twin()) se->mark_ = se->twin_->mark_;
  }

  CGAL_nef3_forall_shalfloops(sl,*this) {
    sl->twin_ = SLM[sl->twin_];
    sl->incident_sface_ = SFM[sl->incident_sface_];
    sl->incident_facet_ = FM[sl->incident_facet_];
  }
  for ( slc = D.shalfloops_begin(), sl = shalfloops_begin();
	slc != D.shalfloops_end(); ++slc, ++sl) {
    /* It is possible that the is_twin() property differs for equivalent 
       sloops on both SNC structures.  So, we need to store the correct
       selection mark in the correct (non-twin) facet of a shalfloop pair. */
    CGAL_nef3_assertion_code( if( slc->is_twin() == sl->is_twin())
			 CGAL_nef3_assertion( slc->mark_ == sl->mark_));
    if( !sl->is_twin() && slc->is_twin()) sl->mark_ = sl->twin_->mark_;
  }

  CGAL_nef3_forall_sfaces(sf,*this) {
    sf->center_vertex_ = VM[sf->center_vertex_];
    sf->incident_volume_ = CM[sf->incident_volume_];
    SFace_cycle_iterator sfc;
    for(sfc = sf->sface_cycles_begin(); 
        sfc != sf->sface_cycles_end(); ++sfc) {
      SVertex_handle sv;
      if ( assign(sv,sfc) ) 
      { *sfc = SObject_handle(EM[sv]); store_sm_boundary_item(sv,sfc); }
      else if ( assign(se,sfc) ) 
      { *sfc = SObject_handle(SEM[se]); store_sm_boundary_item(se,sfc); }
      else if ( assign(sl,sfc) ) 
      { *sfc = SObject_handle(SLM[sl]); store_sm_boundary_item(sl,sfc); }
      else CGAL_nef3_assertion_msg(0,"damn wrong boundary item in sface.");
    }
  }
}

template <typename Items> 
typename SNC_structure<Items>::Object_iterator
SNC_structure<Items>::undef_;


CGAL_END_NAMESPACE
#endif // CGAL_SNC_STRUCTURE_H
