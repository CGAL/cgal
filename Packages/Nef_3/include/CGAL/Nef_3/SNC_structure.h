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
// file          : include/CGAL/Nef_3/SNC_structure.h
// package       : Nef_3
// chapter       : 3D-Nef Polyhedra
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Michael Seel    <seel@mpi-sb.mpg.de>
//                 Miguel Granados <granados@mpi-sb.mpg.de>
//                 Susan Hert      <hert@mpi-sb.mpg.de>
//                 Lutz Kettner    <kettner@mpi-sb.mpg.de>
// maintainer    : Susan Hert      <hert@mpi-sb.mpg.de>
//                 Lutz Kettner    <kettner@mpi-sb.mpg.de>
// coordinator   : MPI Saarbruecken
//
// SNC_structure.h         the SNC class, global structure, allocation and
// ============================================================================
#ifndef CGAL_SNC_STRUCTURE_H
#define CGAL_SNC_STRUCTURE_H

#include <CGAL/basic.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/Nef_S2/Generic_handle_map.h>
#include <CGAL/Nef_S2/Sphere_geometry.h>
#include <CGAL/Nef_3/SNC_iteration.h>
#include <CGAL/Nef_3/SNC_items.h>
#include <CGAL/Nef_3/nef3_assertions.h>
#include <CGAL/Nef_2/iterator_tools.h>
#include <CGAL/Union_find.h>
#include <list>

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
  CGAL_nef3_assertion( hash[o1] != 0 && hash[o2] != 0);
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

  typedef typename Items::Kernel        Kernel;
  typedef typename Kernel::FT           FT;
  typedef typename Kernel::RT           RT;
  typedef typename Items::Sphere_kernel Sphere_kernel;
  

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
    shalfedges_(), shalfloops_(), sfaces_() {}
  ~SNC_structure() { clear(); }

  SNC_structure(const Self& D) : 
    boundary_item_(undef_), sm_boundary_item_(undef_),
    vertices_(D.vertices_), halfedges_(D.halfedges_), 
    halffacets_(D.halffacets_), volumes_(D.volumes_),
    shalfedges_(D.shalfedges_), shalfloops_(D.shalfloops_), sfaces_(D.sfaces_)
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
  bool is_boundary_object(H h) 
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
  
  SNC_io_parser<SNC_structure> *IO;

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
	  TRACEN("UNION of "<<IO->index(fu)<<" & "<<IO->index(ftu));
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
	    SD.delete_vertex_only(src);
	  if( SD.is_isolated(tgt))
	    SD.delete_vertex_only(tgt);
	  /* TO VERIFY: can both svertices be isolated at the same time? */
      }
      }
      else if( assign(l, fc)) {
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
       graph for two antipodal vertices possible connected by a bundle of
       sedges. */
    bool is_part = false;
    SM_decorator SD(v);
    if( !SD.has_loop()) {
      TRACE(SNC_decorator(*this).point(v)<<" is in edge interior? ");
      SVertex_iterator sv(SD.svertices_begin());
      SVertex_handle p1(sv++), p2(sv++); // TODO: is it dangerous?
      if( sv == SD.svertices_end()) {
	TRACE("has two svertices? ");
	Sphere_point sp1(SD.point(p1)), sp2(SD.point(p2));
	if( sp1 == sp2.antipode()) {
	  TRACE("are they antipode? ");
	  SHalfedge_iterator se;
	  CGAL_nef3_forall_sedges_of( se, v) {
	    if( (SD.source(se) != p1 && SD.target(se) != p2) &&
		(SD.source(se) != p2 && SD.target(se) != p1))
	      break;
	  }
	  is_part = (se == SD.shalfedges_end()) ? true: false;
	  TRACE("all sedges conect them? ");
	}
      }
    }
    TRACEN((is_part?"yes":"no"));
    return is_part;
  }
  
  void simplify() {
    TRACEN(">>> simplifying");
    SNC_decorator D(*this);
    SNC_io_parser<SNC_structure> IO_parser(std::cerr, *this);
    IO = &IO_parser;
    
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
	TRACEN("UNION of "<<IO->index(c1)<<" & "<<IO->index(c2));
      }
      f = f_next;
    }

    Halffacet_iterator hf;
    CGAL_nef3_forall_halffacets( hf, *this) {
      hash_facet[hf] = uf_facet.make_set(hf);
      reset_object_list(hf->boundary_entry_objects_);

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
	if(D.mark(e) == D.mark(D.volume(D.source(e)->sfaces_begin()))) {
	  TRACEN("removing pair "<<IO->index(e)<<' '<<IO->index(D.twin(e)));
	  delete_halfedge_pair(e);
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
	    // TRACEN("AFTER "<<PFC(e1)); // e1 not valid after sloop creation
	  }
	}
      }

      e = e_next;
    }
 

    /* 
     * Vertices simplification
     */
    Vertex_iterator v = (*this).vertices_begin();
    while( v != (*this).vertices_end()) {
      SM_decorator SD(v);
      Vertex_iterator v_next(v);
      v_next++;

      CGAL_nef3_assertion( SD.sfaces_begin() != SFace_handle());
      if( is_part_of_volume(v)) {
	TRACEN("mark("<<IO->index(v)<<")="<<D.mark(v)<<", "<<
	       "mark("<<IO->index(D.volume(SD.sfaces_begin()))<<")="<<
	       D.mark(D.volume(SD.sfaces_begin())));
	if(D.mark(v) == D.mark(D.volume(SD.sfaces_begin()))) {
	  TRACEN("removing isolated vertex "<<IO->index(v));
	  delete_vertex(v);
	}
      }
      else if( is_part_of_facet(v)) {
	CGAL_nef3_assertion( D.facet(SD.shalfloop()) != Halffacet_handle());
	if( D.mark(v) == D.mark(D.facet(SD.shalfloop()))) {
	  TRACEN("removing "<<IO->index(v)<<
		 " on facet "<<IO->index(D.facet(SD.shalfloop())));
	  delete_vertex(v);
	}
      }
      else if( is_part_of_edge(v)) {
	SVertex_iterator sv(SD.svertices_begin());
	Halfedge_handle e1(sv++), e2(sv++);
	CGAL_nef3_assertion( sv == SD.svertices_end());
	if( D.mark(e1) == D.mark(v) && D.mark(v) == D.mark(e2)) {
	  TRACEN("merging "<<IO->index(e1)<<" & "<<IO->index(e2)<<
		 " in "<<IO->index(v));
	  merge_halfedge_pairs( e1, e2);
	}
      }
      v = v_next;
    }

    purge_no_find_objects(hash_volume, hash_facet, hash_sface, uf_volume, 
                          uf_facet, uf_sface);
    create_boundary_links_forall_sfaces( hash_sface, uf_sface);
    create_boundary_links_forall_facets( hash_facet, uf_facet);
    create_boundary_links_forall_volumes( hash_volume, uf_volume);

    TRACEN(">>> simplifying done");
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
     CGAL_nef3_forall_sfaces( sf, *this) {
       if( uf_sface.find(hash_sface[sf]) != hash_sface[sf]) {
	 TRACEN("no find object "<<IO->index(sf));
	 SM_decorator SD(D.vertex(sf));
	 SD.delete_face_only(sf);
       }
     }
     Halffacet_iterator f;
     CGAL_nef3_forall_halffacets( f, *this) {
       if( uf_facet.find(hash_facet[f]) != hash_facet[f]) {
	 TRACEN("no find object "<<IO->index(f));
	 delete_halffacet_pair(f);
       }
     }
     Volume_iterator c;
     CGAL_nef3_forall_volumes( c, *this) {
       if( uf_volume.find(hash_volume[c]) != hash_volume[c]) {
	 TRACEN("no find object "<<IO->index(c));
	 delete_volume(c);
       }
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
	SFace_handle sf = *(uf.find(hash[D.source(sv)->sfaces_begin()])); 
	CGAL_nef3_assertion( sf != SFace_handle());
	SD.set_face( sv, sf);
	SD.store_boundary_object( sv, sf);
      }
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
	Point_3 p(D.point(D.vertex(c)));
	if( lexicographically_xyz_smaller( p, D.point(D.vertex(u_min))))
	  u_min = c;
	linked[c] = true;
      }
      /* store the edge use at the lexicographicaly minimum facet vertex, as
	 a cycle entry of f.  The outermost cycle is stored at first
	 on the facet's cycles list. */
      SObject_list f_entries(f->boundary_entry_objects_);
      if( is_empty_range( f_entries.begin(),f_entries.end()))
	D.store_boundary_object( u_min, f);
      else {
	SHalfedge_handle f_sedge;
	CGAL_nef3_assertion( assign( f_sedge, 
				     f->boundary_entry_objects_.front()));
	assign( f_sedge, f->boundary_entry_objects_.front());
	Point_3 p(D.point(D.vertex(f_sedge)));
	if( lexicographically_xyz_smaller( D.point(D.vertex(u_min)), p))
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
    Unique_hash_map< Volume_handle, bool> linked(false);
    SNC_decorator D(*this);
    SFace_iterator sf;
    CGAL_nef3_forall_sfaces(sf, *this) {
      Volume_handle c = *(uf.find(hash[D.volume(sf)]));
      if( linked[c])
	continue;
      linked[c] = true;
      Volume_setter setter(D, c);
      D.visit_shell_objects( sf, setter );
      D.store_boundary_object( sf, c);
    }
  }
  
    // Returns the bounding box of the finite vertices of the polyhedron.
    // Returns $[-1,+1]^3$ as bounding box if no finite vertex exists.
    Bbox_3  bounded_bbox() {
        SNC_decorator deco(*this);
        Vertex_iterator vi = vertices_begin();
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

