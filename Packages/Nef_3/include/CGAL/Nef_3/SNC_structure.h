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
#include <list>

#undef _DEBUG
#define _DEBUG 41
#include <CGAL/Nef_3/debug.h>

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

/*{\Manpage {SNC_structure}{Items}{Selective Nef Complex}{C}}*/

template <typename Items_>
class SNC_structure {
/*{\Mdefinition The extended Wuerzburg structure is the topological
structure of Nef polyhedra. It is programmed around the local
graphs of the vertices of a Nef polyhedron, which describe the
point set completely. All other concepts are either derived from
the local graph or added for the comfort of the user.}*/

public:
  /*{\Mtypes 7}*/
  typedef Items_                 Items;
  typedef SNC_structure<Items>   Self;

  friend class SNC_SM_decorator<Self>;
  friend class SNC_decorator<Self>;
  friend class SNC_io_parser<Self>;

  typedef SNC_SM_decorator<Self> SM_decorator;
  typedef SNC_decorator<Self>    SNC_decorator;

  typedef typename Items::Kernel Kernel;
  typedef typename Items::Sphere_kernel Sphere_kernel;

  typedef typename Kernel::Point_3  Point_3;
  /*{\Mtypemember embedding vertices.}*/
  typedef typename Kernel::Plane_3  Plane_3;
  /*{\Mtypemember supporting facets.}*/
  typedef typename Kernel::Vector_3 Vector_3;
  /*{\Mtypemember normal vectors.}*/
  typedef typename Kernel::Direction_3 Direction_3;
  /*{\Mtypemember normal directions.}*/
  typedef typename Kernel::Segment_3 Segment_3;
  /*{\Mtypemember segments in space.}*/
  typedef typename Kernel::Line_3 Line_3;
  /*{\Mtypemember lines in space.}*/

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
  typedef typename Vertex_list::iterator           Vertex_handle;
  typedef typename Vertex_list::const_iterator     Vertex_const_handle;
  typedef typename Vertex_list::iterator           Vertex_iterator;
  typedef typename Vertex_list::const_iterator     Vertex_const_iterator;

  typedef typename Items::template Halfedge<Self>  Halfedge;
  typedef CGAL::In_place_list<Halfedge,false>      Halfedge_list;
  typedef typename Halfedge_list::iterator         Halfedge_handle;
  typedef typename Halfedge_list::const_iterator   Halfedge_const_handle;
  typedef typename Halfedge_list::iterator         Halfedge_iterator;
  typedef typename Halfedge_list::const_iterator   Halfedge_const_iterator;

  typedef typename Items::template Halffacet<Self> Halffacet;
  typedef CGAL::In_place_list<Halffacet,false>     Halffacet_list;
  typedef typename Halffacet_list::iterator        Halffacet_handle;
  typedef typename Halffacet_list::const_iterator  Halffacet_const_handle;
  typedef typename Halffacet_list::iterator        Halffacet_iterator;
  typedef typename Halffacet_list::const_iterator  Halffacet_const_iterator;

  typedef typename Items::template Volume<Self>    Volume;
  typedef CGAL::In_place_list<Volume,false>        Volume_list;
  typedef typename Volume_list::iterator           Volume_handle;
  typedef typename Volume_list::const_iterator     Volume_const_handle;
  typedef typename Volume_list::iterator           Volume_iterator;
  typedef typename Volume_list::const_iterator     Volume_const_iterator;

  typedef typename Items::template Halfedge<Self>  SVertex;
  typedef typename Halfedge_list::iterator         SVertex_handle;
  typedef typename Halfedge_list::const_iterator   SVertex_const_handle;
  typedef typename Halfedge_list::iterator         SVertex_iterator;
  typedef typename Halfedge_list::const_iterator   SVertex_const_iterator;

  typedef typename Items::template SHalfedge<Self> SHalfedge;
  typedef CGAL::In_place_list<SHalfedge,false>     SHalfedge_list;
  typedef typename SHalfedge_list::iterator        SHalfedge_handle;
  typedef typename SHalfedge_list::const_iterator  SHalfedge_const_handle;
  typedef typename SHalfedge_list::iterator        SHalfedge_iterator;
  typedef typename SHalfedge_list::const_iterator  SHalfedge_const_iterator;

  typedef typename Items::template SHalfloop<Self> SHalfloop;
  typedef CGAL::In_place_list<SHalfloop,false>     SHalfloop_list;
  typedef typename SHalfloop_list::iterator        SHalfloop_handle;
  typedef typename SHalfloop_list::const_iterator  SHalfloop_const_handle;
  typedef typename SHalfloop_list::iterator        SHalfloop_iterator;
  typedef typename SHalfloop_list::const_iterator  SHalfloop_const_iterator;

  typedef typename Items::template SFace<Self>     SFace;
  typedef CGAL::In_place_list<SFace,false>         SFace_list;
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

  SNC_structure() : boundary_item_(undef_),
    vertices_(), halfedges_(), halffacets_(), volumes_(),
    shalfedges_(), shalfloops_(), sfaces_() {}
  ~SNC_structure() { clear(); }

  SNC_structure(const Self& D) : boundary_item_(undef_),
    vertices_(D.vertices_), halfedges_(D.halfedges_), 
    halffacets_(D.halffacets_), volumes_(D.volumes_),
    shalfedges_(D.shalfedges_), shalfloops_(D.shalfloops_), sfaces_(D.sfaces_)
  { pointer_update(D); }

  Self& operator=(const Self& D) 
  { if ( this == &D ) return *this;
    clear();
    boundary_item_.clear(undef_);
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

  void clear()
  { 
    boundary_item_.clear();
    vertices_.destroy();
    halfedges_.destroy();
    halffacets_.destroy();
    volumes_.destroy();
    shalfedges_.destroy();
    shalfloops_.destroy();
    sfaces_.destroy();
  }


  Point_3& point(Vertex_handle v) const
  { return v->point(); }
  Mark& mark(Vertex_handle v) const
  { return v->mark(); }

  Sphere_point& sphere_point(Halfedge_handle e) const
  { return e->sphere_point(); }
  Mark mark(Halfedge_handle e) const
  { return e->mark_; }
  Vertex_handle source(Halfedge_handle e) const
  { return e->center_vertex_; }
  Vertex_handle target(Halfedge_handle e) const
  { return e->twin_->center_vertex_; }
  Halfedge_handle twin(Halfedge_handle e) const
  { return e->twin_; }
  Halffacet_handle twin(Halffacet_handle f) const
  { return f->twin_; }

  Plane_3& supporting_plane(Halffacet_handle f)
  { return f->supporting_plane_; }
  Mark& mark(Halffacet_handle f)
  { return f->mark_; }
  Volume_handle volume(Halffacet_handle f)
  { return f->volume_; }

  Halffacet_cycle_iterator facet_cycles_begin(Halffacet_handle f)
  { return f->facet_cycles_begin(); }
  Halffacet_cycle_iterator facet_cycles_end(Halffacet_handle f)
  { return f->facet_cycles_end(); }


  Shell_entry_iterator shells_begin(Volume_handle c)
  { return c->shells_begin(); }
  Shell_entry_iterator shells_end(Volume_handle c)
  { return c->shells_end(); }


  template <typename H>
  bool is_boundary_object(H h) 
  { return boundary_item_[h]!=undef_; }

  template <typename H>
  Object_iterator& boundary_item(H h)
  { return boundary_item_[h]; }

  template <typename H>
  void store_boundary_item(H h, Object_iterator o)
  { boundary_item_[h] = o; }


  template <typename H>
  void undef_boundary_item(H h)
  { CGAL_nef3_assertion(boundary_item_[h]!=undef_);
    boundary_item_[h] = undef_; }

  void reset_iterator_hash(Object_iterator it)
  { SVertex_handle sv;
    SHalfedge_handle se;
    SHalfloop_handle sl;
    if ( assign(se,*it) ) { undef_boundary_item(se); return; }
    if ( assign(sl,*it) ) { undef_boundary_item(sl); return; }
    if ( assign(sv,*it) ) { undef_boundary_item(sv); return; }
  }

  void reset_object_list(Object_list& L)
  { Object_iterator oit;
    CGAL_forall_iterators(oit,L) reset_iterator_hash(oit);
    L.clear();
  }

  /*{\Moperations 2.5 3}*/

  // The constant iterators and circulators.
  Vertex_const_iterator   vertices_begin()  const { return vertices_.begin();}
  Vertex_const_iterator   vertices_end()    const { return vertices_.end();}
  Halfedge_const_iterator halfedges_begin() const { return halfedges_.begin();}
  Halfedge_const_iterator halfedges_end()   const { return halfedges_.end();}
  Halffacet_const_iterator halffacets_begin() const 
  { return halffacets_.begin();}
  Halffacet_const_iterator halffacets_end()   const 
  { return halffacets_.end();}
  Volume_const_iterator   volumes_begin()   const { return volumes_.begin();}
  Volume_const_iterator   volumes_end()     const { return volumes_.end();}

  SVertex_const_iterator svertices_begin() const { return halfedges_.begin();}
  SVertex_const_iterator svertices_end()   const { return halfedges_.end();}
  SHalfedge_const_iterator shalfedges_begin() const 
    { return shalfedges_.begin();}
  SHalfedge_const_iterator shalfedges_end() const 
    { return shalfedges_.end();}
  SHalfloop_const_iterator shalfloops_begin() const 
    { return shalfloops_.begin();}
  SHalfloop_const_iterator shalfloops_end()   const 
    { return shalfloops_.end();}
  SFace_const_iterator sfaces_begin() const { return sfaces_.begin();}
  SFace_const_iterator sfaces_end()   const { return sfaces_.end();}

  Vertex_iterator    vertices_begin()   { return vertices_.begin();}
  Vertex_iterator    vertices_end()     { return vertices_.end();}
  Halfedge_iterator  halfedges_begin()  { return halfedges_.begin();}
  Halfedge_iterator  halfedges_end()    { return halfedges_.end();}
  Halffacet_iterator halffacets_begin() { return halffacets_.begin();}
  Halffacet_iterator halffacets_end()   { return halffacets_.end();}
  Volume_iterator    volumes_begin()    { return volumes_.begin();}
  Volume_iterator    volumes_end()      { return volumes_.end();}

  SVertex_iterator   svertices_begin() { return halfedges_.begin();}
  SVertex_iterator   svertices_end()   { return halfedges_.end();}
  SHalfedge_iterator shalfedges_begin() { return shalfedges_.begin();}
  SHalfedge_iterator shalfedges_end()   { return shalfedges_.end();}
  SHalfloop_iterator shalfloops_begin() { return shalfloops_.begin();}
  SHalfloop_iterator shalfloops_end()   { return shalfloops_.end();}
  SFace_iterator     sfaces_begin()    { return sfaces_.begin();}
  SFace_iterator     sfaces_end()      { return sfaces_.end();}

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

  bool empty() const 
  { return number_of_vertices() == 0 &&
           number_of_halfedges() == 0 &&
           number_of_halffacets() == 0 &&
           number_of_volumes() == 0 &&
           number_of_shalfedges() == 0 &&
           number_of_shalfloops() == 0 &&
           number_of_sfaces() == 0;
  }

  Vertex_handle new_vertex(const Point_3& p = Point_3(), Mark m = Mark())
  /*{\Mop returns a new vertex at point |p| marked by |m|.}*/
  { 
    Vertex_handle vh = new_vertex_only();
    vh->point_at_center_ = p;
    vh->mark_ = m;
    vh->sncp_ = this;
    vh->svertices_begin_ = vh->svertices_last_ = svertices_end();
    vh->shalfedges_begin_ = vh->shalfedges_last_ = shalfedges_end();
    vh->sfaces_begin_ = vh->sfaces_last_ = sfaces_end();
    vh->shalfloop_ = shalfloops_end();
    TRACEN("new_vertex "<<&*vh);
    return vh;
  }

  template <typename H>
  void make_twins(H h1, H h2)
  { h1->twin_ = h2; h2->twin_ = h1; }

  Halfedge_handle new_halfedge_pair(Vertex_handle v1, Vertex_handle v2,
				    Mark m = Mark())
  /*{\Mop creates a new halfedge pair between the vertices $v_1$
  and $v_2$. The edge is marked by |m|.}*/
  { 
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
  marked with |m|.}*/
  {
    Halffacet_handle f1 = new_halffacet_only();
    Halffacet_handle f2 = new_halffacet_only();
    f1->supporting_plane_ = h; f2->supporting_plane_ = h.opposite();
    f1->mark_ = f2->mark_ = m;
    make_twins(f1,f2);
    return f1;
  }

  Volume_handle new_volume(Mark m = Mark())
  /*{\Mop creates a new volume marked with |m|.}*/
  { Volume_handle vh = new_volume_only();
    vh->mark_ = m;
    return vh;
  }

  void delete_vertex(Vertex_handle v)
  /*{\Mop deletes the vertex including the objects in its local graph.}*/
  { v->clear_local_graph(); vertices_.erase(v); delete &* v; }

  void delete_halfedge_pair(Halfedge_handle e)
  /*{\Mop deletes the halfedge pair of |e,twin(e)|.  Does not care about
  incident objects in the local graph of |source(e)|.}*/ 
  { Halfedge_handle et = e->twin_;
    SM_decorator D1(e->center_vertex_), D2(et->center_vertex_);
    D1.delete_vertex(e);
    D2.delete_vertex(et);
  }

  void delete_halffacet_pair(Halffacet_handle f)
  /*{\Mop deletes the halffacet pair |f,twin(f)|. Does not care about
  boundary cycle objects.}*/ 
  { reset_object_list(f->boundary_entry_objects_);
    reset_object_list(f->twin_->boundary_entry_objects_);
    halffacets_.erase(f->twin_);  delete &* (f->twin_);
    halffacets_.erase(f);         delete &* f;
  }

  void delete_volume(Volume_handle v)
  /*{\Mop deletes the volume |v|. Does not care about shell objects.}*/
  { reset_object_list(v->shell_entry_objects_);
    volumes_.erase(v); delete &* v;
  }


  Vertex_handle new_vertex_only()
  { vertices_.push_back( * new Vertex); return --vertices_end(); }

  Halfedge_handle new_halfedge_only(Halfedge_handle e)
  { return halfedges_.insert(e, * new Halfedge ); }
  Halfedge_handle new_halfedge_only()
  { halfedges_.push_back( * new Halfedge ); return --halfedges_end(); }

  Halffacet_handle new_halffacet_only()
  { halffacets_.push_back( * new Halffacet ); return --halffacets_end(); } 

  Volume_handle new_volume_only()
  { volumes_.push_back( * new Volume ); return --volumes_end(); }

  SHalfedge_handle new_shalfedge_only()
  { shalfedges_.push_back( * new SHalfedge ); return --shalfedges_end(); }

  SHalfloop_handle new_shalfloop_only()
  { shalfloops_.push_back( * new SHalfloop ); return --shalfloops_end(); }

  SFace_handle new_sface_only()
  { sfaces_.push_back( * new SFace ); return --sfaces_end(); }

  void delete_svertex_only(SVertex_handle h)
  { halfedges_.erase(h); delete &* h; }
  void delete_shalfedge_only(SHalfedge_handle h)
  { shalfedges_.erase(h); delete &* h; }
  void delete_shalfloop_only(SHalfloop_handle h)
  { shalfloops_.erase(h); delete &* h; }
  void delete_sface_only(SFace_handle h)
  { sfaces_.erase(h); delete &* h; }
  
  void simplify() {
    Halffacet_iterator f;
    CGAL_forall_facets(f, *this) { // all facets
      if( f->volume_ != 0 ) {
	CGAL_assertion( f->twin_ != 0 );
	CGAL_assertion( f->twin_->volume_ != 0 );
	Volume_handle c1 = f->volume_, c2 = f->twin_->volume_;
	if( c1->mark_ == f->mark_ == c2->mark_
	    // && f not part of bounding box
	    ) {
	  // unify(c1,c2) in a volume partition structure...
	  remove_f_including_all_edge_uses_in_its_boundary_cycles(f);
	}
      }
    }
    Halfedge_iterator e;
    CGAL_forall_edges(e, *this) { // all edges
      // if( e as an svertex is isolated in the local grpah at its source
      //     and mark(e) == mark(volume(incident_sface(e))) )
      SM_decorator D(source(e));
      if( D.is_isolated(e->twin_) 
	  // && mark(e) == mark(volume(e->incident_sface_))
	  ) {
	CGAL_assertion_code( SM_decorator Dt(e->twin_->center_vertex_) );
	CGAL_assertion( D.is_isolated(e->twin_->twin_) );
	// remove e and twin(e) from local graphs and delete edge pairs (*)
	delete_halfedge_pair(e);
      } else { 
	// if( e as an svertex has two incident sedges e1, e2, 
	//     which are part of the same great circle and
	//     mark(e1) == mark(e) == mark(e2) )
	// Precondition: !is_isolated(e->twin_)
	SHalfedge_handle e1, e2;
	e1 = D.first_out_edge(e->twin_);
	e2 = D.cyclic_adj_succ(e1);
	if( (e1 != e2 && e1 == D.cyclic_adj_succ(e2)) &&
	    (D.circle(e1) == D.circle(e2)) &&
	    (D.mark(e1) == e->mark_ == D.mark(e2)) ) {
	  //   merge e1 and e2 and remove e in sphere map (*)
	  if( D.source(e1) != e->twin_ )
	    swap( e1, e2 );
	  SVertex_handle v1 = D.source(e1), v2 = D.source(D.twin(e2));
	  TRACEN("merge_nodes "<<D.point(v1)<<D.point(v2));
	  CGAL_assertion(D.point(v1) == D.point(v2));
	  D.set_source(e1, D.source(e2));
	  D.set_source(D.twin(e2), D.source(D.twin(e1)));
	  D.delete_vertex_only(v1);
	  D.delete_vertex_only(v2);
	  //   unify(facet(e1),facet(e2)) in facet partition structure
	  //   do the same symmetrically at twin(e)
	}
      }
    }
    Vertex_iterator v;
    CGAL_forall_vertices(v, *this) { // all vertices
      // if( v has a trivial local graph consisting of one sface fs and
      //     mark(v) == mark(volume(fs)) )
      //   remove v and merge it into volume(fs)
      // if( v has a local graph consisting only of one sloop ls and
      //     mark(v) == mark(facet(ls)) )
      //   remove v and merge it into facet(ls)
      // if( v has a local graph consisting only of two svertices v1, v2 that
      //     are part of the same line and mark(v1) == mark(v) == mark(v2) )
      //   remove v and merge v1 and v2 into one edge in the global structure
    }
    // purge all sfaces, facets, and volumes that are no find objects
    // create boundary links forall sfaces
    // create boundary links forall facets
    // create boundary links forall volumes
  }
  
  void remove_f_including_all_edge_uses_in_its_boundary_cycles
    (Halffacet_handle& f) {
    Halffacet_cycle_iterator fc;
    // for all edge uses u in the boundary cycles of f
    CGAL_forall_facet_cycles_of(fc, f) {
      SHalfedge_handle e;
      if( assign(e, fc) ) {
	SHalfedge_around_facet_circulator ec(e), ee(e);
	CGAL_For_all(ec, ee) {
	  // if( u as an sedges separates two differente sfaces f1 and f2 )
	  //   unify(f1,f2) is sface partition structure...
	  //   remove u from sphere map
	}
      }
    }
  }
  
  void create_boundary_links_forall_sfaces() {
    SHalfedge_iterator e;
    CGAL_forall_sedges(e, *this) { // all sedges
      // if( e was not yet assigned to an sface cycle )
      //   assign all sedges in sface cycle of e to find(sface(e))
    }
    // forall isolated svertices vs do
    //   assign vs to find(face(vs))
  }

  void create_boundary_links_forall_facets() {
    SHalfedge_iterator e;
    CGAL_forall_edges(u, *this) { // all edge uses
      // if( u was not yet assigned to a facet cycle )
      //   assign all edge uses in facet cycle of u to find(facet(u))
    }
    SHalfloop_iterator l;
    CGAL_forall_shalfloops(l, *this) { // all sloops
      // assign l to find(facet(l))
    }
  }

  void create_boundary_links_forall_volumes() {
    SFace_iterator f;
    CGAL_forall_sfaces(f, *this) { // all sfaces
      // if( f was not yet assigned to a shell )
      //    assign all sfaces a facets in the shell of f to find(volume(f))
    }
  }

  
protected:
  void pointer_update(const Self& D);
  static Object_iterator              undef_;
  Generic_handle_map<Object_iterator> boundary_item_;

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

  CGAL_forall_vertices(v,*this) {
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
  CGAL_forall_halfedges(e,*this) {
    e->center_vertex_ = VM[e->center_vertex_];
    e->twin_ = EM[e->twin_];
    e->out_sedge_ = SEM[e->out_sedge_];
    e->incident_sface_ = SFM[e->incident_sface_];
  }
  // Halffacet update
  CGAL_forall_halffacets(f,*this) {
    f->twin_ = FM[f->twin_];
    f->volume_ = CM[f->volume_];
    Halffacet_cycle_iterator ftc;
    for(ftc = f->boundary_entry_objects_.begin(); 
        ftc !=  f->boundary_entry_objects_.end(); ++ftc) {
      if ( assign(se,SObject_handle(ftc)) ) 
      { *ftc = SObject_handle(SEM[se]); store_boundary_item(se,ftc); }
      else if ( assign(sl,SObject_handle(ftc)) ) 
      { *ftc = SObject_handle(SLM[sl]); store_boundary_item(sl,ftc); }
      else CGAL_nef3_assertion_msg(0,"damn wrong boundary item in facet.");
    }
  }
  // Volume update
  CGAL_forall_volumes(c,*this) {
    Shell_entry_iterator sei;
    CGAL_forall_shells_of(sei,c) {
      sf = sei; // conversion from generic iterator to sface const handle
      *sei = Object_handle(SFM[sf]); 
      store_boundary_item(sf,sei); 
    }
  }

  CGAL_forall_shalfedges(se,*this) {
    se->source_ = EM[se->source_];
    se->sprev_ = SEM[se->sprev_]; se->snext_ = SEM[se->snext_];
    se->incident_sface_ = SFM[se->incident_sface_];
    se->twin_ = SEM[se->twin_];
    se->prev_ = SEM[se->prev_]; se->next_ = SEM[se->next_];
    se->incident_facet_ = FM[se->incident_facet_];
  }
  CGAL_forall_shalfloops(sl,*this) {
    sl->twin_ = SLM[sl->twin_];
    sl->incident_sface_ = SFM[sl->incident_sface_];
    sl->incident_facet_ = FM[sl->incident_facet_];
  }
  CGAL_forall_sfaces(sf,*this) {
    sf->center_vertex_ = VM[sf->center_vertex_];
    sf->incident_volume_ = CM[sf->incident_volume_];
    SFace_cycle_iterator sfc;
    for(sfc = sf->sface_cycles_begin(); 
        sfc != sf->sface_cycles_end(); ++sfc) {
      SVertex_handle sv;
      if ( assign(sv,SObject_handle(sfc)) ) 
      { *sfc = SObject_handle(EM[sv]); store_boundary_item(sv,sfc); }
      else if ( assign(se,SObject_handle(sfc)) ) 
      { *sfc = SObject_handle(SEM[se]); store_boundary_item(se,sfc); }
      else if ( assign(sl,SObject_handle(sfc)) ) 
      { *sfc = SObject_handle(SLM[sl]); store_boundary_item(sl,sfc); }
      else CGAL_nef3_assertion_msg(0,"damn wrong boundary item in sface.");
    }
  }
}

template <typename Items> 
typename SNC_structure<Items>::Object_iterator
SNC_structure<Items>::undef_;


CGAL_END_NAMESPACE
#endif // CGAL_SNC_STRUCTURE_H

