// Copyright (c) 1997-2000  Max-Planck-Institute Saarbruecken (Germany).
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
// Author(s)     : Peter Hachenberger <hachenberger@mpi-sb.mpg.de>

#ifndef CGAL_SNC_LIST_H
#define CGAL_SNC_LIST_H

#include <CGAL/license/Nef_3.h>


#include <CGAL/Nef_S2/SM_list.h>
#include <CGAL/Unique_hash_map.h>

namespace CGAL {

template < class Sphere_map>
class SNC_in_place_list_sm
    : public Sphere_map, 
      public In_place_list_base<SNC_in_place_list_sm<Sphere_map> > {
public:
    typedef SNC_in_place_list_sm<Sphere_map> Self;
    //    typedef typename Vertex::Vertex_handle       Vertex_handle;
    //    typedef typename Vertex::Vertex_const_handle Vertex_const_handle;
    SNC_in_place_list_sm() {}
    SNC_in_place_list_sm(const Sphere_map& sm)   // down cast
        : Sphere_map(sm) {}
    Self& operator=( const Self& sm) {
        // This self written assignment avoids that assigning vertices will
        // overwrite the list linking of the target vertex.
        *((Sphere_map*)this) = ((const Sphere_map&)sm);
        return *this;
    }
};

template < class Halffacet>
class SNC_in_place_list_halffacet
    : public Halffacet, 
      public In_place_list_base<SNC_in_place_list_halffacet<Halffacet> > {
public:
    typedef SNC_in_place_list_halffacet<Halffacet> Self;
    //    typedef typename Halffacet::Halffacet_handle       Halffacet_handle;
    //    typedef typename Halffacet::Halffacet_const_handle Halffacet_const_handle;
    SNC_in_place_list_halffacet() {}
    SNC_in_place_list_halffacet(const Halffacet& v)   // down cast
        : Halffacet(v) {}
    Self& operator=( const Self& v) {
        // This self written assignment avoids that assigning vertices will
        // overwrite the list linking of the target vertex.
        *((Halffacet*)this) = ((const Halffacet&)v);
        return *this;
    }
};

template < class Volume>
class SNC_in_place_list_volume
    : public Volume, 
      public In_place_list_base<SNC_in_place_list_volume<Volume> > {
public:
    typedef SNC_in_place_list_volume<Volume> Self;
    //    typedef typename Volume::Volume_handle       Volume_handle;
    //    typedef typename Volume::Volume_const_handle Volume_const_handle;
    SNC_in_place_list_volume() {}
    SNC_in_place_list_volume(const Volume& v)   // down cast
        : Volume(v) {}
    Self& operator=( const Self& v) {
        // This self written assignment avoids that assigning vertices will
        // overwrite the list linking of the target vertex.
        *((Volume*)this) = ((const Volume&)v);
        return *this;
    }
};

template < class SHalfloop>
class SNC_in_place_list_shalfloop
    : public SHalfloop, 
      public In_place_list_base<SNC_in_place_list_shalfloop<SHalfloop> > {
public:
    typedef SNC_in_place_list_shalfloop<SHalfloop> Self;
    //    typedef typename SHalfloop::SHalfloop_handle       SHalfloop_handle;
    //    typedef typename SHalfloop::SHalfloop_const_handle SHalfloop_const_handle;
    SNC_in_place_list_shalfloop() {}
    SNC_in_place_list_shalfloop(const SHalfloop& v)   // down cast
        : SHalfloop(v) {}
    Self& operator=( const Self& v) {
        // This self written assignment avoids that assigning vertices will
        // overwrite the list linking of the target vertex.
        *((SHalfloop*)this) = ((const SHalfloop&)v);
        return *this;
    }
};

/*
template <typename K, typename I> class SNC_sphere_map;
template <typename S> class SM_decorator;

template<typename Kernel_,typename Items_>
class SNC_list : public SM_list<CGAL::Sphere_geometry<Kernel_>,Items_> {

 public:
  typedef Kernel_                               Kernel;
  typedef Items_                                Items;
  typedef CGAL::Sphere_geometry<Kernel>         Sphere_kernel;
  typedef CGAL::SNC_list<Kernel,Items>          Self;
  typedef CGAL::SM_list<Sphere_kernel,Items>    Base;
  typedef CGAL::SNC_sphere_map<Kernel, Items>   Sphere_map;
  typedef CGAL::SM_decorator<Sphere_map>        SM_decorator;

  typedef typename Base::Mark                   Mark;
  typedef typename Base::Size_type              Size_type;
  
  typedef typename Kernel::Point_3                 Point_3;
  typedef typename Kernel::Plane_3                 Plane_3;
  typedef typename Kernel::Vector_3                Vector_3;
  typedef typename Kernel::Direction_3             Direction_3;
  typedef typename Kernel::Segment_3               Segment_3;
  typedef typename Kernel::Line_3                  Line_3;
  typedef typename Kernel::Ray_3                   Ray_3;
  typedef typename Kernel::Triangle_3              Triangle_3;
  typedef typename Kernel::Aff_transformation_3    Aff_transformation_3;

  typedef typename Sphere_kernel::Sphere_point     Sphere_point;
  typedef typename Sphere_kernel::Sphere_segment   Sphere_segment;
  typedef typename Sphere_kernel::Sphere_circle    Sphere_circle;
  typedef typename Sphere_kernel::Sphere_direction Sphere_direction;

  typedef Sphere_map                                        Vertex_base;
  typedef SNC_in_place_list_sm<Vertex_base>                 Vertex; 
  typedef CGAL::In_place_list<Vertex,false>                 Vertex_list;
  typedef CGAL_ALLOCATOR(Vertex)                            Vertex_alloc;
  typedef typename Vertex_list::iterator                    Vertex_handle;
  typedef typename Vertex_list::const_iterator              Vertex_const_handle;
  typedef typename Vertex_list::iterator                    Vertex_iterator;
  typedef typename Vertex_list::const_iterator              Vertex_const_iterator;

  typedef typename Items::template Halffacet<Self>          Halffacet_base;
  typedef SNC_in_place_list_halffacet<Halffacet_base>       Halffacet;
  typedef CGAL::In_place_list<Halffacet,false>              Halffacet_list;
  typedef CGAL_ALLOCATOR(Halffacet)                         Halffacet_alloc;
  typedef typename Halffacet_list::iterator                 Halffacet_handle;
  typedef typename Halffacet_list::const_iterator           Halffacet_const_handle;
  typedef typename Halffacet_list::iterator                 Halffacet_iterator;
  typedef typename Halffacet_list::const_iterator           Halffacet_const_iterator;

  typedef typename Items::template Volume<Self>             Volume_base;
  typedef SNC_in_place_list_volume<Volume_base>             Volume;
  typedef CGAL::In_place_list<Volume,false>                 Volume_list;
  typedef CGAL_ALLOCATOR(Volume)                            Volume_alloc;
  typedef typename Volume_list::iterator                    Volume_handle;
  typedef typename Volume_list::const_iterator              Volume_const_handle;
  typedef typename Volume_list::iterator                    Volume_iterator;
  typedef typename Volume_list::const_iterator              Volume_const_iterator;

  typedef typename Items::template SVertex<Self>            SVertex_base;
  typedef SNC_in_place_list_svertex<SVertex_base>           SVertex;
  typedef CGAL::In_place_list<SVertex,false>                SVertex_list;
  typedef CGAL_ALLOCATOR(SVertex)                           SVertex_alloc;
  typedef typename SVertex_list::iterator                   SVertex_handle;
  typedef typename SVertex_list::const_iterator             SVertex_const_handle;
  typedef typename SVertex_list::iterator                   SVertex_iterator;
  typedef typename SVertex_list::const_iterator             SVertex_const_iterator;

  typedef typename Items::template SVertex<Self>            Halfedge_base;
  typedef SNC_in_place_list_svertex<SVertex_base>           Halfedge;
  typedef CGAL::In_place_list<SVertex,false>                Halfedge_list;
  typedef CGAL_ALLOCATOR(SVertex)                           Halfedge_alloc;
  typedef typename SVertex_list::iterator                   Halfedge_handle;
  typedef typename SVertex_list::const_iterator             Halfedge_const_handle;
  typedef typename SVertex_list::iterator                   Halfedge_iterator;
  typedef typename SVertex_list::const_iterator             Halfedge_const_iterator;

  typedef typename Items::template SHalfedge<Self>          SHalfedge_base;
  typedef SNC_in_place_list_shalfedge<SHalfedge_base>       SHalfedge;
  typedef CGAL::In_place_list<SHalfedge,false>              SHalfedge_list;
  typedef CGAL_ALLOCATOR(SHalfedge)                         SHalfedge_alloc;
  typedef typename SHalfedge_list::iterator                 SHalfedge_handle;
  typedef typename SHalfedge_list::const_iterator           SHalfedge_const_handle;
  typedef typename SHalfedge_list::iterator                 SHalfedge_iterator;
  typedef typename SHalfedge_list::const_iterator           SHalfedge_const_iterator;

  typedef typename Items::template SHalfloop<Self>          SHalfloop_base;
  typedef SNC_in_place_list_shalfloop<SHalfloop_base>       SHalfloop;
  typedef CGAL::In_place_list<SHalfloop,false>              SHalfloop_list;
  typedef CGAL_ALLOCATOR(SHalfloop)                         SHalfloop_alloc;
  typedef typename SHalfloop_list::iterator                 SHalfloop_handle;
  typedef typename SHalfloop_list::const_iterator           SHalfloop_const_handle;
  typedef typename SHalfloop_list::iterator                 SHalfloop_iterator;
  typedef typename SHalfloop_list::const_iterator           SHalfloop_const_iterator;

  typedef typename Items::template SFace<Self>              SFace_base;
  typedef SNC_in_place_list_sface<SFace_base>               SFace;
  typedef CGAL::In_place_list<SFace,false>                  SFace_list;
  typedef CGAL_ALLOCATOR(SFace)                             SFace_alloc;
  typedef typename SFace_list::iterator                     SFace_handle;
  typedef typename SFace_list::const_iterator               SFace_const_handle;
  typedef typename SFace_list::iterator                     SFace_iterator;
  typedef typename SFace_list::const_iterator               SFace_const_iterator;

  typedef typename Base::SVertex                   Halfedge;
  typedef typename Base::SVertex_handle            Halfedge_handle;
  typedef typename Base::SVertex_iterator          Halfedge_iterator;
  typedef typename Base::SVertex_const_handle      Halfedge_const_handle;
  typedef typename Base::SVertex_const_iterator    Halfedge_const_iterator;
  typedef typename Base::SVertex                   SVertex;
  typedef typename Base::SVertex_handle            SVertex_handle;
  typedef typename Base::SVertex_iterator          SVertex_iterator;
  typedef typename Base::SVertex_const_handle      SVertex_const_handle;
  typedef typename Base::SVertex_const_iterator    SVertex_const_iterator;
  typedef typename Base::SHalfedge                 SHalfedge;
  typedef typename Base::SHalfedge_handle          SHalfedge_handle;
  typedef typename Base::SHalfedge_iterator        SHalfedge_iterator;
  typedef typename Base::SHalfedge_const_handle    SHalfedge_const_handle;
  typedef typename Base::SHalfedge_const_iterator  SHalfedge_const_iterator;
  typedef typename Base::SFace                     SFace;
  typedef typename Base::SFace_handle              SFace_handle;
  typedef typename Base::SFace_iterator            SFace_iterator;
  typedef typename Base::SFace_const_handle        SFace_const_handle;
  typedef typename Base::SFace_const_iterator      SFace_const_iterator;
  typedef typename Base::SHalfloop                 SHalfloop;
  typedef typename Base::SHalfloop_handle          SHalfloop_handle;
  typedef typename Base::SHalfloop_iterator        SHalfloop_iterator;
  typedef typename Base::SHalfloop_const_handle    SHalfloop_const_handle;
  typedef typename Base::SHalfloop_const_iterator  SHalfloop_const_iterator;

  typedef typename Base::Object_handle             Object_handle;
  typedef typename Base::Object_list               Object_list;
  typedef typename Base::Object_iterator           Object_iterator;
  typedef typename Base::Object_const_iterator     Object_const_iterator;

  typedef typename Base::SHalfedge_around_svertex_circulator
                         SHalfedge_around_svertex_circulator;
  typedef typename Base::SHalfedge_around_sface_circulator
                         SHalfedge_around_sface_circulator;
  typedef typename Base::SHalfedge_around_svertex_const_circulator
                         SHalfedge_around_svertex_const_circulator;
  typedef typename Base::SHalfedge_around_sface_const_circulator
                         SHalfedge_around_sface_const_circulator;

  typedef typename Base::SFace_cycle_iterator      
                         SFace_cycle_iterator;
  typedef typename Base::SFace_cycle_const_iterator      
                         SFace_cycle_const_iterator;

  typedef CircFromIt<SHalfedge_const_iterator, 
          move_shalfedge_around_facet<SHalfedge_const_iterator> > 
          SHalfedge_around_facet_const_circulator;

  typedef CircFromIt<SHalfedge_iterator, 
          move_shalfedge_around_facet<SHalfedge_iterator> > 
          SHalfedge_around_facet_circulator;

  class Halffacet_cycle_iterator : public Object_iterator 
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
    { CGAL_error_msg("not impl."); }
  };

  class Halffacet_cycle_const_iterator : public Object_const_iterator 
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
    { CGAL_error_msg("not impl."); }
  };

  class Shell_entry_iterator : public Object_iterator 
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
      CGAL_assertion( CGAL::assign(f,Ibase::operator*()) );
      CGAL::assign(f,Ibase::operator*()); 
      return SFace_const_handle(f); }

    operator Object_handle() const { return Ibase::operator*(); }
    Object_handle& operator*() const { return Ibase::operator*(); }
    Object_handle  operator->() const 
    { CGAL_nef_assertion_msg(0,"not impl."); }
  };

 private:
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

 public:
  SNC_list() : 
    boundary_item_(undef_),
    vertices_(), halffacets_(), volumes_() {}

  ~SNC_list() { clear(); }

  SNC_list(const Self& D) : 
    boundary_item_(undef_),
    vertices_(D.vertices_), 
    halffacets_(D.halffacets_),
    volumes_(D.volumes_) { 
    pointer_update(D); 
  }

  Self& operator=(const Self& D) {
    if ( this == &D )
      return *this;
    clear();
    boundary_item_.clear(undef_);
    vertices_ = D.vertices_;
    halffacets_ = D.halffacets_;
    volumes_ = D.volumes_;
    pointer_update(D);
    return *this;
  }

  void clear() { 
    Base::clear();
    boundary_item_.clear();
    vertices_.destroy();
    halffacets_.destroy();
    volumes_.destroy();
  }

  Vertex_const_iterator vertices_begin() const      {return vertices_.begin();}
  Vertex_const_iterator vertices_end() const        {return vertices_.end();}
  Halfedge_const_iterator halfedges_begin() const   {return svertices_.begin();}
  Halfedge_const_iterator halfedges_end() const     {return svertices_.end();}
  Halffacet_const_iterator halffacets_begin() const {return halffacets_.begin();}
  Halffacet_const_iterator halffacets_end() const   {return halffacets_.end();}
  Volume_const_iterator   volumes_begin() const     {return volumes_.begin();}
  Volume_const_iterator   volumes_end() const       {return volumes_.end();}

  Vertex_iterator    vertices_begin()   { return vertices_.begin();}
  Vertex_iterator    vertices_end()     { return vertices_.end();}
  Halfedge_iterator  halfedges_begin()  { return halfedges_.begin();}
  Halfedge_iterator  halfedges_end()    { return halfedges_.end();}
  Halffacet_iterator halffacets_begin() { return halffacets_.begin();}
  Halffacet_iterator halffacets_end()   { return halffacets_.end();}
  Volume_iterator    volumes_begin()    { return volumes_.begin();}
  Volume_iterator    volumes_end()      { return volumes_.end();}

  Size_type number_of_vertices() const   { return vertices_.size(); }
  Size_type number_of_halfedges() const  { return svertices_.size(); }
  Size_type number_of_edges() const      { return svertices_.size()/2; }
  Size_type number_of_halffacets() const { return halffacets_.size();}
  Size_type number_of_facets() const     { return halffacets_.size()/2;}
  Size_type number_of_volumes() const    { return volumes_.size();}

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
  { CGAL_assertion(boundary_item_[h]!=undef_);
    boundary_item_[h] = undef_; }

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

  void reset_object_list(Object_list& L)
  { Object_iterator oit;
    CGAL_forall_iterators(oit,L) reset_iterator_hash(oit);
    L.clear();
  }

  Vertex_handle new_vertex(const Point_3& p = Point_3(), Mark m = Mark()) { 
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
				    Mark m = Mark()) { 
    SM_decorator D1(&*v1);
    SM_decorator D2(&*v2);
    SVertex_handle e1 = D1.new_vertex();
    SVertex_handle e2 = D2.new_vertex();
    make_twins(e1,e2);
    e1->mark() = m;
    return e1;
  }

  Halffacet_handle new_halffacet_pair(const Plane_3& h = Plane_3(), 
				      Mark m = Mark()) {
    Halffacet_handle f1 = new_halffacet_only();
    Halffacet_handle f2 = new_halffacet_only();
    f1->supporting_plane_ = h; f2->supporting_plane_ = h.opposite();
    make_twins(f1,f2);
    f1->mark() = f2->mark() = m;
    return f1;
  }

  Volume_handle new_volume(Mark m = Mark()) {
    Volume_handle vh = new_volume_only();
    vh->mark() = m;
    return vh;
  }

  void delete_vertex(Vertex_handle v)
    CGAL_NEF_TRACEN("~ deleting vertex "<<&*v<<" from "<<&*this);
    v->clear(true); 
    delete_vertex_only(v);
    CGAL_NEF_TRACEN("~~ vertex deleted"<<&*v);
  }

  void delete_halfedge_pair(Halfedge_handle e)
    CGAL_NEF_TRACEN("~ deleting halfedges pair "<<&*e<<", "<<&*(e->twin())<<
	   " from "<<&*this);
    Halfedge_handle et = e->twin();
    SM_decorator D1(&*e->center_vertex()), D2(&*et->center_vertex());
    D1.delete_vertex(e);
    D2.delete_vertex(et);
  }

  void delete_halffacet_pair(Halffacet_handle f)
    CGAL_NEF_TRACEN("~ deleting halffacets pair "<<&*f<<", "<<&*(f->twin())<<
	   " from "<<&*this);
    reset_object_list(f->boundary_entry_objects());
    reset_object_list(f->twin()->boundary_entry_objects());
    delete_halffacet_only(f->twin());
    delete_halffacet_only(f);
  }

  void delete_volume(Volume_handle c)
    CGAL_NEF_TRACEN("~ deleting volume "<<&*c<<" from "<<&*this);
    reset_object_list(c->shell_entry_objects());
    delete_volume_only(c);
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

 protected:
  Generic_handle_map<Object_iterator> boundary_item_;

  Vertex_list    vertices_;
  Halffacet_list halffacets_;
  Volume_list    volumes_;

  void pointer_update(const Self& D);
};

template <typename Kernel, typename Items>
void SNC_list<Kernel,Items>::
pointer_update(const SNC_list<Kernel,Items>& D)
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
    f->volume_ = CM[f->volume_];
    Halffacet_cycle_iterator ftc;
    for(ftc = f->boundary_entry_objects().begin(); 
        ftc !=  f->boundary_entry_objects().end(); ++ftc) {
      if ( assign( se, ftc) ) 
      { *ftc = Object_handle(SEM[se]); store_boundary_item(se,ftc); }
      else if ( assign( sl, ftc) ) 
      { *ftc = Object_handle(SLM[sl]); store_boundary_item(sl,ftc); }
      else CGAL_error_msg("damn wrong boundary item in facet.");
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
    se->source() = EM[se->source()];
    se->sprev() = SEM[se->sprev()]; se->snext() = SEM[se->snext()];
    se->incident_sface() = SFM[se->incident_sface()];
    se->twin() = SEM[se->twin()];
    se->prev() = SEM[se->prev()]; se->next() = SEM[se->next()];
    se->incident_facet() = FM[se->incident_facet()];
  }

  CGAL_forall_shalfloops(sl,*this) {
    sl->twin() = SLM[sl->twin()];
    sl->incident_sface() = SFM[sl->incident_sface()];
    sl->incident_facet() = FM[sl->incident_facet()];
  }
  for ( slc = D.shalfloops_begin(), sl = shalfloops_begin();
	slc != D.shalfloops_end(); ++slc, ++sl) {
    CGAL_assertion_code( if( slc->is_twin() == sl->is_twin())
			 CGAL_assertion( slc->mark() == sl->mark()));
    if( !sl->is_twin() && slc->is_twin()) sl->mark() = sl->twin()->mark();
  }

  CGAL_forall_sfaces(sf,*this) {
    sf->center_vertex() = VM[sf->center_vertex()];
    sf->incident_volume() = CM[sf->incident_volume()];
    SFace_cycle_iterator sfc;
    for(sfc = sf->sface_cycles_begin(); 
        sfc != sf->sface_cycles_end(); ++sfc) {
      SVertex_handle sv;
      if ( assign(sv,sfc) ) 
      { *sfc = Object_handle(EM[sv]); store_sm_boundary_item(sv,sfc); }
      else if ( assign(se,sfc) ) 
      { *sfc = Object_handle(SEM[se]); store_sm_boundary_item(se,sfc); }
      else if ( assign(sl,sfc) ) 
      { *sfc = Object_handle(SLM[sl]); store_sm_boundary_item(sl,sfc); }
      else CGAL_error_msg("damn wrong boundary item in sface.");
    }
  }
}
*/
} //namespace CGAL

#endif // CGAL_SNC_LIST_H
