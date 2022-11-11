// Copyright (c) 1997-2000  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Seel <seel@mpi-sb.mpg.de>
//                 Peter Hachenberger <hachenberger@mpi-sb.mpg.de>

#ifndef CGAL_SNC_SPHERE_MAP_H
#define CGAL_SNC_SPHERE_MAP_H

#include <CGAL/license/Nef_3.h>


#include <CGAL/basic.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/Nef_2/Object_handle.h>
#include <CGAL/Nef_S2/SM_iteration.h>
#include <CGAL/Nef_S2/Generic_handle_map.h>
#include <CGAL/Nef_2/iterator_tools.h>
#include <CGAL/Nef_3/Infimaximal_box.h>
#include <CGAL/Nef_S2/Sphere_geometry.h>
#include <list>
#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 109
#include <CGAL/Nef_2/debug.h>


namespace CGAL {

template <typename K, typename I, typename M> class SNC_structure;

template <typename HE>
struct move_shalfedge_around_svertex {
  void forward(HE& e) const { e = (e->sprev()->twin()); }
  void backward(HE& e) const { e = (e->twin()->snext()); }
};

template <typename HE>
struct move_shalfedge_around_sface {
  void forward(HE& e)  const { e = (e->snext()); }
  void backward(HE& e) const { e = (e->sprev()); }
};

template <typename Kernel_, typename Items_, typename Mark_>
class SNC_sphere_map : public Items_::template Vertex<SNC_structure<Kernel_, Items_, Mark_> >
{

 public:
  /*{\Mtypes 7}*/
  typedef SNC_sphere_map<Kernel_, Items_,Mark_>     Self;
  typedef Items_                                    Items;
  typedef Kernel_                                   Kernel;
  typedef Mark_                                     Mark;
  typedef CGAL::SNC_structure<Kernel,Items,Mark>    SNC_structure;
  typedef typename Items::template Vertex<SNC_structure>  Base;
  //  typedef bool                                        Mark;
  typedef CGAL::Sphere_geometry<Kernel>             Sphere_kernel;

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

  typedef Infimaximal_box<typename Is_extended_kernel<Kernel>::value_type, Kernel> Infi_box;
  typedef typename Infi_box::Standard_kernel  Standard_kernel;
  typedef Infimaximal_box<typename Is_extended_kernel<Standard_kernel>::value_type, Standard_kernel> No_box;

  typedef Self                                              Vertex_base;
  typedef SNC_in_place_list_sm<Vertex_base>                 Vertex;
  typedef CGAL::In_place_list<Vertex,false>                 Vertex_list;
  typedef CGAL_ALLOCATOR(Vertex)                            Vertex_alloc;
  typedef typename Vertex_list::iterator                    Vertex_handle;
  typedef typename Vertex_list::const_iterator              Vertex_const_handle;
  typedef typename Vertex_list::iterator                    Vertex_iterator;
  typedef typename Vertex_list::const_iterator              Vertex_const_iterator;

  typedef typename Items::template SVertex<SNC_structure>   SVertex_base;
  typedef SNC_in_place_list_svertex<SVertex_base>           SVertex;
  typedef CGAL::In_place_list<SVertex,false>                SVertex_list;
  typedef CGAL_ALLOCATOR(SVertex)                           SVertex_alloc;
  typedef typename SVertex_list::iterator                   SVertex_handle;
  typedef typename SVertex_list::const_iterator             SVertex_const_handle;
  typedef typename SVertex_list::iterator                   SVertex_iterator;
  typedef typename SVertex_list::const_iterator             SVertex_const_iterator;

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

  typedef CGAL::Object_handle Object_handle;
  typedef std::list<Object_handle>    Object_list;
  typedef Object_list::iterator       Object_iterator;
  typedef Object_list::const_iterator Object_const_iterator;
  typedef Object_list::const_iterator Object_const_handle;

  typedef Vertex_handle       Constructor_parameter;
  typedef Vertex_const_handle Constructor_const_parameter;
  using Base::init_range; // AF add CR
 public:
  SNC_sphere_map(bool construct=false) : Base(), destruct(construct) {
    if(!construct) return;
    this->sncp() = new SNC_structure;
    init_range(this->sncp()->svertices_end());
    init_range(this->sncp()->shalfedges_end());
    init_range(this->sncp()->sfaces_end());
    this->shalfloop() = this->sncp()->shalfloops_end();
  }

  ~SNC_sphere_map() { if(destruct) delete this->sncp(); }

  SNC_sphere_map(const Base& v) : Base(v), destruct(false) {}

  SNC_sphere_map(const SNC_sphere_map&)=default;

  Self& operator=(const Self& M) {
    destruct = M.destruct;
    Base* b(this);
    *b = M;
    return *this;
  }

  void clear(bool clear_base=false) {
    if(clear_base) Base::clear();
  }

  template <typename H>
  bool is_sm_boundary_object(H h) const
  { return this->sncp()->is_sm_boundary_object(h); }

  template <typename H>
  Object_iterator& sm_boundary_item(H h)
  { return this->sncp()->sm_boundary_item(h); }

  template <typename H>
  void store_sm_boundary_item(H h, Object_iterator o)
  { this->sncp()->store_sm_boundary_item(h,o); }

  template <typename H>
  void undef_sm_boundary_item(H h)
  { this->sncp()->undef_sm_boundary_item(h); }

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

  void reset_sm_object_list(Object_list& L)
  { Object_iterator oit;
    CGAL_forall_iterators(oit,L) reset_sm_iterator_hash(oit);
    L.clear();
  }

  class SFace_cycle_iterator : public Object_iterator
  /*{\Mtypemember a generic iterator to an object in the boundary
  of a sface. Convertible to |Object_handle|.}*/
  { typedef Object_iterator Ibase;
  public:
    SFace_cycle_iterator() : Ibase() {}
    SFace_cycle_iterator(const Ibase& b) : Ibase(b) {}

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
    { CGAL_assertion_msg(0,"not impl."); return Object_handle(); }
  };

  class SFace_cycle_const_iterator : public Object_const_iterator
  /*{\Mtypemember a generic iterator to an object in the boundary
  of a facet. Convertible to |Object_handle|.}*/
  { typedef Object_const_iterator Ibase;
  public:
    SFace_cycle_const_iterator() : Ibase() {}
    SFace_cycle_const_iterator(const Ibase& b) : Ibase(b) {}

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
    { CGAL_assertion_msg(0,"not impl."); return Object_handle(); }
  };

  SFace_cycle_const_iterator sface_cycles_begin(SFace_const_handle f) const
  { return f->sface_cycles_begin(); }
  SFace_cycle_const_iterator sface_cycles_end(SFace_const_handle f) const
  { return f->sface_cycles_end(); }
  SFace_cycle_iterator sface_cycles_begin(SFace_handle f) const
  { return f->sface_cycles_begin(); }
  SFace_cycle_iterator sface_cycles_end(SFace_handle f) const
  { return f->sface_cycles_end(); }

 typedef CircFromIt<
        SHalfedge_const_iterator,
        move_shalfedge_around_svertex<SHalfedge_const_iterator> >
        SHalfedge_around_svertex_const_circulator;

  typedef CircFromIt<
        SHalfedge_const_iterator,
        move_shalfedge_around_sface<SHalfedge_const_iterator> >
        SHalfedge_around_sface_const_circulator;

  typedef CircFromIt<
        SHalfedge_iterator,
        move_shalfedge_around_svertex<SHalfedge_iterator> >
        SHalfedge_around_svertex_circulator;

  typedef CircFromIt<
        SHalfedge_iterator,
        move_shalfedge_around_sface<SHalfedge_iterator> >
        SHalfedge_around_sface_circulator;

  bool empty() const {
    return this->number_of_svertices() == 0 &&
      this->number_of_shalfedges() == 0 &&
      this->number_of_shalfloops() == 0 &&
      this->number_of_sfaces() == 0;
  }

  /*
  bool has_shalfloop() const {
    return shalfloop_ != 0;
  }

  SHalfloop_handle shalfloop() const {
    return shalfloop_;
  }
  */

  template <typename H>
  void make_twins(H h1, H h2) {
    h1->twin() = h2;
    h2->twin() = h1;
  }

  SVertex_handle new_svertex(const Sphere_point& p,
                           Mark m = Mark()) {
    SVertex_iterator sv;
    if (this->svertices_begin() == this->sncp()->svertices_end()) {
      sv = this->sncp()->new_halfedge_only();
      init_range(sv);
    }
    else {
      SVertex_iterator svn = this->svertices_end();
      sv =  this->sncp()->new_halfedge_only(svn);
      this->svertices_last() = sv;
    }
    sv->point() = p;
    sv->mark() = m;
    sv->center_vertex() = Vertex_handle((SNC_in_place_list_sm<Self>*) this);
    CGAL_NEF_TRACEN("new_svertex "<<&*sv);
    return sv;
  }
  /*
  SFace_handle new_sface() {
    SFace_iterator sf =  this->sncp()->new_sface_only();
    if ( this->sfaces_begin() == this->sncp()->sfaces_end()) init_range(sf);
    else this->sfaces_last() = sf;
    sf->center_vertex() = Vertex_handle((SNC_in_place_list_sm<Self>*) this);
    return sf;
  }
  */
  SFace_handle new_sface() {
    SFace_iterator sf;
    if ( this->sfaces_begin() == this->sncp()->sfaces_end()) {
      sf =  this->sncp()->new_sface_only();
      init_range(sf);
    } else {
      SFace_iterator sfn = this->sfaces_end();
      sf =  this->sncp()->new_sface_only(sfn);
      this->sfaces_last() = sf;
    }
    sf->center_vertex() = Vertex_handle((SNC_in_place_list_sm<Self>*) this);
    return sf;
  }

  /*
  SHalfedge_handle new_shalfedge_pair() {
    SHalfedge_iterator se = this->sncp()->new_shalfedge_only();
    SHalfedge_iterator set = this->sncp()->new_shalfedge_only();
    if(this->shalfedges_begin() == this->sncp()->shalfedges_end())
      init_range(se);
    this->shalfedges_last() = set;
    make_twins(se,set);
    return se;
  }
  */

  SHalfedge_handle new_shalfedge_pair() {
    SHalfedge_iterator se, set;
    if(this->shalfedges_begin() == this->sncp()->shalfedges_end()) {
      se = this->sncp()->new_shalfedge_only();
      set = this->sncp()->new_shalfedge_only();
      init_range(se);
    } else {
      SHalfedge_iterator sen = this->shalfedges_end();
      se = this->sncp()->new_shalfedge_only(sen);
      set = this->sncp()->new_shalfedge_only(sen);
    }
    this->shalfedges_last() = set;
    make_twins(se,set);
    return se;
  }

  SHalfloop_handle new_shalfloop_pair() {
    //    CGAL_assertion( !cv->has_shalfloop() );
    CGAL_assertion( !this->has_shalfloop() );
    SHalfloop_iterator sl =  this->sncp()->new_shalfloop_only();
    SHalfloop_iterator slt = this->sncp()->new_shalfloop_only();
    make_twins(sl,slt);
    //    cv->shalfloop() = sl;
    this->shalfloop() = sl;
    return sl;
  }

  void delete_svertex(SVertex_handle v) {
    if (this->svertices_begin() == this->svertices_last() )
      { CGAL_assertion(v == this->svertices_begin());
      init_range(this->sncp()->svertices_end()); }
    else if (this->svertices_begin() == v ) ++(this->svertices_begin());
    else if (this->svertices_last() == v ) --(this->svertices_last());
    this->sncp()->delete_halfedge_only(v);
  }

  void delete_shalfedge(SHalfedge_handle e) {
    if ( this->shalfedges_begin() == this->shalfedges_last() )
      { CGAL_assertion( e == this->shalfedges_begin() );
      init_range(this->sncp()->shalfedges_end()); }
    else if (this->shalfedges_begin() == e ) ++(this->shalfedges_begin());
    else if (this->shalfedges_last() == e ) --(this->shalfedges_last());
    this->sncp()->delete_shalfedge_only(e);
  }

  void delete_shalfedge_pair(SHalfedge_handle e) {
    delete_shalfedge(e->twin());
    delete_shalfedge(e);
  }

  void delete_sface(SFace_handle f) {
    if (this->sfaces_begin() == this->sfaces_last() )
      { CGAL_assertion( f == this->sfaces_begin() );
      init_range(this->sncp()->sfaces_end()); }
    else if (this->sfaces_begin() == f ) ++(this->sfaces_begin());
    else if (this->sfaces_last() == f )  --(this->sfaces_last());
    this->sncp()->delete_sface_only(f);
  }

  void delete_shalfloop_pair() {
    CGAL_assertion(this->has_shalfloop() );
    this->sncp()->delete_shalfloop_only(this->shalfloop()->twin());
    this->sncp()->delete_shalfloop_only(this->shalfloop());
    this->shalfloop() = this->sncp()->shalfloops_end();
  }

 protected:
  bool            destruct;
};  // SNC_sphere_map


} //namespace CGAL

#endif //  CGAL_SNC_SPHERE_MAP_H
