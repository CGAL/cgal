// Copyright (c) 1997-2000  Max-Planck-Institute Saarbruecken (Germany).
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
// Author(s)     : Michael Seel <seel@mpi-sb.mpg.de>
//                 Peter Hachenberger <hachenberger@mpi-sb.mpg.de>

#ifndef CGAL_SNC_SPHERE_MAP_H
#define CGAL_SNC_SPHERE_MAP_H

#include <CGAL/basic.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/Nef_2/Object_handle.h>
#include <CGAL/Nef_S2/SM_iteration.h>
#include <CGAL/Nef_S2/Generic_handle_map.h>
#include <CGAL/Nef_2/iterator_tools.h>
#include <CGAL/Nef_3/Infimaximal_box.h>
#include <CGAL/Nef_S2/Sphere_geometry.h>
#include <list>
#undef _DEBUG
#define _DEBUG 109
#include <CGAL/Nef_S2/debug.h>


CGAL_BEGIN_NAMESPACE

template <typename K, typename I> class SNC_structure;

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

template <typename Kernel_, typename Items_>
class SNC_sphere_map : public Items_::Vertex<SNC_structure<Kernel_, Items_> >
{

 public:
  /*{\Mtypes 7}*/
  typedef SNC_sphere_map<Kernel_, Items_>           Self;
  typedef Items_                                    Items;
  typedef Kernel_                                   Kernel;
  typedef SNC_structure<Kernel,Items>               SNC_structure;
  typedef typename Items::Vertex<SNC_structure>     Base;
  typedef bool                                      Mark;
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
  //  typedef typename Items::template Vertex<SNC_structure>    Vertex_base;
  //  typedef SNC_in_place_list_vertex<Vertex_base>             Vertex;
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

  //  Constructor_parameter cv;

  //  SNC_structure* sncp() { return sncp(); }
  //  SNC_structure* sncp() const { return sncp(); }

 public:
  SNC_sphere_map() {}
  SNC_sphere_map(const Base& v) : Base(v) {}
  
  //  SNC_sphere_map(const Self& M) : Base(M) {}

  /*
  Self& operator=(const Self& M) {
    cv = M.cv;
    return *this;
  }
  */

  void clear(bool clear_base=false) {
    if(clear_base) Base::clear(); 
  }

  template <typename H>
  bool is_boundary_object(H h) const
  { return sncp()->is_sm_boundary_object(h); }

  template <typename H>
  Object_iterator& boundary_item(H h)
  { return sncp()->sm_boundary_item(h); }

  template <typename H>
  void store_boundary_item(H h, Object_iterator o)
  { sncp()->store_sm_boundary_item(h,o); }

  template <typename H>
  void undef_boundary_item(H h)
  { sncp()->undef_sm_boundary_item(h); }

  void reset_iterator_hash(Object_iterator it)
  { SVertex_handle sv;
    SHalfedge_handle se;
    SHalfloop_handle sl;
    if ( assign(se,*it) ) { undef_boundary_item(se); return; }
    if ( assign(sl,*it) ) { undef_boundary_item(sl); return; }
    if ( assign(sv,*it) ) { undef_boundary_item(sv); return; }
  }

  void reset_sm_iterator_hash(Object_iterator it)
  { SVertex_handle sv;
    SHalfedge_handle se;
    SHalfloop_handle sl;
    if ( assign(se,*it) ) { 
      if( is_boundary_object(se)) 
	undef_sm_boundary_item(se); 
      return; 
    }
    if ( assign(sl,*it) ) { 
      if( is_boundary_object(sl)) 
	undef_sm_boundary_item(sl);
      return; 
    }
    if ( assign(sv,*it) ) { 
      if( is_boundary_object(sv)) 
	undef_sm_boundary_item(sv); 
      return; 
    }
  }

  void reset_object_list(Object_list& L)
  { Object_iterator oit;
    CGAL_forall_iterators(oit,L) reset_iterator_hash(oit);
    L.clear();
  }  

  void reset_sm_object_list(Object_list& L)
  { Object_iterator oit;
    CGAL_forall_iterators(oit,L) reset_sm_iterator_hash(oit);
    L.clear();
  }

  /*
  SVertex_const_iterator   svertices_begin()  const { return svertices_begin();}
  SVertex_const_iterator   svertices_end()    const { return svertices_end();}
  SHalfedge_const_iterator shalfedges_begin() const { return sedges_begin();}
  SHalfedge_const_iterator shalfedges_end()   const { return sedges_end();}
  SHalfloop_const_iterator shalfloops_begin() const { return shalfloop(); }
  SHalfloop_const_iterator shalfloops_end()   const 
    { return shalfloop_ != 0 ? shalfloop_+2 : shalfloop_; }
  SFace_const_iterator     sfaces_begin()     const { return sfaces_begin();}
  SFace_const_iterator     sfaces_end()       const { return sfaces_end();}
  */
  /*
  SHalfloop_const_iterator shalfloops_begin() const { return shalfloop(); }
  SHalfloop_const_iterator shalfloops_end()   const 
    { return shalfloop_ != 0 ? shalfloop_+2 : shalfloop_; }
  SHalfloop_iterator shalfloops_begin()  { return shalfloop(); }
  SHalfloop_iterator shalfloops_end()  
    { return shalfloop_ != 0 ? shalfloop_+2 : shalfloop_; } 
  */
  /*
  SVertex_iterator   svertices_begin()   { return svertices_begin();}
  SVertex_iterator   svertices_end()     { return svertices_end();}
  SHalfedge_iterator shalfedges_begin()  { return sedges_begin();}
  SHalfedge_iterator shalfedges_end()    { return sedges_end();}
  SHalfloop_iterator shalfloops_begin()  { return shalfloop(); }
  SHalfloop_iterator shalfloops_end()  
    { return shalfloop_ != 0 ? shalfloop_+2 : shalfloop_; }
  SFace_iterator     sfaces_begin()      { return sfaces_begin();}
  SFace_iterator     sfaces_end()        { return sfaces_end();}
  */
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

    operator Object_handle() const { return Ibase::operator*(); }
    const Object_handle& operator*() const { return Ibase::operator*(); }
    Object_handle  operator->() const 
    { CGAL_nef_assertion_msg(0,"not impl."); }
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
  /*
  Size_type number_of_svertices() const  { return number_of_svertices(); }
  Size_type number_of_shalfedges() const { return number_of_shalfedges();}
  Size_type number_of_sfaces() const     { return number_of_sfaces();}
  Size_type number_of_shalfloops() const { return number_of_shalfloops(); }
  */
  bool empty() const {
    return number_of_svertices() == 0 &&
      number_of_shalfedges() == 0 &&
      number_of_shalfloops() == 0 &&
      number_of_sfaces() == 0;
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
    SVertex_iterator svn = svertices_end();
    //    SVertex_iterator svn = svertices_end();
    SVertex_iterator sv =  sncp()->new_halfedge_only(svn);
    sv->vector() = p; 
    sv->mark() = m;
    if (svertices_begin() == sncp()->svertices_end()) init_range(sv);
    //    if (svertices_begin() == sncp()->svertices_end()) init_range(sv);
    else svertices_last() = sv;
    //    else svertices_last() = sv;
    //    sv->center_vertex() = cv;
    sv->center_vertex() = Vertex_handle((SNC_in_place_list_sm<Self>*) this);
    TRACEN("new_svertex "<<&*sv);
    return sv;
  }

  /*
  SVertex_handle new_svertex() { 
    svertices_.push_back( * get_vertex_node(SVertex())); 
    return --svertices_end(); 
  }
  */
  SFace_handle new_sface() {
    SFace_iterator sf =  sncp()->new_sface_only();
    //    if ( cv->sfaces_begin() == sncp()->sfaces_end() ) cv->init_range(sf);
    if ( sfaces_begin() == sncp()->sfaces_end() ) init_range(sf);
    //    else cv->sfaces_last() = sf;
    else sfaces_last() = sf;
    //    sf->center_vertex() = cv;
    sf->center_vertex() = Vertex_handle((SNC_in_place_list_sm<Self>*) this);
    return sf; 
  }

  /*
  SFace_handle new_sface() { 
    sfaces_.push_back( * get_face_node(SFace())); 
    return --sfaces_end(); 
  } 
  */

  SHalfedge_handle new_shalfedge_pair() {
    SHalfedge_iterator se = sncp()->new_shalfedge_only();
    SHalfedge_iterator set = sncp()->new_shalfedge_only();
    //  if ( cv->shalfedges_begin() == sncp()->shalfedges_end() ) cv->init_range(se);
    if ( shalfedges_begin() == sncp()->shalfedges_end() ) init_range(se);
    //    cv->shalfedges_last() = set;
    shalfedges_last() = set;
    make_twins(se,set);
    return se; 
  }

  /*
  SHalfedge_handle new_shalfedge_pair() { 
    SHalfedge* ep2 = get_halfedge_node(SHalfedge());
    SHalfedge* ep1 = get_halfedge_node(SHalfedge());
    sedges_.push_back( *ep1 );
    SHalfedge_handle e1 = --shalfedges_end();
    sedges_.push_back( *ep2 );
    SHalfedge_handle e2 = --shalfedges_end();
    make_twins(e1,e2); return e1; }
  */

  SHalfloop_handle new_shalfloop_pair() {
    //    CGAL_assertion( !cv->has_shalfloop() );
    CGAL_assertion( !has_shalfloop() );
    SHalfloop_iterator sl =  sncp()->new_shalfloop_only();
    SHalfloop_iterator slt = sncp()->new_shalfloop_only();
    make_twins(sl,slt);
    //    cv->shalfloop() = sl;
    shalfloop() = sl;
    return sl; 
  }


  /*
  SHalfloop_handle new_shalfloop_pair()
  { SHalfloop_handle ph = new SHalfloop[2];
    SHalfloop* pt(ph); ++pt;
    make_twins(ph,pt);
    shalfloop_=ph; return ph; }
  */

  /*
  SHalfedge_handle new_shalfedge_pair(const SHalfedge& e1)
  { const SHalfedge& e2 = *(e1.twin());
    SHalfedge* ep2 = new SHalfedge[2];
    SHalfedge* ep1 = ep2++;
    *ep1=e1; *ep2=e2;
    sedges_.push_back( *ep1 );
    SHalfedge_handle eh1 = --shalfedges_end();
    sedges_.push_back( *ep2 );
    SHalfedge_handle eh2 = --shalfedges_end();
    make_twins(eh1,eh2); return eh1; }
  */

  /*
  SHalfloop_handle new_shalfloop_pair(const SHalfloop& l1)
  { const SHalfloop& l2 = *(l1.twin());
    SHalfloop* ph = new SHalfloop[2];
    SHalfloop* pt(ph); ++pt;
    *ph=l1; *pt=l2; make_twins(ph,pt);
    shalfloop_=ph; return ph; }
  */
  
void delete_svertex(SVertex_handle v) {
  //  if ( cv->svertices_begin() == cv->svertices_last() ) 
  if (svertices_begin() == svertices_last() ) 
    //    { CGAL_assertion(v == cv->svertices_begin()); 
    { CGAL_assertion(v == svertices_begin()); 
    //    cv->init_range(sncp()->svertices_end()); }
  init_range(sncp()->svertices_end()); }
  //  else if ( cv->svertices_begin() == v ) ++(cv->svertices_begin());
  else if (svertices_begin() == v ) ++(svertices_begin());
  //  else if ( cv->svertices_last() == v ) --(cv->svertices_last());
  else if (svertices_last() == v ) --(svertices_last());
  sncp()->delete_halfedge_only(v);
}

void delete_shalfedge(SHalfedge_handle e) {
  //  if ( cv->shalfedges_begin() == cv->shalfedges_last() ) 
  if ( shalfedges_begin() == shalfedges_last() ) 
    //    { CGAL_assertion( e == cv->shalfedges_begin() ); 
    { CGAL_assertion( e == shalfedges_begin() ); 
    //    cv->init_range(sncp()->shalfedges_end()); }
  init_range(sncp()->shalfedges_end()); }
  //  else if ( cv->shalfedges_begin() == e ) ++(cv->shalfedges_begin());
  else if (shalfedges_begin() == e ) ++(shalfedges_begin());
  //  else if ( cv->shalfedges_last() == e ) --(cv->shalfedges_last());
  else if (shalfedges_last() == e ) --(shalfedges_last());
  sncp()->delete_shalfedge_only(e);
}

void delete_shalfedge_pair(SHalfedge_handle e) { 
  delete_shalfedge(e->twin()); 
  delete_shalfedge(e); 
}

void delete_sface(SFace_handle f) {
  //  if ( cv->sfaces_begin() == cv->sfaces_last() ) 
  if (sfaces_begin() == sfaces_last() ) 
    //    { CGAL_assertion( f == cv->sfaces_begin() ); 
    { CGAL_assertion( f == sfaces_begin() ); 
    //    cv->init_range(sncp()->sfaces_end()); }
  init_range(sncp()->sfaces_end()); }
  //  else if ( cv->sfaces_begin() == f ) ++(cv->sfaces_begin());
  else if (sfaces_begin() == f ) ++(sfaces_begin());
  //  else if ( cv->sfaces_last() == f )  --(cv->sfaces_last());
  else if (sfaces_last() == f )  --(sfaces_last());
  sncp()->delete_sface_only(f); 
}

void delete_shalfloop_pair() {
  //  CGAL_assertion( cv->has_shalfloop() );
  CGAL_assertion(has_shalfloop() );
  //  sncp()->delete_shalfloop_only(cv->shalfloop()->twin());  
  sncp()->delete_shalfloop_only(shalfloop()->twin());  
  //  sncp()->delete_shalfloop_only(cv->shalfloop());  
  sncp()->delete_shalfloop_only(shalfloop());  
  //  cv->shalfloop() = sncp()->shalfloops_end();
  shalfloop() = sncp()->shalfloops_end();
}
  
};  // SNC_sphere_map


CGAL_END_NAMESPACE

#endif //  CGAL_SNC_SPHERE_MAP_H
