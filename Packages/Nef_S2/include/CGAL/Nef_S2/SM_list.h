#ifndef CGAL_SM_LIST_H
#define CGAL_SM_LIST_H

#include <CGAL/In_place_list.h>
#include <CGAL/Nef_S2/SM_items.h>
#include <CGAL/Nef_S2/Sphere_geometry.h>
#include <CGAL/Nef_2/Object_handle.h>
#include <CGAL/Nef_S2/Generic_handle_map.h>
#include <CGAL/Nef_2/iterator_tools.h>
#include <list>

CGAL_BEGIN_NAMESPACE

/*
template <typename HE>
class move_edge_around_svertex {
public:
  void forward(HE& e) const  { e = (e->sprev()->twin()); }
  void backward(HE& e) const { e = (e->twin()->snext()); }
};

template <typename HE>
struct move_edge_around_sface {
  void forward(HE& e)  const { e = (e->snext()); }
  void backward(HE& e) const { e = (e->sprev()); }
};
*/

template < class SVertex>
class SNC_in_place_list_svertex
    : public SVertex, 
      public In_place_list_base<SNC_in_place_list_svertex<SVertex> > {
public:
    typedef SNC_in_place_list_svertex<SVertex> Self;
    //    typedef typename SVertex::SVertex_handle       SVertex_handle;
    //    typedef typename SVertex::SVertex_const_handle SVertex_const_handle;
    SNC_in_place_list_svertex() {}
    SNC_in_place_list_svertex(const SVertex& v)   // down cast
        : SVertex(v) {}
    Self& operator=( const Self& v) {
        // This self written assignment avoids that assigning vertices will
        // overwrite the list linking of the target vertex.
        *((SVertex*)this) = ((const SVertex&)v);
        return *this;
    }
};

template < class SHalfedge>
class SNC_in_place_list_shalfedge
    : public SHalfedge, 
      public In_place_list_base<SNC_in_place_list_shalfedge<SHalfedge> > {
public:
    typedef SNC_in_place_list_shalfedge<SHalfedge> Self;
    //    typedef typename SHalfedge::SHalfedge_handle       SHalfedge_handle;
    //    typedef typename SHalfedge::SHalfedge_const_handle SHalfedge_const_handle;
    SNC_in_place_list_shalfedge() {}
    SNC_in_place_list_shalfedge(const SHalfedge& v)   // down cast
        : SHalfedge(v) {}
    Self& operator=( const Self& v) {
        // This self written assignment avoids that assigning vertices will
        // overwrite the list linking of the target vertex.
        *((SHalfedge*)this) = ((const SHalfedge&)v);
        return *this;
    }
};

/*
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
*/

template < class SFace>
class SNC_in_place_list_sface
    : public SFace, 
      public In_place_list_base<SNC_in_place_list_sface<SFace> > {
public:
    typedef SNC_in_place_list_sface<SFace> Self;
    //    typedef typename SFace::SFace_handle       SFace_handle;
    //    typedef typename SFace::SFace_const_handle SFace_const_handle;
    SNC_in_place_list_sface() {}
    SNC_in_place_list_sface(const SFace& v)   // down cast
        : SFace(v) {}
    Self& operator=( const Self& v) {
        // This self written assignment avoids that assigning vertices will
        // overwrite the list linking of the target vertex.
        *((SFace*)this) = ((const SFace&)v);
        return *this;
    }
};

/*
template<typename Kernel_,typename Items_>
class SM_list {

 public:
  typedef Kernel_                Kernel;
  typedef Items_                 Items;
  typedef SM_list<Kernel,Items>  Self;

  typedef size_t                 Size_type;
  typedef bool                   Mark;

  typedef CGAL::Sphere_geometry<Kernel>                 Sphere_kernel;
  typedef typename Sphere_kernel::Sphere_point          Sphere_point;
  typedef typename Sphere_kernel::Sphere_segment        Sphere_segment;
  typedef typename Sphere_kernel::Sphere_circle         Sphere_circle;
  typedef typename Sphere_kernel::Sphere_direction      Sphere_direction;

  typedef typename Items::template SVertex<Self>        SVertex_base;
  typedef SNC_in_place_list_svertex<SVertex_base>       SVertex;
  typedef CGAL::In_place_list<SVertex,false>            SVertex_list;
  typedef CGAL_ALLOCATOR(SVertex)                       SVertex_alloc;
  typedef typename SVertex_list::iterator               SVertex_handle;
  typedef typename SVertex_list::const_iterator         SVertex_const_handle;
  typedef typename SVertex_list::iterator               SVertex_iterator;
  typedef typename SVertex_list::const_iterator         SVertex_const_iterator;

  typedef typename Items::template SHalfedge<Self>      SHalfedge_base;
  typedef SNC_in_place_list_shalfedge<SHalfedge_base>   SHalfedge;
  typedef CGAL::In_place_list<SHalfedge,false>          SHalfedge_list;
  typedef CGAL_ALLOCATOR(SHalfedge)                     SHalfedge_alloc;
  typedef typename SHalfedge_list::iterator             SHalfedge_handle;
  typedef typename SHalfedge_list::const_iterator       SHalfedge_const_handle;
  typedef typename SHalfedge_list::iterator             SHalfedge_iterator;
  typedef typename SHalfedge_list::const_iterator       SHalfedge_const_iterator;

  typedef typename Items::template SHalfloop<Self>      SHalfloop_base;
  typedef SNC_in_place_list_shalfloop<SHalfloop_base>   SHalfloop;
  typedef CGAL::In_place_list<SHalfloop,false>          SHalfloop_list;
  typedef CGAL_ALLOCATOR(SHalfloop)                     SHalfloop_alloc;
  typedef typename SHalfloop_list::iterator             SHalfloop_handle;
  typedef typename SHalfloop_list::const_iterator       SHalfloop_const_handle;
  typedef typename SHalfloop_list::iterator             SHalfloop_iterator;
  typedef typename SHalfloop_list::const_iterator       SHalfloop_const_iterator;

  typedef typename Items::template SFace<Self>          SFace_base;
  typedef SNC_in_place_list_sface<SFace_base>           SFace;
  typedef CGAL::In_place_list<SFace,false>              SFace_list;
  typedef CGAL_ALLOCATOR(SFace)                         SFace_alloc;
  typedef typename SFace_list::iterator                 SFace_handle;
  typedef typename SFace_list::const_iterator           SFace_const_handle;
  typedef typename SFace_list::iterator                 SFace_iterator;
  typedef typename SFace_list::const_iterator           SFace_const_iterator;

  typedef CGAL::Object_handle                           Object_handle;
  typedef std::list<Object_handle>                      Object_list;
  typedef typename Object_list::iterator                Object_iterator;
  typedef typename Object_list::const_iterator          Object_const_iterator;

  typedef Generic_handle_map<Object_iterator>           Handle_to_iterator_map;

  typedef CircFromIt<
    SHalfedge_const_iterator, 
    move_edge_around_svertex<SHalfedge_const_iterator> > 
    SHalfedge_around_svertex_const_circulator;
  
  typedef CircFromIt<
    SHalfedge_const_iterator, 
    move_edge_around_sface<SHalfedge_const_iterator> > 
    SHalfedge_around_sface_const_circulator;

  typedef CircFromIt<
    SHalfedge_iterator, 
    move_edge_around_svertex<SHalfedge_iterator> > 
    SHalfedge_around_svertex_circulator;
  
  typedef CircFromIt<
    SHalfedge_iterator, 
    move_edge_around_sface<SHalfedge_iterator> > 
    SHalfedge_around_sface_circulator;

  class SFace_cycle_iterator : public Object_iterator 
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
    { CGAL_assertion_msg(0,"not impl."); }
  };

  class SFace_cycle_const_iterator : public Object_const_iterator 
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
    { CGAL_assertion_msg(0,"not impl."); }
  };

  SM_list() : 
    sm_boundary_item_(undef_), 
    svertices_(), 
    shalfedges_(), 
    sfaces_(), 
    shalfloops_() {}

  ~SM_list() { clear(); } 

  SM_list(const Self& L) :
    sm_boundary_item_(undef_),
    svertices_(L.svertices_),
    shalfedges_(L.shalfedges_),
    shalfloops_(L.shalfloops),
    sfaces_(L.sfaces_) {
    pointer_update(D); 
  }

  Self& operator=(const Self& L) {
    if ( this == &L )
      return *this; 
    clear();
    svertices_ = L.svertices_; 
    sfaces_ = L.sfaces_;
    shalfedges_ = L.shalfedges_;
    shalfloops_ = L.shalfloops_;
    pointer_update(L);
    return *this;
  }

  SVertex_const_iterator svertices_begin() const    {return svertices_.begin();}
  SVertex_const_iterator svertices_end() const      {return svertices_.end();}
  SHalfedge_const_iterator shalfedges_begin() const {return shalfedges_.begin();}
  SHalfedge_const_iterator shalfedges_end() const   {return shalfedges_.end();}
  SFace_const_iterator sfaces_begin() const         {return sfaces_.begin();}
  SFace_const_iterator sfaces_end() const           {return sfaces_.end();}
  SHalfloop_const_iterator shalfloops_begin() const {return shalfloops_.begin(); }
  SHalfloop_const_iterator shalfloops_end()   const {return shalfloops_.end(); }

  SVertex_iterator   svertices_begin()  { return svertices_.begin();}
  SVertex_iterator   svertices_end()    { return svertices_.end();}
  SHalfedge_iterator shalfedges_begin() { return shalfedges_.begin();}
  SHalfedge_iterator shalfedges_end()   { return shalfedges_.end();}
  SFace_iterator     sfaces_begin()     { return sfaces_.begin();}
  SFace_iterator     sfaces_end()       { return sfaces_.end();}
  SHalfloop_iterator shalfloops_begin() {return shalfloops_.begin(); }
  SHalfloop_iterator shalfloops_end()   {return shalfloops_.end(); }

  bool has_shalfloop() const { return shalfloops_.begin() != shalfloops_.end(); }
  SHalfloop_handle shalfloop()             { return shalfloops_.begin(); }
  SHalfloop_const_handle shalfloop() const { return shalfloops_.begin(); }

  Size_type number_of_svertices() const  { return svertices_.size();}
  Size_type number_of_shalfedges() const { return shalfedges_.size();}
  Size_type number_of_sedges() const     { return shalfedges_.size()/2;}
  Size_type number_of_sfaces() const     { return sfaces_.size();}
  Size_type number_of_shalfloops() const { return shalfloops_.size();}

 protected:
  SVertex_alloc svertex_allocator;
  SVertex* get_svertex_node( const SVertex& t) {
    SVertex* p = svertex_allocator.allocate(1);
    svertex_allocator.construct( p, SVertex());
    return p;
  }
  void put_svertex_node( SVertex* p) {
    svertex_allocator.destroy(p);
    svertex_allocator.deallocate( p, 1);
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

  template <typename H>
  void make_twins(H h1, H h2) { 
    h1->twin() = h2; 
    h2->twin() = h1; 
  }

 public:
  SVertex_handle new_svertex() { 
    svertices_.push_back( * get_svertex_node(SVertex())); 
    return --svertices_end(); 
  }

  SFace_handle new_sface() { 
    sfaces_.push_back( * get_sface_node(SFace())); 
    return --sfaces_end(); 
  } 

  SHalfedge_handle new_shalfedge_pair() { 
    SHalfedge* ep2 = get_shalfedge_node(SHalfedge());
    SHalfedge* ep1 = get_shalfedge_node(SHalfedge());
    shalfedges_.push_back( *ep1 );
    SHalfedge_handle e1 = --shalfedges_end();
    shalfedges_.push_back( *ep2 );
    SHalfedge_handle e2 = --shalfedges_end();
    make_twins(e1,e2);
    return e1; 
  }

  SHalfloop_handle new_shalfloop_pair() {
    CGAL_assertion(shalfloops_.empty());
    SHalfloop* lp1 = get_shalfloop_node(SHalfloop());
    SHalfloop* lp2 = get_shalfloop_node(SHalfloop());
    shalfloops_.push_back(*lp1);
    SHalfloop_handle l1 = --shalfloops_end();
    shalfloops_.push_back(*lp2);
    SHalfloop_handle l2 = --shalfloops_end();
    make_twins(l1,l2);
    return l1; 
  }

  void delete_svertex(SVertex_handle h) { 
    svertices_.erase(h); 
    put_svertex_node(&*h); 
  }

  void delete_sface(SFace_handle h) { 
    sfaces_.erase(h); 
    put_sface_node(&*h);
  }

  void delete_shalfedge_pair(SHalfedge_handle h) { 
    SHalfedge_handle t = h->twin();
    shalfedges_.erase(h);
    shalfedges_.erase(t);
    put_shalfedge_node(&*h);
    put_shalfedge_node(&*t);
  }

  void delete_shalfloop_pair() { 
    CGAL_assertion(shalfloops_.size()==2);
    SHalfloop_handle h = shalfloops_begin();
    SHalfloop_handle t = h->twin();
    shalfloops_.erase(h);
    shalfloops_.erase(t);
    put_shalfloop_node(&*h);
    put_shalfloop_node(&*t);
  }

  SHalfedge_handle new_shalfedge_pair(const SHalfedge& e1) {
    const SHalfedge& e2 = *(e1.twin());
    SHalfedge* ep2 = new SHalfedge[2];
    SHalfedge* ep1 = ep2++;
    *ep1=e1; *ep2=e2;
    shalfedges_.push_back( *ep1 );
    SHalfedge_handle eh1 = --shalfedges_end();
    shalfedges_.push_back( *ep2 );
    SHalfedge_handle eh2 = --shalfedges_end();
    make_twins(eh1,eh2);
    return eh1;
  }

  SHalfloop_handle new_shalfloop_pair(const SHalfloop& l1) {
    const SHalfedge& l2 = *(l1.twin());
    SHalfloop* lp2 = new SHalfloop[2];
    SHalfloop* lp1 = lp2++;
    *lp1=l1; *lp2=l2;
    shalfloops_.push_back( *lp1 );
    SHalfloop_handle eh1 = --shalfloops_end();
    shalfloops_.push_back( *lp2 );
    SHalfloop_handle eh2 = --shalfloops_end();
    make_twins(lh1,lh2);
    return lh1;
  }

  void clear() { 
    sm_boundary_item_.clear(undef_);
    svertices_.destroy(); 
    sfaces_.destroy();
    while ( shalfedges_begin() != shalfedges_end() )
      delete_shalfedge_pair( shalfedges_begin() );
    if (!shalfloops_.empty())
      delete_shalfloop_pair();
  }

  template <typename H>
  bool is_sm_boundary_object(H h) const
  { return sm_boundary_item_[h]!=undef_; }

  template <typename H>
  Object_iterator& sm_boundary_item(H h)
  { return sm_boundary_item_[h]; }

  template <typename H>
  void store_sm_boundary_item(H h, Object_iterator o)
  { sm_boundary_item_[h] = o; }

  template <typename H>
  void undef_sm_boundary_item(H h)
  { CGAL_assertion(sm_boundary_item_[h]!=undef_);
    sm_boundary_item_[h] = undef_; }

  void reset_sm_iterator_hash(Object_iterator it) { 
    SVertex_handle sv;
    SHalfedge_handle se;
    SHalfloop_handle sl;
    if ( assign(se,*it) ) { undef_sm_boundary_item(se); return; }
    if ( assign(sl,*it) ) { undef_sm_boundary_item(sl); return; }
    if ( assign(sv,*it) ) { undef_sm_boundary_item(sv); return; }
  }

  void reset_sm_object_list(Object_list& L) { 
    Object_iterator oit;
    CGAL_forall_iterators(oit,L)
      reset_sm_iterator_hash(oit);
    L.clear();
  }

 protected:
  static Object_iterator  undef_;

  Handle_to_iterator_map  sm_boundary_item_;
  SVertex_list            svertices_;
  SHalfedge_list          shalfedges_;
  SFace_list              sfaces_;
  SHalfloop_list          shalfloops_; 

  void pointer_update(const Self& L);
};

template <typename Kernel,typename Items> 
typename SM_list<Kernel,Items>::Object_iterator SM_list<Kernel,Items>::undef_;

template <typename K,typename I>
void SM_list<K,I>::
pointer_update(const SM_list<K,I>& L) {

  CGAL::Unique_hash_map<SVertex_const_handle,SVertex_handle>     VM;
  CGAL::Unique_hash_map<SHalfedge_const_handle,SHalfedge_handle> EM;
  CGAL::Unique_hash_map<SHalfloop_const_handle,SHalfloop_handle> LM;
  CGAL::Unique_hash_map<SFace_const_handle,SFace_handle>         FM;

  SVertex_const_iterator vc = L.svertices_begin();
  SVertex_iterator v = svertices_begin();
  for ( ; vc != L.svertices_end(); ++vc,++v) VM[vc] = v;
  VM[L.svertices_end()] = svertices_end();

  SHalfedge_const_iterator ec = L.shalfedges_begin();
  SHalfedge_iterator e = shalfedges_begin();
  for ( ; ec != L.shalfedges_end(); ++ec,++e) { 
    EM[ec] = e;
    e->mark() = ec->mark();
  }
  EM[L.shalfedges_end()] = shalfedges_end();

  SFace_const_iterator fc = L.sfaces_begin();
  SFace_iterator f = sfaces_begin();
  for ( ; fc != L.sfaces_end(); ++fc,++f) FM[fc] = f;
  FM[L.sfaces_end()] = sfaces_end();

  SHalfloop_iterator l;
  if (!L.shalfloops_.empty()) {
    LM[L.shalfloop()] = shalfloop();
    LM[L.shalfloop()->twin()] = shalfloop()->twin();
    l = shalfloop();
    shalfloop()->mark() = L.shalfloop()->mark();
    shalfloop()->twin()->mark() = L.shalfloop()->twin()->mark();    
    if( !l->is_twin() && L.shalfloop()->is_twin()) l->mark() = l->twin()->mark();
  }

  for (v = svertices_begin(); v != svertices_end(); ++v) {
    // Local Graph update: (SVertices are postponed/updated as Edges)
    v->out_sedge() = EM[v->out_sedge()];
    v->incident_sface() = FM[v->incident_sface()];
  }
  // Edge update:
  for (e = shalfedges_begin(); e != shalfedges_end(); ++e) {
    e->twin() = EM[e->twin()];
    e->sprev() = EM[e->sprev()];
    e->snext() = EM[e->snext()];
    e->source() = VM[e->source()];
    e->incident_sface() = FM[e->incident_sface()];
  }
  for ( ec = L.shalfedges_begin(), e = shalfedges_begin();
	ec != L.shalfedges_end(); ++ec, ++e) {
    if( !e->is_twin() && ec->is_twin()) e->mark() = e->twin()->mark();
  }

  for (l = shalfloops_begin(); l != shalfloops_end(); ++l) {
    //    l->twin() = LM[l->twin()];
    l->incident_sface() = FM[l->incident_sface()];
  }

  for (f = sfaces_begin(); f != sfaces_end(); ++f) {
    SFace_cycle_iterator fci; 
    for(fci = f->boundary_entry_objects().begin();
	fci != f->boundary_entry_objects().end(); ++fci) {
      if ( fci.is_svertex() ) 
      { v = SVertex_handle(fci);
	*fci = Object_handle(VM[v]); store_sm_boundary_item(v,fci); }
      else if ( fci.is_shalfedge() ) 
      { e = SHalfedge_handle(fci);
	*fci = Object_handle(EM[e]); store_sm_boundary_item(e,fci); }
      else if ( fci.is_shalfloop() ) 
      { l = SHalfloop_handle(fci);
	*fci = Object_handle(LM[l]); store_sm_boundary_item(l,fci); }
      else CGAL_assertion_msg(0,"damn wrong boundary item in face.");
    }
  }
}
*/

CGAL_END_NAMESPACE

#endif // CGAL_SM_LIST_H
