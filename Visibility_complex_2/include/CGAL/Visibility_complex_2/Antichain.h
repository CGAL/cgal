// Copyright (c) 2001-2004  ENS of Paris (France).
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
// $URL$
// $Id$
//
// Author(s)     : Pierre Angelier, Michel Pocchiola


#ifndef CGAL_VISIBILITY_COMPLEX_2_ANTICHAIN_H
#define CGAL_VISIBILITY_COMPLEX_2_ANTICHAIN_H

#include <CGAL/basic.h>
#include <CGAL/Visibility_complex_2/function_objects.h>
#include <CGAL/In_place_list.h>

#define CGAL_MULTISET_NO_ASSERTIONS
#define CGAL_MULTISET_NO_PRECONDITIONS
#define CGAL_MULTISET_NO_POSTCONDITIONS
#define CGAL_MULTISET_NO_WARNINGS
#include <CGAL/Multiset.h>

#include <CGAL/Visibility_complex_2/Items.h>
#include <CGAL/Visibility_complex_2/Flip_traits.h>
#include <CGAL/Visibility_complex_2/antichain_iterators.h>
#include <CGAL/Visibility_complex_2/ccw_cw_traits.h>
#include <CGAL/Visibility_complex_2/function_objects.h>
#include <CGAL/Visibility_complex_2/Sweep_iterator.h>
#include <CGAL/Visibility_complex_2/Bitangent_2.h>


#include <queue>
#include <list>
#include <set>
#include <vector>
#include <map>
#include <string>
#include <sstream>

CGAL_BEGIN_NAMESPACE

namespace Visibility_complex_2_details {

template<class T> class VC_in_place_list_base: public In_place_list_base<T> {
public:
  VC_in_place_list_base() {
    this->next_link=0;
    this->prev_link=0;
  }    
};

template < class Vertex_base>
class Vertex
  : public Vertex_base, 
    public VC_in_place_list_base< Vertex<Vertex_base> > 
{
public:
  typedef Vertex< Vertex_base> Self;
  typedef typename Vertex_base::Vertex_handle     Vertex_handle;
  typedef typename Vertex_base::Edge_handle       Edge_handle;
  typedef typename Vertex_base::Disk_handle    Disk_handle;
  typedef typename Vertex_base::Bitangent_2       Bitangent_2;
  typedef typename Vertex_base::Type              Type;

  Vertex() {}
  Vertex(Type t , Disk_handle start , 
         Disk_handle finish) 
    : Vertex_base(t,start,finish),final_antichain_(false) {}
  Vertex(Edge_handle start , Edge_handle finish)
    : Vertex_base(start,finish),final_antichain_(false)   {}
  Vertex(const Bitangent_2& b) 
    : Vertex_base(b),final_antichain_(false) {}
  Vertex( const Vertex_base& v)   // down cast
    : Vertex_base(v),final_antichain_(false) {}
  Vertex(const Vertex&sibling,bool reverse,Type t)
    : Vertex_base(sibling,reverse,t),final_antichain_(false) {}
  Self& operator=( const Self& v) {
    this->Vertex_base::operator=(v);
    final_antichain_=v.final_antichain_;
    return *this;
  }
  // Wether this vertex is the pi of a vertex of the original antichain,
  // in which case the sweep should not go beyond.
  bool final_antichain() const {
    return final_antichain_;
  }
  void set_final_antichain(bool b) {
    final_antichain_=b;
  }
private:
  bool final_antichain_;
};

// -----------------------------------------------------------------------------

template <class Edge_base>
class Edge
  : public Edge_base, 
    public VC_in_place_list_base< Edge<Edge_base> > 
{
public:
  typedef Edge<Edge_base> Self;
  typedef typename Edge_base::Vertex_handle  Vertex_handle;
  typedef typename Edge_base::Disk_handle Disk_handle;

  Edge() {}                   
  Edge(bool s, Disk_handle p)
    : Edge_base(s,p) {}
  Edge(Vertex_handle v0 , Vertex_handle v1 , 
       Disk_handle p)
    : Edge_base(v0,v1,p) {}
  Edge( const Edge_base& h)
    : Edge_base(h) {}
  Self& operator=( const Self& e) {
    this->Edge_base::operator=(e);
    return *this;
  }
  bool is_in_antichain() const {

    return this->next_link;
  }
};


template < class Face_base>
class Face
  : public Face_base {
private:
  bool in_initial_antichain;
public:
  typedef Face< Face_base>  Self;
  Face() : in_initial_antichain(false) {}
  bool is_in_initial_antichain() const {
    return in_initial_antichain;
  }
  void set_in_initial_antichain(bool b) {
    in_initial_antichain=b;
  }
};


template<class T, bool b> class VC_in_place_list 
  :public In_place_list<T,b> {
  typedef In_place_list<T,b> Base;
public: 
  void detach() {
    this->node->next_link=&*(this->node);
    this->node->prev_link=&*(this->node);    
  }
  void erase(T* pos) {
    CGAL_precondition(pos->next_link&&pos->prev_link);
    Base::erase(pos);
    pos->next_link=0;
    pos->prev_link=0;
  }
  typename Base::iterator insert(typename Base::iterator pos, T& x) {
    CGAL_precondition(x.next_link==0&&x.prev_link==0);
    return Base::insert(pos, x);
  }
};

template < class Gtr_ , 
	   class It   = Items ,
	   class Flip = Flip_traits >
class Antichain 
  : public VC_in_place_list<
  Edge<
    typename It::template Edge_wrapper<
      Antichain<Gtr_, It, Flip>
      >::Edge>, false>
{
public:
  // -------------------------------------------------------------------------
  typedef Gtr_                                          Gt;
  typedef typename Gtr_::Disk                           Disk; 
  typedef typename Gtr_::Point_2                        Point_2;
  typedef typename Gtr_::Bitangent_2                    Bitangent_2;
  typedef Bitangent_2                                   BT;
  typedef typename BT::Disk_handle                      Disk_handle;
  // -------------------------------------------------------------------------
  typedef Antichain<Gtr_,It,Flip>    Self;
  //  typedef Antichain<Gtr_,It,Flip>    Antichain;
  // -------------------------------------------------------------------------
  typedef Left_ccw_traits_tp<Self>      Left_ccw_traits;
  typedef Right_ccw_traits_tp<Self>     Right_ccw_traits;
  typedef Left_cw_traits_tp<Self>       Left_cw_traits;
  typedef Right_cw_traits_tp<Self>      Right_cw_traits;
  typedef Right_ccw_traits                              Ccw_traits;
  typedef Right_cw_traits                               Cw_traits;
  // -------------------------------------------------------------------------
  /*
    typedef typename Visibility_complex_backward_flip_traits:: Flip_wrapper<Self>
    Flip_wrapper;
  */
  typedef typename  Flip::template Flip_wrapper<Self>             Flip_wrapper;
  typedef  Flip_wrapper FW;
  typedef typename FW::Flip_traits            Flip_traits;
  // -------------------------------------------------------------------------
  typedef typename It::template Vertex_wrapper<Self>             Vertex_wrapper;
  typedef typename Vertex_wrapper::Vertex               Vertex_base;
  typedef Vertex<Vertex_base>       Vertex;
  typedef Vertex*                                       Vertex_handle;
  typedef const Vertex*                                 Vertex_const_handle;
  typedef typename VC_in_place_list<Vertex,false>::iterator Minimals_iterator;
  // -------------------------------------------------------------------------
  typedef typename It::template Edge_wrapper<Self>               Edge_wrapper;
  typedef typename Edge_wrapper::Edge                   Edge_base;
  typedef Edge<Edge_base>           Edge;
  typedef Edge*                           		Edge_handle;
  typedef const Edge*                           	Edge_const_handle;
  typedef VC_in_place_list<Edge,false> Base;
  typedef typename VC_in_place_list<Edge,false>::iterator           Edge_iterator;
  typedef typename VC_in_place_list<Edge,false>::reverse_iterator   Edge_reverse_iterator;
  typedef typename VC_in_place_list<Edge,false>::const_iterator     Edge_const_iterator;
  // -------------------------------------------------------------------------
  typedef typename It::template Face_wrapper<Self>               Face_wrapper;
  typedef typename Face_wrapper::Face                   Face_base;
  typedef Face< Face_base>           Face;
  typedef Face*                           		  Face_handle;
  typedef const Face*                           	  Face_const_handle;
  typedef Antichain_face_iterator<Self,Face,
                                     Face&,Face_handle> Face_iterator;
  typedef Antichain_face_iterator<Self,Face,
                                     const Face&,const Face_handle> 
  Face_const_iterator;
  // -------------------------------------------------------------------------
  typedef Antichain_vertex_iterator<Self,Vertex, Vertex&,
                                       Vertex_handle, typename Ccw_traits::Sup> 
  Vertex_iterator;
  typedef Antichain_vertex_iterator<Self,Vertex, const Vertex&,
                                       const Vertex_handle, 
                                       typename Ccw_traits::Sup> 
  Vertex_const_iterator;
  typedef Antichain_vertex_iterator<Self,Vertex, Vertex&,
                                       Vertex_handle, typename Cw_traits::Sup> 
  Vertex_cw_iterator;
  typedef Antichain_vertex_iterator<Self,Vertex, const Vertex&,
                                       const Vertex_handle,
                                       typename Cw_traits::Sup> 
  Vertex_cw_const_iterator;
  // -------------------------------------------------------------------------
  typedef Visibility_complex_2_details::Linear_sweep_iterator<Antichain,Vertex,
                                                   Vertex&,Vertex_handle,
                                                   typename Gt::Is_upward_directed>
  Linear_sweep_iterator;
  typedef Visibility_complex_2_details::Linear_sweep_iterator<Self,Vertex,
                                                   const Vertex&,const Vertex_handle,
                                                   typename Gt::Is_upward_directed>
  Linear_sweep_const_iterator;
  // -------------------------------------------------------------------------
  typedef Visibility_complex_2_details::Sweep_iterator<Self,Vertex,
                                            Vertex&,Vertex_handle>
  Sweep_iterator;
  typedef Visibility_complex_2_details::Sweep_iterator<Self,Vertex,
                                            const Vertex&,const Vertex_handle>
  Sweep_const_iterator;
  // -------------------------------------------------------------------------
  typedef typename BT::Type_util               Type_util;
  // -------------------------------------------------------------------------
  using Base::destroy;
  using Base::begin;
  using Base::end;
  using Base::rbegin;
  using Base::rend;

  typedef Visibility_complex_2_details::Constraint_input Constraint_input;
  typedef typename std::vector<Disk_handle>::const_iterator
    Disk_const_iterator;
  typedef typename std::vector<Vertex>::const_iterator
    Constraint_const_iterator;
  typedef typename std::vector<Disk_handle>::iterator Disk_iterator;
  typedef typename std::vector<Vertex>::iterator Constraint_iterator;

  Disk_iterator disks_begin() {
    return disks.begin();
  }
  Disk_iterator disks_end() {
    return disks.end();
  }
  Constraint_iterator constraints_begin() {
    return constraints.begin();
  }
  Constraint_iterator constraints_end() {
    return constraints.end();
  }
  Disk_const_iterator disks_begin() const {
    return disks.begin();
  }
  Disk_const_iterator disks_end() const {
    return disks.end();
  }
  Constraint_const_iterator constraints_begin() const {
    return constraints.begin();
  }
  Constraint_const_iterator constraints_end() const {
    return constraints.end();
  }
private:
public:
  bool                        valid_;
  bool                        straight_sweep_;
  bool                        linear_space_;
  Face_handle                 infinite_face_;
  VC_in_place_list<Vertex,false> minimals_ccw_;
  VC_in_place_list<Vertex,false> minimals_cw_;
  std::vector<Disk_handle> disks;
  std::vector<Vertex> constraints;

#if 0
  std::map<long,int>               Key;
#endif

public :

  // -------------------------------------------------------------------------
  Antichain() 
    : valid_(false) , straight_sweep_(false) , linear_space_(true)
  { }
  template < class DiskIterator ,class ConstraintIterator >
  Antichain(bool linear,DiskIterator first, DiskIterator last,
                               ConstraintIterator  firstc,ConstraintIterator lastc); 
  ~Antichain() { 
    this->detach();
    minimals_ccw_.detach();
    minimals_cw_.detach();
  }
  // -------------------------------------------------------------------------
  // Options when sweeping
  bool is_valid()               const { return valid_;          }
  bool is_straight()            const { return straight_sweep_; }
  void make_straight()                { straight_sweep_ = true; }
  bool linear_space()           const { return linear_space_;   }
  void set_linear_space(bool b)       { linear_space_ = b;      }
  // -------------------------------------------------------------------------
  // Identifying convex-hull vertices
  Face_handle infinite_face() const   { return infinite_face_;  }
  void create_infinite_face() { infinite_face_ = new Face;  }
  bool is_on_convex_hull(Vertex_handle v) const;
  // -------------------------------------------------------------------------
  // Iterator pairs for trversing the sink of the faces of the antichain
  // These vertices form the Greedy pseudo-triangulation
  Vertex_iterator       vertices_begin()
  { return Vertex_iterator(this,begin()); } 
  Vertex_const_iterator vertices_begin() const
  { return Vertex_const_iterator(this,begin()); } 
  Vertex_iterator   vertices_end()
  { return Vertex_iterator(this,end()); } 
  Vertex_const_iterator   vertices_end() const
  { return Vertex_const_iterator(this,end()); } 
  // -------------------------------------------------------------------------
  // Iterator pairs for trversing the sources of the faces of the antichain
  // These vertices form the dual Greedy pseudo-triangulation
  Vertex_cw_iterator       cw_vertices_begin()
  { return Vertex_cw_iterator(this,begin()); } 
  Vertex_cw_const_iterator cw_vertices_begin() const
  { return Vertex_cw_const_iterator(this,begin()); } 
  Vertex_cw_iterator   cw_vertices_end()
  { return Vertex_cw_iterator(this,end()); } 
  Vertex_cw_const_iterator   cw_vertices_end() const
  { return Vertex_cw_const_iterator(this,end()); } 
  // -------------------------------------------------------------------------
  Edge_iterator       edges_begin()       { return begin(); } 
  Edge_const_iterator edges_begin() const { return begin(); } 
  Edge_iterator       edges_end()         { return end();   } 
  Edge_const_iterator edges_end()   const { return end();   } 
  Edge_reverse_iterator edges_rbegin()    { return rbegin(); } 
  Edge_reverse_iterator edges_rend()      { return rend(); } 
  // -------------------------------------------------------------------------
  Face_iterator       faces_begin()     {return Face_iterator(this,begin());}
  Face_const_iterator faces_begin()const{return Face_iterator(this,begin());}
  Face_iterator       faces_end()       { return Face_iterator(this,end()); }
  Face_const_iterator faces_end()  const{ return Face_iterator(this,end()); }
  // -------------------------------------------------------------------------
  Minimals_iterator minimals_begin()    { return minimals_ccw_.begin(); }
  Minimals_iterator minimals_end()      { return minimals_ccw_.end(); }
  Minimals_iterator cw_minimals_begin() { return minimals_cw_.begin(); }
  Minimals_iterator cw_minimals_end()   { return minimals_cw_.end();   }
  Minimals_iterator minimals_begin(Ccw_traits) { return minimals_begin();    }
  Minimals_iterator minimals_end  (Ccw_traits) { return minimals_end();      }
  Minimals_iterator minimals_begin(Cw_traits)  { return cw_minimals_begin(); }
  Minimals_iterator minimals_end  (Cw_traits)  { return cw_minimals_end();   }
  // -------------------------------------------------------------------------
  // Iterator pair for linear sweep
  Linear_sweep_iterator       sweep_begin() 
  { return Linear_sweep_iterator(this); }
  Linear_sweep_const_iterator sweep_begin() const 
  { return Linear_sweep_const_iterator(this); }
  Linear_sweep_iterator       sweep_end()
  { return Linear_sweep_iterator(this,0); }
  Linear_sweep_const_iterator sweep_end() const 
  { return Linear_sweep_const_iterator(this,0); }
  // -------------------------------------------------------------------------
  // Testing minimality
  template < class Tr > 
  bool is_minimal(const Vertex_handle& v, Tr tr) const;
  bool is_minimal(const Vertex_handle& v) const
  { return is_minimal(v,Ccw_traits()); }
  template < class Tr > 
  bool is_xx_minimal(const Vertex_handle& v,Tr tr) const;
  bool is_right_minimal(const Vertex_handle& v) const
  { return is_xx_minimal(v,Right_ccw_traits()); }
  bool is_left_minimal(const Vertex_handle& v)  const
  { return is_xx_minimal(v,Left_ccw_traits()); }
  // -------------------------------------------------------------------------
  Vertex_handle pop_minimal(bool finite = false);
  template < class Tr > 
  bool is_swept_regular(const Vertex_handle& v , Tr tr) const;
  template < class Tr > 
  bool is_swept_constraint(const Vertex_handle& v , Tr tr) const;
  template < class Tr > 
  bool is_swept(const Vertex_handle& v , Tr tr) const;
  bool is_swept(const Vertex_handle& v) const 
  { return is_swept(v,Ccw_traits()); }
  // -------------------------------------------------------------------------
  // Adding a minimal
  void push_back_minimal(Vertex_handle v ,Ccw_traits);
  void push_back_minimal(Vertex_handle v ,Cw_traits);
  void push_back_minimal(Vertex_handle v) {push_back_minimal(v,Ccw_traits());}
  // -------------------------------------------------------------------------
  // Removing a minimal
  void erase_minimal(Vertex_handle v, Ccw_traits) { minimals_ccw_.erase(v); }
  void erase_minimal(Vertex_handle v, Cw_traits)  { minimals_cw_.erase(v);  }
  void erase_minimal(Vertex_handle v) { erase_minimal(v, Ccw_traits()); }
  // Depreciated - for backward compatibility
  template < class Tr > 
  void erase_minimal(Minimals_iterator v,Tr tr) { erase_minimal(&(*v),tr); }
  void erase_minimal(Minimals_iterator v) { erase_minimal(v,Ccw_traits()); }
  // -------------------------------------------------------------------------
  // Removing all minimals
  template < class Tr > void clear_minimals(Tr tr);
  void clear_minimals() { clear_minimals(Ccw_traits()); }
  // -------------------------------------------------------------------------
  // Sweeping - flipping operations
  template < class Tr > void sweep(Vertex_handle v, Tr tr);
  template < class Tr > void sweep_good(Vertex_handle v , Tr tr);
  template < class Tr > void sweep_all_minimals(Tr tr);
  template < class Tr > void sweep_good_all_minimals(Tr tr);
  void sweep(Vertex_handle v) { sweep(v,Ccw_traits()); }
  void sweep_all_minimals()   { sweep_all_minimals(Ccw_traits()); }
  // -------------------------------------------------------------------------
  void set_constraint(Vertex_handle v);
  void unset_constraint(Vertex_handle v);
  void add_constraint(Vertex_handle v);
  void remove_constraint(Vertex_handle v);
  // -------------------------------------------------------------------------
protected:
  // -------------------------------------------------------------------------
  // Compute the flipped bitangent phi(v)
  template < class Tr >
  Vertex_handle compute_phi(Vertex_handle v , Tr)  /*const*/;
  // -------------------------------------------------------------------------
  // Update the antichain while adding a constraint
  template < class Tr > void sweep_constraint(const Vertex_handle& v, const Tr&,bool opposite=false) const;
  template < class Tr > void sweep_regular   (const Vertex_handle& v, const Tr&,bool opposite=false) const;
  // -------------------------------------------------------------------------
  // Initialisation methods
  template<class DiskIterator , class ConstraintIterator> 
  void compute_graph(DiskIterator first, DiskIterator last,
                     ConstraintIterator  firstc,ConstraintIterator  lastc);
  struct compute_gr_aux;
  template<class DiskIterator , class ConstraintIterator> 
  void compute_gr(DiskIterator first, DiskIterator last,
                  ConstraintIterator  firstc,ConstraintIterator  lastc);
  template < class Tr > void compute_minimals(Tr tr);
  template < class Tr > void compute_vertices(Tr tr);
  template < class Tr > void fix_extreme_edges(Tr tr) const;
//   template < class Tr > void glue_ccw_cw(Tr tr) const;
  void glue();
  template < class Tr > void initialize_convex_hull(Tr tr) const;
  // -------------------------------------------------------------------------
  // Method used during intialization, to compute candidates during the
  // Bentley-Ottmann rotational sweep.
  //   template<class Tr> bool is_candidate(const Face_handle& fa) const;
  template < class Tr >
  bool
  is_candidate(const Face_handle& f, const Tr&) const
  {
    typename Tr::Dl dl; typename Tr::Ur ur; typename Tr::Sup sup;
    if (f == 0 || f->top_edge() == 0 || f->bottom_edge() == 0) return false;
    if (sup(f) != 0 && sup(f)->is_constraint()) return is_minimal(sup(f),Tr());
    if (f->top_edge()->object() == 0 || f->bottom_edge()->object() == 0)
      return false;
    return (ur(f->bottom_edge()) == f && dl(f->top_edge()) == f);
  }

#if 0

  template<class OutputStream> OutputStream& node(OutputStream &o,Edge_handle e) {
    o<<"edge_node_"<<std::hex<<((size_t) e)<<std::dec;
    return o;
  }
  template<class OutputStream> OutputStream& node(OutputStream &o,Face_handle f,bool b=false) {
    if (f==0) f=infinite_face_;
    if (f!=infinite_face_)
      o<<"face_node_"<<std::hex<<((size_t) f)<<std::dec;
    else 
      o<<"face_node_inf_"<<(b?"below":"above")<<std::hex<<((size_t) f)<<std::dec;

    return o;
  }

  template<class OutputStream> void graphviz_print(OutputStream& o) {
    
    o<<"digraph {ranksep=0.2\n";
    o<<"subgraph{ordering=in;\n";
    for (Edge_iterator i=edges_begin();i!=edges_end();++i) {
      if (i->object()&&i->sign()) node(o,&*i)<<"[label=\""<<std::hex<<((size_t) &*i)<<std::dec<<"\",shape=invtriangle,width=0.3,height=0.2]"<<";\n";
    }
    o<<"\n}\nsubgraph{ordering=out;\n";
    for (Edge_iterator i=edges_begin();i!=edges_end();++i) {
      if (i->object()&&!i->sign()) node(o,&*i)<<"[label=\""<<std::hex<<((size_t) &*i)<<std::dec<<"\",shape=triangle,width=0.3,height=0.2]"<<";\n";
    }
    o<<"\n}\n";


    for (Face_iterator i=faces_begin();i!=faces_end();++i) {
      node(o,&*i)<<"[label=\""<<std::hex<<((size_t) &*i)<<std::dec<<"\",shape=rectangle]"<<";\n";
    }
    node(o,infinite_face_,true)<<"[label=\"inf-"<<std::hex<<((size_t)infinite_face_)<<std::dec<<"\",shape=point]"<<";\n";
    node(o,infinite_face_,false)<<"[label=\"inf+"<<std::hex<<((size_t)infinite_face_)<<std::dec<<"\",shape=point]"<<";\n";
      

    for (Edge_iterator i=edges_begin();i!=edges_end();++i) {
      if (!i->object()) {
        if (i->sign()) node(node(o,i->ul()?i->ul():i->ur())<<"->",&*i)<<"[arrowhead=none];"; else node(node(o,&*i)<<"->",i->dl()?i->dl():i->dr())<<"[arrowhead=none];";
      } else {
        node(node(o,&*i)<<"->",i->dl(),true)<<"[arrowhead=none];";
        if (i->sign()) node(node(o,i->ul(),false)<<"->",&*i)<<"[arrowhead=none];";
        else node(node(o,&*i)<<"->",i->dr(),true)<<"[arrowhead=none];";
        node(node(o,i->ur(),false)<<"->",&*i)<<"[arrowhead=none];";
      }
      o<<"\n";
    }
    o<<"}\n";
  }

public:void graphviz_view() {
  std::ostringstream fn;
  fn<<"antichain"<<getpid();
  std::ostringstream fnd;
  fnd<<fn.str()<<".dot";
  std::ostringstream fnp;
  fnp<<fn.str()<<".ps";
  std::ostringstream cmd;
  cmd<<"(dot -Tps -o "<<fnp.str()<<" "<<fnd.str()<<";gv "<<fnp.str()<<";rm "<<fnd.str()<<" "<<fnp.str()<<")&";
  std::ofstream out(fnd.str().c_str());
  graphviz_print(out);
  out.close();
  system(cmd.str().c_str());
}
#endif
  // -------------------------------------------------------------------------
#if 0
public:


  void print(Vertex_handle v) {
    //std::cout << v << " " << flush;
    if (v == 0) return; 
    std::cout << *v << " , (";
    if (v->is_left_xx()) std::cout << "+"; else std::cout << "-";
    std::cout << Key[long(v->source_object())] << ",";
    if (v->is_xx_left()) std::cout << "+"; else std::cout << "-";
    std::cout << Key[long(v->target_object())] << ")" ;
    /*
      std::cout << " minLR(" << flush 
      << is_xx_minimal(v,Left_ccw_traits()) << "," << flush 
      << is_xx_minimal(v,Right_ccw_traits()) << ")";
      std::cout << " cw_minLR(" << flush 
      << is_xx_minimal(v,Left_cw_traits()) << "," << flush 
      << is_xx_minimal(v,Right_cw_traits()) << ")";
      if (is_swept(v)) std::cout << " swept";
      else std::cout << " not swept";
    */
  }
  void print(Edge_handle e, bool coord = true) {
    std::cout << e << " " << flush;
    if (e == 0) return; 
    if (e->object() == 0) {
      if (e == e->sup()->target_cusp_edge())
        std::cout << "-" << Key[long(e->sup()->target_object())] << "c";
      else std::cout << "+" << Key[long(e->sup()->source_object())] << "c";
    }
    else {
      if (e->sign()) std::cout << "+"; else std::cout << "-";
      std::cout << Key[long(e->object())] ;
    }
    std::cout << " {";
    if (e->sign()) 
      std::cout << e->dl() << "," 
                << e->ur() << "," 
                << e->ul();
    else std::cout << e->dl() << "," 
                   << e->dr() << "," 
                   << e->ul();
    std::cout << "}";
    if (coord) {
      std::cout << " [";
      std::cout << *e->begin() << "," << *--e->end() << "]";
    }
  }
  void print(Face_handle f) {
    std::cout << f << " " << flush;
    if (f == 0) return; 
    std::cout << "[" ; 
    print(f->bottom_edge()); std::cout << ",";
    print(f->top_edge());    std::cout << "]";
    std::cout << " inf = " << f->inf() << " , sup = " << f->sup();
  }
  template < class Tr >
  void print_left_minimal(Vertex_handle v ,Tr tr);
  template < class Tr >
  void print_right_minimal(Vertex_handle v ,Tr tr);
#endif
  // -------------------------------------------------------------------------
public:

#if 0
  // Debugging stuff specialized for circles.
  void view_disks() {
    double ru=900;
    if (fork()==0) {
      int zero=0;
      QApplication app(zero,(char**)0);
      Qt_widget* w;
      w = new CGAL::Qt_widget();
      app.setMainWidget( w );
      w->resize(900, 900);
      w->set_window(0, to_double(ru),0, to_double(ru));
      w->show();
      w->lock();
      {
        int j=0;
        *w<<GREEN;
        for (typename std::vector<Disk_handle>::iterator first=disks.begin();first!=disks.end();first++,j++) {
          *w<<**first;
          typename CGAL::Point_2<CGAL::Simple_cartesian<double> > A(to_double((*first)->center().x()),to_double((*first)->center().y()));
          //           std::cout<<A<<"\n";
          std::ostringstream n;
          n<<j;
          w->get_painter().drawText(w->x_pixel(A.x())-15,w->y_pixel(A.y()),n.str());
        }  }
      w->unlock();;
      app.exec();
      exit(0);
    }
  }
  void view_bit(Vertex_handle v) {
    double ru=900;
    if (fork()==0) {
      int zero=0;
      QApplication app(zero,(char**)0);
      Qt_widget* w;
      w = new CGAL::Qt_widget();
      app.setMainWidget( w );
      w->resize(900, 900);
      w->set_window(0, to_double(ru),0, to_double(ru));
      w->show();
      w->lock();
      {
        int j=0;
        *w<<GREEN;
        for (typename std::vector<Disk_handle>::iterator first=disks.begin();first!=disks.end();first++,j++) {
          *w<<**first;
          typename CGAL::Point_2<CGAL::Simple_cartesian<double> > A(to_double((*first)->center().x()),to_double((*first)->center().y()));
          std::ostringstream n;
          n<<j;
          w->get_painter().drawText(w->x_pixel(A.x())-15,w->y_pixel(A.y()),n.str());
        }  
      }
      *w<<BLACK;
      *w<<static_cast<typename Gt::Segment_2>(*v);
      *w<<typename Gt::R::Circle_2(v->source(),100);
      w->unlock();
      app.exec();
      exit(0);
    }
  }
  template <class ii> void view_pt(ii i,ii ie,Vertex_handle v=0) {
    double ru=900;
    if (fork()==0) {
      int zero=0;
      QApplication app(zero,(char**)0);
      Qt_widget* w;
      w = new CGAL::Qt_widget();
      app.setMainWidget( w );
      w->resize(900, 900);
      w->set_window(0, to_double(ru),0, to_double(ru));
      w->show();
      w->lock();
      {
        *w<<GREEN;
        int j=0;
        for (typename std::vector<Disk_handle>::iterator first=disks.begin();first!=disks.end();first++,j++) {
          *w<<(*(*first));
          typename CGAL::Point_2<CGAL::Simple_cartesian<double> > A(to_double((*first)->center().x()),to_double((*first)->center().y()));
          std::ostringstream n;
          n<<j;
          w->get_painter().drawText(w->x_pixel(A.x())-15,w->y_pixel(A.y()),n.str());
        }  }
      *w<<BLACK;
      for (;i!=ie;++i) {
        *w<<static_cast<typename Gt::Segment_2>(*i);
        *w<<typename Gt::R::Circle_2((i)->source(),100);
      }
      if (v) {
        *w<<BLUE;
        *w<<static_cast<typename Gt::Segment_2>(*v);
        *w<<typename Gt::R::Circle_2(v->source(),100);        
        *w<<RED;
        *w<<static_cast<typename Gt::Segment_2>(*v->sup()->sup());
        *w<<typename Gt::R::Circle_2(v->sup()->sup()->source(),100);        
      }
      w->unlock();
      app.exec();
      exit(0);
    }
  }
#endif
};

// -----------------------------------------------------------------------------
// Computing the Antichain using a Bentley-Ottman sweep

template < class Gtr_ , class It , class Flip >
template < class DiskIterator , class ConstraintIterator >
Antichain<Gtr_,It,Flip>::
Antichain(bool linear,DiskIterator first, DiskIterator last ,
	  ConstraintIterator  firstc,ConstraintIterator  lastc)
{
  // -------------------------------------------------------------------------
  // Default values : linear space topological sweep
  straight_sweep_ = false;
  linear_space_ = linear;

  for (DiskIterator i=first;i!=last;++i) disks.push_back(&*i);
  for (ConstraintIterator i=firstc;i!=lastc;++i) {
    constraints.push_back(Vertex(i->type(),&*disks.begin()[i->source()],
                                 &*disks.begin()[i->target()]));
  }

  if (first == last) { valid_ = false; return; }
  compute_gr(first,last,
             constraints.begin(),constraints.end());

  if (!is_valid()) return;
  // -------------------------------------------------------------------------
  compute_vertices(Ccw_traits());
  compute_vertices(Cw_traits());

  fix_extreme_edges(Ccw_traits());
  fix_extreme_edges(Cw_traits());
  // -------------------------------------------------------------------------
  //     glue_ccw_cw(Ccw_traits());
  //     glue_ccw_cw(Cw_traits());

  glue();
  // -------------------------------------------------------------------------
  initialize_convex_hull(Ccw_traits());
  initialize_convex_hull(Cw_traits());
  // -------------------------------------------------------------------------
  compute_minimals(Ccw_traits());
  compute_minimals(Cw_traits());

  for (Face_iterator f=faces_begin();f!=faces_end();++f) 
    f->set_in_initial_antichain(true);

  // Create symetrical faces for the initial antichain.
  if (!linear_space()) for (Edge_iterator e=begin();e!=end();++e) {
    if (!e->object()) continue;
    Face_handle fs[3]={e->dl(),e->sign()?e->ur():e->dr(),e->ul()};
    for (int i=0;i<3;++i) {
      if (fs[i]==infinite_face()) continue;
      CGAL_assertion(fs[i]->inf()&&fs[i]->sup());
      Vertex_handle v1=fs[i]->inf()->pi();
      Vertex_handle v2=fs[i]->sup()->pi();
      Face_handle f=v1->sup();
      if (!f) {
        if (v1->is_constraint()) {
          if (fs[i]==v1->source_cusp_face()) 
            f=v2->target_cusp_face();
          else
            f=v2->source_cusp_face();
        } else {
          f=new Face();
          v1->set_sup(f);
          f->set_inf(v1);
          v2->set_inf(f);
          f->set_sup(v2);          
        }
      }
      Edge_handle eo;
      if (e.operator->()==e->inf()->ccw_source_edge()) 
        eo=e->inf()->pi()->ccw_target_edge();
      else
        eo=e->inf()->pi()->ccw_source_edge();
      Face_handle f1=eo->dl();
      Face_handle f2=eo->sign()?eo->ur():eo->dr();
      Face_handle f3=eo->ul();
      switch (i) {
      case 0: if (e->sign()) f3=f; else f2=f; break;
      case 1: if (e->sign()) f1=f; else f3=f; break;
      case 2: if (e->sign()) f2=f; else f1=f; break;
      }
      eo->set_adjacent_faces(f1,f2,f3);
    }
  }
}



template < class Gtr_ , class It , class Flip > 
struct Antichain<Gtr_,It,Flip>::compute_gr_aux {


  struct constraint_entry;
  typedef std::vector<constraint_entry*> VV;

  template<class chain> struct wrapper {
    //to x-compare a chain with an object along the sweep-line
    struct Less_chain {
      // position of the sweep line
      Edge_handle edg;
      typename Gtr_::Is_upward_directed is_upward_directed;
      Less_bitangent<Gtr_> lb;
      chain* left_neighbour;
      chain* right_neighbour;
      chain* currently_inserted;
      Less_chain() : edg(0),left_neighbour(0), right_neighbour(0), 
                     currently_inserted(0) {}
      void set_currently_inserted(chain *a) {
        currently_inserted=a;
      }
      void set_left_neighbour(chain *a) {
        left_neighbour=a;
      }
      void set_right_neighbour(chain *a) {
        right_neighbour=a;
      }
      void set_pos(Edge_handle e) {
        edg=e;
      }
      bool leftp(Disk_handle d2) const {
        Bitangent_2 b(Bitangent_2::RL,edg->object(),d2);
        return is_upward_directed(b); 
      }
      Comparison_result auxl (Disk_handle d,Disk_handle d1, 
                              Vertex_handle v2, Bitangent_2* blr) const {
        if (!v2||lb(*v2,*blr)) return SMALLER; else {
          if (v2->is_right_right()) {
            Bitangent_2 b2(Bitangent_2::LR,d,v2->target_object());
            if (is_upward_directed(b2)&&lb(b2,*v2))
              return SMALLER; else return LARGER;
          } else {
            Bitangent_2 b2(Bitangent_2::RL,d,v2->target_object());
            if (is_upward_directed(b2)&&lb(*v2,b2))
              return LARGER; else return SMALLER;
          } 
        }
      }

      Comparison_result auxr (Disk_handle d,Disk_handle d1, 
                              Vertex_handle v2, Bitangent_2* brl) const {
        if (!v2||lb(*brl,*v2)) return LARGER; else {
          if (v2->is_left_left()) {
            Bitangent_2 b2(Bitangent_2::RL,d,v2->target_object());
            if (is_upward_directed(b2)&&lb(*v2,b2))
              return LARGER; else return SMALLER;
          } else {
            Bitangent_2 b2(Bitangent_2::LR,d,v2->target_object());
            if (is_upward_directed(b2)&&lb(b2,*v2))
              return SMALLER; else return LARGER;
          }          
        }
      }
      
      Comparison_result operator()(const chain *a, const chain * b) const {
        CGAL_precondition(a==currently_inserted);
        if (a==b) return EQUAL;
        
        if (a->t==chain::REGULAR_CHAIN) {
          // should not occur, except if insert_before and insert_after
          // perform assertion checking.
          if (b==left_neighbour) return LARGER;
          if (b==right_neighbour) return SMALLER;
          CGAL_assertion(false);
        }
        // a is the chain of the disk being inserted
        if (b->t==chain::OBJECT) {
          if (leftp(b->d)) return LARGER; else return SMALLER;
        } else {
          switch (b->st) {
          case chain::Vtx: {
            Vertex_handle v=b->iplb[0].v->vtx;
            if (v->target_object()==a->d) {
              if (v->is_xx_left()) return SMALLER; else return LARGER;
            }
            if (v->is_xx_left()) {
              Bitangent_2 bt(Bitangent_2::RL,a->d,v->target_object());
              if (is_upward_directed(bt)&&lb(*v,bt)) return LARGER; else
                return SMALLER;
            } else {
              Bitangent_2 bt(Bitangent_2::LR,a->d,v->target_object());
              if (is_upward_directed(bt)&&lb(bt,*v)) return SMALLER; else
                return LARGER;
            }
          }
          case chain::VtxArcVtx: 
            {
              constraint_entry *v1=b->iplb[0].v;
              constraint_entry *v2=b->iplb[1].v;
              if (v2&&v2->vtx->target_object()==a->d) {
                if (v2->vtx->is_xx_left()) return SMALLER; else return LARGER;
              }

              Disk_handle d=
                v1?v1->vtx->target_object():v2->vtx->source_object();

              if ((v1?v1->vtx->is_xx_left():v2->vtx->is_left_xx())) {
                Bitangent_2 brl(Bitangent_2::RL,a->d,d);
                if (!is_upward_directed(brl)) return SMALLER;
                if (!v1||lb(*(v1->vtx),brl)) 
                  return auxr(a->d,d,v2?v2->vtx:0,&brl); 
                else 
                  return SMALLER;
              } else {
                Bitangent_2 blr(Bitangent_2::LR,a->d,d);
                if (!is_upward_directed(blr)) return LARGER;
                if (!v1||lb(blr,*(v1->vtx))) 
                  return auxl(a->d,d,v2?v2->vtx:0,&blr); 
                else 
                  return LARGER;
              }
            }
          case chain::VtxArcVtxArcVtx: 
            {
              constraint_entry *v1=b->iplb[0].v;
              Vertex_handle v2=b->iplb[1].v->vtx;
              constraint_entry *v3=b->iplb[2].v;
              if (v3&&v3->vtx->target_object()==a->d) {
                if (v3->vtx->is_xx_left()) return SMALLER; else return LARGER;
              }
              Disk_handle d1=v2->source_object();
              Disk_handle d2=v2->target_object();

              if (v2->is_left_right()) {
                Bitangent_2 brl1(Bitangent_2::RL,a->d,d1);
                if (is_upward_directed(brl1)) {
                  Bitangent_2 blr2(Bitangent_2::LR,a->d,d2);
                  if (is_upward_directed(blr2)) {
                    if (lb(*v2,blr2)) return LARGER; else {
                      if (lb(brl1,*v2)) {
                        if (v1&&lb(brl1,*(v1->vtx))) return SMALLER; else return LARGER;
                      } else return auxl(a->d,d2,v3?v3->vtx:0,&blr2);
                    }
                  } else return LARGER;                
                } else return SMALLER;
              } else {
                Bitangent_2 blr1(Bitangent_2::LR,a->d,d1);
                if (is_upward_directed(blr1)) {
                  Bitangent_2 brl2(Bitangent_2::RL,a->d,d2);
                  if (is_upward_directed(brl2)) {
                    if (lb(brl2,*v2)) return SMALLER; else {
                      if (lb(*v2,blr1)) {
                        if (v1&&lb(*(v1->vtx),blr1)) return LARGER; else return SMALLER;
                      } else return auxr(a->d,d2,v3?v3->vtx:0,&brl2);
                    }
                  } else return SMALLER;
                } else return LARGER;
              }
            }
          default:
            return EQUAL; // silence silly compiler
          }
        }
      }
    };
    typedef Multiset<chain *, Less_chain>  Ystructure; 
  };
  struct chain {

    typedef typename wrapper<chain>::Ystructure Ystructure;

    enum window_type {
      Vtx,VtxArcVtx,VtxArcVtxArcVtx
    };
    window_type st;

    struct IPLB : VC_in_place_list_base<IPLB> {
      chain * super; // points to the chain
      constraint_entry * v; // points to a vertex
      IPLB(chain * c) : super(c),v(0) {};
      IPLB() {super=0;v=0;}
    };
    
    IPLB iplb[3]; // one IPLB for each vertex in the window

    enum chain_type {
      OBJECT, REGULAR_CHAIN
    };
    
    Face_handle left;
    Face_handle right;

    Disk_handle d;
    chain_type t;
    int status; // -1 not met yet, 0 intersected, 1 past history
    typename Ystructure::iterator position; // position in the ystructure
    chain(Disk_handle d) {
      iplb[0]=iplb[1]=iplb[2]=IPLB(this);
      status=-1;
      st=VtxArcVtx;
      t=OBJECT;
      this->d=d;
    }
    chain(constraint_entry* v){
      iplb[0]=iplb[1]=iplb[2]=IPLB(this);
      status=-1;
      t=REGULAR_CHAIN;
      st=VtxArcVtx;
      left=0;right=0;d=0;
    }
    chain& operator=(const chain& c) {
      st=c.st;
      for (int i=0;i<3;i++) {
        iplb[i]=c.iplb[i];
        iplb[i].super=this;
      }
      left=c.left; right=c.right;
      d=c.d;
      t=c.t;
      status=c.status;
      position=c.position;
      return *this;
    }
    chain() {
      iplb[0]=iplb[1]=iplb[2]=IPLB(this);
    };
  };    

  typedef typename wrapper<chain>::Less_chain Less_chain;
  typedef typename wrapper<chain>::Ystructure Ystructure;

  struct disk_entry {
    Disk_handle d;
    
    struct side {
      VV constraints;
    };
    
    side lhs; // left hand side constraints
    side rhs; // right hand side constraints

    chain ch;

    Edge_handle pos,neg;

    int status;

    disk_entry() {}

    disk_entry(Disk_handle dd) {
      d=dd;
      pos=new Edge(true,d);
      neg=new Edge(false,d);
      //       pos->set_opposite(neg);
      ch=chain(d);
      status=-1;
    }
  };

  struct constraint_entry {
    Vertex_handle vtx;
    // the first constraint that leaves the object this constraint is
    // entering above the entry point
    constraint_entry* target_next_out;
    // positions in the constraint lists of the objects adjacent to this
    // constraint
    typename VV::iterator source_cycle;
    typename VV::iterator target_cycle;
    chain source_chain; // the chain this bitangent is the first
                        // bitangent of
    typedef In_place_list<typename chain::IPLB,false> IPL;
    IPL ipl[3]; // ipl[i] contains the list of the chains whose window has
                // this constraint as ith vertex.
  };

  ;
  struct Less_constraintl {
    Less_bitangent<Gtr_> lb;
    bool operator() (const constraint_entry* a, const constraint_entry* b) {
      return lb(*(b->vtx),*(a->vtx));
    }
  };
  struct Less_constraintr {
    Less_bitangent<Gtr_> lb;
    bool operator() (const constraint_entry* a, const constraint_entry* b) {
      return lb(*(a->vtx),*(b->vtx));
    }
  };
};



template < class Gtr_ , class It , class Flip >
template < class DiskIterator , class ConstraintIterator >
void
Antichain<Gtr_,It,Flip>::
compute_gr(DiskIterator first, DiskIterator last,
           ConstraintIterator  firstc,ConstraintIterator  lastc) {
  valid_=true;

  typedef typename compute_gr_aux::constraint_entry constraint_entry;
  typedef typename compute_gr_aux::VV VV;
  typedef typename compute_gr_aux::chain chain;
  typedef typename compute_gr_aux::Less_chain Less_chain;
  typedef typename compute_gr_aux::Ystructure Ystructure;
  typedef typename compute_gr_aux::disk_entry disk_entry;
  typedef typename compute_gr_aux::Less_constraintl Less_constraintl;
  typedef typename compute_gr_aux::Less_constraintr Less_constraintr;

  

  typedef std::vector<Edge_handle> Xstructure;
  Xstructure XE; // the events


  typedef std::vector<chain> CV;



  Less_constraintl lctl;
  Less_constraintr lctr;


  typename Ccw_traits::Set_adjacent_faces_one_to_one set_adjacent_faces;


  typedef std::map<Disk_handle,disk_entry> DDEM;
  DDEM dem;
  for (DiskIterator i=first;i!= last;++i) {
    typename DDEM::iterator j=
      dem.insert(typename DDEM::value_type(&*i,disk_entry(&*i))).first;
    push_back(*(j->second.neg));
    push_back(*(j->second.pos));
    set_adjacent_faces(j->second.pos,0,new Face,new Face);
    j->second.pos->ul()->set_front_view(j->second.pos);
    j->second.pos->ur()->set_back_view(j->second.pos);
    set_adjacent_faces(j->second.neg,0,0,new Face);
  }




  typename Gtr_::Is_upward_directed is_upward_directed;
  typedef std::map<Vertex_handle,constraint_entry> VCEM;
  VCEM vcem;
  



  for (;firstc!=lastc;++firstc) {
    Vertex_handle c=&*firstc;
    if (!c->is_constraint())       set_constraint(&(*c));
    if (!c->pi()->is_constraint()) set_constraint(c->pi());
    if (!is_upward_directed(*c)) c=c->pi();
    typename DDEM::iterator desi=dem.find(c->source_object());
    CGAL_precondition(desi!=dem.end());
    disk_entry * des=&(desi->second);
    typename DDEM::iterator deti=dem.find(c->target_object());
    CGAL_precondition(deti!=dem.end());
    disk_entry * det=&(deti->second);
    typename VCEM::iterator i=
      vcem.insert(typename VCEM::value_type(c,constraint_entry())).first;
    i->second.vtx=c;

    ((c->is_left_xx())?&des->rhs:&des->lhs)->constraints.push_back(&i->second);
    ((c->is_xx_left())?&det->rhs:&det->lhs)->constraints.push_back(&i->second);

    push_back(*c->target_cusp_edge());
    push_back(*c->source_cusp_edge());
  }




  // sorting bottom-up the constraint lists for each object, and splicing
  // the edges
  for (typename DDEM::iterator dei=dem.begin();dei!=dem.end();++dei) {
    std::sort(dei->second.lhs.constraints.begin(),
              dei->second.lhs.constraints.end(),lctl);
    std::sort(dei->second.rhs.constraints.begin(),
              dei->second.rhs.constraints.end(),lctr);

    typename Ccw_traits::Splice splice_ccw;
    typename Cw_traits::Splice  splice_cw;

    Edge_handle e;
    e=dei->second.pos;
    for (typename VV::iterator i=dei->second.rhs.constraints.begin();
         i!=dei->second.rhs.constraints.end();
         ++i,e=e->sup()->ccw_edge(e->object())) {
      splice_ccw(e,(*i)->vtx);
    }
    e=dei->second.pos;
    for (typename VV::iterator i=dei->second.lhs.constraints.begin();
         i!=dei->second.lhs.constraints.end();
         ++i,e=e->inf()->cw_edge(e->object())) {
      splice_cw(e,(*i)->vtx->pi());
    }
    e=dei->second.neg;
    for (typename VV::reverse_iterator i=dei->second.lhs.constraints.rbegin();
         i!=dei->second.lhs.constraints.rend();
         ++i,e=e->sup()->ccw_edge(e->object())) {
      splice_ccw(e,(*i)->vtx);
    }
    e=dei->second.neg;
    for (typename VV::reverse_iterator i=dei->second.rhs.constraints.rbegin();
         i!=dei->second.rhs.constraints.rend();
         ++i,e=e->inf()->cw_edge(e->object())) {
      splice_cw(e,(*i)->vtx->pi());
    }

    constraint_entry * prev_out_constraint=0;
    chain * prev_out_chain=&dei->second.ch;
    

    if (!dei->second.rhs.constraints.empty())
      for (typename VV::iterator i=--dei->second.rhs.constraints.end();1
             ;--i) {
        if ((*i)->vtx->source_object()==dei->second.d) {
          (*i)->source_chain=chain(*i);
          Edge_handle e=(*i)->vtx->source_cusp_edge();
          Face_handle f=e->ul();
          f->set_back_view((*i)->vtx->ccw_source_edge());
          prev_out_chain->right=f;
          prev_out_chain=&((*i)->source_chain);
          prev_out_chain->left=f;
          prev_out_constraint=*i;
          (*i)->source_cycle=i;
        }
        else {
          (*i)->target_next_out=prev_out_constraint;
          (*i)->target_cycle=i;
        }
        if (i==dei->second.rhs.constraints.begin()) break;
      }
    prev_out_chain->right=dei->second.pos->ur();

    prev_out_constraint=0;
    prev_out_chain=&(dei->second.ch);
    if (!dei->second.lhs.constraints.empty())
      for (typename VV::iterator i=--dei->second.lhs.constraints.end();1
             ;--i) {
        if ((*i)->vtx->source_object()==dei->second.d) {
          (*i)->source_chain=chain(*i);
          Edge_handle e=(*i)->vtx->source_cusp_edge();
          Face_handle f=e->ur();
          f->set_front_view((*i)->vtx->cw_source_edge());
          prev_out_chain->left=f;
          prev_out_chain=&((*i)->source_chain);
          prev_out_chain->right=f;
          prev_out_constraint=*i;
          (*i)->source_cycle=i;
        }
        else {
          (*i)->target_next_out=prev_out_constraint;
          (*i)->target_cycle=i;
        }
        if (i==dei->second.lhs.constraints.begin()) break;
      }
    prev_out_chain->left=dei->second.pos->ul();

    XE.push_back(dei->second.pos);
    XE.push_back(dei->second.neg);
  }
  std::sort(XE.begin(),XE.end(),Less_edge_handle<Gtr_>());
  

  Ystructure YE; // the chains intersected by the sweep line
  Less_chain& LC=YE.key_comp();
 
  struct utils {
    VCEM& vcem;
    DDEM& dem;
    Xstructure& XE;
    Ystructure& YE;
    Less_chain& LC;
    utils(VCEM& v,DDEM& d, Xstructure& xe, Ystructure& ye, Less_chain& lc) :
      vcem(v), dem(d), XE(xe), YE(ye), LC(lc) {};
    typename Ccw_traits::Set_adjacent_faces_one_to_one set_adjacent_faces;

    // installs v as the ith vertex in c's window
    void set_chain_constraint(constraint_entry *v,chain *c,int i) {
      if (v) {
        for (int j=0;j<3;j++) if (c->iplb[j].v==v) {
          if (i==j) return;
          c->iplb[j].v=0;
          v->ipl[j].erase(c->iplb+j);
        }
      }
      if (c->iplb[i].v) c->iplb[i].v->ipl[i].erase(c->iplb+i);
      c->iplb[i].v=v;
      if (v) v->ipl[i].push_front(c->iplb[i]);
    }

    // called when the bottom of the window passes the first cusp
    void remove_chain(chain *ch) {
      if (ch->status!=0) return;
      
      constraint_entry *v=ch->iplb[0].v;
      Edge_handle cusp=v->vtx->target_cusp_edge();
      Face_handle fdead,fcusp;
      if (v->vtx->is_xx_right()) {
        fcusp=cusp->dr();
        fdead=ch->right;
        typename Ystructure::iterator a=ch->position;
        ++a;
        (*a)->left=ch->left;
      } else {
        fcusp=cusp->dl();
        fdead=ch->left;
        typename Ystructure::iterator a=ch->position;
        --a;
        (*a)->right=ch->right;
      }
      Edge_handle edead=fdead->bottom_edge();
      if (edead) {
        if (edead->sign()) {
          set_adjacent_faces(edead,edead->dl(),
                             edead->ur()==fdead?fcusp:edead->ur(),
                             edead->ul()==fdead?fcusp:edead->ul());
        } else {
          set_adjacent_faces(edead,edead->dl(),edead->dr(),fcusp);
        }
        delete fdead;
      }
      YE.erase(ch->position);
      ch->position=YE.end();
      set_chain_constraint(0,ch,0);
      set_chain_constraint(0,ch,1);
      set_chain_constraint(0,ch,2);
    }
    // makes ch a VtxArcVtx with v1 and v2; removes ch if its window passes
    // a cusp; returns true iff ch was not removed
    bool f1(chain * ch,constraint_entry *v1,constraint_entry* v2) {
      if (v1==0||ch==&v1->source_chain) {
        set_chain_constraint(v1,ch,0);
        set_chain_constraint(v2,ch,1);
        set_chain_constraint(0,ch,2);
        ch->st=chain::VtxArcVtx;
        return true;
      } else {
        remove_chain(ch
                     );
        ch->status=1;
        return false;
      }
    }
    // makes ch a Vtx; removes ch if its window passes a cusp; returns true
    // iff ch was not removed
    bool f2(chain *ch,constraint_entry * v) {
      if (ch==&v->source_chain) {
        set_chain_constraint(v,ch,0);
        set_chain_constraint(0,ch,1);
        set_chain_constraint(0,ch,2);
        ch->st=chain::Vtx;
        return true;
      } else {
        remove_chain(ch
                     );
        ch->status=1;
        return false;
      }
    }
    // check whether a VtxArcVtx must be extented into a VtxArcVtxArcVtx,
    // possibly advancing the bottom of the window and removing the chain
    bool f3(chain * ch) {
      constraint_entry * v2=ch->iplb[1].v;
      if (!v2||v2->vtx->is_external()) return true;
      Disk_handle d3=v2->vtx->target_object();
      disk_entry * de3=&(dem[d3]);
      if (de3->status==0) {
        constraint_entry* v3=v2->target_next_out; disk_entry * de4;
        bool b; Edge_handle e;
        if (!v3||v3->vtx->is_external()) goto extend;
        de4=&(dem[v3->vtx->target_object()]);
        if (de4->status!=0) goto extend;
        e=LC.edg;
        LC.set_pos(de4->pos);
        b=LC.leftp(v2->vtx->source_object());        
        LC.set_pos(e);
        if (v2->vtx->is_right_left()) b=!b;
        if (b) {
          if (f1(ch,v2,v3)) {
            set_chain_constraint(v3->target_next_out,ch,3);
            ch->st=chain::VtxArcVtxArcVtx;
            return true;
          } else return false;
        }
      extend:
        set_chain_constraint(v2->target_next_out,ch,2);
        ch->st=chain::VtxArcVtxArcVtx;
        return true;
      }
      return true;
    }
  };

  

  utils U(vcem,dem,XE,YE,LC);
  


  Face lower_outer_face;
  Face_handle upper_outer_face=&lower_outer_face;

  for (typename Xstructure::iterator ei=XE.begin();ei!=XE.end();ei++) {
    Edge_handle e=*ei;

    disk_entry *de;
    {
      typename DDEM::iterator i=dem.find(e->object());
      de=&(i->second);
    }

    LC.set_pos(e);
    
    if (e->sign()) {
      // finding where to insert the new chains
      LC.set_currently_inserted(&de->ch);
      typename Ystructure::iterator rc=YE.upper_bound(&(de->ch));
      LC.set_currently_inserted(0);
      typename Ystructure::iterator lc=YE.end();
      Face_handle f;
      if (rc==YE.end()&&YE.empty()) {
        f=upper_outer_face;
        upper_outer_face=0;
      } else if (rc==YE.begin()) {
        f=(*rc)->left;
      } else {
        lc=rc;--lc;
        f=(*lc)->right;
      }

      if (rc!=YE.end()) {
        LC.set_right_neighbour(*rc);
        (*rc)->left=e->ur();
      } else LC.set_right_neighbour(0);
      if (lc!=YE.end()) {
        LC.set_left_neighbour(*rc);
        (*lc)->right=e->ul();
      } else LC.set_left_neighbour(0);

      e->ur()->set_front_view(f->front_view());
      e->ul()->set_back_view(f->back_view());

      set_adjacent_faces(e,f,e->ur(),e->ul());
      typename Ystructure::iterator k=YE.insert_before(rc,&de->ch);
      LC.set_currently_inserted(0);
      de->ch.position=k;
      de->ch.status=0;
      de->status=0;
      LC.set_right_neighbour(&de->ch);
      for (typename VV::reverse_iterator i=de->lhs.constraints.rbegin();
           i!=de->lhs.constraints.rend();++i) {
        // a new chain
        if ((*i)->vtx->source_object()==de->d) {
          U.f1(&(*i)->source_chain,0,(*i));
          U.f3(&(*i)->source_chain);
          LC.set_currently_inserted(&(*i)->source_chain);
          k=YE.insert_before(k,&(*i)->source_chain);
          LC.set_currently_inserted(0);
          LC.set_right_neighbour(&(*i)->source_chain);
          (*i)->source_chain.position=k;
          (*i)->source_chain.status=0;
        } else { // a chain whose window ends on i
          constraint_entry * ve=(*i);
          Vertex_handle v=ve->vtx;

          typename constraint_entry::IPL::iterator tmp;
          // the window was a single constraint
          for (typename constraint_entry::IPL::iterator j=ve->ipl[0].begin();
               j!=ve->ipl[0].end();j++) {
            U.set_chain_constraint(ve->target_next_out,j->super,1);
            j->super->st=chain::VtxArcVtx;
            U.f3(j->super);
          }
          // the window was a VtxArcVtx
          for (typename constraint_entry::IPL::iterator j=ve->ipl[1].begin();
               j!=ve->ipl[1].end();j=tmp) {
            tmp=j; ++tmp;

            if (v->is_right_right()) {
              if (!LC.leftp(v->source_object())) {
                if (U.f1(j->super,ve,ve->target_next_out)) U.f3(j->super);
              }
            } else U.f3(j->super);
          }
          // the window was a VtxArcVtxArcVtx
          for (typename constraint_entry::IPL::iterator j=ve->ipl[2].begin();
               j!=ve->ipl[2].end();j=tmp) {
            tmp=j; ++tmp;
            if (v->is_right_right()) {
              if (!LC.leftp(v->source_object())) {
                if (U.f1(j->super,j->super->iplb[1].v,j->super->iplb[2].v))
                  if (U.f1(j->super,ve,ve->target_next_out))
                    U.f3(j->super);
              }
            } else {
              if (!LC.leftp(j->super->iplb[1].v->vtx->source_object())) {
                if (U.f1(j->super,j->super->iplb[1].v,j->super->iplb[2].v))
                  U.f3(j->super);
              }
            }
          }
        }
      }
      // same with right hand side constraints
      k=de->ch.position;
      LC.set_left_neighbour(&de->ch);
      LC.set_right_neighbour(rc!=YE.end()?*rc:0);
      for (typename VV::reverse_iterator i=de->rhs.constraints.rbegin();
           i!=de->rhs.constraints.rend();++i) {
        if ((*i)->vtx->source_object()==de->d) {
          U.f1(&(*i)->source_chain,0,(*i));
          U.f3(&(*i)->source_chain);
          LC.set_currently_inserted(&(*i)->source_chain);
          k=YE.insert_after(k,&(*i)->source_chain);
          LC.set_currently_inserted(0);
          LC.set_left_neighbour(&(*i)->source_chain);
          (*i)->source_chain.position=k;
          (*i)->source_chain.status=0;
        } else {
          constraint_entry * ve=(*i);
          Vertex_handle v=ve->vtx;
          typename constraint_entry::IPL::iterator tmp;
          for (typename constraint_entry::IPL::iterator j=ve->ipl[0].begin();
               j!=ve->ipl[0].end();j++) {
            U.set_chain_constraint(ve->target_next_out,j->super,1);
            j->super->st=chain::VtxArcVtx;
            U.f3(j->super);
          }
          for (typename constraint_entry::IPL::iterator j=ve->ipl[1].begin();
               j!=ve->ipl[1].end();j=tmp) {
            tmp=j;++tmp;
            if (v->is_left_left()) {
              if (LC.leftp(v->source_object())) {
                if (U.f1(j->super,ve,ve->target_next_out)) U.f3(j->super);
              }
            } else U.f3(j->super);
          }
          for (typename constraint_entry::IPL::iterator j=ve->ipl[2].begin();
               j!=ve->ipl[2].end();j=tmp) {
            tmp=j;++tmp;
            if (v->is_left_left()) {
              if (LC.leftp(v->source_object())) {
                if (U.f1(j->super,j->super->iplb[1].v,j->super->iplb[2].v))
                  if (U.f1(j->super,ve,ve->target_next_out))
                    U.f3(j->super);
              }
            } else {
              if (LC.leftp(j->super->iplb[1].v->vtx->source_object())) {
                if (U.f1(j->super,j->super->iplb[1].v,j->super->iplb[2].v))
                  U.f3(j->super);
              }
            }
          }
        }
      }
      LC.set_left_neighbour(0);
      LC.set_right_neighbour(0);
    } else {
      // leaving an object
      for (typename VV::reverse_iterator i=de->lhs.constraints.rbegin();
           i!=de->lhs.constraints.rend();++i) {
        constraint_entry * ve=(*i);
        Vertex_handle v=ve->vtx;
        if ((*i)->vtx->source_object()==de->d) {

          typename constraint_entry::IPL::iterator tmp;

          // nothing to be done with ipl[0]: the window is already beyond
          // the disk
          
          for (typename constraint_entry::IPL::iterator j=ve->ipl[1].begin();
               j!=ve->ipl[1].end();j=tmp) {
            tmp=j; ++tmp;
            if (dem[v->target_object()].status==0) {
              if (U.f1(j->super,ve,ve->target_next_out))
                U.f3(j->super);
            } else U.f2(j->super,ve);
          }
          for (typename constraint_entry::IPL::iterator j=ve->ipl[2].begin();
               j!=ve->ipl[2].end();j=tmp) {
            tmp=j; ++tmp;
            if (U.f1(j->super,j->super->iplb[1].v,j->super->iplb[2].v)) {
              if (dem[v->target_object()].status==0) {
                if (U.f1(j->super,ve,ve->target_next_out))
                  U.f3(j->super);
              } else U.f2(j->super,ve);
            }
          }
        } else {
          // entering constraints, all their chains have seen or are seeing
          // their windows pass a cusp, and can therefore be removed
          typename constraint_entry::IPL::iterator tmp;
          
          for (typename constraint_entry::IPL::iterator j=ve->ipl[0].begin();
               j!=ve->ipl[0].end();j=tmp) {
            tmp=j; ++tmp;
            U.remove_chain(j->super);
          }
          for (typename constraint_entry::IPL::iterator j=ve->ipl[1].begin();
               j!=ve->ipl[1].end();j=tmp) {
            tmp=j; ++tmp;
            // pass the cusp
            if (U.f1(j->super,ve,ve->target_next_out))
              U.remove_chain(j->super);
          }
          // ve->ipl[2] must be empty
        }
      }
      // same thing with rhs constraints
      for (typename VV::reverse_iterator i=de->rhs.constraints.rbegin();
           i!=de->rhs.constraints.rend();++i) {
        constraint_entry * ve=(*i);
        Vertex_handle v=ve->vtx;
        if ((*i)->vtx->source_object()==de->d) {

          typename constraint_entry::IPL::iterator tmp;
          
          for (typename constraint_entry::IPL::iterator j=ve->ipl[1].begin();
               j!=ve->ipl[1].end();j=tmp) {
            tmp=j; ++tmp;
            if (dem[v->target_object()].status==0) {
              if (U.f1(j->super,ve,ve->target_next_out))
                U.f3(j->super);              
            } else U.f2(j->super,ve);
          }
          for (typename constraint_entry::IPL::iterator j=ve->ipl[2].begin();
               j!=ve->ipl[2].end();j=tmp) {
            tmp=j; ++tmp;
            if (U.f1(j->super,j->super->iplb[1].v,j->super->iplb[2].v)) {
              if (dem[v->target_object()].status==0) {
                if (U.f1(j->super,ve,ve->target_next_out))
                  U.f3(j->super);
              } else U.f2(j->super,ve);
            }
          }
        } else {
          // entering constraints, all their chains have seen or are seeing
          // their windows pass a cusp, and can therefore be removed
          typename constraint_entry::IPL::iterator tmp;
          
          for (typename constraint_entry::IPL::iterator j=ve->ipl[0].begin();
               j!=ve->ipl[0].end();j=tmp) {
            tmp=j; ++tmp;
            U.remove_chain(j->super);
          }
          for (typename constraint_entry::IPL::iterator j=ve->ipl[1].begin();
               j!=ve->ipl[1].end();j=tmp) {
            tmp=j; ++tmp;
            // pass the cusp
            if (U.f1(j->super,ve,ve->target_next_out))
              U.remove_chain(j->super);
          }
          // ve->ipl[2] must be empty
        }
      }
      de->ch.status=1;
      de->status=1;
      set_adjacent_faces(e,de->ch.left,de->ch.right,e->ul());
      e->ul()->set_back_view(de->ch.left->back_view());
      e->ul()->set_front_view(de->ch.right->front_view());
      if (de->ch.position!=YE.begin()) {
        typename Ystructure::iterator lc=de->ch.position;
        --lc;
        (*lc)->right=e->ul();
      }
      typename Ystructure::iterator rc=de->ch.position;
      ++rc;
      if (rc!=YE.end()) (*rc)->left=e->ul();
      YE.erase(de->ch.position);
      if (YE.empty()) {
        upper_outer_face=e->ul();
      }
    }
  }
  infinite_face_=upper_outer_face;
  infinite_face_->set_top_edge(lower_outer_face.top_edge());
  set_adjacent_faces(infinite_face_->top_edge(),0,
                     infinite_face_->top_edge()->ur(),
                     infinite_face_->top_edge()->ul());
  set_adjacent_faces(infinite_face_->bottom_edge(),
                     infinite_face_->bottom_edge()->dl(),
                     infinite_face_->bottom_edge()->dr(),0);

}


#if 0

template < class Gtr_ , class It , class Flip >
template < class DiskIterator , class ConstraintIterator >
void
Antichain<Gtr_,It,Flip>::
compute_graph(DiskIterator first, DiskIterator last,
	      ConstraintIterator  firstc,ConstraintIterator  lastc)
{
  // -------------------------------------------------------------------------
  // The pseudo-triangulation is valid if no pair of objects intersect
  CGAL_expensive_precondition_code( Do_intersect<Self> do_intersect; );
  valid_ = true;
  // -------------------------------------------------------------------------
  // X-structure containing the nodes of the antichain == edges of 
  // Visibility graph
  typedef std::priority_queue<Edge_handle,
    std::vector<Edge_handle>,
    Less_edge_handle<Gtr_> >      Xstructure;
  Xstructure XE;
  // -------------------------------------------------------------------------

  // -------------------------------------------------------------------------
  // Map to keep track of the edge with opposite sign
  std::map<Edge_handle,Edge_handle> opposite;
  // -------------------------------------------------------------------------

  // -------------------------------------------------------------------------
  // Creating the edges and pushing them in the priority queue
  typename Ccw_traits::Set_adjacent_faces_one_to_one set_adjacent_faces;
  // Regular edges - two per disk
  for (DiskIterator it = first; it != last ; ++it ) {
    Edge_handle neg = new Edge(false,&(*it)); 
    Edge_handle pos = new Edge(true,&(*it));
    push_back(*neg); push_back(*pos); 
    XE.push(neg);  XE.push(pos);
    set_adjacent_faces(neg,new Face,new Face,0);
    set_adjacent_faces(pos,new Face,0,0);
    neg->dl()->set_front_view(neg); 
    neg->dr()->set_back_view(pos);
    opposite[(neg)] = pos;
    opposite[(pos)] = neg;
//     neg->set_opposite(pos);
  }
  // Constraint edges - two per constraint
  // We make sure the constraint is upward directed
  typename Gtr_::Is_upward_directed is_upward_directed;
  for (ConstraintIterator c = firstc; c != lastc ; ++c ) {
    CGAL_precondition(c->sup() == 0);
    if (!c->is_constraint())       set_constraint(&(*c));
    if (!c->pi()->is_constraint()) set_constraint(c->pi());
    Vertex_handle d = (is_upward_directed(*c)) ? &(*c) : c->pi();
    push_back(*d->target_cusp_edge());
    push_back(*d->source_cusp_edge());
    XE.push(d->target_cusp_edge());  
    XE.push(d->source_cusp_edge());
  }
  if (XE.size() == 2) { valid_ = false; return; }
  // -------------------------------------------------------------------------
#if 0
  int index = 1;
  for (Edge_iterator e = edges_begin(); e != edges_end(); ++e) 
    {
      if (!e->sign()) {
        Key[long(e->object())] = index; 
        Key[long(&(*e))] = -index;
      }
      else {
        Key[long(&(*e))] = +index; 
        ++index;
      }
    }
#endif
  // -------------------------------------------------------------------------

  // -------------------------------------------------------------------------
  // Y-structure containing at each moment the ordered list of object
  // intersected by the sweeping line.
  // The faces are ordered by their front view
  // An object is represented by a pair of faces (f1,f2) such that
  // f1->front_object() == f2->back_object()
  typedef std::set<Face_handle,Less_face_handle<Self> >  Ystructure;
  Ystructure YE;
  // -------------------------------------------------------------------------

  // -------------------------------------------------------------------------
  // Splice function object used when inserting a constraint
  typename Ccw_traits::Splice splice_ccw;
  typename Cw_traits::Splice  splice_cw;
  // -------------------------------------------------------------------------
  // The infinite face helps us to identify the convex-hull vertices
  infinite_face_ = new Face;
  infinite_face_->set_bottom_edge(XE.top());
  // -------------------------------------------------------------------------

  // -------------------------------------------------------------------------
  // Start of the Bentley-Ottman sweep
  Face_handle infinite = NULL;
  while ( !XE.empty() && valid_ ) {
    Edge_handle pe = XE.top(); XE.pop();
    // ---------------------------------------------------------------------
    // Treating a node with in_degree 2 and out_degree 1
    if (!pe->sign()) {
      // -----------------------------------------------------------------
      // Inserting a Regular Edge
      if (pe->object() != 0) { // regular Edge
        // -------------------------------------------------------------
        // Insert two new faces in YE
        Face_handle upleft = new Face; upleft->set_front_view(pe);
        typename Ystructure::iterator upleft_it = YE.upper_bound(upleft);
        delete upleft;
        if ( upleft_it == YE.end() ) upleft = infinite;
        else { upleft = *upleft_it; YE.erase(upleft_it); }
        set_adjacent_faces(pe,pe->dl(),pe->dr(),upleft);
        if (upleft != 0) {
          // ---------------------------------------------------------
          CGAL_expensive_precondition_code(
                                           valid_ = !do_intersect(pe->ul()->back_view(),pe);
                                           if (valid_ == false) return;
                                           valid_ = !do_intersect(pe->ul()->front_view(),pe);
                                           if (valid_ == false) return;
                                           );
          // ---------------------------------------------------------
          pe->dr()->set_front_view(pe->ul()->front_view());
          pe->dl()->set_back_view (pe->ul()->back_view());
        }
        YE.insert(pe->dl()); YE.insert(pe->dr());
        // -------------------------------------------------------------
      }
      // -----------------------------------------------------------------
      // Inserting a cusp Edge from an xx-left constraint
      else if (pe->sup()->is_xx_left()) { 
        // -------------------------------------------------------------
        // Find the place where to insert the constraint
        Face_handle right = new Face; right->set_front_view(pe);
        typename Ystructure::iterator right_it = YE.upper_bound(right);
        CGAL_precondition( right_it != YE.end() );
        delete right; right = *right_it;
        // -------------------------------------------------------------
        CGAL_expensive_precondition_code(
                                         valid_ = !do_intersect(right->front_view(),pe);
                                         if (valid_ == false) return;
                                         );
        // -------------------------------------------------------------
        // Splitting the boundary of the target object of pe->sup()
        CGAL_precondition(right->back_view()->object() != 0);
        splice_ccw(right->back_view(),pe->sup());
        // -------------------------------------------------------------
        pe->dl()->set_back_view (pe->sup()->cw_target_edge());
        pe->dl()->set_front_view(pe);
        right->set_back_view(pe);
        YE.insert(pe->dl());
        // -------------------------------------------------------------
      }
      // -----------------------------------------------------------------
      // Inserting a cusp Edge from an xx-right constraint
      else {
        // -------------------------------------------------------------
        // Find the place where to insert the constraint
        Face_handle left = new Face; left->set_front_view(pe);
        typename Ystructure::iterator left_it = YE.upper_bound(left);
        CGAL_precondition( left_it != YE.end() );
        delete left; left = *left_it;
        YE.erase(left_it);
        // -------------------------------------------------------------
        CGAL_expensive_precondition_code(
                                         valid_ = !do_intersect(left->back_view(),pe);
                                         if (valid_ == false) return;
                                         );
        // -------------------------------------------------------------
        CGAL_precondition(left->front_view()->object() != 0);
        splice_ccw(left->front_view(),pe->sup());
        // -------------------------------------------------------------
        pe->dr()->set_back_view (pe);
        pe->dr()->set_front_view(pe->sup()->ccw_target_edge());
        left->set_front_view(pe);
        YE.insert(left); YE.insert(pe->dr());
        // -------------------------------------------------------------
      }
      // -----------------------------------------------------------------
    }
    // ---------------------------------------------------------------------
    // Treating a node with in_degree 1 and out_degree 2
    else {
      // -----------------------------------------------------------------
      // Inserting a regular Edge
      if (pe->object() != 0) {
        // -------------------------------------------------------------
        // Find the two adjacent faces in YE that connect to pe and
        // erase them from YE.
        Face_handle tmpf = new Face; tmpf->set_front_view(pe);
        typename Ystructure::iterator it1 = YE.find(tmpf);
        typename Ystructure::iterator it2 = it1; ++it2;
        Face_handle upleft = *it1; Face_handle upright = *it2;
        delete tmpf; YE.erase(it1); YE.erase(it2);
        // -------------------------------------------------------------
        // Update existing faces
        if ( XE.empty() ) {
          delete pe->dl();
          set_adjacent_faces(pe,0,upright,upleft);
        }
        else { 
          set_adjacent_faces(pe,pe->dl(),upright,upleft);
          pe->dl()->set_front_view(upright->front_view());
          pe->dl()->set_back_view (upleft->back_view());
          // ---------------------------------------------------------
          CGAL_expensive_precondition_code(
                                           valid_ = !do_intersect(pe->dl()->back_view(),
                                                                  pe->dl()->front_view());
                                           if (valid_ == false) return;
                                           );
          // ---------------------------------------------------------
          YE.insert(pe->dl());
        }
        // -------------------------------------------------------------
        infinite = pe->dl();
        // -------------------------------------------------------------
      }
      // -----------------------------------------------------------------
      // Inserting a tail cusp Edge from an left-xx constraint
      else if (pe->sup()->is_left_xx()) {
        // -------------------------------------------------------------
        // Find the face in YE with pe->sup() as a front view.
        Face_handle upleft = new Face; 
        upleft->set_front_view(pe->sup()->target_cusp_edge());
        typename Ystructure::iterator upleft_it = YE.find(upleft);
        typename Ystructure::iterator right = upleft_it; ++right;
        delete upleft; upleft = *upleft_it;
        YE.erase(upleft_it);
        // -------------------------------------------------------------
        // Update degenerate face emanating from edge.
        upleft->set_sup(pe->sup());
        upleft->set_inf(pe->sup()->pi());
        CGAL_precondition(upleft->back_view()->object() != 0);
        splice_ccw(upleft->back_view(),pe->sup());
        // -------------------------------------------------------------
        /* delete pe->ul(); */ set_adjacent_faces(pe,0,0,upleft);
        CGAL_precondition(pe->dl() == 0 && pe->ur() == 0);
        // -------------------------------------------------------------
        // Update the back view of right.
        (*right)->set_back_view(pe->sup()->cw_source_edge());
        // -------------------------------------------------------------
      }
      // -----------------------------------------------------------------
      // Inserting a tail cusp Edge from a right-xx constraint
      else {
        // -------------------------------------------------------------
        // Find the face in YE with pe->sup() as a front view.
        Face_handle left = new Face; 
        left->set_front_view(pe->sup()->target_cusp_edge());
        typename Ystructure::iterator left_it = YE.find(left);
        typename Ystructure::iterator upright = left_it; ++upright;
        delete left; left = *left_it;
        // -------------------------------------------------------------
        // Set upright as the degenerate face of the edge pe.
        (*upright)->set_sup(pe->sup());
        (*upright)->set_inf(pe->sup()->pi());
        CGAL_precondition((*upright)->front_view()->object() != 0);
        splice_ccw((*upright)->front_view(),pe->sup());
        // -------------------------------------------------------------
        /* delete pe->ur(); */ set_adjacent_faces(pe,0,*upright,0);
        CGAL_precondition(pe->dl() == 0 && pe->ul() == 0);
        // -------------------------------------------------------------
        // Update the back view of left and insert the face.
        YE.erase(left_it); YE.erase(upright); 
        left->set_front_view(pe->sup()->ccw_source_edge());
        YE.insert(left); 
        // -------------------------------------------------------------
      }
      // -----------------------------------------------------------------
    }
    if (XE.empty()) infinite_face_->set_top_edge(pe);
  } // end while ( !XE.empty() )
    // -------------------------------------------------------------------------

    // -------------------------------------------------------------------------
    // Splicing the boundary of objects by inserting the pi of the constraints
    // The fact that the edges come in pairs allows us to recover the edge with
    // opposite sign on the same disk.
  for (Edge_iterator e = edges_begin(); e != edges_end() ; ++e) {
    if (e->object() != 0 && e->sup() != 0) {
      Edge_handle eo = opposite[(&(*e))];
//       CGAL_precondition(eo == e->opposite());

      Vertex_handle v = e->sup();
      splice_cw(eo,v->pi());
      while (v->ccw_edge(e->object())->sup() != 0) {
        splice_cw(v->pi()->ccw_edge(e->object()),
                  v->ccw_edge(e->object())->sup()->pi());
        v = v->ccw_edge(e->object())->sup();
      }
    }
  }
  // -------------------------------------------------------------------------
}

#endif

template < class Gtr_ , class It , class Flip >
template < class Tr>
void
Antichain<Gtr_,It,Flip>::compute_vertices(Tr tr)
{
  // -------------------------------------------------------------------------
  // All the operators from this methods are taken from the algorithm traits 
  // class Tr. Two traits classes are given : Ccw_traits and Cw_traits to
  // compute G and G_* respectively
  typename Tr::Sup     sup; typename Tr::Inf     inf;
  typename Tr::Set_sup set_sup;// typename Tr::Set_inf set_inf;
  typename Tr::Vertex_creator vc; 
  typename Tr::Set_adjacent_faces set_adjacent_old_faces;
  typename Tr::Set_adjacent_faces_one_to_one  set_adjacent_faces;
  typename Tr::Dr      dr;  typename Tr::Dl      dl;
  typename Tr::Ur      ur;  typename Tr::Ul      ul;
  typename Tr::Cw_target_edge cw_target_edge;
  typename Tr::Cw_source_edge cw_source_edge;
  typename Tr::Ccw_target_edge ccw_target_edge;
  typename Tr::Ccw_source_edge ccw_source_edge;
  typename Tr::Splice splice;
  typename Tr::Ccw_edge ccw_edge;
  typedef Tr TR;
  // -------------------------------------------------------------------------

  // -------------------------------------------------------------------------
  // Saving the Edge --> Face pointers because we will lose them during the
  // rotational sweep below
  Edge aaa[this->size()];
  {
    size_t i=0;
    for (Edge_iterator e = edges_begin(); e != edges_end() ; ++e,++i) {
      new (aaa+i) Edge(e->sign(),0);
      if (e->sign()) 
        set_adjacent_old_faces(aaa+i,dl(&(*e)),ur(&(*e)),ul(&(*e)));
      else set_adjacent_old_faces(aaa+i,dl(&(*e)),dr(&(*e)),ul(&(*e)));
    }
  }

  // -------------------------------------------------------------------------

  // -------------------------------------------------------------------------
  // death is the event queue. 
  // It contains the vertices which are about to appear.
  // Pushing first candidates to initialize death vertices
  std::set<Vertex_handle,typename Tr::Less_vertex_handle> death;

  for (Face_iterator fit = faces_begin() ; fit != faces_end() ; ++fit) {
    Face_handle f = &(*fit);
    if (is_candidate(f, Tr())) { // af: gcc needs it, for bcc it must not be there
      if (sup(f) == 0) {
        set_sup(f,vc(f->bottom_edge(),f->top_edge()));
      }
      death.insert(sup(f));
    }
  }
  // -------------------------------------------------------------------------

  // -------------------------------------------------------------------------
  // Rotational sweep a la Bentley-Otman
  while ( death.size() != 0 ) {
    Vertex_handle vmin = *death.begin();
    death.erase(death.begin());
    // ---------------------------------------------------------------------
    Edge_handle top = (vmin->is_constraint()) ?  cw_target_edge(vmin) : 
      inf(vmin)->top_edge();
    Edge_handle bot = (vmin->is_constraint()) ?  cw_source_edge(vmin) : 
      inf(vmin)->bottom_edge();
    // ---------------------------------------------------------------------
    // A new Vertex ( == bitangent) has been found 
    // splice the two edges defining it
    if (vmin != sup(top)) splice(top,vmin);
    if (vmin != sup(bot)) splice(bot,vmin);
    // ---------------------------------------------------------------------
    // The two possible adjacent candidates must be erased
    Face_handle ef[6];
    if (bot->sign()) { ef[0] = dl(bot); ef[1] = ur(bot); ef[2] = ul(bot); }
    else             { ef[0] = dl(bot); ef[1] = dr(bot); ef[2] = ul(bot); }
    if (top->sign()) { ef[3] = dl(top); ef[4] = ur(top); ef[5] = ul(top); }
    else             { ef[3] = dl(top); ef[4] = dr(top); ef[5] = ul(top); }
    for (int i = 0 ; i < 6 ; i++) 
      if (is_candidate(ef[i], Tr()) && sup(ef[i]) != vmin) {
        death.erase(sup(ef[i]));
        if (!sup(ef[i])->is_constraint()) {
          delete sup(ef[i]);
          set_sup(ef[i],0);
        }
      }
    // ---------------------------------------------------------------------
    // Update the antichain
    CGAL_precondition(sup(vmin) == 0);
    sweep_constraint(vmin, Tr());
    set_adjacent_old_faces(top,0,0,0);
    set_adjacent_old_faces(bot,0,0,0);
    // ---------------------------------------------------------------------
    // Computing new candidates by looking at the adjacent faces of
    // bot and top.
    bot = ccw_source_edge(vmin);
    top = ccw_target_edge(vmin);
    if (sup(bot) != 0 && sup(bot)->is_constraint() // && 
//         is_minimal(sup(bot),tr)
        ) 
      death.insert(sup(bot));
    if (sup(top) != 0 && sup(top)->is_constraint() // && 
//         is_minimal(sup(top),tr)
        )
      death.insert(sup(top));
    Face_handle af[] = { dl(bot) , ur(bot) , dl(top) , ur(top) };
    for (int i = 0 ; i < 4 ; i++) 
      if (is_candidate(af[i], Tr()) &&
          sup(af[i]) != vmin && sup(af[i]) != vmin->pi()) {
        if (sup(af[i]) == 0) set_sup(af[i],vc(af[i]->bottom_edge(),
                                              af[i]->top_edge()));
        death.insert(sup(af[i]));
      }
    // ---------------------------------------------------------------------
  } // end while (death.size() != 0) 

    // Recovering the initial Edge --> Face pointers with a
  {
    Edge_handle ea = aaa;
    for (Edge_iterator e = edges_begin(); e != edges_end() ; ++e,++ea) {
      if (e->object()) {
        if (e->sign()) 
          set_adjacent_faces(&(*e),dl(ea),ur(ea),ul(ea));
        else set_adjacent_faces(&(*e),dl(ea),dr(ea),ul(ea));
        for (Edge_handle ee=e.operator->(); sup(ee);) {
          ee=ccw_edge(sup(ee),ee->object());
          set_adjacent_faces(ee,0,0,0);     
        }
      }
    }
  }
  // -------------------------------------------------------------------------
}


template < class Gtr_ , class It , class Flip >
template < class Tr>
void
Antichain<Gtr_,It,Flip>::compute_minimals(Tr tr)
{
  // -------------------------------------------------------------------------
  // Operators used in this method.
  typename Tr::Sup sup; 
  typedef typename Tr::Right_traits Right_traits;
  typename Right_traits::Target_cusp_edge target_cusp_edge;
  typename Right_traits::Top_edge         top_edge;
  // -------------------------------------------------------------------------
  // Computing the set of minimal bitangents
  Vertex_handle ll = 0, rr = 0;
  for (Face_iterator fi = faces_begin(); fi != faces_end() ; ++fi) {
    Face_handle f = &(*fi);
    if (sup(f) != 0 && is_minimal(sup(f),tr)) {
      if (!sup(f)->is_constraint() || 
          top_edge(f) == target_cusp_edge(sup(f)))  {
        Vertex_handle v = sup(f);
        if (v != ll && v != rr) push_back_minimal(sup(f),tr);
        if (ll == 0 && v->is_left_left() && 
            is_on_convex_hull(v))      ll = v;
        else if (rr == 0 && v->is_right_right() && 
                 is_on_convex_hull(v)) rr = v;
      }
    }
  }
}


// Find topmost and bottommost edge.
template < class Gtr_ , class It , class Flip >
template < class Tr >
void
Antichain<Gtr_,It,Flip>::
initialize_convex_hull(Tr /*tr*/) const
{
  typename Tr::Set_adjacent_faces            set_adjacent_faces;
  typename Tr::Dr   dr;  typename Tr::Dl  dl;
  typename Tr::Ur   ur;  typename Tr::Ul  ul;
  typename Tr::Ccw_source_edge ccw_source_edge;
  typename Tr::Ccw_target_edge ccw_target_edge;
  typename Tr::Cw_source_edge cw_source_edge;
  typename Tr::Cw_target_edge cw_target_edge;
  typename Tr::Sup sup;
  typename Tr::Inf inf;
 
  Edge_handle e = infinite_face()->bottom_edge();
  do {
    e = ccw_source_edge(sup(e));
    set_adjacent_faces(e,dl(e),dr(e),infinite_face());
  } while (sup(e)&&e != infinite_face()->bottom_edge()); 

  if (!sup(e)) {
    e=infinite_face()->bottom_edge();
    do {
      e = cw_target_edge(inf(e));
      set_adjacent_faces(e,dl(e),dr(e),infinite_face());
    } while (inf(e));
  }

  e = infinite_face()->top_edge();
  do {
    e = ccw_target_edge(sup(e));
    set_adjacent_faces(e,infinite_face(), ur(e),ul(e));
  } while (sup(e)&&e != infinite_face()->top_edge()); 

  if (!sup(e)) {
    e=infinite_face()->top_edge();
    do {
      e = cw_source_edge(inf(e));
      set_adjacent_faces(e,dl(e),dr(e),infinite_face());
    } while (inf(e));
  }
}

// The antichain has 3n-1 faces whereas the pseudo-triangulation has only 3n-3
// bitangents. There are two faces which do not have their sup pointer set at
// this point. These are ul(bottommost) and dr(topmost). This function sets the
// sink of the two supplementary faces.
template < class Gtr_ , class It , class Flip >
template < class Tr >
void
Antichain<Gtr_,It,Flip>::fix_extreme_edges(Tr /*tr*/) const // warning tr is never used 
{
  //--------------------------------------------------------------------------
  // Operators used in this method.
  typename Tr::Sup sup; typename Tr::Inf inf;
  typename Tr::Set_sup set_sup; typename Tr::Set_inf set_inf; 
  typename Tr::Dr   dr;  typename Tr::Dl  dl;
  typename Tr::Ur   ur;  typename Tr::Ul  ul;
  typename Tr::CcL  ccL; typename Tr::CcR ccR;
  typename Tr::Set_adjacent_faces set_adjacent_faces;
//   typename Tr::Ccw_edge ccw_edge;
  typename Tr::Ccw_target_edge ccw_target_edge;
//   typename Tr::Splice splice;
  typename Tr::Source_object source_object;
//   typename Tr::Target_object target_object;
  typename Tr::Ccw_source_edge ccw_source_edge;
  typename Tr::Ccw_edge ccw_edge;
  typename Tr::Sign sign;
  // ------------------------------------------------------------------------- 
  Edge_handle bot = infinite_face()->top_edge();    // Bottommost edge
  Edge_handle top = infinite_face()->bottom_edge(); // Topmost edge
  set_adjacent_faces(bot,infinite_face(),ur(bot),ul(bot));
  set_adjacent_faces(top,dl(top),dr(top),infinite_face());
  // ------------------------------------------------------------------------- 
  // A face may appear twice in the antichain due to our identification 
  // v = v->pi()->pi(). If this is the case, delete one of the two and update
  // pointers. 
  if (ur(bot) == dr(top) && ul(bot) != dl(top)) {
    Face_handle f = ul(bot);
    Edge_handle t = f->top_edge();
    set_inf(dl(top),inf(f));
    set_sup(inf(f),dl(top));
    while (t->object()!=bot->object()) {
      if (sign(t)) {
        set_adjacent_faces(t,dl(top),ur(t),ul(t));
      } else {
        set_adjacent_faces(t,dl(top),dr(t),ul(t));        
      }
      Vertex_handle v=sup(t);
      if (source_object(v)==t->object()) t=ccw_target_edge(v);
      else t=ccw_source_edge(v);
    }
    set_adjacent_faces(t,dl(top),ur(t),ul(t));
    t=ccw_edge(sup(t),t->object());
    set_adjacent_faces(t,dl(top),ur(t),ul(t));
    set_adjacent_faces(bot,dl(bot),ur(bot),dl(top));
    delete f;
  }
  else if (ul(bot) == dl(top) && ur(bot) != dr(top)) {
    Face_handle f = dr(top);
    Edge_handle t = f->bottom_edge();
    set_inf(ur(bot),inf(f));
    set_sup(inf(f),ur(bot));
    while (t->object()!=top->object()) {
      if (sign(t)) {
        set_adjacent_faces(t,dl(t),ur(bot),ul(t));
      } else {
        set_adjacent_faces(t,dl(t),dr(t),ur(bot));        
      }
      Vertex_handle v=sup(t);
      if (source_object(v)==t->object()) t=ccw_target_edge(v);
      else t=ccw_source_edge(v);
    }
    set_adjacent_faces(t,dl(t),ur(bot),ul(t));
    t=ccw_edge(sup(t),t->object());
    set_adjacent_faces(t,dl(t),ur(bot),ul(t));
    set_adjacent_faces(top,dl(top),ur(bot),ul(top));
    delete f;
  }
  // ------------------------------------------------------------------------- 
  // Set the sink of the two faces on the convex-hull.
  if (sup(ur(bot)) != 0 && sup(ul(bot)) == 0) {
    Vertex_handle v = sup(ur(bot));
    Vertex_handle w = inf(dr(top));
    if (*v == *w->pi()) set_sup(ul(bot),w);
    else {
      do { v = ccR(v); } while (v->is_constraint());
      CGAL_precondition(v != 0 && inf(v) != 0);
      set_sup(ul(bot),inf(inf(v)));
    }
    sup(ur(bot))->set_pi(sup(ul(bot)));
  }
  if (sup(dl(top)) != 0 && sup(dr(top)) == 0) {
    Vertex_handle v = sup(dl(top));
    Vertex_handle w = inf(ul(bot));
    if (*v == *w->pi()) set_sup(dr(top),w);
    else {
      do { v = ccL(v); } while (v->is_constraint());
      CGAL_precondition(v != 0 && inf(v) != 0);
      set_sup(dr(top),inf(inf(v)));
    }
    sup(dl(top))->set_pi(sup(dr(top)));
  }
}

template < class Gtr_ , class It , class Flip >
void
Antichain<Gtr_,It,Flip>::glue()
{
  Less_bitangent<Gtr_> lb;
  typename Ccw_traits::Splice splice_ccw;
  typename Cw_traits::Splice splice_cw;
  typename Ccw_traits::Set_adjacent_faces saj;

  for (Edge_iterator e = edges_begin(); e != edges_end() ; ++e) {
    if (e->object() != 0) {
      Edge_iterator eo = e; if (e->sign()) --eo; else ++eo;
      Edge_handle e1=&*e;
      Edge_handle e2=&*eo;
      Edge_handle e20=e2;
      bool b=false;
      Vertex_handle v1=e1->sup();
      bool chiv1e1,chiv1e2;
      while (e2->inf()) {
        e2=e2->inf()->cw_edge(e2->object());
        b=true;
      }
      while (e1->sup()&&e2->sup()) {
        if (e1->sup()==e2->sup()->pi()) {
          e1=e1->sup()->ccw_edge(e1->object());
          e2=e2->sup()->ccw_edge(e2->object());
          if (e1->sup()) chiv1e1=lb(*v1,*e1->sup());
          if (e2->sup()) chiv1e2=lb(*v1,*e2->sup()->pi());
          continue;
        }
        if (*e1->sup()==*e2->sup()->pi()) {
          e1->sup()->set_pi(e2->sup());
          e1=e1->sup()->ccw_edge(e1->object());
          e2=e2->sup()->ccw_edge(e2->object());
          if (e1->sup()) chiv1e1=lb(*v1,*e1->sup());
          if (e2->sup()) chiv1e2=lb(*v1,*e2->sup()->pi());
          continue;
        }
        if (e1->sup()==v1) {
          if (v1->ccw_edge(e1->object())->sup()) chiv1e1=lb(*v1,*v1->ccw_edge(e1->object())->sup());
          chiv1e2=lb(*v1,*e2->sup()->pi());
          splice_cw(e2,e1->sup()->pi());
          saj(e2->inf()->cw_edge(e2->object()),e2->dl(),
              e2->sign()?e2->ur():e2->dr(),e2->ul());
          e1=e1->sup()->ccw_edge(e1->object());
          continue;
        }
        if (chiv1e1? (chiv1e2?lb(*e1->sup(),*e2->sup()->pi()):true)
            : (chiv1e2?false:lb(*e1->sup(),*e2->sup()->pi()))
            ) {
          splice_cw(e2,e1->sup()->pi());
          saj(e2->inf()->cw_edge(e2->object()),e2->dl(),
              e2->sign()?e2->ur():e2->dr(),e2->ul());
          e1=e1->sup()->ccw_edge(e1->object());
          if (e1->sup()) chiv1e1=lb(*v1,*e1->sup());
        } else {
          splice_ccw(e1,e2->sup()->pi());
          saj(e1->sup()->ccw_edge(e1->object()),e1->dl(),
              e1->sign()?e1->ur():e1->dr(),e1->ul());
          e1=e1->sup()->ccw_edge(e1->object());
          e2=e2->sup()->ccw_edge(e2->object());
          if (e2->sup()) chiv1e2=lb(*v1,*e2->sup()->pi());
        }
      }
      for (;e1->sup();e1=e1->sup()->ccw_edge(e1->object()))
        splice_cw(e2,e1->sup()->pi());
      for (;e2!=e20&&e2->sup();e2=e2->sup()->ccw_edge(e2->object())) {
        splice_ccw(e1,e2->sup()->pi());
        e1=e1->sup()->ccw_edge(e1->object());
      }
    }
  }
  for (Edge_iterator e = edges_begin(); e != edges_end() ; ++e) {
    if (e->object()) {
      if (e->sup()) {
        e->sup()->pi()->set_final_antichain(true);
      }
      Edge_handle e1,e2;
      for (e1=&*e;e1->sup();) e1=e1->sup()->ccw_edge(e1->object());
      for (e2=&*e;e2->inf();) e2=e2->inf()->cw_edge(e2->object());
      if (e1!=&*e) {
        e2->set_inf(e1->inf());
        e2->inf()->set_ccw_edge(e2,e2->object());
        saj(e2,e1->dl()?e1->dl():e2->dl(),
            e2->sign()?(e1->ur()?e1->ur():e2->ur()):
            (e1->dr()?e1->dr():e2->dr()),
            e1->ul()?e1->ul():e2->ul());
        delete e1;
      } else {
        e1->set_sup(e2->sup());
        e1->sup()->set_cw_edge(e1,e1->object());
        saj(e1,e1->dl()?e1->dl():e2->dl(),
            e2->sign()?(e1->ur()?e1->ur():e2->ur()):
            (e1->dr()?e1->dr():e2->dr()),
            e1->ul()?e1->ul():e2->ul());
        delete e2;
      }
    }
  }
  
  // ------------------------------------------------------------------------- 
}
    
// Identifying vertices of the convex-hull

template < class Gtr_ , class It , class Flip >
inline bool 
Antichain<Gtr_,It,Flip>::
is_on_convex_hull(Vertex_handle v) const
{
  return ((v->is_left_left()   &&  v->cw_source_edge() != 0 && 
           v->cw_source_edge()->dr() == infinite_face()) ||
          (v->is_right_right() &&  v->cw_target_edge() != 0 &&
           v->cw_target_edge()->ul() == infinite_face()));
}


// Minimality Testing

template < class Gtr_ , class It , class Flip >
template < class Tr > 
inline bool 
Antichain<Gtr_,It,Flip>::
is_minimal(const Vertex_handle& v, Tr /*tr*/)  const
{
  typename  Tr::Left_traits lt;
  typename Tr::Right_traits rt;
  return (is_xx_minimal(v,lt) && 
          is_xx_minimal(v,rt));
}

template < class Gtr_ , class It , class Flip >
template < class Tr > 
inline bool 
Antichain<Gtr_,It,Flip>::
is_xx_minimal(const Vertex_handle& v, Tr /*tr*/ ) const
{ 
  typename Tr::Cw_source_edge cw_source_edge;
  if (v == 0 || cw_source_edge(v) == 0) return false;

  return cw_source_edge(v)->is_in_antichain();
}

// Methods to manage the list of minimals

template < class Gtr_ , class It , class Flip >
template < class Tr > 
void 
Antichain<Gtr_,It,Flip>::clear_minimals(Tr tr) {
  Minimals_iterator first = minimals_begin(tr);
  Minimals_iterator last  = minimals_end(tr);
  while (first != last) erase_minimal(first++,tr);
}

template < class Gtr_ , class It , class Flip >
void          
Antichain<Gtr_,It,Flip>::
push_back_minimal(Vertex_handle v,Ccw_traits)
{
  Less_bitangent<Gtr_> lv;
  Minimals_iterator w = (is_straight()) ? std::lower_bound(minimals_begin(),
                                                           minimals_end(), 
                                                           *v,lv) :
    minimals_end();
  minimals_ccw_.insert(w,*v);
}

template < class Gtr_ , class It , class Flip >
void          
Antichain<Gtr_,It,Flip>::
push_back_minimal(Vertex_handle v,Cw_traits)
{
  Greater_bitangent<Gtr_> lv;
  Minimals_iterator w = (is_straight()) ? std::lower_bound(cw_minimals_begin(),
                                                           cw_minimals_end(), 
                                                           *v,lv) :
    cw_minimals_end();
  minimals_cw_.insert(w,*v);
}

template < class Gtr_ , class It , class Flip >
typename Antichain<Gtr_,It,Flip>::Vertex_handle
Antichain<Gtr_,It,Flip>::pop_minimal(bool finite/* = false*/)
{
  if (minimals_begin() == minimals_end()) return 0;
  if (finite == false) {
    Vertex_handle min = &(*minimals_begin());
    return min;
  }
  else {
    Vertex_handle min = pop_minimal(false);
    if (min != 0 && is_on_convex_hull(min)) return min;
    while (minimals_begin() != minimals_end() && 
           min->sup() != 0 && !is_on_convex_hull(min)) {
      erase_minimal(min);
      min = pop_minimal(false);
    }
    if (min != 0 && is_on_convex_hull(min)) return min;
    return (min == 0 || min->sup() != 0) ? 0 : min;
  }
}

template < class Gtr_ , class It , class Flip >
template < class Tr >
void
Antichain<Gtr_,It,Flip>::sweep(Vertex_handle v , Tr tr)
{
  //--------------------------------------------------------------------------
  // Operators used in this method.
  typename Tr::Inf inf;//  typename Tr::Sup sup; 
  typename Tr::Cw_target_edge  cw_target_edge;
  typename Tr::Cw_source_edge  cw_source_edge;
  typename Tr::Ccw_target_edge ccw_target_edge;
  typename Tr::Ccw_source_edge ccw_source_edge;
  typename Tr::CcR ccR; typename Tr::CcL ccL;
  typename Tr::CwR cwR; typename Tr::CwL cwL;
  typename Tr::Merge  merge(this);
  typename Tr::Ul ul;
  typename Tr::Dl dl;
  typename Tr::Ur ur;
  typename Tr::Dr dr;


  CGAL_precondition(is_minimal(v,tr));
  CGAL_precondition(cwL(v) != 0 && cwR(v) != 0);

  bool dont_delete_inf_v=
    (inf(v)==ul(infinite_face()->top_edge())&&
      inf(v)==dr(infinite_face()->bottom_edge())) ||
    (inf(v)==ur(infinite_face()->top_edge())&&
     inf(v)==dl(infinite_face()->bottom_edge()));
  // -------------------------------------------------------------------------
  // Removing cwR(v) and cwL(v) from the minimal list corresponding to the
  // opposite orientation. Adding v.
  typename Tr::Cw_traits cw_traits;
  erase_minimal(v,tr); push_back_minimal(v,cw_traits);
  CGAL_assertion(v!=cwR(v)&&v!=cwL(v));  
  if (is_minimal(cwL(v),cw_traits)) erase_minimal(cwL(v),cw_traits);
  if (cwR(v) != cwL(v) && is_minimal(cwR(v),cw_traits))
    erase_minimal(cwR(v),cw_traits);
  //--------------------------------------------------------------------------
  // Flipping v and updating the antichain
  if (v->is_constraint()) {
    sweep_constraint(v, Tr());
    if (!linear_space()) sweep_constraint(v->pi(),Tr(),true);
  }
  else {
    compute_phi(v,tr);       // Flip bitangent in pseudo-triangulation
    sweep_regular(v, Tr());   // Modify the antichain
    if (!linear_space()) sweep_regular(v->pi(),Tr(),true);
  }
  erase(cw_source_edge(v)); 
  erase(cw_target_edge(v));
  //--------------------------------------------------------------------------
  // Pushing the two new created edges 
  push_back(*ccw_source_edge(v)); 
  push_back(*ccw_target_edge(v));
  //--------------------------------------------------------------------------
  // Updating minimals
  if (ccL(v)&&is_minimal(ccL(v),tr))     push_back_minimal(ccL(v),tr);
  if (ccR(v)&&ccR(v) != ccL(v) && is_minimal(ccR(v),tr))
    push_back_minimal(ccR(v),tr);

  // Deleting (inf(v)), along with its vertices that are no longer adjacent
  // to any face.
  if (linear_space()&&!dont_delete_inf_v) {
    Face_handle faces[3]={0,0,0};
    if (v->is_constraint()) {
      faces[0]=v->source_cusp_face();
      faces[1]=v->target_cusp_face();
    } else {
      faces[0]=v->inf();
    }
    for (int fi=0;faces[fi];++fi) {
      std::vector<Vertex_handle> vv;
      Edge_handle prev=0;
      for (typename Face::Border_iterator i=--(faces[fi]->top_end());
           i!=faces[fi]->top_end()&&i.operator->()!=prev;--i) {
        vv.push_back(inf(i.operator->()));
        prev=i.operator->();
      }
      prev=0;
      for (typename Face::Border_iterator i=--(faces[fi]->bottom_end());
           i!=faces[fi]->top_end()&&i.operator->()!=prev;--i) {
        if (inf(i.operator->())!=inf(faces[fi])) 
          vv.push_back(inf(i.operator->()));
        prev=i.operator->();
      }
      delete faces[fi];

      for (typename std::vector<Vertex_handle>::iterator i=vv.begin();
           i!=vv.end();++i) {
        if (is_on_convex_hull(*i)||(*i)->is_constraint()||
            (*i)->final_antichain()) 
          continue;
        Face_handle adj[4][3];
        Edge_handle e[4]={cw_source_edge(*i),cw_target_edge(*i),
                          ccw_source_edge(*i),ccw_target_edge(*i)};
        for (int j=0;j<4;j++) {
          adj[j][0]=e[j]->dl();
          adj[j][1]=e[j]->sign()?e[j]->ur():e[j]->dr();
          adj[j][2]=e[j]->ul();
        }
        bool skip=false;
        for (int k=0;k==0;k=1) {
          for (int l=0;l<3;l++) skip=skip||adj[k][l]!=adj[k+2][l];
        }
        if (!skip) { 
          merge(cw_source_edge(*i),ccw_source_edge(*i));
          merge(cw_target_edge(*i),ccw_target_edge(*i));        
          delete *i;
        }
      }
    }
  }
}

template < class Gtr_ , class It , class Flip >
template < class Tr >
void
Antichain<Gtr_,It,Flip>::sweep_all_minimals(Tr tr)
{
  std::list<Vertex_handle> mins;
  for (Minimals_iterator m = minimals_begin(tr); m != minimals_end(tr); ++m)
    mins.push_back(&(*m));
  typename std::list<Vertex_handle>::iterator clic = mins.begin();
  for ( ; clic != mins.end() ; ++clic) 
    if (is_minimal(*clic,tr)) sweep(*clic,tr);
}

template < class Gtr_ , class It , class Flip >
template < class Tr >
void
Antichain<Gtr_,It,Flip>::sweep_good(Vertex_handle v,Tr tr)
{
  typename Tr::CcR ccR; typename Tr::CcL ccL; typename Tr::Sup sup;
  //    typename Tr::CwR cwR; typename Tr::CwL cwL; 
  // -------------------------------------------------------------------------
  sweep(v,tr);
  // -------------------------------------------------------------------------
  // If ccR(v) and/or ccL(v) are geometrically equal to v and if ccR(v) and/or
  // ccL(v) are minimal, we sweep them.
  CGAL_precondition(ccR(v) != 0 && ccL(v) != 0 && sup(sup(v)) != 0);
  Vertex_handle w[] = { ccR(v) , ccL(v) , sup(sup(v)) };
  typename Gtr_::Equal_as_segments equal_as_segments;
  for (int i = 0; i < 3 ; i++)
    if (equal_as_segments(*v,*w[i]) && is_minimal(w[i],tr))
      sweep(w[i],tr);
  // -------------------------------------------------------------------------
}
template < class Gtr_ , class It , class Flip >
template < class Tr >
void
Antichain<Gtr_,It,Flip>::sweep_good_all_minimals(Tr tr)
{
  std::list<Vertex_handle> mins;
  for (Minimals_iterator m = minimals_begin(tr); m != minimals_end(tr); ++m)
    mins.push_back(&(*m));
  typename std::list<Vertex_handle>::iterator clic = mins.begin();
  for ( ; clic != mins.end() ; ++clic) 
    if (is_minimal(*clic,tr)) sweep_good(*clic,tr);
}


template < class Gtr_ , class It , class Flip >
template < class Tr >
typename Antichain<Gtr_,It,Flip>::Vertex_handle
Antichain<Gtr_,It,Flip>::
compute_phi(Vertex_handle v, Tr tr) //const
{
  // -------------------------------------------------------------------------
  typename Tr::Sup sup; 
  typename Tr::Inf inf; 
  typename Tr::Set_sup set_sup; typename Tr::Set_inf set_inf;
  typename Tr::Splice splice;
  typename Tr::Cw_source_edge   cw_source_edge;
  typename Tr::Cw_target_edge   cw_target_edge;
  typename Tr::Ccw_source_edge   ccw_source_edge;
  typename Tr::Ccw_target_edge   ccw_target_edge;
  typename Tr::CcL ccL; typename Tr::CcR ccR; 
  typename Tr::Is_left_xx is_left_xx;
  // -------------------------------------------------------------------------
  Vertex_handle phiv; // The new Vertex that we must compute
  // -------------------------------------------------------------------------
  // Vertex has already been swept.
  if (sup(v) != 0) {
    phiv  = sup(sup(v)); 
  }
  // -------------------------------------------------------------------------
  // The pi of the Vertex has already been swept. We use the formula:
  // phi(v)->pi() = phi(v->pi()).
  else if (!linear_space()&&
           sup(v->pi()) != 0) {
    if (is_on_convex_hull(v)) {
      if (is_left_xx(v)  && ccL(v->pi()) != 0)
        ccR(v)->set_pi(ccL(v->pi()));
      if (!is_left_xx(v) && ccR(v->pi()) != 0) 
        ccL(v)->set_pi(ccR(v->pi()));
    }
    phiv = sup(sup(v->pi()))->pi();
  } else if (linear_space()&&is_on_convex_hull(v)) {
    phiv=v;
    do
      if (is_left_xx(phiv)) phiv=ccR(phiv); else phiv=ccL(phiv);
    while (phiv->is_constraint());
    phiv=phiv->pi();
  }
  // -------------------------------------------------------------------------
  // We flip v in the current pseudo-triangulation. 
  // To compute phi(v) we walk on the incident pseudo-triangles. 
  //else phiv = Chi2_strategy(this)(v,tr);
  else phiv = Flip_traits(this)(v,tr);
  // -------------------------------------------------------------------------
  // We splice the arcs left and right if b is not on the convex hull
  //if (!is_on_convex_hull(v)) { splice(right,phiv); splice(left,phiv); }
  // -------------------------------------------------------------------------
  // Creating a new face with source v and sink phi(v)
  if (sup(v) == 0) set_sup(v,new Face);
  set_inf(sup(v),v);
  set_sup(sup(v),phiv);
  // -------------------------------------------------------------------------
  // Creating a new face with source pi(v) and sink pi(phi(v))
  if (!linear_space()) {
    Vertex_handle vv=inf(cw_source_edge(phiv));
    if (cw_source_edge(phiv)==ccw_target_edge(vv)) {
      splice(ccw_source_edge(vv->pi()),phiv->pi());
    } else {
      CGAL_assertion(cw_source_edge(phiv)==ccw_source_edge(vv));
      splice(ccw_target_edge(vv->pi()),phiv->pi());      
    }
    vv=inf(cw_target_edge(phiv));
    if (cw_target_edge(phiv)==ccw_target_edge(vv)) {
      splice(ccw_source_edge(vv->pi()),phiv->pi());
    } else {
      CGAL_assertion(cw_target_edge(phiv)==ccw_source_edge(vv));
      splice(ccw_target_edge(vv->pi()),phiv->pi());      
    }
    if (sup(v->pi()) == 0)  set_sup(v->pi(),new Face);
    set_inf(sup(v->pi()),v->pi());
    set_sup(sup(v->pi()),phiv?phiv->pi():0);                        
  }
  // -------------------------------------------------------------------------
  CGAL_precondition(sup(v) != 0 && (linear_space()||sup(sup(v)) != 0));
  return phiv;      
}


template < class Gtr_ , class It , class Flip >
void
Antichain<Gtr_,It,Flip>::set_constraint(Vertex_handle v) 
{
  typedef Ccw_traits Tr;
  typename Tr::Set_sup set_sup; typename Tr::Set_inf set_inf;
  typename Tr::Set_adjacent_faces_one_to_one set_adjacent_faces;
  typename Tr::Set_target_cusp_edge set_target_cusp_edge;
  typename Tr::Set_source_cusp_edge set_source_cusp_edge;
  typename Tr::Target_cusp_edge target_cusp_edge;
  typename Tr::Source_cusp_edge source_cusp_edge;
  typename Tr::Is_left_xx is_left_xx;
  typename Tr::Is_xx_left is_xx_left;
  // -------------------------------------------------------------------------
  typename Tr::Sup sup; typename Tr::Inf inf;
  Face_handle fs = sup(v);         Face_handle fi = inf(v);
  Face_handle pifs = sup(v->pi()); Face_handle pifi = inf(v->pi());
  // -------------------------------------------------------------------------
  // New edge for rays emanating from b->target() where b = this
  // with angle \theta such that \theta(b) - \pi <  \theta < \theta(b)
  // New edge for rays emanating from b->source() where b = this
  // with angle \theta such that \theta(b) - \pi <  \theta < \theta(b)
  set_target_cusp_edge(v,new Edge);
  set_source_cusp_edge(v,new Edge); 
  target_cusp_edge(v)->set_sign(false);
  set_sup(target_cusp_edge(v),v);
  set_inf(target_cusp_edge(v),v->pi());
  source_cusp_edge(v)->set_sign(true);
  set_sup(source_cusp_edge(v),v);
  set_inf(source_cusp_edge(v),v->pi());
  // -------------------------------------------------------------------------
  // New face directed from b->target_object() to b, the sink is b
  Face_handle f0 = new Face; 
  set_sup(f0,v);
  set_inf(f0,v->pi());
  f0->set_top_edge(target_cusp_edge(v)); 
  if (is_xx_left(v)) set_adjacent_faces(target_cusp_edge(v),f0,0,0);
  else               set_adjacent_faces(target_cusp_edge(v),0,f0,0);
  // New face directed from b->source_object() to b, the sink is b
  Face_handle f1 = new Face; 
  set_sup(f1,v);
  set_inf(f1,v->pi());
  f1->set_bottom_edge(source_cusp_edge(v)); 
  if (is_left_xx(v)) set_adjacent_faces(source_cusp_edge(v),0,0,f1);
  else               set_adjacent_faces(source_cusp_edge(v),0,f1,0);
  // -------------------------------------------------------------------------
  set_sup(v,fs);         set_inf(v,fi);
  set_sup(v->pi(),pifs); set_inf(v->pi(),pifi);
  // -------------------------------------------------------------------------
  // Do the same with b->pi()
  //if (!v->pi()->is_constraint()) set_constraint(v->pi());
  // -------------------------------------------------------------------------
}


template < class Gtr_ , class It , class Flip >
void
Antichain<Gtr_,It,Flip>::remove_constraint(Vertex_handle v) 
{
  if (!v->is_constraint()) return;
  // -------------------------------------------------------------------------
  //  erase(v->target_cusp_edge());
  //    erase(v->source_cusp_edge());
  unset_constraint(v);
  // -------------------------------------------------------------------------
}


template < class Gtr_ , class It , class Flip >
void
Antichain<Gtr_,It,Flip>::add_constraint(Vertex_handle v) 
{
  if (v->is_constraint()) return;
  // -------------------------------------------------------------------------
  set_constraint(v);
  set_constraint(v->pi());
  //push_back(*v->target_cusp_edge());
  //push_back(*v->source_cusp_edge());
  // -------------------------------------------------------------------------
}


template < class Gtr_ , class It , class Flip >
void
Antichain<Gtr_,It,Flip>::unset_constraint(Vertex_handle v) 
{
  // -------------------------------------------------------------------------
  delete v->target_cusp_face(); delete v->source_cusp_face();
  delete v->target_cusp_edge(); delete v->source_cusp_edge();
  //    if (v->pi()->is_constraint()) unset_constraint(v->pi());
  // -------------------------------------------------------------------------
}


template < class Gtr_ , class It , class Flip >
template < class Tr>
void
Antichain<Gtr_,It,Flip>::
sweep_regular(const Vertex_handle& v, const Tr&,bool opposite) const
{
  // -------------------------------------------------------------------------
  // The operators used by this method
  typename Tr::Sup sup; typename Tr::Inf inf; 
  typename Tr::Set_inf set_inf; 
  typename Tr::Set_adjacent_faces_one_to_one set_adjacent_faces_one_to_one;
  typename Tr::Set_adjacent_faces set_adjacent_faces;
  typename Tr::Dl dl; typename Tr::Dr dr;
  typename Tr::Ul ul; typename Tr::Ur ur;
  typename Tr::Is_left_xx is_left_xx;
  typename Tr::Cw_source_edge cw_source_edge;
  typename Tr::Ccw_source_edge ccw_source_edge;
  typename Tr::Cw_target_edge cw_target_edge;
  typename Tr::Ccw_target_edge ccw_target_edge;
  // -------------------------------------------------------------------------
  // The four edges adjacent to v
  Edge_handle e  = cw_source_edge(v);
  Edge_handle ep = ccw_source_edge(v);
  Edge_handle f  = cw_target_edge(v);
  Edge_handle fp = ccw_target_edge(v);
  // -------------------------------------------------------------------------
  // The six faces adjacent to v.
  Face_handle f0,f1,f2,f3;
  if (e->sign()) { f0 = ul(e); f1 = dl(e); }
  else { 
    f0 = dl(e); if (f0 == 0 || sup(f0) == v) f0 = dl(f); 
    f1 = dr(e); 
  }
  if (f->sign()) { 
    f2 = ur(f); if (f2 == 0 || sup(f2) == v) f2 = ur(e); 
    f3 = ul(f); 
  }
  else           { f2 = dr(f); f3 = ul(f); }
  Face_handle f4 = sup(v); set_inf(f4,v);
  Face_handle f5 = inf(v);
  // -------------------------------------------------------------------------
  // Due to our identification pi^2(v) == v, a face can appear twice in the
  // antichain. As a consequence a face can appear twice in the antichain.
  // This implies that the pointers Edge <---> Face are not necessarily
  // reversible. 
  // Our solution: make the Edge --> Face pointers correct and use the
  // infinite face and the two Edges below.
  Edge_handle botf4 = f4->bottom_edge();
  Edge_handle topf4 = f4->top_edge();
  // -------------------------------------------------------------------------
  // Updating the Edge <--> Face pointers for the two new edges of the
  // antichain.
  if (0&&linear_space() && !is_on_convex_hull(v)) {
    set_adjacent_faces_one_to_one(e,0,0,0);
    set_adjacent_faces_one_to_one(f,0,0,0);
  }
  if (opposite) {
    if (is_left_xx(v)) set_adjacent_faces(ep,f4,f3,f0);
    else               set_adjacent_faces(ep,f0,f4,f3);
    set_adjacent_faces(fp,f1,f2,f4);
    return;
  } else {
    if (is_left_xx(v)) set_adjacent_faces_one_to_one(ep,f4,f3,f0);
    else               set_adjacent_faces_one_to_one(ep,f0,f4,f3);
    set_adjacent_faces_one_to_one(fp,f1,f2,f4);
  }
  // -------------------------------------------------------------------------
  // Fix the pointers due to the problem mentionned above. This only occurs
  // for faces whose sink are on the convex-hull.
  CGAL_precondition(f4 != 0);
  CGAL_precondition(f5 != 0);
  if (is_on_convex_hull(v)) {
    if (f2 == f4) { f4->set_top_edge(ep); f4->set_bottom_edge(fp); }
    else {
      Face_handle w = (is_left_xx(v)) ? f2 : f0;
      if (botf4 != 0 && w->bottom_edge() != botf4) 
        f4->set_bottom_edge(botf4);
      if (topf4 != 0 && w->top_edge()    != topf4) 
        f4->set_top_edge(topf4);
    }
    if (f5->bottom_edge() != e) {
      f5->set_top_edge(infinite_face()->bottom_edge());
      return;
    }
    else if (f5->top_edge() != f) {
      f5->set_bottom_edge(infinite_face()->top_edge());
      return;
    }
  }

//   // The face inf(v) is no longer swept.
//   f5->set_bottom_edge(0);
//   f5->set_top_edge   (0);

}

// -----------------------------------------------------------------------------

template < class Gtr_ , class It , class Flip >
template < class Tr>
void
Antichain<Gtr_,It,Flip>::
sweep_constraint(const Vertex_handle& v, const Tr&,bool opposite) const
{
  // -------------------------------------------------------------------------
  // The operators used by this method
  typename Tr::Sup sup; 
  typename Tr::Set_adjacent_faces set_adjacent_old_faces;
  typename Tr::Set_adjacent_faces_one_to_one set_adjacent_faces;
  typename Tr::Dl dl; typename Tr::Dr dr;
  typename Tr::Ul ul; typename Tr::Ur ur;
  typename Tr::Target_cusp_face target_cusp_face;
  typename Tr::Source_cusp_face source_cusp_face;
  typename Tr::Is_left_xx is_left_xx; typename Tr::Is_xx_left is_xx_left;
  typename Tr::Cw_source_edge cw_source_edge;
  typename Tr::Ccw_source_edge ccw_source_edge;
  typename Tr::Cw_target_edge cw_target_edge;
  typename Tr::Ccw_target_edge ccw_target_edge;
  // -------------------------------------------------------------------------
  // The four edges adjacent to v
  Edge_handle e  = cw_source_edge(v);
  Edge_handle ep = ccw_source_edge(v);
  Edge_handle f  = cw_target_edge(v);
  Edge_handle fp = ccw_target_edge(v);
  // -------------------------------------------------------------------------
  // The four regular faces adjacent to v
  Face_handle f0,f1,f2,f3;
  if (e->sign()) { f0 = ul(e);                               f1 = dl(e); }
  else           { f0 = dl(e); if (sup(f0) == v) f0 = dl(f); f1 = dr(e); }
  if (f->sign()) { f2 = ur(f); if (sup(f2) == v) f2 = ur(e); f3 = ul(f); }
  else           { f2 = dr(f);                               f3 = ul(f); }
  // -------------------------------------------------------------------------
  // The four degenerate faces adjacent to v
  Face_handle a  = source_cusp_face(v);
  Face_handle ap = source_cusp_face(v->pi());
  Face_handle b  = target_cusp_face(v);
  Face_handle bp = target_cusp_face(v->pi());
  // -------------------------------------------------------------------------
  Edge_handle fix = 0;
  if (is_on_convex_hull(v)) {
    if (e->sign() && e != f2->bottom_edge())    fix = f2->bottom_edge();
    else if (!e->sign() && f != f0->top_edge()) fix = f0->top_edge();
  }
  // -------------------------------------------------------------------------
  // Updating the Edge <--> Face pointers for the two new edges of the
  // antichain.
  typename Vertex::Type_util type;
  switch (type(is_left_xx(v),is_xx_left(v))) {
  case Vertex::LL:
    set_adjacent_old_faces(e,f1,f2,f0);
    set_adjacent_old_faces(f,a,b,f3);
    set_adjacent_faces(ep,ap,f3,bp);
    set_adjacent_faces(fp,f1,f2,f0);
    break;
  case Vertex::RR:
    set_adjacent_old_faces(e,a,f1,b);
    set_adjacent_old_faces(f,f0,f2,f3);
    set_adjacent_faces(ep,f0,f2,f3);
    set_adjacent_faces(fp,f1,ap,bp);
    break;
  case Vertex::LR:
    set_adjacent_old_faces(e,f1,b,f0);
    set_adjacent_old_faces(f,a,f2,f3);
    set_adjacent_faces(ep,f2,f3,bp);
    set_adjacent_faces(fp,f1,ap,f0);
    break;
  case Vertex::RL:
    set_adjacent_old_faces(e,a,f1,f2);
    set_adjacent_old_faces(f,f0,b,f3);
    set_adjacent_faces(ep,f0,ap,f3);
    set_adjacent_faces(fp,f1,f2,bp);
    break;
  }
  // -------------------------------------------------------------------------
  if (is_on_convex_hull(v) && fix != 0) {
    if (fix->sign()) set_adjacent_faces(fix,dl(fix),ur(fix),ul(fix));
    else             set_adjacent_faces(fix,dl(fix),dr(fix),ul(fix));
  }
  // -------------------------------------------------------------------------
//   if (a != 0) a->set_top_edge(0);
//   if (b != 0) b->set_bottom_edge(0);
  // -------------------------------------------------------------------------
}



template < class Gtr_ , class It , class Flip >
template < class Tr>
bool
Antichain<Gtr_,It,Flip>::
is_swept_regular(const Vertex_handle& v, Tr /*tr*/) const
{
  // -------------------------------------------------------------------------
  CGAL_precondition(!v->is_constraint());
  // -------------------------------------------------------------------------
  // The operators used by this method
  typename Tr::Sup sup; 
  typename Tr::Dl dl; typename Tr::Dr dr;
  typename Tr::Ul ul; typename Tr::Ur ur;
  typename Tr::Cw_source_edge  cw_source_edge;
  typename Tr::Cw_target_edge  cw_target_edge;
  typename Tr::Ccw_source_edge ccw_source_edge;
  typename Tr::Ccw_target_edge ccw_target_edge;
  typename Tr::Is_left_xx is_left_xx; typename Tr::Is_xx_left is_xx_left;
  // -------------------------------------------------------------------------
  // A vertex v has been swept iff. the two following conditions are
  // satisfied:
  // (1) sup(v) is not 0
  // (2) the faces adjacent to ccw_source_edge(v) and ccw_target_edge(v) are
  // correct. We check this with the help of the faces adjacent to
  // cw_source_edge(v) and cw_target_edge(v) which are assumed to be correct.
  // -------------------------------------------------------------------------
  // Checking (1).
  Face_handle f4 = sup(v); 
  if (f4 == 0) return false;
  // -------------------------------------------------------------------------
  // Checking (2).
  // -------------------------------------------------------------------------
  // The four edges adjacent to v
  Edge_handle e  = cw_source_edge(v);
  Edge_handle ep = ccw_source_edge(v);
  Edge_handle f  = cw_target_edge(v);
  Edge_handle fp = ccw_target_edge(v);
  if (e == 0 || ep == 0 || f == 0 || fp == 0) return false;
  // -------------------------------------------------------------------------
  // The remaining faces adjacent to v different from inf(v).
  Face_handle f0,f1,f2,f3;
  if (e->sign()) { f0 = ul(e); f1 = dl(e); }
  else { 
    f0 = dl(e); if (f0 == 0 || sup(f0) == v) f0 = dl(f); 
    f1 = dr(e); 
  }
  if (f->sign()) { 
    f2 = ur(f); if (f2 == 0 || sup(f2) == v) f2 = ur(e); 
    f3 = ul(f); 
  }
  else           { f2 = dr(f); f3 = ul(f); }
  // -------------------------------------------------------------------------
  if (is_left_xx(v)) {
    if (f4 != dl(ep) || f3 != ur(ep) || f0 != ul(ep)) return false;
  }
  else {
    if (f0 != dl(ep) || f4 != dr(ep) || f3 != ul(ep)) return false;
  }
  if (is_xx_left(v)) {
    if (f1 != dl(fp) || f2 != ur(fp) || f4 != ul(fp)) return false;
  }
  else {
    if (f1 != dl(fp) || f2 != dr(fp) || f4 != ul(fp)) return false;
  }
  // -------------------------------------------------------------------------
  return true;
}


template < class Gtr_ , class It , class Flip >
template < class Tr>
bool
Antichain<Gtr_,It,Flip>::
is_swept_constraint(const Vertex_handle& v, Tr /*tr*/) const
{
  // -------------------------------------------------------------------------
  CGAL_precondition(v->is_constraint());
  // -------------------------------------------------------------------------
  // The operators used by this method
  typename Tr::Sup sup; 
  typename Tr::Dl dl; typename Tr::Dr dr;
  typename Tr::Ul ul; typename Tr::Ur ur;
  typename Tr::Cw_source_edge  cw_source_edge;
  typename Tr::Cw_target_edge  cw_target_edge;
  typename Tr::Ccw_source_edge ccw_source_edge;
  typename Tr::Ccw_target_edge ccw_target_edge;
  typename Tr::Target_cusp_face target_cusp_face;
  typename Tr::Source_cusp_face source_cusp_face;
  typename Tr::Is_left_xx is_left_xx; typename Tr::Is_xx_left is_xx_left;
  // -------------------------------------------------------------------------
  // See is_swept_regular for the criterion. As sup(v) is not defined we check
  // only point (2).
  // -------------------------------------------------------------------------
  // The four edges adjacent to v
  Edge_handle e  = cw_source_edge(v);
  Edge_handle ep = ccw_source_edge(v);
  Edge_handle f  = cw_target_edge(v);
  Edge_handle fp = ccw_target_edge(v);
  if (e == 0 || ep == 0 || f == 0 || fp == 0) return false;
  // -------------------------------------------------------------------------
  // The four regular faces adjacent to v
  Face_handle f0,f1,f2,f3;
  if (e->sign()) { 
    if (ul(e) == 0 || ur(e) == 0) return false;
    f0 = ul(e); f1 = dl(e); 
  }
  else { 
    if (dl(e) == 0 || dr(e) == 0) return false;
    f0 = dl(e); if (sup(f0) == v) f0 = dl(f); f1 = dr(e); 
  }
  if (f->sign()) { 
    if (ul(f) == 0 || ur(f) == 0) return false;
    f2 = ur(f); if (sup(f2) == v) f2 = ur(e); f3 = ul(f); 
  }
  else { 
    if (dl(f) == 0 || dr(f) == 0) return false;
    f2 = dr(f); f3 = ul(f); 
  }
  // -------------------------------------------------------------------------
  // The four degenerate faces adjacent to v
  Face_handle a  = source_cusp_face(v);
  Face_handle b  = target_cusp_face(v);
  Face_handle ap = source_cusp_face(v->pi());
  Face_handle bp = target_cusp_face(v->pi());
  if (a == 0 || b == 0 || ap == 0 || bp == 0) return false;
  // -------------------------------------------------------------------------
  // Updating the Edge <--> Face pointers for the two new edges of the
  // antichain.
  typename Vertex::Type_util type;
  switch (type(is_left_xx(v),is_xx_left(v))) {
  case Vertex::LL:
    if (ap != dl(ep) || f3 != ur(ep) || bp != ul(ep)) return false;
    if (f1 != dl(fp) || f2 != ur(fp) || f0 != ul(fp)) return false;
    break;
  case Vertex::RR:
    if (f0 != dl(ep) || f2 != dr(ep) || f3 != ul(ep)) return false;
    if (f1 != dl(fp) || ap != dr(fp) || bp != ul(fp)) return false;
    break;
  case Vertex::LR:
    if (f2 != dl(ep) || f3 != ur(ep) || bp != ul(ep)) return false;
    if (f1 != dl(fp) || ap != dr(fp) || f0 != ul(fp)) return false;
    break;
  case Vertex::RL:
    if (f0 != dl(ep) || ap != dr(ep) || f3 != ul(ep)) return false;
    if (f1 != dl(fp) || f2 != ur(fp) || bp != ul(fp)) return false;
    break;
  }
  return true;
}


template < class Gtr_ , class It , class Flip >
template < class Tr>
bool
Antichain<Gtr_,It,Flip>::
is_swept(const Vertex_handle& v, Tr tr) const
{
  if (!v->is_constraint()) return is_swept_regular(v,tr);
  return (is_swept_constraint(v,tr) && is_swept_constraint(v->pi(),tr)); 
}

}

CGAL_END_NAMESPACE

#endif // CGAL_VISIBILITY_COMPLEX_2_ANTICHAIN_H
