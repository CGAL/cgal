// Copyright (c) 2012  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Andreas Fabri, Olivier Billet, Mariette Yvinec

#ifndef CGAL_POLYLINE_CONSTRAINT_HIERARCHY_2_H
#define CGAL_POLYLINE_CONSTRAINT_HIERARCHY_2_H

#include <CGAL/basic.h>
#include <utility>
#include <map>
#include <set> 
#include <list> 
#include <CGAL/Skiplist.h>
#include <CGAL/Iterator_project.h>
#include <CGAL/triangulation_assertions.h>

namespace CGAL {

// T               is expected to be Vertex_handle
// Data            is intended to store info on a Vertex
template <class T, class Data>
class Polyline_constraint_hierarchy_2
{
public:
  typedef Data                                    Point;
  typedef std::pair<T, T>                         H_edge;
  typedef T                                       H_vertex;
  typedef T                                       Vertex_handle;
  typedef Polyline_constraint_hierarchy_2<T,Data> Hierarchy;
  typedef std::pair<T, T>                         H_constraint;

  struct Node {
    T vertex;
    Point point;
    int id;

    explicit Node(T vh)
      : vertex(vh), point(vh->point())
    {}
  };

  typedef CGAL::Skiplist<Node>                  H_vertex_list;
  
  // only nodes with a vertex_handle that is still the triangulation
  typedef typename H_vertex_list::skip_iterator H_vertex_it; 
  // all nodes
  typedef typename H_vertex_list::all_iterator  H_all_iterator_it;
  
  typedef std::list<H_constraint>              H_constraint_list;
  typedef typename H_constraint_list::iterator H_constraint_it;
  
  typedef H_vertex_list* Constraint_id;

  class H_context {
    friend class Polyline_constraint_hierarchy_2<T,Data>;
  private:
    H_vertex_list*    enclosing;
    H_vertex_it       pos;
  public:
    H_context()
    {}

    H_context(const H_context& hc)
      : enclosing(hc.enclosing), pos(hc.pos)
    {}

    H_vertex_it    vertices_begin() { return enclosing->begin();}
    H_vertex_it    current() {return pos;}
    H_vertex_it    vertices_end() {return enclosing->end();}
    Constraint_id  id() { return enclosing; }
    std::size_t    number_of_vertices() {return enclosing->size();}
  };                                           

  typedef std::list<H_context>              H_context_list;
  typedef typename H_context_list::iterator H_context_iterator;

  typedef std::set<H_vertex_list*>                      H_constraint_set;
  typedef std::map<H_edge,   H_context_list* >          H_sc_to_c_map;
  typedef typename H_constraint_set::const_iterator     H_c_iterator;
  typedef typename H_sc_to_c_map::const_iterator        H_sc_iterator;
  typedef std::pair<H_constraint, H_vertex_list*>       H_c_value;
  typedef std::pair<H_edge,   H_context_list*>          H_sc_value;
  
private:
  // data for the 1d hierarchy
  H_constraint_set   constraint_set;
  H_sc_to_c_map   sc_to_c_map;

  
public:
  Polyline_constraint_hierarchy_2() { }
  Polyline_constraint_hierarchy_2(const Polyline_constraint_hierarchy_2& ch); 
  ~Polyline_constraint_hierarchy_2(){ clear();}
  void clear();
  Polyline_constraint_hierarchy_2& operator=(const Polyline_constraint_hierarchy_2& ch);

  // Query 
  bool is_subconstrained_edge(T va, T vb) const;
  bool is_constrained_edge(T va, T vb) const;
  bool is_constrained_vertex(T v) const;
  bool vertices_in_constraint(H_constraint hc,  
			      H_vertex_it& v_first,
			      H_vertex_it& v_past) const;
  H_vertex_it vertices_in_constraint_begin(Constraint_id) const;
  H_vertex_it vertices_in_constraint_end(Constraint_id) const;

  bool enclosing_constraint(H_edge he, H_constraint& hc) const;
  bool enclosing_constraint(T  vaa, T  vbb, T& va, T& vb) const;
  bool enclosing_constraints(T vaa, T  vbb,  H_constraint_list& hcl) const;
  bool next_along_sc(T va, T vb, T& w) const;
  void oriented_end(T va, T vb, T& vc) const;

  H_context context(T va, T vb);
  std::size_t number_of_enclosing_constraints(T va, T vb);
  H_context_iterator contexts_begin(T va, T vb);
  H_context_iterator contexts_end(T va, T vb);
  std::size_t number_of_constraints() { return constraint_set.size();}
  std::size_t number_of_subconstraints() {return sc_to_c_map.size();}
  

  // insert/remove
  void add_Steiner(T va, T vb, T vx);
  H_vertex_list* insert_constraint(T va, T vb);
  void append_constraint(H_vertex_list* vl, T va, T vb);
  void swap(Constraint_id first, Constraint_id second);
  void remove_constraint(Constraint_id cid);
  void split_constraint(T va, T vb, T vc);

  void simplify(H_vertex_it u,
                H_vertex_it v,
                H_vertex_it w);

  int remove_points_from_constraint(Constraint_id);
  int remove_points_from_constraints();

  Constraint_id concatenate(Constraint_id first, Constraint_id second);
  Constraint_id concatenate2(Constraint_id first, Constraint_id second);
  Constraint_id split(Constraint_id first, H_vertex_it vcit);
  Constraint_id split2(Constraint_id first, H_vertex_it vcit);

  void constrain_vertex(T v, Data data=Data());
  void unconstrain_vertex(T v);
  void set_data(T v, Data data);
  Data get_data(T v);

  void remove_Steiner(T v, T va, T vb);

  // iterators
  H_sc_iterator sc_begin() const{ return sc_to_c_map.begin(); }
  H_sc_iterator sc_end()   const{ return sc_to_c_map.end();   }
  H_c_iterator  c_begin()  const{ return constraint_set.begin(); }
  H_c_iterator  c_end()    const{ return constraint_set.end();   }
  
  //Helping functions
  void copy(const Polyline_constraint_hierarchy_2& ch);
  void copy(const Polyline_constraint_hierarchy_2& ch, std::map<Node,Node>& vmap);
  void swap(Polyline_constraint_hierarchy_2& ch);

private: 
  H_edge    make_edge(T va, T vb) const;
  H_vertex_it     get_pos(T va, T vb) const;
  bool      get_contexts(T va, T vb, 
			 H_context_iterator& ctxt, 
			 H_context_iterator& past) const;

  bool      get_contexts(T va, T vb, H_context_list*&) const;

  //to_debug
public:
  void   print() const;
};

template <class T, class Data> 
Polyline_constraint_hierarchy_2<T,Data>::
Polyline_constraint_hierarchy_2(const Polyline_constraint_hierarchy_2& ch)
{
  copy(ch);
}

template <class T, class Data> 
Polyline_constraint_hierarchy_2<T,Data>&
Polyline_constraint_hierarchy_2<T,Data>::
operator=(const Polyline_constraint_hierarchy_2& ch){
  copy(ch);
  return *this;
}

template <class T, class Data> 
void
Polyline_constraint_hierarchy_2<T,Data>::
copy(const Polyline_constraint_hierarchy_2& ch1)
{
  // create a identity transfer vertex map
  std::map<Node, Node>  vmap;
  H_c_iterator cit1 = ch1.c_begin();
  for( ; cit1 != ch1.c_end(); ++cit1) {
    H_vertex_it vit = cit1->second->begin();
    for( ; vit != cit1->second->end(); ++vit) {
      vmap[*vit] = *vit;
    }
  }
  copy(ch1, vmap);
}

template <class T, class Data> 
void
Polyline_constraint_hierarchy_2<T,Data>::
copy(const Polyline_constraint_hierarchy_2& ch1, std::map<Node,Node>& vmap)
  // copy with a transfer vertex map
{
  std::cerr << "copy" << std::endl;
  std::map<H_vertex_list*,H_vertex_list*> vlmap;
  clear();
  // copy constraint_set
  H_c_iterator cit1 = ch1.c_begin();
  for( ; cit1 != ch1.c_end(); ++cit1) {
    H_vertex_list* hvl1 = *cit1;
    H_vertex_list* hvl2 = new H_vertex_list;
    vlmap[hvl1] = hvl2;
    H_vertex_it vit = hvl1->begin();
    for( ; vit != hvl1->end(); ++vit) hvl2->push_back(vmap[*vit]);
    constraint_set.insert(hvl2);
  }
  // copy sc_to_c_map
  H_sc_iterator scit1 = ch1.sc_begin();
  for( ; scit1 != ch1.sc_end(); ++scit1) {
    //vertices of the subconstraints
    H_vertex uu2 = vmap[scit1->first.first];
    H_vertex vv2 = vmap[scit1->first.second];
    H_context_list* hcl1  = scit1->second;
    H_context_list* hcl2  = new H_context_list;
    H_context_iterator cit1 = hcl1->begin();
    for( ; cit1 != hcl1->end(); ++cit1){
      // vertices of the enclosing constraints
      H_context ctxt2;
      ctxt2.enclosing = vlmap[cit1->enclosing];
      ctxt2.pos = ctxt2.enclosing->begin();
      H_vertex_it aux = cit1->enclosing->begin();
      while( aux != cit1->pos) {
	++aux;
	++ctxt2.pos;
      }
      hcl2->push_back(ctxt2);
    }
    sc_to_c_map[make_edge(uu2,vv2)] = hcl2;
  }

  return;
}


template <class T, class Data> 
void
Polyline_constraint_hierarchy_2<T,Data>::
swap(Polyline_constraint_hierarchy_2& ch)
{
  std::cerr << "swap" << std::endl;
  constraint_set.swap(ch.constraint_set);
  sc_to_c_map.swap(ch.sc_to_c_map);
}


/*
template <class T, class Data> 
bool Polyline_constraint_hierarchy_2<T,Data>::
is_constrained_edge(T va, T vb) const
{
  return( c_to_sc_map.find(make_edge(va, vb)) != c_to_sc_map.end() );
}
*/

template <class T, class Data> 
bool Polyline_constraint_hierarchy_2<T,Data>::
is_subconstrained_edge(T va, T vb) const
{
  return( sc_to_c_map.find(make_edge(va, vb)) != sc_to_c_map.end() );
}


template <class T, class Data> 
bool Polyline_constraint_hierarchy_2<T,Data>::
vertices_in_constraint(H_constraint hc, 
		       H_vertex_it& v_first,
		       H_vertex_it& v_past ) const
{
  H_sc_iterator sc_iter = sc_to_c_map.find(hc);
  if( sc_iter == sc_to_c_map.end() )
    return false;
  v_first = (*sc_iter).second;
  return true;
}

// af: obsolete
template <class T, class Data>
bool Polyline_constraint_hierarchy_2<T,Data>::
enclosing_constraint(H_edge he, H_constraint& hc) const
{
  H_context_iterator hcit, past;
  if ( !get_contexts(he.first,he.second, hcit ,past)) return false;
  hc = make_edge(hcit->enclosing->front(), hcit->enclosing->back());
  return true;
}


// used by Constrained_triangulation_plus_2::intersect with Exact_intersection_tag
template <class T, class Data>
bool Polyline_constraint_hierarchy_2<T,Data>::
enclosing_constraint(T  vaa, T  vbb, T& va, T& vb) const
{
  H_context_iterator hcit, past;
  if ( !get_contexts(vaa,vbb, hcit ,past)) return false;
  va = hcit->enclosing->front();
  vb = hcit->enclosing->back();
  return true;
}

// af: obsolete
template <class T, class Data>
bool Polyline_constraint_hierarchy_2<T,Data>::
enclosing_constraints(T vaa, T vbb , H_constraint_list& hcl) const
{
  H_context_iterator hcit, past;
  if ( !get_contexts(vaa,vbb, hcit ,past)) return false;
  for (; hcit!=past; hcit++) {
    hcl.push_back(make_edge(hcit->enclosing->front(), 
			    hcit->enclosing->back())); 
  }
  return true;
}

template <class T, class Data>
typename Polyline_constraint_hierarchy_2<T,Data>::H_context
Polyline_constraint_hierarchy_2<T,Data>::
context(T va, T vb)
{
  H_context_iterator hcit, past;
  if(!get_contexts(va,vb, hcit ,past)) CGAL_triangulation_assertion(false);
  return *hcit;
}

template <class T, class Data>
std::size_t 
Polyline_constraint_hierarchy_2<T,Data>::
number_of_enclosing_constraints(T va, T vb)
{
  H_context_list* hcl;
  if (! get_contexts(va,vb,hcl)) CGAL_triangulation_assertion(false);
  return hcl->size();
}

template <class T, class Data>
typename Polyline_constraint_hierarchy_2<T,Data>::H_context_iterator
Polyline_constraint_hierarchy_2<T,Data>::
contexts_begin(T va, T vb)
{
   H_context_iterator first, last;
   if( !get_contexts(va,vb,first,last))  CGAL_triangulation_assertion(false);
   return first;
}

template <class T, class Data>
typename Polyline_constraint_hierarchy_2<T,Data>::H_context_iterator
Polyline_constraint_hierarchy_2<T,Data>::
contexts_end(T va, T vb)
{   
   H_context_iterator first, last;
   if( !get_contexts(va,vb,first,last))  CGAL_triangulation_assertion(false);
   return last;
} 

template <class T, class Data>
typename Polyline_constraint_hierarchy_2<T,Data>::H_vertex_it
Polyline_constraint_hierarchy_2<T,Data>::
vertices_in_constraint_begin(Constraint_id cid) const
{
  return cid->skip_begin();
}
  
template <class T, class Data>
typename Polyline_constraint_hierarchy_2<T,Data>::H_vertex_it
Polyline_constraint_hierarchy_2<T,Data>::
vertices_in_constraint_end(Constraint_id cid) const
{
  return cid->skip_end();
}

template <class T, class Data>
void
Polyline_constraint_hierarchy_2<T,Data>::
swap(Constraint_id first, Constraint_id second){
    // We have to look at all subconstraints
  for(H_vertex_it it = first->begin(), succ = it; 
      ++succ != first->end(); 
      ++it){
    typename H_sc_to_c_map::iterator scit = sc_to_c_map.find(make_edge(*it,*succ));
    CGAL_triangulation_assertion(scit != sc_to_c_map.end());
    H_context_list* hcl = scit->second;

    // and replace the context of the constraint
    for(H_context_iterator ctit=hcl->begin(); ctit != hcl->end(); ctit++) {
      if(ctit->enclosing == first){
	ctit->enclosing = 0;
	break;
      }
    }
  }
    // We have to look at all subconstraints
  for(H_vertex_it it = second->begin(), succ = it; 
      ++succ != second->end(); 
      ++it){
    typename H_sc_to_c_map::iterator scit = sc_to_c_map.find(make_edge(*it,*succ));
    CGAL_triangulation_assertion(scit != sc_to_c_map.end());
    H_context_list* hcl = scit->second;

    // and replace the context of the constraint
    for(H_context_iterator ctit=hcl->begin(); ctit != hcl->end(); ctit++) {
      if(ctit->enclosing == second){
	ctit->enclosing = first;
	break;
      }
    }
  }   
  // We have to look at all subconstraints
  for(H_vertex_it it = first->begin(), succ = it; 
      ++succ != first->end(); 
      ++it){
    typename H_sc_to_c_map::iterator scit = sc_to_c_map.find(make_edge(*it,*succ));
    CGAL_triangulation_assertion(scit != sc_to_c_map.end());
    H_context_list* hcl = scit->second;

    // and replace the context of the constraint
    for(H_context_iterator ctit=hcl->begin(); ctit != hcl->end(); ctit++) {
      if(ctit->enclosing == 0){
	ctit->enclosing = second;
	break;
      }
    }
  }
  first->swap(*second);
}


template <class T, class Data>
void
Polyline_constraint_hierarchy_2<T,Data>::
remove_constraint(Constraint_id hvl){
  std::cerr << "remove_constraint" << std::endl;
  constraint_set.erase(hvl);
  
  // We have to look at all subconstraints
  for(H_vertex_it it = hvl->begin(), succ = it; 
      ++succ != hvl->end(); 
      ++it){
    typename H_sc_to_c_map::iterator scit = sc_to_c_map.find(make_edge(*it,*succ));
    CGAL_triangulation_assertion(scit != sc_to_c_map.end());
    H_context_list* hcl = scit->second;

    // and remove the context of the constraint
    for(H_context_iterator ctit=hcl->begin(); ctit != hcl->end(); ctit++) {
      if(ctit->enclosing == hvl){
	    hcl->erase(ctit);
		break;
      }
    }
    // If the constraint passes several times through the same subconstraint,
    // the above loop maybe removes them in the wrong order

    // If this was the only context in the list, delete the context list
    if(hcl->empty()){
      sc_to_c_map.erase(scit);
      delete hcl;
    }
  }
  delete hvl;
}


// This function removes vertex v from the polyline constraint
// It only works for one polyline passing through v
// and for the case that the constrained edge u,w has no intersections
template <class T, class Data>
void Polyline_constraint_hierarchy_2<T,Data>::simplify(H_vertex_it uc,
                                              H_vertex_it vc,
                                              H_vertex_it wc)

{
  CGAL_assertion(vc->fixed != true);
  CGAL_assertion(vc->removed != true);
  Vertex_handle u = uc->vertex, v = vc->vertex, w = wc->vertex;
  typename H_sc_to_c_map::iterator uv_sc_iter = sc_to_c_map.find(make_edge(u, v));
  CGAL_assertion_msg( uv_sc_iter != sc_to_c_map.end(), "not a subconstraint" );
  H_context_list*  uv_hcl = uv_sc_iter->second;
  CGAL_assertion_msg(uv_hcl->size() == 1, "more than one constraint passing through the subconstraint" );
  if((uv_hcl->front().current())->vertex != u){
    std::swap(u,w);
    uv_sc_iter = sc_to_c_map.find(make_edge(u, v));
    CGAL_assertion_msg( uv_sc_iter != sc_to_c_map.end(), "not a subconstraint" );
    uv_hcl = (*uv_sc_iter).second;
    CGAL_assertion_msg(uv_hcl->size() == 1, "more than one constraint passing through the subconstraint" );
  }
  // now u,v, and w are ordered along the polyline constraint
  typename H_sc_to_c_map::iterator vw_sc_iter = sc_to_c_map.find(make_edge(v, w));
  CGAL_assertion_msg( vw_sc_iter != sc_to_c_map.end(), "not a subconstraint" );
  H_context_list*  vw_hcl = vw_sc_iter->second;
  CGAL_assertion_msg(vw_hcl->size() == 1, "more than one constraint passing through the subconstraint" );
  H_vertex_list* vertex_list = uv_hcl->front().id();
  CGAL_assertion_msg(vertex_list  == vw_hcl->front().id(), "subconstraints from different polyline constraints" );
  // Remove the list item which points to v
  vc->removed = true;
  // Remove the entries for [u,v] and [v,w]
  sc_to_c_map.erase(uv_sc_iter);
  sc_to_c_map.erase(vw_sc_iter);
  delete vw_hcl;
  // reuse other context list
  sc_to_c_map[make_edge(u,w)] = uv_hcl;
}


template <class T, class Data>
int
Polyline_constraint_hierarchy_2<T,Data>::remove_points_from_constraint(Constraint_id cid)
{
  int n=0;
  H_vertex_it b = cid->begin(), e = cid->end();
  while(b!=e){
    if(b->removed){
      H_vertex_it r = b;
      ++b;
      cid->erase(r);
      ++n;
    }else{
      ++b;
    }
  }
  return n;
}


template <class T, class Data>
int
Polyline_constraint_hierarchy_2<T,Data>::remove_points_from_constraints()
{
  int n = 0;
  for(H_c_iterator it = constraint_set.begin(); it!= constraint_set.end(); ++it){
    n+= remove_points_from_constraint(*it);
  }
  std::cerr << "Removed " << n << " points" << std::endl;
  return n;
}


template <class T, class Data>
typename Polyline_constraint_hierarchy_2<T,Data>::Constraint_id
Polyline_constraint_hierarchy_2<T,Data>::concatenate(Constraint_id first, Constraint_id second)
{
  std::cerr << "concatenate" << std::endl;
  constraint_set.erase(first);
  constraint_set.erase(second);
  // We have to look at all subconstraints
  for(H_vertex_it it = second->begin(), succ = it; 
      ++succ != second->end(); 
      ++it){
    typename H_sc_to_c_map::iterator scit = sc_to_c_map.find(make_edge(it->vertex,succ->vertex));
    CGAL_triangulation_assertion(scit != sc_to_c_map.end());
    H_context_list* hcl = scit->second;

    // and replace the context of the constraint
    for(H_context_iterator ctit=hcl->begin(); ctit != hcl->end(); ctit++) {
      if(ctit->enclosing == second){
	ctit->enclosing = first;
	break;
      }
    }
  }
  // now we really concatenate the vertex lists
  // Note that all iterators pointing into second remain valid.
  // This concerns user code, as well as  the data member "pos" of the Context class
  first->pop_back(); // because it is the same as second.front()
  H_vertex_it back_it = first->end();
  --back_it;
  first->splice(first->end(), *second);

  // Note that for VC8 with iterator debugging the iterators pointing into second
  // are NOT valid      So we have to update them
    for(H_vertex_it it = back_it, succ = it; 
      ++succ != first->end(); 
      ++it){
    typename H_sc_to_c_map::iterator scit = sc_to_c_map.find(make_edge(it->vertex,succ->vertex));
    CGAL_triangulation_assertion(scit != sc_to_c_map.end());
    H_context_list* hcl = scit->second;

    // and update pos in the context of the constraint
    for(H_context_iterator ctit=hcl->begin(); ctit != hcl->end(); ctit++) {
      if(ctit->enclosing == first){
	ctit->pos = it;
	break;
      }
    }
    }
  constraint_set.insert(first);

  delete second;
  return first;
}

template <class T, class Data>
typename Polyline_constraint_hierarchy_2<T,Data>::Constraint_id
Polyline_constraint_hierarchy_2<T,Data>::concatenate2(Constraint_id first, Constraint_id second)
{  
  std::cerr << "concatenate2" << std::endl;
  constraint_set.erase(first);
  constraint_set.erase(second);
  // We have to look at all subconstraints
  for(H_vertex_it it = first->begin(), succ = it; 
      ++succ != first->end(); 
      ++it){
    typename H_sc_to_c_map::iterator scit = sc_to_c_map.find(make_edge(it->vertex,succ->vertex));
    CGAL_triangulation_assertion(scit != sc_to_c_map.end());
    H_context_list* hcl = scit->second;

    // and replace the context of the constraint
    for(H_context_iterator ctit=hcl->begin(); ctit != hcl->end(); ctit++) {
      if(ctit->enclosing == first){
	ctit->enclosing = second;
	break;
      }
    }
  }
  // now we really concatenate the vertex lists
  // Note that all iterators pointing into second remain valid.
  first->pop_back(); // because it is the same as second.front()
  H_vertex_it back_it = first->end();
  --back_it;
  second->splice(second->begin(), *first);

  // Note that for VC8 with iterator debugging the iterators pointing into second
  // are NOT valid      So we have to update them
    for(H_vertex_it it = back_it, succ = it; 
      ++succ != first->end(); 
      ++it){
    typename H_sc_to_c_map::iterator scit = sc_to_c_map.find(make_edge(it->vertex,succ->vertex));
    CGAL_triangulation_assertion(scit != sc_to_c_map.end());
    H_context_list* hcl = scit->second;

    // and update pos in the context of the constraint
    for(H_context_iterator ctit=hcl->begin(); ctit != hcl->end(); ctit++) {
      if(ctit->enclosing == second){
	ctit->pos = it;
	break;
      }
    }
  }
  constraint_set.insert(second);

  delete first;
  return second;
}


  // split a constraint in two constraints, so that vcit becomes the first
  // vertex of the new constraint
  // returns the new constraint 
template <class T, class Data>
typename Polyline_constraint_hierarchy_2<T,Data>::Constraint_id
Polyline_constraint_hierarchy_2<T,Data>::split(Constraint_id first, H_vertex_it vcit)
{
  std::cerr << "split" << std::endl;
  constraint_set.erase(first);
  H_vertex_list* second = new H_vertex_list;
  second->splice(second->end(), *first, vcit, first->end());
  first->push_back(second->front()); // Duplicate the common vertex
  constraint_set.insert(first);
  constraint_set.insert(second);
 // We have to look at all subconstraints
  for(H_vertex_it it = second->begin(), succ = it; 
      ++succ != second->end(); 
      ++it){
    typename H_sc_to_c_map::iterator scit = sc_to_c_map.find(make_edge(it->vertex,succ->vertex));
    CGAL_triangulation_assertion(scit != sc_to_c_map.end());
    H_context_list* hcl = scit->second;

    // and replace the context of the constraint
    for(H_context_iterator ctit=hcl->begin(); ctit != hcl->end(); ctit++) {
      if(ctit->enclosing == first){
	ctit->enclosing = second;
	break;
      }
    }
  }
  return second;
}

template <class T, class Data>
typename Polyline_constraint_hierarchy_2<T,Data>::Constraint_id
Polyline_constraint_hierarchy_2<T,Data>::split2(Constraint_id first, H_vertex_it vcit)
{
  std::cerr << "split2" << std::endl;
  constraint_set.erase(first);
  H_vertex_list* second = new H_vertex_list;
  second->splice(second->end(), *first, first->begin(), vcit);
  second->push_back(first->front()); // Duplicate the common vertex
  constraint_set.insert(first);
  constraint_set.insert(second);
 // We have to look at all subconstraints
  for(H_vertex_it it = second->begin(), succ = it; 
      ++succ != second->end(); 
      ++it){
    typename H_sc_to_c_map::iterator scit = sc_to_c_map.find(make_edge(it->vertex,succ->vertex));
    CGAL_triangulation_assertion(scit != sc_to_c_map.end());
    H_context_list* hcl = scit->second;

    // and replace the context of the constraint
    for(H_context_iterator ctit=hcl->begin(); ctit != hcl->end(); ctit++) {
      if(ctit->enclosing == first){
	ctit->enclosing = second;
	break;
      }
    }
  }
  return second;
}

/*
when a constraint is inserted,
it is, at first, both  a constraint and a subconstraint
 */
template <class T, class Data>
typename Polyline_constraint_hierarchy_2<T,Data>::H_vertex_list*
Polyline_constraint_hierarchy_2<T,Data>::
insert_constraint(T va, T vb){
  H_edge        he = make_edge(va, vb);
  H_vertex_list*  children = new H_vertex_list; 
  H_context_list* fathers;

  typename H_sc_to_c_map::iterator scit = sc_to_c_map.find(he);
  if(scit == sc_to_c_map.end()){
    fathers = new H_context_list;
    sc_to_c_map.insert(std::make_pair(he,fathers));
  } else {
    fathers = scit->second;
  }

  children->push_front(Node(va));  // was he.first
  children->push_back(Node(vb));   // was he.second
  constraint_set.insert(children);
  H_context ctxt;
  ctxt.enclosing = children;
  ctxt.pos     = children->skip_begin();
  fathers->push_front(ctxt);

  return children;
}


template <class T, class Data>
void
Polyline_constraint_hierarchy_2<T,Data>::
append_constraint(H_vertex_list* vl, T va, T vb){
  H_edge        he = make_edge(va, vb);
  H_context_list* fathers;

  typename H_sc_to_c_map::iterator scit = sc_to_c_map.find(he);
  if(scit == sc_to_c_map.end()){
    fathers = new H_context_list;
    sc_to_c_map.insert(std::make_pair(he,fathers));
  } else {
    fathers = scit->second;
  }

  typename H_vertex_list::skip_iterator bit = vl->skip_end();
  --bit;
  vl->push_back(Node(vb));
  H_context ctxt;
  ctxt.enclosing = vl;
  ctxt.pos     = bit;
  fathers->push_front(ctxt);
}


template <class T, class Data>
void Polyline_constraint_hierarchy_2<T,Data>::
clear()
{
    std::cerr << "clear" << std::endl;
  H_c_iterator cit;
  H_sc_iterator scit;
  // clean and delete vertices lists
  for(cit=constraint_set.begin(); cit != constraint_set.end();  cit++){
    (*cit)->clear();
    delete (*cit);
  }
  // clean and delete context lists
  for(scit=sc_to_c_map.begin(); scit != sc_to_c_map.end(); scit++){
    (*scit).second->clear();
    delete (*scit).second;
  }
  sc_to_c_map.clear();
  constraint_set.clear();
}


template <class T, class Data>
bool Polyline_constraint_hierarchy_2<T,Data>::
next_along_sc(T va, T vb, T& w) const
{
  // find the next vertex after vb along any enclosing constrained
  // return false if there is no ....
  H_context_iterator  ctxtit, past;
  if(!get_contexts(va, vb, ctxtit, past)) CGAL_triangulation_assertion(false);

  H_vertex_it pos;
  for( ; ctxtit != past; ctxtit++){
    pos = ctxtit->pos;
    if((*pos)==va) {
      ++pos; ++pos;
      if (pos != ctxtit->enclosing->end()) {  w=(*pos); return true;}
    }
    else {
      if (pos != ctxtit->enclosing->begin()) {--pos; w=(*pos); return true;}
    }
  }
  return false;
}



/*
  Attention, le point v DOIT etre un point de Steiner,
  et les segments va,v et v,vb sont des sous contraintes.
*/
template <class T, class Data>
void Polyline_constraint_hierarchy_2<T,Data>::
remove_Steiner(T v, T va, T vb)
{
  // remove a Steiner point
  CGAL_precondition(!is_constrained_vertex(v));
 
  H_context_list*  hcl1;
  H_context_list*  hcl2;
  if(!get_contexts(va,v,hcl1)) CGAL_triangulation_assertion(false);
  if(!get_contexts(v,vb,hcl2)) CGAL_triangulation_assertion(false);

  H_vertex_it      pos;
  for(H_context_iterator ctit=hcl1->begin(); ctit != hcl1->end(); ctit++){
    pos = ctit->pos;
    if((*pos)==va) pos++;
    pos = ctit->enclosing->erase(pos);
    ctit->pos = --pos;
  }

  sc_to_c_map.erase(make_edge(va,v));
  sc_to_c_map.erase(make_edge(v,vb));
  delete hcl2;
  sc_to_c_map.insert(std::make_pair(make_edge(va,vb),hcl1));
}



/*
  same as add_Steiner
  precondition : va,vb est une souscontrainte. 
*/
template <class T, class Data>
void Polyline_constraint_hierarchy_2<T,Data>::
split_constraint(T va, T vb, T vc){
  add_Steiner(va, vb, vc);
}


template <class T, class Data>
void 
Polyline_constraint_hierarchy_2<T,Data>::
add_Steiner(T va, T vb, T vc){
  H_context_list* hcl;
  if(!get_contexts(va,vb,hcl)) CGAL_triangulation_assertion(false);

  H_context_list* hcl2 = new  H_context_list;

  H_vertex_it      pos;
  H_context  ctxt;
  for(H_context_iterator ctit=hcl->begin(); ctit != hcl->end(); ctit++) {
    // insert vc in enclosing constraint
    pos = ctit->current();
    ++pos;
    pos = ctit->enclosing->insert(pos, Node(vc));// fixed == true 
    vc->fixed = true;
    --pos;
    
    // set ctxt to the context of (vc,vb)
    // change *ctit in hcl to the context of (va,vc)
    // add ctxt to hcl2 list
    ctxt.enclosing = ctit->enclosing;  
    if((pos->vertex)==va)    {
      ctit->pos = pos;
      ctxt.pos = ++pos;
    }
    else { //(*pos)==vb
      ctxt.pos = pos;
      ctit->pos = ++pos;
    }
    hcl2->push_back(ctxt);
  }

  H_context_list* hcl3;
  if (get_contexts(va,vc,hcl3)) { // (va,vc) is already a subconstraint
    hcl3->splice(hcl3->end(), *hcl);
    delete hcl;
  }
  else   sc_to_c_map.insert(std::make_pair(make_edge(va,vc), hcl));

  if (get_contexts(vc,vb,hcl3)) {// (vc,vb) is already a subconstraint
    hcl3->splice(hcl3->end(),*hcl2);

    delete hcl2;
  }
  else  sc_to_c_map.insert(std::make_pair(make_edge(vc,vb), hcl2));
    
  
  sc_to_c_map.erase(make_edge(va,vb));
  return;
}


template <class T, class Data>
inline
typename Polyline_constraint_hierarchy_2<T,Data>::H_edge
Polyline_constraint_hierarchy_2<T,Data>::
make_edge(T va, T vb) const
{
  return (va<vb) ? H_edge(va,vb) : H_edge(vb,va);
}

template <class T, class Data>
inline
bool
Polyline_constraint_hierarchy_2<T,Data>::
get_contexts(T va, T vb, H_context_list* & hcl) const
{
  H_sc_iterator sc_iter = sc_to_c_map.find(make_edge(va,vb));
  if( sc_iter == sc_to_c_map.end() )    return(false);
  hcl = (*sc_iter).second;
  return true;
}

template <class T, class Data>
inline
bool
Polyline_constraint_hierarchy_2<T,Data>::
get_contexts(T va, T vb, 
	     H_context_iterator& ctxt, 
	     H_context_iterator& past) const
{
  H_context_list* hcl;
  if (!get_contexts(va,vb,hcl)) return false;
  ctxt = hcl->begin();
  past = hcl->end();
  return true;    
}



template <class T, class Data>
inline
typename Polyline_constraint_hierarchy_2<T,Data>::H_vertex_it
Polyline_constraint_hierarchy_2<T,Data>::
get_pos(T va, T vb) const
  //return pos in the first context
{
    return (*sc_to_c_map.find(make_edge(va,vb))).second->begin().pos;
}

template <class T, class Data>
void
Polyline_constraint_hierarchy_2<T,Data>::
oriented_end(T va, T vb, T& vc) const
{
  H_context_iterator ctxt, past;
  if(!get_contexts(va,vb, ctxt, past) ) CGAL_triangulation_assertion(false);
  if(*(ctxt->pos) == va)
    vc = ctxt->enclosing->back();
  else
    vc = ctxt->enclosing->front();
}


template <class T, class Data>
void
Polyline_constraint_hierarchy_2<T,Data>::
print() const
{ 
  H_c_iterator hcit;
  std::map<T,int>  vertex_num;
  int num = 0;
  for(hcit = c_begin(); hcit != c_end();  hcit++) {
    H_vertex_it vit = (*hcit)->begin();
    for (; vit != (*hcit)->end(); vit++){
      num ++;
      vertex_num.insert(std::make_pair((*vit), num));
    }
  }
//  typename std::map<T,int>::iterator vnit = vertex_num.begin();
//  for(; vnit != vertex_num.end(); vnit++) {
//    vnit->second = ++num;
//    std::cerr << "vertex num " << num  << " " << vnit->first->point()
//	      << std::endl;
//  }

  H_c_iterator cit=c_begin();
  H_sc_iterator scit=sc_begin();

  for(; cit != c_end();  cit++){
    std::cout << std::endl ;
    std::cout << "constraint " ;
    std::cout << (unsigned int)*cit;
    std::cout << "  subconstraints " ;
    H_vertex_it vit = (*cit)->begin();
    for(; vit != (*cit)->end(); vit++){
      std::cout << vertex_num[*vit]  <<" ";
    }
    vit = (*cit)->begin();
    for(; vit != (*cit)->end(); vit++){
      std::cout << (*vit)->point()  <<" ";
    }
  }
  std::cout << std::endl ;
  for(;scit != sc_end(); scit++){
    std::cout << "subconstraint " ;
    std::cout << vertex_num[scit->first.first] << " " 
	      << vertex_num[scit->first.second];
    H_context_iterator cb, ce;
    get_contexts(scit->first.first, scit->first.second, cb, ce);
    
    std::cout << "  enclosing " ;
    for(; cb != ce; cb++) { 
      std::cout << (unsigned int) cb->id();
      std::cout <<  "   " ;
    }
    std::cout << std::endl ;
  }
  return; 
}


} //namespace CGAL
#endif // CGAL_POLYLINE_CONSTRAINT_HIERARCHY_2_H
