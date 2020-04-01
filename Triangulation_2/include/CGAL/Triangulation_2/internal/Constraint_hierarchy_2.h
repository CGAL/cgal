// Copyright (c) 1999  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Olivier Billet, Mariette Yvinec

#ifndef CGAL_CONSTRAINT_HIERARCHY_2_H
#define CGAL_CONSTRAINT_HIERARCHY_2_H

#include <CGAL/license/Triangulation_2.h>


#include <CGAL/basic.h>
#include <utility>
#include <map>
#include <list>
#include <CGAL/triangulation_assertions.h>

namespace CGAL {

// T               is expected to be Vertex_handle
// Compare         is a comparison operator for type T
// Data            is intended to store info on a Vertex
template <class T, class Compare, class Data>
class Constraint_hierarchy_2
{
public:
  typedef std::pair<T, T>                      H_edge;
  typedef T                                    H_vertex;
  typedef Constraint_hierarchy_2<T,
                                 Compare,
                                 Data>         Hierarchy;
  typedef std::pair<T, T>                      H_constraint;
  typedef std::list<T>                         H_vertex_list;
  typedef std::list<H_constraint>              H_constraint_list;
  typedef typename std::list<T>::iterator               H_vertex_it;
  typedef typename std::list<H_constraint>::iterator    H_constraint_it;

  class Pair_compare {
    Compare comp;

  public:
    Pair_compare(const Compare& comp) : comp(comp) {}

    bool operator()(const H_edge& e1, const H_edge& e2) const {
      if(comp(e1.first, e2.first)) {
        return true;
      } else if((! comp(e2.first, e1.first)) && //  !less(e1,e2) && !less(e2,e1) == equal
                comp(e1.second, e2.second)) {
        return true;
      } else {
        return false;
      }
    }
  };

  class H_context {
    friend class Constraint_hierarchy_2<T,Compare,Data>;
  private:
    H_vertex_list*    enclosing;
    H_vertex_it       pos;
  public:
    H_vertex_it    vertices_begin() { return enclosing->begin();}
    H_vertex_it    current() {return pos;}
    H_vertex_it    vertices_end() {return enclosing->end();}
    std::size_t number_of_vertices() {return enclosing->size();}
  };
  typedef std::list<H_context>                 H_context_list;
  typedef typename std::list<H_context>::iterator       H_context_iterator;

  typedef std::map<T, Data>                             H_vertex_map;
  typedef std::map<H_constraint, H_vertex_list*,
                   Pair_compare>                        H_c_to_sc_map;
  typedef std::map<H_edge,   H_context_list*,
                   Pair_compare >                       H_sc_to_c_map;

  typedef typename H_c_to_sc_map::const_iterator        H_c_iterator;
  typedef typename H_sc_to_c_map::const_iterator        H_sc_iterator;
  typedef typename H_vertex_map::const_iterator         H_v_iterator;
  typedef std::pair<H_constraint, H_vertex_list*>       H_c_value;
  typedef std::pair<H_edge,   H_context_list*>          H_sc_value;

private:
  Compare comp;
  // data for the 1d hierarchy
  H_c_to_sc_map   c_to_sc_map;
  H_sc_to_c_map   sc_to_c_map;
  // data for the 0d hierarchy
  H_vertex_map    vertex_map;

public:
  Constraint_hierarchy_2(const Compare& comp_ = Compare())
    : comp(comp_)
    , c_to_sc_map(Pair_compare(comp))
    , sc_to_c_map(Pair_compare(comp))
  { }
  Constraint_hierarchy_2(const Constraint_hierarchy_2& ch);
  ~Constraint_hierarchy_2(){ clear();}
  void clear();
  Constraint_hierarchy_2& operator=(const Constraint_hierarchy_2& ch);

  // Query
  bool is_subconstrained_edge(T va, T vb) const;
  bool is_constrained_edge(T va, T vb) const;
  bool is_constrained_vertex(T v) const;
  bool vertices_in_constraint(H_constraint hc,
                              H_vertex_it& v_first,
                              H_vertex_it& v_past) const;
  H_vertex_it vertices_in_constraint_begin(T va, T vb);
  H_vertex_it vertices_in_constraint_end(T va, T vb);

  bool enclosing_constraint(H_edge he, H_constraint& hc) const;
  bool enclosing_constraint(T  vaa, T  vbb, T& va, T& vb) const;
  bool enclosing_constraints(T vaa, T  vbb,  H_constraint_list& hcl) const;
  bool next_along_sc(T va, T vb, T& w) const;
  void oriented_end(T va, T vb, T& vc) const;

  H_context context(T va, T vb);
  std::size_t number_of_enclosing_constraints(T va, T vb);
  H_context_iterator contexts_begin(T va, T vb);
  H_context_iterator contexts_end(T va, T vb);
  std::size_t number_of_constraints() { return c_to_sc_map.size();}
  std::size_t number_of_subconstraints() {return sc_to_c_map.size();}


  // insert/remove
  void add_Steiner(T va, T vb, T vx);
  bool insert_constraint(T va, T vb);
  void remove_constraint(T va, T vb);
  void split_constraint(T va, T vb, T vc);

  void constrain_vertex(T v, Data data=Data());
  void unconstrain_vertex(T v);
  void set_data(T v, Data data);
  Data get_data(T v);

  void remove_Steiner(T v, T va, T vb);

  // iterators
  H_sc_iterator sc_begin() const{ return sc_to_c_map.begin(); }
  H_sc_iterator sc_end()   const{ return sc_to_c_map.end();   }
  H_c_iterator  c_begin()  const{ return c_to_sc_map.begin(); }
  H_c_iterator  c_end()    const{ return c_to_sc_map.end();   }
  H_v_iterator  v_begin()  const{ return vertex_map.begin(); }
  H_v_iterator  v_end()    const{ return vertex_map.end(); }

  //Helping functions
  void copy(const Constraint_hierarchy_2& ch);
  void copy(const Constraint_hierarchy_2& ch, std::map<T,T>& vmap);
  void swap(Constraint_hierarchy_2& ch);

private:
  H_edge    make_edge(T va, T vb) const;
  H_vertex_it     get_pos(T va, T vb) const;
  bool      get_contexts(T va, T vb,
                         H_context_iterator& ctxt,
                         H_context_iterator& past) const;

  H_context_list * get_contexts(T va, T vb) const;

  //to_debug
public:
  void   print() const;
};

template <class T, class Compare, class Data>
Constraint_hierarchy_2<T,Compare,Data>::
Constraint_hierarchy_2(const Constraint_hierarchy_2& ch)
  : comp(ch.comp)
  , c_to_sc_map(Pair_compare(comp))
  , sc_to_c_map(Pair_compare(comp))
{
  copy(ch);
}

template <class T, class Compare, class Data>
Constraint_hierarchy_2<T,Compare,Data>&
Constraint_hierarchy_2<T,Compare,Data>::
operator=(const Constraint_hierarchy_2& ch){
  copy(ch);
  return *this;
}

template <class T, class Compare, class Data>
void
Constraint_hierarchy_2<T,Compare,Data>::
copy(const Constraint_hierarchy_2& ch1)
{
  // create a identity transfer vertex map
  std::map<T,T>  vmap;
  H_c_iterator cit1 = ch1.c_begin();
  for( ; cit1 != ch1.c_end(); ++cit1) {
    H_vertex_it vit = cit1->second->begin();
    for( ; vit != cit1->second->end(); ++vit) {
      vmap[*vit] = *vit;
    }
  }
  copy(ch1, vmap);
}

template <class T, class Compare, class Data>
void
Constraint_hierarchy_2<T,Compare,Data>::
copy(const Constraint_hierarchy_2& ch1, std::map<T,T>& vmap)
  // copy with a tranfer vertex map
{
  clear();
  // copy c_to_sc_map
  H_c_iterator cit1 = ch1.c_begin();
  for( ; cit1 != ch1.c_end(); ++cit1) {
    H_vertex u2 = vmap[cit1->first.first];
    H_vertex v2 = vmap[cit1->first.second];
    H_vertex_list* hvl1 = cit1->second;
    H_vertex_list* hvl2 = new H_vertex_list;
    H_vertex_it vit = hvl1->begin();
    for( ; vit != hvl1->end(); ++vit) hvl2->push_back(vmap[*vit]);
    c_to_sc_map[make_edge(u2,v2)] = hvl2;
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
      H_vertex u2 = vmap[cit1->enclosing->front()];
      H_vertex v2 = vmap[cit1->enclosing->back()];
      H_context ctxt2;
      ctxt2.enclosing = c_to_sc_map[make_edge(u2,v2)];
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
  // copy of vertex_map
  H_v_iterator hvit1 = ch1.vertex_map.begin();
  for ( ; hvit1 != ch1.vertex_map.end(); ++hvit1){
    vertex_map[vmap[hvit1->first]] = hvit1->second;
  }
  return;
}

template <class T, class Compare, class Data>
void
Constraint_hierarchy_2<T,Compare,Data>::
swap(Constraint_hierarchy_2& ch)
{
  c_to_sc_map.swap(ch.c_to_sc_map);
  sc_to_c_map.swap(ch.sc_to_c_map);
  vertex_map.swap(ch.vertex_map);
}

template <class T, class Compare, class Data>
bool Constraint_hierarchy_2<T,Compare,Data>::
is_constrained_vertex(T v) const
{
  return( vertex_map.find(v) != vertex_map.end() );
}


template <class T, class Compare, class Data>
bool Constraint_hierarchy_2<T,Compare,Data>::
is_constrained_edge(T va, T vb) const
{
  return( c_to_sc_map.find(make_edge(va, vb)) != c_to_sc_map.end() );
}

template <class T, class Compare, class Data>
bool Constraint_hierarchy_2<T,Compare,Data>::
is_subconstrained_edge(T va, T vb) const
{
  return( sc_to_c_map.find(make_edge(va, vb)) != sc_to_c_map.end() );
}

#if 0
template <class T, class Compare, class Data>
bool Constraint_hierarchy_2<T,Compare,Data>::
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
#endif

template <class T, class Compare, class Data>
bool Constraint_hierarchy_2<T,Compare,Data>::
enclosing_constraint(H_edge he, H_constraint& hc) const
{
  H_context_iterator hcit, past;
  if ( !get_contexts(he.first,he.second, hcit ,past)) return false;
  hc = make_edge(hcit->enclosing->front(), hcit->enclosing->back());
  return true;
}



template <class T, class Compare, class Data>
bool Constraint_hierarchy_2<T,Compare,Data>::
enclosing_constraint(T  vaa, T  vbb, T& va, T& vb) const
{
  H_context_iterator hcit, past;
  if ( !get_contexts(vaa,vbb, hcit ,past)) return false;
  va = hcit->enclosing->front();
  vb = hcit->enclosing->back();
  return true;
}


template <class T, class Compare, class Data>
bool Constraint_hierarchy_2<T,Compare,Data>::
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

template <class T, class Compare, class Data>
typename Constraint_hierarchy_2<T,Compare,Data>::H_context
Constraint_hierarchy_2<T,Compare,Data>::
context(T va, T vb)
{
  H_context_iterator hcit, past;
  if(!get_contexts(va,vb, hcit ,past)) CGAL_triangulation_assertion(false);
  return *hcit;
}

template <class T, class Compare, class Data>
std::size_t
Constraint_hierarchy_2<T,Compare,Data>::
number_of_enclosing_constraints(T va, T vb)
{
  H_context_list* hcl = get_contexts(va, vb);
  CGAL_triangulation_assertion(hcl != nullptr);
  return hcl->size();
}

template <class T, class Compare, class Data>
typename Constraint_hierarchy_2<T,Compare,Data>::H_context_iterator
Constraint_hierarchy_2<T,Compare,Data>::
contexts_begin(T va, T vb)
{
   H_context_iterator first, last;
   if( !get_contexts(va,vb,first,last))  CGAL_triangulation_assertion(false);
   return first;
}

template <class T, class Compare, class Data>
typename Constraint_hierarchy_2<T,Compare,Data>::H_context_iterator
Constraint_hierarchy_2<T,Compare,Data>::
contexts_end(T va, T vb)
{
   H_context_iterator first, last;
   if( !get_contexts(va,vb,first,last))  CGAL_triangulation_assertion(false);
   return last;
}

template <class T, class Compare, class Data>
typename Constraint_hierarchy_2<T,Compare,Data>::H_vertex_it
Constraint_hierarchy_2<T,Compare,Data>::
vertices_in_constraint_begin(T va, T vb)
{
  H_c_iterator  cit = c_to_sc_map.find(make_edge(va,vb));
  CGAL_triangulation_assertion( cit != c_to_sc_map.end());
  return cit->second->begin();
}

template <class T, class Compare, class Data>
typename Constraint_hierarchy_2<T,Compare,Data>::H_vertex_it
Constraint_hierarchy_2<T,Compare,Data>::
vertices_in_constraint_end(T va, T vb)
{
  H_c_iterator  cit = c_to_sc_map.find(make_edge(va,vb));
  CGAL_triangulation_assertion( cit != c_to_sc_map.end());
  return cit->second->end();
}


/*
when a constraint is inserted,
it is, at first, both  a constraint and a subconstraint
 */
template <class T, class Compare, class Data>
bool Constraint_hierarchy_2<T,Compare,Data>::
insert_constraint(T va, T vb){
  H_edge        he = make_edge(va, vb);
  H_vertex_list*  children = new H_vertex_list;

  children->push_front(he.first);
  children->push_back(he.second);
  bool insert_ok = (c_to_sc_map.insert(std::make_pair(he,children))).second;
  if (insert_ok) {
    H_context ctxt;
    ctxt.enclosing = children;
    ctxt.pos     = children->begin();

    // It may happen that the sub-constraint he is already in the map.  The
    // following variable 'fathers' is a reference to a pointer. If the
    // sub-constraint he is already in the map, 'fathers' will be the
    // pointer to its fathers list. If the sub-constraint he is new,
    // std::make_pair(he, 0) is inserted in the map (that is what does
    // map::operator[]), 'fathers' will be a default-constructed pointer
    // (that is the nullptr pointer), and it will re-assigned to a newly
    // created context list.  As 'father' is a reference (to a pointer),
    // there is no need to modify the map after that.
    H_context_list*& fathers = sc_to_c_map[he];
    if(fathers == 0) {
      fathers = new H_context_list;
    }
    fathers->push_front(ctxt);
    return true;
  }
  delete children;
  return false; //duplicate constraint - no insertion
}

template <class T, class Compare, class Data>
void
Constraint_hierarchy_2<T,Compare,Data>::
remove_constraint(T va, T vb){
  H_edge   he = make_edge(va, vb);
  typename H_c_to_sc_map::iterator c_to_sc_it = c_to_sc_map.find(he);
  CGAL_triangulation_assertion(c_to_sc_it != c_to_sc_map.end());

  H_vertex_list *hvl = c_to_sc_it->second;

  // We have to look at all subconstraints
  for(H_vertex_it it = hvl->begin(), succ = it;
      ++succ != hvl->end();
      ++it){
    typename H_sc_to_c_map::iterator scit =
                                   sc_to_c_map.find(make_edge(*it,*succ));
    CGAL_triangulation_assertion(scit != sc_to_c_map.end());
    H_context_list* hcl = scit->second;

    // and remove the contraint from the context list of the subcontraint
    for(H_context_iterator ctit=hcl->begin(); ctit != hcl->end(); ctit++) {
      if(ctit->enclosing == hvl){
            hcl->erase(ctit);
                break;
      }
    }
    // If this was the only context in the list, delete the context list
    if(hcl->empty()){
      sc_to_c_map.erase(scit);
      delete hcl;
    }
  }
  c_to_sc_map.erase(c_to_sc_it);
  delete hvl;
}


template <class T, class Compare, class Data>
void Constraint_hierarchy_2<T,Compare,Data>::
constrain_vertex(T v, Data data){
  vertex_map.insert(std::make_pair(v,data));
}


template <class T, class Compare, class Data>
void Constraint_hierarchy_2<T,Compare,Data>::
unconstrain_vertex(T v){
  vertex_map.erase(v);
}


template <class T, class Compare, class Data>
Data Constraint_hierarchy_2<T,Compare,Data>::
get_data(T v){
  CGAL_precondition( is_constrained_vertex(v) );
  return (*vertex_map.find(v)).second;
}


template <class T, class Compare, class Data>
void Constraint_hierarchy_2<T,Compare,Data>::
set_data(T v, Data data){
  vertex_map.erase(v);
  vertex_map.insert(std::make_pair(v,data));
}


template <class T, class Compare, class Data>
void Constraint_hierarchy_2<T,Compare,Data>::
clear()
{
  H_c_iterator cit;
  H_sc_iterator scit;
  // clean and delete vertices lists
  for(cit=c_to_sc_map.begin(); cit != c_to_sc_map.end();  cit++){
    (*cit).second->clear();
    delete (*cit).second;
  }
  // clean and delete context lists
  for(scit=sc_to_c_map.begin(); scit != sc_to_c_map.end(); scit++){
    (*scit).second->clear();
    delete (*scit).second;
  }
  sc_to_c_map.clear();
  c_to_sc_map.clear();
  vertex_map.clear();
}


template <class T, class Compare, class Data>
bool Constraint_hierarchy_2<T,Compare,Data>::
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
template <class T, class Compare, class Data>
void Constraint_hierarchy_2<T,Compare,Data>::
remove_Steiner(T v, T va, T vb)
{
  // remove a Steiner point
  CGAL_precondition(!is_constrained_vertex(v));

  H_context_list*  hcl1 = get_contexts(va, v);
  CGAL_triangulation_assertion(hcl1 != nullptr);
  H_context_list*  hcl2 = get_contexts(v, vb);
  CGAL_triangulation_assertion(hcl2 != nullptr);

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
template <class T, class Compare, class Data>
void Constraint_hierarchy_2<T,Compare,Data>::
split_constraint(T va, T vb, T vc){
  add_Steiner(va, vb, vc);
}


template <class T, class Compare, class Data>
void
Constraint_hierarchy_2<T,Compare,Data>::
add_Steiner(T va, T vb, T vc){
  H_context_list* hcl = get_contexts(va, vb);
  CGAL_triangulation_assertion(hcl != nullptr);

  H_context_list* hcl2 = new  H_context_list;

  H_vertex_it      pos;
  H_context  ctxt;
  for(H_context_iterator ctit=hcl->begin(); ctit != hcl->end(); ctit++) {
    // insert vc in enclosing constraint
    pos = ctit->pos;
    ++pos;
    pos = ctit->enclosing->insert(pos, vc);
    --pos;

    // set ctxt to the context of (vc,vb)
    // change *ctit in hcl to the context of (va,vc)
    // add ctxt to hcl2 list
    ctxt.enclosing = ctit->enclosing;
    if((*pos)==va)    {
      ctit->pos = pos;
      ctxt.pos = ++pos;
    }
    else { //(*pos)==vb
      ctxt.pos = pos;
      ctit->pos = ++pos;
    }
    hcl2->push_back(ctxt);
  }

  if (H_context_list* hcl3 = get_contexts(va,vc)) { // (va,vc) is already a subconstraint
    hcl3->splice(hcl3->end(), *hcl);
    delete hcl;
  }
  else   sc_to_c_map.insert(std::make_pair(make_edge(va,vc), hcl));

  if (H_context_list* hcl3 = get_contexts(vc,vb)) {// (vc,vb) is already a subconstraint
    hcl3->splice(hcl3->end(),*hcl2);
    delete hcl2;
  }
  else  sc_to_c_map.insert(std::make_pair(make_edge(vc,vb), hcl2));


  sc_to_c_map.erase(make_edge(va,vb));
  return;
}


template <class T, class Compare, class Data>
inline
typename Constraint_hierarchy_2<T,Compare,Data>::H_edge
Constraint_hierarchy_2<T,Compare,Data>::
make_edge(T va, T vb) const
{
  return comp(va, vb) ? H_edge(va,vb) : H_edge(vb,va);
}

template <class T, class Compare, class Data>
inline
typename Constraint_hierarchy_2<T,Compare,Data>::H_context_list*
Constraint_hierarchy_2<T,Compare,Data>::
get_contexts(T va, T vb) const
{
  H_sc_iterator sc_iter = sc_to_c_map.find(make_edge(va,vb));
  if ( sc_iter == sc_to_c_map.end() )
    return nullptr;
  return (*sc_iter).second;
}

template <class T, class Compare, class Data>
inline
bool
Constraint_hierarchy_2<T,Compare,Data>::
get_contexts(T va, T vb,
             H_context_iterator& ctxt,
             H_context_iterator& past) const
{
  if (H_context_list * hcl = get_contexts(va, vb)) {
    ctxt = hcl->begin();
    past = hcl->end();
    return true;
  }
  return false;
}



template <class T, class Compare, class Data>
inline
typename Constraint_hierarchy_2<T,Compare,Data>::H_vertex_it
Constraint_hierarchy_2<T,Compare,Data>::
get_pos(T va, T vb) const
  //return pos in the first context
{
    return (*sc_to_c_map.find(make_edge(va,vb))).second->begin().pos;
}

template <class T, class Compare, class Data>
void
Constraint_hierarchy_2<T,Compare,Data>::
oriented_end(T va, T vb, T& vc) const
{
  H_context_iterator ctxt, past;
  if(!get_contexts(va,vb, ctxt, past) ) CGAL_triangulation_assertion(false);
  if(*(ctxt->pos) == va)
    vc = ctxt->enclosing->back();
  else
    vc = ctxt->enclosing->front();
}


template <class T, class Compare, class Data>
void
Constraint_hierarchy_2<T,Compare,Data>::
print() const
{
  H_c_iterator hcit;
  std::map<T,int>  vertex_num;
  int num = 0;
  for(hcit = c_begin(); hcit != c_end();  hcit++) {
    H_vertex_it vit = (*hcit).second->begin();
    for (; vit != (*hcit).second->end(); vit++){
      num ++;
      vertex_num.insert(std::make_pair((*vit), num));
    }
  }
//  typename std::map<T,int>::iterator vnit = vertex_num.begin();
//  for(; vnit != vertex_num.end(); vnit++) {
//    vnit->second = ++num;
//    std::cerr << "vertex num " << num  << " " << vnit->first->point()
//              << std::endl;
//  }

  H_c_iterator cit=c_begin();
  H_sc_iterator scit=sc_begin();

  for(; cit != c_end();  cit++){
    std::cout << std::endl ;
    std::cout << "constraint " ;
    std::cout << vertex_num[cit->first.first]  << " "
              <<  vertex_num[cit->first.second];
    std::cout << "  subconstraints " ;
    H_vertex_it vit = cit->second->begin();
    for(; vit != cit->second->end(); vit++){
      std::cout << vertex_num[*vit]  <<" ";
    }
  }
  std::cout << std::endl ;
  for(;scit != sc_end(); scit++){
    std::cout << "subconstraint " ;
    std::cout << vertex_num[scit->first.first] << " "
              << vertex_num[scit->first.second];
    H_constraint_list  hcl;
    enclosing_constraints( scit->first.first, scit->first.second,
                           hcl);
    H_constraint_it hcit=hcl.begin();
    std::cout << "  enclosing " ;
    for(; hcit != hcl.end(); hcit++) {
      std::cout << vertex_num[hcit->first] <<" "
                << vertex_num[hcit->second] ;
      std::cout <<  "   " ;
    }
    std::cout << std::endl ;
  }
  return;
}


} //namespace CGAL
#endif // CGAL_CONSTRAINT_HIERARCHY_2_H
