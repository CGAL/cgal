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

#include <CGAL/license/Triangulation_2.h>


#include <CGAL/basic.h>
#include <utility>
#include <map>
#include <set> 
#include <list> 
#include <CGAL/Skiplist.h>
#include <CGAL/triangulation_assertions.h>

namespace CGAL {

// T               is expected to be Vertex_handle
// Compare         is a comparison operator for type T
// Data            is intended to store info on a Vertex
template <class T, class Compare, class Data>
class Polyline_constraint_hierarchy_2
{
public:
  typedef Data                                    Point;
  typedef T                                       Vertex_handle;
  typedef std::pair<T, T>                         Edge;
  typedef std::pair<T, T>                         Constraint;
  typedef std::pair<T, T>                         Subconstraint;

private:
  class Node {
  public:
    explicit Node(Vertex_handle vh, bool input = false)
      : vertex_(vh), point_(vh->point()), id(-1), input(input)
    {}
    Point& point() { return point_; }
    const Point& point() const { return point_; }
    Vertex_handle vertex() const { return vertex_; }
  private:
    Vertex_handle vertex_;
    Point point_;
  public:
    int id;
    bool input;
  };

  typedef CGAL::Skiplist<Node>  Vertex_list;
  typedef std::list<Constraint> Constraint_list;

public:
  // the base line is always 
  class Point_it 
    : public boost::iterator_adaptor<
    Point_it
    , typename Vertex_list::all_iterator 
    , Point
    >
  {
  public:
    Point_it() : Vertex_it::iterator_adaptor_() {}
    Point_it(typename Vertex_list::all_iterator it) : Point_it::iterator_adaptor_(it) {}
  private:
    friend class boost::iterator_core_access;
    Point& dereference() const { return this->base()->point(); }
  };

  // only nodes with a vertex_handle that is still in the triangulation
  class Vertex_it 
    : public boost::iterator_adaptor<
    Vertex_it
    , typename Vertex_list::skip_iterator 
    , Vertex_handle
    , boost::use_default
    , Vertex_handle>
  {
  public:
    Vertex_it() : Vertex_it::iterator_adaptor_() {}
    Vertex_it(typename Vertex_list::skip_iterator it) : Vertex_it::iterator_adaptor_(it) {}
    operator Point_it() const { return Point_it(this->base()); }
    bool& input() { return this->base()->input; }
  private:
    friend class boost::iterator_core_access;
    Vertex_handle dereference() const { return this->base()->vertex(); }
  };

  typedef typename Constraint_list::iterator Constraint_it;

  struct Constraint_id {
    Vertex_list* second;

    Constraint_id(): second(NULL) {}

    Constraint_id(Vertex_list* vl)
      : second(vl)
    {}

    Vertex_list* vl_ptr() const {return second;}

    operator std::pair<std::pair<Vertex_handle, Vertex_handle>,Vertex_list*>()
    { 
      if (second!=NULL){
        return std::make_pair(std::make_pair(second->front().vertex(),
                                             second->back().vertex()),second);
      } 
      return std::make_pair(std::make_pair(Vertex_handle(),Vertex_handle()),second);
    }

    bool operator == (const Constraint_id& other) const
    {
      return second == other.second;
    }
    bool operator != (const Constraint_id& other) const
    {
      return second != other.second;
    }

    bool operator<(const Constraint_id& other) const{
      return second < other.second;
    }
  };

  class Pair_compare {
    Compare comp;

  public:
    Pair_compare(const Compare& comp) : comp(comp) {}

    bool operator()(const Edge& e1, const Edge& e2) const {
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

  class Context {
    friend class Polyline_constraint_hierarchy_2<T,Compare,Data>;
  private:
    Vertex_list*    enclosing;
    Vertex_it       pos;
  public:
    Context() : enclosing(NULL) {}

    Context(const Context& hc)
      : enclosing(hc.enclosing), pos(hc.pos)
    {}

    Vertex_it    vertices_begin() { return enclosing->skip_begin();}
    Vertex_it    current() {return pos;}
    Vertex_it    vertices_end() {return enclosing->skip_end();}
    Constraint_id  id() { return enclosing; }
    std::size_t    number_of_vertices() const {return enclosing->skip_size(); }
  };                                           

  typedef std::list<Context>              Context_list;
  typedef typename Context_list::iterator Context_iterator;

  typedef std::set<Constraint_id>           Constraint_set;
  typedef std::map<Edge, Context_list*,
		   Pair_compare>            Sc_to_c_map;
  typedef typename Constraint_set::iterator C_iterator;
  typedef typename Sc_to_c_map::const_iterator    Sc_iterator;
  typedef Sc_iterator Subconstraint_iterator;
  
private:
  // data for the 1d hierarchy
  Compare          comp;
  Constraint_set   constraint_set;
  Sc_to_c_map      sc_to_c_map;
  std::map<std::pair<Vertex_handle, Vertex_handle>,
	   Constraint_id,
	   Pair_compare> constraint_map;
  
public:
  Polyline_constraint_hierarchy_2(const Compare& comp)
    : comp(comp)
    , sc_to_c_map(Pair_compare(comp))
    , constraint_map(Pair_compare(comp))
  { }
  Polyline_constraint_hierarchy_2(const Polyline_constraint_hierarchy_2& ch); 
  ~Polyline_constraint_hierarchy_2(){ clear();}
  void clear();
  Polyline_constraint_hierarchy_2& operator=(const Polyline_constraint_hierarchy_2& ch);

  // Query 
  bool is_subconstrained_edge(T va, T vb) const;
  bool is_constrained_edge(T va, T vb) const;
  bool is_constrained_vertex(T v) const;

  Vertex_it vertices_in_constraint_begin(Constraint_id cid) const
  { return cid.vl_ptr()->skip_begin(); }
  Vertex_it vertices_in_constraint_end(Constraint_id cid) const
  { return cid.vl_ptr()->skip_end(); }

  Vertex_it vertices_in_constraint_begin(Vertex_handle va, Vertex_handle vb) const
  { Constraint_id cid = constraint_map.find(make_edge(va,vb))->second;
    return cid.vl_ptr()->skip_begin(); }

  Vertex_it vertices_in_constraint_end(Vertex_handle va, Vertex_handle vb) const
  { Constraint_id cid = constraint_map.find(make_edge(va,vb))->second;
    return cid.vl_ptr()->skip_end(); }

  Point_it points_in_constraint_begin(Constraint_id cid) const
  { return cid.vl_ptr()->all_begin(); }
  Point_it points_in_constraint_end(Constraint_id cid) const
  { return cid.vl_ptr()->all_end(); }

  bool enclosing_constraint(Edge he, Constraint& hc) const;
  bool enclosing_constraint(T  vaa, T  vbb, T& va, T& vb) const;
  bool enclosing_constraints(T vaa, T  vbb,  Constraint_list& hcl) const;
  bool next_along_sc(T va, T vb, T& w) const;
  void oriented_end(T va, T vb, T& vc) const;

  Context context(T va, T vb);
  std::size_t number_of_enclosing_constraints(T va, T vb) const;
  Context_iterator contexts_begin(T va, T vb) const;
  Context_iterator contexts_end(T va, T vb) const;
  std::size_t number_of_constraints() const  { return constraint_set.size();}
  std::size_t number_of_subconstraints()const {return sc_to_c_map.size();}
  

  // insert/remove
  void add_Steiner(T va, T vb, T vx);
  Vertex_list* insert_constraint(T va, T vb);
  void append_constraint(Constraint_id cid, T va, T vb);
  void swap(Constraint_id first, Constraint_id second);
  void remove_constraint(Constraint_id cid);

  void remove_constraint(Vertex_handle va, Vertex_handle vb)
  {
    remove_constraint(constraint_map[make_edge(va,vb)]);
  }

  void split_constraint(T va, T vb, T vc);

  void simplify(Vertex_it u,
                Vertex_it v,
                Vertex_it w);

  std::size_t remove_points_without_corresponding_vertex(Constraint_id);
  std::size_t remove_points_without_corresponding_vertex();

  Constraint_id concatenate(Constraint_id first, Constraint_id second);
  Constraint_id concatenate2(Constraint_id first, Constraint_id second);
  Constraint_id split(Constraint_id first, Vertex_it vcit);
  Constraint_id split2(Constraint_id first, Vertex_it vcit);

  void remove_Steiner(T v, T va, T vb);

  // iterators

  Subconstraint_iterator subconstraint_begin() const
  { 
    return sc_to_c_map.begin(); 
  }

  Subconstraint_iterator subconstraint_end() const
  { 
    return sc_to_c_map.end();   
  }

  Sc_iterator sc_begin() const{ return sc_to_c_map.begin(); }
  Sc_iterator sc_end()   const{ return sc_to_c_map.end();   }
  C_iterator  c_begin()  const{ return constraint_set.begin(); }
  C_iterator  c_end()    const{ return constraint_set.end();   }
  
  // Helper functions
  void copy(const Polyline_constraint_hierarchy_2& ch);
  void copy(const Polyline_constraint_hierarchy_2& ch, std::map<Vertex_handle,Vertex_handle>& vmap);
  void swap(Polyline_constraint_hierarchy_2& ch);

private: 
  Edge      make_edge(T va, T vb) const;
  Vertex_it get_pos(T va, T vb) const;
  bool      get_contexts(T va, T vb, 
			 Context_iterator& ctxt, 
			 Context_iterator& past) const;

  bool      get_contexts(T va, T vb, Context_list*&) const;

  //to_debug
public:
  void   print() const;
};

template <class T, class Compare, class Data>
Polyline_constraint_hierarchy_2<T,Compare,Data>::
Polyline_constraint_hierarchy_2(const Polyline_constraint_hierarchy_2& ch)
  : comp(ch.comp)
  , sc_to_c_map(Pair_compare(comp))
  , constraint_map(Pair_compare(comp))
{
  copy(ch);
}

template <class T, class Compare, class Data>
Polyline_constraint_hierarchy_2<T,Compare,Data>&
Polyline_constraint_hierarchy_2<T,Compare,Data>::
operator=(const Polyline_constraint_hierarchy_2& ch){
  copy(ch);
  return *this;
}

template <class T, class Compare, class Data>
void
Polyline_constraint_hierarchy_2<T,Compare,Data>::
copy(const Polyline_constraint_hierarchy_2& ch1)
{
  // create a identity transfer vertex map
  std::map<Node, Node>  vmap;
  C_iterator cit1 = ch1.c_begin();
  for( ; cit1 != ch1.c_end(); ++cit1) {
    Vertex_it vit = cit1->second->begin();
    for( ; vit != cit1->second->end(); ++vit) {
      vmap[*vit] = *vit;
    }
  }
  copy(ch1, vmap);
}

template <class T, class Compare, class Data>
void
Polyline_constraint_hierarchy_2<T,Compare,Data>::
copy(const Polyline_constraint_hierarchy_2& ch1, std::map<Vertex_handle,Vertex_handle>& vmap)
  // copy with a transfer vertex map
{
  std::map<Vertex_list*,Vertex_list*> vlmap;
  clear();
  // copy constraint_set
  C_iterator cit1 = ch1.c_begin();
  for( ; cit1 != ch1.c_end(); ++cit1) {
    Vertex_list* hvl1 = cit1->vl_ptr();
    Vertex_list* hvl2 = new Vertex_list;
    vlmap[hvl1] = hvl2;
    Vertex_it vit = hvl1->skip_begin(), end = hvl1->skip_end();
    for( ; vit != end; ++vit) hvl2->push_back(Node(vmap[*vit]));
    constraint_set.insert(hvl2);
  }
  // copy sc_to_c_map
  Sc_iterator scit1 = ch1.sc_begin();
  for( ; scit1 != ch1.sc_end(); ++scit1) {
    //vertices of the subconstraints
    Vertex_handle uu2 = vmap[scit1->first.first];
    Vertex_handle vv2 = vmap[scit1->first.second];
    Context_list* hcl1  = scit1->second;
    Context_list* hcl2  = new Context_list;
    Context_iterator cit1 = hcl1->begin();
    for( ; cit1 != hcl1->end(); ++cit1){
      // vertices of the enclosing constraints
      Context ctxt2;
      ctxt2.enclosing = vlmap[cit1->enclosing];
      ctxt2.pos = ctxt2.enclosing->skip_begin();
      Vertex_it aux = cit1->enclosing->skip_begin();
      while( aux != cit1->pos) {
	++aux;
	++ctxt2.pos;
      }
      hcl2->push_back(ctxt2);
    }
    sc_to_c_map[make_edge(uu2,vv2)] = hcl2;
  }

  comp = ch1.comp;
  return;
}


template <class T, class Compare, class Data>
void
Polyline_constraint_hierarchy_2<T,Compare,Data>::
swap(Polyline_constraint_hierarchy_2& ch)
{
  constraint_set.swap(ch.constraint_set);
  sc_to_c_map.swap(ch.sc_to_c_map);
}


/*
template <class T, class Compare, class Data>
bool Polyline_constraint_hierarchy_2<T,Compare,Data>::
is_constrained_edge(T va, T vb) const
{
  return( c_to_sc_map.find(make_edge(va, vb)) != c_to_sc_map.end() );
}
*/

template <class T, class Compare, class Data>
bool Polyline_constraint_hierarchy_2<T,Compare,Data>::
is_subconstrained_edge(T va, T vb) const
{
  return( sc_to_c_map.find(make_edge(va, vb)) != sc_to_c_map.end() );
}


// af: obsolete
template <class T, class Compare, class Data>
bool Polyline_constraint_hierarchy_2<T,Compare,Data>::
enclosing_constraint(Edge he, Constraint& hc) const
{
  Context_iterator hcit, past;
  if ( !get_contexts(he.first,he.second, hcit ,past)) return false;
  hc = make_edge(hcit->enclosing->front(), hcit->enclosing->back());
  return true;
}


// used by Constrained_triangulation_plus_2::intersect with Exact_intersection_tag
template <class T, class Compare, class Data>
bool Polyline_constraint_hierarchy_2<T,Compare,Data>::
enclosing_constraint(T  vaa, T  vbb, T& va, T& vb) const
{
  Context_iterator hcit, past;
  if ( !get_contexts(vaa,vbb, hcit ,past)) return false;
  // va = hcit->enclosing->front().vertex();
  // vb = hcit->enclosing->back().vertex();
  // Vertex_list* vl = hcit->enclosing;
  Vertex_it pos = hcit->pos;
  if(vaa != *pos){
    std::swap(vaa,vbb);
  }
  while(!pos.input()){
    --pos;
  }
  va = *pos;
  pos = hcit->pos;
  ++pos;
  CGAL_triangulation_assertion(vbb == *pos);
  while(!pos.input()){
    ++pos;
  }
  vb = *pos;
  return true;
}

// af: obsolete
template <class T, class Compare, class Data>
bool Polyline_constraint_hierarchy_2<T,Compare,Data>::
enclosing_constraints(T vaa, T vbb , Constraint_list& hcl) const
{
  Context_iterator hcit, past;
  if ( !get_contexts(vaa,vbb, hcit ,past)) return false;
  for (; hcit!=past; hcit++) {
    hcl.push_back(make_edge(hcit->enclosing->front(), 
			    hcit->enclosing->back())); 
  }
  return true;
}

template <class T, class Compare, class Data>
typename Polyline_constraint_hierarchy_2<T,Compare,Data>::Context
Polyline_constraint_hierarchy_2<T,Compare,Data>::
context(T va, T vb)
{
  Context_iterator hcit, past;
  if(!get_contexts(va,vb, hcit ,past)) CGAL_triangulation_assertion(false);
  return *hcit;
}

template <class T, class Compare, class Data>
std::size_t 
Polyline_constraint_hierarchy_2<T,Compare,Data>::
number_of_enclosing_constraints(T va, T vb) const
{
  Context_list* hcl = NULL;
  CGAL_triangulation_assertion_code( bool found = ) get_contexts(va,vb,hcl);
  CGAL_triangulation_assertion(found);
  return hcl->size();
}

template <class T, class Compare, class Data>
typename Polyline_constraint_hierarchy_2<T,Compare,Data>::Context_iterator
Polyline_constraint_hierarchy_2<T,Compare,Data>::
contexts_begin(T va, T vb) const
{
   Context_iterator first, last;
   if( !get_contexts(va,vb,first,last))  CGAL_triangulation_assertion(false);
   return first;
}

template <class T, class Compare, class Data>
typename Polyline_constraint_hierarchy_2<T,Compare,Data>::Context_iterator
Polyline_constraint_hierarchy_2<T,Compare,Data>::
contexts_end(T va, T vb) const
{   
   Context_iterator first, last;
   if( !get_contexts(va,vb,first,last))  CGAL_triangulation_assertion(false);
   return last;
} 

template <class T, class Compare, class Data>
void
Polyline_constraint_hierarchy_2<T,Compare,Data>::
swap(Constraint_id first, Constraint_id second){
    // We have to look at all subconstraints
  for(Vertex_it it = first.vl_ptr()->skip_begin(), succ = it, end = first.vl_ptr()->skip_end(); 
      ++succ != end; 
      ++it){
    typename Sc_to_c_map::iterator scit = sc_to_c_map.find(make_edge(*it,*succ));
    CGAL_triangulation_assertion(scit != sc_to_c_map.end());
    Context_list* hcl = scit->second;

    // and replace the context of the constraint
    for(Context_iterator ctit=hcl->begin(); ctit != hcl->end(); ctit++) {
      if(ctit->enclosing == first.vl_ptr()){
	ctit->enclosing = 0;
	break;
      }
    }
  }
    // We have to look at all subconstraints
  for(Vertex_it it = second.vl_ptr()->skip_begin(), succ = it, end = second.vl_ptr()->skip_end(); 
      ++succ != end; 
      ++it){
    typename Sc_to_c_map::iterator scit = sc_to_c_map.find(make_edge(*it,*succ));
    CGAL_triangulation_assertion(scit != sc_to_c_map.end());
    Context_list* hcl = scit->second;

    // and replace the context of the constraint
    for(Context_iterator ctit=hcl->begin(); ctit != hcl->end(); ctit++) {
      if(ctit->enclosing == second.vl_ptr()){
	ctit->enclosing = first.vl_ptr();
	break;
      }
    }
  }   
  // We have to look at all subconstraints
  for(Vertex_it it = first.vl_ptr()->skip_begin(), succ = it, end = first.vl_ptr()->skip_end(); 
      ++succ != end; 
      ++it){
    typename Sc_to_c_map::iterator scit = sc_to_c_map.find(make_edge(*it,*succ));
    CGAL_triangulation_assertion(scit != sc_to_c_map.end());
    Context_list* hcl = scit->second;

    // and replace the context of the constraint
    for(Context_iterator ctit=hcl->begin(); ctit != hcl->end(); ctit++) {
      if(ctit->enclosing == 0){
	ctit->enclosing = second.vl_ptr();
	break;
      }
    }
  }
  first.vl_ptr()->swap(*second.vl_ptr());
}


template <class T, class Compare, class Data>
void
Polyline_constraint_hierarchy_2<T,Compare,Data>::
remove_constraint(Constraint_id cid){
  constraint_set.erase(cid);
  
  // We have to look at all subconstraints
  for(Vertex_it it = cid.vl_ptr()->skip_begin(), succ = it, end = cid.vl_ptr()->skip_end(); 
      ++succ != end; 
      ++it){
    typename Sc_to_c_map::iterator scit = sc_to_c_map.find(make_edge(*it,*succ));
    CGAL_triangulation_assertion(scit != sc_to_c_map.end());
    Context_list* hcl = scit->second;

    // and remove the context of the constraint
    for(Context_iterator ctit=hcl->begin(); ctit != hcl->end(); ctit++) {
      if(ctit->enclosing == cid.vl_ptr()){
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
  delete cid.vl_ptr();
}


// This function removes vertex v from the polyline constraint
// It only works for one polyline passing through v
// and for the case that the constrained edge u,w has no intersections
template <class T, class Compare, class Data>
void Polyline_constraint_hierarchy_2<T,Compare,Data>::simplify(Vertex_it uc,
                                                       Vertex_it vc,
                                                       Vertex_it wc)

{
  Vertex_handle u = *uc, v = *vc, w = *wc;
  typename Sc_to_c_map::iterator uv_sc_iter = sc_to_c_map.find(make_edge(u, v));
  CGAL_assertion_msg( uv_sc_iter != sc_to_c_map.end(), "not a subconstraint" );
  Context_list*  uv_hcl = uv_sc_iter->second;
  CGAL_assertion_msg((u == w) || (uv_hcl->size() == 1), "more than one constraint passing through the subconstraint" );

  if(*(uv_hcl->front().current()) != u) {
    std::swap(u,w);
    uv_sc_iter = sc_to_c_map.find(make_edge(u, v));
    CGAL_assertion_msg( uv_sc_iter != sc_to_c_map.end(), "not a subconstraint" );
    uv_hcl = (*uv_sc_iter).second;
    CGAL_assertion_msg((u == w) || (uv_hcl->size() == 1), "more than one constraint passing through the subconstraint" );
  }
  // now u,v, and w are ordered along the polyline constraint
  if(vc.input()){
    uc.input() = true;
    wc.input() = true;
  }
  typename Sc_to_c_map::iterator vw_sc_iter = sc_to_c_map.find(make_edge(v, w));
  CGAL_assertion_msg( vw_sc_iter != sc_to_c_map.end(), "not a subconstraint" );
  Context_list*  vw_hcl = vw_sc_iter->second;
    CGAL_assertion_msg((u == w) || (vw_hcl->size() == 1), "more than one constraint passing through the subconstraint" );
 
  Vertex_list* vertex_list = uv_hcl->front().id().vl_ptr();
  CGAL_assertion_msg(vertex_list  == vw_hcl->front().id().vl_ptr(), "subconstraints from different polyline constraints" );
  // Remove the list item which points to v
  vertex_list->skip(vc.base());
  
  if(u != w){
    // Remove the entries for [u,v] and [v,w]
    sc_to_c_map.erase(uv_sc_iter);
    sc_to_c_map.erase(vw_sc_iter);
    delete vw_hcl;
    // reuse other context list
    sc_to_c_map[make_edge(u,w)] = uv_hcl;
  }else{
    sc_to_c_map.erase(uv_sc_iter);
    delete vw_hcl;
  }
}


template <class T, class Compare, class Data>
std::size_t
Polyline_constraint_hierarchy_2<T,Compare,Data>::remove_points_without_corresponding_vertex(Constraint_id cid)
{
  std::size_t n = 0;
  for(Point_it it = points_in_constraint_begin(cid); 
      it != points_in_constraint_end(cid); ++it) { 
    if(cid.vl_ptr()->is_skipped(it.base())) {
      it = cid.vl_ptr()->erase(it.base());
      ++n;
    }
  }
  return n;
}

template <class T, class Compare, class Data>
std::size_t
Polyline_constraint_hierarchy_2<T,Compare,Data>::remove_points_without_corresponding_vertex()
{
  std::size_t n = 0;
  for(C_iterator it = constraint_set.begin(); it!= constraint_set.end(); ++it){
    n+= remove_points_without_corresponding_vertex(*it);
  }
  return n;
}


template <class T, class Compare, class Data>
typename Polyline_constraint_hierarchy_2<T,Compare,Data>::Constraint_id
Polyline_constraint_hierarchy_2<T,Compare,Data>::concatenate(Constraint_id first, Constraint_id second)
{
  constraint_set.erase(first);
  constraint_set.erase(second);
  // We have to look at all subconstraints
  for(Vertex_it it = second.vl_ptr()->skip_begin(), succ = it, end = second.vl_ptr()->skip_end(); 
      ++succ != end; 
      ++it){
    typename Sc_to_c_map::iterator scit = sc_to_c_map.find(make_edge(*it,*succ));
    CGAL_triangulation_assertion(scit != sc_to_c_map.end());
    Context_list* hcl = scit->second;

    // and replace the context of the constraint
    for(Context_iterator ctit=hcl->begin(); ctit != hcl->end(); ctit++) {
      if(ctit->enclosing == second.vl_ptr()){
	ctit->enclosing = first.vl_ptr();
	break;
      }
    }
  }
  // now we really concatenate the vertex lists
  // Note that all iterators pointing into second remain valid.
  // This concerns user code, as well as  the data member "pos" of the Context class
  first.vl_ptr()->pop_back(); // because it is the same as second.front()
  Vertex_it back_it = first.vl_ptr()->skip_end();
  --back_it;
  first.vl_ptr()->splice(first.vl_ptr()->skip_end(), *(second.vl_ptr()), second.vl_ptr()->skip_begin(), second.vl_ptr()->skip_end());

  // Note that for VC8 with iterator debugging the iterators pointing into second
  // are NOT valid      So we have to update them
  for(Vertex_it it = back_it, succ = it, end = first.vl_ptr()->skip_end(); 
      ++succ != end; 
      ++it){
    typename Sc_to_c_map::iterator scit = sc_to_c_map.find(make_edge(*it,*succ));
    CGAL_triangulation_assertion(scit != sc_to_c_map.end());
    Context_list* hcl = scit->second;

    // and update pos in the context of the constraint
    for(Context_iterator ctit=hcl->begin(); ctit != hcl->end(); ctit++) {
      if(ctit->enclosing == first.vl_ptr()){
	ctit->pos = it;
	break;
      }
    }
    }
  constraint_set.insert(first);

  delete second.vl_ptr();
  return first;
}

template <class T, class Compare, class Data>
typename Polyline_constraint_hierarchy_2<T,Compare,Data>::Constraint_id
Polyline_constraint_hierarchy_2<T,Compare,Data>::concatenate2(Constraint_id first, Constraint_id second)
{  
  constraint_set.erase(first);
  constraint_set.erase(second);
  // We have to look at all subconstraints
  for(Vertex_it it = first.vl_ptr()->skip_begin(), succ = it, end = first.vl_ptr()->skip_end(); 
      ++succ != end; 
      ++it){
    typename Sc_to_c_map::iterator scit = sc_to_c_map.find(make_edge(*it,*succ));
    CGAL_triangulation_assertion(scit != sc_to_c_map.end());
    Context_list* hcl = scit->second;

    // and replace the context of the constraint
    for(Context_iterator ctit=hcl->begin(); ctit != hcl->end(); ctit++) {
      if(ctit->enclosing == first.vl_ptr()){
	ctit->enclosing = second.vl_ptr();
	break;
      }
    }
  }
  // now we really concatenate the vertex lists
  // Note that all iterators pointing into second remain valid.
  first.vl_ptr()->pop_back(); // because it is the same as second.front()
  Vertex_it back_it = first.vl_ptr()->skip_end();
  --back_it;
  second.vl_ptr()->splice(second.vl_ptr()->skip_begin(), *(first.vl_ptr()), first.vl_ptr()->skip_begin(), first.vl_ptr()->skip_end());

  // Note that for VC8 with iterator debugging the iterators pointing into second
  // are NOT valid      So we have to update them
  for(Vertex_it it = back_it, succ = it, end = first.vl_ptr()->skip_end(); 
      ++succ != end; 
      ++it){
    typename Sc_to_c_map::iterator scit = sc_to_c_map.find(make_edge(*it,*succ));
    CGAL_triangulation_assertion(scit != sc_to_c_map.end());
    Context_list* hcl = scit->second;

    // and update pos in the context of the constraint
    for(Context_iterator ctit=hcl->begin(); ctit != hcl->end(); ctit++) {
      if(ctit->enclosing == second.vl_ptr()){
	ctit->pos = it;
	break;
      }
    }
  }
  constraint_set.insert(second);

  delete first.vl_ptr();
  return second.vl_ptr();
}


  // split a constraint in two constraints, so that vcit becomes the first
  // vertex of the new constraint
  // returns the new constraint 
template <class T, class Compare, class Data>
typename Polyline_constraint_hierarchy_2<T,Compare,Data>::Constraint_id
Polyline_constraint_hierarchy_2<T,Compare,Data>::split(Constraint_id first, Vertex_it vcit)
{
  constraint_set.erase(first);
  Vertex_list* second = new Vertex_list;
  second->splice(second->skip_end(), *(first.vl_ptr()), vcit.base(), first.vl_ptr()->skip_end());
  first.vl_ptr()->push_back(second->front()); // Duplicate the common vertex
  Vertex_it vit = second->skip_begin();
  vit.input() = true;
  vit = first.vl_ptr()->skip_end();
  --vit;
  vit.input() = true;
  constraint_set.insert(first);
  constraint_set.insert(second);
 // We have to look at all subconstraints
  for(Vertex_it it = second->skip_begin(), succ = it, end = second->skip_end(); 
      ++succ != end; 
      ++it){
    typename Sc_to_c_map::iterator scit = sc_to_c_map.find(make_edge(*it,*succ));
    CGAL_triangulation_assertion(scit != sc_to_c_map.end());
    Context_list* hcl = scit->second;

    // and replace the context of the constraint
    for(Context_iterator ctit=hcl->begin(); ctit != hcl->end(); ctit++) {
      if(ctit->enclosing == first.vl_ptr()){
	ctit->enclosing = second;
	break;
      }
    }
  }
  return second;
}

template <class T, class Compare, class Data>
typename Polyline_constraint_hierarchy_2<T,Compare,Data>::Constraint_id
Polyline_constraint_hierarchy_2<T,Compare,Data>::split2(Constraint_id first, Vertex_it vcit)
{
  constraint_set.erase(first);
  Vertex_list* second = new Vertex_list;
  second->splice(second->skip_end(), *first.vl_ptr(), first.vl_ptr()->skip_begin(), vcit.base());
  second->push_back(first.vl_ptr()->front()); // Duplicate the common vertex
  Vertex_it vit = second->skip_end();
  --vit;
  vit.input() = true;
  vit = first.vl_ptr()->skip_begin();
  vit.input() = true;
  constraint_set.insert(first);
  constraint_set.insert(second);
 // We have to look at all subconstraints
  for(Vertex_it it = second->skip_begin(), succ = it, end = second->skip_end(); 
      ++succ != end; 
      ++it){
    typename Sc_to_c_map::iterator scit = sc_to_c_map.find(make_edge(*it,*succ));
    CGAL_triangulation_assertion(scit != sc_to_c_map.end());
    Context_list* hcl = scit->second;

    // and replace the context of the constraint
    for(Context_iterator ctit=hcl->begin(); ctit != hcl->end(); ctit++) {
      if(ctit->enclosing == first.vl_ptr()){
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
template <class T, class Compare, class Data>
typename Polyline_constraint_hierarchy_2<T,Compare,Data>::Vertex_list*
Polyline_constraint_hierarchy_2<T,Compare,Data>::
insert_constraint(T va, T vb){
  Edge        he = make_edge(va, vb);
  Vertex_list*  children = new Vertex_list; 
  Context_list* fathers;

  typename Sc_to_c_map::iterator scit = sc_to_c_map.find(he);
  if(scit == sc_to_c_map.end()){
    fathers = new Context_list;
    sc_to_c_map.insert(std::make_pair(he,fathers));
  } else {
    fathers = scit->second;
  }

  children->push_front(Node(va, true));  // was he.first
  children->push_back(Node(vb, true));   // was he.second
  constraint_set.insert(children);
  Context ctxt;
  ctxt.enclosing = children;
  ctxt.pos     = children->skip_begin();
  fathers->push_front(ctxt);

  constraint_map[he] = children;
  return children;
}


template <class T, class Compare, class Data>
void
Polyline_constraint_hierarchy_2<T,Compare,Data>::
append_constraint(Constraint_id cid, T va, T vb){
  Edge        he = make_edge(va, vb);
  Context_list* fathers;

  typename Sc_to_c_map::iterator scit = sc_to_c_map.find(he);
  if(scit == sc_to_c_map.end()){
    fathers = new Context_list;
    sc_to_c_map.insert(std::make_pair(he,fathers));
  } else {
    fathers = scit->second;
  }

  typename Vertex_list::skip_iterator bit = cid.vl_ptr()->skip_end();
  --bit;
  cid.vl_ptr()->push_back(Node(vb, true));
  Context ctxt;
  ctxt.enclosing = cid.vl_ptr();
  ctxt.pos     = bit;
  fathers->push_front(ctxt);
}


template <class T, class Compare, class Data>
void Polyline_constraint_hierarchy_2<T,Compare,Data>::
clear()
{
  C_iterator cit;
  Sc_iterator scit;
  // clean and delete vertices lists
  for(cit=constraint_set.begin(); cit != constraint_set.end();  cit++){
    cit->vl_ptr()->clear();
    delete cit->vl_ptr();
  }
  // clean and delete context lists
  for(scit=sc_to_c_map.begin(); scit != sc_to_c_map.end(); scit++){
    (*scit).second->clear();
    delete (*scit).second;
  }
  sc_to_c_map.clear();
  constraint_set.clear();
  constraint_map.clear();
}


template <class T, class Compare, class Data>
bool Polyline_constraint_hierarchy_2<T,Compare,Data>::
next_along_sc(T va, T vb, T& w) const
{
  // find the next vertex after vb along any enclosing constrained
  // return false if there is no ....
  Context_iterator  ctxtit, past;
  if(!get_contexts(va, vb, ctxtit, past)) CGAL_triangulation_assertion(false);

  Vertex_it pos;
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
void Polyline_constraint_hierarchy_2<T,Compare,Data>::
remove_Steiner(T v, T va, T vb)
{
  // remove a Steiner point
  CGAL_precondition(!is_constrained_vertex(v));
 
  Context_list*  hcl1;
  Context_list*  hcl2;
  if(!get_contexts(va,v,hcl1)) CGAL_triangulation_assertion(false);
  if(!get_contexts(v,vb,hcl2)) CGAL_triangulation_assertion(false);

  Vertex_it      pos;
  for(Context_iterator ctit=hcl1->begin(); ctit != hcl1->end(); ctit++){
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
void Polyline_constraint_hierarchy_2<T,Compare,Data>::
split_constraint(T va, T vb, T vc){
  add_Steiner(va, vb, vc);
}


template <class T, class Compare, class Data>
void 
Polyline_constraint_hierarchy_2<T,Compare,Data>::
add_Steiner(T va, T vb, T vc){
  Context_list* hcl=NULL;
  if(!get_contexts(va,vb,hcl)) CGAL_triangulation_assertion(false);

  Context_list* hcl2 = new  Context_list;

  Vertex_it      pos;
  Context  ctxt;
  for(Context_iterator ctit=hcl->begin(); ctit != hcl->end(); ctit++) {
    // insert vc in enclosing constraint
    pos = ctit->current();
    ++pos;
    pos = ctit->enclosing->insert(pos.base(), Node(vc));
    --pos;
    
    // set ctxt to the context of (vc,vb)
    // change *ctit in hcl to the context of (va,vc)
    // add ctxt to hcl2 list
    ctxt.enclosing = ctit->enclosing;  
    if(*pos == va) {
      ctit->pos = pos;
      ctxt.pos = ++pos;
    }
    else { //(*pos)==vb
      ctxt.pos = pos;
      ctit->pos = ++pos;
    }
    hcl2->push_back(ctxt);
  }

  Context_list* hcl3;
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


template <class T, class Compare, class Data>
inline
typename Polyline_constraint_hierarchy_2<T,Compare,Data>::Edge
Polyline_constraint_hierarchy_2<T,Compare,Data>::
make_edge(T va, T vb) const
{
  return comp(va, vb) ? Edge(va,vb) : Edge(vb,va);
}

template <class T, class Compare, class Data>
inline
bool
Polyline_constraint_hierarchy_2<T,Compare,Data>::
get_contexts(T va, T vb, Context_list* & hcl) const
{
  Sc_iterator sc_iter = sc_to_c_map.find(make_edge(va,vb));
  if( sc_iter == sc_to_c_map.end() )    return(false);
  hcl = (*sc_iter).second;
  return true;
}

template <class T, class Compare, class Data>
inline
bool
Polyline_constraint_hierarchy_2<T,Compare,Data>::
get_contexts(T va, T vb, 
	     Context_iterator& ctxt, 
	     Context_iterator& past) const
{
  Context_list* hcl;
  if (!get_contexts(va,vb,hcl)) return false;
  ctxt = hcl->begin();
  past = hcl->end();
  return true;    
}



template <class T, class Compare, class Data>
inline
typename Polyline_constraint_hierarchy_2<T,Compare,Data>::Vertex_it
Polyline_constraint_hierarchy_2<T,Compare,Data>::
get_pos(T va, T vb) const
  //return pos in the first context
{
    return (*sc_to_c_map.find(make_edge(va,vb))).second->begin().pos;
}

template <class T, class Compare, class Data>
void
Polyline_constraint_hierarchy_2<T,Compare,Data>::
oriented_end(T va, T vb, T& vc) const
{
  Context_iterator ctxt, past;
  if(!get_contexts(va,vb, ctxt, past) ) CGAL_triangulation_assertion(false);
  if(*(ctxt->pos) == va)
    vc = ctxt->enclosing->back();
  else
    vc = ctxt->enclosing->front();
}


template <class T, class Compare, class Data>
void
Polyline_constraint_hierarchy_2<T,Compare,Data>::
print() const
{ 
  C_iterator hcit;
  std::map<T,int>  vertex_num;
  int num = 0;
  for(hcit = c_begin(); hcit != c_end();  hcit++) {
    Constraint_id cid = (*hcit);
    Vertex_it vit =cid.vl_ptr()->skip_begin(), end = cid.vl_ptr()->skip_end();
    for (;vit != end; vit++){
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

  C_iterator cit=c_begin();
  Sc_iterator scit=sc_begin();

  for(; cit != c_end();  cit++){
    std::cout << std::endl ;
    std::cout << "constraint " ;
    std::cout << cit->vl_ptr();
    std::cout << "  subconstraints " ;
    Vertex_it vit = (*cit).vl_ptr()->skip_begin(), end = (*cit).vl_ptr()->skip_end();
    for(; vit != end; vit++){
      std::cout << vertex_num[*vit]  <<" ";
    }
    vit = (*cit).vl_ptr()->skip_begin(), end = (*cit).vl_ptr()->skip_end();
    for(; vit != end; vit++){
      std::cout << (*vit)->point()  <<" ";
    }
  }
  std::cout << std::endl ;
  for(;scit != sc_end(); scit++){
    std::cout << "subconstraint " ;
    std::cout << vertex_num[scit->first.first] << " " 
	      << vertex_num[scit->first.second];
    Context_iterator cb, ce;
    get_contexts(scit->first.first, scit->first.second, cb, ce);
    
    std::cout << "  enclosing " ;
    for(; cb != ce; cb++) { 
      std::cout << cb->id().vl_ptr();
      std::cout <<  "   " ;
    }
    std::cout << std::endl ;
  }
  return; 
}


} //namespace CGAL
#endif // CGAL_POLYLINE_CONSTRAINT_HIERARCHY_2_H
