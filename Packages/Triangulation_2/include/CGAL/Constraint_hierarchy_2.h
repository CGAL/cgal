// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 1.2
// release_date  : 2000/08/10 12:41:17
//
// file          : Constraint_hierarchy_2.h,v
// package       : 
// maintainer    : 
// revision      : 
//
// author(s)     : Olivier Billet
//
// coordinator   : INRIA Sophia Antipolis 
//                 (Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>)
//
// ======================================================================

#ifndef CGAL_CONSTRAINT_HIERARCHY_2_H
#define CGAL_CONSTRAINT_HIERARCHY_2_H

#include <pair.h>
#include <set.h>
#include <map.h> 
#include <list.h> 
#include <CGAL/triangulation_assertions.h>

// T               is expected to be Vertex_handle
// Data            is intended to store info on a Vertex_handle

template <class T, class Data>
class Constraint_hierarchy {
public:
  typedef std::pair<T, T>                      H_edge;
  typedef std::pair<T, Data>                   H_vertex;

  typedef std::pair<T, T>                      H_father;
  typedef std::list<T>                         H_children;
  typedef H_children::iterator                 H_pos;
  typedef struct{
    H_children*    child;
    H_pos          pos;
  }                                            H_context;

  typedef std::map<T, bool>                    H_vertex_map;
  typedef std::map<H_father, H_children*>      H_c_to_sc_map;
  typedef std::map<H_edge,   H_context>        H_sc_to_c_map;

  typedef typename H_c_to_sc_map::iterator     H_c_iterator;
  typedef typename H_sc_to_c_map::iterator     H_sc_iterator;
  typedef std::pair<H_father, H_children*>     H_c_value;
  typedef std::pair<H_edge,   H_context>       H_sc_value;
private:
  // data for the 1d hierarchy
  H_c_to_sc_map   c_to_sc_map;
  H_sc_to_c_map   sc_to_c_map;
  // data for the 0d hierarchy
  H_vertex_map    vertex_map;
  
public:
  Constraint_hierarchy() { }
  ~Constraint_hierarchy(){ }
  Constraint_hierarchy& operator=(const Constraint_hierarchy& ch);

  // Query 
  bool constraint_children(H_edge he,  H_children*& S);
  bool constraint_father(H_edge he, H_edge& hef);
  bool constraint_father(T  vca, T  vcb, T& vfa, T& vfb);
  void next_along_sc(T va, T vb, T& w);
  void oriented_end(T va, T vb, T& vc);

  // insert/remove
  void add_Steiner_in_sc(T va, T vb, T vx);
  void insert_constraint(T va, T vb);
  void remove_constraint(T va, T vb);
  void break_constraint(T va, T vb, T vc);
  
  void constrain_vertex(T v, Data data=Data());
  void unconstrain_vertex(T v);
  void set_data(T v, Data data);
  void get_data(T v, Data& data);

  void remove_of(T v, T va, T vb);
  void clear();

  // iterators
  H_sc_iterator sc_begin(){ return sc_to_c_map.begin(); }
  H_sc_iterator sc_end()  { return sc_to_c_map.end();   }
  H_c_iterator  c_begin() { return c_to_sc_map.begin(); }
  H_c_iterator  c_end()   { return c_to_sc_map.end();   }

  // checks
  bool is_subconstrained_edge(T va, T vb);
  bool is_constrained_edge(T va, T vb);
  bool is_constrained_vertex(T v);
  
private:
  inline
  H_edge    make_edge(T va, T vb);
  inline
  H_pos     get_pos(T va, T vb);public:
  void      get_context(T va, T vb, H_context& ctxt);
  void      replace_children(H_edge he, H_children* S);
};

//--------------------------------------------------------------------------//
template <class T, class Data> 
Constraint_hierarchy<T,Data>&
Constraint_hierarchy<T,Data>::
operator=(const Constraint_hierarchy& ch){
  c_to_sc_map=ch.c_to_sc_map;
  sc_to_c_map=ch.sc_to_c_map;
  vertex_map=ch,vertex_map;
  return *this;
};


template <class T, class Data> 
bool Constraint_hierarchy<T,Data>::
is_constrained_vertex(T v)
{
  return( vertex_map.find(v) != vertex_map.end() );
}


template <class T, class Data> 
bool Constraint_hierarchy<T,Data>::
is_constrained_edge(T va, T vb)
{
  return( c_to_sc_map.find(make_edge(va, vb)) != c_to_sc_map.end() );
}

template <class T, class Data> 
bool Constraint_hierarchy<T,Data>::
is_subconstrained_edge(T va, T vb)
{
  return( sc_to_c_map.find(make_edge(va, vb)) != sc_to_c_map.end() );
}


template <class T, class Data> 
bool Constraint_hierarchy<T,Data>::
constraint_children(H_edge he, H_children*& S) 
{
  H_sc_iterator sc_iter = sc_to_c_map.find(make_edge(va, vb));
  if( sc_iter == sc_to_c_map.end() )
    return false;
  S = (*sc_iter).second.child;
  return true;
}


template <class T, class Data>
bool Constraint_hierarchy<T,Data>::
constraint_father(H_edge he, H_edge& hef)
{
  H_sc_iterator sc_iter = sc_to_c_map.find(he);
  if( sc_iter == sc_to_c_map.end() )
    return(false);
  H_children* child = (*sc_iter).second.child;
  hef = make_edge(child->front(),child->back());
  return(true);
}


template <class T, class Data>
bool Constraint_hierarchy<T,Data>::
constraint_father(T  vca, T  vcb, T& vfa, T& vfb)
{
  H_sc_iterator sc_iter = sc_to_c_map.find(make_edge(vca,vcb));
  if( sc_iter == sc_to_c_map.end() )
    return(false);
  H_children* child = (*sc_iter).second.child;
  vfa = child->front();
  vfb = child->back();
  return(true);
}


/*
  |he| est a la fois sc et c au debut.
  les types de |xA| et |xB| ne sont pas les memes, c'est pour cela que 
  l'on a besoin des deux.
 */
template <class T, class Data>
void Constraint_hierarchy<T,Data>::
insert_constraint(T va, T vb){
  H_edge        he = make_edge(va, vb);
  H_children*   child = new H_children; //est-ce necessaire?
  H_c_value     x_c;
  H_sc_value    x_sc;
  
  x_c.first=he;

  child->clear();
  child->push_front(he.first);
  child->push_back(he.second);
  x_c.second=child;
  c_to_sc_map.insert(x_c);
  
  x_sc.first = he;
  x_sc.second.child= child;
  x_sc.second.pos  = child->begin();
  sc_to_c_map.insert(x_sc);

  //vertex
  constrain_vertex(va);
  constrain_vertex(vb);
}


template <class T, class Data>
void Constraint_hierarchy<T,Data>::
constrain_vertex(T v, Data data){
  H_vertex hv;
  hv.first =v;
  hv.second=data;
  vertex_map.insert(hv);
}


template <class T, class Data>
void Constraint_hierarchy<T,Data>::
unconstrain_vertex(T v){
  vertex_map.erase(v);
}


template <class T, class Data>
void Constraint_hierarchy<T,Data>::
get_data(T v, Data& data){
  //CGAL_precondition( is_constraint(v) );
  data = (*vertex_map.find(v)).second;
}


template <class T, class Data>
void Constraint_hierarchy<T,Data>::
set_data(T v, Data data){
  H_vertex hv;
  vertex_map.erase(v);
  hv.first =v;
  hv.second=data;
  vertex_map.insert(hv);
}


template <class T, class Data>
void Constraint_hierarchy<T,Data>::
clear(){
  sc_to_c_map.clear();
  c_to_sc_map.clear();
  vertex_map.clear();
}



/* 
   Attention, le point vb DOIT etre un point de Steiner,
   et le segment va,vb ne doit pas contenir de sous contrainte.
*/
template <class T, class Data>
void Constraint_hierarchy<T,Data>::
next_along_sc(T va, T vb, T& w){
  /*
  CGAL_precondition(is_subconstrained_edge(va,vb) &&
                    !is_constrained_vertex(vb));
  */ 
  H_pos pos = get_pos(va,vb);
  if((*pos)==va) {pos++; pos++;} else pos--;
  w=(*pos);
}



/*
  Attention, le point v DOIT etre un point de Steiner,
  et les segments va,v et v,vb sont des sous contraintes.
*/
template <class T, class Data>
void Constraint_hierarchy<T,Data>::
remove_of(T v, T va, T vb){
  CGAL_precondition(!is_constrained_vertex(v));
  CGAL_precondition(is_subconstrained_edge(va,v));
  CGAL_precondition(is_subconstrained_edge(vb,v));

  H_sc_value x_sc;
  H_pos      pos;

  get_context(va,v,x_sc.second);
  pos = x_sc.second.pos;

  if((*pos)==va) pos++;
  
  pos = x_sc.second.child->erase(pos);
  pos--;

  sc_to_c_map.erase(make_edge(va,v));
  sc_to_c_map.erase(make_edge(v,vb));

  x_sc.first = make_edge(va,vb);
  x_sc.second.pos= pos;
  sc_to_c_map.insert(x_sc);

}



/*
  precondition : va,vb est une souscontrainte. 
  attention: data not already set !!!
*/
template <class T, class Data>
void Constraint_hierarchy<T,Data>::
break_constraint(T va, T vb, T vc){
  add_Steiner_in_sc(va, vb, vc);
  constrain_vertex(vc,false);
}


template <class T, class Data>
void Constraint_hierarchy<T,Data>::
add_Steiner_in_sc(T va, T vb, T vc){
  H_sc_value x_sc;
  H_pos      pos;
  bool       b=false;

  get_context(va,vb,x_sc.second);
  pos = x_sc.second.pos;
  pos++;  
  pos = x_sc.second.child->insert(pos, vc);
  pos--;
  if((*pos)==va) b=true;

  if(b)
    x_sc.first = make_edge(va,vc);
  else
    x_sc.first = make_edge(vb,vc);

  x_sc.second.pos= pos;
  sc_to_c_map.insert(x_sc);
  
  pos++;
  if(b)
    x_sc.first = make_edge(vb,vc);
  else
    x_sc.first = make_edge(va,vc);
  x_sc.second.pos = pos;
  sc_to_c_map.insert(x_sc); 
  sc_to_c_map.erase(make_edge(va,vb)); 
}


template <class T, class Data>
inline
Constraint_hierarchy<T,Data>::H_edge
Constraint_hierarchy<T,Data>::
make_edge(T va, T vb){
  Constraint_hierarchy<T,Data>::H_edge he;
  if(va<vb)
    { he.first = va; he.second = vb; }
  else
    { he.first = vb; he.second = va; }
  return(he);
}


template <class T, class Data>
inline
void
Constraint_hierarchy<T,Data>::
get_context(T va, T vb, H_context& ctxt){
  ctxt = (*sc_to_c_map.find(make_edge(va,vb))).second;
}

template <class T, class Data>
inline
Constraint_hierarchy<T,Data>::H_pos
Constraint_hierarchy<T,Data>::
get_pos(T va, T vb){
  return (*sc_to_c_map.find(make_edge(va,vb))).second.pos;
}

template <class T, class Data>
void
Constraint_hierarchy<T,Data>::
oriented_end(T va, T vb, T& vc){
  CGAL_precondition(is_subconstrained_edge(va,vb));
  H_context ctxt;
  get_context(va,vb,ctxt);
  if((*ctxt.pos)==va)
    vc = ctxt.child->back();
  else
    vc = ctxt.child->front();
}

#endif // CGAL_CONSTRAINT_HIERARCHY_2_H
