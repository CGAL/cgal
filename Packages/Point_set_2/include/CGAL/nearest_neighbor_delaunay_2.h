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
// release       : 
// release_date  : 
//
// file          : include/CGAL/nearest_neighbor_delaunay_2.h
// package       : Point_set_2 (2.2.1)
// maintainer    : Matthias Baesken <baesken@informatik.uni-trier.de>
// revision      : 2.2.1
// revision_date : 10 July 2001 
// author(s)     : Matthias Baesken
//
// coordinator   : Matthias Baesken, Trier  (<baesken@informatik.uni-trier.de>)
// ======================================================================

#ifndef CGAL_NEAREST_NEIGHBOR_DELAUNAY_2_H
#define CGAL_NEAREST_NEIGHBOR_DELAUNAY_2_H

#include <CGAL/basic.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <list>
#include <queue>
#include <map>
#include <stack>
#include <cmath>

CGAL_BEGIN_NAMESPACE


// compare function objects for the priority queues used in nearest neighbor search
template<class VP, class NT>
class compare_vertices {
 public:
  std::map<VP,NT,std::less<VP> > *pmap;
  
  compare_vertices(std::map<VP,NT,std::less<VP> > *p){ pmap=p; }
  
  bool operator()(VP e1, VP e2)
  // get the priorities from the map and return result of comparison ...
  { NT& v1 = (*pmap)[e1];
    NT& v2 = (*pmap)[e2];
    return (v1 > v2);
  }
};


template<class Dt>
typename Dt::Vertex_handle  nearest_neighbor(const Dt& delau, typename Dt::Vertex_handle v)
{
  typedef typename Dt::Geom_traits                    Gt;
  typedef typename Dt::Point                          Point;  
  typedef typename Dt::Vertex_circulator              Vertex_circulator;
  typedef typename Dt::Vertex_handle                  Vertex_handle;
  typedef typename Gt::Compare_distance_2             Compare_dist_2;
  
  if (delau.number_of_vertices() <= 1) return NULL;    
  Point p = v->point();
     
  Vertex_circulator vc = delau.incident_vertices(v);  
  Vertex_circulator start =vc;
     
  Vertex_handle min_v = (*vc).handle();
  if (delau.is_infinite(min_v)){
     vc++;
     min_v = (*vc).handle();
  }
     
  Vertex_handle act;
     
  // go through the vertices ...
  do {
       act = vc->handle();
 
       if (! delau.is_infinite(act)) {
        if ( Compare_dist_2()(p,act->point(), min_v->point()) == SMALLER ) min_v = act;
       }   
           
       vc++;
  } while (vc != start);     
   
  return min_v;
}


template<class Dt>
typename Dt::Vertex_handle lookup(const Dt& delau, const typename Dt::Point& p)
{ 
 typedef typename Dt::Face                      Face;
 typedef typename Dt::Locate_type               Locate_type;
 typedef typename Dt::Face_handle               Face_handle;
  
 if (delau.number_of_vertices() == 0) return NULL;   
     
 // locate ...
 Locate_type lt;
 int li;
 Face_handle fh = delau.locate(p,lt,li);
     
 if (lt == Dt::VERTEX){
      Face f = *fh;
      return f.vertex(li);
 }
 else return NULL;
}   
   
   
template<class Dt, class OutputIterator>
OutputIterator   nearest_neighbors(Dt& delau, const typename Dt::Point& p, int k, OutputIterator res)
{
  typedef typename Dt::Geom_traits                    Gt;
  typedef typename Dt::Vertex_handle                  Vertex_handle;
  typedef typename Dt::Vertex_iterator                Vertex_iterator;

   int n = delau.number_of_vertices();

   if ( k <= 0 ) return res;
   if ( n <= k ) { // return all finite vertices ...
   
     Vertex_iterator vit = delau.finite_vertices_begin();
     for (; vit != delau.finite_vertices_end(); vit++) { *res= (*vit).handle(); res++; }    
   
     return res;
   }

   // insert p, if necessary
   
   Vertex_handle vh = lookup(delau,p);
   bool old_node = true;
    
   // we have to add a new vertex ...
   if (vh == NULL){
      vh = delau.insert(p);
      old_node = false;
      k++;
   }

   std::list<Vertex_handle> res_list;
   nearest_neighbors_list(delau, vh, k, res_list);
   
   if ( !old_node ) 
   { 
     res_list.pop_front();
     delau.remove(vh);
   }
    
   typename std::list<Vertex_handle>::iterator it = res_list.begin();
    
   for (; it != res_list.end(); it++) { *res= *it; res++; }
   return res;  
}


   
template<class Dt, class OutputIterator>  
OutputIterator  nearest_neighbors(const Dt& delau, typename Dt::Vertex_handle v, int k, OutputIterator res)
{  
  typedef typename Dt::Geom_traits                    Gt;
  typedef typename Dt::Vertex_handle                  Vertex_handle;
  typedef typename Dt::Vertex_iterator                Vertex_iterator;

   int n = delau.number_of_vertices();

   if ( k <= 0 ) return res;
   if ( n <= k ) { // return all (finite) vertices ...
   
     Vertex_iterator vit = delau.finite_vertices_begin();
     for (; vit != delau.finite_vertices_end(); vit++) { *res= (*vit).handle(); res++; }    
   
     return res;
   }
   
   std::list<Vertex_handle> res_list;
   nearest_neighbors_list(delau, v, k, res_list); 
   
   typename std::list<Vertex_handle>::iterator it = res_list.begin();
    
   for (; it != res_list.end(); it++) { *res= *it; res++; }

   return res;     
}


template<class Dt, class OutputIterator>
OutputIterator get_vertices(const Dt& delau, OutputIterator res)
{    
  typedef typename Dt::Vertex_iterator                Vertex_iterator;

  Vertex_iterator vit = delau.finite_vertices_begin();
  for (; vit != delau.finite_vertices_end(); vit++) { *res= (*vit).handle(); res++; }  
  return res;   
}


// second template argument for VC ...

template<class Dt, class T2>
void nearest_neighbors_list(const Dt& delau, typename Dt::Vertex_handle v, int k, std::list<T2>& res) 
{  
  typedef typename Dt::Geom_traits                    Gt;
  typedef typename Dt::Vertex_handle                  Vertex_handle;
  typedef typename Dt::Vertex_iterator                Vertex_iterator;
  typedef typename Dt::Vertex_circulator              Vertex_circulator;
  typedef typename Dt::Vertex                         Vertex;
  typedef typename Dt::Point                          Point;
  typedef typename Gt::FT                             Numb_type;  // field number type ...
  typedef typename Gt::Compute_squared_distance_2     Compute_squared_distance_2;  
  typedef typename std::map<Vertex*,int, std::less<Vertex*> >::iterator map_iterator;

  int n = delau.number_of_vertices();
   
  if ( k <= 0 ) return;
  if ( n <= k ) { 
      get_vertices(delau, std::back_inserter(res)); 
      return; 
  }
     
  Point p = v->point();
     
  std::map<Vertex*,int, std::less<Vertex*> > mark;
  int cur_mark = 1;

  std::map<Vertex*,Numb_type, std::less<Vertex*> > priority_number; // here we save the priorities ...
  compare_vertices<Vertex*,Numb_type> comp(& priority_number);      // comparison object ...
  std::priority_queue<Vertex*, std::vector<Vertex*>, CGAL::compare_vertices<Vertex*,Numb_type> > PQ(comp);

  priority_number[v.ptr()] = 0;
  PQ.push(v.ptr());
     
  //mark_vertex(v);
  Vertex* vptr = v.ptr();
  mark[vptr] = cur_mark;
      
  while ( k > 0 )
  { 
    // find minimum from PQ ...
    Vertex* w = PQ.top(); PQ.pop();
    priority_number.erase(w); // and clean entry in priority map
   
    res.push_back(w->handle()); k--; 

    // get the incident vertices of w ...
    Vertex_circulator vc = delau.incident_vertices(w->handle());
    Vertex_circulator start =vc;
    Vertex* act;
     
    do {
         act = &(*vc);
	 
	 // test, if act is marked ...
	 bool is_marked;
    
         map_iterator mit = mark.find(act);   
         if (mit == mark.end()) is_marked=false;
         else is_marked=true;
	 
         if ( (! is_marked) && (! delau.is_infinite(act->handle())) )
         { 
             priority_number[act] = Compute_squared_distance_2()(p,act->point());
             PQ.push(act);	      
             //mark_vertex(act->handle());
	     Vertex* vptr = act->handle().ptr();
             mark[vptr] = cur_mark;
         }	   
	             
         vc++;
       } while (vc != start);   
        
   }
} 


CGAL_END_NAMESPACE

#endif


