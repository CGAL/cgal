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
// release       : $CGAL_Revision: CGAL-2.4-I-75 $
// release_date  : $CGAL_Date: 2002/04/10 $
//
// file          : include/CGAL/nearest_neighbor_delaunay_2.h
// package       : Point_set_2 (2.3.2)
// maintainer    : Matthias Baesken <baesken@informatik.uni-trier.de>
// revision      : 2.3.2
// revision_date : 11 April 2002 
// author(s)     : Matthias Baesken
//
// coordinator   : Matthias Baesken, Trier  (<baesken@informatik.uni-trier.de>)
// ======================================================================

#ifndef CGAL_NEAREST_NEIGHBOR_DELAUNAY_2_H
#define CGAL_NEAREST_NEIGHBOR_DELAUNAY_2_H

#include <CGAL/basic.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <list>
#include <queue>
#include <map>
#include <stack>
#include <cmath>

CGAL_BEGIN_NAMESPACE


// compare function objects for the priority queues used in nearest neighbor search
template<class VP, class NT,class MAP_TYPE>
class compare_vertices {
 public:
  //std::map<VP,NT,std::less<VP> > *pmap;
  MAP_TYPE* pmap;
  
  compare_vertices(MAP_TYPE *p){ pmap=p; }
  
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
  typedef Unique_hash_map<Vertex_handle, Numb_type>         MAP_TYPE;

  int n = delau.number_of_vertices();
   
  if ( k <= 0 ) return;
  if ( n <= k ) { 
      get_vertices(delau, std::back_inserter(res)); 
      return; 
  }
     
  Point p = v->point();
     
  Unique_hash_map<Vertex_handle, int>  mark;
  int cur_mark = 1;

  MAP_TYPE  priority_number; // here we save the priorities ...
  
  compare_vertices<Vertex_handle,Numb_type,MAP_TYPE> comp(& priority_number);      // comparison object ...
  std::priority_queue<Vertex_handle, std::vector<Vertex_handle>, CGAL::compare_vertices<Vertex_handle,Numb_type,MAP_TYPE> > PQ(comp);

  priority_number[v.ptr()] = 0;
  PQ.push(v.ptr());
     
  // mark vertex v
  mark[v] = cur_mark;
      
  while ( k > 0 )
  { 
    // find minimum from PQ ...
    Vertex_handle w = PQ.top(); 
    PQ.pop();
   
    res.push_back(w);
    k--; 

    // get the incident vertices of w ...
    Vertex_circulator vc = delau.incident_vertices(w->handle());
    Vertex_circulator start =vc;
    Vertex_handle act;
     
    do {
         act = vc->handle();
	 
	 // test, if act is marked ...
	 bool is_marked = mark.is_defined(act);  
	 
         if ( (! is_marked) && (! delau.is_infinite(act)) )
         { 
             priority_number[act] = Compute_squared_distance_2()(p,act->point());
             PQ.push(act);
             mark[act] = cur_mark;
         }	   
	             
         vc++;
       } while (vc != start);   
        
   }
} 


CGAL_END_NAMESPACE

#endif


