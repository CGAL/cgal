// Copyright (c) 1999  
// Max-Planck-Institute Saarbruecken (Germany). All rights reserved.
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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Matthias Baesken

#ifndef CGAL_NEAREST_NEIGHBOR_DELAUNAY_2_H
#define CGAL_NEAREST_NEIGHBOR_DELAUNAY_2_H

#include <CGAL/license/Point_set_2.h>


#include <CGAL/basic.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/squared_distance_2_1.h>
#include <CGAL/compare_vertices.h>
#include <list>
#include <queue>
#include <map>
#include <stack>
#include <cmath>

namespace CGAL {



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
     
  Vertex_handle min_v = vc;
  if (delau.is_infinite(min_v)){
     vc++;
     min_v = vc;
  }
     
  Vertex_handle act;
     
  // go through the vertices ...
  do {
       act = vc;
 
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
OutputIterator   nearest_neighbors(Dt& delau, const typename Dt::Point& p, std::size_t k, OutputIterator res)
{
  typedef typename Dt::size_type                      size_type;
  typedef typename Dt::Vertex_handle                  Vertex_handle;
  typedef typename Dt::Vertex_iterator                Vertex_iterator;

   size_type n = delau.number_of_vertices();

   if ( k <= 0 ) return res;
   if ( n <= k ) { // return all finite vertices ...
   
     Vertex_iterator vit = delau.finite_vertices_begin();
     for (; vit != delau.finite_vertices_end(); vit++) { *res= vit; res++; }    
   
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
OutputIterator  nearest_neighbors(const Dt& delau, typename Dt::Vertex_handle v, std::size_t k, OutputIterator res)
{  
  typedef typename Dt::size_type                      size_type;
  typedef typename Dt::Vertex_handle                  Vertex_handle;
  typedef typename Dt::Vertex_iterator                Vertex_iterator;

   size_type n = delau.number_of_vertices();

   if ( k <= 0 ) return res;
   if ( n <= k ) { // return all (finite) vertices ...
   
     Vertex_iterator vit = delau.finite_vertices_begin();
     for (; vit != delau.finite_vertices_end(); vit++) { *res= vit; res++; }    
   
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
  for (; vit != delau.finite_vertices_end(); vit++) { *res= vit; res++; }  
  return res;   
}


// second template argument for VC ...

template<class Dt, class T2>
void nearest_neighbors_list(const Dt& delau, typename Dt::Vertex_handle v, std::size_t k, std::list<T2>& res) 
{  
  typedef typename Dt::Geom_traits                    Gt;
  typedef typename Dt::size_type                      size_type;
  typedef typename Dt::Vertex_handle                  Vertex_handle;
  typedef typename Dt::Vertex_circulator              Vertex_circulator;
  typedef typename Dt::Point                          Point;
  typedef typename Gt::FT                             Numb_type;  // field number type ...
  typedef typename Gt::Compute_squared_distance_2     Compute_squared_distance_2;   
  typedef Unique_hash_map<Vertex_handle, Numb_type>         MAP_TYPE;

  size_type n = delau.number_of_vertices();
   
  if ( k <= 0 ) return;
  if ( n <= k ) { 
      get_vertices(delau, std::back_inserter(res)); 
      return; 
  }
     
  Point p = v->point();
     
  Unique_hash_map<Vertex_handle, int>  mark;
  int cur_mark = 1;

  MAP_TYPE  priority_number; // here we save the priorities ...
  
  internal::compare_vertices<Vertex_handle,Numb_type,MAP_TYPE> comp(& priority_number);      // comparison object ...
  std::priority_queue<Vertex_handle, std::vector<Vertex_handle>, internal::compare_vertices<Vertex_handle,Numb_type,MAP_TYPE> > PQ(comp);

  priority_number[v] = 0;
  PQ.push(v);
     
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
    Vertex_circulator vc = delau.incident_vertices(w);
    Vertex_circulator start =vc;
    Vertex_handle act;
     
    do {
         act = vc;
	 
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


} //namespace CGAL

#endif
