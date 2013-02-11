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
// 
//
// Author(s)     : Matthias Baesken

#ifndef CGAL_RANGE_SEARCH_DELAUNAY_2_H
#define CGAL_RANGE_SEARCH_DELAUNAY_2_H

#include <CGAL/basic.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/nearest_neighbor_delaunay_2.h>
#include <list>
#include <queue>
#include <map>
#include <stack>
#include <cmath>

namespace CGAL {

// was 
// std::map<typename Dt::Vertex *,int, std::less<typename Dt::Vertex *> >

template<class Dt,class Circle,class OutputIterator,class MAP_TYPE>
void dfs(const Dt& delau, 
         MAP_TYPE& mark, 
         typename Dt::Vertex_handle v, 
	 const Circle& C, 
	 OutputIterator res,
	 bool first_vertex)
{
    typedef typename Dt::Geom_traits               Gt;
    typedef typename Dt::Vertex_circulator         Vertex_circulator;
    typedef typename Dt::Vertex_handle             Vertex_handle;
    typedef typename Gt::Bounded_side_2            Bounded_side_2;

    Bounded_side_2  test;
    
    if (! first_vertex) *res = v;
    
    //mark_vertex v
    mark[v] = true;
    
    // get incident vertices of v ...
    Vertex_circulator vc = delau.incident_vertices(v);
    Vertex_circulator start =vc;
     
    Vertex_handle act;
     
    // go through the vertices ...
    do {
       act = vc;
 
       if (! delau.is_infinite(act)) {
	 // test, if act is marked ...
         bool is_marked = mark[act];     
       
         if ((! is_marked) && !( test(C,act->point()) == ON_UNBOUNDED_SIDE) ) 
           dfs(delau, mark, act, C, res, false);       
       }             
       vc++;
    } while (vc != start);     
}


// second dfs uses test - predicate function object ...

template<class Dt,class Circle,class OutputIterator,class MAP_TYPE,class Pred>
bool dfs(const Dt& delau, 
         MAP_TYPE& mark, 
         typename Dt::Vertex_handle v, 
	 const Circle& C, 
	 OutputIterator res,
	 bool first_vertex,
	 bool return_if_predicate_succeded,
	 Pred& pred)
{
    typedef typename Dt::Geom_traits               Gt;
    typedef typename Dt::Vertex_circulator         Vertex_circulator;
    typedef typename Dt::Vertex_handle             Vertex_handle;
    typedef typename Gt::Bounded_side_2            Bounded_side_2;

    Bounded_side_2  test;
    
    if (! first_vertex) {
      if (pred(v->point())) {
         *res = v;
	 
	 if (return_if_predicate_succeded) return true;
      }
    }
    
    //mark_vertex v
    mark[v] = true;
    
    // get incident vertices of v ...
    Vertex_circulator vc = delau.incident_vertices(v);
    Vertex_circulator start =vc;
     
    Vertex_handle act;
     
    // go through the vertices ...
    do {
       act = vc;
 
       if (! delau.is_infinite(act)) {
	 // test, if act is marked ...
         bool is_marked = mark[act];    
       
         if ((! is_marked) && !( test(C,act->point()) == ON_UNBOUNDED_SIDE) ) 
           if (dfs(delau, mark, act, C, res, false, return_if_predicate_succeded, pred)) return true;       
       }             
       vc++;
    } while (vc != start);   
    
    return false;  
}


template<class Dt,class Circle,class OutputIterator,class MAP_TYPE,class Pred>
void dfs_using_predicate(const Dt& delau, 
         MAP_TYPE& mark, 
         typename Dt::Vertex_handle v, 
	 const Circle& C, 
	 OutputIterator res,
	 bool first_vertex,
	 bool return_if_predicate_succeded,
	 Pred& pred)
{
  bool val = dfs(delau, mark, v, C, res, first_vertex, return_if_predicate_succeded, pred);
  pred.set_result(val);
}

// circular range search ...


template<class Dt, class Circle, class OutputIterator>
OutputIterator range_search(Dt& delau, 
                            const Circle& C, 
			    OutputIterator res)
{ 
  typedef typename Dt::Geom_traits                    Gt;
  typedef typename Dt::Point                          Point;
  typedef typename Dt::Vertex_handle                  Vertex_handle;
  typedef typename Dt::Vertex_iterator                Vertex_iterator;
  typedef typename Gt::Bounded_side_2                 Bounded_side_2;
  typedef typename Gt::Construct_center_2             Construct_center_2;

  Bounded_side_2  test;

  if (delau.number_of_vertices() == 0) return res;
  if (delau.number_of_vertices() == 1) 
  { 
       // get the one vertex ...
       Vertex_iterator vit = delau.finite_vertices_begin();
       Point p = vit->point();
       
       if (! (test(C, p) == ON_UNBOUNDED_SIDE)){
        *res= vit; res++;
       }
       return res;
   }  
     
   // normal case ...
   Point p = Construct_center_2()(C);
   Vertex_handle v = lookup(delau, p);  
   bool new_v = false;     

   // we have to insert the center ...
   if ( v == NULL )
   { 
       new_v = true;
       v = delau.insert(p); 
   }
     
   //std::map<Vertex*,int, std::less<Vertex*> > mark;
   Unique_hash_map<Vertex_handle, bool> mark;
     
   dfs(delau,mark,v,C,res,new_v);
     
   if (new_v)
   { 
     delau.remove(v);
   }    
   
   return res;        
}
   

   
// triangular range search ...   
// Note that the function only works correctly with exact constructions
// because it computes the circumcenter of the points a, b, and c
// and then performs a range query with this circle. 
// When vertices of the trinagulation are on the circle the outcome
// is not deterministic.
// A solution would be to not constuct a circle, but to use the
// function CGAL::side_of_bounded_circle

template<class Dt, class OutputIterator>
OutputIterator range_search(Dt& delau, 
                            const typename Dt::Point& a, 
			    const typename Dt::Point& b, 
                            const typename Dt::Point& c,
			    OutputIterator res)
{
  typedef typename Dt::Geom_traits                    Gt;
  typedef typename Dt::Point                          Point;
  typedef typename Gt::Circle_2                       Circle;
  typedef typename Dt::Vertex_handle                  Vertex_handle;
  typedef typename Gt::Orientation_2                  Orientation_2;
  typedef typename Gt::Construct_circle_2             Construct_circle_2;
  
  Orientation_2 test_ori;
   
  int orient = (int)(test_ori(a,b,c));
  Circle C = Construct_circle_2()(a,b,c);
  std::list<Vertex_handle> L;
  range_search(delau, C, std::back_inserter(L));
      
  typename std::list<Vertex_handle>::const_iterator it = L.begin();
      
  for(;it != L.end();it++)
  { Point p = (*it)->point();
    if ( ((int)(test_ori(a,b,p))) == - orient ||
         ((int)(test_ori(b,c,p))) == - orient ||
         ((int)(test_ori(c,a,p))) == - orient ) { }     
    else { *res = *it; res++; }
  }
  return res;     
}


// rectangular range search ....

template<class Dt,class OutputIterator>
OutputIterator range_search(Dt& delau, 
                            const typename Dt::Point& a1, 
			    const typename Dt::Point& b1, 
                            const typename Dt::Point& c1,
			    const typename Dt::Point& d1,
			    OutputIterator res)
// a1 upper left, b1 lower left , c1 lower right
{
  typedef typename Dt::Geom_traits                    Gt;
  typedef typename Dt::Point                          Point;
  typedef typename Gt::Circle_2                       Circle;
  typedef typename Dt::Vertex_handle                  Vertex_handle;
  typedef typename Gt::Orientation_2                  Orientation_2;
  typedef typename Gt::Construct_circle_2             Construct_circle_2;
  
  Point a=a1,b=b1,c=c1,d=d1;
    
  if (Orientation_2()(a,b,c) == RIGHT_TURN) 
  { Point tmp = b;
    b = d;
    d = tmp;
  }
   
  Circle C = Construct_circle_2()(a,b,c);
     
  std::list<Vertex_handle> L;
  range_search(delau, C, std::back_inserter(L));
  typename std::list<Vertex_handle>::const_iterator it = L.begin();     

  for(;it != L.end();it++)
  { Point p = (*it)->point();
    if ( Orientation_2()(a,b,p) == RIGHT_TURN || Orientation_2()(b,c,p) == RIGHT_TURN ||
         Orientation_2()(c,d,p) == RIGHT_TURN || Orientation_2()(d,a,p) == RIGHT_TURN )  { }
    else { *res = *it; res++; }
  }
  return res;     
}


// ------------------------------------------------------------------------------------------------
// new range search variants using test function object ...

// circular range search ...


template<class Dt, class Circle, class OutputIterator,class Pred>
OutputIterator range_search(Dt& delau, 
                            const Circle& C, 
			    OutputIterator res,
			    Pred& pred,
			    bool return_if_succeded)
{ 
  typedef typename Dt::Geom_traits                    Gt;
  typedef typename Dt::Point                          Point;
  typedef typename Dt::Vertex_handle                  Vertex_handle;
  typedef typename Dt::Vertex_iterator                Vertex_iterator;
  typedef typename Gt::Bounded_side_2                 Bounded_side_2;
  typedef typename Gt::Construct_center_2             Construct_center_2;

  Bounded_side_2  test;

  if (delau.number_of_vertices() == 0) return res;
  if (delau.number_of_vertices() == 1) 
  { 
       // get the one vertex ...
       Vertex_iterator vit = delau.finite_vertices_begin();
       Point p = vit->point();
       
       if (! (test(C, p) == ON_UNBOUNDED_SIDE)){
        *res= vit; res++;
       }
       
       bool val = pred(p);
       pred.set_result(val);
       
       return res;
   }  
     
   // normal case ...
   Point p = Construct_center_2()(C);
   Vertex_handle v = lookup(delau, p);  
   bool new_v = false;     

   // we have to insert the center ...
   if ( v == NULL )
   { 
       new_v = true;
       v = delau.insert(p); 
   }
     
   //std::map<Vertex*,int, std::less<Vertex*> > mark;
   Unique_hash_map<Vertex_handle, bool> mark;
   
   dfs_using_predicate(delau,mark,v,C,res,new_v, return_if_succeded, pred);
     
   if (new_v)
   { 
     delau.remove(v);
   }    
   
   return res;        
}
   
   
// triangular range search ...   

template<class Dt, class OutputIterator,class Pred>
OutputIterator range_search(Dt& delau, 
                            const typename Dt::Point& a, 
			    const typename Dt::Point& b, 
                            const typename Dt::Point& c,
			    OutputIterator res,
			    Pred& pred,
			    bool return_if_succeded)			    
{
  typedef typename Dt::Geom_traits                    Gt;
  typedef typename Dt::Point                          Point;
  typedef typename Gt::Circle_2                       Circle;
  typedef typename Dt::Vertex_handle                  Vertex_handle;
  typedef typename Gt::Orientation_2                  Orientation_2;
  typedef typename Gt::Construct_circle_2             Construct_circle_2;
  
  Orientation_2 test_ori;
   
  int orient = (int)(test_ori(a,b,c));
  Circle C = Construct_circle_2()(a,b,c);
  std::list<Vertex_handle> L;
  
  range_search(delau, C, std::back_inserter(L), pred, return_if_succeded);
  if (return_if_succeded) return res;
      
  typename std::list<Vertex_handle>::const_iterator it = L.begin();
      
  for(;it != L.end();it++)
  { Point p = (*it)->point();
    if ( ((int)(test_ori(a,b,p))) == - orient ||
         ((int)(test_ori(b,c,p))) == - orient ||
         ((int)(test_ori(c,a,p))) == - orient ) { }     
    else { *res = *it; res++; }
  }
  return res;     
}


// rectangular range search ....

template<class Dt,class OutputIterator,class Pred>
OutputIterator range_search(Dt& delau, 
                            const typename Dt::Point& a1, 
			    const typename Dt::Point& b1, 
                            const typename Dt::Point& c1,
			    const typename Dt::Point& d1,
			    OutputIterator res,
			    Pred& pred,
			    bool return_if_succeded)
// a1 upper left, b1 lower left , c1 lower right
{
  typedef typename Dt::Geom_traits                    Gt;
  typedef typename Dt::Point                          Point;
  typedef typename Gt::Circle_2                       Circle;
  typedef typename Dt::Vertex_handle                  Vertex_handle;
  typedef typename Gt::Orientation_2                  Orientation_2;
  typedef typename Gt::Construct_circle_2             Construct_circle_2;
  
  Point a=a1,b=b1,c=c1,d=d1;
    
  if (Orientation_2()(a,b,c) == RIGHT_TURN) 
  { Point tmp = b;
    b = d;
    d = tmp;
  }
   
  Circle C = Construct_circle_2()(a,b,c);
     
  std::list<Vertex_handle> L;
  
  range_search(delau, C, std::back_inserter(L), pred, return_if_succeded);
  if (return_if_succeded) return res;
  
  
  typename std::list<Vertex_handle>::const_iterator it = L.begin();     

  for(;it != L.end();it++)
  { Point p = (*it)->point();
    if ( Orientation_2()(a,b,p) == RIGHT_TURN || Orientation_2()(b,c,p) == RIGHT_TURN ||
         Orientation_2()(c,d,p) == RIGHT_TURN || Orientation_2()(d,a,p) == RIGHT_TURN )  { }
    else { *res = *it; res++; }
  }
  return res;     
}



} //namespace CGAL

#endif
