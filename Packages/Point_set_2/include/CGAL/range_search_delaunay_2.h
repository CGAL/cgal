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
// file          : include/CGAL/range_search_delaunay_2.h
// package       : Point_set_2 (2.2.1)
// maintainer    : Matthias Baesken <baesken@informatik.uni-trier.de>
// revision      : 2.2.1
// revision_date : 10 July 2001 
// author(s)     : Matthias Baesken
//
// coordinator   : Matthias Baesken, Trier  (<baesken@informatik.uni-trier.de>)
// ======================================================================

#ifndef CGAL_RANGE_SEARCH_DELAUNAY_2_H
#define CGAL_RANGE_SEARCH_DELAUNAY_2_H

#include <CGAL/basic.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/nearest_neighbor_delaunay_2.h>
#include <list>
#include <queue>
#include <map>
#include <stack>
#include <cmath>

CGAL_BEGIN_NAMESPACE



template<class Dt,class Circle,class T2>
void dfs(const Dt& delau, std::map<typename Dt::Vertex *,int, std::less<typename Dt::Vertex *> >& mark, 
         typename Dt::Vertex_handle v, const Circle& C, std::list<T2>& L)
{
    typedef typename Dt::Geom_traits               Gt;
    typedef typename Dt::Vertex                    Vertex;
    typedef typename Dt::Vertex_circulator         Vertex_circulator;
    typedef typename Dt::Vertex_handle             Vertex_handle;
    typedef typename Gt::Bounded_side_2            Bounded_side_2;
    typedef typename std::map<Vertex*,int, std::less<Vertex*> >::iterator map_iterator;
 
    Bounded_side_2  test;
    L.push_back(v);
    
    //mark_vertex v
    Vertex* vptr = v.ptr();
    mark[vptr] = 1;
    
    // get incident vertices of v ...
    Vertex_circulator vc = delau.incident_vertices(v);
    Vertex_circulator start =vc;
     
    Vertex_handle act;
     
    // go through the vertices ...
    do {
       act = vc->handle();
 
       if (! delau.is_infinite(act)) {
	 // test, if act is marked ...
	 bool is_marked;
	 Vertex* vtest = act.ptr();
    
         map_iterator mit = mark.find(vtest);   
         if (mit == mark.end()) is_marked=false;
         else is_marked=true;       
       
       
        if ((! is_marked) && !( test(C,act->point()) == ON_UNBOUNDED_SIDE) ) 
           dfs(delau, mark, act, C, L);       
       }             
       vc++;
    } while (vc != start);     
}


template<class Dt, class Circle, class OutputIterator>
OutputIterator range_search(Dt& delau, const Circle& C, OutputIterator res)
{ 
  typedef typename Dt::Geom_traits                    Gt;
  typedef typename Dt::Point                          Point;
  typedef typename Dt::Vertex_handle                  Vertex_handle;
  typedef typename Dt::Vertex                         Vertex;
  typedef typename Dt::Vertex_iterator                Vertex_iterator;
  typedef typename Gt::Bounded_side_2                 Bounded_side_2;
  typedef typename Gt::Construct_center_2             Construct_center_2;

  Bounded_side_2  test;

  if (delau.number_of_vertices() == 0) return res;
  if (delau.number_of_vertices() == 1) 
  { 
       // get the one vertex ...
       Vertex_iterator vit = delau.finite_vertices_begin();
       Vertex v = (*vit);
       Point p = v.point();
       
       if (! (test(C, p) == ON_UNBOUNDED_SIDE)){
        *res= v.handle(); res++;
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
     
   //init_dfs();
   std::map<Vertex*,int, std::less<Vertex*> > mark;
     
   std::list<Vertex_handle> L;
     
   dfs(delau,mark,v,C,L);
     
   if (new_v)
   { L.pop_front();   //first one was inserted in range_search ...
     delau.remove(v);
   }
     
   typename std::list<Vertex_handle>::const_iterator iter = L.begin();
   for(;iter != L.end() ;iter++){ *res= *iter; res++; }
   return res;        
}
   

template<class Dt, class OutputIterator>
OutputIterator range_search(Dt& delau, const typename Dt::Point& a, const typename Dt::Point& b, 
                            const typename Dt::Point& c,OutputIterator res)
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


template<class Dt,class OutputIterator>
OutputIterator range_search(Dt& delau, const typename Dt::Point& a1, const typename Dt::Point& b1, 
                            const typename Dt::Point& c1,const typename Dt::Point& d1,OutputIterator res)
// a1 upper left, b1 lower left , c1 lower right
{
  typedef typename Dt::Geom_traits                    Gt;
  typedef typename Dt::Point                          Point;
  typedef typename Gt::Circle_2                       Circle;
  typedef typename Dt::Vertex_handle                  Vertex_handle;
  typedef typename Gt::Orientation_2                  Orientation_2;
  typedef typename Gt::Construct_circle_2             Construct_circle_2;
  
  Point a=a1,b=b1,c=c1,d=d1;
    
  if (Orientation_2()(a,b,c) == RIGHTTURN) 
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
    if ( Orientation_2()(a,b,p) == RIGHTTURN || Orientation_2()(b,c,p) == RIGHTTURN ||
         Orientation_2()(c,d,p) == RIGHTTURN || Orientation_2()(d,a,p) == RIGHTTURN )  { }
    else { *res = *it; res++; }
  }
  return res;     
}



CGAL_END_NAMESPACE

#endif
