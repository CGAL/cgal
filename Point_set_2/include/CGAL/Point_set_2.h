// Copyright (c) 1999
// Max-Planck-Institute Saarbruecken (Germany). All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Matthias Baesken

#ifndef CGAL_POINT_SET_2_H
#define CGAL_POINT_SET_2_H

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
#include <climits>

namespace CGAL {

template<class Gt, class Tds = Triangulation_data_structure_2 <Triangulation_vertex_base_2<Gt> > >
class  Point_set_2 : public  Delaunay_triangulation_2<Gt,Tds>
{

public:
  typedef Gt Geom_traits;

  typedef typename Geom_traits::Point_2                     Point;
  typedef typename Geom_traits::Segment_2                   Segment;

  typedef typename Geom_traits::Circle_2                    Circle;

  typedef typename Geom_traits::Orientation_2               Orientation_2;
  typedef typename Geom_traits::Side_of_oriented_circle_2   Side_of_oriented_circle_2;
  typedef typename Geom_traits::Construct_circle_2          Construct_circle_2;
  typedef typename Geom_traits::Compute_squared_distance_2  Compute_squared_distance_2;
  typedef typename Geom_traits::FT                          Numb_type;  // field number type ...


  typedef Delaunay_triangulation_2<Gt,Tds>                  Triangulation;
  typedef typename Triangulation::size_type                 size_type;
  typedef typename Triangulation::Locate_type               Locate_type;
  typedef typename Triangulation::Face_handle               Face_handle;
  typedef typename Triangulation::Vertex_handle             Vertex_handle;
  typedef typename Triangulation::Edge                      Edge;
  typedef typename Triangulation::Vertex                    Vertex;
  typedef typename Triangulation::Face                      Face;
  typedef typename Triangulation::Edge_circulator           Edge_circulator;
  typedef typename Triangulation::Finite_edges_iterator     Finite_edges_iterator;
  typedef typename Triangulation::Vertex_iterator           Vertex_iterator;
  typedef typename Triangulation::Vertex_circulator         Vertex_circulator;
  typedef typename Triangulation::Edge_iterator             Edge_iterator;

  typedef typename Geom_traits::Bounded_side_2              Circleptori;
  typedef typename Geom_traits::Compare_distance_2          Comparedist;
  typedef typename Geom_traits::Construct_center_2          Circlecenter;

  typedef Unique_hash_map<Vertex_handle, Numb_type>         MAP_TYPE;
  typedef Delaunay_triangulation_2<Gt,Tds>                  Base;

  using Base::finite_vertices_begin;
  using Base::finite_vertices_end;
  using Base::number_of_vertices;
  using Base::VERTEX;
  using Base::insert;
  using Base::remove;
  using Base::locate;
  using Base::is_infinite;
  using Base::nearest_vertex;
  using Base::incident_vertices;

   Comparedist                   tr_comparedist;
   Orientation_2                 tr_orientation;
   Side_of_oriented_circle_2     tr_so_circle;
   Compute_squared_distance_2    tr_sqrdist;
   Circleptori                   tr_circleptori;
   Circlecenter                  tr_circlecenter;

   //constructions...
   Construct_circle_2            tr_createcircle_3p;

   Point_set_2()
   {
     init_vertex_marks();
   }

   template<class InputIterator>
   Point_set_2(InputIterator first, InputIterator last)
   {
     init_vertex_marks();
     insert(first,last);
   }

   ~Point_set_2() {}

   template<class OutputIterator>
   OutputIterator vertices(OutputIterator res)
   // return vertex handles ...
   {
    Vertex_iterator vit = finite_vertices_begin();
    for (; vit != finite_vertices_end(); vit++) { *res= vit; res++; }
    return res;
   }


   Vertex_handle lookup(Point p) const
   {
     if (number_of_vertices() == 0) return nullptr;

     // locate ...
     Locate_type lt;
     int li;
     Face_handle fh = locate(p,lt,li);

     if (lt == VERTEX){
        Face f = *fh;
        return f.vertex(li);
     }
     else return nullptr;
   }


   Vertex_handle  nearest_neighbor(Point p)
    {
     if (number_of_vertices() == 0) return nullptr;
     return nearest_vertex(p);
   }


   Vertex_handle  nearest_neighbor(Vertex_handle v) const
   {
     if (number_of_vertices() <= 1) return nullptr;
     Point p = v->point();

     Vertex_circulator vc = incident_vertices(v);
     Vertex_circulator start =vc;

     Vertex_handle min_v = vc;
     if (is_infinite(min_v)){
       vc++;
       min_v = vc;
     }

     Vertex_handle act;

     // go through the vertices ...
     do {
       act = vc;

       if (! is_infinite(act)) {
        if ( tr_comparedist(p,act->point(), min_v->point()) == SMALLER ) {
          min_v = act;
        }
       }

       vc++;
     } while (vc != start);

     return min_v;
   }

  template<class OutputIterator>
  OutputIterator   nearest_neighbors(Point p, size_type k, OutputIterator res)
  {
   size_type n = number_of_vertices();

   if ( k <= 0 ) return res;
   if ( n <= k ) { // return all finite vertices ...
     return vertices(res);
   }

   // insert p, if nesessary

    Vertex_handle vh = lookup(p);
    bool old_node = true;

    // we have to add a new vertex ...
    if (vh == nullptr){
      vh = insert(p);
      old_node = false;
      k++;
    }

    std::list<Vertex_handle> res_list;
    nearest_neighbors_list(vh, k, res_list);

    if ( !old_node )
    {
     res_list.pop_front();
     remove(vh);
    }

    typename std::list<Vertex_handle>::iterator it = res_list.begin();

    for (; it != res_list.end(); it++) { *res= *it; res++; }

    return res;
  }

  template<class OutputIterator>
  OutputIterator  nearest_neighbors(Vertex_handle v, size_type k,OutputIterator res)
  {
   size_type n = number_of_vertices();

   if ( k <= 0 ) return res;
   if ( n <= k ) { // return all (finite) vertices ...
     return vertices(res);
   }

   std::list<Vertex_handle> res_list;
   nearest_neighbors_list(v, k, res_list);

   typename std::list<Vertex_handle>::iterator it = res_list.begin();

   for (; it != res_list.end(); it++) { *res= *it; res++; }

   return res;
  }


  void nearest_neighbors_list(Vertex_handle v,
                              size_type k,
                              std::list<Vertex_handle>& res)
  {
     size_type n = number_of_vertices();

     if ( k <= 0 ) return;
     if ( n <= k ) { vertices(std::back_inserter(res)); return; }

     Point p = v->point();

     // "unmark" the vertices ...
     init_search();

     MAP_TYPE                                        priority_number;              // here we save the priorities ...
     internal::compare_vertices<Vertex_handle,Numb_type,MAP_TYPE>
       comp(& priority_number);      // comparison object ...
     std::priority_queue<Vertex_handle, std::vector<Vertex_handle>, internal::compare_vertices<Vertex_handle,Numb_type,MAP_TYPE> > PQ(comp);

     priority_number[v] = 0;
     PQ.push(v);

     mark_vertex(v);

     while ( k > 0 )
     {
       // find minimum from PQ ...
       Vertex_handle w = PQ.top();
       PQ.pop();

       res.push_back(w);
       k--;

       // get the incident vertices of w ...
       Vertex_circulator vc = incident_vertices(w);
       Vertex_circulator start =vc;
       Vertex_handle act;

       do {
         act = vc;

         if ( (!is_marked(act)) && (! is_infinite(act)) )
         {
             priority_number[act] = tr_sqrdist(p,act->point());
             PQ.push(act);
             mark_vertex(act);
         }

         vc++;
       } while (vc != start);

     }
  }


  // for marking nodes in search procedures
  size_type cur_mark;

  Unique_hash_map<Vertex_handle, size_type>  mark;

  void init_vertex_marks()
  {
     cur_mark = 0;
     mark.clear();
  }

  void init_search()
  {
     cur_mark++;
     if (cur_mark == (std::numeric_limits<size_type>::max)()) init_vertex_marks();
  }

  void mark_vertex(Vertex_handle vh)
  // mark vh as visited ...
  {
    mark[vh] = cur_mark;
  }

  bool is_marked(Vertex_handle vh)
  {
    if (! mark.is_defined(vh)) return false;

    return (mark[vh] == cur_mark);
  }

  void search(Vertex_handle v,const Circle& C, std::list<Vertex_handle>& L)
  {
    std::stack<Vertex_handle> todo;
    todo.push(v);

    while (!todo.empty())
    {
      Vertex_handle current = todo.top();
      todo.pop();

      if (is_marked(current))
        continue;

      L.push_back(current);
      mark_vertex(current);

      // get incident vertices of v ...
      Vertex_circulator vc = incident_vertices(current);
      Vertex_circulator start =vc;
      Vertex_handle act;
      // go through the vertices ...
      do {
        act = vc;

        if (! is_infinite(act)) {
          if (!is_marked(act) && ! (tr_circleptori(C,act->point())==ON_UNBOUNDED_SIDE) )
            todo.push(act);
        }
        vc++;
      } while (vc != start);
    }
  }


   template<class OutputIterator>
   OutputIterator range_search(const Circle& C, OutputIterator res)
   {
     if (number_of_vertices() == 0) return res;
     if (number_of_vertices() == 1)
     {
       // get the one vertex ...
       Vertex_iterator vit = finite_vertices_begin();
       Point p = vit->point();

       if (! (tr_circleptori(C, p) == ON_UNBOUNDED_SIDE)){
        *res= vit; res++;
       }
       return res;
     }

     // normal case ...
     Point p = tr_circlecenter(C);
     Vertex_handle v = lookup(p);
     bool new_v = false;

     if ( v == nullptr )
     {
       new_v = true;
       v = insert(p);
     }

     init_search();

     std::list<Vertex_handle> L;
     search(v,C,L);

     if (new_v)
     { L.pop_front();   //first one was inserted in range_search ...
       remove(v);
     }

     typename std::list<Vertex_handle>::const_iterator iter = L.begin();
     for(;iter != L.end() ;iter++){ *res= *iter; res++; }
     return res;
   }


   template<class OutputIterator>
   OutputIterator range_search(const Point& a, const Point& b, const Point& c,OutputIterator res)
   { int orient = (int)(tr_orientation(a,b,c));
     Circle C = tr_createcircle_3p(a,b,c);
     std::list<Vertex_handle> L;
     range_search(C,std::back_inserter(L));

     typename std::list<Vertex_handle>::const_iterator it = L.begin();

     for(;it != L.end();it++)
     { Point p = (*it)->point();
       if ( ((int)(tr_orientation(a,b,p))) == - orient ||
            ((int)(tr_orientation(b,c,p))) == - orient ||
            ((int)(tr_orientation(c,a,p))) == - orient ) { }
        else { *res = *it; res++; }
     }
     return res;
   }


   template<class OutputIterator>
   OutputIterator range_search(const Point& a1, const Point& b1, const Point& c1,const Point&
   d1,OutputIterator res)
   // a1 upper left, b1 lower left , c1 lower right
   {
     //Point b(c.xcoord(),a.ycoord());
     //Point d(a.xcoord(),c.ycoord());
     Point a=a1,b=b1,c=c1,d=d1;

     if (tr_orientation(a,b,c) == RIGHT_TURN)
     { Point tmp = b;
       b = d;
       d = tmp;
      }

     Circle C = tr_createcircle_3p(a,b,c);

     std::list<Vertex_handle> L;
     range_search(C,std::back_inserter(L));
     typename std::list<Vertex_handle>::const_iterator it = L.begin();

     for(;it != L.end();it++)
     { Point p = (*it)->point();
       if ( tr_orientation(a,b,p) == RIGHT_TURN || tr_orientation(b,c,p) == RIGHT_TURN ||
            tr_orientation(c,d,p) == RIGHT_TURN || tr_orientation(d,a,p) == RIGHT_TURN )  { }
        else { *res = *it; res++; }
     }
     return res;
   }

};


} //namespace CGAL


#endif
