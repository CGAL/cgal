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
// file          : include/CGAL/Point_set_2_tb.h
// package       : Point_set_2_tb (0.9)
// maintainer    : Matthias Baesken <baesken@informatik.uni-trier.de>
// revision      : 0.9
// revision_date : 27 March 2001 
// author(s)     : Matthias Baesken
//
// coordinator   : Matthias Baesken, Trier  (<baesken@informatik.uni-trier.de>)
// ======================================================================

#ifndef CGAL_POINT_SET_2_TB_H
#define CGAL_POINT_SET_2_TB_H

#include <CGAL/basic.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <list>
#include <queue>
#include <map>
#include <stack>

#include <cmath>

#include <CGAL/point_set_traits_2.h>


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



template<class Str, class Gt, class Tds>
class  Point_set_2_tb : public  Delaunay_triangulation_2<Gt,Tds>
{

public:  
  typedef Gt Geom_traits;
  
  typedef typename Geom_traits::Rep                         Rep;
  typedef typename Geom_traits::Point_2                     Point;
  typedef typename Geom_traits::Segment_2                   Segment;
  
  typedef typename Geom_traits::Orientation_2               Orientation_2;
  typedef typename Geom_traits::Side_of_oriented_circle_2   Side_of_oriented_circle_2;
  typedef typename Rep::Construct_circle_2                  Construct_circle_2; 
  typedef typename Rep::Compute_squared_distance_2          Compute_squared_distance_2;
  typedef typename Rep::FT                                  Numb_type;  // field number type ...
                        
  
  typedef Triangulation_2<Gt,Tds>                           Triangulation;
  typedef typename Triangulation::Locate_type               Locate_type;
  typedef typename Triangulation::Face_handle               Face_handle;
  typedef typename Triangulation::Vertex_handle             Vertex_handle;
  typedef typename Triangulation::Edge                      Edge;
  typedef typename Triangulation::Vertex                    Vertex;
  typedef typename Triangulation::Face                      Face;
  typedef typename Triangulation::Edge_circulator           Edge_circulator;
  typedef typename Triangulation::Finite_edges_iterator     Finite_edges_iterator;
  typedef typename Triangulation::Vertex_iterator           Vertex_iterator;
  typedef typename Triangulation::Edge_iterator             Edge_iterator;

  // now the types from Str
  typedef typename Str::Circle   Circle;

  typedef typename Str::Compare_dist_2              Comparedist;             
  typedef typename Str::Circle_bounded_side_2       Circleptori;
  typedef typename Str::Circle_center_2             Circlecenter;     
  
   Str                           traits;
   Comparedist                   tr_comparedist;
   Orientation_2                 tr_orientation;  
   Side_of_oriented_circle_2     tr_so_circle;    
   Compute_squared_distance_2    tr_sqrdist;      
   Circleptori                   tr_circleptori;
   Circlecenter                  tr_circlecenter; 
   
   //constructions...
   Construct_circle_2            tr_createcircle_3p;

   Point_set_2_tb()
   { 
     init_vertex_marks();
   }


   Point_set_2_tb(const std::list<Point>& S)
   { 
     init_vertex_marks();
     init(S);
   }
   
   template<class InputIterator>
   Point_set_2_tb(InputIterator first, InputIterator last)
   { 
     init_vertex_marks();    
     init(first,last);
   }
   
   ~Point_set_2_tb() {}


   void  init(const std::list<Point>& L0)
   { 
    clear();
    insert(L0.begin(), L0.end());
    init_vertex_marks();
  }    

   
   template<class InputIterator>
   void init(InputIterator first, InputIterator last)
   { 
    clear();
    insert(first,last);
    init_vertex_marks();
  } 


   std::list<Point>  points() const
   { std::list <Point> L;
     Vertex_iterator vit = finite_vertices_begin();
     for(;vit != finite_vertices_end(); vit++) {
        Vertex v = *vit;
	L.push_back(v.point());
     }
     return L;
   }

   template<class OutputIterator>
   OutputIterator points(OutputIterator out)
   {    
     Vertex_iterator vit = finite_vertices_begin();
     for (; vit != finite_vertices_end(); vit++) { *out= (*vit).point(); out++; }       
     return out;
   }
   
   template<class OutputIterator>
   OutputIterator segments(OutputIterator out)
   {    
     Edge_iterator eit = finite_edges_begin();
     for (; eit != finite_edges_end(); eit++) { 
       *out= segment(eit); 
       out++;
     }       
     return out;
   }     
   
   template<class OutputIterator>
   OutputIterator vertices(OutputIterator res)
   // return vertex handles ...
   {    
    Vertex_iterator vit = finite_vertices_begin();
    for (; vit != finite_vertices_end(); vit++) { *res= (*vit).handle(); res++; }  
    return res;   
   }
   
   template<class OutputIterator>
   OutputIterator edges(OutputIterator out)
   {    
     Edge_iterator eit = finite_edges_begin();
     for (; eit != finite_edges_end(); eit++) { 
       *out = *eit;
     }   
     return out;
   }        
   
   // -----------------------------------------------------
   // d_face_cycle_succ and d_face_cycle_pred removed ...
   // -----------------------------------------------------
   

   bool   empty() { return number_of_vertices() == 0; }
   
   bool   is_empty() { return number_of_vertices() == 0; }

   // -------------------------------------------------------------------------
   // locate removed; use the locate - operations inherited from triangulation
   // -------------------------------------------------------------------------

   // Vertex -> Vertex_handle ...
   Vertex_handle lookup(Point p) const
   { 
     if (number_of_vertices() == 0) return NULL;   
     
     // locate ...
     Locate_type lt;
     int li;
     Face_handle fh = locate(p,lt,li);
     
     if (lt == VERTEX){
        Face f = *fh;
	return f.vertex(li);
     }
     else return NULL;
   }


   void del(Vertex_handle v)
   {  
     remove(v);
   }


   void del(Point p)
   { Vertex_handle v = lookup(p);
     if ( v != NULL ) del(v);
   }


   Vertex_handle  nearest_neighbor(Point p)
    {
     if (number_of_vertices() == 0) return NULL;
     return nearest_vertex(p);
   }
     

   Vertex_handle  nearest_neighbor(Vertex_handle v) const
   {
     if (number_of_vertices() <= 1) return NULL;    
     Point p = v->point();
     
     Vertex_circulator vc = incident_vertices(v);
     Vertex_circulator start =vc;
     
     Vertex_handle min_v = (*vc).handle();
     if (is_infinite(min_v)){
       vc++;
       min_v = (*vc).handle();
     }
     
     Vertex_handle act;
     
     // go through the vertices ...
     do {
       act = vc->handle();
 
       if (! is_infinite(act)) {
        if ( tr_comparedist(p,act->point(), min_v->point()) == SMALLER ) min_v = act;
       }   
           
       vc++;
     } while (vc != start);     
   
     return min_v;
   }
   
  template<class OutputIterator>
  OutputIterator   nearest_neighbors(Point p, int k, OutputIterator res)
  {
   int n = number_of_vertices();

   if ( k <= 0 ) return res;
   if ( n <= k ) { // return all finite vertices ...
     return vertices(res);
   }

   // insert p, if nesessary
   
    Vertex_handle vh = lookup(p);
    bool old_node = true;
    
    // we have to add a new vertex ...
    if (vh == NULL){
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
    
    std::list<Vertex_handle>::iterator it = res_list.begin();
    
    for (; it != res_list.end(); it++) { *res= *it; res++; }

    return res;  
  }
   
  template<class OutputIterator>  
  OutputIterator  nearest_neighbors(Vertex_handle v, int k,OutputIterator res)
  {  
   int n = number_of_vertices();

   if ( k <= 0 ) return res;
   if ( n <= k ) { // return all (finite) vertices ...
     return vertices(res);
   }
   
   std::list<Vertex_handle> res_list;
   nearest_neighbors_list(v, k, res_list); 
   
   std::list<Vertex_handle>::iterator it = res_list.begin();
    
   for (; it != res_list.end(); it++) { *res= *it; res++; }

   return res;     
  }

  void nearest_neighbors_list(Vertex_handle v, int k, std::list<Vertex_handle>& res) 
  {  
     int n = number_of_vertices();
   
     if ( k <= 0 ) return;
     if ( n <= k ) { vertices(std::back_inserter(res)); return; }
     
     Point p = v->point();
     
     // "unmark" the vertices ...
     init_dfs();


     std::map<Vertex*,Numb_type, std::less<Vertex*> > priority_number; // here we save the priorities ...
     compare_vertices<Vertex*,Numb_type> comp(& priority_number);      // comparison object ...
     std::priority_queue<Vertex*, std::vector<Vertex*>, CGAL::compare_vertices<Vertex*,Numb_type> > PQ(comp);

     priority_number[v.ptr()] = 0;
     PQ.push(v.ptr());
     
     mark_vertex(v);
      
     while ( k > 0 )
     { 
       // find minimum from PQ ...
       Vertex* w = PQ.top(); PQ.pop();
       priority_number.erase(w); // and clean entry in priority map
   
       res.push_back(w->handle()); k--; 

       // get the incident vertices of w ...
       Vertex_circulator vc = incident_vertices(w->handle());
       Vertex_circulator start =vc;
       Vertex* act;
     
       do {
         act = &(*vc);
	 
         if ( (!is_marked(act)) && (! is_infinite(act->handle())) )
         { 
             priority_number[act] = tr_sqrdist(p,act->point());
             PQ.push(act);	      
             mark_vertex(act->handle());
         }	   
	             
         vc++;
       } while (vc != start);   
        
     }
  } 
   
   
  // dfs
  // for marking nodes in search procedures
  int cur_mark;
   
  std::map<Vertex*,int, std::less<Vertex*> > mark;  
  typedef typename std::map<Vertex*,int, std::less<Vertex*> >::iterator map_iterator;

  
  void init_vertex_marks()
  {
     cur_mark = 0;
     mark.erase(mark.begin(), mark.end());
  }
  
  void init_dfs()
  {
     cur_mark++; 
     if (cur_mark == INT_MAX) init_vertex_marks();
  }
  
  void mark_vertex(Vertex_handle vh)
  // mark vh as visited ...
  {
    Vertex* v = vh.ptr();
    mark[v] = cur_mark;
  }
  
  bool is_marked(Vertex_handle vh)
  {
    Vertex* v = vh.ptr();
    
    map_iterator mit = mark.find(v);   
    if (mit == mark.end()) return false;

    return (mark[v] == cur_mark);
  }
  
  void dfs(Vertex_handle v,const Circle& C, std::list<Vertex_handle>& L)
  {
    L.push_back(v);
    mark_vertex(v);
    
    // get incident vertices of v ...
    Vertex_circulator vc = incident_vertices(v);
    Vertex_circulator start =vc;
     
    Vertex_handle act;
     
    // go through the vertices ...
    do {
       act = vc->handle();
 
       if (! is_infinite(act)) {
        if (!is_marked(act) && ! (tr_circleptori(C,act->point())==ON_UNBOUNDED_SIDE) ) 
           dfs(act,C,L);       
       }             
       vc++;
    } while (vc != start);     
  }


   template<class OutputIterator>
   OutputIterator range_search(const Circle& C, OutputIterator res)
   { 
     if (number_of_vertices() == 0) return res;
     if (number_of_vertices() == 1) 
     { 
       // get the one vertex ...
       Vertex_iterator vit = finite_vertices_begin();
       Vertex v = (*vit);
       Point p = v.point();
       
       if (! (tr_circleptori(C, p) == ON_UNBOUNDED_SIDE)){
        *res= v.handle(); res++;
       }
       return res;
     }  
     
     // normal case ...
     Point p = tr_circlecenter(C);
     Vertex_handle v = lookup(p);  
     bool new_v = false;     

     if ( v == NULL )
     { 
       new_v = true;
       v = insert(p); 
     }
     
     init_dfs();
     
     std::list<Vertex_handle> L;
     dfs(v,C,L);
     
     if (new_v)
     { L.pop_front();   //first one was inserted in range_search ...
       remove(v);
     }
     
     std::list<Vertex_handle>::const_iterator iter = L.begin();
     for(;iter != L.end() ;iter++){ *res= *iter; res++; }
     return res;        
   }
   

   template<class OutputIterator>
   OutputIterator range_search(const Point& a, const Point& b, const Point& c,OutputIterator res)
   { int orient = (int)(tr_orientation(a,b,c));
     Circle C = tr_createcircle_3p(a,b,c);
     std::list<Vertex_handle> L;
     range_search(C,std::back_inserter(L));
      
     std::list<Vertex_handle>::const_iterator it = L.begin();
      
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
   
     if (tr_orientation(a,b,c) == RIGHTTURN) 
     { Point tmp = b;
       b = d;
       d = tmp;
      }
   
     Circle C = tr_createcircle_3p(a,b,c);
     
     std::list<Vertex_handle> L;
     range_search(C,std::back_inserter(L));
     std::list<Vertex_handle>::const_iterator it = L.begin();     

     for(;it != L.end();it++)
     { Point p = (*it)->point();
       if ( tr_orientation(a,b,p) == RIGHTTURN || tr_orientation(b,c,p) == RIGHTTURN ||
            tr_orientation(c,d,p) == RIGHTTURN || tr_orientation(d,a,p) == RIGHTTURN )  { }
        else { *res = *it; res++; }
     }
     return res;     
   }
 
   bool     IS_NON_DIAGRAM_DART(Edge e) const
   { 
    Face_handle f1 = e.first;
    Face_handle f2 = f1->neighbor(e.second);
    int i = e.second;

    // flip test
    Point a = f1->vertex(cw(i))->point();
    Point b = f1->vertex(i)->point();
    Point c = f1->vertex(ccw(i))->point();
    Point d = f1->mirror_vertex(i)->point();

    if ((tr_orientation(a,b,c)==LEFTTURN) && (tr_orientation(b,c,d)==LEFTTURN) && 
        (tr_orientation(c,d,a)==LEFTTURN) && (tr_orientation(d,a,c)==LEFTTURN) && 
        (tr_so_circle(a,b,c,d)==ON_ORIENTED_BOUNDARY) )
        return true;
    return false;
  }
  
  bool is_non_diagram_edge(Edge e) const
  {
    return IS_NON_DIAGRAM_DART(e); 
  }

   Edge get_cur_dart()  const { return cur_dart; }
   void set_cur_dart(Edge e)  { cur_dart = e;    }
   void set_hull_dart(Edge e) { hull_dart = e;   }

   Point pos(Vertex_handle v) const
   { return v->point(); }
    
   Point pos_source(Edge e) const
   { return segment(e).source(); }
    
   Point pos_target(Edge e) const
   { return segment(e).target(); }
    
   Segment seg(Edge e) const           
   { return segment(e); }

};


CGAL_END_NAMESPACE


#endif

