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
// release_date  : 2000, September 14
//
// file          : include/CGAL/Point_set_2.h
// package       : Point_set_2 (1.2.4)
// maintainer    : Matthias Baesken <baesken@informatik.uni-trier.de>
// revision      : 1.2.4
// revision_date : 14 September 2000 
// author(s)     : Kurt Mehlhorn, Stefan Naeher, Matthias Baesken
//
// coordinator   : Matthias Baesken, Halle  (<baesken@informatik.uni-trier.de>)
// ======================================================================

#ifndef CGAL_POINT_SET_2_H
#define CGAL_POINT_SET_2_H

#if !defined(LEDA_ROOT_INCL_ID)
#define LEDA_ROOT_INCL_ID 400900
#include <LEDA/REDEFINE_NAMES.h>
#endif

#if !defined(LEDA_STL_ITERATORS)
#define LEDA_STL_ITERATORS
#endif


#ifndef POINT_SET_LEDA_ONLY
#include <CGAL/config.h>
#endif

#ifdef POINT_SET_LEDA_ONLY
#  define CGAL_BEGIN_NAMESPACE namespace CGAL {
#  define CGAL_END_NAMESPACE }
#endif

#include <list>
#include <LEDA/basic.h>

#ifndef POINT_SET_LEDA_ONLY

#if (__MSC_VER > 1100)
#if defined(__ALL_MEMBERS_INSTANT__)
#undef __ALL_MEMBERS_INSTANT__
#endif
#endif

#endif

#include <LEDA/p_queue.h>
#include <math.h>

enum cgal_ps_delaunay_edge_info{ CGAL_DIAGRAM_EDGE = 0, CGAL_DIAGRAM_DART = 0,
                         CGAL_NON_DIAGRAM_EDGE = 1, CGAL_NON_DIAGRAM_DART = 1,
                         CGAL_HULL_EDGE = 2, CGAL_HULL_DART = 2
                       };


#include <LEDA/graph.h>

// Traits classes for the 2d point set

#if defined(LEDA_PREFIX)
#include <LEDA/UNDEFINE_NAMES.h>
#endif

#include <CGAL/point_set_traits_2.h>

#if defined(LEDA_PREFIX)
#include <LEDA/REDEFINE_NAMES.h>
#endif


#include <CGAL/point_set_leda_traits_2.h>

#if defined(LEDA_PREFIX)
#include <LEDA/UNDEFINE_NAMES.h>
#endif


// from graph_alg.h
extern leda_list<leda_edge> MIN_SPANNING_TREE(const leda_graph&,int(*cmp)(const leda_edge&,const leda_edge&));


CGAL_BEGIN_NAMESPACE


template<class TR>
class __exportC Point_set_2 : public  GRAPH<__typename TR::Point,int>
{

public:  
  typedef typename TR::FT       Numb_type;
  typedef typename TR::Point    Point;
  typedef typename TR::Circle   Circle;
  typedef typename TR::Segment  Segment;
  typedef typename TR::Line     Line;

  typedef leda_edge Edge;
  typedef leda_node Vertex;

  //functionality on these types...
  typedef typename TR::Compare_xy_2                Comparepoints; 
  typedef typename TR::Compare_dist_2              Comparedist;   
  typedef typename TR::Orientation                 Orientation_2;   
  typedef typename TR::Side_of_oriented_circle_2   Sideofcircle;
  typedef typename TR::Side_of_halfspace_2         Sideofhalfspace;
  typedef typename TR::Segment_has_on_2            Segmentcontains;
  typedef typename TR::Squared_distance            Sqrdist;          
  typedef typename TR::Squared_distance_to_line    Linesqrdist;      
  typedef typename TR::Circle_bounded_side_2       Circleptori;
  typedef typename TR::Circle_center_2             Circlecenter;     
  
  //constructors...
  typedef typename TR::Construct_circle_2          Createcircle_3p;  
  typedef typename TR::Construct_segment_2         Createsegment_2p;  
  typedef typename TR::Construct_line_2            Createline_2p;     
  

private:
   TR               traits;

   Comparepoints    tr_comparepoints;
   Comparedist      tr_comparedist;
   Orientation_2    tr_orientation;  
   Sideofcircle     tr_so_circle;   
   Sideofhalfspace  tr_so_hp;        
   Segmentcontains  tr_seg_contains; 
   Sqrdist          tr_sqrdist;      
   Linesqrdist      tr_linesqrdist; 
   Circleptori      tr_circleptori;
   Circlecenter     tr_circlecenter; 
   
   //constructors...
   Createcircle_3p  tr_createcircle_3p;
   Createsegment_2p tr_createsegment_2p;
   Createline_2p    tr_createline_2p;


   leda_edge cur_dart;
   leda_edge hull_dart;

   bool check; // functions are checked if true

   // for marking nodes in search procedures
   int cur_mark;
   leda_node_map<int> mark;

   void init_node_marks()
   { mark.init(*this,-1);
     cur_mark = 0;
   }

   void mark_node(leda_node v) const 
   { ((leda_node_map<int>&)mark)[v] = cur_mark; }

   void unmark_node(leda_node v) const 
   { ((leda_node_map<int>&)mark)[v] = cur_mark - 1; }

   bool is_marked(leda_node v)  const { return mark[v] == cur_mark; }

   void unmark_all_nodes() const 
   { ((int&)cur_mark)++; 
     if ( cur_mark == MAXINT)
     ((Point_set_2*)this) -> init_node_marks();  
                                   //cast away constness
   }
   
   void mark_edge(leda_edge e, cgal_ps_delaunay_edge_info k) 
   { assign(e,k); }

   leda_node  new_node(Point p) 
   { leda_node v = GRAPH<Point,int>::new_node(p); 
     mark[v] = -1;
     return v;
   }

   leda_edge  new_edge(leda_node v, leda_node w) 
   { return GRAPH<Point,int>::new_edge(v,w,0); }

   leda_edge  new_edge(leda_edge e, leda_node w) 
   { return GRAPH<Point,int>::new_edge(e,w,0,0); }

   void  del_node(leda_node v)  { GRAPH<Point,int>::del_node(v); }
   void  del_edge(leda_edge e)  { GRAPH<Point,int>::del_edge(e); }


   void init_hull();
   
   void make_delaunay(leda_list<leda_edge>& S);
   
   void make_delaunay();

   void check_locate(leda_edge answer,const Point& p) const;

   void dfs(leda_node s, const Circle& C, leda_list<leda_node>& L) const;

   void dfs(leda_node s, const Point& pv, const Point& p, leda_list<leda_node>& L) const;

public:

   Point_set_2()
   { cur_dart = nil; 
     hull_dart = nil;
     init_node_marks();
     check = false;
   }

   Point_set_2(const leda_list<Point>& S)
   { cur_dart = nil; 
     hull_dart = nil;
     init_node_marks();
     check = false;
     
     init(S);
   }

   Point_set_2(const std::list<Point>& S)
   { cur_dart = nil; 
     hull_dart = nil;
     init_node_marks();
     check = false;
     
     init(S);
   }
   
   template<class InputIterator>
   Point_set_2(InputIterator first, InputIterator last)
   { cur_dart = nil; 
     hull_dart = nil;
     init_node_marks();
     check = false;
     
     init(first,last);
   }

   Point_set_2(const GRAPH<Point,int>& G) :GRAPH<Point,int>(G)
   { 
     leda_edge e;
     forall_edges(e,*this) mark_edge(e,CGAL_DIAGRAM_DART);
     init_hull(); 
     make_delaunay(); 
     if (number_of_edges() > 0) 
       cur_dart = reversal(hull_dart); 
     else 
       cur_dart = nil;
     init_node_marks();
     check = false;
   }

   Point_set_2(const Point_set_2& T):GRAPH<Point,int>(T)
   { 
     init_hull(); 
     if (number_of_edges() > 0) 
        cur_dart = reversal(hull_dart); 
     else 
        cur_dart = nil;
     init_node_marks();
     check = false;
   }
   
   Point_set_2& operator=(const Point_set_2& T)
   { 
     GRAPH<Point,int>::operator=(T); 
     init_hull(); 
     if (number_of_edges() > 0) 
       cur_dart = reversal(hull_dart); 
     else 
       cur_dart = nil;
     init_node_marks();
     check = false;
     return *this; 
   } 
   
   ~Point_set_2() {}

   void  init(const leda_list<Point>& L0);
 
   void  init(const std::list<Point>& L0);
   
   template<class InputIterator>
   void init(InputIterator first, InputIterator last)
   { 
    // construct a triangulation for the points in L0

    clear();

    if (first==last) return;

    leda_list<Point> L;
    InputIterator it=first;
    
    for(;it!=last;++it){
      L.append(*it);
    }
    
    L.sort(tr_comparepoints);  

    // initialize graph with a single edge starting at the first point

    Point last_p = L.pop();          // last visited point
    leda_node  last_v = new_node(last_p); // last inserted node

    while (!L.empty() && last_p == L.head()) L.pop();

    if (L.empty()) return;

    last_p = L.pop();
    leda_node v = new_node(last_p);
    leda_edge x = new_edge(last_v,v);
    leda_edge y = new_edge(v,last_v);
    set_reversal(x,y);
    last_v = v;


    // scan remaining points

    Point p=last_p; // changed ....
    forall(p,L) 
    { if (p == last_p) continue; 

      leda_edge e =  last_adj_edge(last_v);

      last_v = new_node(p);
      last_p = p;

      // walk up to upper tangent
      do e = face_cycle_pred(e); while (orientation(e,p) > 0);

      // walk down to lower tangent and triangulate
      do { e = face_cycle_succ(e);
         leda_edge x = new_edge(e,last_v);
         leda_edge y = new_edge(last_v,source(e));
         set_reversal(x,y);
         } while (orientation(e,p) > 0);
     }

    // mark edges of convex hull

    hull_dart = last_edge();
  
    cur_dart = reversal(hull_dart);
 
    leda_edge e0 = hull_dart;
    leda_edge e  = e0;
    do { mark_edge(e,CGAL_HULL_DART);
       e = face_cycle_succ(e);
       } while (e != e0); 

    make_delaunay();

    init_node_marks();
  } 

   leda_list<Point>  points() const;

   template<class OutputIterator>
   OutputIterator points(OutputIterator out)
   {    
     leda_node v;
     forall_nodes(v,*this) { *out= pos(v); out++; } 
     return out;
   }
   
   template<class OutputIterator>
   OutputIterator segments(OutputIterator out)
   {    
     leda_edge e;
     forall_edges(e,*this) { *out= seg(e); out++; } 
     return out;
   }     
   
   template<class OutputIterator>
   OutputIterator vertices(OutputIterator out)
   {    
     leda_node v;
     forall_nodes(v,*this) { *out= v; out++; } 
     return out;
   }
   
   template<class OutputIterator>
   OutputIterator edges(OutputIterator out)
   {    
     leda_edge e;
     forall_edges(e,*this) { *out= e; out++; } 
     return out;
   }        
  
   Edge   d_face_cycle_succ(Edge e) const
   { e = reversal(e);
     do e = cyclic_adj_pred(e); 
     while (!is_diagram_dart(e));
     return e;
   }

   Edge   d_face_cycle_pred(Edge e) const
   { do e = cyclic_adj_succ(e); 
       while (!is_diagram_dart(e));
     return reversal(e);
   }

   bool   empty() { return number_of_nodes() == 0; }
   
   bool   is_empty() { return number_of_nodes() == 0; }

   void   clear() { GRAPH<Point,int>::clear(); cur_dart = hull_dart = nil; }

   Edge   locate(Point p) const;

   Vertex   lookup(Point p) const;

   Vertex   insert(Point p);

   void del(Vertex v);

   void del(Point p);

   Vertex   nearest_neighbor(Point p);

   Vertex   nearest_neighbor(Vertex v) const;   
   
   Vertex   nearest_neighborA(Vertex v) const;
   
   Vertex   nearest_neighborD(Point p) const;

    template<class OutputIterator>
    OutputIterator   nearest_neighbors(Point p, int k,OutputIterator res)
    { leda_list<leda_node> result;
     int n = number_of_nodes();
     leda_node nd;
   
     if ( k <= 0 ) return res;
     if ( n <= k ) {
      result=all_nodes();
      forall(nd,result) { *res=nd; res++; }
      return res;
     }
   
     // insert p and search neighbors graph starting at p
   
     leda_node v = lookup(p);  
     bool old_node = true;
     if ( v == nil ) 
     { v = ((Point_set_2*)this)->insert(p);
       old_node = false;
       k++;
     } 
   
     result = nearest_neighbors(v,k);
      
     if ( !old_node ) 
     { result.pop();
       ((Point_set_2*)this)->del(v);
     }

     forall(nd,result) { *res=nd; res++; }
     return res;   
    }  

   leda_list<Vertex>   nearest_neighbors(Point p, int k);

   template<class OutputIterator>  
   OutputIterator  nearest_neighbors(Vertex v, int k,OutputIterator res) const
   { leda_list<leda_node> result;
     int n = number_of_nodes();
     leda_node nd;
   
     if ( k <= 0 ) return res;
     if ( n <= k ) {
       result=all_nodes();
       forall(nd,result) { *res=nd; res++; }
       return res;
     }
     
     Point p = pos(v);
   
     unmark_all_nodes();
     
     leda_p_queue<Numb_type,leda_node> PQ;
   
     PQ.insert(0,v); mark_node(v);
   
     while ( k > 0 )
     { pq_item it = PQ.find_min();
       leda_node w = PQ.inf(it); PQ.del_item(it);
   
       result.append(w); k--; 
   
       leda_node z;
       forall_adj_nodes(z,w)
       { if ( !is_marked(z) )
         { PQ.insert(tr_sqrdist(p,pos(z)),z);
           mark_node(z);
         }
       }
     }

     forall(nd,result) { *res=nd; res++; }   
     return res;
   } 
   
      
   leda_list<Vertex>   nearest_neighbors(Vertex v, int k) const;

   template<class OutputIterator>
   OutputIterator range_search(const Circle& C, OutputIterator res)
   { 
     leda_list<leda_node> L;
     leda_node nd;
   
     //int orient = C.orientation();  // Achtung !!!!!!!!
     //if (orient == 0)
     //    error_handler(1,"Point_set_2::range_search: circle must be proper");
   
     if (number_of_nodes() == 0) return res;
     if ( number_of_nodes() == 1 && !(tr_circleptori(C,pos(first_node()))==ON_UNBOUNDED_SIDE)) //changed
     { //L.append(first_node());
       //return L;
       *res=first_node(); res++;
       return res;
     }
   
     Point p = tr_circlecenter(C);
     leda_node v = lookup(p);  
     bool new_v = false;
   
     if ( v == nil )
     { new_v = true;
       v = insert(p); 
     }
   
     unmark_all_nodes(); 
   
     dfs(v,C,L);
   
     if (new_v)
     { L.pop();   
       del(v);
     }
     
     forall(nd,L){ *res=nd; res++; }
     return res;
   }
   
   leda_list<Vertex> range_search(const Circle& C);
   
   template<class OutputIterator> 
   OutputIterator range_search(Vertex v,const Point& p, OutputIterator res) const
   { 
    leda_list<leda_node> L;

    Point pv = pos(v);

    unmark_all_nodes(); 

    dfs(v,pv,p,L);

    leda_node nd;
    forall(nd,L){ *res=nd; res++; }
    return res;
  }


   leda_list<Vertex> range_search(Vertex v,const Point& p) const;


   template<class OutputIterator>
   OutputIterator range_search(const Point& a, const Point& b, const Point& c,OutputIterator res)
   { int orient = (int)(tr_orientation(a,b,c));
     Circle C = tr_createcircle_3p(a,b,c);
     leda_list<leda_node> L = range_search(C);
     list_item it = L.first_item();
     while (it != nil)
     { Point p = pos(L[it]);
       list_item next_it = L.succ(it);
       if ( ((int)(tr_orientation(a,b,p))) == - orient ||
            ((int)(tr_orientation(b,c,p))) == - orient ||
            ((int)(tr_orientation(c,a,p))) == - orient )      
          L.del_item(it);
       it = next_it;
     }
     leda_node nd;
     
     forall(nd,L){ *res=nd; res++; }
     return res;     
   }

   leda_list<Vertex> range_search(const Point& a, const Point& b, const Point& c);


   template<class OutputIterator>
   OutputIterator range_search(const Point& a1, const Point& b1, const Point& c1,const Point&
   d1,OutputIterator res)
   // a1 lower left, b1 lower right , c1 upper right
   {
     //Point b(c.xcoord(),a.ycoord());
     //Point d(a.xcoord(),c.ycoord());
     Point a=a1,b=b1,c=c1,d=d1;
   
     if (tr_orientation(a,b,c) == RIGHTTURN) 
     { Point tmp = b;
       b = d;
       d = tmp;
      }
   
     //W.set_color(leda_red);
     //W << a; W << b; W << c;
     Circle C = tr_createcircle_3p(a,b,c);
     //W << C;
   
     leda_list<leda_node> L = range_search(C);
     
     list_item it = L.first_item();
     while (it != nil)
     { Point p = pos(L[it]);
       list_item next_it = L.succ(it);
       if ( tr_orientation(a,b,p) == RIGHTTURN || tr_orientation(b,c,p) == RIGHTTURN ||
            tr_orientation(c,d,p) == RIGHTTURN || tr_orientation(d,a,p) == RIGHTTURN )
          L.del_item(it);
       it = next_it;
     }
     
     leda_node nd; 
     forall(nd,L){ *res=nd; res++; }
     return res;    
   }


   leda_list<Vertex> rectangular_range_search(const Point& a1, const Point& b1, const Point& c1,const Point& d1);

   static const Point_set_2<TR>* T_tmp; 

   static int cmp_edge_length(const leda_edge& e1, const leda_edge& e2)
   { Numb_type l1 = T_tmp->tr_sqrdist(T_tmp->pos_source(e1),T_tmp->pos_target(e1)); 
     Numb_type l2 = T_tmp->tr_sqrdist(T_tmp->pos_source(e2),T_tmp->pos_target(e2)); 
#if defined  POINT_SET_LEDA_ONLY
     return compare(l1,l2);
#else
     return CGAL::compare(l1,l2);
#endif
   }     
   
   leda_list<Edge> minimum_spanning_tree() const
   { T_tmp = this; 
     return MIN_SPANNING_TREE(*this, Point_set_2<TR>::cmp_edge_length);
   }   
   
   template<class OutputIterator>
   OutputIterator minimum_spanning_tree(OutputIterator result) const
   { 
     leda_list<Edge> EL;
     T_tmp = this; 
     EL=MIN_SPANNING_TREE(*this, Point_set_2<TR>::cmp_edge_length);
     Edge e;
     forall(e,EL) {  *result = e;  ++result; }
     return result;
   }      

   void   compute_voronoi(GRAPH<Circle,Point>& VD) const;
 
   void checking_on()
   { check = true; }

   void checking_off()
   { check = false; }

   bool check_state(const leda_string& location) const;
 
   bool     IS_NON_DIAGRAM_DART(Edge e) const;
  
  bool is_non_diagram_edge(Edge e) const
  {
    return IS_NON_DIAGRAM_DART(e); 
  }

   Edge get_cur_dart()  const { return cur_dart; }
   void set_cur_dart(Edge e)  { cur_dart = e;    }
   void set_hull_dart(Edge e) { hull_dart = e;   }

   Point pos(Vertex v) const
   { return inf(v); }
    
   Point pos_source(Edge e) const
   { return inf(source(e)); }
    
   Point pos_target(Edge e) const
   { return inf(target(e)); }
    
   Segment seg(Edge e) const           
   { return tr_createsegment_2p(inf(source(e)),inf(target(e))); }
    
   Line supporting_line(Edge e) const           
   { return tr_createline_2p(inf(source(e)),inf(target(e))); }
    
   Edge  get_hull_dart() const
   { return hull_dart; }
    
   Edge  get_hull_edge() const
   { return hull_dart; }
    
   bool is_diagram_dart(Edge e) const
   { return !(inf(e) & CGAL_NON_DIAGRAM_DART); }
    
   bool  is_diagram_edge(Edge e) const
   { return is_diagram_dart(e); }
    
   bool is_hull_dart(Edge e) const
   { return inf(e) & CGAL_HULL_DART; }
    
   bool is_hull_edge(Edge e) const
   { return is_hull_dart(e); }
    
   int orientation(Edge e, Point p) const
   { return ((int)tr_orientation(pos(source(e)),pos(target(e)),p)); }
    
   int dim() const
   { int n = number_of_nodes();
     if (n <= 1) 
        return n - 1;
     else
        return (is_hull_dart(reversal(hull_dart))) ? 1 : 2 ;
   }

};

template<class TR> const Point_set_2<TR>* Point_set_2<TR>::T_tmp;



// Impl. ...

   template<class TR>
   void Point_set_2<TR>::init_hull()
   { 
    hull_dart = nil;

    leda_edge e;
    forall_edges(e,*this)
    { if ( orientation(e,pos_target(face_cycle_succ(e))) <= 0 ) 
      { hull_dart = e;
        break;
      }
    } 

    if (hull_dart)
    { leda_edge e = hull_dart;
      do { mark_edge(e,CGAL_HULL_DART);
           e = face_cycle_succ(e);
         } 
      while (e != hull_dart);
    }
   }
   
   template<class TR>
   void Point_set_2<TR>::make_delaunay(leda_list<leda_edge>& S)
   {
    // Transforms graph into a Delaunay triangulation by flipping edges.
    // Diagonals of co-circular convex quadrilaterals are marked as 
    // CGAL_NON_DIAGRAM_DART
    // We maintain a stack $S$ of edges containing diagonals which might
    // have to be flipped. 

    if (number_of_nodes() <= 3) return;

    while ( !S.empty() )
    { leda_edge e = S.pop();
      leda_edge r = reversal(e);

    if (is_hull_dart(e) || is_hull_dart(r)) continue;

    mark_edge(e,CGAL_DIAGRAM_DART);
    mark_edge(r,CGAL_DIAGRAM_DART);

    // e1,e2,e3,e4: edges of quadrilateral with diagonal e
    leda_edge e1 = face_cycle_succ(r);
    leda_edge e3 = face_cycle_succ(e);

    // flip test
    Point a = pos_source(e1);
    Point b = pos_target(e1);
    Point c = pos_source(e3);
    Point d = pos_target(e3);

    if ((tr_orientation(b,d,a)==LEFTTURN) && (tr_orientation(b,d,c)==RIGHTTURN) )
    { // the quadrilateral is convex
      int soc = (int)(tr_so_circle(a,b,c,d));

      if (soc == 0) // co-circular quadrilateral(a,b,c,d) 
      { mark_edge(e,CGAL_NON_DIAGRAM_DART);
        mark_edge(r,CGAL_NON_DIAGRAM_DART);
       }
      if (soc > 0) // flip
      { leda_edge e2 = face_cycle_succ(e1);
        leda_edge e4 = face_cycle_succ(e3);
  
        S.push(e1); 
        S.push(e2); 
        S.push(e3); 
        S.push(e4); 
  
        // flip diagonal
        move_edge(e,e2,source(e4));
        move_edge(r,e4,source(e2));
      }
    }
  }

  } 
   
   template<class TR>
   void Point_set_2<TR>::make_delaunay()
   { leda_list<leda_edge> S = all_edges();
     make_delaunay(S);
   }

   template<class TR>
   void Point_set_2<TR>::check_locate(leda_edge answer,const Point& p) const
   { 
    if (answer == nil && dim() < 1) return;

    if ( tr_seg_contains(seg(answer),p) && ( is_hull_dart(answer)
         || (!is_hull_dart(answer) && !is_hull_dart(reversal(answer)) ))) 
    return;

    if (orientation(answer,p) < 0) 
    { std::cerr << "\norientation(" << seg(answer) << "," << p << ") < 0.";
      goto error_in_locate;
    }

    if ( is_hull_dart(answer) && orientation(answer,p) > 0 ) return;

    if (dim() == 1)
    { // orientation = 0 and answer does not contain p; so beyond extreme point
      leda_edge e = face_cycle_succ(answer);
      if ( e == reversal(answer) && !( ((int)(tr_so_hp(pos_source(e),pos_target(e),p))) >= 0 ) )
         return; 
      else 
      { std::cerr << "\n\ndim = 1 error.";
        goto error_in_locate;
      }
    }

    // dim == 2: answer must not be a hull edge and triangle must contain p 
    if ( orientation(answer,p) > 0 &&
         orientation(face_cycle_succ(answer), p) > 0 &&
         orientation(face_cycle_pred(answer), p) > 0 )
       return;
    else
    { std::cerr << "\n\ndim = 2 error";
      goto error_in_locate;
    }

   error_in_locate: 

   std::cerr << "\nAn error occurred in Point_set_2::locate(point).";
 
   exit(1);
  }
   
   template<class TR>
   void Point_set_2<TR>::dfs(leda_node s, const Circle& C, leda_list<leda_node>& L) const
   /* a procedure used in range_search(const Circle& C) */
  { L.append(s);
    mark_node(s);
    leda_node u;
    forall_adj_nodes(u,s)
      if (!is_marked(u) && ! (tr_circleptori(C,pos(u))==ON_UNBOUNDED_SIDE) ) dfs(u,C,L);  //was C.outside
  }
      
   template<class TR>
   void Point_set_2<TR>::dfs(leda_node s, const Point& pv, const Point& p, leda_list<leda_node>& L) const
   /* a procedure used in range_search(leda_node v, const Point& p) */
  { L.append(s);
    mark_node(s);
    leda_node u;
    forall_adj_nodes(u,s)
      if (!is_marked(u)){
       Comparison_result cr = tr_comparedist(pv,pos(u),p);
       if (cr==SMALLER || cr==EQUAL) dfs(u,pv,p,L);
      }
  }  
   
   template<class TR>
   void  Point_set_2<TR>::init(const leda_list<Point>& L0)
   { 
    // construct a triangulation for the points in L0

    clear();

    if (L0.empty()) return;
 
    leda_list<Point> L = L0;

    L.sort(tr_comparepoints);  

    // initialize graph with a single edge starting at the first point

    Point last_p = L.pop();          // last visited point
    leda_node  last_v = new_node(last_p); // last inserted node

    while (!L.empty() && last_p == L.head()) L.pop();

    if (L.empty()) return;

    last_p = L.pop();
    leda_node v = new_node(last_p);
    leda_edge x = new_edge(last_v,v);
    leda_edge y = new_edge(v,last_v);
    set_reversal(x,y);
    last_v = v;


    // scan remaining points

    Point p=last_p; // changed ....
    forall(p,L) 
    { if (p == last_p) continue; 

      leda_edge e =  last_adj_edge(last_v);

      last_v = new_node(p);
      last_p = p;

      // walk up to upper tangent
      do e = face_cycle_pred(e); while (orientation(e,p) > 0);

      // walk down to lower tangent and triangulate
      do { e = face_cycle_succ(e);
         leda_edge x = new_edge(e,last_v);
         leda_edge y = new_edge(last_v,source(e));
         set_reversal(x,y);
         } while (orientation(e,p) > 0);
     }

    // mark edges of convex hull

    hull_dart = last_edge();
  
    cur_dart = reversal(hull_dart);
 
    leda_edge e0 = hull_dart;
    leda_edge e  = e0;
    do { mark_edge(e,CGAL_HULL_DART);
       e = face_cycle_succ(e);
       } while (e != e0); 

    make_delaunay();

    init_node_marks();
  }    
  
   template<class TR>
   void  Point_set_2<TR>::init(const std::list<Point>& L0)
   { 
    // construct a triangulation for the points in L0

    clear();

    if (L0.empty()) return;
 
    leda_list<Point> L;

#if defined(__GNUC__)
    typename std::list<Point>::const_iterator it;
#else    
    std::list<Point>::const_iterator it;
#endif

    it = L0.begin();
    
    for(;it != L0.end(); ++it) L.append(*it);

    L.sort(tr_comparepoints);  

    // initialize graph with a single edge starting at the first point

    Point last_p = L.pop();          // last visited point
    leda_node  last_v = new_node(last_p); // last inserted node

    while (!L.empty() && last_p == L.head()) L.pop();

    if (L.empty()) return;

    last_p = L.pop();
    leda_node v = new_node(last_p);
    leda_edge x = new_edge(last_v,v);
    leda_edge y = new_edge(v,last_v);
    set_reversal(x,y);
    last_v = v;


    // scan remaining points

    Point p=last_p; // changed ....
    forall(p,L) 
    { if (p == last_p) continue; 

      leda_edge e =  last_adj_edge(last_v);

      last_v = new_node(p);
      last_p = p;

      // walk up to upper tangent
      do e = face_cycle_pred(e); while (orientation(e,p) > 0);

      // walk down to lower tangent and triangulate
      do { e = face_cycle_succ(e);
         leda_edge x = new_edge(e,last_v);
         leda_edge y = new_edge(last_v,source(e));
         set_reversal(x,y);
         } while (orientation(e,p) > 0);
     }

    // mark edges of convex hull

    hull_dart = last_edge();
  
    cur_dart = reversal(hull_dart);
 
    leda_edge e0 = hull_dart;
    leda_edge e  = e0;
    do { mark_edge(e,CGAL_HULL_DART);
       e = face_cycle_succ(e);
       } while (e != e0); 

    make_delaunay();

    init_node_marks();
  }    

   template<class TR>
   leda_list<CGAL_TYPENAME_MSVC_NULL Point_set_2<TR>::Point>  Point_set_2<TR>::points() const
   { leda_list <Point> L;
     leda_node v;
     forall_nodes(v,*this) L.append(pos(v));
     return L;
   }      
  
   template<class TR>
   Point_set_2<TR>::Edge   
   Point_set_2<TR>::locate(Point p) const
   { 
    if (number_of_edges() == 0) return nil;

    if (dim() == 1)
    { 
      leda_edge e = hull_dart;
      int orient = orientation(e,p);
      if (orient != 0)
      { if (orient < 0) e = reversal(e);
        if (check) check_locate(e ,p);
        return e;
      }

      // p is collinear with the points in S. We walk 
           
      if ( !(((int) (tr_so_hp(pos_source(e),pos_target(e),p))) >= 0 ) ) e = reversal(e);

      // in the direction of e. We know IN_HALFSPACE(e,p) 
         
      leda_edge e1 = face_cycle_succ(e);
      while ( e1 != reversal(e) && ( ((int)(tr_so_hp(pos_source(e1),pos_target(e1),p))) >= 0 ) ) 
      { e = e1;  
        e1 = face_cycle_succ(e); 
      }
         
      if (check) check_locate(e ,p);

      return e;
    } 

   
    leda_edge e = is_hull_dart(cur_dart) ? reversal(cur_dart) : cur_dart;

    if (p == pos_source(e) ) return reversal(e);
    
    int orient = orientation(e,p);

    if (orient == 0) 
    { e = face_cycle_pred(e);
      orient = orientation(e,p);
    }

    if (orient < 0) e = reversal(e);
   
    //Segment s(pos_source(e),p);
    Point ps1=pos_source(e);
    Point ps2=p;

    while ( true )     
    { 
      if (is_hull_dart(e)) break;
    
      leda_edge e1 = face_cycle_succ(e);
      leda_edge e2 = face_cycle_pred(e);
      Orientation ori = tr_orientation(ps1,ps2,pos_target(e1));
      int d= (int)ori;
      leda_edge e_next = reversal( (d < 0) ? e2 : e1 );
      int orient = orientation(e_next,p);
      if ( orient > 0 )  { e = e_next; continue; }
      if ( orient == 0 ) { e = e_next; break; }
      if ( d == 0 && orient < 0 && orientation(e2,p) == 0 ) 
       e = reversal(e2);
      break;

    }

    if (check) check_locate(e ,p);

    ((leda_edge&)cur_dart) = e;
    return e;
  }

   template<class TR>
   Point_set_2<TR>::Vertex   
   Point_set_2<TR>::lookup(Point p) const
   { 
     if (number_of_nodes() == 0) return nil;   
     
     if (number_of_nodes() == 1)
     { leda_node v = first_node();
       return (pos(v) == p) ? v : nil;
     }
     leda_edge e = locate(p);
     if (pos(source(e)) == p) return source(e);
     if (pos(target(e)) == p) return target(e);
     return nil;
   }

   template<class TR>
   Point_set_2<TR>::Vertex   
   Point_set_2<TR>::insert(Point p)
   { 
     leda_node v; 
  
     if (number_of_nodes() == 0)  
     { v = new_node(p); 
                     if ( check && !check_state("Point_set_2::insert") )   
                     { std::cerr << "The point inserted was " << p; 
                       exit(1);
                     }
                     return v;
     }

     if (number_of_nodes() == 1)
     { leda_node w = first_node();
       if (p == pos(w)) 
       { assign(w,p); 
         v = w; 
      
         if ( check && !check_state("Point_set_2::insert") )   
         { std::cerr << "The point inserted was " << p; 
           exit(1);
         }
         return v;

       }
       else
       { v = new_node(p);
         leda_edge x = new_edge(v,w);
         leda_edge y = new_edge(w,v);
         mark_edge(x,CGAL_HULL_DART);
         mark_edge(y,CGAL_HULL_DART);
         set_reversal(x,y);
         hull_dart = x;
         cur_dart = x;
      
         if ( check && !check_state("Point_set_2::insert") )   
         { std::cerr << "The point inserted was " << p; 
           exit(1);
         }
         return v;

       }
     }


     leda_edge e = locate(p);
     if (p == pos_source(e)) 
       { assign(source(e),p); return source(e); }
     if (p == pos_target(e)) 
       { assign(target(e),p); return target(e); } 

     bool p_on_e = tr_seg_contains(seg(e),p);


     if ( dim() == 1 && orientation(e,p) == 0 )
     { 
       v = new_node(p);
       leda_edge x = new_edge(v,target(e));
       leda_edge y = new_edge(target(e),v);
       mark_edge(x,CGAL_HULL_DART);
       mark_edge(y,CGAL_HULL_DART);
       set_reversal(x,y);
   
       if (p_on_e)
       { x = new_edge(v,source(e));  
         y = new_edge(source(e),v);
         mark_edge(x,CGAL_HULL_DART);
         mark_edge(y,CGAL_HULL_DART);
         set_reversal(x,y);
         hull_dart = cur_dart = x;

         del_edge(reversal(e));
         del_edge(e);
       }

    
       if ( check && !check_state("Point_set_2::insert") )   
       { std::cerr << "The point inserted was " << p; 
         exit(1);
       }
       return v;

    }

    v  = new_node(p);
    leda_edge e1 = e;
    leda_edge e2 = e;
    leda_list<leda_edge> E;
    bool outer_face = is_hull_dart(e);

    if (outer_face) //  move e1/e2 to compute upper/lower tangents
    { do e1 = face_cycle_pred(e1); while (orientation(e1,p) > 0);
      do e2 = face_cycle_succ(e2); while (orientation(e2,p) > 0);
    }

    // insert edges between v and target(e1) ... source(e2)
    e = e1;
    do { e = face_cycle_succ(e);
       leda_edge x = new_edge(e,v);
       leda_edge y = new_edge(v,source(e));
       set_reversal(x,y);
       mark_edge(e,CGAL_DIAGRAM_DART);
       E.append(e);
       E.append(x);
     } while (e != e2);

    if (outer_face) // mark last visited and new edges as hull edges
    { mark_edge(face_cycle_succ(e1),CGAL_HULL_DART);
      mark_edge(face_cycle_pred(e2),CGAL_HULL_DART);
      mark_edge(e2,CGAL_HULL_DART);
      hull_dart = e2;
    }

    make_delaunay(E); // restores Delaunay property
  
    if ( check && !check_state("Point_set_2::insert") )   
    { std::cerr << "The point inserted was " << p; 
      exit(1);
    }
    return v;

   }

   
   template<class TR>
   void Point_set_2<TR>::del(Vertex v)
   {  
     if (v == nil) 
        error_handler(1,"Point_set_2::del: nil argument.");
   
     if (number_of_nodes() == 0) 
        error_handler(1,"Point_set_2::del: graph is empty.");
   
     if ( dim() < 2 )
     { 
       if ( outdeg(v) == 2)
       { leda_node s = target(first_adj_edge(v));
         leda_node t = target(last_adj_edge(v));
         leda_edge x = new_edge(s,t);  
         leda_edge y = new_edge(t,s);
         mark_edge(x,CGAL_HULL_DART);
         mark_edge(y,CGAL_HULL_DART);
         set_reversal(x,y);
       }
   
       del_node(v);
       cur_dart = hull_dart = first_edge();    
     
       
       if ( check && !check_state("Point_set_2::del(node v)") )
       { std::cerr << "deleted the node with position " << pos(v);
         exit(1);
       }
       return;
   
    }
    
    leda_list<leda_edge> E;
   
    int min_deg = 3;
   
    leda_edge e;
    forall_adj_edges(e,v) 
    { E.append(face_cycle_succ(e));
      if (is_hull_dart(e)) min_deg = 2;
    }
   
    int count = 0;
    e = first_adj_edge(v);
   
    while ( outdeg(v) > min_deg && count < outdeg(v) )
    { leda_edge e_pred = cyclic_adj_pred(e);
      leda_edge e_succ = cyclic_adj_succ(e);
      Point a = pos_target(e_pred); Point c = pos_target(e_succ);
   
      if ( ! (tr_orientation(a,c,pos(v))==RIGHTTURN) && (tr_orientation(a,c,pos_target(e))==RIGHTTURN) )
      { // e is flipable
        leda_edge r = reversal(e);
   
        move_edge(e,reversal(e_succ),target(e_pred));

#if (__LEDA__ > 371)
        move_edge(r,reversal(e_pred),target(e_succ),LEDA::before);
#else
        move_edge(r,reversal(e_pred),target(e_succ),leda_before);
#endif
   
        mark_edge(e,CGAL_DIAGRAM_DART);
        mark_edge(r,CGAL_DIAGRAM_DART);
        E.append(e);
   
        e = e_pred;
        count = count - 2;    
        if ( count < 0 ) count = 0;
      }
      else
      { e = e_succ;
        count++;
      }
    }
   
    if ( min_deg == 2 )
    {   
      leda_edge e,x=NULL; // = NULL to supress warning
      forall_adj_edges(e,v)
      { x = face_cycle_succ(e);
        mark_edge(x,CGAL_HULL_DART);
        if ( !is_hull_dart(reversal(x)) ) 
          mark_edge(reversal(x),CGAL_DIAGRAM_DART);
      }
      hull_dart = x;
    }
   
    cur_dart = E.head();
       
    del_node(v);
    make_delaunay(E);
   
   
    
    if ( check && !check_state("Point_set_2::del(node v)") )
    { std::cerr << "deleted the node with position " << pos(v);
      exit(1);
    }
    return;
   
   }

   template<class TR>
   void Point_set_2<TR>::del(Point p)
   { leda_node v = lookup(p);
     if ( v != nil ) del(v);
   }


   template<class TR>
   Point_set_2<TR>::Vertex   
   Point_set_2<TR>::nearest_neighbor(Point p)
    {
     if (number_of_nodes() == 0) return nil;
   
     leda_node v = lookup(p);
     leda_node min_v;
   
     if ( v != nil ) return v;
   
     // insert p and search neighbors of v
   
     leda_node w = insert(p);  
   
     min_v = nearest_neighbor(w);
   
     del(w);
   
     return min_v;
   }
     

   template<class TR>
   Point_set_2<TR>::Vertex   
   Point_set_2<TR>::nearest_neighbor(Vertex v) const
   {
     if (number_of_nodes() <= 1) return nil;
   
     Point p = pos(v);
     leda_edge e = first_adj_edge(v);
   
     leda_node min_v = target(e);
   
     while ((e = adj_succ(e)) != nil)
     { leda_node u = target(e);  
       if ( tr_comparedist(p,pos(u),pos(min_v)) == SMALLER ) min_v = u;
     }
   
     return min_v;
   }
   
   template<class TR>
   Point_set_2<TR>::Vertex   
   Point_set_2<TR>::nearest_neighborA(Vertex v) const
   {
     if (number_of_nodes() <= 1) return nil;
   
     Point p = pos(v);
     leda_edge e = first_adj_edge(v);
   
     
     leda_node min_v = target(e);
     Numb_type min_d = tr_sqrdist(p,pos(min_v));
   
     while ((e = adj_succ(e)) != nil)
     { leda_node u = target(e); 
       Numb_type d_u = tr_sqrdist(p,pos(u)); 
       if ( d_u < min_d ) 
       { min_v = u;
         min_d = d_u;
       }
     }  
   
     return min_v;
   }  
   
   template<class TR>
   Point_set_2<TR>::Vertex   
   Point_set_2<TR>::nearest_neighborD(Point p) const
   { 
     if (number_of_nodes() == 0) return nil;
     if (number_of_nodes() == 1) return first_node();
   
     leda_edge e = locate(p);
   
     if ( is_hull_dart(e) ) 
     { while ( ! ( ((int)( tr_so_hp(pos_source(e),pos_target(e),p))) >= 0 ) )          e = face_cycle_pred(e);
       while ( ! ( ((int)( tr_so_hp(pos_source(reverse(e)),pos_target(reverse(e)),p)) >= 0 ) )) e = face_cycle_succ(e);
     }
   
     unmark_all_nodes(); 
   
     leda_node min_v = source(e);
     Numb_type min_d = tr_sqrdist(p,pos(min_v));
   
     leda_list<leda_node> L;
     L.append(source(e)); 
     L.append(target(e)); 
     mark_node(source(e));
     mark_node(target(e));
   
     while ( !L.empty() )
     { leda_node v = L.pop();
   
       if ( tr_sqrdist(p,pos(v)) < min_d )
       { min_v = v;
         min_d = tr_sqrdist(p,pos(v));
       }
   
       forall_adj_edges(e,v)
       { leda_node w = target(e);
         if ( !is_marked(target(e)) && 
              ( ((int)(tr_so_hp(pos_source(e),pos_target(e),p))) >= 0 ) && \
	      (((int)(tr_so_hp(pos_source(reversal(e)),pos_target(reversal(e)),p)) >= 0 )) && \
              tr_linesqrdist(supporting_line(e),p) < min_d ) 
         { L.append(w); 
           mark_node(w);
         }
       }
     }
   
     return min_v;
   
   }
 
  
   template<class TR>
   leda_list<CGAL_TYPENAME_MSVC_NULL Point_set_2<TR>::Vertex>   
   Point_set_2<TR>::nearest_neighbors(Point p, int k)
    { leda_list<leda_node> result;
     int n = number_of_nodes();
   
     if ( k <= 0 ) return result;
     if ( n <= k ) return all_nodes();
   
     // insert p and search neighbors graph starting at p
   
     leda_node v = lookup(p);  
     bool old_node = true;
     if ( v == nil ) 
     { v = ((Point_set_2*)this)->insert(p);
       old_node = false;
       k++;
     } 
   
     result = nearest_neighbors(v,k);
      
     if ( !old_node ) 
     { result.pop();
       ((Point_set_2*)this)->del(v);
     }
   
     return result;
   } 


   template<class TR>
   leda_list<CGAL_TYPENAME_MSVC_NULL Point_set_2<TR>::Vertex>   
   Point_set_2<TR>::nearest_neighbors(Vertex v, int k) const
   { leda_list<leda_node> result;
     int n = number_of_nodes();
   
     if ( k <= 0 ) return result;
     if ( n <= k ) return all_nodes();
   
     Point p = pos(v);
   
     unmark_all_nodes();
     
     leda_p_queue<Numb_type,leda_node> PQ;
   
     PQ.insert(0,v); mark_node(v);
   
     while ( k > 0 )
     { pq_item it = PQ.find_min();
       leda_node w = PQ.inf(it); PQ.del_item(it);
   
       result.append(w); k--; 
   
       leda_node z;
       forall_adj_nodes(z,w)
       { if ( !is_marked(z) )
         { PQ.insert(tr_sqrdist(p,pos(z)),z);
           mark_node(z);
         }
       }
     }
   
     return result;
   } 


   template<class TR>
   leda_list<CGAL_TYPENAME_MSVC_NULL Point_set_2<TR>::Vertex> 
   Point_set_2<TR>::range_search(const Circle& C)
   { 
     leda_list<leda_node> L;
   
     //int orient = C.orientation();  // Achtung !!!!!!!!
     //if (orient == 0)
     //    error_handler(1,"Point_set_2::range_search: circle must be proper");
   
     if (number_of_nodes() == 0) return L;
     if ( number_of_nodes() == 1 && !(tr_circleptori(C,pos(first_node()))==ON_UNBOUNDED_SIDE )) //changed
     { L.append(first_node());
       return L;
     }
   
     Point p = tr_circlecenter(C);
     leda_node v = lookup(p);  
     bool new_v = false;
   
     if ( v == nil )
     { new_v = true;
       v = insert(p); 
     }
   
     unmark_all_nodes(); 
   
     dfs(v,C,L);
   
     if (new_v)
     { L.pop();   
       del(v);
     }
     return L;
   }

   template<class TR>
   leda_list<CGAL_TYPENAME_MSVC_NULL Point_set_2<TR>::Vertex> 
   Point_set_2<TR>::range_search(Vertex v,const Point& p) const
   { 
    leda_list<leda_node> L;

    Point pv = pos(v);

    unmark_all_nodes(); 

    dfs(v,pv,p,L);

    return L;
  }

   template<class TR>
   leda_list<CGAL_TYPENAME_MSVC_NULL Point_set_2<TR>::Vertex> 
   Point_set_2<TR>::range_search(const Point& a, const Point& b, const Point& c)
   { int orient = (int)(tr_orientation(a,b,c));
     Circle C = tr_createcircle_3p(a,b,c);
     leda_list<leda_node> L = range_search(C);
     list_item it = L.first_item();
     while (it != nil)
     { Point p = pos(L[it]);
       list_item next_it = L.succ(it);
       if ( ((int)(tr_orientation(a,b,p))) == - orient ||
            ((int)(tr_orientation(b,c,p))) == - orient ||
            ((int)(tr_orientation(c,a,p))) == - orient )      
          L.del_item(it);
       it = next_it;
     }
     return L;
   }

   template<class TR>
   leda_list<CGAL_TYPENAME_MSVC_NULL Point_set_2<TR>::Vertex> 
   Point_set_2<TR>::rectangular_range_search(const Point& a1, const Point& b1, const Point& c1,const Point& d1)
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
   
     leda_list<leda_node> L = range_search(C);
     list_item it = L.first_item();
     while (it != nil)
     { Point p = pos(L[it]);
       list_item next_it = L.succ(it);
       if ( tr_orientation(a,b,p) == RIGHTTURN || tr_orientation(b,c,p) == RIGHTTURN ||
            tr_orientation(c,d,p) == RIGHTTURN || tr_orientation(d,a,p) == RIGHTTURN )
          L.del_item(it);
       it = next_it;
     }
     return L;
   }
           
   template<class TR>
   void  Point_set_2<TR>::compute_voronoi(GRAPH<Circle,Point>& VD) const
   { 
     VD.clear();
   
     if (number_of_nodes() < 2) return;
   
     // create Voronoi nodes
   
     leda_edge_array<leda_node> vnode(*this,nil);
   
     // for outer face
   
     leda_edge  e = hull_dart;
     Point a = pos_source(e);
     leda_edge  x = e;
     do { Point b = pos_target(x);
          vnode[x] =  VD.new_node(tr_createcircle_3p(a,a,b)); // midpoint has to be changed (was midpoint(a,b) instead of a
          a = b;
          x = face_cycle_succ(x);
        } while ( x != e);
   
     // for all other faces
   
     forall_edges(e,*this)
     { 
       if (vnode[e] || !is_diagram_dart(e)) continue;
   
       leda_edge  x = face_cycle_succ(e);
       Point a = pos_source(e);
       Point b = pos_target(e);
       Point c = pos_target(x);
       leda_node  v = VD.new_node(tr_createcircle_3p(a,b,c));
   
       x = e;
       do { vnode[x] = v;
            x = d_face_cycle_succ(x);
           } while( x != e);
      }
   
     // construct Voronoi edges
   
     // for outer face
   
     e = hull_dart;
     x = e;
     do { leda_edge r = reversal(x); 
          Point p = pos_target(x);
          VD.new_edge(vnode[x],vnode[r],p);
          x = cyclic_adj_pred(r);
        } while ( x != e);
   
     // for all other faces
   
     forall_edges(e,*this)
     { leda_node v = vnode[e]; 
   
       if (!is_diagram_dart(e) || VD.outdeg(v) > 0) continue;
   
       Point p = pos_target(e);
       leda_edge x = e;
       do { leda_edge r = reversal(x); 
            VD.new_edge(v,vnode[r],p);
            x = d_face_cycle_succ(x);
          } while (x != e);
      }
      // add a check for correctness
   }

   template<class TR>
   bool Point_set_2<TR>::check_state(const leda_string& location) const
   { 
    if ( !check ) return true;

    write("Point_set_2_error.graph_after_op");

    // check hull_dart and cur_dar
    if ( (hull_dart == nil || cur_dart == nil)
       && (number_of_nodes() >= 2) )
    { std::cerr << "\nhull_dart or cur_dart contradicts number of nodes\n";
      
      std::cerr << "\n\nCheck_state was called at " << location;
      std::cerr << "\n\nThe situation before the call of " << location;
      std::cerr << "\nwas saved into files Point_set_2_error.graph";
      std::cerr << "\nand Point_set_2_error.aux.";
      std::cerr << "\n\nPlease send these files to ledares@mpi-sb.mpg.de.\n\n";

      return false;

    }     
    if ( hull_dart == nil ) return true;
  
    // check edge labels
    leda_edge e;
    forall_edges(e,*this)
    { switch ( inf(e)) {

      case CGAL_HULL_DART:
      { leda_edge next = face_cycle_succ(e);
        int orient = orientation(e,pos_target(next));
        if ( orient > 0 )
        { std::cerr << "\n\nwrongly labeled hull_dart.";
        
          std::cerr << "\n\nCheck_state was called at " << location;
          std::cerr << "\n\nThe situation before the call of " << location;
          std::cerr << "\nwas saved into files Point_set_2_error.graph";
          std::cerr << "\nand Point_set_2_error.aux.";
          std::cerr << "\n\nPlease send these files to ledares@mpi-sb.mpg.de.\n\n";

          return false;

        }
        break;
      }
          
      case CGAL_DIAGRAM_DART:
      { if ( inf(reversal(e)) != CGAL_HULL_DART &&
                          IS_NON_DIAGRAM_DART(e) )
        { std::cerr << "\n\nwrongly labeled diagram_dart.";
        
          std::cerr << "\n\nCheck_state was called at " << location;
          std::cerr << "\n\nThe situation before the call of " << location;
          std::cerr << "\nwas saved into files Point_set_2_error.graph";
          std::cerr << "\nand Point_set_2_error.aux.";
          std::cerr << "\n\nPlease send these files to ledares@mpi-sb.mpg.de.\n\n";

          return false;

        }
        break;
      }
          
      case CGAL_NON_DIAGRAM_DART:
      { if ( inf(reversal(e)) != CGAL_NON_DIAGRAM_DART ||
                          !IS_NON_DIAGRAM_DART(e) )
        { std::cerr << "\n\nwrongly labeled non_diagram_dart.";
        
          std::cerr << "\n\nCheck_state was called at " << location;
          std::cerr << "\n\nThe situation before the call of " << location;
          std::cerr << "\nwas saved into files Point_set_2_error.graph";
          std::cerr << "\nand Point_set_2_error.aux.";
          std::cerr << "\n\nPlease send these files to ledares@mpi-sb.mpg.de.\n\n";

          return false;

        }
        break;
      }


      default:  std::cerr << "\n\nillegal edge label.";
            
            std::cerr << "\n\nCheck_state was called at " << location;
            std::cerr << "\n\nThe situation before the call of " << location;
            std::cerr << "\nwas saved into files Point_set_2_error.graph";
            std::cerr << "\nand Point_set_2_error.aux.";
            std::cerr << "\n\nPlease send these files to ledares@mpi-sb.mpg.de.\n\n";

            return false;

      }
   }
                     
   if ( inf(hull_dart) != CGAL_HULL_DART )
   { std::cerr << "\n\nis_hull_dart gives wrong information.";
   
     std::cerr << "\n\nCheck_state was called at " << location;
     std::cerr << "\n\nThe situation before the call of " << location;
     std::cerr << "\nwas saved into files Point_set_2_error.graph";
     std::cerr << "\nand Point_set_2_error.aux.";
     std::cerr << "\n\nPlease send these files to ledares@mpi-sb.mpg.de.\n\n";

     return false;

   }
   return true;
  } 
   
   template<class TR>
   bool     Point_set_2<TR>::IS_NON_DIAGRAM_DART(Edge e) const
   { leda_edge r = reversal(e);
  
    // e1,e2,e3,e4: edges of quadrilateral with diagonal e
    leda_edge e1 = face_cycle_succ(r);
    leda_edge e3 = face_cycle_succ(e);

    // flip test
    Point a = pos_source(e1);
    Point b = pos_target(e1);
    Point c = pos_source(e3);
    Point d = pos_target(e3);

    if ((tr_orientation(a,b,c)==LEFTTURN) && (tr_orientation(b,c,d)==LEFTTURN) && 
        (tr_orientation(c,d,a)==LEFTTURN) && (tr_orientation(d,a,c)==LEFTTURN) && 
        (tr_so_circle(a,b,c,d)==ON_ORIENTED_BOUNDARY) )
        return true;
    return false;
  }


CGAL_END_NAMESPACE



#if LEDA_ROOT_INCL_ID == 400900
#undef LEDA_ROOT_INCL_ID
#include <LEDA/UNDEFINE_NAMES.h>
#endif


#endif

