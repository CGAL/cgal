// Copyright (c) 1997-2000  Max-Planck-Institute Saarbruecken (Germany).
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
// Author(s)     : Michael Seel <seel@mpi-sb.mpg.de>

#ifndef CGAL_PM_PERSISTENT_PL_H
#define CGAL_PM_PERSISTENT_PL_H

#include <CGAL/Nef_2/gen_point_location.h>

template <typename PMPL>
struct PM_persistent_PL_traits 
{
  typedef PMPL  Graph;
  typedef typename PMPL::Vertex_const_handle   Node;
  typedef typename PMPL::Halfedge_const_handle Edge;
  typedef typename PMPL::Face_const_handle     Face;
  typedef typename PMPL::Object_handle         Object_handle;

  typedef typename PMPL::Geometry  Geometry;
  typedef typename PMPL::Point     Point;
  typedef typename PMPL::Segment   Segment;
  const Geometry* pK;

  typedef typename PMPL::Vertex_const_iterator NodeIterator;
  NodeIterator Nodes_begin(const Graph& G) const { return G.vertices_begin(); }
  NodeIterator Nodes_end(const Graph& G) const { return G.vertices_end(); }
  Node toNode(const NodeIterator& nit) const { return nit; }

  typedef typename PMPL::Halfedge_around_vertex_const_circulator HAVC;
  struct IncEdgeIterator {
    HAVC _start, _curr;
    bool met;
    IncEdgeIterator() {}
    IncEdgeIterator(HAVC c) : 
      _start(c), _curr(c), met(false) {}
    IncEdgeIterator& operator++()
    { if (_curr==_start)
        if (!met)  { met=true; ++_curr; }
        else       { _curr=HAVC(); }
      else ++_curr;
      return *this; 
    }
    bool operator==(const IncEdgeIterator& it2) const
    { return _curr==it2._curr; }
    bool operator!=(const IncEdgeIterator& it2) const
    { return !(*this==it2); }
  };
  Edge toEdge(const IncEdgeIterator& eit) const { return eit._curr; }

  IncEdgeIterator IncEdges_begin(const Graph& G, const Node& n) 
  { return IncEdgeIterator(HAVC(G.first_out_edge(n))); }
  IncEdgeIterator IncEdges_end(const Graph& G, const Node& n)   
  { return IncEdgeIterator(); }

  enum EdgeCategory 
  { StartingNonVertical, StartingVertical, EndingNonVertical, EndingVertical };

  Node opposite(const Graph& G, const Edge& e, const Node& u)
  { if ( G.source(e) == u ) return G.target(e);
    else                    return G.source(e); }

  EdgeCategory ClassifyEdge(const Graph& G, const Edge& e, const Node& u)
  {
    Point p_u = G.point(u);
    Point p_v = G.point(opposite(G,e,u));

    int cmpX = pK->compare_x(p_u, p_v);
    if ( cmpX < 0 ) return StartingNonVertical;
    if ( cmpX > 0 ) return EndingNonVertical;

    int cmpY = pK->compare_y(p_u, p_v); 
    CGAL_assertion(cmpY != 0);
    if ( cmpY < 0 ) return StartingVertical;
    return EndingVertical;
  }    

  typedef Point XCoord;
  const XCoord getXCoord(const Point& p) const 
  { return p; }
  const XCoord getXCoord(const Graph& G, const Node& n) const 
  { return G.point(n); }

  class PredLessThanX {
    const Geometry* pK;
  public:
    PredLessThanX() : pK(0) {}
    PredLessThanX(const Geometry* pKi) : pK(pKi) {}
    PredLessThanX(const PredLessThanX& P) : pK(P.pK) 
    { CGAL_NEF_TRACEN("copy PredLessThanX"); }
    int operator() (const XCoord& x1, const XCoord& x2) const
    { return pK->compare_x(x1,x2) < 0; }
  };
  PredLessThanX getLessThanX() const { return PredLessThanX(pK); }

  // Curve connected functionality:
  typedef Segment  Curve;

  Curve makeCurve(const Point& p) const 
  { return pK->construct_segment(p,p); }
  Curve makeCurve(const Graph& G, const Node& n) const
  { return makeCurve(G.point(n)); }
  Curve makeCurve(const Graph& G, const Edge& e) const
  { Point ps = G.point(G.source(e)), pt = G.point(G.target(e));
    Curve res(G.point(G.source(e)),G.point(G.target(e)));
    if ( pK->compare_xy(ps,pt) < 0 ) res = pK->construct_segment(ps,pt);
    else                             res = pK->construct_segment(pt,ps);
    return res; 
  }

  struct PredCompareCurves {
   const Geometry* pK;
   PredCompareCurves() : pK(0) {}
   PredCompareCurves(const Geometry* pKi) : pK(pKi) {}
   PredCompareCurves(const PredCompareCurves& P) : pK(P.pK) {}

   int cmppntseg(const Point& p, const Curve& s) const
   { 
     if ( pK->compare_x(pK->source(s),pK->target(s)) != 0 ) // !vertical
       return pK->orientation(pK->source(s),pK->target(s), p); 
     if ( pK->compare_y(p,pK->source(s)) <= 0 ) return -1;
     if ( pK->compare_y(p,pK->target(s)) >= 0 ) return +1;
     return 0;
   }

   int operator()(const Curve& s1, const Curve& s2) const
   { 
     Point a = pK->source(s1); 
     Point b = pK->target(s1);
     Point c = pK->source(s2);
     Point d = pK->target(s2);
     if ( a==b ) 
       if ( c==d ) return pK->compare_y(a,c);
       else        return  cmppntseg(a, s2);
     if ( c==d )   return -cmppntseg(c, s1);
     // now both are non-trivial:
     int cmpX = pK->compare_x(a, c);
     if ( cmpX < 0 ) 
       return - pK->orientation(a,b,c);
     if ( cmpX > 0 ) 
       return   pK->orientation(c,d,a);

     int cmpY = pK->compare_y(a, c);
     if ( cmpY < 0 ) return -1;
     if ( cmpY > 0 )  return +1;
      
     // cmpX == cmpY == 0 => a == c
     return pK->orientation(c,d,b);
   }
  };

  PredCompareCurves getCompareCurves() const
  { return PredCompareCurves(pK); }

  typedef GenericLocation<Node, Edge> Location;
  typedef Object_handle QueryResult;

  virtual Object_handle 
  PostProcess(const Location& L, const Location& L_plus, 
    const Point& p) const 
  { /* we only get an L_plus (non-nil) if L is ABOVE a vertex
       in which case we want to extract the face from the edge
       below (p+epsilon) available via L_plus. */
    if (!L_plus.is_nil()) { CGAL_assertion(L_plus.is_edge());
      return Object_handle(Edge(L_plus));
    } else { 
      if ( L.is_edge() ) {
        return Object_handle(Edge(L));
      }
      if ( L.is_node() ) {
        Node v(L); CGAL_assertion( v != Node() );
        return Object_handle(v);
      }
      return Object_handle();
    }
  }

  PM_persistent_PL_traits() : pK(0) {}
  PM_persistent_PL_traits(const Geometry& k) : pK(&k) {}
  virtual ~PM_persistent_PL_traits() {}
  virtual void sweep_begin(const Graph&) {}
  virtual void sweep_moveto(const XCoord&) {}
  virtual void sweep_end() {}
  virtual void clear() {}

};


#endif // CGAL_PM_PM_PERSISTENT_PL_H
