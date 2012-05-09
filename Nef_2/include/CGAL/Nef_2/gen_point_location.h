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


#ifndef CGAL_NEF_2_GEN_POINT_LOCATION_H
#define CGAL_NEF_2_GEN_POINT_LOCATION_H

#include <CGAL/LEDA_basic.h>
#if CGAL_LEDA_VERSION < 500
#include <LEDA/pp_dictionary.h>
#else
#include <LEDA/core/pp_dictionary.h>
#endif
#include <sstream>
#include <string>
#include <list>
#include <vector>
#include <map>
#ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
#include <CGAL/Nef_2/geninfo.h>
#else
#include <boost/any.hpp>
#endif


#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 17
#include <CGAL/Nef_2/debug.h>

// #define CHECKING_OFF

// for dictionary
template <class Node>
inline std::ostream& operator<<(std::ostream& o, const std::list<Node>& n) 
{ return o; }

/*{\Manpage {GenericLocation}{Node, Edge}
{Return Type for Planar Point Location}{L}}*/

template <class Node, class Edge>
class GenericLocation {
/*{\Mdefinition
  An instance of the data type |\Mtype| is used as return value for planar 
  point location. It can store a node or an edge of a graph or the special
  value |nil| which is used to signal that no node or edge could be found.
}*/
  #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
  typedef void*  GenPtr;
  #else
  typedef boost::any GenPtr;
  #endif
public:
/*{\Mtypes}*/ 
  enum Type { NIL, NODE, EDGE };
  /*{\Menum This enumeration allows to specify the 3 basic types of the 
     values that a |\Mtype| can represent.}*/

  /*{\Mcreation}*/ 
    GenericLocation()       { init(); }
    /*{\Mcreate creates a |\Mtype| and initializes with the value |nil|.}*/
    GenericLocation(Node n) { init(n); }
    /*{\Mcreate creates a |\Mtype| and initializes with the node |n|.}*/
    GenericLocation(Edge e) { init(e); }
    /*{\Mcreate creates a |\Mtype| and initializes with the edge |e|.}*/

    ~GenericLocation() { clear(); }

    GenericLocation(const GenericLocation<Node, Edge>& L) { assign(L); }
    GenericLocation<Node, Edge>& operator=(
      const GenericLocation<Node, Edge>& L)
    { clear(); assign(L); return *this; }

  /*{\Moperations}*/
    operator const Node&() const
    {
  #if !defined(CHECKING_OFF)
      if (type != NODE) 
        CGAL_LEDA_SCOPE::error_handler(1, "Location: not convertible to node");
  #endif
      #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
      return geninfo<Node>::const_access(value);
      #else
      return 
        *boost::any_cast<Node>(&value);
      #endif
    }
    /*{\Mconversion converts |\Mvar| into a node.\\
       \precond |\Mvar| represents a node.}*/

    operator const Edge&() const
    { 
  #if !defined(CHECKING_OFF)
      if (type != EDGE) 
        CGAL_LEDA_SCOPE::error_handler(1, "Location: not convertible to edge");
  #endif
      #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
      return geninfo<Edge>::const_access(value);
      #else
      return 
        *boost::any_cast<Edge>(&value);
      #endif
    }
    /*{\Mconversion converts |\Mvar| into an edge.\\
       \precond |\Mvar| represents an edge.}*/

    GenericLocation<Node, Edge>& operator=(Node n)
    { clear(); init(n); return *this; }
    /*{\Mbinop makes |\Mvar| represent the node |n|.}*/

    GenericLocation<Node, Edge>& operator=(Edge e)
    { clear(); init(e); return *this; }
    /*{\Mbinop makes |\Mvar| represent the edge |e|.}*/

    Type get_type() const { return type; }
    /*{\Mop returns the type of the value contained in |\Mvar|.}*/

    bool is_nil()  const { return type == NIL; }
    /*{\Mop returns |true| iff |\Mvar| represents the value |nil|.}*/

    bool is_node() const { return type == NODE; }
    /*{\Mop returns |true| iff |\Mvar| represents a node.}*/

    bool is_edge() const { return type == EDGE; }
    /*{\Mop returns |true| iff |\Mvar| represents an edge.}*/

  private:
    void init() { type = NIL; }
    void init(Node n) 
    { type = NODE; 
      #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
      geninfo<Node>::create(value); 
      geninfo<Node>::access(value) = n;
      #else
      value=n;
      #endif
    }
    void init(Edge e) 
    { type = EDGE; 
      #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
      geninfo<Edge>::create(value); 
      geninfo<Edge>::access(value) = e;
      #else
      value=e;
      #endif
    }

    void clear()
    { 
      switch(type) {
      #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
       case NODE: geninfo<Node>::clear(value); break;
       case EDGE: geninfo<Edge>::clear(value); break;
      #else
       case NODE: value=boost::any(); break;
       case EDGE: value=boost::any(); break;
      #endif
       case NIL: break;
      }
    }

    void assign(const GenericLocation<Node, Edge>& L)
    { 
      type = L.type;
      switch(type) {
       case NODE: 
        #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
         geninfo<Node>::access(value) = geninfo<Node>::const_access(L.value);
        #else
         *boost_any_cast<Node>(&value) = boost::any_cast<Node>(L.value);
        #endif
         break;
       case EDGE: 
        #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
         geninfo<Edge>::access(value) = geninfo<Edge>::const_access(L.value);
        #else
         *boost::any_cast<Edge>(&value) = boost::any_cast<Edge>(L.value);
        #endif
         break;
       case NIL: break;
      }
    }

  public:
    typedef GenericLocation<Node, Edge> self;


private:
  Type   type;
  GenPtr value;
};

/*{\Mimplementation
  The data type |\Mtype| is implemented as a union of the types |Node| and 
  |Edge|. There is only constant time and space overhead.
}*/
template <class Node, class Edge>
inline
bool 
operator==(const GenericLocation<Node, Edge>& L1, 
           const GenericLocation<Node, Edge>& L2)
{ 
  if (L1.get_type() != L2.get_type()) return false;
  switch (L1.get_type()) {
   case GenericLocation<Node, Edge>::NIL:  return true;
   case GenericLocation<Node, Edge>::NODE: return Node(L1) == Node(L2);
   case GenericLocation<Node, Edge>::EDGE: return Edge(L1) == Edge(L2);
  }
}

template <class Node, class Edge>
inline
bool 
operator!=(const GenericLocation<Node, Edge>& L1, 
           const GenericLocation<Node, Edge>& L2)
{ return ! (L1==L2); }

template <class Node, class Edge>
std::ostream& operator<<(std::ostream& o, 
                         const GenericLocation<Node, Edge>& L)
{
  switch (L.get_type()) {
   case GenericLocation<Node, Edge>::NIL:  return o<<"nil";
   case GenericLocation<Node, Edge>::NODE: return o<<"node("<<&*Node(L)<<')';
   case GenericLocation<Node, Edge>::EDGE: return o<<"edge("<<&*Edge(L)<<')';
  }
  return o;
}

template <class Node, class Edge>
std::istream& operator>>(std::istream& i, GenericLocation<Node, Edge>&)
{ return i; }
template <class XCoord, class PredLessThanX, class Sweepline>
class GenericXStructure {
public:
  typedef std::vector<XCoord>    Array_Coordinates;
  typedef std::vector<Sweepline> Array_Sweeplines;
  typedef typename Array_Coordinates::const_iterator Coord_iterator;
  typedef typename Array_Sweeplines::const_iterator  Sweepline_iterator;
private:
  int stops;
  PredLessThanX     LtX;
  Array_Coordinates Coordinates;
  Array_Sweeplines  SweepLines;
  // SweepLines[0] is EmptyLine;
public:
  GenericXStructure() { clear(); }
  GenericXStructure(int n, const PredLessThanX& cmp) { init(n, cmp); }
  ~GenericXStructure() { clear(); }

  void init(int n, const PredLessThanX& cmp)
  { CGAL_NEF_TRACEN("XSinit "<<n); 
    LtX = cmp;
    Coordinates = Array_Coordinates(n);
    SweepLines =  Array_Sweeplines(2*n+1);
    stops = 0;
  }

  void clear()
  { Coordinates.clear();
    SweepLines.clear();
    stops = 0;
  }

  void insertLines(const XCoord& X, 
                   const Sweepline& atX, const Sweepline& inXplus)
  { CGAL_NEF_TRACEN("XSinsert "<<X); 
    Coordinates[stops]    = X;
    SweepLines[2*stops+1]   = atX;
    SweepLines[2*stops+2] = inXplus;
    ++stops;
  }

  Sweepline_iterator getLineAt(const XCoord& X) const
  { CGAL_NEF_TRACEN("XSgetLineAt "<<X);
    Sweepline_iterator sit = SweepLines.begin(); // EmptyLine
    if ( LtX(X,*Coordinates.begin()) ) {
      CGAL_NEF_TRACEN("infinity first");
      return sit; // ]-infinity, x0[
    }
    Coord_iterator stopit = std::lower_bound (
        Coordinates.begin(),Coordinates.end(),X,LtX);
    /* determines stopit maximal such that 
       \forall j \in [begin,stopit) : *j < X
       as a consequence now: *stopit >= X */
    bool found_exact = false;
    if ( LtX(X,*stopit) ) --stopit;  // X <  *stopit
    else found_exact = true;         // X >= *stopit
      
    CGAL_NEF_TRACEN("stopit "<<*stopit);
    int offset = stopit-Coordinates.begin();
    return found_exact ? 
      SweepLines.begin() + (2*offset+1) :
      SweepLines.begin() + (2*offset+2);
  }

  Sweepline_iterator begin() const { return SweepLines.begin(); }
  Sweepline_iterator end() const { return SweepLines.end();}

};
/*{\Manpage {PointLocator} {PLocTraits} {Planar Point Location} {PL}}*/

template <class PLocTraits>
class PointLocator {
/*{\Mdefinition
   An instance |\Mvar| of the parameterized data type |\Mtype| can be used
   to perform point location queries in the two-dimensional plane.
   Every non-empty instance |\Mvar| is associated with an embedded planar 
   graph |G|, which has to remain unchanged while it is referenced by |PL|.\\
   A location query for a point |p| returns the first object (node or edge)
   of |G| which is intersected by the straight ray starting in |p| and going
   vertically downwards/upwards.
   If the ray does not intersect any node or edge of |G|, then |nil| is 
   returned.\\
   The class |\Mtype| is generic, it is parameterized with a traits class 
   |PLocTraits| which widely controls its behaviour.
   The traits may even change the return type of a query and its semantics.
   There are predined traits classes for the LEDA graph types, which are
   described below in a seperate section.
}*/
public:
  // copied types from PLocTraits
  typedef typename PLocTraits::Point             Point;
  typedef typename PLocTraits::XCoord            XCoord;
  typedef typename PLocTraits::PredLessThanX     PredLessThanX;

  typedef typename PLocTraits::Graph             Graph;
  typedef typename PLocTraits::Node              Node;
  typedef typename PLocTraits::Edge              Edge;
  typedef typename PLocTraits::NodeIterator      NodeIterator;
  typedef typename PLocTraits::IncEdgeIterator   IncEdgeIterator;

  typedef typename PLocTraits::Curve             Curve;
  typedef typename PLocTraits::PredCompareCurves PredCompareCurves;

  typedef typename PLocTraits::QueryResult       QueryResult;

/*{\Mtypes}*/
  // define additional types 
  typedef GenericLocation<Node, Edge> Location;
  /*{\Mtypedef usual return value for the point loaction.}*/

  enum Direction { downwards, upwards};
  /*{\Menum used to specify the direction for the point location.}*/

  typedef CGAL_LEDA_SCOPE::pp_dictionary<Curve, Location, PredCompareCurves>   
                                                              Sweepline;
  typedef GenericXStructure<XCoord, PredLessThanX, Sweepline> XStructure;
  typedef typename Sweepline::item                            SL_item;
  typedef typename XStructure::Sweepline_iterator Sweepline_iterator;
  /*{\Mcreation}*/ 
  PointLocator() { clear(); }
  /*{\Mcreate creates an empty |\Mtype|.}*/

  PointLocator(const Graph& G, const PLocTraits& PLT = PLocTraits()) : 
   traits(PLT) { init(G); }
  /*{\Mcreate creates a |\Mtype| for the graph |G| and the traits |PLT|.}*/
     
  /*{\Moperations}*/
  void clear() { X_Structure.clear(); traits.clear(); }
  /*{\Mop makes |\Mvar| empty.}*/


  void init(const Graph& G, const PLocTraits& PLT) { traits = PLT; init(G); }
  /*{\Mop makes |\Mvar| a |\Mtype| for the graph |G| and the traits |PLT|.}*/

  void init(const Graph& G);
  /*{\Mop makes |\Mvar| a |\Mtype| for the graph |G|.}*/


  QueryResult locate(const Point& p, const Direction dir) const
  { return dir == downwards ? locate_down(p) : locate_up(p); }
  /*{\Mop locates the point |p| in the direction |dir|.}*/

  QueryResult locate_down(const Point& p) const;
  /*{\Mop locates the point |p| vertically downwards.}*/

  QueryResult locate_up(const Point& p) const;
  /*{\Mop locates the point |p| vertically upwards.}*/

  Location location(Sweepline_iterator S, SL_item it) const
  { return (it == nil ? Location() : S->inf(it)); }

  std::string str(const Sweepline& S) const
  { std::ostringstream os; os << "Sweepline:\n";
    SL_item it;
    forall_items(it,S) {  os << "  " << S.key(it) << std::endl; }
    return os.str();
  }


private:
  PLocTraits traits;
  XStructure X_Structure;

};

template <class PLocTraits>
void
PointLocator<PLocTraits>::init(const Graph& G)
{
  traits.sweep_begin(G);
  PredLessThanX LtX = traits.getLessThanX();
  typedef std::map<XCoord, std::list<Node>, PredLessThanX> dictionary;
  typedef typename dictionary::iterator dic_iterator;
  dictionary stops(LtX);
  // Note: X_Structure, Sweepline, and stops copy compare object

  NodeIterator ni = traits.Nodes_begin(G), beyond = traits.Nodes_end(G);
  for(; ni != beyond; ++ni) {
    XCoord currentX = traits.getXCoord(G, ni);
    stops[currentX].push_front(traits.toNode(ni));
  }

  Sweepline SL(traits.getCompareCurves());
  X_Structure.init(stops.size(), LtX);
  dic_iterator stop;
  for(stop = stops.begin(); stop != stops.end(); ++stop) {
    std::list<Node>& NodesOnSL = stop->second;
    traits.sweep_moveto(traits.getXCoord(G, *NodesOnSL.begin()));
      std::list<Edge> EmergingEdges, VerticalEdges;

      // explore the nodes on SL
      typename std::list<Node>::iterator cur_node;
      for(cur_node = NodesOnSL.begin(); 
          cur_node != NodesOnSL.end(); ++cur_node) {
        IncEdgeIterator ei     = traits.IncEdges_begin(G, *cur_node);
        IncEdgeIterator beyond = traits.IncEdges_end(G, *cur_node);
          CGAL_NEF_TRACEN("NODE: "<<(*cur_node)->point());
        for(; ei != beyond; ++ei) { 
          switch (traits.ClassifyEdge(G, traits.toEdge(ei), *cur_node)) {
            case PLocTraits::StartingNonVertical: 
              EmergingEdges.push_front(traits.toEdge(ei)); break;
            case PLocTraits::StartingVertical:    
              VerticalEdges.push_front(traits.toEdge(ei)); break;
            case PLocTraits::EndingNonVertical:
              SL.del(traits.makeCurve(G, traits.toEdge(ei))); break;
            case PLocTraits::EndingVertical: break;
          }
        }
      }

      // compute SL_at_X

      typename std::list<Edge>::iterator cur_edge;
      for(cur_edge=VerticalEdges.begin(); 
          cur_edge!=VerticalEdges.end(); ++cur_edge) 
        SL.insert(traits.makeCurve(G, *cur_edge), Location(*cur_edge));
      for(cur_node=NodesOnSL.begin();
          cur_node!=NodesOnSL.end(); ++cur_node)
        SL.insert(traits.makeCurve(G, *cur_node), Location(*cur_node));
      Sweepline SL_at_X = SL;

      // compute SL_in_X_plus

      for(cur_edge=VerticalEdges.begin(); 
          cur_edge!=VerticalEdges.end(); ++cur_edge)
        SL.del(traits.makeCurve(G, *cur_edge));
      for(cur_node=NodesOnSL.begin();
          cur_node!=NodesOnSL.end(); ++cur_node)
        SL.del(traits.makeCurve(G, *cur_node));

      for(cur_edge=EmergingEdges.begin();
          cur_edge!=EmergingEdges.end(); ++cur_edge) 
        SL.insert(traits.makeCurve(G, *cur_edge), Location(*cur_edge));

    X_Structure.insertLines(traits.getXCoord(G, *NodesOnSL.begin()), 
                            SL_at_X, SL);
  }

  traits.sweep_end();
}

template <class PLocTraits>
typename PointLocator<PLocTraits>::QueryResult
PointLocator<PLocTraits>::
locate_down(const typename PLocTraits::Point& p) const
{
  Sweepline_iterator line_at_x = X_Structure.getLineAt(traits.getXCoord(p)),
                     line_plus = line_at_x;
  CGAL_NEF_TRACEN("locate_down "<<str(*line_at_x));
  Curve p_curve = traits.makeCurve(p);
  PredCompareCurves cmp = traits.getCompareCurves();
  SL_item it = line_at_x->locate_pred(p_curve), it_plus(0);
  if ( it && line_at_x->inf(it).is_node() &&
       cmp(p_curve, line_at_x->key(it))!=0 ) {
    // p hit a feature exactly
    line_plus = line_at_x+1;
    if ( line_plus != X_Structure.end() )
      it_plus = line_plus->locate_pred(p_curve);
  }
  return traits.PostProcess(location(line_at_x,it),
                            location(line_plus,it_plus),p);
}

template <class PLocTraits>
typename PointLocator<PLocTraits>::QueryResult
PointLocator<PLocTraits>::locate_up(const typename PLocTraits::Point& p) const
{
  Sweepline_iterator line_at_x = 
    X_Structure.getLineAt(traits.getXCoord(p)), line_plus;
  Curve p_curve = traits.makeCurve(p);
  PredCompareCurves cmp = traits.getCompareCurves();
  SL_item it = line_at_x->locate_succ(p_curve), it_plus(0);
  if ( it && line_at_x->inf(it).is_node() &&
       cmp(p_curve, line_at_x->key(it))!=0 ) {
    // p hit a feature exactly
    line_plus = line_at_x+1;
    if ( line_plus != X_Structure.end() )
      it_plus = line_plus->locate_succ(p_curve);
  }
  return traits.PostProcess(location(line_at_x,it),
                            location(line_plus,it_plus), p);
}

/*{\Mimplementation
  The implementation of the data type |\Mtype| is based on partially 
  persistent binary search trees.
  The expected space requirement is $O(k)$ where $k$ is the sum of the number 
  of nodes and the number of edges in the graph $G$.
  The expected time needed for construction and the operation |init| is 
  $O(k \cdot \log k)$, for the |locate|-operations it is $O(\log k)$. The
  operation |clear| runs in $O(k)$.
}*/
/*{\Mtext
  \headerline{\arabic{manctr}. Predefined traits classes}
  \stepcounter{manctr}
  All predefined traits classes have in common that the return type of a query 
  is the type |Location|.
  The embedding of the given graph |G| is a straight-line embedding, so that
  it is totally determined by the position of the nodes of |G|.
  Such a position is specified by a |Point| which can be one of the LEDA point
  types |point| or |rat_point|. The positions can be specified implicitly by 
  the node attribute of a parameterized graph (e.g. |GRAPH<Point,...>|) or 
  explicitly by a |node_array<Point>|. In case of explicit specification a
  |node_array| with the positions of the nodes can be passed to the constructor
  of the traits class.
  Further, the point location processes for maps and for standard graphs differ
  slightly. As a map is a bidirected graph where each edge knows its reversal,
  the traits classes for maps can ensure the following property:
  If the result of a query for point |p| is an edge |e| (not containing |p|),
  then |e| bounds the face of |G| which contains |p|, i.e. |p| lies to the
  left of |e|.\\
  Here comes a list of the predefined traits classes:\\[-5.5ex]
  \begin{itemize}
  \item |PLocTraits<Graph>|: standard traits for implicitly specified node 
    positions\\
    |Graph| can be |GRAPH<Point,...>| (standard graph) or 
    |PLANAR_MAP<Point,...,...>| (map).
  \item |PLocTraits_NodeArray<Graph,Point>|: std. traits for explicitly 
    specified node positions\\
    |Graph| can be |graph| (standard graph) or |planar_map| (map).
  \item |PLocTraits_Map<Graph>| and |PLocTraits_Map_NodeArray<Graph,Point>|:\\
    The parameter |Graph| can be |GRAPH<Point,...>| and |graph| respectively.
    These traits classes assume that the given graphs are maps.
  \item |PLocTraits< GRAPH<Circle,Point> >|: traits class for closest-site 
    voronoi diagrams
  \end{itemize}
  Note that a traits class instantiated with |Graph| can also handle graph 
  types which are derived from |Graph|. Thus |PLocTraits< graph<Point,T> >| 
  can be used for graphs of type |ugraph<Point,T>| for example.
}*/
/*{\Mexample
First we show an example where the node positions are given implicitly
as node attributes:
\begin{verbatim}
  typedef PointLocator< PLocTraits< GRAPH<Point, int> > > PLocator1;
  typedef PLocator1::Location Location;

  UGRAPH<Point, int> G;
  ... // construct G

  PLocator1 PL1(G);
  Point p = ...; // compute p
  Location L1 = PL1.locate_down(p);
\end{verbatim}

The second example shows how a |node_array| can be used to determine the
node positions:
\begin{verbatim}
  typedef PLocTraits_NodeArray<planar_map,Point> PLocTraits2;
  typedef PointLocator<PLocTraits2> PLocator2;

  planar_map pm;
  node_array<Point> na;
  ... // construct pm and na

  PLocator2 PL2(pm, PLocTraits2(na));
  Point q = ...; // compute q
  Location L2 = PL2.locate_up(q);
\end{verbatim}
}*/

#endif // GEN_POINT_LOCATION_H
