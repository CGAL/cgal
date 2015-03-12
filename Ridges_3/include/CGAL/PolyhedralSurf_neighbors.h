// Copyright (c) 2007  INRIA Sophia-Antipolis (France), INRIA Lorraine LORIA.
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
// Author(s)     : Marc Pouget and Frédéric Cazals
#ifndef CGAL_POLYHEDRALSURF_NEIGHBORS_H_
#define CGAL_POLYHEDRALSURF_NEIGHBORS_H_

#include <queue>
#include <algorithm>
#include <list>
#include <CGAL/basic.h>

namespace CGAL {

//---------------------------------------------------------------------------
//T_Gate : element of the priority queue. A gate is a halfedge and a
//number giving the max distance from v to the vertices of the
//triangle incident to the halfedge.
//---------------------------------------------------------------------------
template < class TriangleMesh > class T_Gate
{
  typedef typename boost::property_map<TriangleMesh,CGAL::vertex_point_t>::type VPM;
  typedef typename boost::property_traits<VPM>::value_type Point_3;
  typedef typename Kernel_traits<Point_3>::Kernel Kernel;

public:
  typedef typename Kernel::FT       FT;
  typedef typename Kernel::Vector_3 Vector_3;

  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor  halfedge_descriptor;
 
  T_Gate(FT d, const halfedge_descriptor he);
  FT& d() { return m_d;}
  const FT d() const { return m_d;}            
  const halfedge_descriptor he() { return m_he;}

private:
  FT m_d;
  halfedge_descriptor m_he;
};

//////////////IMPLEMENTATION//////////////////////////
template < class TriangleMesh > 
T_Gate<TriangleMesh>::T_Gate(FT d, 
					    const halfedge_descriptor he)
  : m_d(d), m_he(he)
{}

//---------------------------------------------------------------------------
// functor for priority queue
// order so that the top element is the smallest in the queue
//---------------------------------------------------------------------------
template<class g>
struct compare_gates 
{       
        bool operator()(const g& g1, 
                        const g& g2) const
        {       
                return g1.d() > g2.d();
        }
};

//---------------------------------------------------------------------------
//T_PolyhedralSurf_neighbors : MAIN class for computation, it uses the
//class Gate and the functor compare_gates for the definition of a
//priority queue
//---------------------------------------------------------------------------
template < class TriangleMesh >
class T_PolyhedralSurf_neighbors
{

  typedef typename boost::property_map<TriangleMesh,CGAL::vertex_point_t>::const_type VPM;
  typedef typename boost::property_traits<VPM>::value_type Point_3;
  typedef typename Kernel_traits<Point_3>::Kernel Kernel;

public:
  typedef typename Kernel::FT        FT;
  typedef typename Kernel::Vector_3  Vector_3;

  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor     vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor   halfedge_descriptor;
  typedef CGAL::Halfedge_around_target_circulator<TriangleMesh>
  Halfedge_around_vertex_const_circulator;
  typedef typename boost::graph_traits<TriangleMesh>::vertex_iterator   Vertex_const_iterator;
  typedef T_Gate<TriangleMesh> Gate;

  T_PolyhedralSurf_neighbors(const TriangleMesh& P);
  // vertex_neigh stores the vertex v and its 1Ring neighbors contour
  // stores halfedges, oriented CW, following the 1Ring disk border
  // OneRingSize is the max distance from v to its OneRing
  // neighbors. (the tag is_visited is not mofified)
  void compute_one_ring(const vertex_descriptor v,
			std::vector<vertex_descriptor> &vertex_neigh,
			std::list<halfedge_descriptor> &contour,
			FT &OneRingSize);
  // call compute_one_ring and expand the contour (circle of halfedges
  // CW), vertex_neigh are vertices on and inside the contour (there
  // tag is_visited is set to true, but reset to false at the end),
  // size is such that gates with distance less than size*OneRingSize
  // are processed
  void compute_neighbors(const vertex_descriptor v,
			 std::vector<vertex_descriptor> &vertex_neigh,
			 std::list<halfedge_descriptor> &contour,
			 const FT size); 
  //vertex tags is_visited are set to false
  void reset_is_visited_map(std::vector<vertex_descriptor> &vces);


  Gate make_gate(const vertex_descriptor v,
                 const halfedge_descriptor he)
  {
    Point_3 p0 = get(vpm, v),
      p1 = get(vpm, target(he,P)),
      p2 = get(vpm, target(next(he,P),P)),
      p3 = get(vpm, target(prev(he,P),P));
    Vector_3 p0p1 = p0 - p1,
      p0p2 = p0 - p2,
      p0p3 = p0 - p3;
    FT d1 = p0p1*p0p1,
      d2 = p0p2*p0p2,
      d3 = p0p3*p0p3;
    FT d = CGAL::sqrt( (std::max)( (std::max)(d1,d2), d3) );
    return Gate(d,he);
  }

 protected:
  /*
  //tag to visit vertices
  struct Vertex_cmp{//comparison is wrt vertex addresses
    bool operator()(const vertex_descriptor a, const vertex_descriptor b) const{
      return &*a < &*b;
    }
  };
  */
  const TriangleMesh& P;
  VPM vpm;
  typedef std::map<vertex_descriptor, bool/*, Vertex_cmp*/> Vertex2bool_map;
  Vertex2bool_map is_visited_map;
};

//////////////IMPLEMENTATION//////////////////////////
template < class TriangleMesh >
T_PolyhedralSurf_neighbors < TriangleMesh >::
T_PolyhedralSurf_neighbors(const TriangleMesh& P)
  :P(P), vpm(get(vertex_point,P))
{
  //init the is_visited_map
  Vertex_const_iterator itb, ite;
  boost::tie(itb,ite) = vertices(P);
  for(;itb!=ite;itb++) is_visited_map[*itb] = false; 
}

template < class TriangleMesh >
void T_PolyhedralSurf_neighbors < TriangleMesh >::
compute_one_ring(const vertex_descriptor v,
		 std::vector<vertex_descriptor> &vertex_neigh,
		 std::list<halfedge_descriptor> &contour,
		 FT &OneRingSize)
{
  vertex_neigh.push_back(v);
  Halfedge_around_vertex_const_circulator he_circ(halfedge(v,P),P),
                                    he_end = he_circ;
  do {
    if ( is_border(*he_circ,P) )//then he and he->next follow the contour CW
	{contour.push_back(*he_circ);
          contour.push_back(next(*he_circ,P));}
    else contour.push_back(opposite(prev(*he_circ,P),P));//not border, he->prev->opp on contour CW
    vertex_neigh.push_back(target(opposite(*he_circ,P),P));
      he_circ++;
  } while (he_circ != he_end);

  //compute OneRingSize = distance(v, 1Ring)
  OneRingSize = 0;
  typename std::vector<vertex_descriptor>::const_iterator itb = vertex_neigh.begin(),
    ite = vertex_neigh.end();
  itb++;//the first vertex v is the center to which distances are
	//computed from, for other 1ring neighbors
  Point_3 p0 = get(vpm, v), p;
  Vector_3 p0p;
  FT d = OneRingSize;
  for (; itb != ite; itb++){

    p = get(vpm, *itb);
    p0p = p0 - p;
    d =  CGAL::sqrt(p0p*p0p);
    if (d > OneRingSize) OneRingSize = d;
  }
}

template < class TriangleMesh >
void T_PolyhedralSurf_neighbors < TriangleMesh >::
compute_neighbors(const vertex_descriptor v,
		  std::vector<vertex_descriptor> &vertex_neigh,
		  std::list<halfedge_descriptor> &contour,
		  const FT size)  
{
  FT OneRingSize;
  compute_one_ring(v, vertex_neigh, contour, OneRingSize);
  const FT d_max = OneRingSize*size;
  std::priority_queue< Gate, std::vector< Gate >, compare_gates< Gate > > GatePQ;
  // tag neighbors 
  typename std::vector<vertex_descriptor>::const_iterator itbv = vertex_neigh.begin(),
    itev = vertex_neigh.end();
  for (; itbv != itev; itbv++) is_visited_map.find(*itbv)->second = true;

  // init GatePQ
  typename std::list<halfedge_descriptor>::const_iterator itb = contour.begin(),
                                       ite = contour.end();
  for (; itb != ite; itb++) {
    if (!( is_border(*itb,P) )) GatePQ.push(make_gate(v, *itb));
  }
  // init d_current
  Gate firstGate = GatePQ.top();
  FT d_current = firstGate.d();
  // main loop
  while ( !GatePQ.empty() && d_current <= d_max ) {
    Gate gate = GatePQ.top();
    GatePQ.pop();
    d_current = gate.d();
    halfedge_descriptor he = gate.he(), he1, he2;
    vertex_descriptor v1;
    // find the gate on the contour
    typename std::list<halfedge_descriptor>::iterator pos_he, pos_prev, pos_next, iter;
   
    pos_he = find(contour.begin(), contour.end(), he);
    iter = pos_he;
    /**
       there are different cases to expand the contour : 
       (case 3) he is not on the contour, nothing to do
       (case 2) he is on the contour and either the previous or the next 
       following edge in the triangle is also on the contour, then delete
       these 2 he from the contour and add the third one to the contour 
       and the PQ.
       (case1) the vertex opposite to he is not visited, then the he is removed
       from the contour, the two others are added to the contour and PQ, the 
       vertex is set visited.
    */

    // if the gate is not encountered on the contour (case 3)
    if ( pos_he == contour.end() ) continue;
    // simulate a circulator on the contour: 
    // find the prev and next pos on coutour
    if ( ite != (++iter) ) pos_next = iter;
    else pos_next = contour.begin();
    iter = pos_he;
    if ( iter != contour.begin() ) pos_prev = --iter;
    else pos_prev = --contour.end();

    if ( next(he,P) == *pos_next )
      {  // case 2a
	//contour
	he1 = opposite(prev(he,P),P);
	contour.insert(pos_he, he1);
	contour.erase(pos_he);
	contour.erase(pos_next);
	//GatePQ
	if ( !is_border(he1,P) ) GatePQ.push(make_gate(v, he1));
	continue;
      }
    else if ( prev(he,P) == *pos_prev )
      {  // case 2b
	//contour
	he1 = opposite(next(he,P),P);
	contour.insert(pos_prev, he1);
	contour.erase(pos_prev);
	contour.erase(pos_he);
	//GatePQ
	if ( ! is_border(he1,P) ) GatePQ.push(make_gate(v, he1));
	continue;
      }
    v1 = target(next(he,P),P);
    if ( !is_visited_map.find(v1)->second )
      {  // case 1
	//vertex
	is_visited_map.find(v1)->second = true;
	vertex_neigh.push_back(v1);
	//contour
	he1 = opposite(prev(he,P),P);
	he2 = opposite(next(he,P),P);
	contour.insert(pos_he, he1);
	contour.insert(pos_he, he2);
	contour.erase(pos_he);
	//GatePQ
	if ( ! is_border(he1,P) ) GatePQ.push(make_gate(v, he1));
	if ( ! is_border(he2,P) ) GatePQ.push(make_gate(v, he2));
	continue;
      }
    //else do nothing (keep the he on the contour, and continue) to
    //prevent a change of the topology.
  }// end while
  
  reset_is_visited_map(vertex_neigh);
}

template < class TriangleMesh >
void T_PolyhedralSurf_neighbors < TriangleMesh >::
reset_is_visited_map(std::vector<vertex_descriptor> &vces)
{
  typename std::vector<vertex_descriptor>::const_iterator 
    itb = vces.begin(), ite = vces.end();
  for (;itb != ite; itb++) is_visited_map[*itb] = false;
}

} //namespace CGAL

#endif
