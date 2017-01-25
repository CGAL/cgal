// Copyright (c) 2005 INRIA Sophia-Antipolis (France).
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
// Author(s)     : Nico Kruithof <Nico@cs.rug.nl>
//                 Developed at Rijksuniversiteit Groningen (Netherlands)

#ifndef CGAL_TRIANGULATION_SIMPLEX_3_H
#define CGAL_TRIANGULATION_SIMPLEX_3_H

#include <CGAL/license/Triangulation_3.h>


#include <CGAL/assertions.h>
#include <algorithm>

namespace CGAL {

template < class TriangulationDataStructure_3 >
class Triangulation_simplex_3 {
  typedef TriangulationDataStructure_3         TDS;
  typedef Triangulation_simplex_3<TDS>           Self;
public:
  typedef Self                                 Simplex;

  typedef typename TDS::Vertex_handle            Vertex_handle;
  typedef typename TDS::Edge                     Edge;
  typedef typename TDS::Facet                    Facet;
  typedef typename TDS::Cell_handle              Cell_handle;

  typedef typename TDS::Cell_circulator          Cell_circulator;
  typedef typename TDS::Facet_circulator         Facet_circulator;

  typedef typename TDS::Edge_iterator            Edge_iterator;
  typedef typename TDS::Facet_iterator           Facet_iterator;

  // Constructors

  // Default constructor initialises to undefined simplex:
  Triangulation_simplex_3() : ref(-1), ch() { }
	
  Triangulation_simplex_3(Vertex_handle vh) {
    set_vertex(vh);
  }
  Triangulation_simplex_3(const Edge &e) {
    set_edge(e);
  }
  Triangulation_simplex_3(const Facet &f) {
    set_facet(f);
  }
  Triangulation_simplex_3(Cell_handle ch_) {
    set_cell(ch_);
  }

  Triangulation_simplex_3(Cell_circulator ccir) {
    set_cell(ccir);
  }
  Triangulation_simplex_3(Facet_circulator fcir) {
    set_facet(*fcir);
  }

  Triangulation_simplex_3(Edge_iterator eit) {
    set_edge(*eit);
  }
  Triangulation_simplex_3(Facet_iterator fit) {
    set_facet(*fit);
  }

  // Conversions:
  operator Vertex_handle () const
  {
    CGAL_assertion(dimension() == 0);
    return ch->vertex(index(0));
  }
  operator Edge () const
  {
    CGAL_assertion(dimension() == 1);
    return Edge(ch,index(0),index(1));
  }
  operator Facet () const
  {
    CGAL_assertion(dimension() == 2);
    return Facet(ch,index(0));
  }
  operator Cell_handle () const
  {
    CGAL_assertion(dimension() == 3);
    return ch;
  }

  // returns the dimension of the simplex
  int dimension () const {
    return (ref & 3);
  }
  // returns an incident cell:
  Cell_handle incident_cell() {
    return ch;
  }

  template < class TDS2 >
  friend bool operator==(Triangulation_simplex_3<TDS2> s0,
			 Triangulation_simplex_3<TDS2> s1);
  template < class TDS2 >
  friend bool operator< (Triangulation_simplex_3<TDS2> s0,
			 Triangulation_simplex_3<TDS2> s1);
	
private:
  void set_vertex(const Vertex_handle vh) {
    ch = vh->cell();
    ref = (ch->index(vh) << 2); /* dim == 0 */
    CGAL_assertion (ch != Cell_handle());
  }
  void set_edge(const Edge &e) {
    ch = e.first;
    ref = (((e.third<< 2) + e.second) << 2) + 1; /* dim */
    CGAL_assertion (ch != Cell_handle());
  }
  void set_facet(const Facet &f) {
    ch = f.first;
    ref = (f.second << 2) + 2; /* dim */
    CGAL_assertion (ch != Cell_handle());
  }
  void set_cell(Cell_handle ch_) {
    ch = ch_;
    ref = 3; /* dim */
    CGAL_assertion (ch != Cell_handle());
  }

  inline int index(int i) const {
    CGAL_assertion (i==0 || ((i==1) && (dimension()==1)));
    return (ref >> (2*(i+1))) & 3;
  }

  int ref; // storage iijjdd (index i, index j, dimension of simplex)
  Cell_handle ch; // Corresponding cell handle
};

///////////////////////////////
// Simplex functions
///////////////////////////////
template < class TriangulationDataStructure_3 >
bool
operator!=(Triangulation_simplex_3<TriangulationDataStructure_3> s0,
	   Triangulation_simplex_3<TriangulationDataStructure_3> s1) {
  return !(s0==s1);
}

template < class TriangulationDataStructure_3 >
bool
operator==(Triangulation_simplex_3<TriangulationDataStructure_3> s0,
	   Triangulation_simplex_3<TriangulationDataStructure_3> s1) {
  typedef Triangulation_simplex_3<TriangulationDataStructure_3> Sim;
  if (s0.dimension() != s1.dimension()) return false;
	
  typename Sim::Cell_handle neighbor;
	
  switch (s0.dimension()) {
  case (0): // Vertex
    return (s0.ch->vertex(s0.index(0)) == s1.ch->vertex(s1.index(0)));
  case (1): // Edge
    return ((s0.ch->vertex(s0.index(0)) == s1.ch->vertex(s1.index(0)) &&
	     s0.ch->vertex(s0.index(1)) == s1.ch->vertex(s1.index(1))) ||
	    (s0.ch->vertex(s0.index(1)) == s1.ch->vertex(s1.index(0)) &&
	     s0.ch->vertex(s0.index(0)) == s1.ch->vertex(s1.index(1))));
  case (2):
    if (s0.ch == s1.ch && s0.index(0) == s1.index(0)) {
      return true;
    }
			
    neighbor = s0.ch->neighbor(s0.index(0));
    if (neighbor == s1.ch &&
	neighbor->index(s0.ch) == s1.index(0)) {
      return true;
    }
    return false;
  case (3):
    return (&(*s0.ch) == &(*s1.ch));
  }
  CGAL_error();
  return false;
}

template < class TriangulationDataStructure_3 >
bool
operator<(Triangulation_simplex_3<TriangulationDataStructure_3> s0,
	  Triangulation_simplex_3<TriangulationDataStructure_3> s1) {
  typedef Triangulation_simplex_3<TriangulationDataStructure_3> Sim;

  if (s0 == s1) return false;
  if (s0.dimension() < s1.dimension()) return true;
  if (s0.dimension() > s1.dimension()) return false;
	
  // Dimensions are equal, compare the memory addresses of the simplices
  typename Sim::Cell_handle ch1, ch2;
  typename Sim::Vertex_handle vh1, vh2, vh3, vh4;
  switch (s0.dimension()) {
  case (0): // Vertex
    // Vertextices are not equal
    return (&(*s0.ch->vertex(s0.index(0))) <
	    &(*s1.ch->vertex(s1.index(0))));
  case (1): // Edge
    vh1 = s0.ch->vertex(s0.index(0));
    vh2 = s0.ch->vertex(s0.index(1));
    vh3 = s1.ch->vertex(s1.index(0));
    vh4 = s1.ch->vertex(s1.index(1));
			
    if ((std::min)(&(*vh1), &(*vh2)) < (std::min)(&(*vh3), &(*vh4)))
      return true;
			
    if ((std::min)(&(*vh1), &(*vh2)) > (std::min)(&(*vh3), &(*vh4)))
      return false;
			
    if ((std::max)(&(*vh1), &(*vh2)) < (std::max)(&(*vh3), &(*vh4)))
      return true;
			
    return false;
  case (2): // Facet
    ch1 = s0.ch->neighbor(s0.index(0));
    ch2 = s1.ch->neighbor(s1.index(0));
			
    if ((std::min)(&(*s0.ch), &(*ch1)) < (std::min)(&(*s1.ch), &(*ch2)))
      return true;
			
    if ((std::min)(&(*s0.ch), &(*ch1)) > (std::min)(&(*s1.ch), &(*ch2)))
      return false;
			
    if ((std::max)(&(*s0.ch), &(*ch1)) < (std::max)(&(*s1.ch), &(*ch2)))
      return true;
			
    return false;
  case (3): // Cell
    return (&(*s0.ch) < &(*s1.ch));
  }
  CGAL_error();
  return false;
}

template < class TriangulationDataStructure_3 >
std::ostream &
operator<< (std::ostream& os,
	    const Triangulation_simplex_3<TriangulationDataStructure_3> &s)
{
  typename TriangulationDataStructure_3::Vertex_handle vh;
  typename TriangulationDataStructure_3::Edge e;
  typename TriangulationDataStructure_3::Facet f;
  typename TriangulationDataStructure_3::Cell_handle ch;
  switch (s.dimension()) {
    case 0:
      vh = s;
      os << &*vh;
      break;
    case 1:
      e = s;
      os << &*(e.first->vertex(e.second)) << " "
	 << &*(e.first->vertex(e.third));
      break;
    case 2:
      f = s;
      os << &*(f.first->vertex((f.second+1)&3)) << " "
	 << &*(f.first->vertex((f.second+2)&3)) << " "
	 << &*(f.first->vertex((f.second+3)&3));
      break;
    case 3:
      ch = s;
      os << &*(ch->vertex(0)) << " "
	 << &*(ch->vertex(1)) << " "
	 << &*(ch->vertex(2)) << " "
	 << &*(ch->vertex(3));
      break;
  }
  return os;
}


} //namespace CGAL

#endif // CGAL_TRIANGULATION_SIMPLEX_3_H
