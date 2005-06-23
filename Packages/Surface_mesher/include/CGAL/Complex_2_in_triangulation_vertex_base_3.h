// Copyright (c) 2003-2005  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source: 
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Steve OUDOT


#ifndef CGAL_COMPLEX_2_IN_TRIANGULATION_VERTEX_BASE_3_H
#define CGAL_COMPLEX_2_IN_TRIANGULATION_VERTEX_BASE_3_H

#include <CGAL/Triangulation_vertex_base_3.h>
#include <CGAL/Node.h>
#include <CGAL/Graph.h>

namespace CGAL {

  template < class GT, class Vb = Triangulation_vertex_base_3 <GT> > 
    class Complex_2_in_triangulation_vertex_base_3 : public Vb {    
    
  public:
    typedef Complex_2_in_triangulation_vertex_base_3 <GT, Vb> Self;
    
    template < class TDS3 >
    struct Rebind_TDS {
      typedef typename Vb::template Rebind_TDS<TDS3>::Other  Vb3;
      typedef Complex_2_in_triangulation_vertex_base_3 <GT, Vb3> Other;
    };
    
    typedef typename Vb::Triangulation_data_structure Tds;
    typedef typename Tds::Vertex_handle Vertex_handle;
    typedef typename Tds::Cell_handle Cell_handle;
    typedef typename Tds::Facet Facet;
    typedef std::list<Facet> Facets;
    typedef typename Facets::iterator Facets_iterator;
    typedef CGAL::Node<Facet> Node;
    typedef CGAL::Graph<Facet> Graph;
    typedef std::map<Facet, Node*> Nodes_map;
    typedef typename Nodes_map::iterator Nodes_map_iterator;


  private:
    bool visited;
    Graph umbrellas_dual;

  private:
    // computes and return the semi-facet with the smallest cell_handle
    Facet facet_with_smallest_cell_handle(const Facet& f) const {
      Cell_handle c = f.first;
      int i = f.second;
      Cell_handle c2 = c->neighbor(i);
      int i2 = c2->index(c);
      
      Cell_handle cmin = c;
      int imin = i;
      
      if (c2 < cmin) {
	cmin = c2;
	imin = i2;
      }
      
      return std::make_pair(cmin, imin);
    }

  public:  
    // Constructors

    Complex_2_in_triangulation_vertex_base_3() : Vb() {
      
      visited = false;
    }
    
  Graph& get_umbrellas_dual() {
    return umbrellas_dual;
  }

    bool is_visited() const {
      return visited;
    }
    
    void set_visited(const bool b) {
      visited = b;
    }

    // parameters are : a facet to be added as a vertex of the graph
    // the list of incident facets of this facet in the complex
    void add_in_graph(const Facet& f, const Facets& lof) {
      umbrellas_dual.add_node(f, lof);
    }

    void remove_from_graph(const Facet& f) {
      umbrellas_dual.remove_node(f);
    }

    bool is_graph_connected() {
      return ( umbrellas_dual.is_graph_connected ());
    }
  
  };  // end Complex_2_in_triangulation_vertex_base_3


}  // namespace CGAL


#endif  // CGAL_COMPLEX_2_IN_TRIANGULATION_CELL_BASE_3_H
