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
#ifdef USE_GRAPH
#include <CGAL/Node.h>
#include <CGAL/Graph.h>
#endif

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
#ifdef USE_GRAPH
    typedef std::list<Facet> Facets;
    typedef typename Facets::iterator Facets_iterator;

    typedef CGAL::Node<Facet> Node;
    typedef CGAL::Graph<Facet> Graph;
    typedef std::map<Facet, Node*> Nodes_map;
    typedef typename Nodes_map::iterator Nodes_map_iterator;
#endif

  private:
    bool visited;
#ifdef USE_GRAPH
    Graph umbrellas_dual;
#endif
    public: // AF: todo: make private and wrap in functions
      bool regular_is_cached;
      bool regular;

  public:  
    // Constructors

    Complex_2_in_triangulation_vertex_base_3()
      : Vb(), visited(false), regular_is_cached(false), regular(false)
      {}

#ifdef USE_GRAPH
  Graph& get_umbrellas_dual() {
    return umbrellas_dual;
  }
#endif

    bool is_visited() const {
      return visited;
    }
    
    void set_visited(const bool b) {
      visited = b;
    }

#ifdef USE_GRAPH
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
#endif

  };  // end Complex_2_in_triangulation_vertex_base_3


}  // namespace CGAL


#endif  // CGAL_COMPLEX_2_IN_TRIANGULATION_CELL_BASE_3_H
