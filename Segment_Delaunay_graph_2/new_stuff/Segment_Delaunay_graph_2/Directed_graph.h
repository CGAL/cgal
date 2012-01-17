// Copyright (c) 2006 Foundation for Research and Technology-Hellas (Greece).
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
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>

#ifndef CGAL_DIRECTED_GRAPH_H
#define CGAL_DIRECTED_GRAPH_H 1

#include <CGAL/basic.h>

namespace CGAL {

template<typename Data_t = int>
class Directed_graph_node
{
public:
  typedef Data_t                               Data;

private:
  typedef Directed_graph_node<Data>            Self;

public:
  typedef Self*                                Node_handle;
  typedef std::set<Node_handle>                Container;
  typedef typename Container::const_iterator   neighbor_iterator;
  typedef typename Container::size_type        size_type;

private:
  Data data_;
  unsigned long id_;
  Container out_neighbors_;
  Container in_neighbors_;

public:
  Directed_graph_node(const Data& data = Data(),
		      unsigned long id = 0)
    : data_(data), id_(id) {}

public:
  bool add_out_neighbor(Node_handle out_neighbor)
  {
    std::pair<neighbor_iterator,bool>
      res = out_neighbors_.insert(out_neighbor);
    return res.second;
  }

  bool add_in_neighbor(Node_handle in_neighbor)
  {
    std::pair<neighbor_iterator,bool>
      res = in_neighbors_.insert(in_neighbor);
    return res.second;
  }

  inline neighbor_iterator out_neighbors_begin() const {
    return out_neighbors_.begin();
  }

  inline neighbor_iterator out_neighbors_end() const {
    return out_neighbors_.end();
  }

  inline neighbor_iterator in_neighbors_begin() const {
    return in_neighbors_.begin();
  }

  inline neighbor_iterator in_neighbors_end() const {
    return in_neighbors_.end();
  }

  //  unsigned long id() const { return id_; }

  inline size_type out_degree() const { return out_neighbors_.size(); }
  inline size_type in_degree() const { return in_neighbors_.size(); }

  const Data& data() const { return data_; }

  static
  Node_handle make_node(const Data& data = Data(), unsigned long id = 0) {
    return new Self(data, id);
  }

  static void delete_node(Node_handle n) {
    delete n;
  }
};


//=====================================================================
//=====================================================================


template<typename Data_t = void>
class Directed_graph
{
public:
  typedef Data_t  Data;
  typedef Directed_graph_node<Data>            Node;

  typedef typename Node::Node_handle           Node_handle;
  typedef typename Node::size_type             size_type;
  typedef std::set<Node_handle>                Container;
  typedef typename Container::const_iterator   Node_iterator;

public:
  Directed_graph() {}
  ~Directed_graph() { clear(); }

  Node_handle add_node(const Data& data) {
#if 1
    Node_handle np = Node::make_node(data);
    nodes_.insert(np);
    return np;
#else
    Node_handle np = Node::make_node(data);
    Node_iterator np_it = nodes_.find(np);
    if ( np_it == nodes_.end() ) {
      std::pair<Node_iterator,bool> it = nodes_.insert(np);
      return *(it.first);
    } else {
      delete np;
      return *np_it;
    }
#endif
  }

  void add_directed_edge(const Node_handle& h1, const Node_handle& h2)
  {
    h1->add_out_neighbor(h2);
    h2->add_in_neighbor(h1);
  }

  inline Node_iterator nodes_begin() const {
    return nodes_.begin();
  }

  inline Node_iterator nodes_end() const {
    return nodes_.end();
  }

  void clear()
  {
    for(Node_iterator it = nodes_begin(); it != nodes_end(); ++it) {
      Node::delete_node(*it);
    }
  }

private:
  Container nodes_;
};




} //namespace CGAL

#endif // CGAL_DIRECTED_GRAPH_H
