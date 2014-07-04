// Copyright (c) 2013 CNRS and LIRIS' Establishments (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s) : Pierre Talbot

#ifndef CGAL_BOOST_GRAPH_GRAPH_TRAITS_CMAP_H
#define CGAL_BOOST_GRAPH_GRAPH_TRAITS_CMAP_H

#include <utility>
#include <iterator>

#include <boost/config.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

#include <CGAL/Combinatorial_map.h>
#include <CGAL/Dart_iterators.h>

#define CGAL_CMAP_BASE_TEMPLATE_ARGS template<unsigned int d, class Refs, class Items, class Alloc>
#define CGAL_CMAP_BASE_TYPE CGAL::Combinatorial_map_base<d, Refs, Items, Alloc>

#define CGAL_CMAP_TEMPLATE_ARGS template<unsigned int d, class Items, class Alloc>
#define CGAL_CMAP_TYPE CGAL::Combinatorial_map<d, Items, Alloc>

#define CGAL_LCC_TEMPLATE_ARGS template < unsigned int d_, unsigned int ambient_dim, \
             class Traits_, \
             class Items_, \
             class Alloc_, \
             template<unsigned int, class,class,class>\
             class CMap>

#define CGAL_LCC_TYPE CGAL::Linear_cell_complex<d_, ambient_dim, Traits_, Items_, Alloc_, CMap> 

namespace CGAL {


template <class CMap, typename Dart_Iterator>
class CMap_dart_handle_iterator 
{
public:
  typedef Dart_Iterator Iterator;    

  typedef typename CMap::Dart_handle Dart_handle;

  typedef CMap_dart_handle_iterator<CMap, Dart_Iterator> Self; 

  typedef typename std::iterator_traits<Iterator>::iterator_category iterator_category;
  typedef typename std::iterator_traits<Iterator>::difference_type   difference_type;
  typedef Dart_handle                                                value_type;
  typedef value_type                                                 reference;
  typedef value_type                                                 pointer;

public:

// OPERATIONS Forward Category
// ---------------------------

  bool operator==( const Self& i) const { return ( nt == i.nt); }
  bool operator!=( const Self& i) const { return !(nt == i.nt );}
  value_type operator*() const { return nt; }
  value_type operator->() { return nt; }

  Self& operator++() {
    ++nt;
    return *this;
  }

  Self  operator++(int) {
    Self tmp = *this;
    ++*this;
    return tmp;
  }

  CMap_dart_handle_iterator(Iterator iter):
    nt(iter)
  {}

  // Default constructor
  CMap_dart_handle_iterator():
    nt(get_default(), Dart_handle())
  {}

  CMap_dart_handle_iterator(const CMap_dart_handle_iterator& it)
  : nt(it.nt)
  {}
  
  CMap_dart_handle_iterator& operator=(const CMap_dart_handle_iterator& it)
  {
    nt = const_cast<CMap_dart_handle_iterator&>(it).nt;
  }

private:
  Iterator nt;
  static CMap& get_default()
  {
    static CMap* m = new CMap();
    return *m;
  }
};

template <typename Dart_handle>
struct EdgeHandle : Dart_handle
{
  EdgeHandle() : Dart_handle(NULL){}
  EdgeHandle(const Dart_handle& h): Dart_handle(h)
  {}
};

template <class CMap>
struct CMap_Base_graph_traits
{

public :  
  struct CMap_graph_traversal_category : public virtual boost::bidirectional_graph_tag,
                                        public virtual boost::vertex_list_graph_tag,
                                        public virtual boost::edge_list_graph_tag
  {};

  // Expose types required by the boost::Graph concept.
  typedef typename CMap::Dart_handle vertex_descriptor;
  typedef EdgeHandle<typename CMap::Dart_handle> edge_descriptor;
  typedef boost::directed_tag directed_category;
  typedef boost::allow_parallel_edge_tag edge_parallel_category; 
  typedef CMap_graph_traversal_category traversal_category;

  // Expose types required by the boost::IncidenceGraph concept.
  typedef CMap_dart_handle_iterator<CMap, CMap_dart_iterator_of_cell<CMap, 0> > out_edge_iterator;
  typedef typename CMap::size_type degree_size_type;

  // Expose types required by the boost::BidirectionalGraph concept.
  typedef CMap_dart_handle_iterator<CMap, CMap_dart_iterator_of_second_vertex<CMap> > in_edge_iterator;
  typedef typename CMap::size_type edges_size_type;

  // Expose types required by the boost::EdgeListGraph concept.
  typedef CMap_dart_handle_iterator<CMap, typename CMap::Dart_range::iterator> edge_iterator;

  // Expose types required by the boost::VertexListGraph concept.
  typedef typename CMap::size_type vertices_size_type;
  typedef CMap_dart_handle_iterator<CMap, typename CMap::template One_dart_per_cell_range<0>::iterator> vertex_iterator;
};

} //namespace CGAL

namespace boost{

// Specialization of graph_traits for Combinatorial map.
CGAL_CMAP_TEMPLATE_ARGS
struct graph_traits<CGAL_CMAP_TYPE >
: CGAL::CMap_Base_graph_traits<typename CGAL_CMAP_TYPE::Base >
{};

CGAL_CMAP_TEMPLATE_ARGS
struct graph_traits<CGAL_CMAP_TYPE const>
: CGAL::CMap_Base_graph_traits<typename CGAL_CMAP_TYPE::Base >
{};

// Specialization of graph_traits for Combinatorial map base.
CGAL_CMAP_BASE_TEMPLATE_ARGS
struct graph_traits<CGAL_CMAP_BASE_TYPE >
: CGAL::CMap_Base_graph_traits<CGAL_CMAP_BASE_TYPE >
{};

CGAL_CMAP_BASE_TEMPLATE_ARGS
struct graph_traits<CGAL_CMAP_BASE_TYPE const>
: CGAL::CMap_Base_graph_traits<CGAL_CMAP_BASE_TYPE >
{};

// Specialization of graph_traits for Linear Cell Complex.
CGAL_LCC_TEMPLATE_ARGS
struct graph_traits<CGAL_LCC_TYPE >
: CGAL::CMap_Base_graph_traits<typename CGAL_LCC_TYPE::Base >
{};

CGAL_LCC_TEMPLATE_ARGS
struct graph_traits<CGAL_LCC_TYPE const>
: CGAL::CMap_Base_graph_traits<typename CGAL_LCC_TYPE::Base >
{};

// Expression required by the boost::IncidenceGraph concept.

CGAL_CMAP_BASE_TEMPLATE_ARGS
typename boost::graph_traits<CGAL_CMAP_BASE_TYPE >::vertex_descriptor 
source(typename boost::graph_traits<CGAL_CMAP_BASE_TYPE >::edge_descriptor e, const CGAL_CMAP_BASE_TYPE&)
{
  return e;
}

CGAL_CMAP_BASE_TEMPLATE_ARGS
typename boost::graph_traits<CGAL_CMAP_BASE_TYPE >::vertex_descriptor 
target(typename boost::graph_traits<CGAL_CMAP_BASE_TYPE >::edge_descriptor e, const CGAL_CMAP_BASE_TYPE&)
{
  return e->opposite();
}

CGAL_CMAP_BASE_TEMPLATE_ARGS
std::pair<typename boost::graph_traits<CGAL_CMAP_BASE_TYPE >::out_edge_iterator, 
          typename boost::graph_traits<CGAL_CMAP_BASE_TYPE >::out_edge_iterator>
out_edges(typename boost::graph_traits<CGAL_CMAP_BASE_TYPE >::vertex_descriptor u, const CGAL_CMAP_BASE_TYPE& cm)
{
  typedef typename boost::graph_traits<CGAL_CMAP_BASE_TYPE >::out_edge_iterator iter_type;

  CGAL_CMAP_BASE_TYPE& cmap = const_cast<CGAL_CMAP_BASE_TYPE&>(cm);

  return std::make_pair(
            cmap.template darts_of_cell<0>(u).begin(), 
            cmap.template darts_of_cell<0>(u).end());
}

CGAL_CMAP_BASE_TEMPLATE_ARGS
typename boost::graph_traits<CGAL_CMAP_BASE_TYPE >::degree_size_type
out_degree(typename boost::graph_traits<CGAL_CMAP_BASE_TYPE >::vertex_descriptor u, const CGAL_CMAP_BASE_TYPE& cm)
{
  typedef typename boost::graph_traits<CGAL_CMAP_BASE_TYPE >::out_edge_iterator iter_type;
  std::pair<iter_type, iter_type> iter = out_edges(u, cm);

  typename boost::graph_traits<CGAL_CMAP_BASE_TYPE >::degree_size_type degree=0;
  for(;iter.first != iter.second; ++(iter.first))
    ++degree;
  return degree;
}

// Expression required by the boost::BidirectionalGraph concept.

CGAL_CMAP_BASE_TEMPLATE_ARGS
std::pair<typename boost::graph_traits<CGAL_CMAP_BASE_TYPE >::in_edge_iterator, typename boost::graph_traits<CGAL_CMAP_BASE_TYPE >::in_edge_iterator>
in_edges(typename boost::graph_traits<CGAL_CMAP_BASE_TYPE >::vertex_descriptor v, const CGAL_CMAP_BASE_TYPE& cm)
{
  typedef typename boost::graph_traits<CGAL_CMAP_BASE_TYPE >::in_edge_iterator iter_type;

  CGAL_CMAP_BASE_TYPE& cmap = const_cast<CGAL_CMAP_BASE_TYPE&>(cm);

  return std::make_pair<iter_type, iter_type>
    (cmap.darts_of_second_vertex(v).begin(), cmap.darts_of_second_vertex(v).end());
}

CGAL_CMAP_BASE_TEMPLATE_ARGS
typename boost::graph_traits<CGAL_CMAP_BASE_TYPE >::degree_size_type
in_degree(typename boost::graph_traits<CGAL_CMAP_BASE_TYPE >::vertex_descriptor v, const CGAL_CMAP_BASE_TYPE& cm)
{
  typedef typename boost::graph_traits<CGAL_CMAP_BASE_TYPE >::in_edge_iterator iter_type;
  std::pair<iter_type, iter_type> iter = in_edges(v, cm);

  typename boost::graph_traits<CGAL_CMAP_BASE_TYPE >::degree_size_type degree=0;
  for(;iter.first != iter.second; ++(iter.first))
    ++degree;
  return degree;
}

// We suppose there are no loops.
CGAL_CMAP_BASE_TEMPLATE_ARGS
typename boost::graph_traits<CGAL_CMAP_BASE_TYPE >::degree_size_type
degree(typename boost::graph_traits<CGAL_CMAP_BASE_TYPE >::vertex_descriptor v, const CGAL_CMAP_BASE_TYPE& cm)
{
  return in_degree(v, cm) + out_degree(v, cm);
}

// Expression required by the boost::VertexListGraph concept.

CGAL_CMAP_BASE_TEMPLATE_ARGS
std::pair<typename boost::graph_traits<CGAL_CMAP_BASE_TYPE >::vertex_iterator, typename boost::graph_traits<CGAL_CMAP_BASE_TYPE >::vertex_iterator>
vertices(const CGAL_CMAP_BASE_TYPE& cm)
{
  typedef typename boost::graph_traits<CGAL_CMAP_BASE_TYPE >::vertex_iterator iter_type;

  CGAL_CMAP_BASE_TYPE& cmap = const_cast<CGAL_CMAP_BASE_TYPE&>(cm);

  return std::make_pair<iter_type, iter_type>
    (iter_type(cmap.template one_dart_per_cell<0>().begin()),
     iter_type(cmap.template one_dart_per_cell<0>().end()));
}

CGAL_CMAP_BASE_TEMPLATE_ARGS
typename boost::graph_traits<CGAL_CMAP_BASE_TYPE >::vertices_size_type
num_vertices(const CGAL_CMAP_BASE_TYPE& cm)
{
  CGAL_CMAP_BASE_TYPE& cmap = const_cast<CGAL_CMAP_BASE_TYPE&>(cm);
  return cmap.template one_dart_per_cell<0>().size();
}

// Expression required by the boost::EdgeListGraph concept.

CGAL_CMAP_BASE_TEMPLATE_ARGS
std::pair<typename boost::graph_traits<CGAL_CMAP_BASE_TYPE >::edge_iterator, typename boost::graph_traits<CGAL_CMAP_BASE_TYPE >::edge_iterator>
edges(const CGAL_CMAP_BASE_TYPE& cm)
{
  typedef typename boost::graph_traits<CGAL_CMAP_BASE_TYPE >::edge_iterator iter_type;

  CGAL_CMAP_BASE_TYPE& cmap = const_cast<CGAL_CMAP_BASE_TYPE&>(cm);

  return std::make_pair<iter_type, iter_type>
    (iter_type(cmap.darts().begin()),
     iter_type(cmap.darts().end()));
}

CGAL_CMAP_BASE_TEMPLATE_ARGS
typename boost::graph_traits<CGAL_CMAP_BASE_TYPE >::edges_size_type
num_edges(const CGAL_CMAP_BASE_TYPE& cm)
{
  typedef typename boost::graph_traits<CGAL_CMAP_BASE_TYPE >::edge_iterator iter_type;
  std::pair<iter_type, iter_type> iter = edges(cm);

  typename boost::graph_traits<CGAL_CMAP_BASE_TYPE >::edges_size_type degree=0;
  for(;iter.first != iter.second; ++(iter.first))
    ++degree;
  return degree;
}


}// namespace boost

#undef CGAL_CMAP_BASE_TEMPLATE_ARGS
#undef CGAL_CMAP_TEMPLATE_ARGS
#undef CGAL_CMAP_TYPE
#undef CGAL_CMAP_BASE_TYPE
#undef CGAL_LCC_TEMPLATE_ARGS
#undef CGAL_LCC_TYPE

#endif // CGAL_BOOST_GRAPH_GRAPH_TRAITS_CMAP_H
