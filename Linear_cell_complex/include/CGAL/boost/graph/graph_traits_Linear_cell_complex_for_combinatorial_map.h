// Copyright (c) 2017 CNRS and LIRIS' Establishments (France).
// All rights reserved.
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
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef CGAL_BOOST_GRAPH_TRAITS_LINEAR_CELL_COMPLEX_FOR_COMBINATORIAL_MAP_H
#define CGAL_BOOST_GRAPH_TRAITS_LINEAR_CELL_COMPLEX_FOR_COMBINATORIAL_MAP_H

#include <utility>
#include <iterator>

#include <boost/config.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/properties_Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/boost/graph/graph_traits_HalfedgeDS.h>

#include <CGAL/Linear_cell_complex_for_bgl_combinatorial_map_helper.h>
#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/boost/graph/io.h>
#include <CGAL/Dart_iterators.h>

#include <CGAL/boost/graph/helpers.h>
#include <CGAL/iterator.h>

#define CGAL_LCC_TEMPLATE_ARGS template<unsigned int d_, unsigned int ambient_dim, \
                                        class Traits_,                \
                                        class Items_,                 \
                                        class Alloc_,                 \
                                        template<unsigned int,class,class,class,class> \
                                        class CMap,                   \
                                        class Storage_>
#define CGAL_LCC_TEMPLATE_ARGS_NOTEND template<unsigned int d_, unsigned int ambient_dim, \
                                        class Traits_,                \
                                        class Items_,                 \
                                        class Alloc_,                 \
                                        template<unsigned int,class,class,class,class> \
                                        class CMap,                   \
                                        class Storage_,

#define CGAL_LCC_TYPE CGAL::Linear_cell_complex_for_combinatorial_map \
                <d_, ambient_dim, Traits_, Items_, Alloc_, CMap, Storage_>

namespace CGAL {

namespace internal
{

/// A struct to define edge based on halfedge.
template <typename Dart_handle>
struct EdgeHandle : Dart_handle
{
  EdgeHandle() : Dart_handle(NULL){}
  explicit EdgeHandle(Dart_handle h): Dart_handle(h)
  {}

  Dart_handle first_halfedge() const
  { return *this; }

  Dart_handle second_halfedge() const
  {
    assert(*this!=NULL);
    return (*this)->get_f(2);
  }
  
  bool operator==(const EdgeHandle& h) const
  {
    return first_halfedge()==h.first_halfedge() ||
        (h.first_halfedge()!=NULL &&
        first_halfedge()==h.second_halfedge());
  }

  bool operator!=(const EdgeHandle& other) const
  { return !(*this==other); }

  friend bool operator<(const EdgeHandle& a, const EdgeHandle& b)
  { return a.first_halfedge()<b.first_halfedge(); }

  friend bool operator>(const EdgeHandle& a, const EdgeHandle& b)
  { return b<a; }

  friend bool operator<=(const EdgeHandle& a, const EdgeHandle& b)
  { return !(a>b); }

  friend bool operator>=(const EdgeHandle& a, const EdgeHandle& b)
  { return !(a<b); }

  std::size_t id() const
  { return first_halfedge()->id()/2; }

  friend std::size_t hash_value(const EdgeHandle& i)
  {
    if (i.first_halfedge()==NULL) return 0;
    return hash_value(i.first_halfedge()<i.second_halfedge()?
                      i.first_halfedge():i.second_halfedge());
  }
};

// make edge_descriptor hashable by default in Unique_hash_map
namespace handle{
  template<typename Dart_handle>
  struct Hash_functor< EdgeHandle<Dart_handle> >
  {
    std::size_t
    operator()(const EdgeHandle<Dart_handle>& edge)
    { return hash_value(edge); }
  };
} //end of namespace handle

template <class CMap, typename Dart_Iterator>
class CMap_dart_handle_edge_iterator 
{
public:
  CMap_dart_handle_edge_iterator(){}

  typedef Dart_Iterator Iterator;    

  typedef typename CMap::Dart_handle Dart_handle;

  typedef CMap_dart_handle_edge_iterator<CMap, Dart_Iterator> Self; 

  typedef typename std::iterator_traits<Iterator>::iterator_category iterator_category;
  typedef typename std::iterator_traits<Iterator>::difference_type   difference_type;
  typedef EdgeHandle<Dart_handle>                                    value_type;
  typedef value_type                                                 reference;
  typedef value_type                                                 pointer;

public:

// OPERATIONS Forward Category
// ---------------------------

  bool operator==( const Self& i) const { return ( nt == i.nt); }
  bool operator!=( const Self& i) const { return !(nt == i.nt );}
  value_type operator*() const { return value_type(nt); }
  value_type operator->() { return value_type(nt); }

  Self& operator++()
  {
    typedef typename Dart_handle::CC CC;
    ++nt;
    // We need to test if we are at the end of the compact container.
    // (case where we were previously on the last element)
    if (!CC::is_begin_or_end(&*value_type(nt)))
      ++nt; // halfedges are always created by pair (two consecutive elements)
    return *this;
  }

  Self  operator++(int)
  {
    Self tmp = *this;
    operator++();
    return tmp;
  }

  CMap_dart_handle_edge_iterator(const Iterator& iter) :
    nt(iter)
  {}

private:
  Iterator nt;
};

} // internal

template <class CMap>
struct CMap_Base_graph_traits
{

public :  
  struct CMap_graph_traversal_category : public virtual boost::bidirectional_graph_tag,
                                         public virtual boost::vertex_list_graph_tag,
                                         public virtual boost::edge_list_graph_tag
  {};

  // Expose types required by the boost::Graph concept.
  typedef typename CMap::template Attribute_handle<0>::type vertex_descriptor;
  typedef internal::EdgeHandle<typename CMap::Dart_handle>  edge_descriptor;
  typedef typename CMap::template Attribute_handle<2>::type face_descriptor;
  typedef typename CMap::Dart_handle                        halfedge_descriptor;

  typedef boost::directed_tag directed_category;
  typedef boost::allow_parallel_edge_tag edge_parallel_category; 
  typedef CMap_graph_traversal_category traversal_category;

  typedef Prevent_deref<typename CMap::template Attribute_range<0>::type::iterator> vertex_iterator;
  typedef Prevent_deref<typename CMap::template Attribute_range<2>::type::iterator> face_iterator;
  typedef Prevent_deref<typename CMap::Dart_range::iterator> halfedge_iterator;
  
  typedef internal::CMap_dart_handle_edge_iterator<CMap, typename CMap::Dart_range::iterator> edge_iterator;
  
  typedef typename CMap::size_type degree_size_type;
  typedef typename CMap::size_type halfedges_size_type;
  typedef typename CMap::size_type vertices_size_type;
  typedef typename CMap::size_type edges_size_type;
  typedef typename CMap::size_type faces_size_type;
  
  typedef CGAL::In_edge_iterator<CMap>  in_edge_iterator;
  typedef CGAL::Out_edge_iterator<CMap> out_edge_iterator;

  // nulls
  static vertex_descriptor   null_vertex()   { return NULL; }
  static face_descriptor     null_face()     { return NULL; }
  static halfedge_descriptor null_halfedge() { return NULL; }
};

} //namespace CGAL

namespace boost
{
  // Specialization of graph_traits for Linear Cell Complex.
  CGAL_LCC_TEMPLATE_ARGS
  struct graph_traits<CGAL_LCC_TYPE >
    : public CGAL::CMap_Base_graph_traits<typename CGAL_LCC_TYPE >
  {
    typedef typename CGAL_LCC_TYPE::Point vertex_property_type;
  };
  
  CGAL_LCC_TEMPLATE_ARGS
  struct graph_traits<CGAL_LCC_TYPE const>
    : public CGAL::CMap_Base_graph_traits<typename CGAL_LCC_TYPE >
  {
    typedef typename CGAL_LCC_TYPE::Point vertex_property_type;
  };

}// namespace boost

namespace CGAL 
{

//  
// 1) boost::Graph
//
  
CGAL_LCC_TEMPLATE_ARGS
typename boost::graph_traits<CGAL_LCC_TYPE >::vertices_size_type
num_vertices(const CGAL_LCC_TYPE& lcc)
{ return lcc.template attributes<0>().size(); }

CGAL_LCC_TEMPLATE_ARGS
typename boost::graph_traits<CGAL_LCC_TYPE>::edges_size_type
num_edges(const CGAL_LCC_TYPE& lcc)
{ return lcc.number_of_darts()/2; }
  
// We suppose there are no loops.
CGAL_LCC_TEMPLATE_ARGS
typename boost::graph_traits<CGAL_LCC_TYPE >::degree_size_type
degree(typename boost::graph_traits<CGAL_LCC_TYPE >::vertex_descriptor v,
       const CGAL_LCC_TYPE& lcc)
{
  typename boost::graph_traits<CGAL_LCC_TYPE >::degree_size_type degree=0;
  for (typename CGAL_LCC_TYPE::template Dart_of_cell_range<0>::const_iterator
         it=lcc.template darts_of_cell<0>(v->dart()).begin(),
         itend=lcc.template darts_of_cell<0>(v->dart()).end();
       it!=itend; ++it)
  { ++degree; }
  return degree;
}

CGAL_LCC_TEMPLATE_ARGS
typename boost::graph_traits<CGAL_LCC_TYPE >::degree_size_type
out_degree(typename boost::graph_traits<CGAL_LCC_TYPE >::vertex_descriptor u,
           const CGAL_LCC_TYPE& lcc)
{ return degree(u, lcc); }

CGAL_LCC_TEMPLATE_ARGS
typename boost::graph_traits<CGAL_LCC_TYPE >::degree_size_type
in_degree(typename boost::graph_traits<CGAL_LCC_TYPE >::vertex_descriptor v,
          const CGAL_LCC_TYPE& lcc)
{ return degree(v, lcc); }

CGAL_LCC_TEMPLATE_ARGS
typename boost::graph_traits<CGAL_LCC_TYPE >::degree_size_type
degree(typename boost::graph_traits<CGAL_LCC_TYPE >::face_descriptor f,
       const CGAL_LCC_TYPE& lcc)

{
  typename boost::graph_traits<CGAL_LCC_TYPE >::degree_size_type degree=0;
  for (typename CGAL_LCC_TYPE::template Dart_of_cell_range<2>::const_iterator
         it=lcc.template darts_of_cell<2>(f->dart()).begin(),
         itend=lcc.template darts_of_cell<2>(f->dart()).end();
       it!=itend; ++it)
  { ++degree; }
  return degree;
}

CGAL_LCC_TEMPLATE_ARGS
typename boost::graph_traits<CGAL_LCC_TYPE >::vertex_descriptor 
source(typename boost::graph_traits<CGAL_LCC_TYPE >::halfedge_descriptor h,
       const CGAL_LCC_TYPE& lcc)
{ return const_cast<CGAL_LCC_TYPE&>(lcc).template attribute<0>(h); }

CGAL_LCC_TEMPLATE_ARGS
typename boost::graph_traits<CGAL_LCC_TYPE >::vertex_descriptor 
target(typename boost::graph_traits<CGAL_LCC_TYPE >::halfedge_descriptor h,
       const CGAL_LCC_TYPE& lcc)
{
  return const_cast<CGAL_LCC_TYPE&>(lcc).template attribute<0>
    (const_cast<CGAL_LCC_TYPE&>(lcc).template beta<2>(h));
}

CGAL_LCC_TEMPLATE_ARGS
typename boost::graph_traits<CGAL_LCC_TYPE >::vertex_descriptor 
source(typename boost::graph_traits<CGAL_LCC_TYPE >::edge_descriptor e,
       const CGAL_LCC_TYPE& lcc)
{ return source(e.first_halfedge(), lcc); }

CGAL_LCC_TEMPLATE_ARGS
typename boost::graph_traits<CGAL_LCC_TYPE >::vertex_descriptor 
target(typename boost::graph_traits<CGAL_LCC_TYPE >::edge_descriptor e,
       const CGAL_LCC_TYPE& lcc)
{ return target(e.first_halfedge(), lcc); }

CGAL_LCC_TEMPLATE_ARGS
std::pair<
  typename boost::graph_traits<CGAL_LCC_TYPE>::halfedge_descriptor,
  bool>
halfedge(typename boost::graph_traits<CGAL_LCC_TYPE>::vertex_descriptor u,
         typename boost::graph_traits<CGAL_LCC_TYPE>::vertex_descriptor v,
         const CGAL_LCC_TYPE& lcc)
{
  for (typename CGAL_LCC_TYPE::template Dart_of_cell_range<0>::iterator
         it=const_cast<CGAL_LCC_TYPE&>(lcc).template
       darts_of_cell<0>(u->dart()).begin(),
         itend=const_cast<CGAL_LCC_TYPE&>(lcc).template
       darts_of_cell<0>(u->dart()).end();
       it!=itend; ++it)
  {
    if (lcc.template attribute<0>(lcc.template beta<2>(it))==v)
    { return std::make_pair(it, true); }
  }

  return std::make_pair(lcc.null_handle, false);
}

CGAL_LCC_TEMPLATE_ARGS
std::pair<typename boost::graph_traits<CGAL_LCC_TYPE>::edge_descriptor, bool>
edge(typename boost::graph_traits<CGAL_LCC_TYPE>::vertex_descriptor u, 
     typename boost::graph_traits<CGAL_LCC_TYPE>::vertex_descriptor v, 
     const CGAL_LCC_TYPE& lcc)
{
  std::pair<typename boost::graph_traits<CGAL_LCC_TYPE>::halfedge_descriptor,
            bool> res=halfedge(u,v,lcc);
  return std::make_pair(internal::EdgeHandle<typename CGAL_LCC_TYPE::Dart_handle>(res.first),
                        res.second);
}

CGAL_LCC_TEMPLATE_ARGS
CGAL::Iterator_range<typename boost::graph_traits<CGAL_LCC_TYPE >::vertex_iterator>
vertices(const CGAL_LCC_TYPE& lcc)
{
  typedef typename boost::graph_traits<CGAL_LCC_TYPE >::vertex_iterator Iter;
  CGAL_LCC_TYPE& lccnc = const_cast<CGAL_LCC_TYPE&>(lcc);
  return CGAL::make_range(Iter(lccnc.template attributes<0>().begin()),
                          Iter(lccnc.template attributes<0>().end()));
}

CGAL_LCC_TEMPLATE_ARGS
CGAL::Iterator_range<typename boost::graph_traits<CGAL_LCC_TYPE>::edge_iterator>
edges(const CGAL_LCC_TYPE& lcc)
{
  typedef typename boost::graph_traits<CGAL_LCC_TYPE >::edge_iterator Iter;
  CGAL_LCC_TYPE& lccnc = const_cast<CGAL_LCC_TYPE&>(lcc);
  return CGAL::make_range(Iter(lccnc.darts().begin()),
                          Iter(lccnc.darts().end()));
}

CGAL_LCC_TEMPLATE_ARGS
CGAL::Iterator_range<typename boost::graph_traits<CGAL_LCC_TYPE>::in_edge_iterator>
in_edges(typename boost::graph_traits<CGAL_LCC_TYPE >::vertex_descriptor v,
         const CGAL_LCC_TYPE& lcc)
{
  typedef typename boost::graph_traits<CGAL_LCC_TYPE >::in_edge_iterator Iter;
  return make_range(Iter(halfedge(v, lcc), lcc), Iter(halfedge(v, lcc), lcc, 1));
}

CGAL_LCC_TEMPLATE_ARGS
CGAL::Iterator_range<typename boost::graph_traits<CGAL_LCC_TYPE>::out_edge_iterator>
out_edges(typename boost::graph_traits<CGAL_LCC_TYPE>::vertex_descriptor v,
          const CGAL_LCC_TYPE& lcc)
{
  typedef typename boost::graph_traits<CGAL_LCC_TYPE>::out_edge_iterator Iter;
  return make_range(Iter(halfedge(v, lcc), lcc), Iter(halfedge(v, lcc), lcc, 1));
}

//
// 2) HalfedgeGraph
//

CGAL_LCC_TEMPLATE_ARGS
typename boost::graph_traits<CGAL_LCC_TYPE>::edge_descriptor
edge(typename boost::graph_traits<CGAL_LCC_TYPE>::halfedge_descriptor h,
     const CGAL_LCC_TYPE&/* lcc*/)
{ return internal::EdgeHandle
      <typename boost::graph_traits<CGAL_LCC_TYPE>::halfedge_descriptor>(h); }

CGAL_LCC_TEMPLATE_ARGS
typename boost::graph_traits<CGAL_LCC_TYPE>::halfedge_descriptor
halfedge(typename boost::graph_traits<CGAL_LCC_TYPE>::edge_descriptor e,
         const CGAL_LCC_TYPE&)
{ return e.first_halfedge(); }
  
CGAL_LCC_TEMPLATE_ARGS
typename boost::graph_traits<CGAL_LCC_TYPE>::halfedge_descriptor
halfedge(typename boost::graph_traits<CGAL_LCC_TYPE>::vertex_descriptor v,
         const CGAL_LCC_TYPE& lcc)
{
  if (v->dart()==NULL) return NULL;
  return const_cast<CGAL_LCC_TYPE&>(lcc).template beta<2>(v->dart());
}

CGAL_LCC_TEMPLATE_ARGS
typename boost::graph_traits<CGAL_LCC_TYPE>::halfedge_descriptor
opposite(typename boost::graph_traits<CGAL_LCC_TYPE>::halfedge_descriptor dh,
         const CGAL_LCC_TYPE& lcc)
{ return const_cast<CGAL_LCC_TYPE&>(lcc).template beta<2>(dh); }

CGAL_LCC_TEMPLATE_ARGS
typename boost::graph_traits<CGAL_LCC_TYPE>::halfedge_descriptor
next(typename boost::graph_traits<CGAL_LCC_TYPE>::halfedge_descriptor dh,
         const CGAL_LCC_TYPE& lcc)
{ return const_cast<CGAL_LCC_TYPE&>(lcc).template beta<1>(dh); }

CGAL_LCC_TEMPLATE_ARGS
typename boost::graph_traits<CGAL_LCC_TYPE>::halfedge_descriptor
prev(typename boost::graph_traits<CGAL_LCC_TYPE>::halfedge_descriptor dh,
         const CGAL_LCC_TYPE& lcc)
{ return const_cast<CGAL_LCC_TYPE&>(lcc).template beta<0>(dh); }

//
// 3) HalfedgeListGraph
//

CGAL_LCC_TEMPLATE_ARGS
Iterator_range<typename boost::graph_traits<CGAL_LCC_TYPE>::halfedge_iterator>
halfedges(const CGAL_LCC_TYPE& lcc)
{
  typedef typename boost::graph_traits<CGAL_LCC_TYPE>::halfedge_iterator Iter;
  CGAL_LCC_TYPE& lccnc = const_cast<CGAL_LCC_TYPE&>(lcc);
  return CGAL::make_range(Iter(lccnc.darts().begin()),
                          Iter(lccnc.darts().end()));
}

CGAL_LCC_TEMPLATE_ARGS
typename boost::graph_traits<CGAL_LCC_TYPE>::halfedges_size_type
num_halfedges(const CGAL_LCC_TYPE& lcc)
{ return lcc.number_of_darts(); }

//
// 4) MutableHalfedgeGraph
// 

CGAL_LCC_TEMPLATE_ARGS
typename boost::graph_traits<CGAL_LCC_TYPE>::vertex_descriptor
add_vertex(CGAL_LCC_TYPE& lcc)
{ return lcc.template create_attribute<0>(); }

CGAL_LCC_TEMPLATE_ARGS
typename boost::graph_traits<CGAL_LCC_TYPE>::vertex_descriptor
add_vertex(const typename boost::graph_traits<CGAL_LCC_TYPE>::vertex_property_type& p,
           CGAL_LCC_TYPE& lcc)
{ return lcc.template create_attribute<0>(p); }

CGAL_LCC_TEMPLATE_ARGS
void
remove_vertex(typename boost::graph_traits<CGAL_LCC_TYPE>::vertex_descriptor v,
              CGAL_LCC_TYPE& lcc)
{ lcc.template erase_attribute<0>(v); }

CGAL_LCC_TEMPLATE_ARGS
void set_halfedge(typename boost::graph_traits<CGAL_LCC_TYPE>::vertex_descriptor v,
                  typename boost::graph_traits<CGAL_LCC_TYPE>::halfedge_descriptor h,
                  CGAL_LCC_TYPE& lcc)
{
  if (h!=NULL)
    lcc.template set_dart_of_attribute<0>(v, lcc.template beta<2>(h));
  else
    lcc.template set_dart_of_attribute<0>(v, h);
}

CGAL_LCC_TEMPLATE_ARGS
typename boost::graph_traits<CGAL_LCC_TYPE>::edge_descriptor
add_edge(CGAL_LCC_TYPE& lcc)
{
  typename CGAL_LCC_TYPE::Dart_handle actu = lcc.create_dart();
  lcc.template link_beta<2>(actu, lcc.create_dart());
  return typename boost::graph_traits<CGAL_LCC_TYPE>::edge_descriptor(actu);
}
  
CGAL_LCC_TEMPLATE_ARGS
void
remove_edge(typename boost::graph_traits<CGAL_LCC_TYPE>::edge_descriptor e, 
            CGAL_LCC_TYPE& lcc)
{
  assert ( !lcc.template is_free<2>(e.first_halfedge()) );
  lcc.restricted_erase_dart(lcc.template beta<2>(e.first_halfedge()));
  lcc.restricted_erase_dart(e);
}

CGAL_LCC_TEMPLATE_ARGS
void set_target(typename boost::graph_traits<CGAL_LCC_TYPE>::halfedge_descriptor h1,
                typename boost::graph_traits<CGAL_LCC_TYPE>::vertex_descriptor v,
                CGAL_LCC_TYPE& lcc)
{ lcc.template restricted_set_dart_attribute<0>(lcc.template beta<2>(h1), v); }

CGAL_LCC_TEMPLATE_ARGS
void set_next(typename boost::graph_traits<CGAL_LCC_TYPE>::halfedge_descriptor h1,
              typename boost::graph_traits<CGAL_LCC_TYPE>::halfedge_descriptor h2,
              CGAL_LCC_TYPE& lcc)
{ lcc.basic_link_beta_1(h1, h2); }

//
// 5) FaceGraph
//

CGAL_LCC_TEMPLATE_ARGS
typename boost::graph_traits<CGAL_LCC_TYPE>::face_descriptor
face(typename boost::graph_traits<CGAL_LCC_TYPE>::halfedge_descriptor h,
     const CGAL_LCC_TYPE& lcc)
{ return const_cast<CGAL_LCC_TYPE&>(lcc).template attribute<2>(h); }

CGAL_LCC_TEMPLATE_ARGS
typename boost::graph_traits<CGAL_LCC_TYPE>::halfedge_descriptor
halfedge(typename boost::graph_traits<CGAL_LCC_TYPE>::face_descriptor f,
         const CGAL_LCC_TYPE&) 
{ return f->dart(); }
  
CGAL_LCC_TEMPLATE_ARGS
CGAL::Iterator_range<typename boost::graph_traits<CGAL_LCC_TYPE>::face_iterator>
faces(const CGAL_LCC_TYPE& lcc)
{
  typedef typename boost::graph_traits<CGAL_LCC_TYPE >::face_iterator iter_type;
  CGAL_LCC_TYPE& lccnc = const_cast<CGAL_LCC_TYPE&>(lcc);
  return CGAL::make_range(iter_type(lccnc.template attributes<2>().begin()),
                         iter_type(lccnc.template attributes<2>().end()));
}

CGAL_LCC_TEMPLATE_ARGS
typename boost::graph_traits<CGAL_LCC_TYPE>::faces_size_type
num_faces(const CGAL_LCC_TYPE& lcc)
{ return lcc.template attributes<2>().size(); }
  
CGAL_LCC_TEMPLATE_ARGS
void reserve(CGAL_LCC_TYPE& g,
             typename boost::graph_traits<CGAL_LCC_TYPE>::vertices_size_type nv,
             typename boost::graph_traits<CGAL_LCC_TYPE>::edges_size_type ne,
             typename boost::graph_traits<CGAL_LCC_TYPE>::faces_size_type nf)
{
  g.template attributes<0>().reserve(nv);
  g.darts().reserve(2*ne);
  g.template attributes<2>().reserve(nf);
}

CGAL_LCC_TEMPLATE_ARGS
bool is_valid(const CGAL_LCC_TYPE& lcc, bool verbose=false)
{
  if (verbose)
    lcc.display_characteristics(std::cout)<<std::endl;
  return lcc.is_valid();
}

//
// 6) MutableFaceGraph 
//

CGAL_LCC_TEMPLATE_ARGS
typename boost::graph_traits<CGAL_LCC_TYPE>::face_descriptor
add_face(CGAL_LCC_TYPE& lcc)
{ return lcc.template create_attribute<2>(); }

CGAL_LCC_TEMPLATE_ARGS_NOTEND typename InputIterator>
typename boost::graph_traits<CGAL_LCC_TYPE>::face_descriptor
add_face(InputIterator begin, InputIterator end, CGAL_LCC_TYPE& lcc)
{
  while(begin!=end)
  {
    lcc.template create_attribute<2>();
    ++begin;
  }
}

CGAL_LCC_TEMPLATE_ARGS
void remove_face(typename boost::graph_traits<CGAL_LCC_TYPE>::face_descriptor f,
                 CGAL_LCC_TYPE& lcc)
{ lcc.template erase_attribute<2>(f); }

CGAL_LCC_TEMPLATE_ARGS
void set_face(typename boost::graph_traits<CGAL_LCC_TYPE>::halfedge_descriptor h,              
              typename boost::graph_traits<CGAL_LCC_TYPE>::face_descriptor f,
              CGAL_LCC_TYPE& lcc)
{ lcc.template restricted_set_dart_attribute<2>(h, f); }

CGAL_LCC_TEMPLATE_ARGS
void set_halfedge(typename boost::graph_traits<CGAL_LCC_TYPE>::face_descriptor f,
                  typename boost::graph_traits<CGAL_LCC_TYPE>::halfedge_descriptor h,
                  const CGAL_LCC_TYPE&)
{ f->set_dart(h); }

} // namespace CGAL

#undef CGAL_LCC_TEMPLATE_ARGS
#undef CGAL_LCC_TYPE

#endif // CGAL_BOOST_GRAPH_TRAITS_LINEAR_CELL_COMPLEX_FOR_COMBINATORIAL_MAP_H
