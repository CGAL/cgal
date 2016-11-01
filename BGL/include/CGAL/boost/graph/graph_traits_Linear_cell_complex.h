// Copyright (c) 2016  GeometryFactory (France).  All rights reserved.
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
//
// Author(s)     : Guillaume Damiand

#ifndef CGAL_BOOST_GRAPH_TRAITS_LINEAR_CELL_COMPLEX_H
#define CGAL_BOOST_GRAPH_TRAITS_LINEAR_CELL_COMPLEX_H

#include <utility>
#include <iterator>

#include <boost/config.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/properties_Linear_cell_complex.h>
#include <CGAL/boost/graph/graph_traits_HalfedgeDS.h>

#include <CGAL/Linear_cell_complex.h>
#include <CGAL/Dart_iterators.h>
#include <CGAL/boost/graph/helpers.h>

#define CGAL_LCC_TEMPLATE_ARGS template < unsigned int d_, unsigned int ambient_dim, \
                                          class Traits_,                \
                                          class Items_,                 \
                                          class Alloc_,                 \
                                          template<unsigned int, class,class,class,class> \
                                          class CMap,                   \
                                          class Storage_>
#define CGAL_LCC_TEMPLATE_ARGS_NOTEND template < unsigned int d_, unsigned int ambient_dim, \
                                          class Traits_,                \
                                          class Items_,                 \
                                          class Alloc_,                 \
                                          template<unsigned int, class,class,class,class> \
                                          class CMap,                   \
                                          class Storage_,

#define CGAL_LCC_TYPE CGAL::Linear_cell_complex<d_, ambient_dim, Traits_, Items_, Alloc_, CMap, Storage_> 

namespace CGAL {

template<typename CMap>
typename CMap::Dart_handle opposite(typename CMap::Dart_handle dh, const CMap& cmap)
{ return const_cast<CMap&>(cmap).template beta<2>(dh); }

template<typename CMap>
typename CMap::Dart_handle prev(typename CMap::Dart_handle dh, const CMap& cmap)
{ return const_cast<CMap&>(cmap).template beta<0>(dh); }
  
template<typename CMap>
typename CMap::Dart_handle next(typename CMap::Dart_handle dh, const CMap& cmap)
{ return const_cast<CMap&>(cmap).template beta<1>(dh); }

template <typename Dart_handle>
struct EdgeHandle : Dart_handle
{
  EdgeHandle() : Dart_handle(NULL){}
  EdgeHandle(Dart_handle& h): Dart_handle(h)
  {}
  EdgeHandle(const Dart_handle& h): Dart_handle(h)
  {}

  Dart_handle first_halfedge()
  { return *this; }

  Dart_handle second_halfedge()
  { return this->beta(2); }
  
  bool operator==(const EdgeHandle& h)
  { return (*this)==h || h->beta(2)==*this; }
};

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

  Self& operator++()
  {
    ++nt; ++nt; // halfedges are always created by pair (two consecutive elements)
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
  typedef EdgeHandle<typename CMap::Dart_handle>            edge_descriptor;
  typedef typename CMap::template Attribute_handle<2>::type face_descriptor;
  typedef typename CMap::Dart_handle                        halfedge_descriptor;

  typedef boost::directed_tag directed_category;
  typedef boost::allow_parallel_edge_tag edge_parallel_category; 
  typedef CMap_graph_traversal_category traversal_category;

  typedef internal::Prevent_deref<typename CMap::template Attribute_range<0>::type::iterator> vertex_iterator;
  typedef internal::Prevent_deref<typename CMap::template Attribute_range<2>::type::iterator> face_iterator;
  typedef internal::Prevent_deref<typename CMap::Dart_range::iterator> halfedge_iterator;
  
  typedef CMap_dart_handle_edge_iterator<CMap, typename CMap::Dart_range::iterator> edge_iterator;
  
  typedef typename CMap::size_type degree_size_type;
  typedef typename CMap::size_type halfedges_size_type;
  typedef typename CMap::size_type vertices_size_type;
  typedef typename CMap::size_type edges_size_type;
  typedef typename CMap::size_type faces_size_type;
  
  typedef CGAL::In_edge_iterator<CMap> in_edge_iterator;
  typedef CGAL::Out_edge_iterator<CMap> out_edge_iterator;

  // nulls
  static vertex_descriptor   null_vertex()   { return NULL; } // vertex_descriptor(); }
  static face_descriptor     null_face()     { return NULL; } // face_descriptor(); }
  static halfedge_descriptor null_halfedge() { return NULL; } // halfedge_descriptor(); }
};

} //namespace CGAL

namespace boost{
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
// Expression required by the boost::IncidenceGraph concept.

CGAL_LCC_TEMPLATE_ARGS
typename boost::graph_traits<CGAL_LCC_TYPE >::vertex_descriptor 
source(typename boost::graph_traits<CGAL_LCC_TYPE >::edge_descriptor e,
       const CGAL_LCC_TYPE& amap)
{
  return const_cast<CGAL_LCC_TYPE&>(amap).template attribute<0>(e);
  /*  return const_cast<CGAL_LCC_TYPE&>(amap).template beta<2>(e)->
      template attribute<0>(); */
}

CGAL_LCC_TEMPLATE_ARGS
typename boost::graph_traits<CGAL_LCC_TYPE >::vertex_descriptor 
target(typename boost::graph_traits<CGAL_LCC_TYPE >::edge_descriptor e,
       const CGAL_LCC_TYPE& amap)
{
  return const_cast<CGAL_LCC_TYPE&>(amap).template attribute<0>
    (const_cast<CGAL_LCC_TYPE&>(amap).template beta<2>(e));
  //  return const_cast<CGAL_LCC_TYPE&>(amap).template attribute<0>(e);
}

CGAL_LCC_TEMPLATE_ARGS
std::pair<typename boost::graph_traits<CGAL_LCC_TYPE >::out_edge_iterator, 
          typename boost::graph_traits<CGAL_LCC_TYPE >::out_edge_iterator>
out_edges(typename boost::graph_traits<CGAL_LCC_TYPE >::vertex_descriptor u,
          const CGAL_LCC_TYPE& cm)
{
  typedef typename boost::graph_traits<CGAL_LCC_TYPE >::out_edge_iterator
      iter_type;

  //  CGAL_LCC_TYPE& cmap = const_cast<CGAL_LCC_TYPE&>(cm);

  return std::make_pair<iter_type, iter_type>(iter_type(u->dart(), cm),
                                              iter_type(u->dart(), cm, 1));
}

CGAL_LCC_TEMPLATE_ARGS
typename boost::graph_traits<CGAL_LCC_TYPE >::degree_size_type
out_degree(typename boost::graph_traits<CGAL_LCC_TYPE >::vertex_descriptor u,
           const CGAL_LCC_TYPE& cm)
{ return degree(u, cm); }

// Expression required by the boost::BidirectionalGraph concept.

CGAL_LCC_TEMPLATE_ARGS
std::pair<typename boost::graph_traits<CGAL_LCC_TYPE >::in_edge_iterator,
          typename boost::graph_traits<CGAL_LCC_TYPE >::in_edge_iterator>
in_edges(typename boost::graph_traits<CGAL_LCC_TYPE >::vertex_descriptor v,
         const CGAL_LCC_TYPE& cm)
{
  typedef typename boost::graph_traits<CGAL_LCC_TYPE >::in_edge_iterator
      iter_type;

  // CGAL_LCC_TYPE& cmap = const_cast<CGAL_LCC_TYPE&>(cm);

  return std::make_pair<iter_type, iter_type>(iter_type(v->dart(), cm),
                                              iter_type(v->dart(), cm, 1));
}

CGAL_LCC_TEMPLATE_ARGS
typename boost::graph_traits<CGAL_LCC_TYPE >::degree_size_type
in_degree(typename boost::graph_traits<CGAL_LCC_TYPE >::vertex_descriptor v,
          const CGAL_LCC_TYPE& cm)
{ return degree(v, cm); }

// We suppose there are no loops.
CGAL_LCC_TEMPLATE_ARGS
typename boost::graph_traits<CGAL_LCC_TYPE >::degree_size_type
degree(typename boost::graph_traits<CGAL_LCC_TYPE >::vertex_descriptor v,
       const CGAL_LCC_TYPE& cm)
{
  typename boost::graph_traits<CGAL_LCC_TYPE >::degree_size_type degree=0;

  for (typename CGAL_LCC_TYPE::template Dart_of_cell_range<0>::const_iterator
         it=cm.template darts_of_cell<0>(v->dart()).begin(),
         itend=cm.template darts_of_cell<0>(v->dart()).end();
         //         it=cm.template darts_of_cell<0>(const_cast<CGAL_LCC_TYPE&>(cm).template beta<2>(v->dart())).begin(),
         //         itend=cm.template darts_of_cell<0>(const_cast<CGAL_LCC_TYPE&>(cm).template beta<2>(v->dart())).end();
       it!=itend; ++it)
  {
    ++degree;
  }
  return degree;
}

// Expression required by the boost::VertexListGraph concept.

CGAL_LCC_TEMPLATE_ARGS
CGAL::Iterator_range<typename boost::graph_traits<CGAL_LCC_TYPE >::vertex_iterator>
vertices(const CGAL_LCC_TYPE& cm)
{
  typedef typename boost::graph_traits<CGAL_LCC_TYPE >::vertex_iterator
      iter_type;

  CGAL_LCC_TYPE& cmap = const_cast<CGAL_LCC_TYPE&>(cm);

  return CGAL::make_range(iter_type(cmap.template attributes<0>().begin()),
                          iter_type(cmap.template attributes<0>().end()));
}

CGAL_LCC_TEMPLATE_ARGS
typename boost::graph_traits<CGAL_LCC_TYPE >::vertices_size_type
num_vertices(const CGAL_LCC_TYPE& cm)
{
  CGAL_LCC_TYPE& cmap = const_cast<CGAL_LCC_TYPE&>(cm);
  return cmap.template attributes<0>().size();
}

// Expression required by the boost::EdgeListGraph concept.

CGAL_LCC_TEMPLATE_ARGS
std::pair<typename boost::graph_traits<CGAL_LCC_TYPE >::edge_iterator,
          typename boost::graph_traits<CGAL_LCC_TYPE >::edge_iterator>
edges(const CGAL_LCC_TYPE& cm)
{
  typedef typename boost::graph_traits<CGAL_LCC_TYPE >::edge_iterator iter_type;
  CGAL_LCC_TYPE& cmap = const_cast<CGAL_LCC_TYPE&>(cm);

  return std::make_pair(iter_type(cmap.darts().begin()),
                        iter_type(cmap.darts().end()));
}

CGAL_LCC_TEMPLATE_ARGS
typename boost::graph_traits<CGAL_LCC_TYPE>::edges_size_type
num_edges(const CGAL_LCC_TYPE& cm)
{ return cm.number_of_darts()/2; }
  
CGAL_LCC_TEMPLATE_ARGS
std::pair<typename boost::graph_traits<CGAL_LCC_TYPE>::edge_descriptor, bool>
edge(typename boost::graph_traits<CGAL_LCC_TYPE>::vertex_descriptor u, 
     typename boost::graph_traits<CGAL_LCC_TYPE>::vertex_descriptor v, 
     const CGAL_LCC_TYPE& cm)
{
  std::pair<typename boost::graph_traits<CGAL_LCC_TYPE>::edge_descriptor,
            bool> res=halfedge(u,v,cm);
  return std::make_pair(EdgeHandle<typename CGAL_LCC_TYPE::Dart_handle>(res.first),
                        res.second);
}

// Expression required by the boost::MutableGraph concept.

CGAL_LCC_TEMPLATE_ARGS
std::pair<typename boost::graph_traits<CGAL_LCC_TYPE>::edge_descriptor,
          bool>
add_edge(typename boost::graph_traits<CGAL_LCC_TYPE>::vertex_descriptor u, 
         typename boost::graph_traits<CGAL_LCC_TYPE>::vertex_descriptor v, 
         CGAL_LCC_TYPE& cm)
{
  typename CGAL_LCC_TYPE::Dart_handle actu = cm.create_dart(u);
  cm.template link_beta<2>(actu, cm.create_dart(v));
  
  return std::make_pair(EdgeHandle<typename CGAL_LCC_TYPE::Dart_handle>(actu),
                        true);
}

CGAL_LCC_TEMPLATE_ARGS
void remove_edge(typename boost::graph_traits<CGAL_LCC_TYPE>::vertex_descriptor u, 
                 typename boost::graph_traits<CGAL_LCC_TYPE>::vertex_descriptor v, 
                 CGAL_LCC_TYPE& cm)
{
  std::pair<typename boost::graph_traits<CGAL_LCC_TYPE>::edge_descriptor, bool>
    e = edge(u,v,cm);
  if ( e.second )
  {
    assert ( !cm.template is_free<2>(e.first.first_halfedge()) );
    cm.restricted_erase_dart(cm.template beta<2>(e.first.first_halfedge()));
    cm.restricted_erase_dart(e.first.first_halfedge());
  }
}

CGAL_LCC_TEMPLATE_ARGS
void
remove_edge(typename boost::graph_traits<CGAL_LCC_TYPE>::edge_descriptor e, 
            CGAL_LCC_TYPE& cm)
{
  assert ( !cm.template is_free<2>(e.first_halfedge()) );
  cm.restricted_erase_dart(cm.template beta<2>(e.first_halfedge()));
  cm.restricted_erase_dart(e);
}

CGAL_LCC_TEMPLATE_ARGS
typename boost::graph_traits<CGAL_LCC_TYPE>::vertex_descriptor
add_vertex(const typename boost::graph_traits<CGAL_LCC_TYPE>::vertex_property_type& p,
           CGAL_LCC_TYPE& cm)
{
  return cm.template create_attribute<0>(p);
}

CGAL_LCC_TEMPLATE_ARGS
typename boost::graph_traits<CGAL_LCC_TYPE>::vertex_descriptor
add_vertex(CGAL_LCC_TYPE& cm)
{
  return cm.template create_attribute<0>();
}

CGAL_LCC_TEMPLATE_ARGS
void
remove_vertex(typename boost::graph_traits<CGAL_LCC_TYPE>::vertex_descriptor v, 
              CGAL_LCC_TYPE& cm)
{
  cm.template erase_attribute<0>(v);
  // Useless because in CMap, attributes are automatically deleted thanks
  // to ref counting
}

CGAL_LCC_TEMPLATE_ARGS
void
clear_vertex(typename boost::graph_traits<CGAL_LCC_TYPE>::vertex_descriptor/* v */, 
             CGAL_LCC_TYPE&/* cm*/)
{
 assert(false); // TODO remove all the edges incident to v
}

//
// FaceGraph
//
CGAL_LCC_TEMPLATE_ARGS
typename boost::graph_traits<CGAL_LCC_TYPE>::halfedge_descriptor
halfedge(typename boost::graph_traits<CGAL_LCC_TYPE>::edge_descriptor e,
         const CGAL_LCC_TYPE&)
{ return e.first_halfedge(); }
  
CGAL_LCC_TEMPLATE_ARGS
typename boost::graph_traits<CGAL_LCC_TYPE>::halfedge_descriptor
halfedge(typename boost::graph_traits<CGAL_LCC_TYPE>::vertex_descriptor v,
         const CGAL_LCC_TYPE&)
{ return v->dart(); }

CGAL_LCC_TEMPLATE_ARGS
std::pair<
  typename boost::graph_traits<CGAL_LCC_TYPE>::halfedge_descriptor,
  bool>
halfedge(typename boost::graph_traits<CGAL_LCC_TYPE>::vertex_descriptor u,
         typename boost::graph_traits<CGAL_LCC_TYPE>::vertex_descriptor v,
         const CGAL_LCC_TYPE& g) 
{
  for (typename CGAL_LCC_TYPE::template Dart_of_cell_range<0>::iterator
         it=const_cast<CGAL_LCC_TYPE&>(g).template
       darts_of_cell<0>(u->dart()).begin(),
         itend=const_cast<CGAL_LCC_TYPE&>(g).template
       darts_of_cell<0>(u->dart()).end();
       /*       darts_of_cell<0>(const_cast<CGAL_LCC_TYPE&>(g).template beta<2>(u->dart())).begin(),
         itend=const_cast<CGAL_LCC_TYPE&>(g).template
         darts_of_cell<0>(const_cast<CGAL_LCC_TYPE&>(g).template beta<2>(u->dart())).end();*/
       it!=itend; ++it)
  {
    if (it->template attribute<0>()==v)
    {
      return std::make_pair(it, true);
      // return std::make_pair(const_cast<CGAL_LCC_TYPE&>(g).template beta<2>(it), true);
    }
  }

  return std::make_pair(g.null_handle, false);
}

CGAL_LCC_TEMPLATE_ARGS
typename boost::graph_traits<CGAL_LCC_TYPE>::halfedge_descriptor
halfedge(typename boost::graph_traits<CGAL_LCC_TYPE>::face_descriptor f,
         const CGAL_LCC_TYPE&) 
{ return f->dart(); }
  
CGAL_LCC_TEMPLATE_ARGS
typename boost::graph_traits<CGAL_LCC_TYPE>::face_descriptor
face(typename boost::graph_traits<CGAL_LCC_TYPE>::halfedge_descriptor h,
     const CGAL_LCC_TYPE& cm)
{ return const_cast<CGAL_LCC_TYPE&>(cm).template attribute<2>(h); }

CGAL_LCC_TEMPLATE_ARGS
typename boost::graph_traits<CGAL_LCC_TYPE>::edge_descriptor
edge(typename boost::graph_traits<CGAL_LCC_TYPE>::halfedge_descriptor h,
     const CGAL_LCC_TYPE&/* cm*/)
{ return EdgeHandle
      <typename boost::graph_traits<CGAL_LCC_TYPE>::halfedge_descriptor>(h); }

CGAL_LCC_TEMPLATE_ARGS
CGAL::Iterator_range<typename boost::graph_traits<CGAL_LCC_TYPE>::face_iterator>
faces(const CGAL_LCC_TYPE& cm)
{
  typedef typename boost::graph_traits<CGAL_LCC_TYPE >::face_iterator iter_type;
  CGAL_LCC_TYPE& cmap = const_cast<CGAL_LCC_TYPE&>(cm);
  return CGAL::make_range(iter_type(cmap.template attributes<2>().begin()),
                         iter_type(cmap.template attributes<2>().end()));
}

CGAL_LCC_TEMPLATE_ARGS
typename boost::graph_traits<CGAL_LCC_TYPE>::edges_size_type
num_faces(const CGAL_LCC_TYPE& cm)
{ return cm.template attributes<2>().size(); }

CGAL_LCC_TEMPLATE_ARGS_NOTEND typename InputIterator>
typename boost::graph_traits<CGAL_LCC_TYPE>::face_descriptor
add_face(InputIterator begin, InputIterator end, CGAL_LCC_TYPE& cm)
{
  while(begin!=end)
  {
    cm.template create_attribute<2>();
    ++begin;
  }
}

CGAL_LCC_TEMPLATE_ARGS
bool is_valid(const CGAL_LCC_TYPE& cm, bool = false)
{
  // cm.display_darts(std::cout,true);
  return cm.is_valid(); // true to inverse the convention between darts and 0-attributes
}

CGAL_LCC_TEMPLATE_ARGS
Iterator_range<typename boost::graph_traits<CGAL_LCC_TYPE>::halfedge_iterator>
halfedges(const CGAL_LCC_TYPE& cm)
{
  typedef typename boost::graph_traits<CGAL_LCC_TYPE>::halfedge_iterator Iter;
  CGAL_LCC_TYPE& cmap = const_cast<CGAL_LCC_TYPE&>(cm);
  return CGAL::make_range(Iter(cmap.darts().begin()), Iter(cmap.darts().end()));
}

CGAL_LCC_TEMPLATE_ARGS
typename boost::graph_traits<CGAL_LCC_TYPE>::halfedges_size_type
num_halfedges(const CGAL_LCC_TYPE& cm)
{ return cm.number_of_darts(); }

CGAL_LCC_TEMPLATE_ARGS
typename boost::graph_traits<CGAL_LCC_TYPE>::edge_descriptor
add_edge(CGAL_LCC_TYPE& cm)
{
  typename CGAL_LCC_TYPE::Dart_handle actu = cm.create_dart();
  cm.template link_beta<2>(actu, cm.create_dart());
  return actu;
}
  
CGAL_LCC_TEMPLATE_ARGS
void set_target(typename boost::graph_traits<CGAL_LCC_TYPE>::halfedge_descriptor h1,
                typename boost::graph_traits<CGAL_LCC_TYPE>::vertex_descriptor v,
                CGAL_LCC_TYPE& cm)
{
  cm.template restricted_set_dart_attribute<0>(cm.template beta<2>(h1), v);
  // cm.template restricted_set_dart_attribute<0>(h1, v);
  // cm.template set_dart_attribute<0>(h1, v);
}

CGAL_LCC_TEMPLATE_ARGS
void set_next(typename boost::graph_traits<CGAL_LCC_TYPE>::halfedge_descriptor h1,
              typename boost::graph_traits<CGAL_LCC_TYPE>::halfedge_descriptor h2,
              CGAL_LCC_TYPE& cm)
{ cm.basic_link_beta_1(h1, h2); }

CGAL_LCC_TEMPLATE_ARGS
void set_halfedge(typename boost::graph_traits<CGAL_LCC_TYPE>::vertex_descriptor v,
                  typename boost::graph_traits<CGAL_LCC_TYPE>::halfedge_descriptor h,
                  CGAL_LCC_TYPE& cm)
{
  // v->set_dart(cm.template beta<2>(h));
  cm.template set_dart_of_attribute<0>(v, h);
}

CGAL_LCC_TEMPLATE_ARGS
typename boost::graph_traits<CGAL_LCC_TYPE>::face_descriptor
add_face(CGAL_LCC_TYPE& cm)
{ return cm.template create_attribute<2>(); }

CGAL_LCC_TEMPLATE_ARGS
void remove_face(typename boost::graph_traits<CGAL_LCC_TYPE>::face_descriptor f,
                 CGAL_LCC_TYPE& cm)
{
  cm.template erase_attribute<2>(f);
  // Useled because in CMap, attributes are automatically deleted thanks
  // to ref counting
}

CGAL_LCC_TEMPLATE_ARGS
void set_face(typename boost::graph_traits<CGAL_LCC_TYPE>::halfedge_descriptor h,              
              typename boost::graph_traits<CGAL_LCC_TYPE>::face_descriptor f,
              CGAL_LCC_TYPE& cm)
{
  cm.template restricted_set_dart_attribute<2>(h, f);
  // cm.template set_dart_attribute<2>(h, f);
}

CGAL_LCC_TEMPLATE_ARGS
void set_halfedge(typename boost::graph_traits<CGAL_LCC_TYPE>::face_descriptor f,
                  typename boost::graph_traits<CGAL_LCC_TYPE>::halfedge_descriptor h,
                  const CGAL_LCC_TYPE&)
{ f->set_dart(h); }

} // namespace CGAL

//#undef CGAL_CMAP_BASE_TEMPLATE_ARGS
//#undef CGAL_CMAP_TEMPLATE_ARGS
//#undef CGAL_CMAP_TYPE
//#undef CGAL_CMAP_BASE_TYPE
#undef CGAL_LCC_TEMPLATE_ARGS
#undef CGAL_LCC_TYPE

#endif // CGAL_BOOST_GRAPH_TRAITS_LINEAR_CELL_COMPLEX_H
