// Copyright (c) 2007  GeometryFactory (France).  All rights reserved.
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
//
// Author(s)     : Andreas Fabri, Fernando Cacciola

#ifndef CGAL_BOOST_GRAPH_GRAPH_TRAITS_HALFEDGEDS_H
#define CGAL_BOOST_GRAPH_GRAPH_TRAITS_HALFEDGEDS_H

#include <functional>

// include this to avoid a VC15 warning
#include <CGAL/boost/graph/named_function_params.h>

#include <boost/config.hpp>
#include <boost/iterator/iterator_adaptor.hpp>
#include <CGAL/boost/iterator/transform_iterator.hpp>
#include <boost/type_traits/remove_const.hpp>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

#include <CGAL/basic.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/Handle_hash_function.h>
#include <CGAL/iterator.h>

#ifndef CGAL_NO_DEPRECATED_CODE
#include <CGAL/boost/graph/halfedge_graph_traits.h>
#endif

namespace CGAL {

namespace internal {

// a HDS_halfedge pretending to be an Edge
template<typename Halfedge_handle>
struct HDS_edge {
  HDS_edge() : halfedge_() { }

  explicit HDS_edge(const Halfedge_handle& h) : halfedge_(h) {}

  bool operator==(const HDS_edge& other) const {
    // equality is tricky, we are equal when both halfedges are the
    // same or when the opposite halfedge of this is the same as
    // halfedge_.other but we need to be careful not to apply this
    // should this be the invalid halfedge or opposite is going to
    // blow up
    if(halfedge_ == other.halfedge_) {
      return true;
    } else if(halfedge_ != Halfedge_handle()) { // not default constructed
      return halfedge_->opposite() == other.halfedge_;
    } else {
      // this is the invalid halfedge, it can only be equal to the
      // invalid halfedge and this is covered by the first case
      return false;
    }
  }

  bool operator!=(const HDS_edge& other) const {
    return !(*this == other);
  }

  friend bool operator<(const HDS_edge& a,const HDS_edge& b)
  {
    if(a==b) return false;
    return a.halfedge_ < b.halfedge_;
  }

  // forward some function to avoid boilerplate and typedefs inside
  // the free functions

  HDS_edge next() { return HDS_edge(halfedge_->next()); }

  // this is potentially broken as it does not use the decorator to
  // find prev, but we cannot instantiate the entire decorator
  // without the full polyhedron type and taking all necessary
  // template parameters seems overkill
  HDS_edge prev() { return HDS_edge(halfedge_->prev()); }

  HDS_edge opposite() { return HDS_edge(halfedge_->opposite()); }

  // this is hacky, we don't know the actual type of the id and if we
  // start adding decltype special cases we have to do it consistently
  // up to the property map and maybe back down to Polyhedron.
  std::size_t id() const { return halfedge_->id() / 2; }

  Halfedge_handle halfedge() const { return halfedge_; }

  // save us some work to do function chaining
  HDS_edge
  opposite_next() { return HDS_edge(halfedge_->opposite()->next()); }

  HDS_edge
  next_opposite() { return HDS_edge(halfedge_->next()->opposite()); }

  HDS_edge
  prev_opposite() { return HDS_edge(halfedge_->prev()->opposite()); }

  HDS_edge
  opposite_prev() { return HDS_edge(halfedge_->opposite()->prev()); }

  friend  std::size_t hash_value(const HDS_edge&  i)
  {
    if (i.halfedge()==Halfedge_handle()) return 0;
    return hash_value(i.halfedge()<i.halfedge()->opposite()?
                      i.halfedge():i.halfedge()->opposite());
  }

private:
  Halfedge_handle halfedge_;
};

// make edge_descriptor hashable by default in Unique_hash_map
namespace handle{
  template<typename Halfedge_handle>
  struct Hash_functor< HDS_edge<Halfedge_handle> >
  {
    std::size_t
    operator()(const HDS_edge<Halfedge_handle>& edge)
    {
      Halfedge_handle he = edge.halfedge();
      if (he==Halfedge_handle()) return 0;
      if ( he < he->opposite() )
        return Hash_functor<Halfedge_handle>()(he);
      return Hash_functor<Halfedge_handle>()(he->opposite());
    }
  };
} //end of namespace handle

template<typename Halfedge_handle>
struct Construct_edge {
  typedef HDS_edge<Halfedge_handle> result_type;
  HDS_edge<Halfedge_handle> operator()(const Halfedge_handle& he) const
  { return HDS_edge<Halfedge_handle>(he); }
};

template<typename Halfedge_handle>
struct Construct_edge_opposite {
  typedef HDS_edge<Halfedge_handle> result_type;
  HDS_edge<Halfedge_handle> operator()(const Halfedge_handle& he) const
  { return HDS_edge<Halfedge_handle>(he->opposite()); }
};

} // internal

template <class HDS>
struct HDS_graph_traits
{
private:
  struct HDS_graph_traversal_category : public virtual boost::bidirectional_graph_tag,
                                        public virtual boost::vertex_list_graph_tag,
                                        public virtual boost::edge_list_graph_tag
  {};

public:
  typedef typename HDS::Vertex_handle                                vertex_descriptor;
  typedef typename internal::HDS_edge<typename HDS::Halfedge_handle> edge_descriptor;
  typedef typename HDS::Face_handle                                  face_descriptor;
  typedef typename HDS::Halfedge_handle                              halfedge_descriptor;

  typedef CGAL::Prevent_deref<typename HDS::Vertex_iterator>     vertex_iterator;
  typedef CGAL::Prevent_deref<typename HDS::Face_iterator>       face_iterator;
  typedef CGAL::Prevent_deref<typename HDS::Edge_iterator>       edge_iterator_i;
  typedef CGAL::Prevent_deref<typename HDS::Halfedge_iterator>   halfedge_iterator;


  typedef boost::transform_iterator<
    internal::Construct_edge<typename HDS::Halfedge_handle>,
    edge_iterator_i,
    edge_descriptor> edge_iterator;

  typedef Out_edge_iterator<HDS> out_edge_iterator;

  typedef In_edge_iterator<HDS> in_edge_iterator;

  typedef boost::undirected_tag             directed_category;
  typedef boost::disallow_parallel_edge_tag edge_parallel_category;
  typedef HDS_graph_traversal_category      traversal_category;

  typedef typename HDS::size_type vertices_size_type;
  typedef vertices_size_type      edges_size_type;
  typedef vertices_size_type      halfedges_size_type;
  typedef vertices_size_type      degree_size_type;
  typedef vertices_size_type      faces_size_type;

  static vertex_descriptor null_vertex() { return vertex_descriptor(); }
  static halfedge_descriptor null_halfedge() { return halfedge_descriptor(); }
  static face_descriptor null_face() { return face_descriptor(); }
};


} //namespace CGAL


namespace std {

#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4099) // For VC10 it is class hash 
#endif

#ifndef CGAL_CFG_NO_STD_HASH

  template <typename H>
  struct hash<CGAL::internal::HDS_edge<H> > {
    std::size_t operator()(const CGAL::internal::HDS_edge<H>& e) const
    {
      std::hash<H> fct;
      if (e.halfedge()==H()) return 0;
      return fct(e.halfedge()<e.halfedge()->opposite()?
                 e.halfedge():e.halfedge()->opposite());
    }
  };

#endif // CGAL_CFG_NO_STD_HASH

#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif

} // namespace std

#endif // CGAL_BOOST_GRAPH_GRAPH_TRAITS_HALFEDGEDS_H
