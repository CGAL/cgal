// Copyright (c) 2015  GeometryFactory (France).  All rights reserved.
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
// Author(s)     : Andreas Fabri

#ifndef CGAL_BGL_DUAL_H
#define CGAL_BGL_DUAL_H

#include <CGAL/boost/graph/properties.h>
#include <boost/range/distance.hpp>
#include <CGAL/boost/graph/iterator.h>

namespace CGAL {

/*!
\ingroup PkgBGLAdaptors

The class template `Dual` is an adaptor that creates the dual view of
a `FaceGraph`. Faces of the original graph correspond to vertices in
the `Dual` and vice versa.

Note that border edges in a `Dual` have the `null_face` of the
original graph as either source or target. This is unusual and might
break other algorithms since edges are always assumed to have non-null
vertices as a source and target. It is possible to filter border edges
using `boost::filtered_graph` as shown in example
\ref BGL_surface_mesh/surface_mesh_dual.cpp

Property forwarding
-------------------
\cgalAdvancedBegin
Edge properties of the underlying graph are forwarded directly. For
faces and vertices only the `face_index` and `vertex_index` properties
are forwarded. Accessing other properties will lead to a compilation
error.
\cgalAdvancedEnd

\tparam Primal_ must be a model of `FaceGraph`

\cgalModels `FaceGraph`

*/
template <typename Primal_>
class Dual
{
  const Primal_& primal_;

public:
  /*! The underlying primal type. */
  typedef Primal_ Primal;

  /*! Construct a Dual from a given primal. */
  Dual(const Primal& primal)
    : primal_(primal) {}

  /*! Returns the underlying primal. */
  const Primal& primal() const
  { return primal_; }
};


/*!
  Construct a `Dual` from a given `primal`.
  \relates CGAL::Dual
 */
template<typename Primal>
Dual<Primal> dual(const Primal& primal)
{ return Dual<Primal>(primal); }

} // namespace CGAL

namespace boost {
  /*!
\ingroup PkgBGLTraits

  */
template <typename Primal>
class graph_traits<CGAL::Dual<Primal> >
{
public:
  typedef boost::graph_traits<Primal> GTP;
  typedef typename GTP::face_descriptor     vertex_descriptor;
  typedef typename GTP::vertex_descriptor   face_descriptor;
  typedef typename GTP::halfedge_descriptor halfedge_descriptor;
  typedef typename GTP::edge_descriptor     edge_descriptor;
  typedef typename GTP::directed_category   directed_category;
  typedef boost::allow_parallel_edge_tag    edge_parallel_category; 
  typedef typename GTP::traversal_category  traversal_category;

  typedef typename GTP::faces_size_type          vertices_size_type;
  typedef typename GTP::vertices_size_type       faces_size_type;
  typedef typename GTP::edges_size_type          edges_size_type;
  typedef typename GTP::halfedges_size_type      halfedges_size_type;
  typedef typename GTP::degree_size_type         degree_size_type;

  typedef typename GTP::face_iterator     vertex_iterator;
  typedef typename GTP::vertex_iterator   face_iterator;
  typedef typename GTP::halfedge_iterator halfedge_iterator;
  typedef typename GTP::edge_iterator     edge_iterator;

  typedef CGAL::Edge_around_face_iterator<Primal>      out_edge_iterator;
  typedef CGAL::Opposite_edge_around_face_iterator<Primal> in_edge_iterator;

  static vertex_descriptor   null_vertex()   { return vertex_descriptor(); }
  static face_descriptor     null_face()     { return face_descriptor(); }
  static halfedge_descriptor null_halfedge() { return halfedge_descriptor(); }
};
 
template<typename P>
struct graph_traits< const CGAL::Dual<P> >  
  : public graph_traits< CGAL::Dual<P> >
{}; 

namespace internal{

template <class G>
struct Dual_vertex_index_pmap{
  typedef typename boost::property_map<G, boost::face_index_t>::type  Property_map;
  Property_map m_pmap;

  typedef typename boost::graph_traits<G>::face_descriptor key_type;
  typedef typename Property_map::value_type value_type;
  typedef typename Property_map::reference reference;
  typedef typename Property_map::category category;

  Dual_vertex_index_pmap(const G& g)
    : m_pmap( get(boost::face_index, g) )
  {}

  friend reference get(const Dual_vertex_index_pmap<G>& pmap, key_type fd) {
    return get(pmap.m_pmap, fd);
  }
};

template <class G>
struct Dual_face_index_pmap{
  typedef typename boost::property_map<G, boost::vertex_index_t>::type  Property_map;
  Property_map m_pmap;

  typedef typename boost::graph_traits<G>::vertex_descriptor key_type;
  typedef typename Property_map::value_type value_type;
  typedef typename Property_map::reference reference;
  typedef typename Property_map::category category;

  Dual_face_index_pmap(const G& g)
    : m_pmap( get(boost::vertex_index, g) )
  {}

  friend reference get(const Dual_face_index_pmap<G>& pmap, key_type vd) {
    return get(pmap.m_pmap, vd);
  }
};

template<typename P, typename Property,
         bool is_edge = boost::is_same<boost::edge_property_tag,
                                       typename boost::property_kind<Property>::type>::value>
struct Dual_property_maps : boost::property_map<P, Property> {};

template< typename P, typename Property>
struct Dual_property_maps<P, Property, false> {};

} //end of namespace internal

template <typename P, typename Property>
struct property_map<CGAL::Dual<P>, Property>
  : internal::Dual_property_maps<P, Property> {};

template <typename P>
struct property_map<CGAL::Dual<P>, boost::vertex_index_t>
{
  typedef internal::Dual_vertex_index_pmap<P> type;
  typedef internal::Dual_vertex_index_pmap<P> const_type;
};

template <typename P>
struct property_map<CGAL::Dual<P>, boost::face_index_t>
{
  typedef internal::Dual_face_index_pmap<P> type;
  typedef internal::Dual_face_index_pmap<P> const_type;
};

} // namespace boost

namespace CGAL {

template <typename P, typename Property>
typename boost::property_map<P, Property>::type
get(Property p, Dual<P>& dual)
{
  return get(p, dual.primal());
}

template <typename P, typename Property>
typename boost::property_map<P, Property>::const_type
get(Property p, const Dual<P>& dual)
{
  return get(p, dual.primal());
}

template <typename P, typename Property, typename Key >
typename boost::property_map_value<P, Property>::type
get(Property p, const Dual<P>& dual, const Key& k)
{
  return get(p, dual.primal(), k);
}

template<typename G, typename P, typename>
struct Property_map_value_dummy {
  typedef typename boost::property_map_value<G, P>::type type;
};

template <typename P, typename Key>
typename Property_map_value_dummy<Dual<P>, boost::vertex_index_t, Key>::type
get(boost::vertex_index_t, const Dual<P>& dual, const Key& k)
{
  return get(typename boost::internal::Dual_vertex_index_pmap<P>(dual.primal()), k);
}

template <typename P, typename Key>
typename Property_map_value_dummy<Dual<P>, boost::face_index_t, Key>::type
get(boost::face_index_t, const Dual<P>& dual, const Key& k)
{
  return get(typename boost::internal::Dual_face_index_pmap<P>(dual.primal()), k);
}

template <typename P, typename Property, typename Key, typename Value>
void
put(Property p, const Dual<P>& dual, const Key& k, const Value& val)
{
  put(p, dual.primal(), k, val);
}

template <typename P>
typename boost::internal::Dual_vertex_index_pmap<P>
get(boost::vertex_index_t, const Dual<P>& dual)
{
  return typename boost::internal::Dual_vertex_index_pmap<P>(dual.primal());
}

template <typename P>
typename boost::internal::Dual_face_index_pmap<P>
get(boost::face_index_t, const Dual<P>& dual)
{
  return typename boost::internal::Dual_face_index_pmap<P>(dual.primal());
}

template <typename P>
typename boost::graph_traits<CGAL::Dual<P> >::vertices_size_type
num_vertices(const CGAL::Dual<P>& dual)
{
  return num_faces(dual.primal());
}

template <typename P>
typename boost::graph_traits<CGAL::Dual<P> >::edges_size_type
num_edges(const CGAL::Dual<P>& dual)
{
  return num_edges(dual.primal());
}

template <typename P>
typename boost::graph_traits<CGAL::Dual<P> >::halfedges_size_type
num_halfedges(const CGAL::Dual<P>& dual)
{
  return num_halfedges(dual.primal());
}
     
template <typename P>
typename boost::graph_traits<CGAL::Dual<P> >::faces_size_type
num_faces(const CGAL::Dual<P>& dual)
{
  return num_vertices(dual.primal());
}

template <typename P>
Iterator_range<typename boost::graph_traits<Dual<P> >::vertex_iterator>
vertices(const CGAL::Dual<P>& dual)
{
  return faces(dual.primal()); 
}

template <typename P>
Iterator_range<typename boost::graph_traits<Dual<P> >::face_iterator>
faces(const CGAL::Dual<P>& dual)
{
  return vertices(dual.primal()); 
}
    
template <typename P>
Iterator_range<typename boost::graph_traits<Dual<P> >::halfedge_iterator>
halfedges(const CGAL::Dual<P>& dual)
{
  return halfedges(dual.primal()); 
}
  
template <typename P>
Iterator_range<typename boost::graph_traits<Dual<P> >::edge_iterator>
edges(const CGAL::Dual<P>& dual)
{
  return edges(dual.primal()); 
}

template <typename P>
std::pair<typename boost::graph_traits<Dual<P> >::edge_descriptor, bool>
edge(typename boost::graph_traits<Dual<P> >::vertex_descriptor u, 
     typename boost::graph_traits<Dual<P> >::vertex_descriptor v, 
     const Dual<P>& dual)
{
  typename boost::graph_traits<Dual<P> >::out_edge_iterator e, e_end;
  for(boost::tie(e, e_end) = out_edges(u, dual); e != e_end; ++e) {
    if(target(*e, dual) == v)
      return std::make_pair(*e, true);
  }
  
  return std::make_pair(typename boost::graph_traits<Dual<P> >::edge_descriptor(), false);
}

template <typename P>
typename boost::graph_traits<Dual<P> >::edge_descriptor
edge(typename boost::graph_traits<Dual<P> >::halfedge_descriptor h,
     const Dual<P>& dual)
{
  return edge(h, dual.primal());
}

template <typename P>
typename boost::graph_traits<Dual<P> >::vertex_descriptor
source(typename boost::graph_traits<Dual<P> >::halfedge_descriptor h,
       const Dual<P>& dual)
{
  const typename Dual<P>::Primal& primal = dual.primal();
  return face(h,primal);
}
 
template <typename P>
typename boost::graph_traits<Dual<P> >::vertex_descriptor
target(typename boost::graph_traits<Dual<P> >::halfedge_descriptor h,
       const Dual<P>& dual)
{
  const typename Dual<P>::Primal& primal = dual.primal();
  return face(opposite(h,primal),primal);
}
 

template <typename P>
typename boost::graph_traits<Dual<P> >::vertex_descriptor
source(typename boost::graph_traits<Dual<P> >::edge_descriptor h,
       const Dual<P>& dual)
{
  const typename Dual<P>::Primal& primal = dual.primal();
  return face(halfedge(h,primal),primal);
}
 
template <typename P>
typename boost::graph_traits<Dual<P> >::vertex_descriptor
target(typename boost::graph_traits<Dual<P> >::edge_descriptor h,
       const Dual<P>& dual)
{
  const typename Dual<P>::Primal& primal = dual.primal();
  return face(opposite(halfedge(h,primal),primal),primal);
}

template <typename P>
typename boost::graph_traits<Dual<P> >::halfedge_descriptor
halfedge(typename boost::graph_traits<Dual<P> >::vertex_descriptor v,
         const Dual<P>& dual)
{
  const typename Dual<P>::Primal& primal = dual.primal();
  return opposite(halfedge(v, primal),primal);
}

template <typename P>
typename boost::graph_traits<Dual<P> >::halfedge_descriptor
halfedge(typename boost::graph_traits<Dual<P> >::face_descriptor f,
         const Dual<P>& dual)
{
  return halfedge(f, dual.primal());
}

template <typename P>
typename boost::graph_traits<Dual<P> >::halfedge_descriptor
halfedge(typename boost::graph_traits<Dual<P> >::edge_descriptor e,
         const Dual<P>& dual)
{
  return halfedge(e, dual.primal());
}

template <typename P>
std::pair<typename boost::graph_traits<Dual<P> >::halfedge_descriptor, bool>
halfedge(typename boost::graph_traits<Dual<P> >::vertex_descriptor u,
         typename boost::graph_traits<Dual<P> >::vertex_descriptor v,
         const Dual<P>& dual)
{
  typename boost::graph_traits<Dual<P> >::out_edge_iterator e, e_end;
  for(boost::tie(e, e_end) = out_edges(u, dual); e != e_end; ++e) {
    if(target(*e, dual) == v)
      return std::make_pair(halfedge(*e, dual), true);
  }
  
  return std::make_pair(boost::graph_traits<Dual<P> >::null_halfedge(), false);
}

template <typename P>
typename boost::graph_traits<Dual<P> >::face_descriptor
face(typename boost::graph_traits<Dual<P> >::halfedge_descriptor h,
       const Dual<P>& dual)
{
  const typename Dual<P>::Primal& primal = dual.primal();
  return target(h,primal);
}

template <typename P>
typename boost::graph_traits<Dual<P> >::halfedge_descriptor
opposite(typename boost::graph_traits<Dual<P> >::halfedge_descriptor h,
         const Dual<P>& dual)
{
  return opposite(h, dual.primal());
}

template <typename P>
typename boost::graph_traits<Dual<P> >::halfedge_descriptor
next(typename boost::graph_traits<Dual<P> >::halfedge_descriptor h,
         const Dual<P>& dual)
{
  const typename Dual<P>::Primal& primal = dual.primal();
  return prev(opposite(h,primal),primal);
}  

template <typename P>
typename boost::graph_traits<Dual<P> >::halfedge_descriptor
prev(typename boost::graph_traits<Dual<P> >::halfedge_descriptor h,
         const Dual<P>& dual)
{
  const typename Dual<P>::Primal& primal = dual.primal();
  return opposite(next(h,primal),primal);
}  

template <typename P>
Iterator_range<typename boost::graph_traits<Dual<P> >::out_edge_iterator>
out_edges(typename boost::graph_traits<Dual<P> >::vertex_descriptor v,
          const Dual<P>& dual)
{
  const typename Dual<P>::Primal& primal = dual.primal();
  return edges_around_face(halfedge(v,primal),primal);
}

template <typename P>
Iterator_range<typename boost::graph_traits<Dual<P> >::in_edge_iterator>
in_edges(typename boost::graph_traits<Dual<P> >::vertex_descriptor v,
         const Dual<P>& dual)
{
  const typename Dual<P>::Primal& primal = dual.primal();
  return opposite_edges_around_face(halfedge(v,primal),primal);
}
       
template <typename P>
typename boost::graph_traits<Dual<P> >::degree_size_type
out_degree(typename boost::graph_traits<Dual<P> >::vertex_descriptor v,
           const Dual<P>& dual)
{
  const typename Dual<P>::Primal& primal = dual.primal();
  return boost::distance(halfedges_around_face(halfedge(v,primal),primal));
}

 template <typename P>
typename boost::graph_traits<Dual<P> >::degree_size_type
in_degree(typename boost::graph_traits<Dual<P> >::vertex_descriptor v,
           const Dual<P>& dual)
{
  return out_degree(v,dual);
}
         
        
} // namespace CGAL

#endif
