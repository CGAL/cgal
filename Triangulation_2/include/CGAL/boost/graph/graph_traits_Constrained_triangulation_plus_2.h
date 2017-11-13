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
// 
//
// Author(s)     : Andreas Fabri, Fernando Cacciola

#ifndef CGAL_GRAPH_TRAITS_CONSTRAINED_TRIANGULATION_PLUS_2_H
#define CGAL_GRAPH_TRAITS_CONSTRAINED_TRIANGULATION_PLUS_2_H

#include <CGAL/license/Triangulation_2.h>


// include this to avoid a VC15 warning
#include <CGAL/boost/graph/named_function_params.h>

#include <boost/config.hpp>
#include <boost/iterator_adaptors.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>
#include <CGAL/boost/graph/graph_traits_Constrained_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>

// The functions and classes in this file allows the user to
// treat a CGAL Constrained_triangulation_plus_2 object as a boost graph "as is". No
// wrapper is needed for the Constrained_triangulation_plus_2 object.



namespace boost { 

  template <class Tr>
  struct graph_traits< CGAL::Constrained_triangulation_plus_2<Tr> > {

    struct DT2_graph_traversal_category : 
      public virtual bidirectional_graph_tag,
      public virtual adjacency_graph_tag,        
      public virtual edge_list_graph_tag,
      public virtual vertex_list_graph_tag { };

    typedef typename boost::graph_traits<Tr>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<Tr>::face_descriptor face_descriptor;
    typedef typename boost::graph_traits<Tr>::edge_descriptor  edge_descriptor;
    typedef typename boost::graph_traits<Tr>::edge_iterator edge_iterator;

    typedef typename boost::graph_traits<Tr>::halfedge_descriptor halfedge_descriptor;

    typedef typename boost::graph_traits<Tr>::halfedge_iterator halfedge_iterator;
    typedef typename boost::graph_traits<Tr>::vertex_iterator vertex_iterator;
    typedef typename boost::graph_traits<Tr>::face_iterator face_iterator;
    typedef typename boost::graph_traits<Tr>::out_edge_iterator out_edge_iterator;
    typedef typename boost::graph_traits<Tr>::in_edge_iterator in_edge_iterator;
    typedef typename boost::graph_traits<Tr>::Incident_vertices_iterator Incident_vertices_iterator;
    typedef undirected_tag directed_category;
    typedef disallow_parallel_edge_tag edge_parallel_category; 
    typedef DT2_graph_traversal_category traversal_category;
    typedef typename boost::graph_traits<Tr>::size_type size_type;
    typedef size_type vertices_size_type;
    typedef size_type edges_size_type;
    typedef size_type degree_size_type;
  };


} // namespace boost

namespace CGAL {



  template <class Tr>
  typename boost::graph_traits< CGAL::Constrained_triangulation_plus_2<Tr> >::vertex_descriptor
  source(typename boost::graph_traits< CGAL::Constrained_triangulation_plus_2<Tr> >::edge_descriptor e,
         const CGAL::Constrained_triangulation_plus_2<Tr>& g)
  {
    return e.first->vertex(g.ccw(e.second));
  }

  template <class Tr>
  typename boost::graph_traits< CGAL::Constrained_triangulation_plus_2<Tr> >::vertex_descriptor
  target(typename boost::graph_traits< CGAL::Constrained_triangulation_plus_2<Tr> >::edge_descriptor e,
         const CGAL::Constrained_triangulation_plus_2<Tr>& g)
  {
    return e.first->vertex(g.cw(e.second));
  }


  template <class Tr>
  inline std::pair<
    typename boost::graph_traits< CGAL::Constrained_triangulation_plus_2<Tr> >::out_edge_iterator,
    typename boost::graph_traits< CGAL::Constrained_triangulation_plus_2<Tr> >::out_edge_iterator >  
  out_edges(
    typename boost::graph_traits< CGAL::Constrained_triangulation_plus_2<Tr> >::vertex_descriptor u, 
    const CGAL::Constrained_triangulation_plus_2<Tr>& g)
  {
    typename CGAL::Constrained_triangulation_plus_2<Tr>::Edge_circulator ec(u,u->face());
    typename boost::graph_traits< CGAL::Constrained_triangulation_plus_2<Tr> >::degree_size_type out_deg = out_degree(u,g);
    typedef typename boost::graph_traits< CGAL::Constrained_triangulation_plus_2<Tr> >
      ::out_edge_iterator Iter;
    
    return std::make_pair( Iter(ec), Iter(ec,out_deg) );
  }

  template <class Tr>
  inline std::pair<
    typename boost::graph_traits< CGAL::Constrained_triangulation_plus_2<Tr> >::in_edge_iterator,
    typename boost::graph_traits< CGAL::Constrained_triangulation_plus_2<Tr> >::in_edge_iterator >  
  in_edges(
    typename boost::graph_traits< CGAL::Constrained_triangulation_plus_2<Tr> >::vertex_descriptor u, 
    const CGAL::Constrained_triangulation_plus_2<Tr>& g)
  {
    typename CGAL::Constrained_triangulation_plus_2<Tr>::Edge_circulator ec(u,u->face());
    typename boost::graph_traits< CGAL::Constrained_triangulation_plus_2<Tr> >::degree_size_type out_deg = out_degree(u,g);
    typedef typename boost::graph_traits< CGAL::Constrained_triangulation_plus_2<Tr> >
      ::in_edge_iterator Iter;
    return std::make_pair( Iter(ec), Iter(ec,out_deg) );
  }



  // property maps
  template <class Tr>
  class CTP2_vertex_id_map
    : public boost::put_get_helper<int, CTP2_vertex_id_map<Tr> >
  {
  public:
    typedef boost::readable_property_map_tag category;
    typedef int value_type;
    typedef int reference;
    typedef typename CGAL::Constrained_triangulation_plus_2<Tr>::Vertex_handle key_type;
    
    CTP2_vertex_id_map()
    {}
    
    long operator[](key_type vh) const {
      return vh->id(); 
    }
  };

 template <class Tr>
  class CTP2_vertex_point_map
  {
  public:
    typedef boost::lvalue_property_map_tag category;
    typedef typename CGAL::Constrained_triangulation_plus_2<Tr>::Point value_type;
    typedef value_type& reference;
    typedef typename CGAL::Constrained_triangulation_plus_2<Tr>::Vertex_handle key_type;

    friend reference get(CTP2_vertex_point_map<Tr>, key_type vh)
    { 
      return vh->point(); 
    }
    friend void put(CTP2_vertex_point_map<Tr>, key_type vh, reference v)
    {
      vh->point()=v; 
    }
    reference operator[](key_type vh) const {
      return vh->point();
    }
  };

  template <class Tr>
  class CTP2_edge_id_map
    : public boost::put_get_helper<int, CTP2_edge_id_map<Tr> >
  {
  public:
    typedef boost::readable_property_map_tag category;
    typedef int value_type;
    typedef int reference;
    typedef typename CGAL::Constrained_triangulation_plus_2<Tr>::Edge key_type;
    
    CTP2_edge_id_map()
    {}
    
    long operator[](key_type e) const {
      return (3 * e.first.id()) + e.second; 
    }
  };

  template <class Tr>
  class CTP2_edge_weight_map
    : public boost::put_get_helper<typename Tr::Geom_traits::FT, CTP2_edge_weight_map<Tr> >
  {
  private:
    const CGAL::Constrained_triangulation_plus_2<Tr>& tr;
  public:
    typedef boost::readable_property_map_tag category;
    typedef typename Tr::Geom_traits::FT value_type;
    typedef value_type reference;
    typedef typename CGAL::Constrained_triangulation_plus_2<Tr>::Edge key_type;

    CTP2_edge_weight_map(const CGAL::Constrained_triangulation_plus_2<Tr>& tr_) 
      : tr(tr_) 
    { }

    typename Tr::Geom_traits::FT operator[](key_type e) const {
      return approximate_sqrt(tr.segment(e).squared_length());
    }
  };


  template <class Tr>
  inline CTP2_vertex_id_map<Tr>
  get(boost::vertex_index_t, const CGAL::Constrained_triangulation_plus_2<Tr>& ) {
    CTP2_vertex_id_map<Tr> m;
    return m;
  }

  template <class Tr>
  inline CTP2_vertex_point_map<Tr>
  get(boost::vertex_point_t, const CGAL::Constrained_triangulation_plus_2<Tr>& ) {
    CTP2_vertex_point_map<Tr> m;
    return m;
  }
  
  template <class Tr>
  inline CTP2_edge_id_map<Tr>
  get(boost::edge_index_t, const CGAL::Constrained_triangulation_plus_2<Tr>& ) {
    CTP2_edge_id_map<Tr> m;
    return m;
  }

  template <class Tr>
  inline CTP2_edge_weight_map<Tr>
  get(boost::edge_weight_t, const CGAL::Constrained_triangulation_plus_2<Tr>& g) {
    CTP2_edge_weight_map<Tr> m(g);
    return m;
  }

  template <class Tag>
  struct CTP2_property_map { };

  template <>
  struct CTP2_property_map<boost::vertex_index_t> {
    template <class Tr>
    struct bind_ {
      typedef CTP2_vertex_id_map<Tr> type;
      typedef CTP2_vertex_id_map<Tr> const_type;
    };
  };

  template <>
  struct CTP2_property_map<boost::vertex_point_t> {
    template <class Tr>
    struct bind_ {
      typedef CTP2_vertex_point_map<Tr> type;
      typedef CTP2_vertex_point_map<Tr> const_type;
    };
  };



  template <>
  struct CTP2_property_map<boost::edge_index_t> {
    template <class Tr>
    struct bind_ {
      typedef CTP2_edge_id_map<Tr> type;
      typedef CTP2_edge_id_map<Tr> const_type;
    };
  };


  template <>
  struct CTP2_property_map<boost::edge_weight_t> {
    template <class Tr>
    struct bind_ {
      typedef CTP2_edge_weight_map<Tr> type;
      typedef CTP2_edge_weight_map<Tr> const_type;
    };
  };

} // namespace CGAL

namespace boost {

  // g++ 'enumeral_type' in template unification not implemented workaround
  template <class Tr, class Tag>
  struct property_map<CGAL::Constrained_triangulation_plus_2<Tr>, Tag> {
    typedef typename 
    CGAL::CTP2_property_map<Tag>::template bind_<Tr> map_gen;
    typedef typename map_gen::type type;
    typedef typename map_gen::const_type const_type;
  };

  // see struct property_map in Polyhedron for an explanation
  template <class Tr, class Tag>
  struct property_map<const CGAL::Constrained_triangulation_plus_2<Tr>, Tag> {
    typedef typename 
    CGAL::CTP2_property_map<Tag>::template bind_<Tr> map_gen;
    typedef typename map_gen::type type;
    typedef typename map_gen::const_type const_type;
  };


  // What are those needed for ???
  template <typename Tr>
  struct edge_property_type<CGAL::Constrained_triangulation_plus_2<Tr> > {
    typedef void type;
  };  

  template <typename Tr>
  struct vertex_property_type<CGAL::Constrained_triangulation_plus_2<Tr> > {
    typedef void type;
  };
} // namespace boost


#endif // CGAL_GRAPH_TRAITS_CONSTRAINED_TRIANGULATION_PLUS_2_H
