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

#ifndef CGAL_GRAPH_TRAITS_DELAUNAY_TRIANGULATION_2_H
#define CGAL_GRAPH_TRAITS_DELAUNAY_TRIANGULATION_2_H

#include <boost/config.hpp>
#include <boost/iterator_adaptors.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>
#include <CGAL/boost/graph/graph_traits_Triangulation_2.h>

// The functions and classes in this file allows the user to
// treat a CGAL Delaunay_triangulation_2 object as a boost graph "as is". No
// wrapper is needed for the Delaunay_triangulation_2 object.



namespace boost { 

  template <class GT, class TDS>
  struct graph_traits< CGAL::Delaunay_triangulation_2<GT,TDS> > {

    struct DT2_graph_traversal_category : 
      public virtual bidirectional_graph_tag,
      public virtual adjacency_graph_tag,        
      public virtual edge_list_graph_tag,
      public virtual vertex_list_graph_tag { };

    typedef CGAL::Delaunay_triangulation_2<GT,TDS> Delaunay_triangulation;

    typedef typename CGAL::Delaunay_triangulation_2<GT,TDS>::Vertex_handle vertex_descriptor;
    typedef CGAL::detail::Edge<CGAL::Delaunay_triangulation_2<GT,TDS>, typename CGAL::Delaunay_triangulation_2<GT,TDS>::Edge>  edge_descriptor;
    typedef typename CGAL::Delaunay_triangulation_2<GT,TDS>::All_edges_iterator  edge_iterator;

    typedef CGAL::detail::boost_all_vertices_iterator<Delaunay_triangulation> vertex_iterator;
    typedef CGAL::Counting_iterator<CGAL::detail::Out_edge_circulator<typename Delaunay_triangulation::Edge_circulator, edge_descriptor>, edge_descriptor > out_edge_iterator;
    typedef CGAL::Counting_iterator<CGAL::detail::In_edge_circulator<typename Delaunay_triangulation::Edge_circulator, edge_descriptor>, edge_descriptor > in_edge_iterator;
    typedef CGAL::Counting_iterator<typename Delaunay_triangulation::Vertex_circulator> Incident_vertices_iterator;
    typedef Incident_vertices_iterator adjacency_iterator;

    typedef undirected_tag directed_category;
    typedef disallow_parallel_edge_tag edge_parallel_category; 
    typedef DT2_graph_traversal_category traversal_category;
    typedef typename Delaunay_triangulation::size_type size_type;
    typedef size_type vertices_size_type;
    typedef size_type edges_size_type;
    typedef size_type degree_size_type;
  };


} // namespace boost


namespace CGAL {

  template <class Gt, class Tds>
  typename boost::graph_traits< CGAL::Delaunay_triangulation_2<Gt,Tds> >::vertex_descriptor
  source(typename boost::graph_traits< CGAL::Delaunay_triangulation_2<Gt,Tds> >::edge_descriptor e,
         const CGAL::Delaunay_triangulation_2<Gt,Tds>& g)
  {
    return e.first->vertex(g.ccw(e.second));
  }

  template <class Gt, class Tds>
  typename boost::graph_traits< CGAL::Delaunay_triangulation_2<Gt,Tds> >::vertex_descriptor
  target(typename boost::graph_traits< CGAL::Delaunay_triangulation_2<Gt,Tds> >::edge_descriptor e,
         const CGAL::Delaunay_triangulation_2<Gt,Tds>& g)
  {
    return e.first->vertex(g.cw(e.second));
  }

  template <class Gt, class Tds>
  inline std::pair<
    typename boost::graph_traits< CGAL::Delaunay_triangulation_2<Gt,Tds> >::vertex_iterator,
    typename boost::graph_traits< CGAL::Delaunay_triangulation_2<Gt,Tds> >::vertex_iterator >  
  vertices(const CGAL::Delaunay_triangulation_2<Gt,Tds>& g)
  {
    typedef typename boost::graph_traits< CGAL::Delaunay_triangulation_2<Gt,Tds> >::vertex_iterator
      Iter;
    return std::make_pair( Iter(g.all_vertices_begin()), Iter(g.all_vertices_end()) );
  }


  template <class Gt, class Tds>
  inline std::pair<
    typename boost::graph_traits< CGAL::Delaunay_triangulation_2<Gt,Tds> >::edge_iterator,
    typename boost::graph_traits< CGAL::Delaunay_triangulation_2<Gt,Tds> >::edge_iterator >  
  edges(const CGAL::Delaunay_triangulation_2<Gt,Tds>& g)
  {    
    return std::make_pair(g.all_edges_begin(), g.all_edges_end());
  }

  template <class Gt, class Tds>
  inline std::pair<
    typename boost::graph_traits< CGAL::Delaunay_triangulation_2<Gt,Tds> >::out_edge_iterator,
    typename boost::graph_traits< CGAL::Delaunay_triangulation_2<Gt,Tds> >::out_edge_iterator >  
  out_edges(
    typename boost::graph_traits< CGAL::Delaunay_triangulation_2<Gt,Tds> >::vertex_descriptor u, 
    const CGAL::Delaunay_triangulation_2<Gt,Tds>& g)
  {
    typename CGAL::Delaunay_triangulation_2<Gt,Tds>::Edge_circulator ec(u,u->face());
    typename boost::graph_traits< CGAL::Delaunay_triangulation_2<Gt,Tds> >::degree_size_type out_deg = out_degree(u,g);
    typedef typename boost::graph_traits< CGAL::Delaunay_triangulation_2<Gt,Tds> >
      ::out_edge_iterator Iter;
    
    return std::make_pair( Iter(ec), Iter(ec,out_deg) );
  }

  template <class Gt, class Tds>
  inline std::pair<
    typename boost::graph_traits< CGAL::Delaunay_triangulation_2<Gt,Tds> >::in_edge_iterator,
    typename boost::graph_traits< CGAL::Delaunay_triangulation_2<Gt,Tds> >::in_edge_iterator >  
  in_edges(
    typename boost::graph_traits< CGAL::Delaunay_triangulation_2<Gt,Tds> >::vertex_descriptor u, 
    const CGAL::Delaunay_triangulation_2<Gt,Tds>& g)
  {
    typename CGAL::Delaunay_triangulation_2<Gt,Tds>::Edge_circulator ec(u,u->face());
    typename boost::graph_traits< CGAL::Delaunay_triangulation_2<Gt,Tds> >::degree_size_type out_deg = out_degree(u,g);
    typedef typename boost::graph_traits< CGAL::Delaunay_triangulation_2<Gt,Tds> >
      ::in_edge_iterator Iter;
    return std::make_pair( Iter(ec), Iter(ec,out_deg) );
  }

  template <class Gt, class Tds>
  inline std::pair<
    typename boost::graph_traits< CGAL::Delaunay_triangulation_2<Gt,Tds> >::adjacency_iterator,
    typename boost::graph_traits< CGAL::Delaunay_triangulation_2<Gt,Tds> >::adjacency_iterator >  
  adjacent_vertices(
    typename boost::graph_traits< CGAL::Delaunay_triangulation_2<Gt,Tds> >::vertex_descriptor u, 
    const CGAL::Delaunay_triangulation_2<Gt,Tds>& g)
  {
    typename CGAL::Delaunay_triangulation_2<Gt,Tds>::Vertex_circulator vc = out_edge_iterator(u,u.face());
    typename boost::graph_traits< CGAL::Delaunay_triangulation_2<Gt,Tds> >::degree_size_type out_deg = out_degree(u,g);
    typedef typename boost::graph_traits< CGAL::Delaunay_triangulation_2<Gt,Tds> >
      ::adjacency_iterator Iter;
    return std::make_pair( Iter(vc), Iter(vc,out_deg) );
  }

  template <class Gt, class Tds>
  typename boost::graph_traits< CGAL::Delaunay_triangulation_2<Gt,Tds> >::vertices_size_type
  num_vertices(const CGAL::Delaunay_triangulation_2<Gt,Tds>& g)
  {
    return g.number_of_vertices()+1;
  }  

  template <class Gt, class Tds>
  typename boost::graph_traits< CGAL::Delaunay_triangulation_2<Gt,Tds> >::edges_size_type
  num_edges(const CGAL::Delaunay_triangulation_2<Gt,Tds>& g)
  {
    return  g.number_of_vertices() + 1 + g.number_of_faces() + degree(g.infinite_vertex(), g) - 2;
  }  

  template <class Gt, class Tds>
  typename boost::graph_traits< CGAL::Delaunay_triangulation_2<Gt,Tds> >::degree_size_type
  out_degree(
    typename boost::graph_traits< CGAL::Delaunay_triangulation_2<Gt,Tds> >::vertex_descriptor u, 
    const CGAL::Delaunay_triangulation_2<Gt,Tds>& g)
  {
    typename boost::graph_traits< CGAL::Delaunay_triangulation_2<Gt,Tds> >::degree_size_type deg = 0;
    typename CGAL::Delaunay_triangulation_2<Gt,Tds>::Edge_circulator c = g.incident_edges(u), done(c);
    if ( c != 0) {
        do {
            ++deg;
        } while (++c != done);
    }
    return deg;
  }

  template <class Gt, class Tds>
  typename boost::graph_traits< CGAL::Delaunay_triangulation_2<Gt,Tds> >::degree_size_type
  in_degree(
    typename boost::graph_traits< CGAL::Delaunay_triangulation_2<Gt,Tds> >::vertex_descriptor u, 
    const CGAL::Delaunay_triangulation_2<Gt,Tds>& g)
  {
    typename boost::graph_traits< CGAL::Delaunay_triangulation_2<Gt,Tds> >::degree_size_type deg = 0;
    typename CGAL::Delaunay_triangulation_2<Gt,Tds>::Edge_circulator c = g.incident_edges(u), done(c);
    if ( c != 0) {
        do {
            ++deg;
        } while (++c != done);
    }
    return deg;
  }

  template <class Gt, class Tds>
  typename boost::graph_traits< CGAL::Delaunay_triangulation_2<Gt,Tds> >::degree_size_type
  degree(
    typename boost::graph_traits< CGAL::Delaunay_triangulation_2<Gt,Tds> >::vertex_descriptor u, 
    const CGAL::Delaunay_triangulation_2<Gt,Tds>& g)
  {
    typename boost::graph_traits< CGAL::Delaunay_triangulation_2<Gt,Tds> >::degree_size_type deg = 0;
    typename CGAL::Delaunay_triangulation_2<Gt,Tds>::Edge_circulator c = g.incident_edges(u), done(c);
    if ( c != 0) {
        do {
            ++deg;
        } while (++c != done);
    }
    return deg;
  }


  // property maps
  template <class Gt, class Tds>
  class DT2_vertex_id_map
    : public boost::put_get_helper<int, DT2_vertex_id_map<Gt,Tds> >
  {
  public:
    typedef boost::readable_property_map_tag category;
    typedef int value_type;
    typedef int reference;
    typedef typename CGAL::Delaunay_triangulation_2<Gt,Tds>::Vertex_handle key_type;
    
    DT2_vertex_id_map()
    {}
    
    long operator[](key_type vh) const {
      return vh->id(); 
    }
  };

  template <class Gt, class Tds>
  class DT2_edge_id_map
    : public boost::put_get_helper<int, DT2_edge_id_map<Gt,Tds> >
  {
  public:
    typedef boost::readable_property_map_tag category;
    typedef int value_type;
    typedef int reference;
    typedef typename CGAL::Delaunay_triangulation_2<Gt,Tds>::Edge key_type;
    
    DT2_edge_id_map()
    {}
    
    long operator[](key_type e) const {
      return (3 * e.first.id()) + e.second; 
    }
  };

  template <class Gt, class Tds>
  class DT2_edge_weight_map
    : public boost::put_get_helper<typename Gt::FT, DT2_edge_weight_map<Gt, Tds> >
  {
  private:
    const CGAL::Delaunay_triangulation_2<Gt,Tds>& tr;
  public:
    typedef boost::readable_property_map_tag category;
    typedef typename Gt::FT value_type;
    typedef value_type reference;
    typedef typename CGAL::Delaunay_triangulation_2<Gt,Tds>::Edge key_type;

    DT2_edge_weight_map(const CGAL::Delaunay_triangulation_2<Gt,Tds>& tr_) 
      : tr(tr_) 
    { }

    typename Gt::FT operator[](key_type e) const {
      return tr.segment(e).squared_length();
    }
  };


  template <class Gt, class Tds>
  inline DT2_vertex_id_map<Gt,Tds>
  get(boost::vertex_index_t, const CGAL::Delaunay_triangulation_2<Gt,Tds>& ) {
    DT2_vertex_id_map<Gt,Tds> m;
    return m;
  }

  template <class Gt, class Tds>
  inline DT2_edge_id_map<Gt,Tds>
  get(boost::edge_index_t, const CGAL::Delaunay_triangulation_2<Gt,Tds>& ) {
    DT2_edge_id_map<Gt,Tds> m;
    return m;
  }

  template <class Gt, class Tds>
  inline DT2_edge_weight_map<Gt,Tds>
  get(boost::edge_weight_t, const CGAL::Delaunay_triangulation_2<Gt,Tds>& g) {
    DT2_edge_weight_map<Gt,Tds> m(g);
    return m;
  }

  template <class Tag>
  struct DT2_property_map { };

  template <>
  struct DT2_property_map<boost::vertex_index_t> {
    template <class Gt, class Tds>
    struct bind_ {
      typedef DT2_vertex_id_map<Gt,Tds> type;
      typedef DT2_vertex_id_map<Gt,Tds> const_type;
    };
  };


  template <>
  struct DT2_property_map<boost::edge_index_t> {
    template <class Gt, class Tds>
    struct bind_ {
      typedef DT2_edge_id_map<Gt,Tds> type;
      typedef DT2_edge_id_map<Gt,Tds> const_type;
    };
  };


  template <>
  struct DT2_property_map<boost::edge_weight_t> {
    template <class Gt, class Tds>
    struct bind_ {
      typedef DT2_edge_weight_map<Gt,Tds> type;
      typedef DT2_edge_weight_map<Gt,Tds> const_type;
    };
  };

} // namespace CGAL

namespace boost {

  // g++ 'enumeral_type' in template unification not implemented workaround
  template <class Gt, class Tds, class Tag>
  struct property_map<CGAL::Delaunay_triangulation_2<Gt,Tds>, Tag> {
    typedef typename 
    CGAL::DT2_property_map<Tag>::template bind_<Gt,Tds> map_gen;
    typedef typename map_gen::type type;
    typedef typename map_gen::const_type const_type;
  };

  // see struct property_map in Polyehdron for an explanation
  template <class Gt, class Tds, class Tag>
  struct property_map<const CGAL::Delaunay_triangulation_2<Gt,Tds>, Tag> {
    typedef typename 
    CGAL::DT2_property_map<Tag>::template bind_<Gt,Tds> map_gen;
    typedef typename map_gen::type type;
    typedef typename map_gen::const_type const_type;
  };

} // namespace boost

namespace CGAL {
  template <class Gt, class Tds, class PropertyTag, class Key>
  inline
  typename boost::property_traits<
    typename boost::property_map<CGAL::Delaunay_triangulation_2<Gt,Tds>,PropertyTag>::const_type>::value_type
  get(PropertyTag p, const CGAL::Delaunay_triangulation_2<Gt,Tds>& g, const Key& key) {
    return get(get(p, g), key);
  }
  
  template <class Gt, class Tds, class PropertyTag, class Key,class Value>
  inline void
  put(PropertyTag p, CGAL::Delaunay_triangulation_2<Gt,Tds>& g, 
      const Key& key, const Value& value)
  {
    typedef typename boost::property_map<CGAL::Delaunay_triangulation_2<Gt,Tds>, PropertyTag>::type Map;
    Map pmap = get(p, g);
    put(pmap, key, value);
  }

} // namespace CGAL

namespace boost {

  // What are those needed for ???
  template <typename Gt, typename Tds>
  struct edge_property_type<CGAL::Delaunay_triangulation_2<Gt,Tds> > {
    typedef void type;
  };  

  template <typename Gt, typename Tds>
  struct vertex_property_type<CGAL::Delaunay_triangulation_2<Gt,Tds> > {
    typedef void type;
  };
} // namespace boost

//#include <CGAL/graph_traits_Delaunay_delaunay_triangulation_2.h>

#endif // CGAL_GRAPH_TRAITS_DELAUNAY_TRIANGULATION_2_H
