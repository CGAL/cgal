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

#ifndef CGAL_GRAPH_TRAITS_CONSTRAINED_TRIANGULATION_2_H
#define CGAL_GRAPH_TRAITS_CONSTRAINED_TRIANGULATION_2_H

#include <CGAL/license/Triangulation_2.h>


// include this to avoid a VC15 warning
#include <CGAL/boost/graph/named_function_params.h>

#include <boost/config.hpp>
#include <boost/iterator_adaptors.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>
#include <CGAL/boost/graph/graph_traits_Triangulation_2.h>
#include <CGAL/Constrained_triangulation_2.h>

// The functions and classes in this file allows the user to
// treat a CGAL Constrained_triangulation_2 object as a boost graph "as is". No
// wrapper is needed for the Constrained_triangulation_2 object.



namespace boost { 

  template <class GT, class TDS, class ITAG>
  struct graph_traits< CGAL::Constrained_triangulation_2<GT,TDS,ITAG> > {

    struct DT2_graph_traversal_category : 
      public virtual bidirectional_graph_tag,
      public virtual adjacency_graph_tag,        
      public virtual edge_list_graph_tag,
      public virtual vertex_list_graph_tag { };

    typedef CGAL::Constrained_triangulation_2<GT,TDS,ITAG> Constrained_triangulation;

    typedef typename CGAL::Constrained_triangulation_2<GT,TDS,ITAG>::Vertex_handle vertex_descriptor;
    typedef typename CGAL::Constrained_triangulation_2<GT,TDS,ITAG>::Face_handle face_descriptor;
    typedef CGAL::detail::Edge<CGAL::Constrained_triangulation_2<GT,TDS,ITAG>, typename CGAL::Constrained_triangulation_2<GT,TDS,ITAG>::Edge>  edge_descriptor;
    typedef typename CGAL::Constrained_triangulation_2<GT,TDS,ITAG>::All_edges_iterator  edge_iterator;

    typedef CGAL::detail::T2_halfedge_descriptor<typename Constrained_triangulation::Triangulation> halfedge_descriptor;

    typedef typename Constrained_triangulation::All_halfedges_iterator  halfedge_iterator;
    typedef CGAL::Prevent_deref<typename Constrained_triangulation::All_vertices_iterator> vertex_iterator;
    typedef CGAL::Prevent_deref<typename Constrained_triangulation::All_faces_iterator> face_iterator;
    typedef CGAL::Counting_iterator<CGAL::detail::Out_edge_circulator<typename Constrained_triangulation::Edge_circulator, edge_descriptor>, edge_descriptor > out_edge_iterator;
    typedef CGAL::Counting_iterator<CGAL::detail::In_edge_circulator<typename Constrained_triangulation::Edge_circulator, edge_descriptor>, edge_descriptor > in_edge_iterator;
    typedef CGAL::Counting_iterator<typename Constrained_triangulation::Vertex_circulator> Incident_vertices_iterator;
    typedef Incident_vertices_iterator adjacency_iterator;

    typedef undirected_tag directed_category;
    typedef disallow_parallel_edge_tag edge_parallel_category; 
    typedef DT2_graph_traversal_category traversal_category;
    typedef typename Constrained_triangulation::size_type size_type;
    typedef size_type vertices_size_type;
    typedef size_type edges_size_type;
    typedef size_type degree_size_type;
  };


} // namespace boost

namespace CGAL {



  template <class Gt, class Tds, class ITAG>
  typename boost::graph_traits< CGAL::Constrained_triangulation_2<Gt,Tds,ITAG> >::vertex_descriptor
  source(typename boost::graph_traits< CGAL::Constrained_triangulation_2<Gt,Tds,ITAG> >::edge_descriptor e,
         const CGAL::Constrained_triangulation_2<Gt,Tds,ITAG>& g)
  {
    return e.first->vertex(g.ccw(e.second));
  }

  template <class Gt, class Tds, class ITAG>
  typename boost::graph_traits< CGAL::Constrained_triangulation_2<Gt,Tds,ITAG> >::vertex_descriptor
  target(typename boost::graph_traits< CGAL::Constrained_triangulation_2<Gt,Tds,ITAG> >::edge_descriptor e,
         const CGAL::Constrained_triangulation_2<Gt,Tds,ITAG>& g)
  {
    return e.first->vertex(g.cw(e.second));
  }



  template <class Gt, class Tds, class ITAG>
  inline std::pair<
    typename boost::graph_traits< CGAL::Constrained_triangulation_2<Gt,Tds,ITAG> >::out_edge_iterator,
    typename boost::graph_traits< CGAL::Constrained_triangulation_2<Gt,Tds,ITAG> >::out_edge_iterator >  
  out_edges(
    typename boost::graph_traits< CGAL::Constrained_triangulation_2<Gt,Tds,ITAG> >::vertex_descriptor u, 
    const CGAL::Constrained_triangulation_2<Gt,Tds,ITAG>& g)
  {
    typename CGAL::Constrained_triangulation_2<Gt,Tds,ITAG>::Edge_circulator ec(u,u->face());
    typename boost::graph_traits< CGAL::Constrained_triangulation_2<Gt,Tds,ITAG> >::degree_size_type out_deg = out_degree(u,g);
    typedef typename boost::graph_traits< CGAL::Constrained_triangulation_2<Gt,Tds,ITAG> >
      ::out_edge_iterator Iter;
    
    return std::make_pair( Iter(ec), Iter(ec,out_deg) );
  }

  template <class Gt, class Tds, class ITAG>
  inline std::pair<
    typename boost::graph_traits< CGAL::Constrained_triangulation_2<Gt,Tds,ITAG> >::in_edge_iterator,
    typename boost::graph_traits< CGAL::Constrained_triangulation_2<Gt,Tds,ITAG> >::in_edge_iterator >  
  in_edges(
    typename boost::graph_traits< CGAL::Constrained_triangulation_2<Gt,Tds,ITAG> >::vertex_descriptor u, 
    const CGAL::Constrained_triangulation_2<Gt,Tds,ITAG>& g)
  {
    typename CGAL::Constrained_triangulation_2<Gt,Tds,ITAG>::Edge_circulator ec(u,u->face());
    typename boost::graph_traits< CGAL::Constrained_triangulation_2<Gt,Tds,ITAG> >::degree_size_type out_deg = out_degree(u,g);
    typedef typename boost::graph_traits< CGAL::Constrained_triangulation_2<Gt,Tds,ITAG> >
      ::in_edge_iterator Iter;
    return std::make_pair( Iter(ec), Iter(ec,out_deg) );
  }



  // property maps
  template <class Gt, class Tds, class ITAG>
  class CT2_vertex_id_map
    : public boost::put_get_helper<int, CT2_vertex_id_map<Gt,Tds,ITAG> >
  {
  public:
    typedef boost::readable_property_map_tag category;
    typedef int value_type;
    typedef int reference;
    typedef typename CGAL::Constrained_triangulation_2<Gt,Tds,ITAG>::Vertex_handle key_type;
    
    CT2_vertex_id_map()
    {}
    
    long operator[](key_type vh) const {
      return vh->id(); 
    }
  };

 template <class Gt, class Tds, class ITAG>
  class CT2_vertex_point_map
  {
  public:
    typedef boost::lvalue_property_map_tag category;
    typedef typename CGAL::Constrained_triangulation_2<Gt,Tds,ITAG>::Point value_type;
    typedef value_type& reference;
    typedef typename CGAL::Constrained_triangulation_2<Gt,Tds,ITAG>::Vertex_handle key_type;

    friend reference get(CT2_vertex_point_map<Gt,Tds,ITAG>, key_type vh)
    { 
      return vh->point(); 
    }
    friend void put(CT2_vertex_point_map<Gt,Tds,ITAG>, key_type vh, reference v)
    {
      vh->point()=v; 
    }
    reference operator[](key_type vh) const {
      return vh->point();
    }
  };

  template <class Gt, class Tds, class ITAG>
  class CT2_edge_id_map
    : public boost::put_get_helper<int, CT2_edge_id_map<Gt,Tds,ITAG> >
  {
  public:
    typedef boost::readable_property_map_tag category;
    typedef int value_type;
    typedef int reference;
    typedef typename CGAL::Constrained_triangulation_2<Gt,Tds,ITAG>::Edge key_type;
    
    CT2_edge_id_map()
    {}
    
    long operator[](key_type e) const {
      return (3 * e.first.id()) + e.second; 
    }
  };

  template <class Gt, class Tds, class ITAG>
  class CT2_edge_weight_map
    : public boost::put_get_helper<typename Gt::FT, CT2_edge_weight_map<Gt, Tds,ITAG> >
  {
  private:
    const CGAL::Constrained_triangulation_2<Gt,Tds,ITAG>& tr;
  public:
    typedef boost::readable_property_map_tag category;
    typedef typename Gt::FT value_type;
    typedef value_type reference;
    typedef typename CGAL::Constrained_triangulation_2<Gt,Tds,ITAG>::Edge key_type;

    CT2_edge_weight_map(const CGAL::Constrained_triangulation_2<Gt,Tds,ITAG>& tr_) 
      : tr(tr_) 
    { }

    typename Gt::FT operator[](key_type e) const {
      return approximate_sqrt(tr.segment(e).squared_length());
    }
  };


  template <class Gt, class Tds, class ITAG>
  inline CT2_vertex_id_map<Gt,Tds,ITAG>
  get(boost::vertex_index_t, const CGAL::Constrained_triangulation_2<Gt,Tds,ITAG>& ) {
    CT2_vertex_id_map<Gt,Tds,ITAG> m;
    return m;
  }

  template <class Gt, class Tds, class ITAG>
  inline CT2_vertex_point_map<Gt,Tds,ITAG>
  get(boost::vertex_point_t, const CGAL::Constrained_triangulation_2<Gt,Tds,ITAG>& ) {
    CT2_vertex_point_map<Gt,Tds,ITAG> m;
    return m;
  }
  
  template <class Gt, class Tds, class ITAG>
  inline CT2_edge_id_map<Gt,Tds,ITAG>
  get(boost::edge_index_t, const CGAL::Constrained_triangulation_2<Gt,Tds,ITAG>& ) {
    CT2_edge_id_map<Gt,Tds,ITAG> m;
    return m;
  }

  template <class Gt, class Tds, class ITAG>
  inline CT2_edge_weight_map<Gt,Tds,ITAG>
  get(boost::edge_weight_t, const CGAL::Constrained_triangulation_2<Gt,Tds,ITAG>& g) {
    CT2_edge_weight_map<Gt,Tds,ITAG> m(g);
    return m;
  }

  template <class Tag>
  struct CT2_property_map { };

  template <>
  struct CT2_property_map<boost::vertex_index_t> {
    template <class Gt, class Tds, class ITAG>
    struct bind_ {
      typedef CT2_vertex_id_map<Gt,Tds,ITAG> type;
      typedef CT2_vertex_id_map<Gt,Tds,ITAG> const_type;
    };
  };

  template <>
  struct CT2_property_map<boost::vertex_point_t> {
    template <class Gt, class Tds, class ITAG>
    struct bind_ {
      typedef CT2_vertex_point_map<Gt,Tds,ITAG> type;
      typedef CT2_vertex_point_map<Gt,Tds,ITAG> const_type;
    };
  };



  template <>
  struct CT2_property_map<boost::edge_index_t> {
    template <class Gt, class Tds, class ITAG>
    struct bind_ {
      typedef CT2_edge_id_map<Gt,Tds,ITAG> type;
      typedef CT2_edge_id_map<Gt,Tds,ITAG> const_type;
    };
  };


  template <>
  struct CT2_property_map<boost::edge_weight_t> {
    template <class Gt, class Tds, class ITAG>
    struct bind_ {
      typedef CT2_edge_weight_map<Gt,Tds,ITAG> type;
      typedef CT2_edge_weight_map<Gt,Tds,ITAG> const_type;
    };
  };

} // namespace CGAL

namespace boost {

  // g++ 'enumeral_type' in template unification not implemented workaround
  template <class Gt, class Tds, class ITAG, class Tag>
  struct property_map<CGAL::Constrained_triangulation_2<Gt,Tds,ITAG>, Tag> {
    typedef typename 
    CGAL::CT2_property_map<Tag>::template bind_<Gt,Tds,ITAG> map_gen;
    typedef typename map_gen::type type;
    typedef typename map_gen::const_type const_type;
  };

  // see struct property_map in Polyehdron for an explanation
  template <class Gt, class Tds, class ITAG, class Tag>
  struct property_map<const CGAL::Constrained_triangulation_2<Gt,Tds,ITAG>, Tag> {
    typedef typename 
    CGAL::CT2_property_map<Tag>::template bind_<Gt,Tds,ITAG> map_gen;
    typedef typename map_gen::type type;
    typedef typename map_gen::const_type const_type;
  };


  // What are those needed for ???
  template <typename Gt, typename Tds, typename ITAG>
  struct edge_property_type<CGAL::Constrained_triangulation_2<Gt,Tds,ITAG> > {
    typedef void type;
  };  

  template <typename Gt, typename Tds, typename ITAG>
  struct vertex_property_type<CGAL::Constrained_triangulation_2<Gt,Tds,ITAG> > {
    typedef void type;
  };
} // namespace boost


#endif // CGAL_GRAPH_TRAITS_CONSTRAINED_TRIANGULATION_2_H
