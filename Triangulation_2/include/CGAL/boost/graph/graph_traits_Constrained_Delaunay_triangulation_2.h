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

#ifndef CGAL_GRAPH_TRAITS_CONSTRAINED_DELAUNAY_TRIANGULATION_2_H
#define CGAL_GRAPH_TRAITS_CONSTRAINED_DELAUNAY_TRIANGULATION_2_H

#include <CGAL/license/Triangulation_2.h>


// include this to avoid a VC15 warning
#include <CGAL/boost/graph/named_function_params.h>

#include <boost/config.hpp>
#include <boost/iterator_adaptors.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>
#include <CGAL/boost/graph/graph_traits_Constrained_triangulation_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>

// The functions and classes in this file allows the user to
// treat a CGAL Constrained_triangulation_2 object as a boost graph "as is". No
// wrapper is needed for the Constrained_triangulation_2 object.



namespace boost { 

  template <class GT, class TDS, class ITAG>
  struct graph_traits< CGAL::Constrained_Delaunay_triangulation_2<GT,TDS,ITAG> > {

    struct DT2_graph_traversal_category : 
      public virtual bidirectional_graph_tag,
      public virtual adjacency_graph_tag,        
      public virtual edge_list_graph_tag,
      public virtual vertex_list_graph_tag { };

    typedef CGAL::Constrained_Delaunay_triangulation_2<GT,TDS,ITAG> Constrained_Delaunay_triangulation;

    typedef typename CGAL::Constrained_Delaunay_triangulation_2<GT,TDS,ITAG>::Vertex_handle vertex_descriptor;
    typedef typename CGAL::Constrained_Delaunay_triangulation_2<GT,TDS,ITAG>::Face_handle face_descriptor;
    typedef CGAL::detail::Edge<CGAL::Constrained_Delaunay_triangulation_2<GT,TDS,ITAG>, typename CGAL::Constrained_Delaunay_triangulation_2<GT,TDS,ITAG>::Edge>  edge_descriptor;
    typedef typename CGAL::Constrained_Delaunay_triangulation_2<GT,TDS,ITAG>::All_edges_iterator  edge_iterator;

    typedef CGAL::detail::T2_halfedge_descriptor<typename Constrained_Delaunay_triangulation::Triangulation> halfedge_descriptor;
    typedef typename Constrained_Delaunay_triangulation::All_halfedges_iterator  halfedge_iterator;
    typedef CGAL::Prevent_deref<typename Constrained_Delaunay_triangulation::All_vertices_iterator> vertex_iterator;
    typedef CGAL::Prevent_deref<typename Constrained_Delaunay_triangulation::All_faces_iterator> face_iterator;
    typedef CGAL::Counting_iterator<CGAL::detail::Out_edge_circulator<typename Constrained_Delaunay_triangulation::Edge_circulator, edge_descriptor>, edge_descriptor > out_edge_iterator;
    typedef CGAL::Counting_iterator<CGAL::detail::In_edge_circulator<typename Constrained_Delaunay_triangulation::Edge_circulator, edge_descriptor>, edge_descriptor > in_edge_iterator;
    typedef CGAL::Counting_iterator<typename Constrained_Delaunay_triangulation::Vertex_circulator> Incident_vertices_iterator;
    typedef Incident_vertices_iterator adjacency_iterator;

    typedef undirected_tag directed_category;
    typedef disallow_parallel_edge_tag edge_parallel_category; 
    typedef DT2_graph_traversal_category traversal_category;
    typedef typename Constrained_Delaunay_triangulation::size_type size_type;
    typedef size_type vertices_size_type;
    typedef size_type edges_size_type;
    typedef size_type degree_size_type;
  };


} // namespace boost

namespace CGAL {



  template <class Gt, class Tds, class ITAG>
  typename boost::graph_traits< CGAL::Constrained_Delaunay_triangulation_2<Gt,Tds,ITAG> >::vertex_descriptor
  source(typename boost::graph_traits< CGAL::Constrained_Delaunay_triangulation_2<Gt,Tds,ITAG> >::edge_descriptor e,
         const CGAL::Constrained_Delaunay_triangulation_2<Gt,Tds,ITAG>& g)
  {
    return e.first->vertex(g.ccw(e.second));
  }

  template <class Gt, class Tds, class ITAG>
  typename boost::graph_traits< CGAL::Constrained_Delaunay_triangulation_2<Gt,Tds,ITAG> >::vertex_descriptor
  target(typename boost::graph_traits< CGAL::Constrained_Delaunay_triangulation_2<Gt,Tds,ITAG> >::edge_descriptor e,
         const CGAL::Constrained_Delaunay_triangulation_2<Gt,Tds,ITAG>& g)
  {
    return e.first->vertex(g.cw(e.second));
  }



  template <class Gt, class Tds, class ITAG>
  inline std::pair<
    typename boost::graph_traits< CGAL::Constrained_Delaunay_triangulation_2<Gt,Tds,ITAG> >::out_edge_iterator,
    typename boost::graph_traits< CGAL::Constrained_Delaunay_triangulation_2<Gt,Tds,ITAG> >::out_edge_iterator >  
  out_edges(
    typename boost::graph_traits< CGAL::Constrained_Delaunay_triangulation_2<Gt,Tds,ITAG> >::vertex_descriptor u, 
    const CGAL::Constrained_Delaunay_triangulation_2<Gt,Tds,ITAG>& g)
  {
    typename CGAL::Constrained_Delaunay_triangulation_2<Gt,Tds,ITAG>::Edge_circulator ec(u,u->face());
    typename boost::graph_traits< CGAL::Constrained_Delaunay_triangulation_2<Gt,Tds,ITAG> >::degree_size_type out_deg = out_degree(u,g);
    typedef typename boost::graph_traits< CGAL::Constrained_Delaunay_triangulation_2<Gt,Tds,ITAG> >
      ::out_edge_iterator Iter;
    
    return std::make_pair( Iter(ec), Iter(ec,out_deg) );
  }

  template <class Gt, class Tds, class ITAG>
  inline std::pair<
    typename boost::graph_traits< CGAL::Constrained_Delaunay_triangulation_2<Gt,Tds,ITAG> >::in_edge_iterator,
    typename boost::graph_traits< CGAL::Constrained_Delaunay_triangulation_2<Gt,Tds,ITAG> >::in_edge_iterator >  
  in_edges(
    typename boost::graph_traits< CGAL::Constrained_Delaunay_triangulation_2<Gt,Tds,ITAG> >::vertex_descriptor u, 
    const CGAL::Constrained_Delaunay_triangulation_2<Gt,Tds,ITAG>& g)
  {
    typename CGAL::Constrained_Delaunay_triangulation_2<Gt,Tds,ITAG>::Edge_circulator ec(u,u->face());
    typename boost::graph_traits< CGAL::Constrained_Delaunay_triangulation_2<Gt,Tds,ITAG> >::degree_size_type out_deg = out_degree(u,g);
    typedef typename boost::graph_traits< CGAL::Constrained_Delaunay_triangulation_2<Gt,Tds,ITAG> >
      ::in_edge_iterator Iter;
    return std::make_pair( Iter(ec), Iter(ec,out_deg) );
  }



  template <class Gt, class Tds, class ITAG>
  inline CT2_vertex_id_map<Gt,Tds,ITAG>
  get(boost::vertex_index_t, const CGAL::Constrained_Delaunay_triangulation_2<Gt,Tds,ITAG>& ) {
    CT2_vertex_id_map<Gt,Tds,ITAG> m;
    return m;
  }

  template <class Gt, class Tds, class ITAG>
  inline CT2_vertex_point_map<Gt,Tds,ITAG>
  get(boost::vertex_point_t, const CGAL::Constrained_Delaunay_triangulation_2<Gt,Tds,ITAG>& ) {
    CT2_vertex_point_map<Gt,Tds,ITAG> m;
    return m;
  }
  template <class Gt, class Tds, class ITAG>
  inline CT2_edge_id_map<Gt,Tds,ITAG>
  get(boost::edge_index_t, const CGAL::Constrained_Delaunay_triangulation_2<Gt,Tds,ITAG>& ) {
    CT2_edge_id_map<Gt,Tds,ITAG> m;
    return m;
  }

  template <class Gt, class Tds, class ITAG>
  inline CT2_edge_weight_map<Gt,Tds,ITAG>
  get(boost::edge_weight_t, const CGAL::Constrained_Delaunay_triangulation_2<Gt,Tds,ITAG>& g) {
    CT2_edge_weight_map<Gt,Tds,ITAG> m(g);
    return m;
  }

} // namespace CGAL

namespace boost {

  // g++ 'enumeral_type' in template unification not implemented workaround
  template <class Gt, class Tds, class ITAG, class Tag>
  struct property_map<CGAL::Constrained_Delaunay_triangulation_2<Gt,Tds,ITAG>, Tag> {
    typedef typename 
    CGAL::CT2_property_map<Tag>::template bind_<Gt,Tds,ITAG> map_gen;
    typedef typename map_gen::type type;
    typedef typename map_gen::const_type const_type;
  };

  // see struct property_map in Polyehdron for an explanation
  template <class Gt, class Tds, class ITAG, class Tag>
  struct property_map<const CGAL::Constrained_Delaunay_triangulation_2<Gt,Tds,ITAG>, Tag> {
    typedef typename 
    CGAL::CT2_property_map<Tag>::template bind_<Gt,Tds,ITAG> map_gen;
    typedef typename map_gen::type type;
    typedef typename map_gen::const_type const_type;
  };



  // What are those needed for ???
  template <typename Gt, typename Tds, typename ITAG>
  struct edge_property_type<CGAL::Constrained_Delaunay_triangulation_2<Gt,Tds,ITAG> > {
    typedef void type;
  };  

  template <typename Gt, typename Tds, typename ITAG>
  struct vertex_property_type<CGAL::Constrained_Delaunay_triangulation_2<Gt,Tds,ITAG> > {
    typedef void type;
  };
} // namespace boost


#endif // CGAL_GRAPH_TRAITS_CONSTRAINED_DELAUNAY_TRIANGULATION_2_H
