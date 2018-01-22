// Copyright (c) 2009-2010 INRIA Sophia-Antipolis (France).
// Copyright (c) 2014-2017 GeometryFactory Sarl (France)
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
// Author(s)     : Laurent Rineau, Stéphane Tayeb
//
//******************************************************************************
// File Description :
//
//******************************************************************************

#ifndef CGAL_POLYHEDRAL_MESH_DOMAIN_WITH_FEATURES_3_H
#define CGAL_POLYHEDRAL_MESH_DOMAIN_WITH_FEATURES_3_H

#include <CGAL/license/Mesh_3.h>


#include <CGAL/Mesh_3/config.h>

#include <CGAL/Random.h>
#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/Mesh_domain_with_polyline_features_3.h>
#include <CGAL/Mesh_polyhedron_3.h>

#include <CGAL/Mesh_3/Detect_polylines_in_polyhedra.h>
#include <CGAL/Mesh_3/Polyline_with_context.h>
#include <CGAL/Polygon_mesh_processing/Detect_features_in_polyhedra.h>
#include <CGAL/Mesh_3/properties_Polyhedron_3.h>

#include <CGAL/enum.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/helpers.h>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <CGAL/boost/graph/split_graph_into_polylines.h>

#include <boost/iterator/transform_iterator.hpp>
#include <boost/foreach.hpp>

#include <string>
#include <vector>
#include <fstream>


namespace CGAL {

namespace internal {
namespace Mesh_3 {

template <typename Kernel>
struct Angle_tester
{
  template <typename vertex_descriptor, typename Graph>
  bool operator()(vertex_descriptor& v, const Graph& g) const
  {
    typedef typename boost::graph_traits<Graph>::out_edge_iterator out_edge_iterator;
    if (out_degree(v, g) != 2)
      return true;
    else
    {
      out_edge_iterator out_edge_it, out_edges_end;
      boost::tie(out_edge_it, out_edges_end) = out_edges(v, g);

      vertex_descriptor v1 = target(*out_edge_it++, g);
      vertex_descriptor v2 = target(*out_edge_it++, g);
      CGAL_assertion(out_edge_it == out_edges_end);

      const typename Kernel::Point_3& p = g[v];
      const typename Kernel::Point_3& p1 = g[v1];
      const typename Kernel::Point_3& p2 = g[v2];

      return (CGAL::angle(p1, p, p2) == CGAL::ACUTE);
    }
  }
};

template <typename Polyhedron>
struct Is_featured_edge {
  const Polyhedron* polyhedron;
  typename boost::property_map<Polyhedron, halfedge_is_feature_t>::type hifm;
  Is_featured_edge() 
    : polyhedron(0) 
  {} // required by boost::filtered_graph
  
  Is_featured_edge(const Polyhedron& polyhedron)
    : polyhedron(&polyhedron), hifm(get(halfedge_is_feature,polyhedron))
  {}

  bool operator()(typename boost::graph_traits<Polyhedron>::edge_descriptor e) const {
    return get(hifm, halfedge(e, *polyhedron));
  }
}; // end Is_featured_edge<Polyhedron>

template <typename Polyhedron>
struct Is_border_edge {
  const Polyhedron* polyhedron;
  Is_border_edge() : polyhedron(0) {} // required by boost::filtered_graph
  Is_border_edge(const Polyhedron& polyhedron) : polyhedron(&polyhedron) {}

  bool operator()(typename boost::graph_traits<Polyhedron>::edge_descriptor e) const {
    return is_border(halfedge(e, *polyhedron), *polyhedron) ||
      is_border(opposite(halfedge(e, *polyhedron), *polyhedron), *polyhedron);
  }
}; // end Is_featured_edge<Polyhedron>

template<typename Polyhedral_mesh_domain,
         typename Polyline_with_context,
         typename Graph>
struct Extract_polyline_with_context_visitor
{
  typedef typename Polyhedral_mesh_domain::Polyhedron Polyhedron;
  std::vector<Polyline_with_context>& polylines;
  const Graph& graph;

  Extract_polyline_with_context_visitor
  (const Graph& graph,
   typename std::vector<Polyline_with_context>& polylines)
    : polylines(polylines), graph(graph)
  {}

  void start_new_polyline()
  {
    polylines.push_back(Polyline_with_context());
  }

  void add_node(typename boost::graph_traits<Graph>::vertex_descriptor vd)
  {
    if(polylines.back().polyline_content.empty()) {
      polylines.back().polyline_content.push_back(graph[vd]);
    }
  }

  void add_edge(typename boost::graph_traits<Graph>::edge_descriptor ed)
  {
    typename boost::graph_traits<Graph>::vertex_descriptor
      s = source(ed, graph),
      t = target(ed, graph);
    Polyline_with_context& polyline = polylines.back();
    CGAL_assertion(!polyline.polyline_content.empty());
    if(polyline.polyline_content.back() != graph[s]) {
      polyline.polyline_content.push_back(graph[s]);
    } else if(polyline.polyline_content.back() != graph[t]) {
      // if the edge is zero-length, it is ignored
      polyline.polyline_content.push_back(graph[t]);
    }
    const typename boost::edge_bundle_type<Graph>::type &
      set_of_indices = graph[ed];
    polyline.context.adjacent_patches_ids.insert(set_of_indices.begin(),
                                                 set_of_indices.end());
  }

  void end_polyline()
  {
    // ignore degenerated polylines
    if(polylines.back().polyline_content.size() < 2)
      polylines.resize(polylines.size() - 1);
    // else {
    //   std::cerr << "Polyline with " << polylines.back().polyline_content.size()
    //             << " vertices, incident to "
    //             << polylines.back().context.adjacent_patches_ids.size()
    //             << " patches:\n ";
    //   for(auto p: polylines.back().polyline_content)
    //     std::cerr << " " << p;
    //   std::cerr << "\n";
    // }
  }
};


} // end CGAL::internal::Mesh_3
} // end CGAL::internal

/**
 * @class Polyhedral_mesh_domain_with_features_3
 *
 *
 */
template < class IGT_,
           class Polyhedron_ = typename Mesh_polyhedron_3<IGT_>::type,
           class TriangleAccessor=Triangle_accessor_3<Polyhedron_,IGT_>,
           class Patch_id=int,
           class Use_exact_intersection_construction_tag = Tag_true >
class Polyhedral_mesh_domain_with_features_3
  : public Mesh_domain_with_polyline_features_3<
      Polyhedral_mesh_domain_3< Polyhedron_,
                                IGT_,
                                TriangleAccessor,
                                Patch_id,
                                Use_exact_intersection_construction_tag > >
{
  typedef Mesh_domain_with_polyline_features_3<
    Polyhedral_mesh_domain_3<
      Polyhedron_, IGT_, TriangleAccessor,
      Patch_id, Use_exact_intersection_construction_tag > > Base;

  typedef boost::adjacency_list<
    boost::setS, // this avoids parallel edges
    boost::vecS,
    boost::undirectedS,
    typename IGT_::Point_3,
    std::set<typename Base::Surface_patch_index> > Featured_edges_copy_graph;
public:
  typedef Polyhedron_ Polyhedron;

  // Index types
  typedef typename Base::Index                Index;
  typedef typename Base::Corner_index         Corner_index;
  typedef typename Base::Curve_segment_index  Curve_segment_index;
  typedef typename Base::Surface_patch_index  Surface_patch_index;
  typedef typename Base::Subdomain_index      Subdomain_index;

  typedef typename CGAL::Mesh_3::details::Surface_patch_index_generator
  <Subdomain_index,
   Polyhedron,
   Patch_id>::Patch_id P_id;

  // Backward compatibility
#ifndef CGAL_MESH_3_NO_DEPRECATED_SURFACE_INDEX
  typedef Surface_patch_index                 Surface_index;
#endif // CGAL_MESH_3_NO_DEPRECATED_SURFACE_INDEX

  typedef typename Base::R         R;
  typedef typename Base::Point_3   Point_3;
  typedef typename Base::FT        FT;

  typedef CGAL::Tag_true           Has_features;

  typedef std::vector<Point_3> Bare_polyline;
  typedef Mesh_3::Polyline_with_context<Surface_patch_index, Curve_segment_index,
                                        Bare_polyline > Polyline_with_context;
  /// Constructors
  Polyhedral_mesh_domain_with_features_3(const Polyhedron& p,
                                         CGAL::Random* p_rng = NULL)
    : Base(p_rng) , borders_detected_(false)
  {
    stored_polyhedra.resize(1);
    stored_polyhedra[0] = p;
    this->add_primitives(stored_polyhedra[0]);
    this->build();
  }

  CGAL_DEPRECATED
  Polyhedral_mesh_domain_with_features_3(const std::string& filename,
                                         CGAL::Random* p_rng = NULL)
    : Base(p_rng) , borders_detected_(false)
  {
    load_from_file(filename.c_str());
  }

  // The following is needed, because otherwise, when a "const char*" is
  // passed, the constructors templates are a better match, than the
  // constructor with `std::string`.
  CGAL_DEPRECATED
  Polyhedral_mesh_domain_with_features_3(const char* filename,
                                         CGAL::Random* p_rng = NULL)
    : Base(p_rng) , borders_detected_(false)
  {
    load_from_file(filename);
  }

  Polyhedral_mesh_domain_with_features_3(const Polyhedron& p,
                                         const Polyhedron& bounding_p,
                                         CGAL::Random* p_rng = NULL)
    : Base(p_rng) , borders_detected_(false)
  {
    stored_polyhedra.resize(2);
    stored_polyhedra[0] = p;
    stored_polyhedra[1] = bounding_p;
    this->add_primitives(stored_polyhedra[0]);
    this->add_primitives(stored_polyhedra[1]);
    if(CGAL::is_empty(bounding_p)) {
      this->set_surface_only();
    } else {
      this->add_primitives_to_bounding_tree(stored_polyhedra[1]);
    }
  }

  template <typename InputPolyhedraPtrIterator>
  Polyhedral_mesh_domain_with_features_3(InputPolyhedraPtrIterator begin,
                                         InputPolyhedraPtrIterator end,
                                         CGAL::Random* p_rng = NULL)
    : Base(p_rng) , borders_detected_(false)
  {
    stored_polyhedra.reserve(std::distance(begin, end));
    for (; begin != end; ++begin) {
      stored_polyhedra.push_back(**begin);
      this->add_primitives(stored_polyhedra.back());
    }
    this->set_surface_only();
    this->build();
  }

  template <typename InputPolyhedraPtrIterator>
  Polyhedral_mesh_domain_with_features_3(InputPolyhedraPtrIterator begin,
                                         InputPolyhedraPtrIterator end,
                                         const Polyhedron& bounding_polyhedron,
                                         CGAL::Random* p_rng = NULL)
    : Base(p_rng) , borders_detected_(false)
  {
    stored_polyhedra.reserve(std::distance(begin, end)+1);
    if(begin != end) {
      for (; begin != end; ++begin) {
        stored_polyhedra.push_back(**begin);
        this->add_primitives(stored_polyhedra.back());
      }
      stored_polyhedra.push_back(bounding_polyhedron);
      this->add_primitives(stored_polyhedra.back());
    }
    if(bounding_polyhedron.empty()) {
      this->set_surface_only();
    } else {
      this->add_primitives_to_bounding_tree(stored_polyhedra.back());
    }
    this->build();
  }

  /// Destructor
  ~Polyhedral_mesh_domain_with_features_3() {}

  /// Detect features
  void initialize_ts(Polyhedron& p);

  void detect_features(FT angle_in_degree, std::vector<Polyhedron>& p);
  void detect_features(FT angle_in_degree = FT(60))
  {
    detect_features(angle_in_degree, stored_polyhedra);
  }

  void detect_borders(std::vector<Polyhedron>& p);
  void detect_borders() { detect_borders(stored_polyhedra); };

private:
  void load_from_file(const char* filename) {
    // Create input polyhedron
    std::ifstream input(filename);
    stored_polyhedra.resize(1);
    input >> stored_polyhedra[0];
    this->add_primitives(stored_polyhedra[0]);
    this->build();
  }

  void add_features_from_split_graph_into_polylines(Featured_edges_copy_graph& graph);

  template <typename Edge_predicate, typename P2Vmap>
  void add_featured_edges_to_graph(const Polyhedron& poly,
                                   const Edge_predicate& pred,
                                   Featured_edges_copy_graph& graph,
                                   P2Vmap& p2vmap);

  std::vector<Polyhedron> stored_polyhedra;
  bool borders_detected_;

private:
  // Disabled copy constructor & assignment operator
  typedef Polyhedral_mesh_domain_with_features_3 Self;
  Polyhedral_mesh_domain_with_features_3(const Self& src);
  Self& operator=(const Self& src);

};  // end class Polyhedral_mesh_domain_with_features_3


template < typename GT_, typename P_, typename TA_,
           typename Tag_, typename E_tag_>
void
Polyhedral_mesh_domain_with_features_3<GT_,P_,TA_,Tag_,E_tag_>::
initialize_ts(Polyhedron& p)
{
  typedef typename boost::property_map<Polyhedron,vertex_time_stamp_t>::type Vtmap;
  typedef typename boost::property_map<Polyhedron,halfedge_time_stamp_t>::type Htmap;
  typedef typename boost::property_map<Polyhedron,face_time_stamp_t>::type Ftmap;
  Vtmap vtm = get(vertex_time_stamp,p);
  Htmap htm = get(halfedge_time_stamp,p);
  Ftmap ftm = get(face_time_stamp,p);

  std::size_t ts = 0;
  BOOST_FOREACH(typename boost::graph_traits<Polyhedron>::vertex_descriptor vd, vertices(p))
  {
    put(vtm,vd,ts++);
  }

  BOOST_FOREACH(typename boost::graph_traits<Polyhedron>::face_descriptor fd, faces(p))
  {
    put(ftm,fd,ts++);
  }

  BOOST_FOREACH(typename boost::graph_traits<Polyhedron>::halfedge_descriptor hd, halfedges(p))
  {
    put(htm,hd,ts++);
  }
}




template <typename Graph>
void dump_graph_edges(std::ostream& out, const Graph& g)
{
  typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<Graph>::edge_descriptor edge_descriptor;

  out.precision(17);
  BOOST_FOREACH(edge_descriptor e, edges(g))
  {
    vertex_descriptor s = source(e, g);
    vertex_descriptor t = target(e, g);
    out << "2 " << g[s] << " " << g[t] << "\n";
  }
}

template <typename Graph>
void dump_graph_edges(const char* filename, const Graph& g)
{
  std::ofstream out(filename);
  dump_graph_edges(out, g);
}

template < typename GT_, typename P_, typename TA_,
           typename Tag_, typename E_tag_>
void
Polyhedral_mesh_domain_with_features_3<GT_,P_,TA_,Tag_,E_tag_>::
detect_features(FT angle_in_degree, std::vector<Polyhedron>& poly)
{
  CGAL_assertion(!borders_detected_);
  if (borders_detected_)
    return;//prevent from not-terminating

  typedef Featured_edges_copy_graph G_copy;
  G_copy g_copy;
  typedef typename boost::graph_traits<G_copy>::vertex_descriptor vertex_descriptor;
  typedef typename boost::property_traits<typename boost::property_map<Polyhedron,vertex_point_t>::type>::value_type Point;
  typedef std::map<Point,
                   vertex_descriptor> P2vmap;
  // TODO: replace this map by and unordered_map
  P2vmap p2vmap;

  CGAL::Polygon_mesh_processing::Detect_features_in_polyhedra<Polyhedron,Surface_patch_index> detect_features;
  BOOST_FOREACH(Polyhedron& p, poly)
  {
    initialize_ts(p);

    // Get sharp features
    detect_features.detect_sharp_edges(p, angle_in_degree);
    detect_features.detect_surface_patches(p);
    detect_features.detect_vertices_incident_patches(p);

    internal::Mesh_3::Is_featured_edge<Polyhedron> is_featured_edge(p);

    add_featured_edges_to_graph(p, is_featured_edge, g_copy, p2vmap);
  }
  add_features_from_split_graph_into_polylines(g_copy);
  borders_detected_ = true;/*done by Mesh_3::detect_features*/
}

template < typename GT_, typename P_, typename TA_,
           typename Tag_, typename E_tag_>
void
Polyhedral_mesh_domain_with_features_3<GT_,P_,TA_,Tag_,E_tag_>::
detect_borders(std::vector<Polyhedron>& poly)
{
  if (borders_detected_)
    return;//border detection has already been done

  detect_features(180, poly);

  borders_detected_ = true;
}

template < typename GT_, typename P_, typename TA_,
           typename Tag_, typename E_tag_>
void
Polyhedral_mesh_domain_with_features_3<GT_,P_,TA_,Tag_,E_tag_>::
add_features_from_split_graph_into_polylines(Featured_edges_copy_graph& g_copy)
{
  std::vector<Polyline_with_context> polylines;

  internal::Mesh_3::Extract_polyline_with_context_visitor<
    Polyhedral_mesh_domain_with_features_3,
    Polyline_with_context,
    Featured_edges_copy_graph
    > visitor(g_copy, polylines);
  internal::Mesh_3::Angle_tester<GT_> angle_tester;
  split_graph_into_polylines(g_copy, visitor, angle_tester);

  this->add_features_with_context(polylines.begin(),
                                  polylines.end());

#if CGAL_MESH_3_PROTECTION_DEBUG > 1
  {//DEBUG
    std::ofstream og("polylines_graph.polylines.txt");
    og.precision(17);
    BOOST_FOREACH(const Polyline_with_context& poly, polylines)
    {
      og << poly.polyline_content.size() << " ";
      BOOST_FOREACH(const Point_3& p, poly.polyline_content)
        og << p << " ";
      og << std::endl;
    }
    og.close();
  }
#endif // CGAL_MESH_3_PROTECTION_DEBUG > 1

}

template < typename GT_, typename P_, typename TA_,
           typename Tag_, typename E_tag_>
template <typename Edge_predicate, typename P2vmap>
void
Polyhedral_mesh_domain_with_features_3<GT_,P_,TA_,Tag_,E_tag_>::
add_featured_edges_to_graph(const Polyhedron& p,
                            const Edge_predicate& pred,
                            Featured_edges_copy_graph& g_copy,
                            P2vmap& p2vmap)
{
  typedef boost::filtered_graph<Polyhedron,
                                Edge_predicate > Featured_edges_graph;
  Featured_edges_graph orig_graph(p, pred);

  typedef Featured_edges_graph Graph;
  typedef typename boost::graph_traits<Graph>::vertex_descriptor Graph_vertex_descriptor;
  typedef typename boost::graph_traits<Graph>::edge_descriptor Graph_edge_descriptor;
  typedef Featured_edges_copy_graph G_copy;
  typedef typename boost::graph_traits<G_copy>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<G_copy>::edge_descriptor edge_descriptor;

  const Featured_edges_graph& graph = orig_graph;

  typedef typename boost::property_map<Polyhedron,vertex_point_t>::const_type Vpm;
  Vpm vpm = get(vertex_point, p);
  BOOST_FOREACH(Graph_vertex_descriptor v, vertices(graph)){
    vertex_descriptor vc;
    typename P2vmap::iterator it = p2vmap.find(get(vpm,v));
    if(it == p2vmap.end()) {
      vc = add_vertex(g_copy);
      g_copy[vc] = get(vpm, v);
      p2vmap[get(vpm,v)] = vc;
    }
  }

  typedef typename boost::property_map<Polyhedron,face_patch_id_t<P_id> >::type Face_patch_id_pmap;
  Face_patch_id_pmap fpm = get(face_patch_id_t<P_id>(),p);

  BOOST_FOREACH(Graph_edge_descriptor e, edges(graph)){
    vertex_descriptor vs = p2vmap[get(vpm,source(e,graph))];
    vertex_descriptor vt = p2vmap[get(vpm,target(e,graph))];
    CGAL_warning_msg(vs != vt, "ignore self loop");
    if(vs != vt) {
      const std::pair<edge_descriptor, bool> pair = add_edge(vs,vt,g_copy);
      typename boost::graph_traits<Polyhedron>::halfedge_descriptor he = halfedge(e, p);
      if(!is_border(he, p)) {
        g_copy[pair.first].insert(get(fpm, face(he,p)));
      }
      he = opposite(he,p);
      if(!is_border(he, p)) {
        g_copy[pair.first].insert(get(fpm, face(he,p)));
      }
    }
  }

#if CGAL_MESH_3_PROTECTION_DEBUG > 1
  {// DEBUG
    dump_graph_edges("edges-graph.polylines.txt", g_copy);
  }
#endif
}


} //namespace CGAL


#endif // CGAL_POLYHEDRAL_MESH_DOMAIN_WITH_FEATURES_3_H
