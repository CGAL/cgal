// Copyright (c) 2009-2010 INRIA Sophia-Antipolis (France).
// Copyright (c) 2014-2017 GeometryFactory Sarl (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
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

#include <CGAL/disable_warnings.h>

#include <CGAL/Mesh_3/config.h>

#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/Mesh_domain_with_polyline_features_3.h>
#include <CGAL/Mesh_polyhedron_3.h>
#include <CGAL/Mesh_3/Polyline_with_context.h>

#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/internal/Mesh_3/helpers.h>

#include <CGAL/boost/graph/split_graph_into_polylines.h>
#include <CGAL/boost/iterator/transform_iterator.hpp>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/enum.h>
#include <CGAL/Default.h>
#include <CGAL/Polygon_mesh_processing/detect_features.h>
#include <CGAL/Random.h>

#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/dynamic_bitset.hpp>

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <set>
#include <vector>
#include <utility>

namespace CGAL {

/**
 * @class Polyhedral_mesh_domain_with_features_3
 *
 *
 */
template < class IGT_,
           class Polyhedron_ = typename Mesh_polyhedron_3<IGT_>::type,
           class TriangleAccessor= CGAL::Default,
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
  typedef Polyhedron Polyhedron_type;

  // Index types
  typedef typename Base::Index                Index;
  typedef typename Base::Corner_index         Corner_index;
  typedef typename Base::Curve_index          Curve_index;
  typedef typename Base::Surface_patch_index  Surface_patch_index;
  typedef typename Base::Subdomain_index      Subdomain_index;

#ifndef CGAL_NO_DEPRECATED_CODE
  typedef Curve_index Curve_segment_index; ///< Backward-compatibility
#endif

  typedef typename boost::property_map<Polyhedron,
                                       face_patch_id_t<Patch_id>
                                       >::type            Face_patch_id_pmap;
  typedef typename boost::property_traits<
    Face_patch_id_pmap>::value_type                       P_id;

  // Backward compatibility
#ifndef CGAL_MESH_3_NO_DEPRECATED_SURFACE_INDEX
  typedef Surface_patch_index                 Surface_index;
#endif // CGAL_MESH_3_NO_DEPRECATED_SURFACE_INDEX

  typedef typename Base::R         R;
  typedef typename Base::Point_3   Point_3;
  typedef typename Base::FT        FT;

  typedef CGAL::Tag_true           Has_features;

  typedef std::vector<Point_3> Bare_polyline;
  typedef Mesh_3::Polyline_with_context<Surface_patch_index, Curve_index,
                                        Bare_polyline > Polyline_with_context;
  /// Constructors
  Polyhedral_mesh_domain_with_features_3(const Polyhedron& p,
                                         CGAL::Random* p_rng = nullptr)
    : Base(p_rng) , borders_detected_(false)
  {
    stored_polyhedra.resize(1);
    stored_polyhedra[0] = p;
    get(face_patch_id_t<Patch_id>(), stored_polyhedra[0]);
    this->add_primitives(stored_polyhedra[0]);
    this->build();
  }

#ifndef CGAL_NO_DEPRECATED_CODE

  CGAL_DEPRECATED
  Polyhedral_mesh_domain_with_features_3(const std::string& filename,
                                         CGAL::Random* p_rng = nullptr)
    : Base(p_rng) , borders_detected_(false)
  {
    load_from_file(filename.c_str());
  }

  // The following is needed, because otherwise, when a "const char*" is
  // passed, the constructors templates are a better match, than the
  // constructor with `std::string`.
  CGAL_DEPRECATED
  Polyhedral_mesh_domain_with_features_3(const char* filename,
                                         CGAL::Random* p_rng = nullptr)
    : Base(p_rng) , borders_detected_(false)
  {
    load_from_file(filename);
  }
#endif // not CGAL_NO_DEPRECATED_CODE

  Polyhedral_mesh_domain_with_features_3(const Polyhedron& p,
                                         const Polyhedron& bounding_p,
                                         CGAL::Random* p_rng = nullptr)
    : Base(p_rng) , borders_detected_(false)
  {
    stored_polyhedra.resize(2);
    stored_polyhedra[0] = p;
    stored_polyhedra[1] = bounding_p;
    get(face_patch_id_t<Patch_id>(), stored_polyhedra[0]);
    get(face_patch_id_t<Patch_id>(), stored_polyhedra[1]);
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
                                         CGAL::Random* p_rng = nullptr)
    : Base(p_rng) , borders_detected_(false)
  {
    stored_polyhedra.reserve(std::distance(begin, end));
    for (; begin != end; ++begin) {
      stored_polyhedra.push_back(**begin);
      get(face_patch_id_t<Patch_id>(), stored_polyhedra.back());
      this->add_primitives(stored_polyhedra.back());
    }
    this->set_surface_only();
    this->build();
  }

  template <typename InputPolyhedraPtrIterator>
  Polyhedral_mesh_domain_with_features_3(InputPolyhedraPtrIterator begin,
                                         InputPolyhedraPtrIterator end,
                                         const Polyhedron& bounding_polyhedron,
                                         CGAL::Random* p_rng = nullptr)
    : Base(p_rng) , borders_detected_(false)
  {
    stored_polyhedra.reserve(std::distance(begin, end)+1);
    if(begin != end) {
      for (; begin != end; ++begin) {
        stored_polyhedra.push_back(**begin);
        get(face_patch_id_t<Patch_id>(), stored_polyhedra.back());
        this->add_primitives(stored_polyhedra.back());
      }
    }
    stored_polyhedra.push_back(bounding_polyhedron);
    get(face_patch_id_t<Patch_id>(), stored_polyhedra.back());
    this->add_primitives(stored_polyhedra.back());
    if(CGAL::is_empty(bounding_polyhedron)) {
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

  template <typename InputIterator>
  void
  add_features(InputIterator first, InputIterator end)
  {
    auto max = 0;
    auto min = (std::numeric_limits<int>::max)();
    for(const auto& polyhedron: stored_polyhedra) {
      auto f_pid = get(face_patch_id_t<Patch_id>(), polyhedron);
      for(auto fd : faces(polyhedron)) {
        const auto patch_id = get(f_pid, fd);
        min = (std::min)(patch_id, min);
        max = (std::max)(patch_id, max);
      }
    }
    boost::dynamic_bitset<> patch_ids_bitset;
    patch_ids_bitset.resize(max - min + 1);
    for(const auto& polyhedron: stored_polyhedra) {
      auto f_pid = get(face_patch_id_t<Patch_id>(), polyhedron);
      for(auto fd : faces(polyhedron)) {
        const auto patch_id = get(f_pid, fd);
        patch_ids_bitset.set(patch_id - min);
      }
    }
    using Patch_ids_container = std::vector<int>;
    Patch_ids_container all_patch_ids;
    all_patch_ids.reserve(patch_ids_bitset.count());
    for(auto i = patch_ids_bitset.find_first();
        i != patch_ids_bitset.npos;
        i = patch_ids_bitset.find_next(i))
    {
      all_patch_ids.push_back(static_cast<int>(i + min));
    }
    using Polyline = typename std::iterator_traits<InputIterator>::value_type;
    auto identity_property_map = boost::typed_identity_property_map<Polyline>();
    auto all_patch_ids_pmap =
      boost::static_property_map<Patch_ids_container>(all_patch_ids);
    Base::add_features_and_incidences(first, end,
                                      identity_property_map,
                                      all_patch_ids_pmap);
  }

  // non-documented, provided to the FEniCS project
  const std::vector<Polyhedron>& polyhedra()const
  {
    return stored_polyhedra;
  }

private:
  void load_from_file(const char* filename) {
    // Create input polyhedron
    std::ifstream input(filename);
    stored_polyhedra.resize(1);
    input >> stored_polyhedra[0];
    get(face_patch_id_t<Patch_id>(), stored_polyhedra[0]);
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
  for(typename boost::graph_traits<Polyhedron>::vertex_descriptor vd : vertices(p))
  {
    put(vtm,vd,ts++);
  }

  for(typename boost::graph_traits<Polyhedron>::face_descriptor fd : faces(p))
  {
    put(ftm,fd,ts++);
  }

  for(typename boost::graph_traits<Polyhedron>::halfedge_descriptor hd : halfedges(p))
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
  for(edge_descriptor e : make_range(edges(g)))
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
  namespace PMP = CGAL::Polygon_mesh_processing;
  std::size_t nb_of_patch_plus_one = 1;
  for(Polyhedron& p : poly)
  {
    initialize_ts(p);
    typedef typename boost::property_map<Polyhedron,CGAL::face_patch_id_t<Tag_> >::type PIDMap;
    typedef typename boost::property_map<Polyhedron,CGAL::vertex_incident_patches_t<P_id> >::type VIPMap;
    typedef typename boost::property_map<Polyhedron, CGAL::edge_is_feature_t>::type EIFMap;

    PIDMap pid_map = get(face_patch_id_t<Tag_>(), p);
    VIPMap vip_map = get(vertex_incident_patches_t<P_id>(), p);
    EIFMap eif_map = get(CGAL::edge_is_feature, p);

    // Get sharp features
    nb_of_patch_plus_one += PMP::sharp_edges_segmentation(p, angle_in_degree
      , eif_map
      , pid_map
      , PMP::parameters::first_index(nb_of_patch_plus_one)
                        .face_index_map(get_initialized_face_index_map(p))
                        .vertex_incident_patches_map(vip_map));

    Mesh_3::internal::Is_featured_edge<Polyhedron> is_featured_edge(p);

    add_featured_edges_to_graph(p, is_featured_edge, g_copy, p2vmap);
  }
  add_features_from_split_graph_into_polylines(g_copy);
  borders_detected_ = true;/*done by PMP::detect_features*/
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

  Mesh_3::internal::Extract_polyline_with_context_visitor<
    Polyline_with_context,
    Featured_edges_copy_graph
    > visitor(g_copy, polylines);
  Mesh_3::internal::Angle_tester<GT_> angle_tester;
  split_graph_into_polylines(g_copy, visitor, angle_tester);

  this->add_features_with_context(polylines.begin(),
                                  polylines.end());

#if CGAL_MESH_3_PROTECTION_DEBUG & 2
  {//DEBUG
    std::ofstream og("polylines_graph.polylines.txt");
    og.precision(17);
    for(const Polyline_with_context& poly : polylines)
    {
      og << poly.polyline_content.size() << " ";
      for(const Point_3& p : poly.polyline_content)
        og << p << " ";
      og << std::endl;
    }
    og.close();
  }
#endif // CGAL_MESH_3_PROTECTION_DEBUG & 2

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
  for(Graph_vertex_descriptor v : make_range(vertices(graph))){
    vertex_descriptor vc;
    typename P2vmap::iterator it = p2vmap.find(get(vpm,v));
    if(it == p2vmap.end()) {
      vc = add_vertex(g_copy);
      g_copy[vc] = get(vpm, v);
      p2vmap[get(vpm,v)] = vc;
    }
  }

  typedef typename boost::property_map<Polyhedron,face_patch_id_t<Tag_> >::type Face_patch_id_pmap;
  Face_patch_id_pmap fpm = get(face_patch_id_t<Tag_>(),p);

  for(Graph_edge_descriptor e : make_range(edges(graph))){
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

#if CGAL_MESH_3_PROTECTION_DEBUG & 2
  {// DEBUG
    dump_graph_edges("edges-graph.polylines.txt", g_copy);
  }
#endif
}


} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_POLYHEDRAL_MESH_DOMAIN_WITH_FEATURES_3_H
