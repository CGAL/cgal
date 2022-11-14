// Copyright (c) 2016 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sebastien Loriot

#ifndef CGAL_POLYGON_MESH_PROCESSING_INTERNAL_COREFINEMENT_VISITOR_H
#define CGAL_POLYGON_MESH_PROCESSING_INTERNAL_COREFINEMENT_VISITOR_H

#include <CGAL/license/Polygon_mesh_processing/corefinement.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Polygon_mesh_processing/internal/Corefinement/predicates.h>
#include <CGAL/Polygon_mesh_processing/internal/Corefinement/face_graph_utils.h>
#include <CGAL/utility.h>
#include <CGAL/Default.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Projection_traits_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

#include <boost/container/flat_map.hpp>
#include <boost/container/small_vector.hpp>

namespace CGAL{
namespace Polygon_mesh_processing {
namespace Corefinement{

// TODO option to ignore internal edges for patches of coplanar faces

//binds two edge constrained pmaps
template <class G, class Ecm1, class Ecm2=Ecm1>
struct Ecm_bind{
  G& g1;
  Ecm1& ecm1;
  G& g2;
  Ecm2& ecm2;

  Ecm_bind(G& g1, G& g2, Ecm1& ecm1, Ecm2& ecm2)
  : g1(g1), ecm1(ecm1), g2(g2), ecm2(ecm2)
  {}

  typedef typename boost::graph_traits<G>::edge_descriptor edge_descriptor;

  void call_put(G& g, edge_descriptor e, bool b) const
  {
    if ( &g==&g1 )
      put(ecm1,e,b);
    else
    {
      CGAL_assertion( &g==&g2 );
      put(ecm2,e,b);
    }
  }

  bool call_get(G& g, edge_descriptor e) const
  {
    if ( &g==&g1 )
      return get(ecm1,e);
    CGAL_assertion( &g==&g2 );
    return get(ecm2,e);
  }
};

template <class G>
struct Ecm_bind<G, No_mark<G>, No_mark<G> >
{
  No_mark<G> ecm1, ecm2;
  Ecm_bind(G&, G&, const No_mark<G>&, const No_mark<G>&){}
  typedef typename boost::graph_traits<G>::edge_descriptor edge_descriptor;
  void call_put(G&, edge_descriptor, bool) const {}
  bool call_get(G&, edge_descriptor) const {
    return false;
  }
};

template<class G>
struct No_extra_output_from_corefinement
{
  void start_new_polyline(std::size_t, std::size_t) {}
  void add_node_to_polyline(std::size_t){}
  template<class Node_id_pair, class halfedge_descriptor>
  void set_edge_per_polyline(G& /*tm*/,
                             Node_id_pair /*indices*/,
                             halfedge_descriptor /*hedge*/){}
  template <class vertex_descriptor, class Node_id>
  void set_vertex_id(vertex_descriptor, Node_id, const G&){}
  template <class Node_vector,
            class Mesh_to_map_node>
  void operator()(
    const Node_vector& /*nodes*/,
    bool /*input_have_coplanar_faces*/,
    const boost::dynamic_bitset<>& /* is_node_of_degree_one */,
    const Mesh_to_map_node& /*mesh_to_node_id_to_vertex*/) const
  {}
};

template <class TriangleMesh, bool doing_autorefinement /* = false */>
class Graph_node_classifier
{
  typedef std::size_t Node_id;
  typedef boost::graph_traits<TriangleMesh> Graph_traits;
  typedef typename Graph_traits::vertex_descriptor vertex_descriptor;
  typedef typename Graph_traits::halfedge_descriptor halfedge_descriptor;
  boost::dynamic_bitset<> is_node_on_boundary; // indicate if a vertex is a border vertex in tm1 or tm2
  boost::container::flat_map<TriangleMesh*, std::vector< vertex_descriptor> > m_node_on_vertex_map;
  boost::container::flat_map<TriangleMesh*, std::vector< halfedge_descriptor> > m_node_on_edge_map;
// variables filled by preprocessing
  TriangleMesh* m_tm1_ptr = nullptr;
  const std::vector<vertex_descriptor>* m_node_on_vertex_1_ptr = nullptr;
  const std::vector<halfedge_descriptor>* m_node_on_edge_1_ptr = nullptr;
  TriangleMesh* m_tm2_ptr = nullptr;
  const std::vector<vertex_descriptor>* m_node_on_vertex_2_ptr = nullptr;
  const std::vector<halfedge_descriptor>* m_node_on_edge_2_ptr = nullptr;

  bool is_on_border(std::size_t node_id1, std::size_t node_id2,
                    const std::vector<vertex_descriptor>* node_on_vertex_ptr,
                    const std::vector<halfedge_descriptor>* node_on_edge_ptr,
                    TriangleMesh* tm_ptr)
  {
    if (tm_ptr == nullptr) return false;

    if (node_on_vertex_ptr!=nullptr)
    {
      vertex_descriptor v1 = (*node_on_vertex_ptr)[node_id1];
      if ( v1 != Graph_traits::null_vertex() )
      {
        vertex_descriptor v2 = (*node_on_vertex_ptr)[node_id2];
        if ( v2 != Graph_traits::null_vertex() )
        {
          std::pair< halfedge_descriptor, bool > res =
            halfedge(v1, v2, *tm_ptr);
          CGAL_assertion(res.second);
          return res.second && is_border_edge(res.first, *tm_ptr);
        }
        if (node_on_edge_ptr!=nullptr)
        {
          halfedge_descriptor h = (*node_on_edge_ptr)[node_id2];
          if (h != Graph_traits::null_halfedge() && is_border_edge(h, *tm_ptr))
            return source(h, *tm_ptr)==v1 || target(h, *tm_ptr)==v1;
          return false;
        }
      }
    }

    if (node_on_edge_ptr!=nullptr)
    {
      halfedge_descriptor h = (*node_on_edge_ptr)[node_id1];
      if (h != Graph_traits::null_halfedge() && is_border_edge(h, *tm_ptr))
      {
        if ( node_on_vertex_ptr!=nullptr )
        {
          vertex_descriptor v2 = (*node_on_vertex_ptr)[node_id2];
          if (v2 != Graph_traits::null_vertex() )
            return source(h, *tm_ptr)==v2 || target(h, *tm_ptr)==v2;
        }
        halfedge_descriptor h_bis = (*node_on_edge_ptr)[node_id2];
        if (h_bis != Graph_traits::null_halfedge())
          return h==h_bis || h==opposite(h_bis, *tm_ptr);
      }
    }
    return false;
  }
public:

  void preprocessing()
  {
    boost::container::flat_set<TriangleMesh*> mesh_ptrs;
    mesh_ptrs.reserve(2);
    for (const auto& k_v : m_node_on_vertex_map) mesh_ptrs.insert(k_v.first);
    for (const auto& k_v : m_node_on_edge_map) mesh_ptrs.insert(k_v.first);

    if (!mesh_ptrs.empty())
    {
      m_tm1_ptr=*mesh_ptrs.begin();
      auto itv = m_node_on_vertex_map.find(m_tm1_ptr);
      if (itv != m_node_on_vertex_map.end()) m_node_on_vertex_1_ptr= &(itv->second);
      auto ite = m_node_on_edge_map.find(m_tm1_ptr);
      if (ite != m_node_on_edge_map.end()) m_node_on_edge_1_ptr= &(ite->second);
      if (mesh_ptrs.size()==2)
      {
        m_tm2_ptr=*std::next(mesh_ptrs.begin());
        itv = m_node_on_vertex_map.find(m_tm2_ptr);
        if (itv != m_node_on_vertex_map.end()) m_node_on_vertex_2_ptr= &(itv->second);
        ite = m_node_on_edge_map.find(m_tm2_ptr);
        if (ite != m_node_on_edge_map.end()) m_node_on_edge_2_ptr= &(ite->second);
      }
    }
  }

  void node_on_vertex(Node_id node_id, vertex_descriptor v, const TriangleMesh& tm)
  {
    m_node_on_vertex_map[const_cast<TriangleMesh*>(&tm)][node_id] = v;

    //we turn around the hedge and check no halfedge is a border halfedge
    for(halfedge_descriptor hc :halfedges_around_target(halfedge(v,tm),tm))
      if ( is_border_edge(hc,tm) )
      {
        is_node_on_boundary.set(node_id);
        return;
      }
  }

  void node_on_edge(Node_id node_id, halfedge_descriptor h, const TriangleMesh& tm)
  {
    if ( is_border_edge(h,tm) )
      is_node_on_boundary.set(node_id);
    m_node_on_edge_map[const_cast<TriangleMesh*>(&tm)][node_id] = h;

  }

  void new_node(Node_id node_id, const TriangleMesh& tm)
  {
    is_node_on_boundary.resize(node_id+1, false);
    TriangleMesh* tm_ptr = const_cast<TriangleMesh*>(&tm);
    m_node_on_edge_map[tm_ptr].resize(node_id+1, Graph_traits::null_halfedge());
    m_node_on_vertex_map[tm_ptr].resize(node_id+1, Graph_traits::null_vertex());
  }

  bool is_terminal(Node_id node_id,  const std::vector<Node_id>& neighbor_nodes)
  {
    if ( is_node_on_boundary.test(node_id) && neighbor_nodes.size()==2)
    {
      std::size_t nn1 = neighbor_nodes[0], nn2 = neighbor_nodes[1];

      return  is_on_border(node_id, nn1, m_node_on_vertex_1_ptr, m_node_on_edge_1_ptr, m_tm1_ptr) !=
              is_on_border(node_id, nn2, m_node_on_vertex_1_ptr, m_node_on_edge_1_ptr, m_tm1_ptr)
          ||
              is_on_border(node_id, nn1, m_node_on_vertex_2_ptr, m_node_on_edge_2_ptr, m_tm2_ptr) !=
              is_on_border(node_id, nn2, m_node_on_vertex_2_ptr, m_node_on_edge_2_ptr, m_tm2_ptr);
    }
    return false;
  }
};

template <class TriangleMesh>
class Graph_node_classifier<TriangleMesh, /* doing_autorefinement = */ true >
{
  typedef std::size_t Node_id;
  typedef boost::graph_traits<TriangleMesh> Graph_traits;
  typedef typename Graph_traits::vertex_descriptor vertex_descriptor;
  typedef typename Graph_traits::halfedge_descriptor halfedge_descriptor;
  boost::dynamic_bitset<> m_is_node_on_boundary; // indicate if a vertex is a border vertex in tm1 or tm2
  std::vector<boost::container::small_vector<vertex_descriptor,2>> m_node_on_vertex;
  std::vector<boost::container::small_vector<halfedge_descriptor,2>> m_node_on_edge;
  TriangleMesh* m_tm_ptr = nullptr;

  bool is_on_border(std::size_t node_id1, std::size_t node_id2)
  {
    if (m_tm_ptr == nullptr) return false;

    for (vertex_descriptor v1 : m_node_on_vertex[node_id1])
    {
      for (vertex_descriptor v2 : m_node_on_vertex[node_id2])
      {
        //vertex-vertex case
        std::pair< halfedge_descriptor, bool > res =
          halfedge(v1, v2, *m_tm_ptr);
        if (res.second && is_border_edge(res.first, *m_tm_ptr))
          return true;
      }
      for (halfedge_descriptor h2 : m_node_on_edge[node_id2])
      {
        // vertex-edge case
        if ( (source(h2, *m_tm_ptr)==v1 || target(h2, *m_tm_ptr)==v1)
              && is_border_edge(h2, *m_tm_ptr) )
        {
          return true;
        }
      }
    }

    for (halfedge_descriptor h1 : m_node_on_edge[node_id1])
    {
      if (!is_border_edge(h1, *m_tm_ptr)) continue;
      for (vertex_descriptor v2 : m_node_on_vertex[node_id2])
      {
        // edge-vertex case
        if (source(h1, *m_tm_ptr)==v2 || target(h1, *m_tm_ptr)==v2)
          return true;
      }
      for (halfedge_descriptor h2 : m_node_on_edge[node_id2])
      {
        if (h1==h2 || h1==opposite(h2, *m_tm_ptr))
          return true;
      }
    }

    return false;
  }

public:

  void preprocessing(){}

  void node_on_vertex(Node_id node_id, vertex_descriptor v, const TriangleMesh& tm)
  {
    m_node_on_vertex[node_id].push_back(v);

    //we turn around the hedge and check no halfedge is a border halfedge
    for(halfedge_descriptor hc :halfedges_around_target(halfedge(v,tm),tm))
      if ( is_border_edge(hc,tm) )
      {
        m_is_node_on_boundary.set(node_id);
        return;
      }
  }

  void node_on_edge(Node_id node_id, halfedge_descriptor h, const TriangleMesh& tm)
  {
     if ( is_border_edge(h,tm) )
       m_is_node_on_boundary.set(node_id);
     m_node_on_edge[node_id].push_back(h);

  }

  void new_node(Node_id node_id, const TriangleMesh& tm)
  {
    m_is_node_on_boundary.resize(node_id+1, false);
    m_tm_ptr = const_cast<TriangleMesh*>(&tm);
    m_node_on_edge.resize(node_id+1);
    m_node_on_vertex.resize(node_id+1);
  }

  bool is_terminal(Node_id node_id,  const std::vector<Node_id>& neighbor_nodes)
  {
    if ( m_is_node_on_boundary.test(node_id) && neighbor_nodes.size()==2)
    {
      std::size_t nn1 = neighbor_nodes[0], nn2 = neighbor_nodes[1];

      return  is_on_border(node_id, nn1) != is_on_border(node_id, nn2);
    }
    return false;
  }
};

namespace internal{
//version for corefinement, only one vertex per node_id
template <class TriangleMesh, bool doing_autorefinement=false>
struct Node_id_to_vertex
{
  typedef boost::graph_traits<TriangleMesh>                        Graph_traits;
  typedef typename Graph_traits::vertex_descriptor            vertex_descriptor;
  std::vector<vertex_descriptor> data;

  vertex_descriptor get_vertex(std::size_t i) const
  {
    return data[i];
  }
  void register_vertex(std::size_t i, vertex_descriptor v)
  {
    data[i] = v;
  }
  void set_vertex_for_retriangulation(std::size_t i, vertex_descriptor v)
  {
    data[i] = v;
  }
  void set_temporary_vertex_for_retriangulation(std::size_t i, vertex_descriptor v)
  {
    data[i] = v;
  }
  void resize(std::size_t n)
  {
    data.resize(n,Graph_traits::null_vertex());
  }
  std::size_t size() const
  {
    return data.size();
  }
  template <class VPM, class Point_3>
  void update_vertex_point(std::size_t i, const Point_3& p, const VPM& vpm) const
  {
    if (data[i]!=Graph_traits::null_vertex())
      put(vpm, data[i], p);
  }
};

//version for autorefinement and non-manifold corefinement, several vertices per node_id
template <class TriangleMesh>
struct Node_id_to_vertex<TriangleMesh, true>
{
  typedef boost::graph_traits<TriangleMesh>                        Graph_traits;
  typedef typename Graph_traits::vertex_descriptor            vertex_descriptor;
  std::vector< std::vector<vertex_descriptor> > data;

  vertex_descriptor get_vertex(std::size_t i) const
  {
    if (data[i].empty())
      return Graph_traits::null_vertex();
    return data[i].back();
  }
  void register_vertex(std::size_t i, vertex_descriptor v)
  {
    data[i].push_back(v);
  }
  void set_temporary_vertex_for_retriangulation(std::size_t i, vertex_descriptor v)
  {
    data[i].assign(1,v);
  }
  // warning: data[i] might then contains several times the same vertex
  //          but it is probably still a better option than look for the
  //          vertex and remove it
  void set_vertex_for_retriangulation(std::size_t i, vertex_descriptor v)
  {
    assert(!data[i].empty());
    if (data[i].back()!=v)
      data[i].push_back(v);
  }
  void resize(std::size_t n)
  {
    data.resize(n);
  }
  std::size_t size() const
  {
    return data.size();
  }
  template <class VPM, class Point_3>
  void update_vertex_point(std::size_t i, const Point_3& p, const VPM& vpm) const
  {
    for (vertex_descriptor v : data[i])
      put(vpm, v, p);
  }
};
}


// A visitor for Intersection_of_triangle_meshes that can be used to corefine
// two meshes
template< class TriangleMesh,
          class VertexPointMap1,
          class VertexPointMap2,
          class OutputBuilder_ = Default,
          class EdgeMarkMapBind_ = Default,
          class UserVisitor_ = Default,
          bool doing_autorefinement = false,
          bool handle_non_manifold_features = false >
class Surface_intersection_visitor_for_corefinement{
//default template parameters
  typedef typename Default::Get<EdgeMarkMapBind_,
    Ecm_bind<TriangleMesh, No_mark<TriangleMesh> > >::type      EdgeMarkMapBind;
  typedef typename Default::Get<OutputBuilder_,
    No_extra_output_from_corefinement<TriangleMesh> >::type       OutputBuilder;
  typedef typename Default::Get<
    UserVisitor_, Default_visitor<TriangleMesh> >::type  UserVisitor;

// config flags
public:
  static const bool Predicates_on_constructions_needed = true;
  static const bool do_need_vertex_graph = true;
// typdefs
private:
  typedef std::size_t                                                   Node_id;
  typedef boost::graph_traits<TriangleMesh>                        Graph_traits;
  typedef typename Graph_traits::edge_descriptor                edge_descriptor;
  typedef typename Graph_traits::face_descriptor                face_descriptor;
  typedef typename Graph_traits::vertex_descriptor            vertex_descriptor;
  typedef typename Graph_traits::halfedge_descriptor        halfedge_descriptor;
   typedef std::vector<Node_id>                                        Node_ids;
   typedef std::unordered_map<face_descriptor,Node_ids>             On_face_map;
   typedef std::unordered_map<edge_descriptor,Node_ids>             On_edge_map;
   //to keep the correspondance between node_id and vertex_handle in each mesh
   typedef internal::Node_id_to_vertex<TriangleMesh,
                                       doing_autorefinement||
                                       handle_non_manifold_features>
                                                              Node_id_to_vertex;
   typedef std::map<const TriangleMesh*, Node_id_to_vertex >   Mesh_to_map_node;
   //to handle coplanar halfedge of polyhedra that are full in the intersection
   typedef std::multimap<Node_id,halfedge_descriptor> Node_to_target_of_hedge_map;
   typedef std::map<TriangleMesh*,Node_to_target_of_hedge_map>
                                           Mesh_to_vertices_on_intersection_map;
   typedef std::unordered_map<vertex_descriptor,Node_id>      Vertex_to_node_id;
   typedef std::map<TriangleMesh*, Vertex_to_node_id> Mesh_to_vertex_to_node_id;
   typedef Non_manifold_feature_map<TriangleMesh>               NM_features_map;
// typedef for the CDT
   typedef Intersection_nodes<TriangleMesh, VertexPointMap1, VertexPointMap2,
            Predicates_on_constructions_needed>                              INodes;
   typedef typename INodes::Exact_kernel                                     EK;
    typedef Projection_traits_3<EK>                                          CDT_traits;
    typedef Triangulation_vertex_base_with_info_2<Node_id,CDT_traits>        Vb;
    typedef Constrained_triangulation_face_base_2<CDT_traits>                Fb;
    typedef Triangulation_data_structure_2<Vb,Fb>                         TDS_2;
    typedef Constrained_Delaunay_triangulation_2<CDT_traits,TDS_2>          CDT;
    typedef typename CDT::Vertex_handle                       CDT_Vertex_handle;
// data members
private:

  Graph_node_classifier<TriangleMesh, doing_autorefinement> graph_node_classifier;
  std::vector< std::vector<Node_id> > graph_of_constraints;
  boost::dynamic_bitset<> is_node_of_degree_one;
  //nb of intersection points between coplanar faces, see fixes XSL_TAG_CPL_VERT
  std::size_t number_coplanar_vertices;
  std::map<TriangleMesh*,On_face_map> on_face;
  std::map<TriangleMesh*,On_edge_map> on_edge;
  Mesh_to_vertices_on_intersection_map mesh_to_vertices_on_inter;
  Mesh_to_map_node mesh_to_node_id_to_vertex;
  // map an input vertex to a node id (input vertex on the intersection)
  Mesh_to_vertex_to_node_id mesh_to_vertex_to_node_id;

  std::map< Node_id,std::set<Node_id> > coplanar_constraints;

// optional data members to handle non-manifold issues
  std::map<const TriangleMesh*, const NM_features_map*> non_manifold_feature_maps;

//data members that require initialization in the constructor
  UserVisitor& user_visitor;
  OutputBuilder& output_builder;
  EdgeMarkMapBind marks_on_edges;
  bool input_with_coplanar_faces;
  TriangleMesh* const_mesh_ptr;

  template <class Ecm1, class Ecm2>
  void call_put(Ecm_bind<TriangleMesh, Ecm1, Ecm2>& ecm,
                TriangleMesh& tm, edge_descriptor ed, bool v)
  {
    ecm.call_put(tm, ed, v);
  }
  template <class Ecm>
  void call_put(Ecm& ecm,
                TriangleMesh&, edge_descriptor ed, bool v)
  {
    put(ecm, ed, v);
  }

  template <class Ecm1, class Ecm2>
  bool call_get(const Ecm_bind<TriangleMesh, Ecm1, Ecm2>& ecm,
                TriangleMesh& tm, edge_descriptor ed)
  {
    return ecm.call_get(tm, ed);
  }
  template <class Ecm>
  bool call_get(const Ecm& ecm,
                TriangleMesh&, edge_descriptor ed)
  {
    return get(ecm, ed);
  }
// visitor public functions
public:
  Surface_intersection_visitor_for_corefinement(
    UserVisitor& uv, OutputBuilder& o, const EdgeMarkMapBind& emm, TriangleMesh* const_mesh_ptr=nullptr)
    : number_coplanar_vertices(0)
    , user_visitor(uv)
    , output_builder(o)
    , marks_on_edges(emm)
    , input_with_coplanar_faces(false)
    , const_mesh_ptr(const_mesh_ptr)
  {}


  void start_filtering_intersections() const
  {
    user_visitor.start_filtering_intersections();
  }


  void progress_filtering_intersections(double d) const
  {
    user_visitor.progress_filtering_intersections(d);
  }

  void end_filtering_intersections() const
  {
    user_visitor.end_filtering_intersections();
  }


  void start_handling_edge_face_intersections(std::size_t i) const
  {
    user_visitor.start_handling_edge_face_intersections(i);
  }

  void edge_face_intersections_step() const
  {
    user_visitor.edge_face_intersections_step();
  }

  void end_handling_edge_face_intersections() const
  {
    user_visitor.end_handling_edge_face_intersections();
  }

  void start_handling_intersection_of_coplanar_faces(std::size_t i) const
  {
    user_visitor.start_handling_intersection_of_coplanar_faces(i);
  }

  void intersection_of_coplanar_faces_step() const
  {
    user_visitor.intersection_of_coplanar_faces_step();
  }

  void end_handling_intersection_of_coplanar_faces() const
  {
    user_visitor.end_handling_intersection_of_coplanar_faces();
  }

  void start_building_output() const
  {
    user_visitor.start_building_output();
  }

  void build_output_step() const
  {
    user_visitor.build_output_step();
  }

  void end_building_output() const
  {
    user_visitor.end_building_output();
  }

  void
  set_non_manifold_feature_map(
    const TriangleMesh& tm,
    const NM_features_map& nm)
  {
    non_manifold_feature_maps[&tm] = &nm;
  }

  void copy_nodes_ids_for_non_manifold_features()
  {
    static const constexpr std::size_t NM_NID((std::numeric_limits<std::size_t>::max)());

    for(const std::pair<const TriangleMesh* const, const NM_features_map*>& tm_and_nm :
        non_manifold_feature_maps)
    {
      TriangleMesh* tm_ptr = const_cast<TriangleMesh*>(tm_and_nm.first);
      // update nodes on edges
      On_edge_map& on_edge_map = on_edge[tm_ptr];
      std::vector< std::pair<std::size_t, const Node_ids*> > edges_to_copy;
      for (const std::pair<const edge_descriptor, Node_ids>& ed_and_ids : on_edge_map)
      {
        std::size_t eid = get(tm_and_nm.second->e_nm_id, ed_and_ids.first);
        if (eid!=NM_NID)
          edges_to_copy.push_back(std::make_pair(eid,&(ed_and_ids.second)));
      }
      for(const std::pair<std::size_t, const Node_ids*>& id_and_nodes : edges_to_copy)
      {
        const std::vector<edge_descriptor>& nm_edges =
          tm_and_nm.second->non_manifold_edges[id_and_nodes.first];
        CGAL_assertion( on_edge_map.count(nm_edges.front())==1 );

        for (std::size_t i=1; i<nm_edges.size(); ++i)
          on_edge_map[nm_edges[i]] = *id_and_nodes.second;
      }

      // update map vertex -> node_id
      Vertex_to_node_id& vertex_to_node_id = mesh_to_vertex_to_node_id[tm_ptr];
      Node_to_target_of_hedge_map& vertices_on_inter = mesh_to_vertices_on_inter[tm_ptr];

      std::vector< std::pair<vertex_descriptor, Node_id> > vertices_to_add;
      for (const typename std::pair<const vertex_descriptor, Node_id>& vd_and_id
          : vertex_to_node_id)
      {
        std::size_t vid = get(tm_and_nm.second->v_nm_id, vd_and_id.first);
        if (vid!=NM_NID)
          vertices_to_add.push_back(std::make_pair(vd_and_id.first,vd_and_id.second));
      }

      for(const std::pair<vertex_descriptor, Node_id>& vd_and_nid : vertices_to_add)
      {
        std::size_t vid = get(tm_and_nm.second->v_nm_id, vd_and_nid.first);
        for(vertex_descriptor vd : tm_and_nm.second->non_manifold_vertices[vid])
        {
          if (vd != vd_and_nid.first)
          {
            vertex_to_node_id.insert(std::make_pair(vd,vd_and_nid.second));
            output_builder.set_vertex_id(vd, vd_and_nid.second, *tm_ptr);
            vertices_on_inter.insert(std::make_pair(vd_and_nid.second,halfedge(vd,*tm_ptr)));
          }
        }
      }
    }
  }

  template<class Graph_node>
  void annotate_graph(std::vector<Graph_node>& graph)
  {
    std::size_t nb_nodes=graph.size();
    graph_of_constraints.resize(nb_nodes);
    is_node_of_degree_one.resize(nb_nodes);
//TODO: pas bon avec autoref et la collecte des infos aussi...

    graph_node_classifier.preprocessing();
    for(std::size_t node_id=0;node_id<nb_nodes;++node_id)
    {
      graph_of_constraints[node_id].assign(
        graph[node_id].neighbors.begin(),
        graph[node_id].neighbors.end());

      if (graph_of_constraints[node_id].size()==1)
        is_node_of_degree_one.set(node_id);

      if (handle_non_manifold_features) continue; // skip the rest of the function meant to prepare
                                                  // the intersection polylines for the output builder
                                                  // (which is not called with non-manifold features)

      // mark every vertex contained by the intersection polyline, incident to a border and a non-border edge.
      // The logic is somehow equivalent to what was done with the container `non_manifold_nodes`
      // that was used to split polylines at certains points. Non-manifold was used in the context
      // of the combinatorial map where the import inside the combinatorial map was possible (see broken_bound-[12].off)
      if (  graph_node_classifier.is_terminal(node_id, graph_of_constraints[node_id]) )
        graph[node_id].make_terminal();
    }
  }

  void start_new_polyline(Node_id i, Node_id j)
  {
    if ( i==j ) //case of a single point
    {
      // TODO shall we insert the point?????
      //TAG SL001
      //nothing is done
      return;
    }
    output_builder.start_new_polyline(i,j);
  }

  void add_node_to_polyline(Node_id i)
  {
    output_builder.add_node_to_polyline(i);
  }

  void set_number_of_intersection_points_from_coplanar_faces(std::size_t n)
  {
      number_coplanar_vertices=n;
  }

  void input_have_coplanar_faces()
  {
    input_with_coplanar_faces=true;
  }

  void update_terminal_nodes(std::vector<bool>&)
  {
    CGAL_error_msg("This function should not be called");
  }

  void check_node_on_boundary_edge_case(std::size_t node_id,
                                      halfedge_descriptor h,
                                      const TriangleMesh& tm)
  {
    graph_node_classifier.node_on_edge(node_id, h, tm);
  }

  void check_node_on_boundary_vertex_case(std::size_t node_id,
                                        halfedge_descriptor h,
                                        const TriangleMesh& tm)
  {
    graph_node_classifier.node_on_vertex(node_id, target(h,tm), tm);
  }

  //keep track of the fact that a polyhedron original vertex is a node
  void all_incident_faces_got_a_node_as_vertex(
      halfedge_descriptor h,
      Node_id node_id,
      TriangleMesh& tm)
  {
    CGAL_assertion_code(bool insert_ok = )
    mesh_to_vertex_to_node_id[&tm].insert(std::make_pair(target(h,tm),node_id))
    CGAL_assertion_code(.second);
    CGAL_assertion(insert_ok || mesh_to_vertex_to_node_id[&tm][target(h,tm)]==node_id);
  }

  void new_node_added_triple_face(std::size_t node_id,
                                  face_descriptor f1,
                                  face_descriptor f2,
                                  face_descriptor f3,
                                  const TriangleMesh& tm) // TODO check if we need a special case if the endpoint of the intersect edge is on the third face
  {
    CGAL_assertion(f1!=f2 && f1!=f3 && f2!=f3);
    graph_node_classifier.new_node(node_id, tm);

//    user_visitor.new_node_added_triple_face(node_id, f1, f2, f3, tm); // NODE_VISITOR_TAG
#ifdef CGAL_DEBUG_AUTOREFINEMENT
    std::cout << "adding node " << node_id << " " << f1 << " " << f2 << " " << f3 << "\n";
#endif
    TriangleMesh* tm_ptr = const_cast<TriangleMesh*>(&tm);
    on_face[tm_ptr][f1].push_back(node_id);
    on_face[tm_ptr][f2].push_back(node_id);
    on_face[tm_ptr][f3].push_back(node_id);
  }

  void new_node_added(std::size_t node_id,
                      Intersection_type type,
                      halfedge_descriptor h_1,
                      halfedge_descriptor h_2,
                      const TriangleMesh& tm1,
                      const TriangleMesh& tm2,
                      bool is_target_coplanar,
                      bool is_source_coplanar)
  {
    TriangleMesh* tm1_ptr = const_cast<TriangleMesh*>(&tm1);
    TriangleMesh* tm2_ptr = const_cast<TriangleMesh*>(&tm2);
    graph_node_classifier.new_node(node_id, *tm1_ptr);
    graph_node_classifier.new_node(node_id, *tm2_ptr);

    //forward to the visitor
    user_visitor.intersection_point_detected(node_id, type, h_1, h_2, tm1, tm2, is_target_coplanar, is_source_coplanar);
    if (tm2_ptr!=const_mesh_ptr)
    {
      switch(type)
      {
        case ON_FACE: //Face intersected by an edge
          on_face[tm2_ptr][face(h_2,tm2)].push_back(node_id);
        break;
        case ON_EDGE: //Edge intersected by an edge
        {
          on_edge[tm2_ptr][edge(h_2,tm2)].push_back(node_id);
          check_node_on_boundary_edge_case(node_id,h_2,tm2);
        }
        break;
        case ON_VERTEX:
        {
          //grab original vertex that is on commom intersection
          mesh_to_vertices_on_inter[tm2_ptr].insert(std::make_pair(node_id,h_2));
          Node_id_to_vertex& node_id_to_vertex=mesh_to_node_id_to_vertex[tm2_ptr];
          if (node_id_to_vertex.size()<=node_id)
            node_id_to_vertex.resize(node_id+1);
          node_id_to_vertex.register_vertex(node_id, target(h_2,tm2));
          all_incident_faces_got_a_node_as_vertex(h_2,node_id,*tm2_ptr);
          check_node_on_boundary_vertex_case(node_id,h_2,tm2);
          output_builder.set_vertex_id(target(h_2, tm2), node_id, tm2);
        }
        break;
        default:
        return;
      }
    }

    if (tm1_ptr==const_mesh_ptr)
      return;

    CGAL_assertion(!is_target_coplanar || !is_source_coplanar); //coplanar edge are not forwarded

    if ( is_target_coplanar )
    {
      //grab original vertex that is on commom intersection
      mesh_to_vertices_on_inter[tm1_ptr].insert(std::make_pair(node_id,h_1));
      Node_id_to_vertex& node_id_to_vertex=mesh_to_node_id_to_vertex[tm1_ptr];
      if (node_id_to_vertex.size()<=node_id)
        node_id_to_vertex.resize(node_id+1);
      node_id_to_vertex.register_vertex(node_id, target(h_1,tm1));
      all_incident_faces_got_a_node_as_vertex(h_1,node_id, *tm1_ptr);
      // register the vertex in the output builder
      output_builder.set_vertex_id(target(h_1, tm1), node_id, tm1);
      check_node_on_boundary_vertex_case(node_id,h_1,tm1);
    }
    else{
      if ( is_source_coplanar ){
        //grab original vertex that is on commom intersection
        halfedge_descriptor h_1_opp=opposite(h_1,tm1);
        mesh_to_vertices_on_inter[tm1_ptr].insert(std::make_pair(node_id,h_1_opp));
        Node_id_to_vertex& node_id_to_vertex=mesh_to_node_id_to_vertex[tm1_ptr];
        if(node_id_to_vertex.size()<=node_id)
          node_id_to_vertex.resize(node_id+1);
        node_id_to_vertex.register_vertex(node_id, source(h_1,tm1));
        all_incident_faces_got_a_node_as_vertex(h_1_opp,node_id, *tm1_ptr);
        // register the vertex in the output builder
        output_builder.set_vertex_id(source(h_1, tm1), node_id, tm1);
        check_node_on_boundary_vertex_case(node_id,h_1_opp,tm1);
      }
      else{
        //handle intersection on principal edge
        typename std::map<const TriangleMesh*, const NM_features_map*>::iterator it_find =
          non_manifold_feature_maps.find(&tm1);
        if ( it_find != non_manifold_feature_maps.end() )
        {
          // update h_1 if it is not the canonical non-manifold edge
          // This is important to make sure intersection points on non-manifold
          // edges are all connected for the same edge so that the redistribution
          // on other edges does not overwrite some nodes.
          // This update might be required in case of EDGE-EDGE intersection or
          // COPLANAR intersection.
          const NM_features_map& nm_features_map_1 = *it_find->second;
          std::size_t eid1 = nm_features_map_1.non_manifold_edges.empty()
                           ? std::size_t(-1)
                           : get(nm_features_map_1.e_nm_id, edge(h_1, tm1));

          if (eid1 != std::size_t(-1))
          {
            if ( edge(h_1, tm1) != nm_features_map_1.non_manifold_edges[eid1].front() )
              h_1 = halfedge(nm_features_map_1.non_manifold_edges[eid1].front(), tm1);
          }
        }

        on_edge[tm1_ptr][edge(h_1,tm1)].push_back(node_id);
        check_node_on_boundary_edge_case(node_id,h_1,tm1);
      }
    }
  }

  //sort node ids so that we can split the hedge
  //consecutively
  template <class VPM, class Node_vector>
  void sort_vertices_along_hedge(std::vector<std::size_t>& node_ids,
                                 halfedge_descriptor hedge,
                                 const TriangleMesh& tm,
                                 const VPM& vpm,
                                 const Node_vector& nodes)
  {
    std::sort(node_ids.begin(),
              node_ids.end(),
              Less_along_a_halfedge<TriangleMesh, VPM, Node_vector>
                (hedge, tm, vpm, nodes)
    );
  }

  struct Face_boundary{
    std::vector<std::size_t> node_ids_array[3]; // the node_ids on each halfedges
    std::map<halfedge_descriptor,int> hedges_ids;
    halfedge_descriptor halfedges[3]; //the three halfedges of the original face
    vertex_descriptor   vertices[3];  //the three vertices  of the original face
    //node_ids_array[0] corresponds to the original edge vertices[0],vertices[1] = halfedges[0]
    //node_ids_array[1] corresponds to the original edge vertices[1],vertices[2] = halfedges[1]
    //node_ids_array[2] corresponds to the original edge vertices[2],vertices[0] = halfedges[2]
    Face_boundary(halfedge_descriptor first, TriangleMesh& tm)
    {
      CGAL_assertion(is_triangle(first,tm));
      halfedges[0]=first;
      halfedges[1]=next(first,tm);
      halfedges[2]=next(halfedges[1],tm);

      vertices[0]=source(halfedges[0],tm);
      vertices[1]=source(halfedges[1],tm);
      vertices[2]=source(halfedges[2],tm);

      hedges_ids.insert(std::make_pair(halfedges[0],0));
      hedges_ids.insert(std::make_pair(halfedges[1],1));
      hedges_ids.insert(std::make_pair(halfedges[2],2));
    }

    //used when object was created with hedge but opposite was used to split the original face
    void update_original_halfedge(halfedge_descriptor original,
                                  halfedge_descriptor new_hedge,
                                  TriangleMesh& /*tm*/)
    {
      typename std::map<halfedge_descriptor,int>::iterator it_id =
        hedges_ids.find(original);
      CGAL_assertion(it_id!=hedges_ids.end());
      int index=it_id->second;
      CGAL_assertion(halfedges[index]==original);
      hedges_ids.erase(it_id);
      hedges_ids.insert(std::make_pair(new_hedge,index));
      halfedges[index]=new_hedge;
    }

    template <class Iterator>
    void copy_node_ids(halfedge_descriptor hedge,Iterator begin,Iterator end)
    {
      typename std::map<halfedge_descriptor,int>::iterator it_id =
        hedges_ids.find(hedge);
      CGAL_assertion(it_id!=hedges_ids.end());
      std::copy(begin,end,std::back_inserter(node_ids_array[it_id->second]));
    }

    // Used by the autorefinement and non-manifold edge handling to re-set
    // the id of nodes on the boundary of a face since another vertex
    // (inside a face or on another edge) might have
    // overwritten the vertex in node_id_to_vertex
    template <class Node_id_to_vertex>
    void update_node_id_to_vertex_map(Node_id_to_vertex& node_id_to_vertex,
                                      TriangleMesh& tm)
    {
      for (int i=0; i<3; ++i)
      {
        halfedge_descriptor h = halfedges[(i+2)%3];
        h = next(h, tm);
        for(std::size_t id : node_ids_array[i])
        {
          // needed when we triangulate a face --> need to pick the right vertex
          node_id_to_vertex.set_vertex_for_retriangulation(id, target(h, tm));
          h = next(h, tm);
        }
        CGAL_assertion(h ==  halfedges[i]);
      }
    }

  };

  typedef std::unordered_map<face_descriptor,Face_boundary>  Face_boundaries;

  //update the id of input mesh vertex that are also a node
  void update_face_indices(
    std::array<vertex_descriptor,3>& f_vertices,
    std::array<Node_id,3>& f_indices,
    Vertex_to_node_id& vertex_to_node_id)
  {
    for (int k=0;k<3;++k){
      typename std::unordered_map<vertex_descriptor,Node_id>::iterator it =
        vertex_to_node_id.find(f_vertices[k]);
      if (it!=vertex_to_node_id.end())
        f_indices[k]=it->second;
    }
  }

  //insert intersection edge as constrained edges in a CDT triangulation
  void insert_constrained_edges_coplanar_case(
    Node_id node_id,
    CDT& cdt,
    std::map<Node_id,CDT_Vertex_handle>& id_to_CDT_vh)
  {
    if (node_id < number_coplanar_vertices){
      //XSL_TAG_CPL_VERT
      // Insert constrained edges from coplanar faces that have been
      // retriangulated. This ensure that triangulations are compatible.
      // This edges were not constrained in the first mesh but are in the
      // second (ensuring compatibility)
      typename std::map< Node_id,std::set<Node_id> >::iterator it_neighbors =
        coplanar_constraints.find(node_id);
      if (it_neighbors!=coplanar_constraints.end())
      {
        CDT_Vertex_handle vh=id_to_CDT_vh[node_id];
        for(Node_id id :it_neighbors->second)
        {
          typename std::map<Node_id,CDT_Vertex_handle>
            ::iterator it_vh=id_to_CDT_vh.find(id);
          // this condition ensures to consider only graph edges that are in
          // the same triangle (not in a neighbor one when involving node on
          // a triangle edge) here we can't make the difference between a point
          // on the interior or the boundary, so points_on_triangle is not used.
          if ( it_vh!=id_to_CDT_vh.end() )
            cdt.insert_constraint(vh,it_vh->second);
        }
      }
    }
  }

  //insert intersection edges as constrained edges in a CDT triangulation
  void insert_constrained_edges(
    Node_ids& node_ids,
    CDT& cdt,
    std::map<Node_id, CDT_Vertex_handle>& id_to_CDT_vh,
    std::vector<std::pair<Node_id,Node_id> >& constrained_edges,
    bool points_on_triangle=false)
  {
    for(Node_id id : node_ids)
    {
      CGAL_assertion(id < graph_of_constraints.size());
      std::vector<Node_id>& neighbors=graph_of_constraints[id];
      if (!neighbors.empty())
      {
        CDT_Vertex_handle vh=id_to_CDT_vh.find(id)->second;
        for(Node_id id_n :neighbors)
        {
        //   if (id_n < id) continue; //no need to do it twice
          typename std::map<Node_id,CDT_Vertex_handle>
            ::iterator it_vh=id_to_CDT_vh.find(id_n);
          // this condition ensures to consider only graph edges that are in
          // the same triangle
          if ( !points_on_triangle || it_vh!=id_to_CDT_vh.end() ){
            CGAL_assertion(doing_autorefinement || handle_non_manifold_features || it_vh!=id_to_CDT_vh.end());
            if (it_vh==id_to_CDT_vh.end()) continue; // needed for autorefinement (interior nodes)
            cdt.insert_constraint(vh,it_vh->second);
            constrained_edges.push_back(std::make_pair(id,id_n));
            constrained_edges.push_back(std::make_pair(id_n,id));
          }
        }
      }
      #ifdef CGAL_COREFINEMENT_DEBUG
      else
        std::cout << "X0: Found an isolated point" << std::endl;
      #endif
      insert_constrained_edges_coplanar_case(id,cdt,id_to_CDT_vh);
    }
  }

  // insert a point on the convex_hull using the infinite incident face to the edge the point is inserted into
  // warning: fh is updated
  typename CDT::Vertex_handle
  insert_point_on_ch_edge(CDT& cdt, typename CDT::Face_handle& fh, const typename CDT::Point& p)
  {
    CGAL_assertion(cdt.is_infinite(fh));
    int fi = fh->index(cdt.infinite_vertex());
    typename CDT::Vertex_handle vh = cdt.insert(p, CDT::EDGE, fh, fi);
    typename CDT::Edge_circulator ec = cdt.incident_edges(vh);
    while(ec->first->vertex(CDT::ccw(ec->second)) != cdt.infinite_vertex()){
      ++ec;
    }
    fh = ec->first->neighbor(ec->second);
    CGAL_assertion( cdt.is_valid() );
    return vh;
  }

  template <class OnEdgeMapIterator, class VPM>
  void split_halfedges(OnEdgeMapIterator it,
                       const VPM& vpm,
                       INodes& nodes,
                       std::map<TriangleMesh*, Face_boundaries>& mesh_to_face_boundaries)
  {
    TriangleMesh& tm=*it->first;
    CGAL_assertion(&tm!=const_mesh_ptr);

    On_edge_map& on_edge_map=it->second;
    On_face_map& on_face_map=on_face[&tm];
    Face_boundaries& face_boundaries=mesh_to_face_boundaries[&tm];

    for(typename On_edge_map::iterator it2=on_edge_map.begin();
                                       it2!=on_edge_map.end();
                                       ++it2)
    {
      //the edge to be split
      halfedge_descriptor hedge=halfedge(it2->first,tm);
      //indices of the nodes to be inserted
      Node_ids& node_ids=it2->second;
      CGAL_assertion( std::set<Node_id>(node_ids.begin(), node_ids.end())
                        .size()==node_ids.size() );
      //sort nodes along the egde to allow consecutive splits
      sort_vertices_along_hedge(node_ids,hedge,tm,vpm,nodes);

      //save original face and nodes for face of hedge (1)
      if ( !is_border(hedge,tm) ){
        face_descriptor f=face(hedge,tm);
        typename Face_boundaries::iterator it_face = face_boundaries.find(f);
        if (it_face==face_boundaries.end())
          it_face=face_boundaries.insert(std::make_pair(f,Face_boundary(hedge,tm))).first;
        it_face->second.copy_node_ids(hedge,node_ids.begin(),node_ids.end());
      }

      //save original face and nodes for face of hedge->opposite (2)
      typename Face_boundaries::iterator opposite_original_info=face_boundaries.end();
      halfedge_descriptor hedge_opp = opposite(hedge,tm);
      if ( !is_border(hedge_opp,tm) ){
        face_descriptor f=face(hedge_opp,tm);
        opposite_original_info=face_boundaries.find(f);
        if (opposite_original_info==face_boundaries.end())
          opposite_original_info=face_boundaries.insert(std::make_pair(f,Face_boundary(hedge_opp,tm))).first;
        opposite_original_info->second.copy_node_ids(hedge_opp,node_ids.rbegin(),node_ids.rend());
      }

      typename Mesh_to_map_node::iterator it_map=mesh_to_node_id_to_vertex.find(&tm);
      CGAL_assertion(it_map!=mesh_to_node_id_to_vertex.end());
      //a map to identify the vertex in the polyhedron corresponding to an intersection point
      Node_id_to_vertex& node_id_to_vertex=it_map->second;

      CGAL_assertion_code(vertex_descriptor original_vertex=source(hedge,tm);)

      //We need an edge incident to the source vertex of hedge. This is the first opposite edge created.
      bool first=true;
      halfedge_descriptor hedge_incident_to_src=Graph_traits::null_halfedge();
      bool hedge_is_marked = call_get(marks_on_edges,tm,edge(hedge,tm));
      //do split the edges
      CGAL_assertion_code(vertex_descriptor expected_src=source(hedge,tm));
      user_visitor.before_edge_split(hedge, tm);
      for(std::size_t node_id : node_ids)
      {
        halfedge_descriptor hnew = Euler::split_edge(hedge, tm);
        CGAL_assertion(expected_src==source(hnew,tm));
        vertex_descriptor vnew=target(hnew,tm);
        nodes.call_put(vpm, vnew, node_id, tm);
        // register the new vertex in the output builder
        output_builder.set_vertex_id(vnew, node_id, tm);
        node_id_to_vertex.register_vertex(node_id, vnew);
        if (first){
          first=false;
          hedge_incident_to_src=next(opposite(hedge,tm),tm);
        }

        //update marker tags. If the edge was marked, then the resulting edges in the split must be marked
        if ( hedge_is_marked )
          call_put(marks_on_edges,tm,edge(hnew,tm),true);
        user_visitor.new_vertex_added(node_id, target(hnew, tm), tm);
        user_visitor.edge_split(hnew, tm);

        CGAL_assertion_code(expected_src=vnew);
      }
      user_visitor.after_edge_split();

      CGAL_assertion(target(hedge_incident_to_src,tm)==original_vertex);
      CGAL_assertion(face(hedge_incident_to_src,tm)==face(hedge_opp,tm));

      //save original face and nodes for face of hedge->opposite (2)
      if ( !is_border(hedge_opp,tm) ){
        CGAL_assertion(opposite_original_info!=face_boundaries.end());
        opposite_original_info->second.update_original_halfedge(
          hedge_opp,hedge_incident_to_src,tm);
      }

      //insert the two incident faces in on_face map so that they will be triangulated.
      if (!is_border(hedge,tm)) on_face_map[face(hedge,tm)];
      if (!is_border(hedge_opp,tm)) on_face_map[face(hedge_opp,tm)];
    }
  }

  template <class OnFaceMapIterator, class VPM>
  void triangulate_intersected_faces(OnFaceMapIterator it,
                                     const VPM& vpm,
                                     INodes& nodes,
                                     std::map<TriangleMesh*, Face_boundaries>& mesh_to_face_boundaries)
  {
    TriangleMesh& tm=*it->first;
    CGAL_assertion(&tm!=const_mesh_ptr);

    On_face_map& on_face_map=it->second;
    Face_boundaries& face_boundaries=mesh_to_face_boundaries[&tm];
    Node_id_to_vertex& node_id_to_vertex=mesh_to_node_id_to_vertex[&tm];
    Vertex_to_node_id& vertex_to_node_id=mesh_to_vertex_to_node_id[&tm];

    const Node_id nb_nodes = nodes.size();

    for (typename On_face_map::iterator it=on_face_map.begin();
          it!=on_face_map.end();++it)
    {
      user_visitor.triangulating_faces_step();
      face_descriptor f = it->first; //the face to be triangulated
      Node_ids& node_ids  = it->second; // ids of nodes in the interior of f
      typename Face_boundaries::iterator it_fb=face_boundaries.find(f);

      std::map<Node_id,typename CDT::Vertex_handle> id_to_CDT_vh;

      //associate an edge of the triangulation to a halfedge in a given polyhedron
      std::map<std::pair<Node_id,Node_id>,halfedge_descriptor> edge_to_hedge;

      // the vertices of f
      std::array<vertex_descriptor,3> f_vertices;
      // the node_id of an input vertex or a fake id (>=nb_nodes)
      std::array<Node_id,3> f_indices = {{nb_nodes,nb_nodes+1,nb_nodes+2}};
      if (it_fb!=face_boundaries.end()){ //the boundary of the triangle face was refined
        f_vertices[0]=it_fb->second.vertices[0];
        f_vertices[1]=it_fb->second.vertices[1];
        f_vertices[2]=it_fb->second.vertices[2];
        update_face_indices(f_vertices,f_indices,vertex_to_node_id);
        if (doing_autorefinement || handle_non_manifold_features)
          it_fb->second.update_node_id_to_vertex_map(node_id_to_vertex, tm);
      }
      else{
        CGAL_assertion( is_triangle(halfedge(f,tm),tm) );
        halfedge_descriptor h0=halfedge(f,tm), h1=next(h0,tm), h2=next(h1,tm);
        f_vertices[0]=target(h0,tm); //nb_nodes
        f_vertices[1]=target(h1,tm); //nb_nodes+1
        f_vertices[2]=target(h2,tm); //nb_nodes+2

        update_face_indices(f_vertices,f_indices,vertex_to_node_id);
        edge_to_hedge[std::make_pair( f_indices[2],f_indices[0] )] = h0;
        edge_to_hedge[std::make_pair( f_indices[0],f_indices[1] )] = h1;
        edge_to_hedge[std::make_pair( f_indices[1],f_indices[2] )] = h2;
      }

      // handle possible presence of degenerate faces
      if (const_mesh_ptr && collinear( get(vpm,f_vertices[0]), get(vpm,f_vertices[1]), get(vpm,f_vertices[2]) ) )
      {
        Node_ids face_vertex_nids;

        //check if one of the triangle input vertex is also a node
        for (int ik=0;ik<3;++ik)
          if ( f_indices[ik]<nb_nodes )
            face_vertex_nids.push_back(f_indices[ik]);

        // collect nodes on edges (if any)
        if (it_fb != face_boundaries.end())
        {
          Face_boundary& f_boundary=it_fb->second;
          for (int i=0;i<3;++i)
            std::copy(f_boundary.node_ids_array[i].begin(),
                      f_boundary.node_ids_array[i].end(),
                      std::back_inserter(face_vertex_nids));
        }

        std::sort(face_vertex_nids.begin(), face_vertex_nids.end());
        std::vector<std::array<std::pair<halfedge_descriptor,Node_id>,2>> constraints;
        for(Node_id id : face_vertex_nids)
        {
          CGAL_assertion(id < graph_of_constraints.size());
          const std::vector<Node_id>& neighbors=graph_of_constraints[id];
          if (!neighbors.empty())
          {
            for(Node_id id_n :neighbors)
            {
              if (id_n<id) continue;
              if (std::binary_search(face_vertex_nids.begin(), face_vertex_nids.end(), id_n))
              {
                vertex_descriptor vi = node_id_to_vertex.get_vertex(id),
                                  vn = node_id_to_vertex.get_vertex(id_n);
                bool is_face_border = false;
                halfedge_descriptor h;

                std::tie(h, is_face_border) = halfedge(vi,vn, tm);
                if (is_face_border)
                {
                  call_put(marks_on_edges,tm,edge(h,tm),true);
                  output_builder.set_edge_per_polyline(tm,std::make_pair(id, id_n),h);
                }
                else
                {
                  halfedge_descriptor hi=halfedge(vi, tm);
                  while(face(hi, tm) != f)
                    hi=opposite(next(hi, tm), tm);

                  halfedge_descriptor hn=halfedge(vn, tm);
                  while(face(hn, tm) != f)
                    hn=opposite(next(hn, tm), tm);
                  constraints.emplace_back(make_array(std::make_pair(hi,id),std::make_pair(hn, id_n)));
                }
              }
            }
          }
          #ifdef CGAL_COREFINEMENT_DEBUG
          else
            std::cout << "X0bis: Found an isolated point" << std::endl;
          #endif
        }

        CGAL_assertion(constraints.empty() || it_fb != face_boundaries.end());
        std::vector<face_descriptor> new_faces;
        for (const std::array<std::pair<halfedge_descriptor, Node_id>, 2>& a : constraints)
        {
          halfedge_descriptor nh = Euler::split_face(a[0].first, a[1].first, tm);
          new_faces.push_back(face(opposite(nh, tm), tm));

          call_put(marks_on_edges,tm,edge(nh,tm),true);
          output_builder.set_edge_per_polyline(tm,std::make_pair(a[0].second, a[1].second),nh);
        }

        // now triangulate new faces
        if (!new_faces.empty())
        {
          new_faces.push_back(f);
          for(face_descriptor nf : new_faces)
          {
            halfedge_descriptor h = halfedge(nf, tm),
                                nh = next(next(h,tm),tm);
            while(next(nh, tm)!=h)
              nh=next(Euler::split_face(h, nh, tm), tm);
          }
        }

        continue;
      }

      typename EK::Point_3 p = nodes.to_exact(get(vpm,f_vertices[0])),
                           q = nodes.to_exact(get(vpm,f_vertices[1])),
                           r = nodes.to_exact(get(vpm,f_vertices[2]));
///TODO use a positive normal and remove all workaround to guarantee that triangulation of coplanar patches are compatible
      CDT_traits traits(typename EK::Construct_normal_3()(p,q,r));
      CDT cdt(traits);

      // insert triangle points
      std::array<CDT_Vertex_handle,3> triangle_vertices;
      //we can do this to_exact because these are supposed to be input points.
      triangle_vertices[0]=cdt.insert_outside_affine_hull(p);
      triangle_vertices[1]=cdt.insert_outside_affine_hull(q);
      triangle_vertices[2]=cdt.tds().insert_dim_up(cdt.infinite_vertex(), false);
      triangle_vertices[2]->set_point(r);

      triangle_vertices[0]->info()=f_indices[0];
      triangle_vertices[1]->info()=f_indices[1];
      triangle_vertices[2]->info()=f_indices[2];

      node_id_to_vertex.set_temporary_vertex_for_retriangulation(nb_nodes, f_vertices[0]);
      node_id_to_vertex.set_temporary_vertex_for_retriangulation(nb_nodes+1, f_vertices[1]);
      node_id_to_vertex.set_temporary_vertex_for_retriangulation(nb_nodes+2, f_vertices[2]);

      //if one of the triangle input vertex is also a node
      for (int ik=0;ik<3;++ik){
        if ( f_indices[ik]<nb_nodes )
        {
          id_to_CDT_vh.insert(
              std::make_pair(f_indices[ik],triangle_vertices[ik]));
          if (doing_autorefinement || handle_non_manifold_features)
            // update the current vertex in node_id_to_vertex
            // to match the one of the face
            node_id_to_vertex.set_temporary_vertex_for_retriangulation(f_indices[ik], f_vertices[ik]);
            // Note on set_temporary_vertex instead of set_vertex: here since the point is an input point
            // it is OK not to store all vertices corresponding to this id as the approximate version
            // is already tight and the call in Intersection_nodes::finalize() will not fix anything
        }
      }
      //insert points on edges
      if (it_fb!=face_boundaries.end()) //if f not a triangle?
      {
        // collect infinite faces incident to the initial triangle
        typename CDT::Face_handle infinite_faces[3];
        for (int i=0;i<3;++i)
        {
          int oi=-1;
          CGAL_assertion_code(bool is_edge = )
          cdt.is_edge(triangle_vertices[i], triangle_vertices[(i+1)%3], infinite_faces[i], oi);
          CGAL_assertion(is_edge);
          CGAL_assertion( cdt.is_infinite( infinite_faces[i]->vertex(oi) ) );
        }

        // In this loop, for each original edge of the triangle, we insert
        // the constrained edges and we recover the halfedge_descriptor
        // corresponding to these constrained (they are already in tm)
        Face_boundary& f_boundary=it_fb->second;
        for (int i=0;i<3;++i){
          //handle case of halfedge starting at triangle_vertices[i]
          // and ending at triangle_vertices[(i+1)%3]

          const Node_ids& ids_on_edge=f_boundary.node_ids_array[i];
          CDT_Vertex_handle previous=triangle_vertices[i];
          Node_id prev_index=f_indices[i];// node-id of the mesh vertex
          halfedge_descriptor hedge = next(f_boundary.halfedges[(i+2)%3],tm);
          CGAL_assertion( source(hedge,tm)==f_boundary.vertices[i] );
          if (!ids_on_edge.empty()){ //is there at least one node on this edge?
            // fh must be an infinite face
            // The points must be ordered from fh->vertex(cw(infinite_vertex)) to fh->vertex(ccw(infinite_vertex))
            for(Node_id id : ids_on_edge)
            {
              CDT_Vertex_handle vh=insert_point_on_ch_edge(cdt,infinite_faces[i],nodes.exact_node(id));
              vh->info()=id;
              id_to_CDT_vh.insert(std::make_pair(id,vh));
              edge_to_hedge[std::make_pair(prev_index,id)]=hedge;
              previous=vh;
              hedge=next(hedge,tm);
              prev_index=id;
            }
          }
          else{
          CGAL_assertion_code(halfedge_descriptor hd=f_boundary.halfedges[i]);
            CGAL_assertion( target(hd,tm) == f_boundary.vertices[(i+1)%3] );
            CGAL_assertion( source(hd,tm) == f_boundary.vertices[ i ] );
          }
          CGAL_assertion(hedge==f_boundary.halfedges[i]);
          edge_to_hedge[std::make_pair(prev_index,f_indices[(i+1)%3])] =
            it_fb->second.halfedges[i];
        }
      }

      //insert point inside face
      for(Node_id node_id : node_ids)
      {
        CDT_Vertex_handle vh=cdt.insert(nodes.exact_node(node_id));
        vh->info()=node_id;
        id_to_CDT_vh.insert(std::make_pair(node_id,vh));
      }

      std::vector<std::pair<Node_id,Node_id> > constrained_edges;

      // insert constraints that are interior to the triangle (in the case
      // no edges are collinear in the meshes)
      insert_constrained_edges(node_ids,cdt,id_to_CDT_vh,constrained_edges);

      // insert constraints between points that are on the boundary
      // (not a contrained on the triangle boundary)
      if (it_fb!=face_boundaries.end()) //is f not a triangle ?
      {
        for (int i=0;i<3;++i)
        {
          Node_ids& ids=it_fb->second.node_ids_array[i];
          insert_constrained_edges(ids,cdt,id_to_CDT_vh,constrained_edges,1);
        }
      }

      //insert coplanar edges for endpoints of triangles
      for (int i=0;i<3;++i){
        Node_id nindex=triangle_vertices[i]->info();
        if ( nindex < nb_nodes )
          insert_constrained_edges_coplanar_case(nindex,cdt,id_to_CDT_vh);
      }

      //XSL_TAG_CPL_VERT
      //collect edges incident to a point that is the intersection of two
      // coplanar faces. This ensure that triangulations are compatible.
      if (it_fb!=face_boundaries.end()) //is f not a triangle ?
      {
        for (typename CDT::Finite_vertices_iterator
              vit=cdt.finite_vertices_begin(),
              vit_end=cdt.finite_vertices_end();vit_end!=vit;++vit)
        {
          //skip original vertices (that are not nodes) and non-coplanar face
          // issued vertices (this is working because intersection points
          // between coplanar facets are the first inserted)
          if (vit->info() >= nb_nodes ||
              vit->info() >= number_coplanar_vertices) continue;
          // \todo no need to insert constrained edges (they also are constrained
          // in the other mesh)!!
          typename std::map< Node_id,std::set<Node_id> >::iterator res =
              coplanar_constraints.insert(
                  std::make_pair(vit->info(),std::set<Node_id>())).first;
          //turn around the vertex and get incident edge
          typename CDT::Edge_circulator  start=cdt.incident_edges(vit);
          typename CDT::Edge_circulator  curr=start;
          do{
            if (cdt.is_infinite(*curr) ) continue;
            typename CDT::Edge mirror=cdt.mirror_edge(*curr);
            if ( cdt.is_infinite( curr->first->vertex(curr->second) ) ||
                 cdt.is_infinite( mirror.first->vertex(mirror.second) ) )
              continue; // skip edges that are on the boundary of the triangle
                        // (these are already constrained)
            //insert edges in the set of constraints
            CDT_Vertex_handle vh=vit;
            int nindex = curr->first->vertex((curr->second+1)%3)==vh
                           ? (curr->second+2)%3
                           : (curr->second+1)%3;
            CDT_Vertex_handle vn=curr->first->vertex(nindex);
            if ( vit->info() > vn->info() || vn->info()>=nb_nodes)
              continue; //take only one out of the two edges + skip input
            CGAL_assertion(vn->info()<nb_nodes);
            res->second.insert( vn->info() );
          }while(start!=++curr);
        }
      }

      // import the triangle in `cdt` in the face `f` of `tm`
      triangulate_a_face(f, tm, nodes, node_ids, node_id_to_vertex,
        edge_to_hedge, cdt, vpm, output_builder, user_visitor);

      // TODO Here we do the update only for internal edges.
      // Update for border halfedges could be done during the split

      //3) mark halfedges that are common to two polyhedral surfaces
      //recover halfedges inserted that are on the intersection
      typedef std::pair<Node_id,Node_id> Node_id_pair;
      for(const Node_id_pair& node_id_pair : constrained_edges)
      {
        typename std::map<Node_id_pair,halfedge_descriptor>
          ::iterator it_poly_hedge=edge_to_hedge.find(node_id_pair);
        //we cannot have an assertion here in case an edge or part of an edge is a constraints.
        //Indeed, the graph_of_constraints report an edge 0,1 and 1,0 for example while only one of the two
        //is defined as one of them defines an adjacent face
        //CGAL_assertion(it_poly_hedge!=edge_to_hedge.end());
        if( it_poly_hedge!=edge_to_hedge.end() ){
          call_put(marks_on_edges,tm,edge(it_poly_hedge->second,tm),true);
          output_builder.set_edge_per_polyline(tm,node_id_pair,it_poly_hedge->second);
        }
        else{
          //WARNING: in few case this is needed if the marked edge is on the border
          //to optimize it might be better to only use sorted pair. TAG_SLXX1
          Node_id_pair opposite_pair(node_id_pair.second,node_id_pair.first);
          it_poly_hedge=edge_to_hedge.find(opposite_pair);
          CGAL_assertion( it_poly_hedge!=edge_to_hedge.end() );

          call_put(marks_on_edges,tm,edge(it_poly_hedge->second,tm),true);
          output_builder.set_edge_per_polyline(tm,opposite_pair,it_poly_hedge->second);
        }
      }
    }
  }

  void check_no_duplicates(const INodes& nodes) const
  {
    if (const_mesh_ptr == nullptr) // actually only needed for clip
      nodes.check_no_duplicates();
  }

  void finalize(INodes& nodes,
                const TriangleMesh& tm1,
                const TriangleMesh& tm2,
                const VertexPointMap1& vpm1,
                const VertexPointMap2& vpm2)
  {
    copy_nodes_ids_for_non_manifold_features();

    nodes.all_nodes_created();

    TriangleMesh* tm1_ptr = const_cast<TriangleMesh*>(&tm1);
    TriangleMesh* tm2_ptr = const_cast<TriangleMesh*>(&tm2);

    const Node_id nb_nodes = nodes.size();
    // we reserve nb_nodes+3 because we use the last three entries for the
    // face triangulation
    mesh_to_node_id_to_vertex[tm1_ptr].resize(nb_nodes+3);
    mesh_to_node_id_to_vertex[tm2_ptr].resize(nb_nodes+3);

    //store for each triangle face which boundary is intersected by the other surface,
    //original vertices (and halfedges in the refined mesh pointing on these vertices)
    std::map<TriangleMesh*,Face_boundaries> mesh_to_face_boundaries;

    //0) For each triangle mesh, collect original vertices that belongs to the intersection.
    //   From the graph of constraints, extract intersection edges that are incident to such vertices. In case
    //   there exists another original vertex adjacent to the first one found, this halfedge must be
    //   marked on the boundary (and possibly update an_edge_per_polyline).
    //   This is done first to avoid halfedges stored to be modified in the steps following.
    for (typename Mesh_to_vertices_on_intersection_map::iterator
          it=mesh_to_vertices_on_inter.begin();
          it!=mesh_to_vertices_on_inter.end();
          ++it)
    {
      TriangleMesh& tm=*it->first;
      CGAL_assertion(&tm!=const_mesh_ptr);

    //   Face_boundaries& face_boundaries=mesh_to_face_boundaries[&tm];

      Node_to_target_of_hedge_map& nodes_to_hedge=it->second;

      // iterate on the vertices that are on the intersection between the input meshes
      for(typename Node_to_target_of_hedge_map::iterator
            it_node_2_hedge=nodes_to_hedge.begin();
            it_node_2_hedge!=nodes_to_hedge.end();
            ++it_node_2_hedge)
      {
        Node_id node_id_of_first=it_node_2_hedge->first;
        // look for neighbors of the current node in the intersection graph
        std::vector<Node_id>& neighbors=graph_of_constraints[node_id_of_first];
        if ( !neighbors.empty() )
        {
          // for all neighbors look for input vertices that are also on the intersection
          for(Node_id node_id : neighbors)
          {
            //if already done for the opposite
            if (node_id >= node_id_of_first) continue;

            typename Node_to_target_of_hedge_map::iterator
              it_node_2_hedge_two = nodes_to_hedge.find(node_id);
            if ( it_node_2_hedge_two!=nodes_to_hedge.end() ) //a full edge is on intersection
            {
              //get the corresponding halfedge with vertex corresponding to node_id_of_first
              halfedge_descriptor hedge=it_node_2_hedge->second;
              halfedge_descriptor start=hedge;
              bool did_break=false;
              while ( source(hedge,tm) !=
                      target(it_node_2_hedge_two->second,tm) )
              {
                hedge=opposite(next(hedge,tm),tm);
                if ((doing_autorefinement || handle_non_manifold_features) && hedge==start)
                {
                  ++it_node_2_hedge_two; // we are using a multimap and
                                         // the halfedge we are looking for
                                         // might be on another sheet
                  if (it_node_2_hedge_two==nodes_to_hedge.end() ||
                      node_id!=it_node_2_hedge_two->first)
                  {
                    did_break=true;
                    break;
                  }
                  CGAL_assertion(it_node_2_hedge_two!=nodes_to_hedge.end());
                  CGAL_assertion(it_node_2_hedge->first==node_id_of_first);
                }
                else
                {
                  CGAL_assertion(hedge!=start);
                }
              }
              if (did_break) continue;
              std::pair<Node_id,Node_id> edge_pair(node_id,node_id_of_first);
              call_put(marks_on_edges,tm,edge(hedge,tm),true);
              output_builder.set_edge_per_polyline(tm,edge_pair,hedge);
            }
          }
        }
        #ifdef CGAL_COREFINEMENT_DEBUG
        else
        {
          std::cout << "X1: Found an isolated point" << std::endl;
        }
        #endif
      }
    }

    //1) First split halfedges cut by the intersection polyline(s)
    for (typename std::map<TriangleMesh*,On_edge_map>::iterator
      it=on_edge.begin(); it!=on_edge.end(); ++it)
    {
      if(it->first == tm1_ptr)
        split_halfedges(it, vpm1, nodes, mesh_to_face_boundaries);
      else
        split_halfedges(it, vpm2, nodes, mesh_to_face_boundaries);
    }

    //2)triangulation of the triangle faces containing intersection point in their interior
    //  and also those with intersection points only on the boundary.
    std::size_t total_size = 0;
    for (typename std::map<TriangleMesh*,On_face_map>::iterator
           it=on_face.begin(); it!=on_face.end(); ++it)
    {
      total_size += it->second.size();
    }

    user_visitor.start_triangulating_faces(total_size);

    for (typename std::map<TriangleMesh*,On_face_map>::iterator
           it=on_face.begin(); it!=on_face.end(); ++it)
    {
      if(it->first == tm1_ptr)
        triangulate_intersected_faces(it, vpm1, nodes, mesh_to_face_boundaries);
      else
        triangulate_intersected_faces(it, vpm2, nodes, mesh_to_face_boundaries);
    }

    user_visitor.end_triangulating_faces();

    nodes.finalize(mesh_to_node_id_to_vertex);
    user_visitor.start_building_output();
    // additional operations
    output_builder(nodes,
                   input_with_coplanar_faces,
                   is_node_of_degree_one,
                   mesh_to_node_id_to_vertex);

    user_visitor.end_building_output();
  }
};

} } } // CGAL::Polygon_mesh_processing::Corefinement

#include <CGAL/enable_warnings.h>

#endif //CGAL_POLYGON_MESH_PROCESSING_INTERNAL_COREFINEMENT_VISITOR_H
