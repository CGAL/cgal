// Copyright (c) 2016 GeometryFactory (France).
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
// Author(s)     : Sebastien Loriot

#ifndef CGAL_POLYGON_MESH_PROCESSING_INTERNAL_COREFINEMENT_VISITOR_H
#define CGAL_POLYGON_MESH_PROCESSING_INTERNAL_COREFINEMENT_VISITOR_H

#include <CGAL/license/Polygon_mesh_processing/corefinement.h>


#include <CGAL/Polygon_mesh_processing/internal/Corefinement/predicates.h>
#include <CGAL/Polygon_mesh_processing/internal/Corefinement/face_graph_utils.h>
#include <CGAL/utility.h>
#include <CGAL/Default.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_2_projection_traits_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

namespace CGAL{
namespace Corefinement{

// TODO option to ignore internal edges for patches of coplanar faces

template <class TriangleMesh>
struct Default_node_visitor{
  typedef boost::graph_traits<TriangleMesh> GT;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;
  typedef typename GT::vertex_descriptor vertex_descriptor;

  void new_node_added(  std::size_t /* node_id */,
                        Intersection_type /* type */,
                        halfedge_descriptor /* principal_edge */,
                        halfedge_descriptor /* additional_edge */,
                        bool /* is_target_coplanar */,
                        bool /* is_source_coplanar */ )
  {}

  void new_vertex_added(std::size_t /* node_id */,
                        vertex_descriptor /* vh */,
                        TriangleMesh& /*tm*/){}
};

template <class TriangleMesh>
struct Default_face_visitor{
  typedef boost::graph_traits<TriangleMesh> GT;
  typedef typename GT::face_descriptor face_descriptor;

  void before_subface_creations(face_descriptor /*f_old*/,TriangleMesh&){}
  void after_subface_created(face_descriptor /*f_new*/,TriangleMesh&){}
};

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
  template <class Mesh_to_intersection_edge_map,
            class Node_vector,
            class An_edge_per_polyline_map,
            class Mesh_to_map_node>
  void operator()(
    const std::map<const G*,Mesh_to_intersection_edge_map>& /*mesh_to_intersection_edges*/,
    const Node_vector& /*nodes*/,
    const An_edge_per_polyline_map& /*an_edge_per_polyline*/,
    bool /*input_have_coplanar_faces*/,
    const boost::dynamic_bitset<>& /* is_node_of_degree_one */,
    const Mesh_to_map_node& /*mesh_to_node_id_to_vertex*/) const
  {}
};

// A visitor for Intersection_of_triangle_meshes that can be used to corefine
// two meshes
template< class TriangleMesh,
          class VertexPointMap,
          class OutputBuilder_ = Default,
          class EdgeMarkMapBind_ = Default,
          class NewNodeVisitor_ = Default,
          class NewFaceVisitor_ = Default >
class Visitor{
//default template parameters
  typedef typename Default::Get<EdgeMarkMapBind_,
    Ecm_bind<TriangleMesh, No_mark<TriangleMesh> > >::type      EdgeMarkMapBind;
  typedef typename Default::Get<OutputBuilder_,
    No_extra_output_from_corefinement<TriangleMesh> >::type       OutputBuilder;
  typedef typename Default::Get<
    NewNodeVisitor_, Default_node_visitor<TriangleMesh> >::type  NewNodeVisitor;
  typedef typename Default::Get<
    NewFaceVisitor_, Default_face_visitor<TriangleMesh> >::type  NewFaceVisitor;

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
  // to maintain a halfedge on each polyline per TriangleMesh + pair<bool,size_t>
  // with first = "is the key (pair<Node_id,Node_id>) was reversed?" and
  // second is the number of edges +1 in the polyline
  typedef std::map< std::pair<Node_id,Node_id>,
                    std::pair< std::map<TriangleMesh*,halfedge_descriptor>,
                               std::pair<bool,std::size_t> > >
                                                       An_edge_per_polyline_map;

   typedef std::vector<Node_id>                                        Node_ids;
   typedef boost::unordered_map<face_descriptor,Node_ids>           On_face_map;
   typedef boost::unordered_map<edge_descriptor,Node_ids>           On_edge_map;
   //to keep the correspondance between node_id and vertex_handle in each mesh
   typedef std::vector<vertex_descriptor>                     Node_id_to_vertex;
   typedef std::map<TriangleMesh*, Node_id_to_vertex >         Mesh_to_map_node;
   //to handle coplanar halfedge of polyhedra that are full in the intersection
   typedef std::map<Node_id,halfedge_descriptor>    Node_to_target_of_hedge_map;
   typedef std::map<TriangleMesh*,Node_to_target_of_hedge_map>
                                           Mesh_to_vertices_on_intersection_map;
   typedef boost::unordered_map<vertex_descriptor,Node_id>    Vertex_to_node_id;
   typedef std::map<TriangleMesh*, Vertex_to_node_id> Mesh_to_vertex_to_node_id;
// typedef for the CDT
   typedef typename Intersection_nodes<TriangleMesh,
        VertexPointMap, Predicates_on_constructions_needed>::Exact_kernel    EK;
    typedef Triangulation_2_projection_traits_3<EK>                  CDT_traits;
    typedef Triangulation_vertex_base_with_info_2<Node_id,CDT_traits>        Vb;
    typedef Constrained_triangulation_face_base_2<CDT_traits>                Fb;
    typedef Triangulation_data_structure_2<Vb,Fb>                         TDS_2;
    typedef Constrained_Delaunay_triangulation_2<CDT_traits,TDS_2>          CDT;
    typedef typename CDT::Vertex_handle                       CDT_Vertex_handle;
// data members
private:
  // boost::dynamic_bitset<> non_manifold_nodes;
  std::vector< std::vector<Node_id> > graph_of_constraints;
  boost::dynamic_bitset<> is_node_of_degree_one;
  An_edge_per_polyline_map an_edge_per_polyline;
  //nb of intersection points between coplanar faces, see fixes XSL_TAG_CPL_VERT
  std::size_t number_coplanar_vertices;
  typename An_edge_per_polyline_map::iterator last_polyline;
  std::map<TriangleMesh*,On_face_map> on_face;
  std::map<TriangleMesh*,On_edge_map> on_edge;
  Mesh_to_vertices_on_intersection_map mesh_to_vertices_on_inter;
  Mesh_to_map_node mesh_to_node_id_to_vertex;
  // map an input vertex to a node id (input vertex on the intersection)
  Mesh_to_vertex_to_node_id mesh_to_vertex_to_node_id;

  std::map< Node_id,std::set<Node_id> > coplanar_constraints;

//data members that require initialization in the constructor
  NewNodeVisitor& new_node_visitor;
  NewFaceVisitor& new_face_visitor;
  OutputBuilder& output_builder;
  const EdgeMarkMapBind& marks_on_edges;
  bool input_with_coplanar_faces;
// visitor public functions
public:
  Visitor(NewNodeVisitor& v, NewFaceVisitor& f,
          OutputBuilder& o, const EdgeMarkMapBind& emm)
    : new_node_visitor(v)
    , new_face_visitor(f)
    , output_builder(o)
    , marks_on_edges(emm)
    , input_with_coplanar_faces(false)
  {}

  template<class Graph_node>
  void annotate_graph(std::vector<Graph_node>& graph)
  {
    std::size_t nb_nodes=graph.size();
    graph_of_constraints.resize(nb_nodes);
    is_node_of_degree_one.resize(nb_nodes);
    for(std::size_t node_id=0;node_id<nb_nodes;++node_id)
    {
    //   if (non_manifold_nodes.test(node_id))
    //     graph[node_id].make_terminal();
      graph_of_constraints[node_id].assign(
        graph[node_id].neighbors.begin(),
        graph[node_id].neighbors.end());

      if (graph_of_constraints[node_id].size()==1)
        is_node_of_degree_one.set(node_id);
    }
  }

  void start_new_polyline(std::size_t i, std::size_t j)
  {
    if ( i==j ) //case of a single point
    {
      // TODO shall we insert the point?????
      //TAG SL001
      //nothing is done
      return;
    }
    std::pair<typename An_edge_per_polyline_map::iterator,bool> res=
      an_edge_per_polyline.insert(
        std::make_pair( make_sorted_pair(i,j),
          std::make_pair( std::map<TriangleMesh*,halfedge_descriptor>(),std::make_pair(false,0))  )
      );
    CGAL_assertion(res.second);
    last_polyline=res.first;
    if ( i !=last_polyline->first.first )
      last_polyline->second.second.first=true;
  }

  void add_node_to_polyline(std::size_t)
  {
    ++(last_polyline->second.second.second);
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
    CGAL_assertion(!"This function should not be called");
  }

// The following code was used to split polylines at certains points.
// Here manifold was used in the context of the combinatorial map where
// the import inside the combinatorial map was possible (see broken_bound-[12].off)
// I keep the code here for now as it could be use to detect non-manifold
// situation of surfaces
// void check_node_on_non_manifold_edge(
//     std::size_t node_id,
//     halfedge_descriptor h,
//     const TriangleMesh& tm)
// {
//   if ( is_border_edge(h,tm) )
//    non_manifold_nodes.set(node_id);
// }
//
// void check_node_on_non_manifold_vertex(
//   std::size_t node_id,
//   halfedge_descriptor h,
//   const TriangleMesh& tm)
// {
//   //we turn around the hedge and check no halfedge is a border halfedge
//   BOOST_FOREACH(halfedge_descriptor hc,halfedges_around_target(h,tm))
//     if ( is_border_edge(hc,tm) )
//     {
//       non_manifold_nodes.set(node_id);
//       return;
//     }
// }

  //keep track of the fact that a polyhedron original vertex is a node
  void all_incident_faces_got_a_node_as_vertex(
      halfedge_descriptor h,
      Node_id node_id,
      TriangleMesh& tm)
  {
    mesh_to_vertex_to_node_id[&tm].insert(std::make_pair(target(h,tm),node_id));
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
    // non_manifold_nodes.resize(node_id+1);

    TriangleMesh* tm1_ptr = const_cast<TriangleMesh*>(&tm1);
    TriangleMesh* tm2_ptr = const_cast<TriangleMesh*>(&tm2);

    //forward to the visitor
    new_node_visitor.new_node_added(node_id, type, h_1, h_2, is_target_coplanar, is_source_coplanar);
    switch(type)
    {
      case ON_FACE: //Face intersected by an edge
        on_face[tm2_ptr][face(h_2,tm2)].push_back(node_id);
      break;
      case ON_EDGE: //Edge intersected by an edge
      {
        on_edge[tm2_ptr][edge(h_2,tm2)].push_back(node_id);
      //   check_node_on_non_manifold_edge(node_id,h_2,tm2);
      }
      break;
      case ON_VERTEX:
      {
        //grab original vertex that is on commom intersection
        mesh_to_vertices_on_inter[tm2_ptr].insert(std::make_pair(node_id,h_2));
        Node_id_to_vertex& node_id_to_vertex=mesh_to_node_id_to_vertex[tm2_ptr];
        if (node_id_to_vertex.size()<=node_id)
          node_id_to_vertex.resize(node_id+1,Graph_traits::null_vertex());
        node_id_to_vertex[node_id]=target(h_2,tm2);
        all_incident_faces_got_a_node_as_vertex(h_2,node_id,*tm2_ptr);
      //   check_node_on_non_manifold_vertex(node_id,h_2,tm2);
      }
      break;
      default:
      return;
    }

    CGAL_assertion(!is_target_coplanar || !is_source_coplanar); //coplanar edge are not forwarded

    if ( is_target_coplanar )
    {
      //grab original vertex that is on commom intersection
      mesh_to_vertices_on_inter[tm1_ptr].insert(std::make_pair(node_id,h_1));
      Node_id_to_vertex& node_id_to_vertex=mesh_to_node_id_to_vertex[tm1_ptr];
      if (node_id_to_vertex.size()<=node_id)
        node_id_to_vertex.resize(node_id+1,Graph_traits::null_vertex());
      node_id_to_vertex[node_id]=target(h_1,tm1);
      all_incident_faces_got_a_node_as_vertex(h_1,node_id, *tm1_ptr);
      // check_node_on_non_manifold_vertex(node_id,h_1,tm1);
    }
    else{
      if ( is_source_coplanar ){
        //grab original vertex that is on commom intersection
        halfedge_descriptor h_1_opp=opposite(h_1,tm1);
        mesh_to_vertices_on_inter[tm1_ptr].insert(std::make_pair(node_id,h_1_opp));
        Node_id_to_vertex& node_id_to_vertex=mesh_to_node_id_to_vertex[tm1_ptr];
        if(node_id_to_vertex.size()<=node_id)
          node_id_to_vertex.resize(node_id+1,Graph_traits::null_vertex());
        node_id_to_vertex[node_id]=source(h_1,tm1);
        all_incident_faces_got_a_node_as_vertex(h_1_opp,node_id, *tm1_ptr);
      //   check_node_on_non_manifold_vertex(node_id,h_1_opp,tm1);
      }
      else{
        //handle intersection on principal edge
        on_edge[tm1_ptr][edge(h_1,tm1)].push_back(node_id);
      //   check_node_on_non_manifold_edge(node_id,h_1,tm1);
      }
    }
  }

  //sort node ids so that we can split the hedge
  //consecutively
  template <class Node_vector>
  void sort_vertices_along_hedge(std::vector<std::size_t>& node_ids,
                                 halfedge_descriptor hedge,
                                 const TriangleMesh& tm,
                                 const VertexPointMap& vpm,
                                 const Node_vector& nodes)
  {
    std::sort(node_ids.begin(),
              node_ids.end(),
              Less_along_a_halfedge<TriangleMesh,VertexPointMap,Node_vector>
                (hedge, tm, vpm, nodes)
    );
  }

  void set_edge_per_polyline(TriangleMesh& tm,
                             std::pair<std::size_t,std::size_t> indices,
                             halfedge_descriptor hedge)
  {
    if (indices.first>indices.second)
    {
      std::swap(indices.first,indices.second);
      hedge=opposite(hedge,tm);
    }
    typename An_edge_per_polyline_map::iterator it =
      an_edge_per_polyline.find(indices);

    if (it!=an_edge_per_polyline.end()){
      CGAL_assertion(it->second.first.count(&tm) == 0 ||
                     it->second.first[&tm]==hedge);
      it->second.first.insert( std::make_pair( &tm,hedge) );
    }
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
  };

  //update the id of input mesh vertex that are also a node
  void update_face_indices(
    cpp11::array<vertex_descriptor,3>& f_vertices,
    cpp11::array<Node_id,3>& f_indices,
    Vertex_to_node_id& vertex_to_node_id)
  {
    for (int k=0;k<3;++k){
      typename boost::unordered_map<vertex_descriptor,Node_id>::iterator it =
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
      std::map< Node_id,std::set<Node_id> >::iterator it_neighbors =
        coplanar_constraints.find(node_id);
      if (it_neighbors!=coplanar_constraints.end())
      {
        CDT_Vertex_handle vh=id_to_CDT_vh[node_id];
        BOOST_FOREACH(Node_id id,it_neighbors->second)
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
    BOOST_FOREACH(Node_id id, node_ids)
    {
      CGAL_assertion(id < graph_of_constraints.size());
      std::vector<Node_id>& neighbors=graph_of_constraints[id];
      if (!neighbors.empty())
      {
        CDT_Vertex_handle vh=id_to_CDT_vh.find(id)->second;
        BOOST_FOREACH(Node_id id_n,neighbors)
        {
        //   if (id_n < id) continue; //no need to do it twice
          typename std::map<Node_id,CDT_Vertex_handle>
            ::iterator it_vh=id_to_CDT_vh.find(id_n);
          // this condition ensures to consider only graph edges that are in
          // the same triangle
          if ( !points_on_triangle || it_vh!=id_to_CDT_vh.end() ){
            CGAL_assertion(it_vh!=id_to_CDT_vh.end());
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

  void finalize(Intersection_nodes<TriangleMesh,
                 VertexPointMap, Predicates_on_constructions_needed>& nodes,
                 const TriangleMesh& tm1,
                 const TriangleMesh& tm2,
                 const VertexPointMap& vpm1,
                 const VertexPointMap& vpm2)
  {
    nodes.all_nodes_created();

    TriangleMesh* tm1_ptr = const_cast<TriangleMesh*>(&tm1);
    TriangleMesh* tm2_ptr = const_cast<TriangleMesh*>(&tm2);

    std::map<TriangleMesh*, VertexPointMap> vpms;
    vpms[tm1_ptr] = vpm1;
    vpms[tm2_ptr] = vpm2;

    vertex_descriptor null_vertex = Graph_traits::null_vertex();
    const Node_id nb_nodes = nodes.size();
    // we reserve nb_nodes+3 because we use the last three entries for the
    // face triangulation
    mesh_to_node_id_to_vertex[tm1_ptr].resize(nb_nodes+3, null_vertex);
    mesh_to_node_id_to_vertex[tm2_ptr].resize(nb_nodes+3, null_vertex);

    ///TODO check these comments
    //mark halfedge that are on the intersection
    //SL: I needed to use a map because to get the orientation around the edge,
    //    I need to know in the case the third vertex is a node its index (for exact construction)
    typedef boost::unordered_map<edge_descriptor,
                                 std::pair<Node_id,Node_id> > Intersection_edge_map;
    std::map<const TriangleMesh*,Intersection_edge_map> mesh_to_intersection_edges;

    //store for each triangle face which boundary is intersected by the other surface,
    //original vertices (and halfedges in the refined mesh pointing on these vertices)
    typedef boost::unordered_map<face_descriptor,Face_boundary> Face_boundaries;
    std::map<TriangleMesh*,Face_boundaries> mesh_to_face_boundaries;

    //0) For each triangle mesh, collect original vertices that belongs to the intersection.
    //   From the graph of constaints, extract intersection edges that are incident to such vertices. In case
    //   there exists another original vertex adjacent to the first one found, this halfedge must be
    //   marked on the boundary (and possibly update an_edge_per_polyline).
    //   This is done first to avoid halfedges stored to be modified in the steps following.
    for (typename Mesh_to_vertices_on_intersection_map::iterator
          it=mesh_to_vertices_on_inter.begin();
          it!=mesh_to_vertices_on_inter.end();
          ++it)
    {
      TriangleMesh& tm=*it->first;
      Intersection_edge_map& intersection_edges = mesh_to_intersection_edges[&tm];
    //   Face_boundaries& face_boundaries=mesh_to_face_boundaries[&tm];

      std::set<std::pair<Node_id,Node_id> > already_done;
      Node_to_target_of_hedge_map& nodes_to_hedge=it->second;
      for(typename Node_to_target_of_hedge_map::iterator
            it_node_2_hedge=nodes_to_hedge.begin();
            it_node_2_hedge!=nodes_to_hedge.end();
            ++it_node_2_hedge)
      {
        Node_id node_id_of_first=it_node_2_hedge->first;
        std::vector<Node_id>& neighbors=graph_of_constraints[node_id_of_first];
        if ( !neighbors.empty() )
        {
          BOOST_FOREACH(Node_id node_id, neighbors)
          {
            //if already done for the opposite
            if ( !already_done.insert(
              make_sorted_pair(node_id,node_id_of_first)).second ) continue;

            typename Node_to_target_of_hedge_map::iterator
              it_node_2_hedge_two = nodes_to_hedge.find(node_id);
            if ( it_node_2_hedge_two!=nodes_to_hedge.end() ) //a full edge is on intersection
            {
              //get the corresponding halfedge with vertex corresponding to node_id_of_first
              halfedge_descriptor hedge=it_node_2_hedge->second;
              CGAL_assertion_code(halfedge_descriptor start=hedge;)
              while ( source(hedge,tm) !=
                      target(it_node_2_hedge_two->second,tm) )
              {
                hedge=opposite(next(hedge,tm),tm);
                CGAL_assertion(hedge!=start);
              }
              std::pair<Node_id,Node_id> edge_pair(node_id,node_id_of_first);
              if ( intersection_edges.insert( std::make_pair(edge(hedge,tm),edge_pair) ).second)
                marks_on_edges.call_put(tm,edge(hedge,tm),true);
              set_edge_per_polyline(tm,edge_pair,hedge);
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
      TriangleMesh& tm=*it->first;
      const VertexPointMap& vpm=vpms[&tm];
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
        bool hedge_is_marked = marks_on_edges.call_get(tm,edge(hedge,tm));
        //do split the edges
        CGAL_assertion_code(vertex_descriptor expected_src=source(hedge,tm));
        BOOST_FOREACH(std::size_t node_id, node_ids)
        {
          halfedge_descriptor hnew = Euler::split_edge(hedge, tm);
          CGAL_assertion(expected_src==source(hnew,tm));
          vertex_descriptor vnew=target(hnew,tm);
          new_node_visitor.new_vertex_added(node_id, vnew, tm);
          nodes.call_put(vpm, vnew, node_id, tm);

          node_id_to_vertex[node_id]=vnew;
          if (first){
            first=false;
            hedge_incident_to_src=next(opposite(hedge,tm),tm);
          }

          //update marker tags. If the edge was marked, then the resulting edges in the split must be marked
          if ( hedge_is_marked )
            marks_on_edges.call_put(tm,edge(hnew,tm),true);

          CGAL_assertion_code(expected_src=vnew);
        }

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

    //2)triangulation of the triangle faces containing intersection point in their interior
    //  and also those with intersection points only on the boundary.
    for (typename std::map<TriangleMesh*,On_face_map>::iterator
      it=on_face.begin(); it!=on_face.end(); ++it)
    {
      TriangleMesh& tm=*it->first;
      const VertexPointMap& vpm=vpms[&tm];
      On_face_map& on_face_map=it->second;
      Face_boundaries& face_boundaries=mesh_to_face_boundaries[&tm];
      Node_id_to_vertex& node_id_to_vertex=mesh_to_node_id_to_vertex[&tm];
      Vertex_to_node_id& vertex_to_node_id=mesh_to_vertex_to_node_id[&tm];
      Intersection_edge_map& intersection_edges = mesh_to_intersection_edges[&tm];

      for (typename On_face_map::iterator it=on_face_map.begin();
            it!=on_face_map.end();++it)
      {
        face_descriptor f = it->first; //the face to be triangulated
        Node_ids& node_ids  = it->second; // ids of node in the interior of f
        typename Face_boundaries::iterator it_fb=face_boundaries.find(f);

        std::map<Node_id,typename CDT::Vertex_handle> id_to_CDT_vh;

        //associate an edge of the triangulation to a halfedge in a given polyhedron
        std::map<std::pair<Node_id,Node_id>,halfedge_descriptor> edge_to_hedge;

        // the vertices of f
        cpp11::array<vertex_descriptor,3> f_vertices;
        // the node_id of an input vertex or a fake id (>=nb_nodes)
        cpp11::array<Node_id,3> f_indices = {{nb_nodes,nb_nodes+1,nb_nodes+2}};
        if (it_fb!=face_boundaries.end()){ //the boundary of the triangle face was refined
          f_vertices[0]=it_fb->second.vertices[0];
          f_vertices[1]=it_fb->second.vertices[1];
          f_vertices[2]=it_fb->second.vertices[2];
          update_face_indices(f_vertices,f_indices,vertex_to_node_id);
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

        typename EK::Point_3 p = nodes.to_exact(get(vpm,f_vertices[0])),
                             q = nodes.to_exact(get(vpm,f_vertices[1])),
                             r = nodes.to_exact(get(vpm,f_vertices[2]));

        CDT_traits traits(typename EK::Construct_normal_3()(p,q,r));
        CDT cdt(traits);

        // insert triangle points
        cpp11::array<CDT_Vertex_handle,3> triangle_vertices;
        //we can do this to_exact because these are supposed to be input points.
        triangle_vertices[0]=cdt.insert_outside_affine_hull(p);
        triangle_vertices[1]=cdt.insert_outside_affine_hull(q);
        triangle_vertices[2]=cdt.tds().insert_dim_up(cdt.infinite_vertex(), false);
        triangle_vertices[2]->set_point(r);


        triangle_vertices[0]->info()=f_indices[0];
        triangle_vertices[1]->info()=f_indices[1];
        triangle_vertices[2]->info()=f_indices[2];

        node_id_to_vertex[nb_nodes  ]=f_vertices[0];
        node_id_to_vertex[nb_nodes+1]=f_vertices[1];
        node_id_to_vertex[nb_nodes+2]=f_vertices[2];

        //if one of the triangle input vertex is also a node
        for (int ik=0;ik<3;++ik){
          if ( f_indices[ik]<nb_nodes )
            id_to_CDT_vh.insert(
                std::make_pair(f_indices[ik],triangle_vertices[ik]));
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
              BOOST_FOREACH(Node_id id, ids_on_edge)
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
        BOOST_FOREACH(Node_id node_id, node_ids)
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
            std::map< Node_id,std::set<Node_id> >::iterator res =
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
          edge_to_hedge, cdt, vpm, new_node_visitor, new_face_visitor);

        // TODO Here we do the update only for internal edges.
        // Update for border halfedges could be done during the split

        //3) mark halfedges that are common to two polyhedral surfaces
        //recover halfedges inserted that are on the intersection
        typedef std::pair<Node_id,Node_id> Node_id_pair;
        BOOST_FOREACH(const Node_id_pair& node_id_pair, constrained_edges)
        {
          typename std::map<Node_id_pair,halfedge_descriptor>
            ::iterator it_poly_hedge=edge_to_hedge.find(node_id_pair);
          //we cannot have an assertion here in case an edge or part of an edge is a constraints.
          //Indeed, the graph_of_constraints report an edge 0,1 and 1,0 for example while only one of the two
          //is defined as one of them defines an adjacent face
          //CGAL_assertion(it_poly_hedge!=edge_to_hedge.end());
          if( it_poly_hedge!=edge_to_hedge.end() ){
            if ( intersection_edges.insert(
                    std::make_pair(edge(it_poly_hedge->second,tm),node_id_pair) ).second)
              marks_on_edges.call_put(tm,edge(it_poly_hedge->second,tm),true);
            set_edge_per_polyline(tm,node_id_pair,it_poly_hedge->second);
          }
          else{
            //WARNING: in few case this is needed if the marked edge is on the border
            //to optimize it might be better to only use sorted pair. TAG_SLXX1
            Node_id_pair opposite_pair(node_id_pair.second,node_id_pair.first);
            it_poly_hedge=edge_to_hedge.find(opposite_pair);
            CGAL_assertion( it_poly_hedge!=edge_to_hedge.end() );

            if ( intersection_edges.insert(
                std::make_pair(edge(it_poly_hedge->second,tm),opposite_pair) ).second )
              marks_on_edges.call_put(tm,edge(it_poly_hedge->second,tm),true);
            set_edge_per_polyline(tm,opposite_pair,it_poly_hedge->second);
          }
        }
      }
    }

    nodes.finalize();

    // additional operations
    output_builder(mesh_to_intersection_edges,
                   nodes,
                   an_edge_per_polyline,
                   input_with_coplanar_faces,
                   is_node_of_degree_one,
                   mesh_to_node_id_to_vertex);
  }
};

} } //end of namespace CGAL::Corefinement

#endif //CGAL_POLYGON_MESH_PROCESSING_INTERNAL_COREFINEMENT_VISITOR_H
