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

#ifndef CGAL_POLYGON_MESH_PROCESSING_INTERNAL_FACE_GRAPH_UTILS_H
#define CGAL_POLYGON_MESH_PROCESSING_INTERNAL_FACE_GRAPH_UTILS_H

#include <CGAL/license/Polygon_mesh_processing/corefinement.h>

#include <CGAL/AABB_trees/intersection.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>

#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/property_map.h>
#include <fstream>
#include <sstream>
#include <set>
#include <type_traits>

namespace CGAL {
namespace Polygon_mesh_processing {
namespace Corefinement {

template<class T>
struct is_surface_mesh : std::false_type {};

enum Boolean_operation_type {UNION = 0, INTERSECTION,
                             TM1_MINUS_TM2, TM2_MINUS_TM1, NONE };

template <typename G>
struct No_mark
{
  friend
  bool get(No_mark<G>,
                  typename boost::graph_traits<G>::edge_descriptor)
  {
    return false;
  }
  friend void put(No_mark<G>,
                  typename boost::graph_traits<G>::edge_descriptor, bool)
  {}
};

template<class G,
         class EdgeMarkMap>
void mark_all_edges(G& tm, EdgeMarkMap& edge_mark_map)
{
  for(typename boost::graph_traits<G>::edge_descriptor ed :
                edges(tm))
  {
    put(edge_mark_map, ed, true);
  }
}

template<class G>
void mark_all_edges(G&, No_mark<G>&)
{} //nothing to do

template<class G,
         class EdgeMarkMap,
         class HalfedgeRange>
void unmark_edges(      G& tm,
                        EdgeMarkMap& edge_mark_map,
                  const HalfedgeRange& hedges)
{
  for(typename boost::graph_traits<G>::halfedge_descriptor hd : hedges)
    put(edge_mark_map, edge(hd, tm), false);
}

template <class G,
          class HalfedgeRange>
void unmark_edges(      G&,
                        No_mark<G>&,
                  const HalfedgeRange&)
{} //nothing to do

template <class G,
          class edge_descriptor,
          class EdgeMarkMapIn,
          class EdgeMarkMapOut>
void copy_edge_mark(edge_descriptor ed_in,
                    edge_descriptor ed_out,
                    const EdgeMarkMapIn& edge_mark_map_in,
                          EdgeMarkMapOut& edge_mark_map_out)
{
  if(get(edge_mark_map_in, ed_in))
    put(edge_mark_map_out, ed_out, true);
}

template <class G,
          class edge_descriptor,
          class EdgeMarkMapOut>
void copy_edge_mark(edge_descriptor,
                    edge_descriptor,
                    No_mark<G>,
                    EdgeMarkMapOut)
{} // nothing to do

template <class G,
          class edge_descriptor,
          class EdgeMarkMapIn>
void copy_edge_mark(edge_descriptor,
                    edge_descriptor,
                    EdgeMarkMapIn,
                    No_mark<G>)
{} // nothing to do

template <class G,
          class edge_descriptor>
void copy_edge_mark(edge_descriptor,
                    edge_descriptor,
                    No_mark<G>,
                    No_mark<G>)
{} // nothing to do

template <class G,
          class EdgeMarkMapIn,
          class EdgeMarkMapOut>
void copy_edge_mark(G& g,
                    const EdgeMarkMapIn & edge_mark_map_in,
                          EdgeMarkMapOut& edge_mark_map_out)
{
  for(typename boost::graph_traits<G>::edge_descriptor ed : edges(g))
    if(get(edge_mark_map_in, ed))
      put(edge_mark_map_out, ed, true);
}

template <class G,
          class EdgeMarkMapOut>
void copy_edge_mark(G&,
                    const No_mark<G> &,
                          EdgeMarkMapOut&)
{} // nothing to do

template <class G,
          class EdgeMarkMapIn>
void copy_edge_mark(G&,
                    const EdgeMarkMapIn&,
                          No_mark<G>&)
{} // nothing to do

template <class G>
void copy_edge_mark(G&,
                    const No_mark<G>&,
                          No_mark<G>&)
{} // nothing to do


// For exact side_of_triangle_mesh
template <class Node_id_map,
          class VertexPointMap,
          class NodeVector>
struct Node_vector_exact_vertex_point_map
{
// map type definitions
  typedef typename boost::property_traits<VertexPointMap>::key_type key_type;
  typedef typename NodeVector::Exact_kernel Exact_kernel;
  typedef typename Exact_kernel::Point_3 value_type;
  typedef value_type reference;
  typedef boost::readable_property_map_tag category;
// internal type definitions
  typedef std::size_t Node_id;

  Node_vector_exact_vertex_point_map(){}

  Node_vector_exact_vertex_point_map(const Node_id_map& node_ids,
                                     const VertexPointMap& vpm,
                                     const NodeVector& node_vector)
    : node_ids(&node_ids)
    , vpm(&vpm)
    , node_vector(&node_vector)
  {}

  friend value_type get(Node_vector_exact_vertex_point_map m, key_type k)
  {
    typename Node_id_map::const_iterator it = m.node_ids->find(k);
    if (it == m.node_ids->end())
      return m.node_vector->to_exact( get(*(m.vpm), k) );
    return m.node_vector->exact_node(it->second);
  }

  const Node_id_map* node_ids;
  const VertexPointMap* vpm;
  const NodeVector* node_vector;
};

// For exact side_of_triangle_mesh
template <class TriangleMesh, class PPM, class TreeTraits>
struct Split_primitives
{
  Split_primitives(TriangleMesh& tm, PPM ppm)
    : tm(tm)
    , ppm(ppm)
  {}

  template<typename PrimitiveIterator>
  void operator()(PrimitiveIterator first,
                  PrimitiveIterator beyond,
                  const CGAL::Bbox_3& bbox) const
    {
      typedef typename std::iterator_traits<PrimitiveIterator>::value_type Prmtv;
      PrimitiveIterator middle = first + (beyond - first)/2;
      typedef typename std::iterator_traits<PrimitiveIterator>::value_type Prmtv;
      switch(TreeTraits::longest_axis(bbox))
      {
      case TreeTraits::CGAL_AXIS_X: // sort along x
        std::nth_element(first, middle, beyond, [this](const Prmtv& p1, const Prmtv& p2){ return get(ppm, p1.id()).x() < get(ppm,p2.id()).x(); });
        break;
      case TreeTraits::CGAL_AXIS_Y: // sort along y
        std::nth_element(first, middle, beyond, [this](const Prmtv& p1, const Prmtv& p2){ return get(ppm, p1.id()).y() < get(ppm,p2.id()).y(); });
        break;
      case TreeTraits::CGAL_AXIS_Z: // sort along z
        std::nth_element(first, middle, beyond, [this](const Prmtv& p1, const Prmtv& p2){ return get(ppm, p1.id()).z() < get(ppm,p2.id()).z(); });
        break;
      default:
        CGAL_error();
      }
    }
  TriangleMesh& tm;
  PPM ppm;
};

// For exact side_of_triangle_mesh
template <class BPM>
struct Compute_bbox {
  Compute_bbox(const BPM& bpm)
    : bpm(bpm)
  {}

  template<typename ConstPrimitiveIterator>
  CGAL::Bbox_3 operator()(ConstPrimitiveIterator first,
                          ConstPrimitiveIterator beyond) const
  {
    CGAL::Bbox_3 bbox = get(bpm, first->id());
    for(++first; first != beyond; ++first)
    {
      bbox += get(bpm, first->id());
    }
    return bbox;
  }
  BPM bpm;
};

// For exact side_of_triangle_mesh
template <class TriangleMesh,
          class Node_id_map,
          class VertexPointMap,
          class NodeVector,
          class Input_Kernel>
struct Side_of_helper
{
  typedef Node_vector_exact_vertex_point_map<Node_id_map, VertexPointMap, NodeVector> VPM;

  typedef CGAL::AABB_face_graph_triangle_primitive<TriangleMesh, VPM> Primitive;
  typedef CGAL::AABB_traits_3<typename NodeVector::Exact_kernel, Primitive> Traits;
  typedef CGAL::AABB_tree<Traits> Tree_type;

  static
  VPM get_vpm(const Node_id_map& node_ids,
               const VertexPointMap& vpm,
               const NodeVector& node_vector)
  {
    return VPM(node_ids, vpm, node_vector);
  }

  template <class FaceIdMap>
  static
  void build_tree(TriangleMesh& tm,
                  Tree_type& tree,
                  const Node_id_map& node_ids,
                  FaceIdMap fid,
                  const VertexPointMap& vpm,
                  const NodeVector& node_vector)
  {
    typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;
    typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;
    typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;

    // add primitives
    tree.insert(faces(tm).begin(), faces(tm).end(), tm, get_vpm(node_ids, vpm, node_vector));

    // pre-build bboxes (using approximation)
    std::vector<Bbox_3> face_bboxes(num_faces(tm));

    auto get_v_box = [&node_ids, &node_vector, &vpm](vertex_descriptor v)
    {
      typename Node_id_map::const_iterator it = node_ids.find(v);
      if (it == node_ids.end())
        return get(vpm, v).bbox();
      return approx(node_vector.exact_node(it->second)).bbox();
    };

    for (face_descriptor f : faces(tm))
    {
      halfedge_descriptor h = halfedge(f, tm);
      face_bboxes[get(fid, f)] = get_v_box( source(h, tm) ) +
                                 get_v_box( target(h, tm) ) +
                                 get_v_box( target(next(h, tm), tm) );
    }

    typedef CGAL::Pointer_property_map<CGAL::Bbox_3>::type Id_to_box;
    Id_to_box id_to_box = CGAL::make_property_map(face_bboxes);
    typedef Compose_property_map<FaceIdMap, Id_to_box> BPM;
    BPM bpm(fid, id_to_box);
    Compute_bbox<BPM> compute_bbox(bpm);

    typedef One_point_from_face_descriptor_map<TriangleMesh, VertexPointMap> PPM;
    PPM ppm(&tm, vpm);

    Split_primitives<TriangleMesh, PPM, Traits> split_primitives(tm, ppm);
    tree.custom_build(compute_bbox, split_primitives);
  }
};

template <class TriangleMesh,
          class Node_id_map,
          class VertexPointMap,
          class NodeVector>
struct Side_of_helper<TriangleMesh, Node_id_map, VertexPointMap, NodeVector, typename NodeVector::Exact_kernel>
{
  typedef VertexPointMap VPM;
  static
  VPM get_vpm(const Node_id_map&,
              const VertexPointMap& vpm,
              const NodeVector&)
  {
    return vpm;
  }

  typedef CGAL::AABB_face_graph_triangle_primitive<TriangleMesh, VPM> Primitive;
  typedef CGAL::AABB_traits_3<typename NodeVector::Exact_kernel, Primitive> Traits;
  typedef CGAL::AABB_tree<Traits> Tree_type;


  template <class FaceIdMap>
  static
  void build_tree(TriangleMesh& tm,
                  Tree_type& tree,
                  const Node_id_map& /* node_ids */,
                  FaceIdMap /* fid */,
                  const VertexPointMap& vpm,
                  const NodeVector& /* node_vector */)
  {
    tree.insert(faces(tm).begin(), faces(tm).end(), tm, vpm);
    tree.build();
  }
};

// Parts to get default property maps for output meshes based on the value type
// of input vertex point maps.
template <typename Point_3, typename vertex_descriptor>
struct Dummy_default_vertex_point_map
{
  typedef vertex_descriptor key_type;
  typedef Point_3 value_type;
  typedef Point_3 reference;
  typedef boost::read_write_property_map_tag category;

  inline friend
  value_type
  get(Dummy_default_vertex_point_map, key_type)
  {
    CGAL_assertion(false ||
      !"This property map should not be used."
       "Check the value type of the output vpm vs that of input");
    return Point_3();
  }

  inline friend
  void
  put(Dummy_default_vertex_point_map, key_type, value_type)
  {
    CGAL_assertion(false ||
      !"This property map should not be used."
       "Check the value type of the output vpm vs that of input");
  }
};

template <class Point_3, class NamedParameters, class PolygonMesh>
struct TweakedGetVertexPointMap
{
  typedef typename GetVertexPointMap<PolygonMesh,
                                     NamedParameters>::type Default_map;
  typedef typename std::is_same<Point_3,
    typename boost::property_traits<Default_map>::value_type>::type Use_default_tag;

  typedef std::conditional_t<
    Use_default_tag::value,
    Default_map,
    Dummy_default_vertex_point_map<Point_3,
      typename boost::graph_traits<PolygonMesh>::vertex_descriptor >
  > type;
};

template <class PT, class NP, class PM>
std::optional< typename TweakedGetVertexPointMap<PT, NP, PM>::type >
get_vpm(const NP& np, std::optional<PM*> opm, std::true_type)
{
  if (std::nullopt == opm) return std::nullopt;
  return parameters::choose_parameter(
           parameters::get_parameter(np, internal_np::vertex_point),
           get_property_map(boost::vertex_point, *(*opm)) );
}

template <class PT, class NP, class PM>
std::optional< typename TweakedGetVertexPointMap<PT, NP, PM>::type >
get_vpm(const NP&, std::optional<PM*> opm, std::false_type)
{
  if (std::nullopt == opm) return std::nullopt;
  return typename TweakedGetVertexPointMap<PT, NP, PM>::type();
}
//

template <class TriangleMesh>
struct Default_visitor{
  typedef boost::graph_traits<TriangleMesh> GT;
  typedef typename GT::face_descriptor face_descriptor;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;
  typedef typename GT::vertex_descriptor vertex_descriptor;
// face visitor functions
  void before_subface_creations(face_descriptor /*f_old*/, const TriangleMesh&){}
  void after_subface_creations(const TriangleMesh&){}
  void before_subface_created(const TriangleMesh&){}
  void after_subface_created(face_descriptor /*f_new*/,const TriangleMesh&){}
  void before_face_copy(face_descriptor /*f_old*/, const TriangleMesh&, const TriangleMesh&){}
  void after_face_copy(face_descriptor /*f_old*/, const TriangleMesh&,
                       face_descriptor /* f_new */, const TriangleMesh&){}
  void subface_of_coplanar_faces_intersection(face_descriptor /*f*/, const TriangleMesh&) {}
// edge visitor functions
  void before_edge_split(halfedge_descriptor /* h */, const TriangleMesh& /* tm */){}
  void edge_split(halfedge_descriptor /* hnew */, const TriangleMesh& /* tm */){}
  void after_edge_split(){}

  void add_retriangulation_edge(halfedge_descriptor /* h */ , const TriangleMesh& /* tm */) {} // edges added during split face retriangulation

  void before_edge_copy(halfedge_descriptor /*h_old*/, const TriangleMesh&, const TriangleMesh&){}
  void after_edge_copy(halfedge_descriptor /*h_old*/, const TriangleMesh&,
                       halfedge_descriptor /* f_new */, const TriangleMesh&){}

  void before_edge_duplicated(halfedge_descriptor /*h_old*/, const TriangleMesh&){} // called before a patch border edge is duplicated
  void after_edge_duplicated(halfedge_descriptor /*h_old*/,
                             halfedge_descriptor /* f_new */, const TriangleMesh&){} // called after a patch border edge is duplicated

  void intersection_edge_copy(halfedge_descriptor /* h_old1 */, const TriangleMesh& /* tm1 */,
                              halfedge_descriptor /* h_old2 */, const TriangleMesh& /* tm2 */,
                              halfedge_descriptor /* h_new */,  const TriangleMesh& /* tm_new */){}
// vertex visitor functions
  void new_vertex_added(std::size_t /* node_id */,
                        vertex_descriptor /* vh */,
                        const TriangleMesh& /* tm */){}

  void intersection_point_detected(std::size_t /* node_id */,
                                   int /* sdim */,
                                   halfedge_descriptor /* principal_edge */,
                                   halfedge_descriptor /* additional_edge */,
                                   const TriangleMesh& /* tm1 */,
                                   const TriangleMesh& /* tm2 */,
                                   bool /* is_target_coplanar */,
                                   bool /* is_source_coplanar */) {}
  void before_vertex_copy(vertex_descriptor /*v_src*/, const TriangleMesh& /*tm_src*/, const TriangleMesh& /*tm_tgt*/){}
  void after_vertex_copy(vertex_descriptor /*v_src*/, const TriangleMesh& /*tm_src*/,
                         vertex_descriptor /* v_tgt */, const TriangleMesh& /*tm_tgt*/){}

  // progress tracking
  void start_filtering_intersections() const {}
  void progress_filtering_intersections(double ) const {}
  void end_filtering_intersections() const {}

  void start_triangulating_faces(std::size_t) const {}
  void triangulating_faces_step() const {}
  void end_triangulating_faces() const {}

  void start_handling_intersection_of_coplanar_faces(std::size_t) const {}
  void intersection_of_coplanar_faces_step() const {}
  void end_handling_intersection_of_coplanar_faces() const {}

  void start_handling_edge_face_intersections(std::size_t) const {}
  void edge_face_intersections_step() const {}
  void end_handling_edge_face_intersections() const {}

  void start_building_output() const {}
  void end_building_output() const {}

// Required by Face_graph_output_builder
  void filter_coplanar_edges() const {}
  void detect_patches() const {}
  void classify_patches() const {}
  void classify_intersection_free_patches(const TriangleMesh&) const {}
  void out_of_place_operation(Boolean_operation_type) const {}
  void in_place_operation(Boolean_operation_type) const {}
  void in_place_operations(Boolean_operation_type,Boolean_operation_type) const {}

  // calls commented in the code and probably incomplete due to the migration
  // see NODE_VISITOR_TAG
  /*
  // autorefinement only
  void new_node_added_triple_face(std::size_t node_id,
                                  face_descriptor f1,
                                  face_descriptor f2,
                                  face_descriptor f3,
                                  const TriangleMesh& tm)
  {}
  */
};

template < class TriangleMesh,
           class VertexPointMap,
           class Node_id,
           class Node_vector,
           class Node_id_to_vertex,
           class CDT,
           class OutputBuilder,
           class UserVisitor>
void
triangulate_a_face(
  typename boost::graph_traits<TriangleMesh>::face_descriptor current_face,
  TriangleMesh& tm,
  Node_vector& nodes,
  const std::vector<Node_id>& node_ids,
  Node_id_to_vertex& node_id_to_vertex,
  std::map<std::pair<Node_id,Node_id>,
           typename boost::graph_traits<TriangleMesh>
                ::halfedge_descriptor>& edge_to_hedge,
  const CDT& cdt,
  const VertexPointMap& vpm,
  OutputBuilder& output_builder,
  UserVisitor& user_visitor)
{
  typedef boost::graph_traits<TriangleMesh> GT;
  typedef typename GT::vertex_descriptor vertex_descriptor;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;
  typedef typename GT::edge_descriptor edge_descriptor;

  //insert the intersection point interior to the face inside the mesh and
  //save their vertex_descriptor
  CGAL_assertion( node_ids.size()== std::set<Node_id>(node_ids.begin(), node_ids.end()).size() );
  for(Node_id node_id : node_ids)
  {
    vertex_descriptor v=add_vertex(tm);
    nodes.call_put(vpm, v, node_id, tm);
    // register the new vertex in the output builder
    output_builder.set_vertex_id(v, node_id, tm);
    user_visitor.new_vertex_added(node_id, v, tm);

    CGAL_assertion(node_id_to_vertex.size()>node_id);
    node_id_to_vertex.register_vertex(node_id, v);
  }

  //insert the new halfedge and set their incident vertex

  //grab edges that are not on the convex hull (these have already been created)
  for (typename CDT::Finite_edges_iterator it=cdt.finite_edges_begin();
                                           it!=cdt.finite_edges_end(); ++it)
  {
    typename CDT::Vertex_handle cdt_v0=it->first->vertex( cdt.ccw(it->second) );
    typename CDT::Vertex_handle cdt_v1=it->first->vertex( cdt.cw(it->second) );

    // consider edges not on the convex hull (not on the boundary of the face)
    // and create the corresponding halfedges
    if ( !cdt.is_infinite(it->first->vertex(it->second)) &&
         !cdt.is_infinite(cdt.mirror_vertex(it->first,it->second)) )
    {
      edge_descriptor e=add_edge(tm);
      user_visitor.add_retriangulation_edge(halfedge(e, tm), tm);
      halfedge_descriptor h=halfedge(e,tm), h_opp=opposite(h,tm);

      Node_id i0=cdt_v0->info(), i1=cdt_v1->info();
      CGAL_assertion( node_id_to_vertex.get_vertex(i0)!=GT::null_vertex());
      CGAL_assertion( node_id_to_vertex.get_vertex(i1)!=GT::null_vertex());

      vertex_descriptor v0=node_id_to_vertex.get_vertex(i0), v1=node_id_to_vertex.get_vertex(i1);

      set_target(h,v0,tm);
      set_target(h_opp,v1,tm);
      set_halfedge(v0,h,tm);
      set_halfedge(v1,h_opp,tm);

      edge_to_hedge[std::make_pair(i0,i1)]=h_opp;
      edge_to_hedge[std::make_pair(i1,i0)]=h;
    }
  }

  //grab triangles.
  user_visitor.before_subface_creations(current_face,tm);
  for (typename CDT::Finite_faces_iterator it=cdt.finite_faces_begin(),
                                           it_end=cdt.finite_faces_end();;)
  {
    typename CDT::Vertex_handle cdt_v0=it->vertex(0);
    typename CDT::Vertex_handle cdt_v1=it->vertex(1);
    typename CDT::Vertex_handle cdt_v2=it->vertex(2);

    Node_id i0=cdt_v0->info(), i1=cdt_v1->info(), i2=cdt_v2->info();

    CGAL_assertion(edge_to_hedge.count(std::make_pair(i0,i1))!= 0);
    CGAL_assertion(edge_to_hedge.count(std::make_pair(i1,i2))!= 0);
    CGAL_assertion(edge_to_hedge.count(std::make_pair(i2,i0))!= 0);

    halfedge_descriptor h01=edge_to_hedge[std::make_pair(i0,i1)];
    halfedge_descriptor h12=edge_to_hedge[std::make_pair(i1,i2)];
    halfedge_descriptor h20=edge_to_hedge[std::make_pair(i2,i0)];

    CGAL_assertion(target(h01,tm)==node_id_to_vertex.get_vertex(i1));
    CGAL_assertion(target(h12,tm)==node_id_to_vertex.get_vertex(i2));
    CGAL_assertion(target(h20,tm)==node_id_to_vertex.get_vertex(i0));

    set_next(h01,h12,tm);
    set_next(h12,h20,tm);
    set_next(h20,h01,tm);

    //update face halfedge
    set_halfedge(current_face,h01,tm);

    //update face of halfedges
    set_face(h01,current_face,tm);
    set_face(h12,current_face,tm);
    set_face(h20,current_face,tm);

    if ( ++it!=it_end )
    {
      user_visitor.before_subface_created(tm);
      current_face=add_face(tm);
      user_visitor.after_subface_created(current_face,tm);
    }
    else
      break;
  }
  user_visitor.after_subface_creations(tm);
}

template <class PolygonMesh>
class Border_edge_map {
  typedef boost::graph_traits<PolygonMesh> GT;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;
  typedef typename GT::edge_descriptor edge_descriptor;
  typedef std::unordered_set<edge_descriptor> Intersection_edge_map;
  const Intersection_edge_map* intersection_edges;
  const PolygonMesh* tm;
public:
  // required by the property forwarder of the filtered_graph in boost
  Border_edge_map()
   : intersection_edges(NULL)
   , tm(NULL)
  {}

  Border_edge_map(const Intersection_edge_map& intersection_edges,
                 const PolygonMesh& tm)
   : intersection_edges(&intersection_edges)
   , tm(&tm)
  {}

  friend
  bool get(const Border_edge_map& map, edge_descriptor e)
  {
    if ( is_border(e,*map.tm) ) return false;
    return map.intersection_edges->count(e)!=0;
  }
};

template<class PolygonMesh>
struct Patch_description{
  typedef boost::graph_traits<PolygonMesh> GT;
  typedef typename GT::face_descriptor face_descriptor;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;
  typedef typename GT::vertex_descriptor vertex_descriptor;

  std::vector<face_descriptor> faces;
  std::vector<vertex_descriptor> interior_vertices;
  std::vector<halfedge_descriptor> interior_edges;
  std::vector<halfedge_descriptor> shared_edges;

  // Three sets used only for patch with borders not shared, their are empty if not
  std::set<halfedge_descriptor> border_edges;
  std::set<halfedge_descriptor> border_with_shared_target;
  std::set<halfedge_descriptor> border_with_shared_source;

  bool is_initialized;

  Patch_description(): is_initialized(false) {};
};

// shared_edges will be filled by halfedges pointing in the patch
// that are inside `is_intersection_edge`, thus mesh boundary halfedges
// are not necessarily inside.
template <class PolygonMesh, class IsIntersectionEdge>
void extract_patch_simplices(
  PolygonMesh& pm,
  Patch_description<PolygonMesh>& patch,
  const IsIntersectionEdge& is_intersection_edge)
{
  typedef boost::graph_traits<PolygonMesh> GT;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;
  typedef typename GT::vertex_descriptor vertex_descriptor;
  typedef typename GT::face_descriptor face_descriptor;

  for(face_descriptor f : patch.faces)  {
    for(halfedge_descriptor h : halfedges_around_face(halfedge(f, pm),pm))    {
      if ( !is_intersection_edge.count(edge(h, pm)) )
      {
        halfedge_descriptor h_opp = opposite(h, pm);
        if(is_border(h_opp,pm)){
          patch.border_edges.insert(h_opp);
          patch.interior_edges.push_back(h);
        }
        else if ( h < h_opp)
          patch.interior_edges.push_back(h);
      }
      else
        patch.shared_edges.push_back(h);
    }
  }

  std::set<vertex_descriptor> border_vertices;
  for(halfedge_descriptor h : patch.shared_edges){
    border_vertices.insert( target(h,pm) );
    // if the model is not closed i.e. patch_border_halfedge is not cycle-only
    if(!patch.interior_edges.empty())
      border_vertices.insert( source(h,pm) );
  }

  // a patch is a disk v = 2 + e - f
  // v_in = v - v_border, e = e_in + e_shared, v_border = e_shared ( +1 if not closed)
  // v_in = 2 + e_in - f (-1 if not closed)
  patch.interior_vertices.reserve(2 + patch.interior_edges.size() - patch.faces.size() );
  for(halfedge_descriptor h : patch.interior_edges){
    if ( !border_vertices.count( target(h,pm) ) &&
         halfedge(target(h,pm),pm) == h ) // only add the vertex once
      patch.interior_vertices.push_back( target(h,pm) );
    if ( !border_vertices.count( source(h,pm) ) &&
         halfedge(source(h,pm),pm) == opposite(h,pm) ) // only add the vertex once
      patch.interior_vertices.push_back( source(h,pm) );
  }

  for(halfedge_descriptor h : patch.border_edges){
    if ( border_vertices.count( target(h,pm) ) )
      patch.border_with_shared_target.emplace(h);
    if ( border_vertices.count( source(h,pm) ) )
      patch.border_with_shared_source.emplace(h);
  }
}

template <class PolygonMesh, class FaceIndexMap, class IsIntersectionEdge>
struct Patch_container{
// data members
  std::vector< Patch_description<PolygonMesh> > patches;
// external data members
  PolygonMesh& pm;
  const std::vector<std::size_t>& patch_ids;
  const FaceIndexMap fids;
  const IsIntersectionEdge& is_intersection_edge;
// constructor
  Patch_container(
    PolygonMesh& pm,
    const std::vector<std::size_t>& patch_ids,
    const FaceIndexMap fids,
    const IsIntersectionEdge& is_intersection_edge,
    std::size_t nb_patches
  ) : patches(nb_patches)
    , pm(pm)
    , patch_ids(patch_ids)
    , fids(fids)
    , is_intersection_edge(is_intersection_edge)
  {
    typedef boost::graph_traits<PolygonMesh> GT;
    typedef typename GT::face_descriptor face_descriptor;

    for(face_descriptor f : faces(pm))
      patches[patch_ids[ get(fids, f) ]].faces.push_back( f );
  }

  Patch_description<PolygonMesh>& operator[](std::size_t i) {
    if ( !patches[i].is_initialized )
    {
      extract_patch_simplices(pm, patches[i], is_intersection_edge);
      patches[i].is_initialized=true;
    }
    return patches[i];
  }

  /// debug
  std::ostream& dump_patch(std::size_t i, std::ostream& out)
  {
    typedef boost::graph_traits<PolygonMesh> GT;
    typedef typename GT::vertex_descriptor vertex_descriptor;
    typedef typename GT::halfedge_descriptor halfedge_descriptor;
    typedef typename GT::face_descriptor face_descriptor;

    Patch_description<PolygonMesh>& patch=this->operator[](i);

    std::stringstream ss;
    ss.precision(17);
    std::map<vertex_descriptor, int> vertexid;
    int id=0;
    for(vertex_descriptor vh : patch.interior_vertices)
    {
      vertexid[vh]=id++;
      ss << get(vertex_point, pm, vh) << "\n";
    }

    for(halfedge_descriptor hh : patch.shared_edges)
    {
      if( vertexid.insert( std::make_pair(target(hh,pm), id) ).second )
      {
        ss << get(vertex_point, pm, target(hh, pm)) << "\n";
        ++id;
      }
      if( vertexid.insert( std::make_pair(source(hh,pm), id) ).second )
      {
        ss << get(vertex_point, pm, source(hh, pm)) << "\n";
        ++id;
      }
    }

    out << "OFF\n" << id << " " << patch.faces.size() << " 0\n";
    out << ss.str();

    for(face_descriptor f : patch.faces)
    {
      out << "3 " << vertexid.at(source(halfedge(f,pm),pm)) <<
             " "  << vertexid.at(target(halfedge(f,pm),pm)) <<
             " "  << vertexid.at(target(next(halfedge(f,pm),pm),pm)) << "\n";
    }

    return out;
  }

  void dump_patches(const boost::dynamic_bitset<>& selection, std::string prefix)
  {
    for (std::size_t i=selection.find_first();
                     i < selection.npos; i = selection.find_next(i))
    {
      std::stringstream ss;
      ss << prefix << "-" << i << ".off";
      std::ofstream output(ss.str().c_str());
      dump_patch(i, output);
    }
  }
  void dump_patches(std::string prefix)
  {
    for (std::size_t i=0;i <patches.size() ; ++i)
    {
      std::stringstream ss;
      ss << prefix << "-" << i << ".off";
      std::ofstream output(ss.str().c_str());
      dump_patch(i, output);
    }
  }
};

template <class PolygonMesh,
          class MarkedEdgeSet>
typename boost::graph_traits<PolygonMesh>::halfedge_descriptor
next_marked_halfedge_around_target_vertex(
   typename boost::graph_traits<PolygonMesh>::halfedge_descriptor h,
   const PolygonMesh& pm,
   const MarkedEdgeSet& marked_edges)
{
  CGAL_assertion( marked_edges.count(edge(h,pm))!= 0 );
  typename boost::graph_traits<PolygonMesh>::halfedge_descriptor nxt =
    next(h, pm);
  while( !marked_edges.count(edge(nxt,pm)) )
  {
    nxt=next(opposite(nxt,pm),pm);
  }
  CGAL_assertion(nxt!=h);
  return nxt;
}

template <class PolygonMesh,
          class EdgeMap,
          class VertexMap,
          class VertexPointMap1,
          class VertexPointMap2,
          class VertexPointMapOut,
          class IntersectionEdgeMap,
          class UserVisitor>
void import_polyline(
  PolygonMesh& output,
  typename boost::graph_traits<PolygonMesh>::halfedge_descriptor h1,
  typename boost::graph_traits<PolygonMesh>::halfedge_descriptor h2,
  const PolygonMesh& pm1,
  const PolygonMesh& pm2,
  std::size_t nb_segments,
  EdgeMap& pm1_to_output_edges,
  EdgeMap& pm2_to_output_edges,
  VertexMap& pm1_to_output_vertices,
  VertexMap& pm2_to_output_vertices,
  const IntersectionEdgeMap& intersection_edges1,
  const IntersectionEdgeMap& intersection_edges2,
  const VertexPointMap1& vpm1,
  const VertexPointMap2& /*vpm2*/,
  const VertexPointMapOut& vpm_out,
  std::vector<typename boost::graph_traits<PolygonMesh>
                ::edge_descriptor>& output_shared_edges,
  UserVisitor& user_visitor)
{
  typedef boost::graph_traits<PolygonMesh> GT;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;
  typedef typename GT::vertex_descriptor vertex_descriptor;

  output_shared_edges.push_back(add_edge(output));
  halfedge_descriptor h_out = halfedge(output_shared_edges.back(),output);

  auto set_output_vertex = [&](vertex_descriptor v, halfedge_descriptor h_out){
    user_visitor.before_vertex_copy(v, pm1, output);
    vertex_descriptor new_v = add_vertex(output);
    set_halfedge(new_v, h_out, output);
    put(vpm_out, new_v, get(vpm1, v));
    user_visitor.after_vertex_copy(v, pm1, new_v, output);
    return new_v;
  };

  //make sure the first vertex does not already exist
  vertex_descriptor src = GT::null_vertex();
  if( get(pm1_to_output_vertices, source(h1,pm1)) == GT::null_vertex() )
  {
    src = set_output_vertex(source(h1, pm1), opposite(h_out, output));
    put(pm1_to_output_vertices, source(h1,pm1), src);
    put(pm2_to_output_vertices, source(h2,pm2), src);
  }
  else
    src = get(pm1_to_output_vertices, source(h1,pm1));

  //make sure the target vertex does not already exist if it is a polyline endpoint
  vertex_descriptor tgt = GT::null_vertex();
  if ( nb_segments==1 )
  {
    if( get(pm1_to_output_vertices, target(h1,pm1)) == GT::null_vertex() )
    {
      tgt = set_output_vertex(target(h1, pm1), h_out);
      put(pm1_to_output_vertices, target(h1,pm1), tgt);
      put(pm2_to_output_vertices, target(h2,pm2), tgt);
    }
    else
      tgt = get(pm1_to_output_vertices, target(h1,pm1));
  }
  else
  {
    tgt = set_output_vertex(target(h1, pm1), h_out);
    put(pm1_to_output_vertices, target(h1,pm1), tgt);
    put(pm2_to_output_vertices, target(h2,pm2), tgt);
  }

  //update source and target vertex of the edge created
  set_target(h_out, tgt, output);
  set_target(opposite(h_out,output), src, output);

  halfedge_descriptor prev_out=h_out;
  halfedge_descriptor prev1=h1;
  halfedge_descriptor prev2=h2;

  //set the correspondence
  put(pm1_to_output_edges, edge(prev1, pm1), edge(prev_out, output));
  put(pm2_to_output_edges, edge(prev2, pm2), edge(prev_out, output));

  user_visitor.intersection_edge_copy(prev1, pm1, prev2, pm2, h_out, output);

  src=tgt;
  for (std::size_t i=1; i<nb_segments; ++i)
  {
    //create new edge
    output_shared_edges.push_back(add_edge(output));
    h_out = halfedge(output_shared_edges.back(),output);
    //get the new edge
    h1 = next_marked_halfedge_around_target_vertex(prev1, pm1, intersection_edges1);
    h2 = next_marked_halfedge_around_target_vertex(prev2, pm2, intersection_edges2);

    user_visitor.intersection_edge_copy(h1, pm1, h2, pm2, h_out, output);

    //if this is the final segment, only create a target vertex if it does not exist
    if (i+1!=nb_segments)
    {
      tgt = set_output_vertex(target(h1, pm1), h_out);
      put( pm1_to_output_vertices, target(h1,pm1), tgt );
      put( pm2_to_output_vertices, target(h2,pm2), tgt );
    }
    else{
      if( get(pm1_to_output_vertices, target(h1,pm1)) == GT::null_vertex() )
      {
        tgt = set_output_vertex(target(h1, pm1), h_out);
        put( pm1_to_output_vertices, target(h1,pm1), tgt );
        put( pm2_to_output_vertices, target(h2,pm2), tgt );
      }
      else
        tgt = get(pm1_to_output_vertices, target(h1,pm1));
    }

    set_target(h_out, tgt, output);
    set_target(opposite(h_out, output), src, output);

    prev_out=h_out;
    prev1 = h1;
    prev2 = h2;
    src = tgt;

    put(pm1_to_output_edges, edge(prev1, pm1), edge(prev_out, output));
    put(pm2_to_output_edges, edge(prev2, pm2), edge(prev_out, output));
  }
}

template <class TriangleMesh, class EdgeMap, class VertexMap, bool reverse_patch_orientation>
struct Triangle_mesh_extension_helper
{
  typedef boost::graph_traits<TriangleMesh> GT;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;
  typedef typename GT::edge_descriptor edge_descriptor;
  typedef typename GT::face_descriptor face_descriptor;

  EdgeMap& tm_to_output_edges;
  VertexMap& tm_to_output_vertices;
  const TriangleMesh& tm;
  TriangleMesh& output;

  Triangle_mesh_extension_helper(EdgeMap& tm_to_output_edges,
                                 VertexMap& tm_to_output_vertices,
                                 const TriangleMesh& tm,
                                 TriangleMesh& output)
    : tm_to_output_edges(tm_to_output_edges)
    , tm_to_output_vertices(tm_to_output_vertices)
    , tm(tm)
    , output(output)
  {}

  halfedge_descriptor get_hedge(halfedge_descriptor h_tm)
  {
    edge_descriptor key = edge(h_tm, tm);
    edge_descriptor value = get(tm_to_output_edges, key);
    CGAL_assertion( value != edge(GT::null_halfedge(), tm) );
    if constexpr(reverse_patch_orientation)
      return target(halfedge(value, output), output) == get(tm_to_output_vertices, target(h_tm, tm))
            ? opposite(halfedge(value, output), output)
            : halfedge(value, output);
    else
      return target(halfedge(value, output), output) == get(tm_to_output_vertices, target(h_tm, tm))
            ? halfedge(value, output)
            : opposite(halfedge(value, output), output);
  }

  std::array<halfedge_descriptor,3>
  halfedges(face_descriptor f)
  {
     halfedge_descriptor h=halfedge(f,tm);
     if constexpr(reverse_patch_orientation)
      return make_array( get_hedge( h ),
                         get_hedge( prev(h,tm) ),
                         get_hedge( next(h,tm) ) );
    else
      return make_array( get_hedge( h ),
                         get_hedge( next(h,tm) ),
                         get_hedge( prev(h,tm) ) );
  }
};

template < bool reverse_patch_orientation,
           class TriangleMesh,
           class PatchContainer,
           class VertexPointMap,
           class VertexPointMapOut,
           class EdgeMarkMapOut,
           class EdgeMarkMapIn ,
           class EdgetoEdgeMap,
           class VertextoVertexMap,
           class UserVisitor>
void append_patches_to_triangle_mesh(
  TriangleMesh& output,
  const boost::dynamic_bitset<>& patches_to_append,
  PatchContainer& patches,
  const VertexPointMapOut& vpm_out,
  const VertexPointMap& vpm_tm,
  EdgeMarkMapOut& edge_mark_map_out,
  const EdgeMarkMapIn& edge_mark_map_in,
  EdgetoEdgeMap& tm_to_output_edges,
  VertextoVertexMap& tm_to_output_vertices,
  UserVisitor& user_visitor)
{
  typedef boost::graph_traits<TriangleMesh> GT;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;
  typedef typename GT::edge_descriptor edge_descriptor;
  typedef typename GT::vertex_descriptor vertex_descriptor;
  typedef typename GT::face_descriptor face_descriptor;

  CGAL::Real_timer t;
  t.start();

  const TriangleMesh& tm = patches.pm;
  Triangle_mesh_extension_helper<TriangleMesh, EdgetoEdgeMap, VertextoVertexMap, reverse_patch_orientation> helper(tm_to_output_edges, tm_to_output_vertices, tm, output);

  std::vector<std::size_t> ids_of_patches_to_append;
  ids_of_patches_to_append.reserve(patches_to_append.count());
  for (std::size_t i=patches_to_append.find_first();
                   i < patches_to_append.npos;
                   i = patches_to_append.find_next(i))
  {
    ids_of_patches_to_append.push_back(i);
  }

  for (std::size_t i : ids_of_patches_to_append)
  {
    #ifdef CGAL_COREFINEMENT_POLYHEDRA_DEBUG
    #warning the size of tm_to_output_edges will increase at each step \
             when adding new patches  by the size of internal edges. \
             Maybe the use of a copy would be better since we do not need \
             the internal edges added?
    #endif

    Patch_description<TriangleMesh>& patch = patches[i];

    // copy the vertices in the output without connecting them
    for(vertex_descriptor v: patch.interior_vertices)
    {
      vertex_descriptor new_v = add_vertex(output);
      set_halfedge(new_v, GT::null_halfedge(), output);

      put(tm_to_output_vertices, v, new_v);
      put(vpm_out, new_v, get(vpm_tm, v));
    }

    // insert interior halfedges and connect vertices
    for(halfedge_descriptor h : patch.interior_edges)
    {
      user_visitor.before_edge_copy(h, tm, output);
      edge_descriptor new_edge = add_edge(output), ed = edge(h,tm);
      user_visitor.after_edge_copy(h, tm, halfedge(new_edge, output), output);

      // copy the mark on input edge to the output edge
      copy_edge_mark<TriangleMesh>(ed, new_edge,
                                   edge_mark_map_in, edge_mark_map_out);

      halfedge_descriptor new_h = halfedge(new_edge, output);
      put(tm_to_output_edges, ed, new_edge);

      set_face(new_h, GT::null_face(), output);
      set_face(opposite(new_h, output), GT::null_face(), output);

      CGAL_assertion(is_border(new_h, output));
      CGAL_assertion(is_border(opposite(new_h,output), output));

      //create a copy of interior vertices only once
      auto fill_vertex = [&](vertex_descriptor v, vertex_descriptor new_v, halfedge_descriptor new_h){
        user_visitor.before_vertex_copy(v, tm, output);
        set_halfedge(new_v, new_h, output);
        user_visitor.after_vertex_copy(v, tm, new_v, output);
      };

      vertex_descriptor tgt = target(h, tm);
      vertex_descriptor new_tgt = get(tm_to_output_vertices, tgt);
      CGAL_assertion( new_tgt != GT::null_vertex() );
      set_target(new_h, new_tgt, output);
      if (  halfedge(tgt,tm)==h &&
            halfedge(new_tgt, output) == GT::null_halfedge())
      {
        fill_vertex(tgt, new_tgt, new_h);
      }

      vertex_descriptor src = source(h, tm);
      vertex_descriptor new_src = get(tm_to_output_vertices, src);
      CGAL_assertion( new_src != GT::null_vertex() );
      halfedge_descriptor new_h_opp = opposite(new_h, output);
      set_target(new_h_opp, new_src, output);
      if (  halfedge(src,tm)==opposite(h, tm) &&
            halfedge(new_src, output) == GT::null_halfedge())
      {
        fill_vertex(src, new_src, new_h_opp);
      }
    }
  }

  for (std::size_t i : ids_of_patches_to_append)
  {
    Patch_description<TriangleMesh>& patch = patches[i];

    //create faces and connect halfedges
    for(face_descriptor f : patch.faces)
    {
      std::array<halfedge_descriptor,3> hedges = helper.halfedges(f);

      user_visitor.before_face_copy(f, patches.pm, output);
      face_descriptor new_f = add_face(output);
      user_visitor.after_face_copy(f, patches.pm, new_f, output);
      set_halfedge(new_f, hedges[0], output);

      for (int i=0;i<3;++i)
      {
        CGAL_assertion(hedges[i]!=GT::null_halfedge());
        set_next(hedges[i], hedges[(i+1)%3], output);
        set_face(hedges[i], new_f, output);
      }
    }
  }

  // handle interior edges that are on the border of the mesh:
  // they do not have a prev/next pointer set since only the pointers
  // of patch interior halfedges part a face have been. In the following
  // (i) we set the next/prev pointer around interior vertices on the mesh
  // boundary and (ii) we collect interior mesh border halfedges incident to
  // a patch border vertex and set their next/prev pointer (possibly of
  // another patch)

  // Containers used for step (ii) for collecting mesh border halfedges
  // with source/target on an intersection polyline that needs it prev/next
  // pointer to be set
  std::vector<halfedge_descriptor> border_halfedges_source_to_link;
  std::vector<halfedge_descriptor> border_halfedges_target_to_link;
  for (std::size_t i : ids_of_patches_to_append)
  {
    Patch_description<TriangleMesh>& patch = patches[i];
    for(halfedge_descriptor h : patch.border_edges){
      if( (!reverse_patch_orientation && patch.border_with_shared_target.count(h)) ||
          ( reverse_patch_orientation && patch.border_with_shared_source.count(h)) )
        continue; // since the next halfedge should not be in the same patch
      halfedge_descriptor h_out=helper.get_hedge(h);
      halfedge_descriptor h_out_next = reverse_patch_orientation
                                        ? helper.get_hedge(prev(h,tm))
                                        : helper.get_hedge(next(h,tm));
      CGAL_assertion(is_border(h_out,output) && is_border(h_out_next,output));
      set_next(h_out, h_out_next, output);
    }
    if(reverse_patch_orientation){
      border_halfedges_target_to_link.insert(border_halfedges_target_to_link.begin(),
                                             patch.border_with_shared_source.begin(),
                                             patch.border_with_shared_source.end());
      border_halfedges_source_to_link.insert(border_halfedges_target_to_link.begin(),
                                             patch.border_with_shared_target.begin(),
                                             patch.border_with_shared_target.end());
    } else {
      border_halfedges_target_to_link.insert(border_halfedges_target_to_link.begin(),
                                             patch.border_with_shared_target.begin(),
                                             patch.border_with_shared_target.end());
      border_halfedges_source_to_link.insert(border_halfedges_target_to_link.begin(),
                                             patch.border_with_shared_source.begin(),
                                             patch.border_with_shared_source.end());
    }
  }

  // now the step (ii) we look for the candidate halfedge by turning around
  // the vertex in the direction of the interior of the patch
  for(halfedge_descriptor h_out : border_halfedges_target_to_link)
  {
    halfedge_descriptor candidate =
      opposite(prev(opposite(h_out, output), output), output);
    CGAL_assertion_code(halfedge_descriptor start=candidate);
    while (!is_border(candidate, output)){
      candidate=opposite(prev(candidate, output), output);
      CGAL_assertion(candidate!=start);
    }
    set_next(h_out, candidate, output);
  }

  for(halfedge_descriptor h_out : border_halfedges_source_to_link)
  {
    halfedge_descriptor candidate =
    opposite(next(opposite(h_out, output), output), output);
    while (!is_border(candidate, output))
      candidate = opposite(next(candidate, output), output);
    set_next(candidate, h_out, output);
  }

  for (std::size_t i : ids_of_patches_to_append)
  {
    Patch_description<TriangleMesh>& patch = patches[i];
    // For all patch boundary vertices, update the vertex pointer
    // of all but the vertex halfedge
    for(halfedge_descriptor h : patch.shared_edges)
    {
      //check for a halfedge pointing inside an already imported patch
      halfedge_descriptor h_out = helper.get_hedge(h);
      CGAL_assertion( next(h_out, output)!=GT::null_halfedge() );
      // update the pointers on the target
      halfedge_descriptor next_around_target=h_out;
      vertex_descriptor v=target(h_out, output);
      do{
        next_around_target = opposite(next(next_around_target, output), output);
        set_target(next_around_target, v, output);
      }while(next(next_around_target, output)!=GT::null_halfedge() &&
             next_around_target!=h_out &&
             !is_border(next_around_target, output));
      // update the pointers on the source
      halfedge_descriptor next_around_source=prev(h_out, output);
      CGAL_assertion(next_around_source!=GT::null_halfedge());
      v = source(h_out, output);
      do{
        set_target(next_around_source, v, output);
        next_around_source = prev(opposite(next_around_source, output), output);
      }while( next_around_source!=GT::null_halfedge() &&
              next_around_source!=opposite(h_out, output) &&
              !is_border(next_around_source, output));
    }
  }
}

template < class TriangleMesh,
           class IntersectionEdgeMap,
           class VertexPointMap1,
           class VertexPointMap2,
           class VertexPointMapOut,
           class EdgeMarkMap1,
           class EdgeMarkMap2,
           class EdgeMarkMapOut,
           class IntersectionPolylines,
           class PatchContainer1,
           class PatchContainer2,
           class UserVisitor>
void fill_new_triangle_mesh(
  TriangleMesh& output,
  const boost::dynamic_bitset<>& patches_of_tm1_to_import,
  const boost::dynamic_bitset<>& patches_of_tm2_to_import,
  PatchContainer1& patches_of_tm1,
  PatchContainer2& patches_of_tm2,
  bool reverse_orientation_of_patches_from_tm1,
  bool reverse_orientation_of_patches_from_tm2,
  const IntersectionPolylines& polylines,
  const IntersectionEdgeMap& intersection_edges1,
  const IntersectionEdgeMap& intersection_edges2,
  const VertexPointMap1& vpm1,
  const VertexPointMap2& vpm2,
  const VertexPointMapOut& vpm_out,
  const EdgeMarkMap1& edge_mark_map1,
  const EdgeMarkMap2& edge_mark_map2,
        EdgeMarkMapOut& edge_mark_map_out,
  std::vector< typename boost::graph_traits<TriangleMesh>::edge_descriptor>&
                                                            output_shared_edges,
  UserVisitor& user_visitor)
{
  using GT = boost::graph_traits<TriangleMesh>;
  using vertex_descriptor = typename GT::vertex_descriptor;
  using edge_descriptor = typename GT::edge_descriptor;

  using V2V_tag = typename CGAL::dynamic_vertex_property_t<vertex_descriptor>;
  using Vertex_to_vertex_map = typename boost::property_map<TriangleMesh, V2V_tag>::const_type;

  using E2E_tag = typename CGAL::dynamic_edge_property_t<SM_Edge_index>;
  using Edge_to_edge_map = typename boost::property_map<TriangleMesh, E2E_tag>::const_type;

  const TriangleMesh& tm1 = patches_of_tm1.pm;
  const TriangleMesh& tm2 = patches_of_tm2.pm;

  Vertex_to_vertex_map tm1_to_output_vertices = get(V2V_tag(), tm1, GT::null_vertex()),
                       tm2_to_output_vertices = get(V2V_tag(), tm2, GT::null_vertex());
  Edge_to_edge_map tm1_to_output_edges = get(E2E_tag(), tm1, edge(GT::null_halfedge(), tm1)),
                   tm2_to_output_edges = get(E2E_tag(), tm2, edge(GT::null_halfedge(), tm2));

  // this is the minimal number of edges that will be marked (intersection edge).
  // We cannot easily have the total number since some patch interior edges might be marked
  output_shared_edges.reserve(std::accumulate(polylines.lengths.begin(),polylines.lengths.end(),std::size_t(0)) );

  //add a polyline inside O for each intersection polyline
  std::size_t nb_polylines = polylines.lengths.size();

  for (std::size_t i=0; i < nb_polylines; ++i)
    if (!polylines.to_skip.test(i))
      import_polyline(output,
                      polylines.tm1[i], polylines.tm2[i],
                      tm1, tm2,
                      polylines.lengths[i],
                      tm1_to_output_edges, tm2_to_output_edges,
                      tm1_to_output_vertices, tm2_to_output_vertices,
                      intersection_edges1, intersection_edges2,
                      vpm1, vpm2, vpm_out,
                      output_shared_edges,
                      user_visitor);

  //import patches of tm1
  if (reverse_orientation_of_patches_from_tm1)
    append_patches_to_triangle_mesh<true>(output,
                                          patches_of_tm1_to_import,
                                          patches_of_tm1,
                                          vpm_out,
                                          vpm1,
                                          edge_mark_map_out,
                                          edge_mark_map1,
                                          tm1_to_output_edges,
                                          tm1_to_output_vertices,
                                          user_visitor);
  else
    append_patches_to_triangle_mesh<false>(output,
                                           patches_of_tm1_to_import,
                                           patches_of_tm1,
                                           vpm_out,
                                           vpm1,
                                           edge_mark_map_out,
                                           edge_mark_map1,
                                           tm1_to_output_edges,
                                           tm1_to_output_vertices,
                                           user_visitor);

  //import patches from tm2
  if (reverse_orientation_of_patches_from_tm2)
    append_patches_to_triangle_mesh<true>(output,
                                          patches_of_tm2_to_import,
                                          patches_of_tm2,
                                          vpm_out,
                                          vpm2,
                                          edge_mark_map_out,
                                          edge_mark_map2,
                                          tm2_to_output_edges,
                                          tm2_to_output_vertices,
                                          user_visitor);
  else
    append_patches_to_triangle_mesh<false>(output,
                                           patches_of_tm2_to_import,
                                           patches_of_tm2,
                                           vpm_out,
                                           vpm2,
                                           edge_mark_map_out,
                                           edge_mark_map2,
                                           tm2_to_output_edges,
                                           tm2_to_output_vertices,
                                           user_visitor);
}

template <class TriangleMesh,
          class PatchContainer,
          class EdgeMap,
          class VertexMap,
          class UserVisitor>
void disconnect_patches(
  TriangleMesh& tm1,
  const boost::dynamic_bitset<>& patches_to_remove,
  PatchContainer& patches_of_tm1,
  const EdgeMap& tm1_edge_to_tm2_edge, //map intersection edges of tm1 to the equivalent in tm2
        EdgeMap& new_tm1_edge_to_tm2_edge, //map the new intersection edges of tm1 to the equivalent in tm2
        VertexMap& tm1_vertex_to_tm2_vertex, //map intersection vertices of tm1 to the equivalent in tm2
        VertexMap& new_tm1_vertex_to_tm2_vertex, //map the new intersection vertices
        UserVisitor& user_visitor)
{
  typedef boost::graph_traits<TriangleMesh> GT;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;
  typedef typename GT::face_descriptor face_descriptor;

  // disconnect each patch one by one
  for (std::size_t i=patches_to_remove.find_first();
                   i < patches_to_remove.npos;
                   i = patches_to_remove.find_next(i))
  {
    Patch_description<TriangleMesh>& patch=patches_of_tm1[i];
    std::vector<halfedge_descriptor> new_patch_border;

    // start the duplicate step
    std::size_t nb_shared_edges = patch.shared_edges.size();
    new_patch_border.reserve( nb_shared_edges );
    std::unordered_map<halfedge_descriptor, halfedge_descriptor> old_to_new;

    // save faces inside the patch and set the halfedge
    // to be duplicated on the boundary
    std::vector<face_descriptor> face_backup;
    face_backup.reserve( nb_shared_edges );
    for(halfedge_descriptor h : patch.shared_edges)
    {
      face_backup.push_back( face(h, tm1) );
      set_face(h, GT::null_face(), tm1);
    }

    // look for the prev/next halfedges after the patch will be disconnected
    // for the halfedges to be duplicated
    std::vector<halfedge_descriptor> shared_next, shared_prev;
    shared_next.reserve( nb_shared_edges );
    shared_prev.reserve( nb_shared_edges );
    for(halfedge_descriptor h : patch.shared_edges)
    {
      halfedge_descriptor nxt=next(h, tm1);
      while(!is_border(nxt, tm1))
        nxt=next(opposite(nxt, tm1), tm1);
      shared_next.push_back(nxt);

      // setting prev is only needed in case tm1 is open
      // and the intersection polyline intersects its boundary
      halfedge_descriptor prv=prev(h, tm1);
      while( !is_border(prv, tm1) )
        prv=prev(opposite(prv, tm1), tm1);
      shared_prev.push_back(prv);

      set_halfedge(target(h, tm1), h, tm1);
      set_halfedge(source(h, tm1), opposite(h, tm1), tm1); //only needed if tm1 is open
    }

    // now duplicate the edge and set its pointers
    for(std::size_t k=0; k<nb_shared_edges; ++k)
    {
      halfedge_descriptor h = patch.shared_edges[k];
      user_visitor.before_edge_duplicated(h, tm1);
      halfedge_descriptor new_hedge = halfedge(add_edge(tm1), tm1);
      user_visitor.after_edge_duplicated(h, new_hedge, tm1);
      set_next(new_hedge, next(h, tm1), tm1);
      set_next(prev(h, tm1), new_hedge, tm1);
      set_face(new_hedge, face_backup[k], tm1);
      set_target(new_hedge, target(h, tm1), tm1);
      set_target(opposite(new_hedge,tm1), source(h, tm1), tm1);
      set_halfedge(face_backup[k], new_hedge, tm1);

      new_patch_border.push_back(new_hedge);
      set_face(opposite(new_hedge, tm1), GT::null_face(), tm1);
      old_to_new.insert( std::make_pair(h, new_hedge) );

      CGAL_assertion( next(opposite(new_hedge, tm1), tm1)==GT::null_halfedge() );
      CGAL_assertion( prev(opposite(new_hedge, tm1), tm1)==GT::null_halfedge() );
      CGAL_assertion( next(prev(new_hedge, tm1), tm1) == new_hedge );
      CGAL_assertion( prev(next(new_hedge, tm1), tm1) == new_hedge );
    }

    // update next/prev pointer of new hedges in case it is one of the new hedges
    for(halfedge_descriptor h : new_patch_border)
      if (is_border(next(h, tm1), tm1))
        set_next(h, old_to_new[next(h,tm1)], tm1);

    // set next/prev pointers on the border of the neighbor patch
    for(std::size_t k=0; k<nb_shared_edges; ++k)
    {
      halfedge_descriptor h = patch.shared_edges[k];
      set_next(h, shared_next[k], tm1);
      set_next(shared_prev[k], h, tm1);
    }

    // update next/prev pointers on the border of the patch
    for(halfedge_descriptor h : new_patch_border)
    {
      halfedge_descriptor h_opp = opposite(h, tm1);
      //set next pointer if not already set
      if ( next(h_opp, tm1)==GT::null_halfedge() )
      {
        // we visit faces inside the patch we consider
        halfedge_descriptor candidate = opposite(prev(h, tm1), tm1);
        while ( !is_border(candidate, tm1) )
          candidate = opposite(prev(candidate, tm1), tm1);
        set_next(h_opp, candidate, tm1);
      }
      CGAL_assertion( prev(next(h_opp, tm1), tm1)==h_opp );

      // set prev pointer if not already set
      if ( prev(h_opp, tm1) == GT::null_halfedge() )
      {
        halfedge_descriptor candidate = opposite(next(h, tm1), tm1);
        while ( !is_border(candidate, tm1) )
          candidate = opposite(next(candidate, tm1), tm1);
        set_next(candidate, h_opp, tm1);
      }

      CGAL_assertion( prev(next(h_opp, tm1), tm1) == h_opp );
      CGAL_assertion( is_border(prev(h_opp, tm1), tm1) );
      CGAL_assertion( is_border(next(h_opp, tm1), tm1) );
    }
    // end of the duplicate step

    CGAL_assertion( new_patch_border.size() == nb_shared_edges );

    for (std::size_t k=0; k < nb_shared_edges; ++k){
      CGAL_assertion( target(patch.shared_edges[k], tm1) == target(new_patch_border[k], tm1) );
      CGAL_assertion( source(patch.shared_edges[k], tm1) == source(new_patch_border[k], tm1) );
      CGAL_assertion( is_border_edge(new_patch_border[k], tm1) );
      CGAL_assertion( !is_border(new_patch_border[k], tm1) );
      CGAL_assertion( is_border(next(opposite(new_patch_border[k], tm1), tm1), tm1) );
      CGAL_assertion( is_border(prev(opposite(new_patch_border[k], tm1), tm1), tm1) );

      auto e = get(tm1_edge_to_tm2_edge, edge(patch.shared_edges[k], tm1));
      auto src = get(tm1_vertex_to_tm2_vertex, source(patch.shared_edges[k], tm1));
      auto tgt = get(tm1_vertex_to_tm2_vertex, target(patch.shared_edges[k], tm1));
      CGAL_assertion( src != GT::null_vertex() );
      CGAL_assertion( tgt != GT::null_vertex() );

      /* halfedge(edge(h, g)) == h is a requirement of the concept and so the if is useless, TODO Sebastien can you confirm
      put(new_tm1_edge_to_tm2_edge,
        patch.shared_edges[k]==halfedge(edge(patch.shared_edges[k], tm1), tm1)
        ? edge(new_patch_border[k], tm1)
        : edge(opposite(new_patch_border[k], tm1), tm1),
      e);
      */
      put(new_tm1_edge_to_tm2_edge, edge(new_patch_border[k], tm1), e);
      put(new_tm1_vertex_to_tm2_vertex, target(new_patch_border[k], tm1), tgt);
      put(new_tm1_vertex_to_tm2_vertex, source(new_patch_border[k], tm1), src);
    }

    patch.shared_edges.swap(new_patch_border);
  }
}

template <class TriangleMesh,
          class PatchContainer1,
          class PatchContainer2,
          class IntersectionPolylines,
          class EdgeMap,
          class VertexMap,
          class VertexPointMap1,
          class VertexPointMap2,
          class EdgeMarkMapIn1,
          class EdgeMarkMapIn2,
          class EdgeMarkMapOut,
          class UserVisitor>
void compute_inplace_operation_delay_removal_and_insideout(
  TriangleMesh& tm1,
  TriangleMesh& tm2,
  const boost::dynamic_bitset<>& patches_of_tm1_to_keep,
  const boost::dynamic_bitset<>& patches_of_tm2_to_import,
  PatchContainer1& patches_of_tm1,
  PatchContainer2& patches_of_tm2,
  bool reverse_patch_orientation_tm2,
  const IntersectionPolylines& polylines,
  const VertexPointMap1& vpm1,
  const VertexPointMap2& vpm2,
        EdgeMarkMapIn1&,
  const EdgeMarkMapIn2& edge_mark_map2,
  const EdgeMarkMapOut& edge_mark_map_out1,
  EdgeMap& disconnected_patches_edge_to_tm2_edge,
  VertexMap& disconnected_patches_vertex_to_tm2_vertex,
  UserVisitor& user_visitor)
{
  typedef boost::graph_traits<TriangleMesh> GT;
  typedef typename GT::vertex_descriptor vertex_descriptor;
  typedef typename GT::edge_descriptor edge_descriptor;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;

  using V2V_tag = typename CGAL::dynamic_vertex_property_t<vertex_descriptor>;
  using E2E_tag = typename CGAL::dynamic_edge_property_t<edge_descriptor>;

  EdgeMap tm2_edge_to_tm1_edge = get(E2E_tag(), tm2, edge(GT::null_halfedge(), tm2)),
          tm1_edge_to_tm2_edge = get(E2E_tag(), tm1, edge(GT::null_halfedge(), tm1));
  VertexMap tm2_vertex_to_tm1_vertex = get(V2V_tag(), tm2, GT::null_vertex()),
            tm1_vertex_to_tm2_vertex = get(V2V_tag(), tm1, GT::null_vertex());
  //maps intersection edges from tm2 to tm1
  std::size_t nb_polylines = polylines.lengths.size();
  for(std::size_t i=0; i<nb_polylines; ++i)
  {
    halfedge_descriptor h1 = polylines.tm1[i];
    halfedge_descriptor h2 = polylines.tm2[i];
    std::size_t nb_segments = polylines.lengths[i];

    for (std::size_t k=0;;)
    {
      put(tm2_edge_to_tm1_edge, edge(h2, tm2), edge(h1, tm1));
      put(tm1_edge_to_tm2_edge, edge(h1, tm1), edge(h2, tm2));
      put(tm2_vertex_to_tm1_vertex, target(h2, tm2), target(h1, tm1));
      put(tm2_vertex_to_tm1_vertex, source(h2, tm2), source(h1, tm1));
      put(tm1_vertex_to_tm2_vertex, target(h1, tm1), target(h2, tm2));
      put(tm1_vertex_to_tm2_vertex, source(h1, tm1), source(h2, tm2));
      if (++k==nb_segments) break;
      h2 = next_marked_halfedge_around_target_vertex(h2, tm2,
             patches_of_tm2.is_intersection_edge);
      h1 = next_marked_halfedge_around_target_vertex(h1, tm1,
             patches_of_tm1.is_intersection_edge);
    }
  }

#ifdef CGAL_COREFINEMENT_POLYHEDRA_DEBUG
  #warning do not try to disconnect if the patch is isolated? i.e opposite(border_edge_of_patch)->is_border()
#endif
  // disconnect patches inside tm2
  // For the patches scheduled to be removed, their patch descriptions
  // in patches_of_tm1 will be updated so that is_intersection_edge are
  // the newly created halfedges within disconnect_patches.
  // Note that disconnected_patches_edge_to_tm2_edge also refers to those halfedges
  //init the map with the previously filled one (needed when reusing patches in two operations)
  disconnected_patches_edge_to_tm2_edge=tm1_edge_to_tm2_edge;
  disconnected_patches_vertex_to_tm2_vertex=tm1_vertex_to_tm2_vertex;
  disconnect_patches(tm1, ~patches_of_tm1_to_keep, patches_of_tm1,
                     tm1_edge_to_tm2_edge, disconnected_patches_edge_to_tm2_edge,
                     tm1_vertex_to_tm2_vertex, disconnected_patches_vertex_to_tm2_vertex, user_visitor);

  //we import patches from tm2
  if (reverse_patch_orientation_tm2)
    append_patches_to_triangle_mesh<true>(tm1,
                                          patches_of_tm2_to_import,
                                          patches_of_tm2,
                                          vpm1,
                                          vpm2,
                                          edge_mark_map_out1,
                                          edge_mark_map2,
                                          tm2_edge_to_tm1_edge,
                                          tm2_vertex_to_tm1_vertex,
                                          user_visitor);
  else
    append_patches_to_triangle_mesh<false>(tm1,
                                           patches_of_tm2_to_import,
                                           patches_of_tm2,
                                           vpm1,
                                           vpm2,
                                           edge_mark_map_out1,
                                           edge_mark_map2,
                                           tm2_edge_to_tm1_edge,
                                           tm2_vertex_to_tm1_vertex,
                                           user_visitor);
}

template <class TriangleMesh,
          class PatchContainer,
          class EdgeMarkMap>
void
remove_patches(TriangleMesh& tm,
               const boost::dynamic_bitset<>& patches_to_remove,
               PatchContainer& patches,
               const EdgeMarkMap& edge_mark_map)
{
  typedef boost::graph_traits<TriangleMesh> GT;
  typedef typename GT::face_descriptor face_descriptor;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;
  typedef typename GT::vertex_descriptor vertex_descriptor;

  for (std::size_t i=patches_to_remove.find_first();
                   i < patches_to_remove.npos;
                   i = patches_to_remove.find_next(i))
  {
    Patch_description<TriangleMesh>& patch=patches[i];

    // put the halfedges on the boundary of the patch on the boundary of tm
    for(halfedge_descriptor h : patch.shared_edges)
      set_face(h, GT::null_face(), tm);

    // set next/prev relationship of border halfedges
    for(halfedge_descriptor h : patch.shared_edges)
    {
      halfedge_descriptor nxt=next(h, tm);
      while(!is_border(nxt,tm))
        nxt=next(opposite(nxt, tm), tm);
      set_next(h, nxt, tm);
      set_halfedge(target(h, tm), h, tm);
    }

    // edges removed must be unmarked to avoid issues when adding new elements
    // that could be marked because they retrieve a previously set property
    unmark_edges(tm, edge_mark_map, patch.interior_edges);

    // In case a ccb of the patch is not a cycle (the source and target vertices
    // are border vertices), the first halfedge of that ccb will not have its
    // prev pointer set correctly. To fix that, we consider all interior edges
    // and check for one that is on the border of the patch and that is incident
    // to a border vertex and use it to get the missing prev pointer.
    for(halfedge_descriptor h : patch.border_with_shared_target)
    {
      // look for the halfedge belonging to shared_edges
      // having the prev pointer not correctly set
      halfedge_descriptor nxt=next(h, tm);
      while(!is_border(nxt, tm))
        nxt=next(opposite(nxt, tm), tm);
      CGAL_assertion( is_border(nxt, tm) );//we marked it above!
      // now update the prev pointer
      halfedge_descriptor prv=prev(opposite(h, tm), tm);
      set_next(prv, nxt, tm);
      set_halfedge(target(prv, tm), prv, tm);
    }

     //now remove the simplices
    for(halfedge_descriptor h : patch.interior_edges)
      remove_edge(edge(h, tm), tm);
    for(face_descriptor f : patch.faces)
      remove_face(f, tm);
    for(vertex_descriptor v : patch.interior_vertices)
      remove_vertex(v, tm);
  }
}

template <class TriangleMesh,
          class PatchContainer1,
          class PatchContainer2,
          class VertexPointMap1,
          class VertexPointMap2,
          class EdgeMarkMapIn1,
          class EdgeMarkMapIn2,
          class EdgeMarkMapOut1,
          class EdgeToEdgeMap,
          class VertexToVertexMap,
          class UserVisitor>
void compute_inplace_operation(
        TriangleMesh& tm1,
  const TriangleMesh& tm2,
  const boost::dynamic_bitset<>& patches_of_tm1_to_keep,
  const boost::dynamic_bitset<>& patches_of_tm2_to_import,
  PatchContainer1& patches_of_tm1,
  PatchContainer2& patches_of_tm2,
  bool reverse_patch_orientation_tm1,
  bool reverse_patch_orientation_tm2,
  const VertexPointMap1& vpm1,
  const VertexPointMap2& vpm2,
        EdgeMarkMapIn1& edge_mark_map_in1,
  const EdgeMarkMapIn2& edge_mark_map_in2,
        EdgeMarkMapOut1& edge_mark_map_out1,
        EdgeToEdgeMap& tm2_edge_to_tm1_edge,
        VertexToVertexMap& tm2_vertex_to_tm1_vertex,
        UserVisitor& user_visitor)
{
  using edge_descriptor = typename boost::graph_traits<TriangleMesh>::edge_descriptor;
  //clean up patches not kept
  remove_patches(tm1, ~patches_of_tm1_to_keep, patches_of_tm1, edge_mark_map_in1);

  // transfer marks of edges of patches kept to the output edge mark property
  copy_edge_mark<TriangleMesh>(tm1, edge_mark_map_in1, edge_mark_map_out1);

  if (reverse_patch_orientation_tm1){
    Polygon_mesh_processing::reverse_face_orientations_of_mesh_with_polylines(tm1);
    // here we need to update the mapping to use the correct border
    // halfedges while appending the patches from tm2
    for(edge_descriptor e2: edges(tm2)){
      edge_descriptor e1 = get(tm2_edge_to_tm1_edge, e2);
      if (e1 != edge(boost::graph_traits<TriangleMesh>::null_halfedge(), tm1))
        put(tm2_edge_to_tm1_edge, e2, edge(opposite(halfedge(e1, tm1), tm1), tm1));
    }
  }

  //we import patches from tm2
  if ( reverse_patch_orientation_tm2 )
    append_patches_to_triangle_mesh<true>(tm1,
                                          patches_of_tm2_to_import,
                                          patches_of_tm2,
                                          vpm1,
                                          vpm2,
                                          edge_mark_map_out1,
                                          edge_mark_map_in2,
                                          tm2_edge_to_tm1_edge,
                                          tm2_vertex_to_tm1_vertex,
                                          user_visitor);
  else
    append_patches_to_triangle_mesh<false>(tm1,
                                           patches_of_tm2_to_import,
                                           patches_of_tm2,
                                           vpm1,
                                           vpm2,
                                           edge_mark_map_out1,
                                           edge_mark_map_in2,
                                           tm2_edge_to_tm1_edge,
                                           tm2_vertex_to_tm1_vertex,
                                           user_visitor);
}

template <class TriangleMesh,
          class IntersectionPolylines,
          class PatchContainer1,
          class PatchContainer2,
          class EdgeMap,
          class VertexMap>
void compute_border_edge_map(
  const TriangleMesh& tm1,
  const TriangleMesh& tm2,
  const IntersectionPolylines& polylines,
  PatchContainer1& patches_of_tm1,
  PatchContainer2& patches_of_tm2,
  EdgeMap& tm2_edge_to_tm1_edge,
  VertexMap& tm2_vertex_to_tm1_vertex)
{
  typedef boost::graph_traits<TriangleMesh> GT;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;

  std::size_t nb_polylines = polylines.lengths.size();
  for( std::size_t i=0; i<nb_polylines; ++i)
  {
    if (polylines.to_skip.test(i)) continue;
    halfedge_descriptor h1 = polylines.tm1[i];
    halfedge_descriptor h2 = polylines.tm2[i];
    std::size_t nb_segments = polylines.lengths[i];

    for (std::size_t k=0;;)
    {
      put(tm2_edge_to_tm1_edge, edge(h2, tm2), edge(h1, tm1));
      put(tm2_vertex_to_tm1_vertex, target(h2, tm2), target(h1, tm1));
      put(tm2_vertex_to_tm1_vertex, source(h2, tm2), source(h1, tm1));
      if (++k==nb_segments) break;
      h2 = next_marked_halfedge_around_target_vertex(
            h2, tm2, patches_of_tm2.is_intersection_edge);
      h1 = next_marked_halfedge_around_target_vertex(
            h1, tm1, patches_of_tm1.is_intersection_edge);
    }
  }
}


template <class TriangleMesh,
          class PatchContainer1,
          class PatchContainer2,
          class IntersectionPolylines,
          class VertexPointMap1,
          class VertexPointMap2,
          class EdgeMarkMapIn1,
          class EdgeMarkMapIn2,
          class EdgeMarkMapOut1,
          class UserVisitor>
void compute_inplace_operation(
        TriangleMesh& tm1,
  const TriangleMesh& tm2,
  const boost::dynamic_bitset<>& patches_of_tm1_to_keep,
  const boost::dynamic_bitset<>& patches_of_tm2_to_import,
  PatchContainer1& patches_of_tm1,
  PatchContainer2& patches_of_tm2,
  bool reverse_patch_orientation_tm1,
  bool reverse_patch_orientation_tm2,
  const VertexPointMap1& vpm1,
  const VertexPointMap2& vpm2,
  const EdgeMarkMapIn1& edge_mark_map_in1,
  const EdgeMarkMapIn2& edge_mark_map_in2,
  const EdgeMarkMapOut1& edge_mark_map_out1,
  const IntersectionPolylines& polylines,
        UserVisitor& user_visitor)
{
  typedef boost::graph_traits<TriangleMesh> GT;
  typedef typename GT::edge_descriptor edge_descriptor;
  typedef typename GT::vertex_descriptor vertex_descriptor;

  using V2V_tag = typename CGAL::dynamic_vertex_property_t<vertex_descriptor>;
  using Vertex_to_vertex_map = typename boost::property_map<TriangleMesh, V2V_tag>::const_type;

  using E2E_tag = typename CGAL::dynamic_edge_property_t<edge_descriptor>;
  using Edge_to_edge_map = typename boost::property_map<TriangleMesh, E2E_tag>::const_type;

  Edge_to_edge_map tm2_edge_to_tm1_edge = get(E2E_tag(), tm2);
  Vertex_to_vertex_map tm2_vertex_to_tm1_vertex = get(V2V_tag(), tm2, GT::null_vertex());

  //maps intersection edges from tm2 to the equivalent in tm1
  compute_border_edge_map(tm1, tm2,
                          polylines,
                          patches_of_tm1, patches_of_tm2,
                          tm2_edge_to_tm1_edge,
                          tm2_vertex_to_tm1_vertex);

  compute_inplace_operation(tm1, tm2,
                            patches_of_tm1_to_keep, patches_of_tm2_to_import,
                            patches_of_tm1, patches_of_tm2,
                            reverse_patch_orientation_tm1,
                            reverse_patch_orientation_tm2,
                            vpm1,
                            vpm2,
                            edge_mark_map_in1,
                            edge_mark_map_in2,
                            edge_mark_map_out1,
                            tm2_edge_to_tm1_edge,
                            tm2_vertex_to_tm1_vertex,
                            user_visitor);
}

// function used to remove polylines imported or kept that are incident only
// to patches not kept for the operation P_ptr is used for storing
// the result. We look for edges with halfedges both on the border of
// the mesh. The vertices incident only to such edges should be removed.
// Here to detect vertices that should be kept, we abuse the fact that
// the halfedge to be removed and incident to a vertex that should not be
// removed will still have its next pointer set to a halfedge part of
// the result.
template <class TriangleMesh, class PatchContainer>
void remove_unused_polylines(
  TriangleMesh& tm,
  const boost::dynamic_bitset<>& patches_to_remove,
  PatchContainer& patches)
{
  typedef boost::graph_traits<TriangleMesh> GT;
  typedef typename GT::vertex_descriptor vertex_descriptor;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;
  typedef typename GT::edge_descriptor edge_descriptor;

  std::set<vertex_descriptor> vertices_to_remove;
  std::set<edge_descriptor> edges_to_remove;

  for (std::size_t i = patches_to_remove.find_first();
                   i < patches_to_remove.npos;
                   i = patches_to_remove.find_next(i))
  {
    Patch_description<TriangleMesh>& patch=patches[i];
    for(halfedge_descriptor h : patch.shared_edges)
    {
      if (is_border(h, tm) && is_border(opposite(h, tm), tm)){
        vertices_to_remove.insert(target(h, tm));
        vertices_to_remove.insert(source(h, tm));
        edges_to_remove.insert(edge(h, tm));
      }
    }
  }

  std::vector<vertex_descriptor> vertices_kept;
  for(vertex_descriptor v : vertices_to_remove)
  {
    bool to_remove=true;
    for(halfedge_descriptor h : halfedges_around_target(v,tm))
      if (!is_border(h, tm) || !is_border(opposite(h,tm),tm))
      {
        to_remove=false;
        // in case the vertex halfedge was one that is going to remove,
        // update it
        set_halfedge(v, h, tm);
        break;
      }
    if (to_remove)
      remove_vertex(v,tm);
    else
      vertices_kept.push_back(v);
  }

  // update next/prev pointers around vertices in vertices_kept
  for(vertex_descriptor v: vertices_kept)
  {
    halfedge_descriptor h = halfedge(v, tm), start=GT::null_halfedge();

    do{

      halfedge_descriptor starter = h;
      while ( !is_border(h, tm) || is_border(opposite(h, tm), tm) )
      {
        h = opposite(next(h, tm), tm);
        if (starter==h) break; // it might happen that no update is required because
                               // it is already closed thanks to another patch from tm'
                               // (usually incident to coplanar patches)

      }
      if( !is_border(h, tm) || is_border(opposite(h, tm), tm) )
      {
        // nothing to do: the vertex has already been updated and is now in the middle of a patch kept.
        // This function can be called after the stitching of the patches kept, the vertex halfedge
        // can have been updated and no border halfedge might be found
        break;
      }
      halfedge_descriptor in = h;

      if (start==GT::null_halfedge())
        start=in;
      else
        if (start==in)
          break;
      while ( is_border(h, tm) )
        h = opposite(next(h, tm), tm);
      set_next(in, opposite(h, tm), tm);
    }
    while(true);//this loop handles non-manifold vertices
  }

  for(edge_descriptor e : edges_to_remove)
    remove_edge(e,tm);
}

template <class TriangleMesh, class PatchContainer, class EdgeMarkMap>
void remove_disconnected_patches(
  TriangleMesh& tm,
  PatchContainer& patches,
  const boost::dynamic_bitset<>& patches_to_remove,
  EdgeMarkMap& edge_mark_map)
{
  typedef boost::graph_traits<TriangleMesh> GT;
  typedef typename GT::vertex_descriptor vertex_descriptor;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;
  typedef typename GT::face_descriptor face_descriptor;

  for (std::size_t i=patches_to_remove.find_first();
                   i < patches_to_remove.npos;
                   i = patches_to_remove.find_next(i))
  {
    Patch_description<TriangleMesh>& patch = patches[i];

    // edges removed must be unmarked to avoid issues when adding new elements
    // that could be marked because they retrieve a previously set property
    unmark_edges(tm, edge_mark_map, patch.interior_edges);

    for(halfedge_descriptor h : patch.interior_edges)
      remove_edge(edge(h, tm), tm);
    // There is no shared halfedge between duplicated patches even
    // if they were before the duplication. Thus the erase that follows is safe.
    // However remember that vertices were not duplicated which is why their
    // removal is not handled here (still in use or to be removed in
    // remove_unused_polylines())
    for(halfedge_descriptor h : patch.shared_edges)
      remove_edge(edge(h, tm), tm);
    for(face_descriptor f : patch.faces)
      remove_face(f, tm);
    for(vertex_descriptor v : patch.interior_vertices)
      remove_vertex(v, tm);
  }
}

} } } // CGAL::Polygon_mesh_processing::Corefinement

#endif // CGAL_POLYGON_MESH_PROCESSING_INTERNAL_FACE_GRAPH_UTILS_H
