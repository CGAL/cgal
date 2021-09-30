// Copyright (c) 2018 GeometryFactory (France).
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

#ifndef CGAL_POLYGON_MESH_PROCESSING_PLANAR_SEGMENTATION_H
#define CGAL_POLYGON_MESH_PROCESSING_PLANAR_SEGMENTATION_H

#include <CGAL/license/Polygon_mesh_processing/meshing_hole_filling.h>

#include <CGAL/plane_mesh_segmentation.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Projection_traits_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>

#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_polygon_mesh.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing_on_segment_set.h>

#include <CGAL/boost/graph/properties.h>
#include <boost/unordered_map.hpp>

#ifdef DEBUG_PCA
#include <fstream>
#endif
#include <algorithm>

/// @TODO remove Kernel_traits
/// @TODO function to move in segmentation package: segment_via_plane_fitting(in, vci, ecm, fccid, np) (pca is a np option)
///       + version with a range of meshes
/// @TODO function to move in PMP: retriangulate_planar_patches(in, out, vci, ecm, fccid, np) (pca is a np option)
/// @TODO expose mark_corner_vertices (or make it part of segment_via_plane_fitting?)
/// @TODO check we can always retriangulate the mesh without creating non-manifoldness issue.

namespace CGAL{

/// @TODO these predicates should probably end-up in the Kernel
namespace Predicates {

// predicates
template <typename K>
struct Is_angle_close_to_Pi_impl{
  typedef bool result_type;

  /// Compares the value of the angle between the vectors `pq` and `pr`
  /// with the angles \$f\theta\f$ and \f$\pi + \theta\f$.
  /// `false` is returned if the angle is less or equal to \f$\pi\f$.
  /// If the cosinus of the angle is smaller or equal to the cosinus of \f$\theta\f$
  /// `true` is returned and `false` otherwise. \f$theta\f$ is provided using the square of its
  /// cosinus.

  bool
  operator()(const typename K::Point_3& p,
             const typename K::Point_3& q,
             const typename K::Point_3& r,
             const typename K::FT min_cosinus_squared) const
  {
    const typename K::Vector_3 pq(p,q);
    const typename K::Vector_3 pr(p,r);

    typename K::Compute_scalar_product_3 scalar_product = K().compute_scalar_product_3_object();
    typename K::FT pqpr = scalar_product(pq, pr);
    // this assumes that cos_theta_bound was meant to be negative
    if (sign(pqpr) != NEGATIVE) return false;
    return compare( square(pqpr),
                    min_cosinus_squared
                      * scalar_product(pq, pq)
                      * scalar_product(pr, pr) ) != SMALLER;
  }
};

template<typename K, bool has_filtered_predicates = K::Has_filtered_predicates>
struct Is_angle_close_to_Pi
  : public Is_angle_close_to_Pi_impl<K>
{
  using Is_angle_close_to_Pi_impl<K>::operator();
};

template<typename K>
struct Is_angle_close_to_Pi<K, true>
  : public Get_filtered_predicate_RT<K, Is_angle_close_to_Pi_impl>::type
{
  using Get_filtered_predicate_RT<K, Is_angle_close_to_Pi_impl>::type::operator();
};

} //end of Predicates namespace

namespace Polygon_mesh_processing {

namespace Planar_segmentation{

inline std::size_t init_id()
{
  return std::size_t(-1);
}

inline std::size_t default_id()
{
  return std::size_t(-2);
}

inline bool is_init_id(std::size_t i)
{
  return i == init_id();
}

inline bool is_corner_id(std::size_t i)
{
  return i < default_id();
}

template <typename Vector_3>
bool is_vector_positive(const Vector_3& normal)
{
  if (normal.x()==0)
  {
    if (normal.y()==0)
      return normal.z() > 0;
    else
      return normal.y() > 0;
  }
  else
    return normal.x() > 0;
}

struct FaceInfo2
{
  FaceInfo2():m_in_domain(-1){}
  int m_in_domain;

  void set_in_domain()   { m_in_domain=1; }
  void set_out_domain()  { m_in_domain=0; }
  bool visited() const   { return m_in_domain!=-1; }
  bool in_domain() const { return m_in_domain==1; }
};

template <typename TriangleMesh,
          typename VertexPointMap,
          typename halfedge_descriptor,
          typename EdgeIsConstrainedMap>
bool is_target_vertex_a_corner(halfedge_descriptor h,
                               EdgeIsConstrainedMap edge_is_constrained,
                               const TriangleMesh& tm,
                               double min_cosinus_squared,
                               const VertexPointMap& vpm)
{
  typedef typename boost::property_traits<VertexPointMap>::value_type Point_3;
  typedef typename Kernel_traits<Point_3>::type K;
  typedef typename boost::graph_traits<TriangleMesh> graph_traits;

  halfedge_descriptor h2 = graph_traits::null_halfedge();
  for(halfedge_descriptor h_loop : halfedges_around_target(h, tm))
  {
    if (h_loop==h) continue;
    if (get(edge_is_constrained, edge(h_loop, tm)))
    {
      if (h2 != graph_traits::null_halfedge()) return true;
      h2=h_loop;
    }
  }

  // handle case when the graph of constraints does not contains only cycle
  // (for example when there is a tangency between surfaces and is shared)
  if (h2 == graph_traits::null_halfedge()) return true;

  const Point_3& p = get(vpm, target(h, tm));
  const Point_3& q = get(vpm, source(h, tm));
  const Point_3& r = get(vpm, source(h2, tm));

  if (min_cosinus_squared==1)
    return !collinear(p, q, r);
  else
  {
    Predicates::Is_angle_close_to_Pi<K> pred;
    return !pred(p, q, r, min_cosinus_squared);
  }
}

template <typename TriangleMesh,
          typename VertexPointMap,
          typename FaceCCIdMap,
          typename VertexCornerIdMap>
std::size_t
mark_corner_vertices_with_region_growing(
  const TriangleMesh& triangle_mesh,
  const FaceCCIdMap& face_cc_ids,
  VertexCornerIdMap& vertex_corner_id,
  const double max_distance,
  const double cos_value_squared,
  const VertexPointMap& vpm)
{
  using Kernel = typename Kernel_traits<typename boost::property_traits<VertexPointMap>::value_type>::Kernel;
  using Vertex_type = typename boost::property_traits<VertexPointMap>::key_type;
  using FT = typename Kernel::FT;

  // Get face/edge range and types.
  const auto edge_range = edges(triangle_mesh);
  using Face_range          = decltype(faces(triangle_mesh));
  using Edge_range          = decltype(edge_range);
  using Triangle_mesh       = TriangleMesh;
  using Vertex_to_point_map = VertexPointMap;
  using Face_to_region_map  = FaceCCIdMap;

  using Polyline_graph = CGAL::Shape_detection::Polygon_mesh::
    Polyline_graph<Kernel, Triangle_mesh, Face_to_region_map, Face_range, Edge_range, Vertex_to_point_map>;

  using Segment_range = typename Polyline_graph::Segment_range;
  using Segment_map   = typename Polyline_graph::Segment_map;

  using Line_region  = CGAL::Shape_detection::Segment_set::
    Least_squares_line_fit_region<Kernel, Segment_range, Segment_map>;
  using Line_sorting = CGAL::Shape_detection::Segment_set::
    Least_squares_line_fit_sorting<Kernel, Segment_range, Polyline_graph, Segment_map>;

  using Region_growing = CGAL::Shape_detection::
    Region_growing<Segment_range, Polyline_graph, Line_region, typename Line_sorting::Seed_map>;

  // Create a graph of line segments.
  // TODO: Can we use edge_is_constrained directly inside Line_sorting and RG?
  Polyline_graph pgraph(triangle_mesh, CGAL::parameters::
    face_index_map(face_cc_ids).vertex_point_map(vpm));
  const auto& segment_range = pgraph.segment_range();

  // Create instances of the classes Region_type and Sorting.
  Line_region line_region(
    segment_range, CGAL::parameters::
    segment_map(pgraph.segment_map()).
    maximum_distance(static_cast<FT>(max_distance)).
    cosine_value(static_cast<FT>(CGAL::sqrt(cos_value_squared))));

  Line_sorting line_sorting(
    segment_range, pgraph, CGAL::parameters::segment_map(pgraph.segment_map()));
  line_sorting.sort();

  // Create an instance of the region growing class.
  Region_growing region_growing(
    segment_range, pgraph, line_region, line_sorting.seed_map());

  // Run the algorithm.
  std::vector< std::vector<std::size_t> > regions;
  region_growing.detect(std::back_inserter(regions));
  // std::cout << "- found linear regions: " << regions.size() << std::endl;

  std::set<Vertex_type> unique;
  for (const auto& region : regions) {

    // TODO: Can we do it faster?
    // We should not have more than 2 corners per region!
    for (const std::size_t segment_index : region) {

      const std::size_t ei = pgraph.edge_index(segment_index);
      const auto edge = *(edge_range.begin() + ei);
      const auto he = halfedge(edge, triangle_mesh);
      const auto svertex = source(he, triangle_mesh);
      const auto tvertex = target(he, triangle_mesh);

      const auto& sneighbors = pgraph.source_neighbors(segment_index);
      const auto& tneighbors = pgraph.target_neighbors(segment_index);
      CGAL_assertion(sneighbors.size() > 0);
      CGAL_assertion(tneighbors.size() > 0);

      if (sneighbors.size() == 1) {
        put(vertex_corner_id, svertex, default_id());
      } else {
        unique.insert(svertex);
      }

      if (tneighbors.size() == 1) {
        put(vertex_corner_id, tvertex, default_id());
      } else {
        unique.insert(tvertex);
      }
    }
  }

  std::size_t corner_id = 0;
  for (const auto& vertex : unique) {
    put(vertex_corner_id, vertex, corner_id++);
  }
  // std::cout << "- found corners: " << corner_id << std::endl;
  return corner_id;
}

template <typename TriangleMesh,
          typename VertexPointMap,
          typename EdgeIsConstrainedMap,
          typename VertexCornerIdMap>
std::size_t
mark_corner_vertices(
  TriangleMesh& tm,
  EdgeIsConstrainedMap& edge_is_constrained,
  VertexCornerIdMap& vertex_corner_id,
  double min_cosinus_squared,
  const VertexPointMap& vpm)
{
  typedef boost::graph_traits<TriangleMesh> graph_traits;
  std::size_t corner_id = 0;
  for(typename graph_traits::edge_descriptor e : edges(tm))
  {
    if (!get(edge_is_constrained, e)) continue;
    typename graph_traits::halfedge_descriptor h = halfedge(e, tm);

    if (is_init_id(get(vertex_corner_id, target(h, tm))))
    {
      if (is_target_vertex_a_corner(h, edge_is_constrained, tm, min_cosinus_squared, vpm))
        put(vertex_corner_id, target(h, tm), corner_id++);
      else
        put(vertex_corner_id, target(h, tm), default_id());
    }
    if (is_init_id(get(vertex_corner_id, source(h, tm))))
    {
      if (is_target_vertex_a_corner(opposite(h, tm), edge_is_constrained, tm, min_cosinus_squared, vpm))
        put(vertex_corner_id, source(h, tm), corner_id++);
      else
        put(vertex_corner_id, source(h, tm), default_id());
    }
  }

  return corner_id;
}

template <typename CDT>
void mark_face_triangles(CDT& cdt)
{
  //look for a triangle inside the domain of the face
  typename CDT::Face_handle fh = cdt.infinite_face();
  fh->info().set_out_domain();
  std::vector<typename CDT::Edge> queue;
  for (int i=0; i<3; ++i)
    queue.push_back(typename CDT::Edge(fh, i) );
  while(true)
  {
    typename CDT::Edge e = queue.back();
    queue.pop_back();
    e=cdt.mirror_edge(e);
    if (e.first->info().visited()) continue;
    if (cdt.is_constrained(e))
    {
      queue.clear();
      queue.push_back(e);
      break;
    }
    else
    {
      for(int i=1; i<3; ++i)
      {
        typename CDT::Edge candidate(e.first, (e.second+i)%3);
        if (!candidate.first->neighbor(candidate.second)->info().visited())
          queue.push_back( candidate );
      }
      e.first->info().set_out_domain();
    }
  }
  // now extract triangles inside the face
  while(!queue.empty())
  {
    typename CDT::Edge e = queue.back();
    queue.pop_back();
    if (e.first->info().visited()) continue;
    e.first->info().set_in_domain();

    for(int i=1; i<3; ++i)
    {
      typename CDT::Edge candidate(e.first, (e.second+i)%3);
      if (!cdt.is_constrained(candidate) &&
          !candidate.first->neighbor(candidate.second)->info().visited())
      {
        queue.push_back( cdt.mirror_edge(candidate) );
      }
    }
  }
}

template <typename Kernel>
bool add_triangle_faces(const std::vector< std::pair<std::size_t, std::size_t> >& csts,
                        typename Kernel::Vector_3 normal,
                        const std::vector<typename Kernel::Point_3>& corners,
                        std::vector<cpp11::array<std::size_t, 3> >& triangles)
{
  typedef Projection_traits_3<Kernel>                            P_traits;
  typedef Triangulation_vertex_base_with_info_2<std::size_t, P_traits> Vb;
  typedef Triangulation_face_base_with_info_2<FaceInfo2,P_traits>     Fbb;
  typedef Constrained_triangulation_face_base_2<P_traits,Fbb>          Fb;
  typedef Triangulation_data_structure_2<Vb,Fb>                       TDS;
  typedef Exact_predicates_tag                                       Itag;
  typedef Constrained_Delaunay_triangulation_2<P_traits, TDS, Itag>   CDT;
  typedef typename Kernel::Point_3 Point_3;

  std::vector<std::pair<Point_3, std::size_t> > points;
  points.reserve(csts.size()/2);

  typedef std::pair<std::size_t, std::size_t> Id_pair;
  for(const Id_pair& p : csts)
  {
    CGAL_assertion(p.first<corners.size());
    points.push_back( std::make_pair(corners[p.first], p.first) );
  }

  bool reverse_face_orientation = is_vector_positive(normal);
  if (reverse_face_orientation)
    normal=-normal;

  // create cdt and insert points
  P_traits p_traits(normal);
  CDT cdt(p_traits);
  cdt.insert(points.begin(), points.end());

  // note that nbv might be different from points.size() in case of hole
  // tangent to the principal CCB
  std::size_t nbv=cdt.number_of_vertices();

  // insert constrained edges
  boost::unordered_map<std::size_t, typename CDT::Vertex_handle> vertex_map;
  for(typename CDT::Finite_vertices_iterator vit = cdt.finite_vertices_begin(),
                                             end = cdt.finite_vertices_end(); vit!=end; ++vit)
  {
    vertex_map[vit->info()]=vit;
  }

  std::vector< std::pair<typename CDT::Vertex_handle, typename CDT::Vertex_handle> > local_csts;
  local_csts.reserve(csts.size());
  for(const Id_pair& p : csts)
    cdt.insert_constraint(vertex_map[p.first], vertex_map[p.second]);

  if (cdt.number_of_vertices() != nbv)
    return false;

  mark_face_triangles(cdt);

  for (typename CDT::Finite_faces_iterator fit=cdt.finite_faces_begin(),
                                           end=cdt.finite_faces_end(); fit!=end; ++fit)
  {
    if (!fit->info().in_domain()) continue;
    if (cdt.is_infinite(fit)) return false;

    if (reverse_face_orientation)
      triangles.push_back( make_array(fit->vertex(1)->info(),
                                      fit->vertex(0)->info(),
                                      fit->vertex(2)->info()) );
    else
      triangles.push_back( make_array(fit->vertex(0)->info(),
                                      fit->vertex(1)->info(),
                                      fit->vertex(2)->info()) );
  }

  return true;
}

template <typename TriangleMesh,
          typename VertexCornerIdMap,
          typename EdgeIsConstrainedMap,
          typename FaceCCIdMap,
          typename VertexPointMap>
std::pair<std::size_t, std::size_t>
tag_corners_and_constrained_edges(TriangleMesh& tm,
                                  double min_cosinus_squared,
                                  double max_frechet_distance,
                                  bool use_region_growing,
                                  VertexCornerIdMap& vertex_corner_id,
                                  EdgeIsConstrainedMap& edge_is_constrained,
                                  FaceCCIdMap& face_cc_ids,
                                  const VertexPointMap& vpm)
{
  std::size_t nb_cc =
    segment_via_plane_fitting(tm, face_cc_ids,
                              parameters::edge_is_constrained_map(edge_is_constrained)
                                         .minimum_cosinus_squared(min_cosinus_squared)
                                         .maximum_Frechet_distance(max_frechet_distance)
                                         .use_region_growing(use_region_growing));

  std::size_t nb_corners = 0;
  if (!use_region_growing) {
    nb_corners =
      mark_corner_vertices(tm, edge_is_constrained, vertex_corner_id, min_cosinus_squared, vpm);
  } else {
    nb_corners =
      mark_corner_vertices_with_region_growing(tm, face_cc_ids, vertex_corner_id, max_frechet_distance, min_cosinus_squared, vpm);
  }

  return std::make_pair(nb_corners, nb_cc);
}

template <typename TriangleMesh,
          typename VertexCornerIdMap,
          typename EdgeIsConstrainedMap,
          typename FaceCCIdMap,
          typename VertexPointMap>
bool decimate_impl(TriangleMesh& tm,
                   const std::pair<std::size_t, std::size_t>& nb_corners_and_nb_cc,
                   VertexCornerIdMap& vertex_corner_id,
                   EdgeIsConstrainedMap& edge_is_constrained,
                   FaceCCIdMap& face_cc_ids,
                   const VertexPointMap& vpm)
{
  typedef typename boost::property_traits<VertexPointMap>::value_type Point_3;
  typedef typename Kernel_traits<Point_3>::type K;
  typedef typename boost::graph_traits<TriangleMesh> graph_traits;
  typedef typename graph_traits::halfedge_descriptor halfedge_descriptor;
  typedef typename graph_traits::vertex_descriptor vertex_descriptor;
  typedef std::pair<std::size_t, std::size_t> Id_pair;
  std::vector< std::vector<Id_pair> > face_boundaries(nb_corners_and_nb_cc.second);
  std::vector< typename K::Vector_3 > face_normals(nb_corners_and_nb_cc.second, NULL_VECTOR);
  //collect corners
  std::vector< Point_3 > corners(nb_corners_and_nb_cc.first);
  for(vertex_descriptor v : vertices(tm))
  {
    std::size_t i = get(vertex_corner_id, v);
    if ( is_corner_id(i) )
      corners[i]=get(vpm, v);
  }

  // collect maximal constrained edges per cc
  for(halfedge_descriptor h : halfedges(tm))
  {
    if (!get(edge_is_constrained, edge(h, tm)) || is_border(h, tm)) continue;

    std::size_t i1 = get(vertex_corner_id, source(h, tm));
    if ( is_corner_id(i1) )
    {
      halfedge_descriptor h_init = h;
      do{
        std::size_t i2 = get(vertex_corner_id, target(h_init, tm));
        if ( is_corner_id(i2) )
        {
          std::size_t cc_id = get(face_cc_ids, face(h_init, tm));
          face_boundaries[ cc_id ].push_back( Id_pair(i1,i2) );
          if (face_normals[ cc_id ] == NULL_VECTOR)
          {
            face_normals[ cc_id ] = normal(get(vpm, source(h, tm)),
                                           get(vpm, target(h, tm)),
                                           get(vpm, target(next(h, tm), tm)));
          }
          break;
        }

        do{
          h_init=opposite(next(h_init, tm), tm);
        } while( !get(edge_is_constrained, edge(h_init, tm)) );
        h_init=opposite(h_init, tm);
      }
      while(true);
    }
  }

  /// @TODO this is rather drastic in particular if the mesh has almost none simplified faces
/// TODO use add_faces?

  // compute the new mesh
  std::vector< cpp11::array<std::size_t, 3> > triangles;
  triangles.reserve(nb_corners_and_nb_cc.second);

  for(std::size_t cc_id=0; cc_id<nb_corners_and_nb_cc.second; ++cc_id)
  {
    const std::vector< Id_pair >& csts = face_boundaries[cc_id];
    if (csts.size() < 3)
      return false;
    if (csts.size()==3)
    {
      triangles.push_back( make_array(csts[0].first,
                                      csts[0].second,
                                      csts[0].first==csts[1].first ||
                                      csts[0].second==csts[1].first ?
                                      csts[1].second:csts[1].first) );
    }
    else
      if (!add_triangle_faces<K>(csts, face_normals[cc_id], corners, triangles))
        return false;
  }

  if (!is_polygon_soup_a_polygon_mesh(triangles))
    return false;

  //clear(tm);
  tm.clear_without_removing_property_maps();
  polygon_soup_to_polygon_mesh(corners, triangles, tm, parameters::all_default(), parameters::vertex_point_map(vpm));
  return true;
}

template <typename vertex_descriptor,
          typename Point_3,
          typename OutputIterator>
void extract_meshes_containing_a_point(
  const Point_3& pt,
  const std::map<Point_3, std::map<std::size_t, vertex_descriptor> >& point_to_vertex_maps,
  OutputIterator out)
{
  typedef std::pair<const std::size_t, vertex_descriptor> Pair_type;
  for(const Pair_type& p : point_to_vertex_maps.find(pt)->second)
    *out++=p.first;
}

template <typename TriangleMesh,
          typename Point_3,
          typename vertex_descriptor,
          typename EdgeIsConstrainedMap,
          typename VertexIsSharedMap,
          typename VertexPointMap>
void mark_boundary_of_shared_patches_as_constrained_edges(
  std::vector<TriangleMesh*>& mesh_ptrs,
  std::map<Point_3, std::map<std::size_t, vertex_descriptor> >& point_to_vertex_maps,
  std::vector<EdgeIsConstrainedMap>& edge_is_constrained_maps,
  std::vector<VertexIsSharedMap>& vertex_shared_maps,
  const std::vector<VertexPointMap>& vpms)
{
  typedef boost::graph_traits<TriangleMesh> graph_traits;
  typedef typename graph_traits::edge_descriptor edge_descriptor;
  typedef typename graph_traits::halfedge_descriptor halfedge_descriptor;

  std::size_t mesh_id = 0;
  for(TriangleMesh* tm_ptr : mesh_ptrs)
  {
    TriangleMesh& tm=*tm_ptr;
    EdgeIsConstrainedMap& edge_is_constrained = edge_is_constrained_maps[mesh_id];
    VertexIsSharedMap& is_vertex_shared = vertex_shared_maps[mesh_id];

    for(edge_descriptor e : edges(tm))
    {
      if (is_border(e, tm)) continue; //border edges will be automatically marked as constrained

      halfedge_descriptor h = halfedge(e, tm);
      vertex_descriptor src = source(h, tm), tgt = target(h, tm);
      if (get(is_vertex_shared, src) && get(is_vertex_shared, tgt))
      {
        //extract the set of meshes having both vertices
        std::set<std::size_t> src_set, tgt_set, inter_set;
        extract_meshes_containing_a_point(get(vpms[mesh_id], src),
                                          point_to_vertex_maps,
                                          std::inserter(src_set, src_set.begin()));
        extract_meshes_containing_a_point(get(vpms[mesh_id], tgt),
                                          point_to_vertex_maps,
                                          std::inserter(tgt_set, tgt_set.begin()));

        std::set_intersection(src_set.begin(), src_set.end(),
                              tgt_set.begin(), tgt_set.end(),
                              std::inserter(inter_set, inter_set.begin()));

        std::map<std::size_t, vertex_descriptor>& mesh_to_vertex_src =
          point_to_vertex_maps[get(vpms[mesh_id], src)];
        std::map<std::size_t, vertex_descriptor>& mesh_to_vertex_tgt =
          point_to_vertex_maps[get(vpms[mesh_id], tgt)];

        std::set<Point_3> incident_face_points;
        incident_face_points.insert(get(vpms[mesh_id], target(next(h, tm), tm)));
        h=opposite(h, tm);
        incident_face_points.insert(get(vpms[mesh_id], target(next(h, tm), tm)));

        // we mark as constrained edge, any edge that is shared between more than 2 meshes
        // such that at least one of the two incident faces to the edge are not present in
        // all the meshes containing the edge
        for(std::size_t other_mesh_id : inter_set)
        {
          TriangleMesh* other_tm_ptr = mesh_ptrs[other_mesh_id];
          if (other_tm_ptr==&tm) continue;
          vertex_descriptor other_src=mesh_to_vertex_src[other_mesh_id];
          vertex_descriptor other_tgt=mesh_to_vertex_tgt[other_mesh_id];
          std::pair<halfedge_descriptor, bool> hres = halfedge(other_src, other_tgt, *other_tm_ptr);
          if (hres.second)
          {
            if (is_border_edge(hres.first, *other_tm_ptr))
            {
              put(edge_is_constrained, e, true);
              break;
            }
            if (incident_face_points.count(
                  get(vpms[other_mesh_id], target(next(hres.first, *other_tm_ptr), *other_tm_ptr)))==0)
            {
              put(edge_is_constrained, e, true);
              break;
            }
            hres.first=opposite(hres.first, *other_tm_ptr);
            if (incident_face_points.count(
                  get(vpms[other_mesh_id], target(next(hres.first, *other_tm_ptr), *other_tm_ptr)))==0)
            {
              put(edge_is_constrained, e, true);
              break;
            }
          }
        }
      }
    }
    ++mesh_id;
  }
}

template<typename Point_3,
         typename vertex_descriptor,
         typename VertexIsSharedMap >
void propagate_corner_status(
  std::vector<VertexIsSharedMap>& vertex_corner_id_maps,
  std::map<Point_3, std::map<std::size_t, vertex_descriptor> >& point_to_vertex_maps,
  std::vector< std::pair<std::size_t, std::size_t> >& nb_corners_and_nb_cc_all)
{
  typedef std::pair<const Point_3, std::map<std::size_t, vertex_descriptor> > Pair_type;
  for(Pair_type& p : point_to_vertex_maps)
  {
    // if one vertex is a corner, all should be
    typedef std::pair<const std::size_t, vertex_descriptor> Map_pair_type;
    bool is_corner=false;
    for (Map_pair_type& mp : p.second)
    {
      std::size_t mesh_id = mp.first;
      if ( is_corner_id( get(vertex_corner_id_maps[mesh_id], mp.second) ))
      {
        is_corner=true;
        break;
      }
    }
    if (is_corner)
    {
      for(Map_pair_type& mp : p.second)
      {
        std::size_t mesh_id = mp.first;
        if ( !is_corner_id(get(vertex_corner_id_maps[mesh_id], mp.second)) )
          put(vertex_corner_id_maps[mesh_id], mp.second,
            nb_corners_and_nb_cc_all[mesh_id].first++);
      }
    }
  }
}

#ifndef CGAL_DO_NOT_USE_PCA
template <typename TriangleMesh,
          typename EdgeIsConstrainedMap,
          typename CornerIdMap,
          typename VertexPointMap>
std::size_t
mark_extra_corners_with_pca(TriangleMesh& tm,
                            double max_frechet_distance,
                            std::size_t nb_corners,
                            EdgeIsConstrainedMap& edge_is_constrained,
                            CornerIdMap& vertex_corner_id,
                            const VertexPointMap& vpm)
{
  typedef typename Kernel_traits<typename boost::property_traits<VertexPointMap>::value_type>::Kernel IK; // input kernel
  typedef CGAL::Exact_predicates_inexact_constructions_kernel PCA_K;
  typedef typename boost::graph_traits<TriangleMesh> graph_traits;
  typedef typename graph_traits::halfedge_descriptor halfedge_descriptor;
  typedef typename graph_traits::vertex_descriptor vertex_descriptor;

  CGAL::Cartesian_converter<IK, PCA_K> to_pca_k;

  const double max_squared_frechet_distance = max_frechet_distance * max_frechet_distance;

  for(halfedge_descriptor h : halfedges(tm))
  {
    if (!get(edge_is_constrained, edge(h, tm)) || is_border(h, tm)) continue;

    std::size_t i1 = get(vertex_corner_id, source(h, tm));
    if ( is_corner_id(i1) )
    {
      std::size_t i2 = get(vertex_corner_id, target(h, tm));
      if ( is_corner_id(i2)) continue;

      std::vector<vertex_descriptor> edge_boundary_vertices;

      edge_boundary_vertices.push_back(source(h, tm));

      halfedge_descriptor h_init = h;
      do{
        if ( is_corner_id(i2))
          break;
        do{
          h_init=opposite(next(h_init, tm), tm);
        } while( !get(edge_is_constrained, edge(h_init, tm)) );
        h_init=opposite(h_init, tm);
        edge_boundary_vertices.push_back(target(h_init, tm));
        i2 = get(vertex_corner_id, target(h_init, tm));
      }
      while(true);

      //create the set of segments from the chain of vertices
      std::vector<typename PCA_K::Segment_3> edge_boundary_segments;
      std::size_t nb_segments=edge_boundary_vertices.size()-1;
      CGAL_assertion(nb_segments>0);
      edge_boundary_segments.reserve(nb_segments);
      for (std::size_t i=0; i<nb_segments; ++i)
      {
        edge_boundary_segments.push_back(
          typename PCA_K::Segment_3(
            to_pca_k(get(vpm, edge_boundary_vertices[i])),
            to_pca_k(get(vpm, edge_boundary_vertices[i+1]) )) );
      }
      typename PCA_K::Line_3 line;
      typename PCA_K::Point_3 centroid;

      auto does_fitting_respect_distance_bound = [&to_pca_k, &vpm, max_squared_frechet_distance](
        typename std::vector<vertex_descriptor>::const_iterator begin,
        typename std::vector<vertex_descriptor>::const_iterator end,
        const PCA_K::Line_3& line)
      {
        typename PCA_K::Compare_squared_distance_3 compare_squared_distance;
        while ( begin != end)
        {
          if (compare_squared_distance(to_pca_k(get(vpm, *begin++)), line, max_squared_frechet_distance) == LARGER)
            return false;
        }
        return true;
      };

      // first look if the whole vertex chain is a good fit
      linear_least_squares_fitting_3(edge_boundary_segments.begin(),
                                     edge_boundary_segments.end(),
                                     line,
                                     centroid,
                                     Dimension_tag<0>()); /// @TODO: use dimension 1 when the bug in CGAL is fixed

      if ( !does_fitting_respect_distance_bound(edge_boundary_vertices.begin(),
                                                edge_boundary_vertices.end(),
                                                line)) continue;
#ifdef DEBUG_PCA
      std::cout << "  The whole chain cannot be fit, nb_segments=" << nb_segments << " line: " << line << "\n";
#endif
      // iteratively increase the boundary edge length while it is a good fit, then
      // continue with the next part
      std::size_t b=0, e=2;
      while(e<=nb_segments)
      {
#ifdef DEBUG_PCA
      std::cout << "  b=" << b << " and e="<< e << "\n";
#endif
        linear_least_squares_fitting_3(edge_boundary_segments.begin()+b,
                                       edge_boundary_segments.begin()+e,
                                       line,
                                       centroid,
                                       Dimension_tag<0>());  /// @TODO: use dimension 1 when the bug in CGAL is fixed
        CGAL_assertion(edge_boundary_vertices.size() >= e+1 );

        if (does_fitting_respect_distance_bound(edge_boundary_vertices.begin()+b,
                                                edge_boundary_vertices.begin()+e+1,
                                                line))
        {
          ++e;
        }
        else
        {
          CGAL_assertion( !is_corner_id(get(vertex_corner_id, edge_boundary_vertices[e-1])) );
          put(vertex_corner_id, edge_boundary_vertices[e-1], nb_corners++);
          b=e;
          e+=2;
        }
      }
    }
  }

  return nb_corners;
}
#endif

template <typename TriangleMeshRange,
          typename MeshMap,
          typename VertexPointMap>
bool decimate_meshes_with_common_interfaces_impl(TriangleMeshRange& meshes,
                                                 MeshMap mesh_map,
                                                 double max_frechet_distance, // !=0 if pca should be used
                                                 double min_cosinus_squared,
                                                 const std::vector<VertexPointMap>& vpms)
{
  typedef typename boost::property_traits<MeshMap>::value_type Triangle_mesh;
  typedef typename std::iterator_traits<typename TriangleMeshRange::iterator>::value_type Mesh_descriptor;
  typedef typename boost::property_traits<VertexPointMap>::value_type Point_3;
  typedef typename boost::graph_traits<Triangle_mesh> graph_traits;
  typedef typename graph_traits::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits::edge_descriptor edge_descriptor;
  typedef typename graph_traits::face_descriptor face_descriptor;
  CGAL_assertion(min_cosinus_squared>=0);
  const bool use_PCA = max_frechet_distance!=0;
#ifdef CGAL_DO_NOT_USE_PCA
  if (use_PCA)
    std::cerr << "Warning: ask for using PCA while it was disabled at compile time!\n";
#endif

  // declare and init all property maps
  typedef typename boost::property_map<Triangle_mesh, CGAL::dynamic_face_property_t<std::size_t> >::type Face_cc_map;
  typedef typename boost::property_map<Triangle_mesh, CGAL::dynamic_edge_property_t<bool> >::type Edge_is_constrained_map;
  typedef typename boost::property_map<Triangle_mesh, CGAL::dynamic_vertex_property_t<std::size_t> >::type Vertex_corner_id_map;
  typedef typename boost::property_map<Triangle_mesh, CGAL::dynamic_vertex_property_t<bool> >::type Vertex_is_shared_map;

  std::vector<Vertex_is_shared_map > vertex_shared_maps;
  std::vector<Edge_is_constrained_map> edge_is_constrained_maps;
  std::vector<Vertex_corner_id_map> vertex_corner_id_maps;
  std::vector<Face_cc_map> face_cc_ids_maps;
  const std::size_t nb_meshes = meshes.size();
  vertex_corner_id_maps.reserve(nb_meshes);
  vertex_shared_maps.reserve(nb_meshes);
  edge_is_constrained_maps.reserve(nb_meshes);
  face_cc_ids_maps.reserve(nb_meshes);

  std::vector<Triangle_mesh*> mesh_ptrs;
  mesh_ptrs.reserve(nb_meshes);
  for(Mesh_descriptor& md : meshes)
    mesh_ptrs.push_back( &(mesh_map[md]) );


  for(Triangle_mesh* tm_ptr : mesh_ptrs)
  {
    Triangle_mesh& tm = *tm_ptr;

    vertex_shared_maps.push_back( get(CGAL::dynamic_vertex_property_t<bool>(), tm) );
    for(vertex_descriptor v : vertices(tm))
      put(vertex_shared_maps.back(), v, false);

    edge_is_constrained_maps.push_back( get(CGAL::dynamic_edge_property_t<bool>(), tm) );
    for(edge_descriptor e : edges(tm))
      put(edge_is_constrained_maps.back(), e, false);

    vertex_corner_id_maps.push_back( get(CGAL::dynamic_vertex_property_t<std::size_t>(), tm) );
    for(vertex_descriptor v : vertices(tm))
      put(vertex_corner_id_maps.back(), v, Planar_segmentation::init_id());

    face_cc_ids_maps.push_back( get(CGAL::dynamic_face_property_t<std::size_t>(), tm) );
    for(face_descriptor f : faces(tm))
      put(face_cc_ids_maps.back(), f, -1);
  }

  std::map<Point_3, std::map<std::size_t, vertex_descriptor> > point_to_vertex_maps;

  //start by detecting and marking all shared vertices
  std::size_t mesh_id = 0;
  for(Triangle_mesh* tm_ptr : mesh_ptrs)
  {
    Triangle_mesh& tm = *tm_ptr;

    for(vertex_descriptor v : vertices(tm))
    {
      std::map<std::size_t, vertex_descriptor>& mesh_id_to_vertex =
        point_to_vertex_maps[get(vpms[mesh_id], v)];
      if (!mesh_id_to_vertex.empty())
        put(vertex_shared_maps[mesh_id], v, true);
      if (mesh_id_to_vertex.size()==1)
      {
        std::pair<std::size_t, vertex_descriptor> other=*mesh_id_to_vertex.begin();
        put(vertex_shared_maps[other.first], other.second, true);
      }
      mesh_id_to_vertex.insert( std::make_pair(mesh_id, v) );
    }

    ++mesh_id;
  }
#ifndef CGAL_DO_NOT_USE_PCA
  if (use_PCA)
  {
    mesh_id=0;
    for(Triangle_mesh* tm_ptr : mesh_ptrs)
    {
      Triangle_mesh& tm = *tm_ptr;

      // mark constrained edges of coplanar regions detected with PCA
      segment_via_plane_fitting(tm, face_cc_ids_maps[mesh_id], parameters::vertex_point_map(vpms[mesh_id])
                                                                          .minimum_cosinus_squared(min_cosinus_squared)
                                                                          .maximum_Frechet_distance(max_frechet_distance)
                                                                          .edge_is_constrained_map(edge_is_constrained_maps[mesh_id]));
      ++mesh_id;
    }
    /// @TODO in this version there is no guarantee that an edge internal to a shared patch will
    ///       be constrained in all the meshes sharing the patch. I think this is a bug!
  }
#endif

  //then detect edge on the boundary of shared patches and mark them as constrained
  mark_boundary_of_shared_patches_as_constrained_edges(mesh_ptrs, point_to_vertex_maps, edge_is_constrained_maps, vertex_shared_maps, vpms);

  // first tag corners and constrained edges
  std::vector< std::pair<std::size_t, std::size_t> > nb_corners_and_nb_cc_all(nb_meshes);
  mesh_id=0;
  for(Triangle_mesh* tm_ptr : mesh_ptrs)
  {
    Triangle_mesh& tm = *tm_ptr;

    //reset face cc ids as it was set by segment_via_plane_fitting
    for(face_descriptor f : faces(tm))
      put(face_cc_ids_maps[mesh_id], f, -1);

    nb_corners_and_nb_cc_all[mesh_id] =
      tag_corners_and_constrained_edges(tm,
                                        min_cosinus_squared,
                                        33.0,
                                        false,
                                        vertex_corner_id_maps[mesh_id],
                                        edge_is_constrained_maps[mesh_id],
                                        face_cc_ids_maps[mesh_id],
                                        vpms[mesh_id]);
    ++mesh_id;
  }
#ifndef CGAL_DO_NOT_USE_PCA
  if (use_PCA)
  {
    mesh_id=0;
    for(Triangle_mesh* tm_ptr : mesh_ptrs)
    {
      Triangle_mesh& tm = *tm_ptr;

      nb_corners_and_nb_cc_all[mesh_id].first = mark_extra_corners_with_pca(tm,
                                                  max_frechet_distance,
                                                  nb_corners_and_nb_cc_all[mesh_id].first,
                                                  edge_is_constrained_maps[mesh_id],
                                                  vertex_corner_id_maps[mesh_id],
                                                  vpms[mesh_id]);
      ++mesh_id;
    }
  }
#endif

  // extra step to propagate is_corner to all meshes to make sure shared vertices are kept
  propagate_corner_status(vertex_corner_id_maps, point_to_vertex_maps, nb_corners_and_nb_cc_all);

  /// @TODO: make identical patches normal identical (up to the sign). Needed only in the approximate case

  // now call the decimation
  bool res=false;
  mesh_id=0;
  for(Triangle_mesh* tm_ptr : mesh_ptrs)
  {
    Triangle_mesh& tm = *tm_ptr;

    res = decimate_impl(tm,
                        nb_corners_and_nb_cc_all[mesh_id],
                        vertex_corner_id_maps[mesh_id],
                        edge_is_constrained_maps[mesh_id],
                        face_cc_ids_maps[mesh_id],
                        vpms[mesh_id]);
    if (!res) break;
    ++mesh_id;
  }

  return res;
}

} //end of namespace Planar_segmentation

template <typename TriangleMesh, typename NamedParameters>
bool decimate(TriangleMesh& tm, const NamedParameters& np)
{
  const double min_cosinus_squared =
    parameters::choose_parameter<double>(parameters::get_parameter(np, internal_np::min_cosinus_squared), 1.0);
  CGAL_assertion(min_cosinus_squared >= 0.0);
  const double max_frechet_distance =
    parameters::choose_parameter<double>(parameters::get_parameter(np, internal_np::max_Frechet_distance), 33.0);
  CGAL_assertion(max_frechet_distance >= 0.0);
  const bool use_region_growing =
    parameters::choose_parameter<bool>(parameters::get_parameter(np, internal_np::use_region_growing), true);

  typedef typename boost::graph_traits<TriangleMesh> graph_traits;
  typedef typename graph_traits::edge_descriptor edge_descriptor;
  typedef typename graph_traits::vertex_descriptor vertex_descriptor;
  typedef typename graph_traits::face_descriptor face_descriptor;

  // initialize property maps
  typename boost::property_map<TriangleMesh, CGAL::dynamic_edge_property_t<bool> >::type edge_is_constrained = get(CGAL::dynamic_edge_property_t<bool>(), tm);
  for(edge_descriptor e : edges(tm))
    put(edge_is_constrained, e, false);

  typename boost::property_map<TriangleMesh, CGAL::dynamic_vertex_property_t<std::size_t> >::type vertex_corner_id = get(CGAL::dynamic_vertex_property_t<std::size_t>(), tm);
  for(vertex_descriptor v : vertices(tm))
    put(vertex_corner_id, v, Planar_segmentation::init_id());

  typename boost::property_map<TriangleMesh, CGAL::dynamic_face_property_t<std::size_t> >::type face_cc_ids = get(CGAL::dynamic_face_property_t<std::size_t>(), tm);
  for(face_descriptor f : faces(tm))
    put(face_cc_ids, f, -1);

  /// @TODO turn into a named parameter
  typename boost::property_map<TriangleMesh, boost::vertex_point_t>::type vpm = get(boost::vertex_point, tm);
  std::pair<std::size_t, std::size_t> nb_corners_and_nb_cc =
    Planar_segmentation::tag_corners_and_constrained_edges(tm, min_cosinus_squared, max_frechet_distance, use_region_growing, vertex_corner_id, edge_is_constrained, face_cc_ids, vpm);
  bool res=Planar_segmentation::decimate_impl(tm,
                                              nb_corners_and_nb_cc,
                                              vertex_corner_id,
                                              edge_is_constrained,
                                              face_cc_ids,
                                              vpm);

  return res;
}

template <typename TriangleMesh>
bool decimate(TriangleMesh& tm, const double min_cosinus_squared = 1.0) {
  return decimate(tm, CGAL::parameters::minimum_cosinus_squared(min_cosinus_squared));
}

// MeshMap must be a mutable lvalue pmap with Triangle_mesh as value_type
template <typename TriangleMeshRange, typename MeshMap>
bool decimate_meshes_with_common_interfaces(TriangleMeshRange& meshes, double min_cosinus_squared, MeshMap mesh_map)
{
  CGAL_assertion(min_cosinus_squared>=0);
  typedef typename boost::property_traits<MeshMap>::value_type Triangle_mesh;
  typedef typename std::iterator_traits<typename TriangleMeshRange::iterator>::value_type Mesh_descriptor;

  /// @TODO turn into a range of named parameter
  std::vector<typename boost::property_map<Triangle_mesh, boost::vertex_point_t>::type > vpms;
  vpms.reserve(meshes.size());

  for(Mesh_descriptor& md : meshes)
    vpms.push_back( get(boost::vertex_point, mesh_map[md]) );
  return Planar_segmentation::decimate_meshes_with_common_interfaces_impl(meshes, mesh_map, 0, min_cosinus_squared, vpms);
}

template <class TriangleMesh>
bool decimate_meshes_with_common_interfaces(std::vector<TriangleMesh>& meshes, double min_cosinus_squared=1)
{
  return decimate_meshes_with_common_interfaces(meshes, min_cosinus_squared, CGAL::Identity_property_map<TriangleMesh>());
}

#ifndef CGAL_DO_NOT_USE_PCA
template <typename TriangleMeshRange, typename MeshMap>
bool decimate_meshes_with_common_interfaces_and_pca_for_coplanarity(
  TriangleMeshRange& meshes,
  double max_frechet_distance,
  double min_cosinus_squared,
  MeshMap mesh_map)
{
  CGAL_assertion(min_cosinus_squared>=0);
  typedef typename boost::property_traits<MeshMap>::value_type Triangle_mesh;
  typedef typename std::iterator_traits<typename TriangleMeshRange::iterator>::value_type Mesh_descriptor;

  /// @TODO turn into a range of named parameter
  std::vector<typename boost::property_map<Triangle_mesh, boost::vertex_point_t>::type > vpms;
  vpms.reserve(meshes.size());
  for(Mesh_descriptor& md : meshes)
    vpms.push_back( get(boost::vertex_point, mesh_map[md]) );
  return Planar_segmentation::decimate_meshes_with_common_interfaces_impl(meshes, mesh_map, max_frechet_distance, min_cosinus_squared, vpms);
}

template <typename TriangleMesh>
bool decimate_meshes_with_common_interfaces_and_pca_for_coplanarity(
  std::vector<TriangleMesh>& meshes,
  double max_frechet_distance,
  double min_cosinus_squared)
{
  return decimate_meshes_with_common_interfaces_and_pca_for_coplanarity(
    meshes, max_frechet_distance, min_cosinus_squared, CGAL::Identity_property_map<TriangleMesh>());
}

///@TODO use something better than a fitting score
///@TODO remove debug
template <typename TriangleMesh>
bool decimate_with_pca_for_coplanarity(TriangleMesh& tm,
                                       double max_frechet_distance,
                                       double min_cosinus_squared)
{
  typedef typename boost::graph_traits<TriangleMesh> graph_traits;
  typedef typename graph_traits::edge_descriptor edge_descriptor;
  typedef typename graph_traits::face_descriptor face_descriptor;
  typedef typename graph_traits::vertex_descriptor vertex_descriptor;

  ///@TODO turn it into a named parameter XXX
  typename boost::property_map<TriangleMesh, boost::vertex_point_t>::type vpm =
    get(boost::vertex_point, tm);

  CGAL_assertion(min_cosinus_squared>=0);
  // initialize property maps
  typename boost::property_map<TriangleMesh, CGAL::dynamic_edge_property_t<bool> >::type edge_is_constrained = get(CGAL::dynamic_edge_property_t<bool>(), tm);
  for(edge_descriptor e : edges(tm))
    put(edge_is_constrained, e, false);

  typename boost::property_map<TriangleMesh, CGAL::dynamic_vertex_property_t<std::size_t> >::type vertex_corner_id = get(CGAL::dynamic_vertex_property_t<std::size_t>(), tm);
  for(vertex_descriptor v : vertices(tm))
    put(vertex_corner_id, v, Planar_segmentation::init_id());

  typename boost::property_map<TriangleMesh, CGAL::dynamic_face_property_t<std::size_t> >::type face_cc_ids = get(CGAL::dynamic_face_property_t<std::size_t>(), tm);
  for(face_descriptor f : faces(tm))
    put(face_cc_ids, f, -1);

  std::size_t nb_cc = segment_via_plane_fitting(tm, face_cc_ids, parameters::vertex_point_map(vpm)
                                                                            .minimum_cosinus_squared(min_cosinus_squared)
                                                                            .maximum_Frechet_distance(max_frechet_distance)
                                                                            .edge_is_constrained_map(edge_is_constrained));

  // initial set of corner vertices
  std::size_t nb_corners =
    Planar_segmentation::mark_corner_vertices(tm, edge_is_constrained, vertex_corner_id, min_cosinus_squared, vpm);

#ifdef DEBUG_PCA
  std::cout << "found " << nb_cc << " components\n";
  std::ofstream tmp_out("/tmp/csts.cgal");
  for(edge_descriptor e : edges(tm))
    if (get(edge_is_constrained, e))
      tmp_out << "2 " << get(vpm, source(halfedge(e, tm), tm)) << " "
              << get(vpm, target(halfedge(e, tm), tm)) << "\n";
  tmp_out.close();
  std::cout << "  initial nb_corners: " << nb_corners << "\n";
#endif

  // apply pca also for patch boundaries: This will lead to the tagging of new corner vertices
  nb_corners = Planar_segmentation::mark_extra_corners_with_pca(tm, max_frechet_distance, nb_corners, edge_is_constrained, vertex_corner_id, vpm);

#ifdef DEBUG_PCA
  std::cout << "  nb_corners after constraint graph simplification: " << nb_corners << "\n";
#endif

  // now run the main decimation function
  bool res=Planar_segmentation::decimate_impl(tm, std::make_pair(nb_corners, nb_cc), vertex_corner_id, edge_is_constrained, face_cc_ids, vpm);

  return res;
}
#endif

} } // end of CGAL::Polygon_mesh_processing

#endif // CGAL_POLYGON_MESH_PROCESSING_PLANAR_SEGMENTATION_H
