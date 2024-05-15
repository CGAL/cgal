// Copyright (c) 2023 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sven Oesau, Florent Lafarge, Dmitry Anisimov, Simon Giraudot

#ifndef CGAL_KSP_3_DATA_STRUCTURE_H
#define CGAL_KSP_3_DATA_STRUCTURE_H

#include <CGAL/license/Kinetic_space_partition.h>

#include <boost/iterator/transform_iterator.hpp>
#include <boost/function.hpp>

// Internal includes.
#include <CGAL/KSP/utils.h>
#include <CGAL/KSP/debug.h>
#include <CGAL/KSP/parameters.h>

#include <CGAL/KSP_3/Support_plane.h>
#include <CGAL/KSP_3/Intersection_graph.h>

namespace CGAL {
namespace KSP_3 {
namespace internal {

#ifdef DOXYGEN_RUNNING
#else

template<typename GeomTraits, typename IntersectionKernel>
class Data_structure {

public:
  using Kernel = GeomTraits;
  using Intersection_kernel = IntersectionKernel;

  using Support_plane = CGAL::KSP_3::internal::Support_plane<Kernel, Intersection_kernel>;
  using Intersection_graph = CGAL::KSP_3::internal::Intersection_graph<Kernel, Intersection_kernel>;
  using Face_event = typename Support_plane::Face_event;

  using FT = typename Kernel::FT;
  using Point_2 = typename Kernel::Point_2;
  using IkPoint_2 = typename Intersection_kernel::Point_2;
  using Point_3 = typename Kernel::Point_3;
  using IkPoint_3 = typename Intersection_kernel::Point_3;
  using Segment_2 = typename Kernel::Segment_2;
  using IkSegment_2 = typename Intersection_kernel::Segment_2;
  using Segment_3 = typename Kernel::Segment_3;
  using IkSegment_3 = typename Intersection_kernel::Segment_3;
  using Vector_2 = typename Kernel::Vector_2;
  using Direction_2 = typename Kernel::Direction_2;
  using IkDirection_2 = typename Intersection_kernel::Direction_2;
  using Triangle_2  = typename Kernel::Triangle_2;
  using Line_2      = typename Kernel::Line_2;
  using IkLine_2    = typename Intersection_kernel::Line_2;
  using Plane_3     = typename Kernel::Plane_3;

  using Polygon_2  = CGAL::Polygon_2<Kernel>;
  using Parameters = KSP::internal::Parameters_3<FT>;

  using To_exact = CGAL::Cartesian_converter<Kernel, Intersection_kernel>;
  using From_exact = CGAL::Cartesian_converter<Intersection_kernel, Kernel>;

public:
  using Mesh           = typename Support_plane::Mesh;
  using Vertex_index   = typename Mesh::Vertex_index;
  using Face_index     = typename Mesh::Face_index;
  using Edge_index     = typename Mesh::Edge_index;
  using Halfedge_index = typename Mesh::Halfedge_index;

  using PVertex = std::pair<std::size_t, Vertex_index>;
  using PFace   = std::pair<std::size_t, Face_index>;
  using PEdge   = std::pair<std::size_t, Edge_index>;

  template<typename PSimplex>
  struct Make_PSimplex {
    using argument_type = typename PSimplex::second_type;
    using result_type   = PSimplex;

    const std::size_t support_plane_idx;
    Make_PSimplex(const std::size_t sp_idx) :
    support_plane_idx(sp_idx)
    { }

    const result_type operator()(const argument_type& arg) const {
      return result_type(support_plane_idx, arg);
    }
  };

  // ToDo:: check all kind of iterators/circulators
  using PVertex_iterator = boost::transform_iterator<Make_PSimplex<PVertex>, typename Mesh::Vertex_range::iterator>;
  using PVertices = CGAL::Iterator_range<PVertex_iterator>;

  using PFace_iterator = boost::transform_iterator<Make_PSimplex<PFace>, typename Mesh::Face_range::iterator>;
  using PFaces = CGAL::Iterator_range<PFace_iterator>;

  using PEdge_iterator = boost::transform_iterator<Make_PSimplex<PEdge>, typename Mesh::Edge_range::iterator>;
  using PEdges = CGAL::Iterator_range<PEdge_iterator>;

  struct Halfedge_to_pvertex {
    using argument_type = Halfedge_index;
    using result_type   = PVertex;

    const std::size_t support_plane_idx;
    const Mesh& mesh;

    Halfedge_to_pvertex(const std::size_t sp_idx, const Mesh& m) :
    support_plane_idx(sp_idx),
    mesh(m)
    { }

    const result_type operator()(const argument_type& arg) const {
      return result_type(support_plane_idx, mesh.target(arg));
    }
  };

  using PVertex_of_pface_iterator = boost::transform_iterator<Halfedge_to_pvertex, CGAL::Halfedge_around_face_iterator<Mesh> >;
  using PVertices_of_pface = CGAL::Iterator_range<PVertex_of_pface_iterator>;

  struct Halfedge_to_pedge {
    using argument_type = Halfedge_index;
    using result_type = PEdge;

    const std::size_t support_plane_idx;
    const Mesh& mesh;

    Halfedge_to_pedge(const std::size_t sp_idx, const Mesh& m) :
    support_plane_idx(sp_idx),
    mesh(m)
    { }

    const result_type operator()(const argument_type& arg) const {
      return result_type(support_plane_idx, mesh.edge(arg));
    }
  };

  struct Halfedge_to_pface {
    using argument_type = Halfedge_index;
    using result_type = PFace;

    const std::size_t support_plane_idx;
    const Mesh& mesh;

    Halfedge_to_pface(const std::size_t sp_idx, const Mesh& m) :
    support_plane_idx(sp_idx),
    mesh(m)
    { }

    const result_type operator()(const argument_type& arg) const {
      return result_type(support_plane_idx, mesh.face(arg));
    }
  };

  // ToDo:: check all kind of iterators/circulators
  using PEdge_around_pvertex_iterator = boost::transform_iterator<Halfedge_to_pedge, CGAL::Halfedge_around_target_iterator<Mesh> >;
  using PEdges_around_pvertex = CGAL::Iterator_range<PEdge_around_pvertex_iterator>;

  using PEdge_of_pface_iterator = boost::transform_iterator<Halfedge_to_pedge, CGAL::Halfedge_around_face_iterator<Mesh> >;
  using PEdges_of_pface = CGAL::Iterator_range<PEdge_of_pface_iterator>;

  using PFace_around_pvertex_iterator = boost::transform_iterator<Halfedge_to_pface, CGAL::Halfedge_around_target_iterator<Mesh> >;
  using PFaces_around_pvertex = CGAL::Iterator_range<PFace_around_pvertex_iterator>;

  using IVertex = typename Intersection_graph::Vertex_descriptor;
  using IEdge   = typename Intersection_graph::Edge_descriptor;
  using IFace   = typename Intersection_graph::Face_descriptor;

  using IEdge_set = typename Intersection_graph::IEdge_set;

  using Index = std::pair<std::size_t, std::size_t>; // first is partition index, second is volume index. Partition -1 indicates outside. 0-5 indicate which side of the bbox

  struct Volume_cell {
    std::vector<PFace> pfaces;// index compatible to faces and neighbors
    std::vector<size_t> faces;// Indices into m_face2vertices in m_data.
    std::vector<bool> pface_oriented_outwards;
    std::vector<int> neighbors;
    std::set<PVertex> pvertices;
    std::size_t index = std::size_t(-1);
    Point_3 centroid;

    FT inside  = FT(1);
    FT outside = FT(0);
    FT weight  = FT(0);

    FT inside_count = FT(0);
    FT outside_count = FT(0);

    void add_pface(const PFace& pface, const int neighbor) {
      pfaces.push_back(pface);
      neighbors.push_back(neighbor);
    }
    void set_index(const std::size_t idx) {
      index = idx;
    }
    void set_centroid(const Point_3& point) {
      centroid = point;
    }
  };

  struct Reconstructed_model {
    std::vector<PFace> pfaces;
    void clear() {
      pfaces.clear();
    }
  };

private:
  std::vector<Support_plane> m_support_planes;
  std::vector<typename Support_plane::Data> m_initial_support_planes;
  Intersection_graph m_intersection_graph;

  std::vector<std::vector<Point_3> > m_input_polygons;

  To_exact to_exact;
  From_exact from_exact;

  const Parameters& m_parameters;

  std::string m_prefix;

  std::vector<Volume_cell> m_volumes;
  std::vector<Point_3> m_vertices;
  std::vector<typename Intersection_kernel::Point_3> m_exact_vertices;
  std::vector<int> m_ivertex2vertex; // Used to map ivertices to m_vertices which only contain vertices of the finalized kinetic partition.
  std::vector<std::pair<int, int> > m_face2volumes;

  std::map<PFace, std::size_t> m_face2index;
  std::vector<std::vector<std::size_t> > m_face2vertices;
  std::vector<bool> m_part_of_initial;
  std::map<PFace, std::pair<int, int> > m_pface_neighbors;
  std::vector<std::size_t> m_face2sp;
  std::vector<std::set<std::size_t> > m_sp2input_polygon;
  std::map<std::size_t, std::size_t> m_input_polygon_map; // Maps index of input polygon onto support plane indices.
  Reconstructed_model m_reconstructed_model;

public:
  Data_structure(const Parameters& parameters, const std::string &prefix) : to_exact(), from_exact(), m_parameters(parameters), m_prefix(prefix) { }

  template<typename Type1, typename Type2, typename ResultType>
  static bool intersection(const Type1& t1, const Type2& t2, ResultType& result) {
    const auto inter = CGAL::intersection(t1, t2);
    if (!inter) return false;
    if (CGAL::assign(result, inter))
      return true;

    return false;
  }

  /*******************************
  **      INITIALIZATION        **
  ********************************/

  void clear() {
    m_support_planes.clear();
    m_intersection_graph.clear();

    m_volumes.clear();
    m_pface_neighbors.clear();
    m_input_polygon_map.clear();
    m_reconstructed_model.clear();
  }

  void precompute_iedge_data() {
    for (std::size_t i = 0; i < number_of_support_planes(); ++i) {
      auto& unique_iedges = support_plane(i).unique_iedges();
      CGAL_assertion(unique_iedges.size() > 0);

      auto& iedges = this->iedges(i);

      iedges.clear();
      iedges.reserve(unique_iedges.size());
      std::copy(unique_iedges.begin(), unique_iedges.end(), std::back_inserter(iedges));
      unique_iedges.clear();
    }
  }

  void initialization_done() {
    m_intersection_graph.initialization_done();
    m_initial_support_planes.resize(m_support_planes.size());
    for (std::size_t i = 0; i < m_support_planes.size(); i++)
      m_initial_support_planes[i] = m_support_planes[i].data();

  }
  void reset_to_initialization() {
    m_intersection_graph.reset_to_initialization();

    CGAL_assertion(m_support_planes.size() == m_initial_support_planes.size());
    for (std::size_t i = 0; i < m_support_planes.size(); i++) {
      m_support_planes[i].data() = m_initial_support_planes[i];

      m_support_planes[i].link_property_maps();
    }

    m_volumes.clear();
    m_vertices.clear();
    m_exact_vertices.clear();
    m_ivertex2vertex.clear();
    m_face2index.clear();
    m_face2vertices.clear();
    m_pface_neighbors.clear();
    m_face2sp.clear();
  }

  /*******************************
  **           ACCESS           **
  ********************************/

  const std::vector<std::vector<Point_3> >& input_polygons() const {
    return m_input_polygons;
  }

  int support_plane_index(const std::size_t polygon_index) const {
    CGAL_assertion(m_input_polygon_map.find(polygon_index) != m_input_polygon_map.end());
    const std::size_t sp_idx = m_input_polygon_map.at(polygon_index);
    return static_cast<int>(sp_idx);
  }

  std::size_t number_of_volumes() const {
    return m_volumes.size();
  }

  /*******************************
  **          GENERAL           **
  ********************************/

  std::map<PFace, std::pair<int, int> >& pface_neighbors() { return m_pface_neighbors; }
  const std::map<PFace, std::pair<int, int> >& pface_neighbors() const { return m_pface_neighbors; }
  std::map<PFace, std::size_t> &face_to_index() { return m_face2index; }
  std::vector<int>& ivertex_to_index() {
    if (m_ivertex2vertex.size() == 0)
      m_ivertex2vertex.resize(m_intersection_graph.number_of_vertices(), -1);

    return m_ivertex2vertex;
  }
  std::vector<std::pair<int, int> >& face_to_volumes() { return m_face2volumes; }
  const std::vector<std::pair<int, int> >& face_to_volumes() const { return m_face2volumes; }
  std::vector<Point_3>& vertices() { return m_vertices; }
  const std::vector<Point_3>& vertices() const { return m_vertices; }
  std::vector<typename Intersection_kernel::Point_3>& exact_vertices() { return m_exact_vertices; }
  const std::vector<typename Intersection_kernel::Point_3>& exact_vertices() const { return m_exact_vertices; }
  std::vector<std::vector<std::size_t> >& face_to_vertices() { return m_face2vertices; }
  const std::vector<std::vector<std::size_t> >& face_to_vertices() const { return m_face2vertices; }

  std::vector<bool>& face_is_part_of_input_polygon() { return m_part_of_initial; }
  const std::vector<bool>& face_is_part_of_input_polygon() const { return m_part_of_initial; }

  std::vector<std::size_t>& face_to_support_plane() { return m_face2sp; }
  std::vector<std::set<std::size_t> >& support_plane_to_input_polygon() { return m_sp2input_polygon; }

  const std::vector<Support_plane>& support_planes() const { return m_support_planes; }
  std::vector<Support_plane>& support_planes() { return m_support_planes; }

  const Intersection_graph& igraph() const { return m_intersection_graph; }
  Intersection_graph& igraph() { return m_intersection_graph; }

  const std::string& prefix() const { return m_prefix; }

  void resize(const std::size_t number_of_items) {
    m_support_planes.resize(number_of_items);
  }

  FT calculate_edge_intersection_time(std::size_t sp_idx, IEdge edge, Face_event &event) {
    // Not need to calculate for border edges.
    if (m_intersection_graph.iedge_is_on_bbox(edge))
      return 0;

    Support_plane& sp = m_support_planes[sp_idx];

    Point_2 centroid = sp.data().centroid;

    typename Intersection_graph::Kinetic_interval& kinetic_interval = m_intersection_graph.kinetic_interval(edge, sp_idx);

    if (kinetic_interval.size() != 0) {
      int a;
      a = 3;
    }

    Point_2 s = sp.to_2d(from_exact(point_3(m_intersection_graph.source(edge))));
    Point_2 t = sp.to_2d(from_exact(point_3(m_intersection_graph.target(edge))));
    Vector_2 segment = t - s;
    FT segment_length = sqrt(segment * segment);
    CGAL_assertion(segment_length > 0);
    segment = segment / segment_length;
    Direction_2 to_source(s - centroid);
    Direction_2 to_target(t - centroid);

    std::size_t source_idx = -1;
    std::size_t target_idx = -1;

    event.crossed_edge = edge;
    event.support_plane = sp_idx;

    std::pair<IFace, IFace> faces;
    m_intersection_graph.get_faces(sp_idx, edge, faces);

    if (m_intersection_graph.face(faces.first).part_of_partition)
      event.face = faces.second;
    else
      event.face = faces.first;

    for (std::size_t i = 0; i < sp.data().original_directions.size(); i++) {
      if (source_idx == -1 && sp.data().original_directions[i] > to_source)
        source_idx = i;

      if (target_idx == -1 && sp.data().original_directions[i] > to_target)
        target_idx = i;
    }

    source_idx = (source_idx == -1) ? 0 : source_idx;
    target_idx = (target_idx == -1) ? 0 : target_idx;

    std::size_t num;

    Vector_2 tt = to_target.vector(), ts = to_source.vector();
    bool ccw = (tt.x() * ts.y() - tt.y() * ts.x()) < 0;

    // Check whether the segment is cw or ccw oriented.
    if (!ccw) {
      std::size_t tmp = source_idx;
      source_idx = target_idx;
      target_idx = tmp;

      Point_2 tmp_p = s;
      s = t;
      t = tmp_p;
    }

    if (source_idx <= target_idx)
      num = target_idx - source_idx;
    else
      num = (sp.data().original_directions.size() + target_idx - source_idx);

    std::vector<FT> time(num);
    std::vector<Point_2> intersections(num);
    std::vector<FT> intersections_bary(num);

    // Shooting rays to find intersection with line of IEdge
    typename Intersection_kernel::Line_2 l = sp.to_2d(m_intersection_graph.line_3(edge));
    for (std::size_t i = 0; i < num; i++) {
      std::size_t idx = (i + source_idx) % sp.data().original_directions.size();
      const auto result = CGAL::intersection(l, sp.data().original_rays[idx]);
      if (!result) {
        time[i] = (std::numeric_limits<double>::max)();
        continue;
      }
      IkPoint_2 p;
      if (CGAL::assign(p, result)) {
        FT l = CGAL::sqrt(sp.data().original_vectors[idx].squared_length());

        double l2 = CGAL::to_double((p - sp.data().original_rays[idx].point(0)).squared_length());
        time[i] = l2 / l;
        CGAL_assertion(0 <= time[i]);
        intersections[i] = from_exact(p);
        intersections_bary[i] = abs(((from_exact(p) - s) * segment)) / segment_length;
        if (!ccw)
          intersections_bary[i] = 1.0 - intersections_bary[i];
      }
      // If the intersection is a segment, it can be safely ignored as there are also two intersections with the adjacent edges.
    }

    // Calculate pedge vs ivertex collision
    FT edge_time[2];

    // Source edge time
    std::size_t adjacent = (source_idx + sp.data().original_vertices.size() - 1) % sp.data().original_vertices.size();
    Vector_2 dir = sp.data().original_vertices[source_idx] - sp.data().original_vertices[adjacent];
    dir = dir / CGAL::sqrt(dir * dir);

    // Orthogonal direction matching the direction of the adjacent vertices
    dir = Vector_2(dir.y(), -dir.x());

    // Moving speed matches the speed of adjacent vertices
    FT speed = (dir * sp.data().original_vectors[source_idx]);

    if (speed < 0)
      speed = -speed;

    // Distance from edge to endpoint of iedge
    FT dist = (s - sp.data().original_vertices[source_idx]) * dir;

    edge_time[0] = dist / speed;

    // Target edge time
    adjacent = (target_idx + sp.data().original_vertices.size() - 1) % sp.data().original_vertices.size();
    dir = sp.data().original_vertices[target_idx] - sp.data().original_vertices[adjacent];
    dir = dir / CGAL::sqrt(dir * dir);

    // Orthogonal direction matching the direction of the adjacent vertices
    dir = Vector_2(dir.y(), -dir.x());

    // Moving speed matches the speed of adjacent vertices
    speed = (dir * sp.data().original_vectors[target_idx]);

    if (speed < 0)
      speed = -speed;

    // Distance from edge to endpoint of iedge
    dist = (t - sp.data().original_vertices[target_idx]) * dir;

    edge_time[1] = dist / speed;

    // Fill event structure and kinetic intervals.
    if (ccw)
      kinetic_interval.push_back(std::pair<FT, FT>(0, edge_time[0]));
    else
      kinetic_interval.push_back(std::pair<FT, FT>(0, edge_time[1]));

    event.time = kinetic_interval.back().second;
    event.intersection_bary = 0;

    for (std::size_t i = 0; i < num; i++) {
      kinetic_interval.push_back(std::pair<FT, FT>(intersections_bary[i], time[i]));
      if (event.time > time[i]) {
        event.time = time[i];
        event.intersection_bary = intersections_bary[i];
      }
    }

    if (ccw)
      kinetic_interval.push_back(std::pair<FT, FT>(1, edge_time[1]));
    else
      kinetic_interval.push_back(std::pair<FT, FT>(1, edge_time[0]));

    if (event.time > kinetic_interval.back().second) {
      event.time = kinetic_interval.back().second;
      event.intersection_bary = 1;
    }

    if (kinetic_interval.size() > 4) {
      if (kinetic_interval[2].first > kinetic_interval[1].first) {
        int a;
        a = 2;
      }
    }

    CGAL_assertion(0 <= event.intersection_bary && event.intersection_bary <= 1);

    return event.time;
  }

  template<typename Queue>
  void fill_event_queue(Queue& queue) {
    for (std::size_t sp_idx = 6; sp_idx < m_support_planes.size(); sp_idx++) {
      std::vector<IEdge> border;
      m_support_planes[sp_idx].get_border(m_intersection_graph, border);

      for (IEdge edge : border) {
        if (m_intersection_graph.has_crossed(edge, sp_idx))
          continue;

        Face_event fe;
        FT t = calculate_edge_intersection_time(sp_idx, edge, fe);
        if (t > 0)
          queue.push(fe);
      }
    }
  }

  std::vector<Volume_cell>& volumes() { return m_volumes; }
  const std::vector<Volume_cell>& volumes() const { return m_volumes; }

  const std::vector<std::size_t>& volume(std::size_t volume_index) {
    return m_volumes[volume_index].faces;
  }

  const std::vector<std::size_t> &face(std::size_t face_index) const { return m_face2vertices[face_index]; }
  const Point_3& vertex(std::size_t vertex_index) const { return m_vertices[vertex_index]; }

  Reconstructed_model& reconstructed_model() { return m_reconstructed_model; }
  const Reconstructed_model& reconstructed_model() const { return m_reconstructed_model; }

  /*******************************
  **      SUPPORT PLANES        **
  ********************************/

  template<typename PSimplex>
  const Support_plane& support_plane(const PSimplex& psimplex) const { return support_plane(psimplex.first); }
  const Support_plane& support_plane(const std::size_t idx) const { return m_support_planes[idx]; }

  template<typename PSimplex>
  Support_plane& support_plane(const PSimplex& psimplex) { return support_plane(psimplex.first); }
  Support_plane& support_plane(const std::size_t idx) { return m_support_planes[idx]; }

  template<typename PSimplex>
  const Mesh& mesh(const PSimplex& psimplex) const { return mesh(psimplex.first); }
  const Mesh& mesh(const std::size_t support_plane_idx) const { return support_plane(support_plane_idx).mesh(); }

  template<typename PSimplex>
  Mesh& mesh(const PSimplex& psimplex) { return mesh(psimplex.first); }
  Mesh& mesh(const std::size_t support_plane_idx) { return support_plane(support_plane_idx).mesh(); }

  std::size_t number_of_support_planes() const {
    return m_support_planes.size();
  }

  bool is_bbox_support_plane(const std::size_t support_plane_idx) const {
    return (support_plane_idx < 6);
  }

  template<typename PointRange>
  std::pair<std::size_t, bool> add_support_plane(const PointRange& polygon, const bool is_bbox, const typename Intersection_kernel::Plane_3& plane) {
    const Support_plane new_support_plane(polygon, is_bbox, plane);
    std::size_t support_plane_idx = std::size_t(-1);

    for (std::size_t i = 0; i < number_of_support_planes(); ++i) {
      if (new_support_plane == support_plane(i)) {
        support_plane_idx = i;
        return std::make_pair(support_plane_idx, false);
      }
    }

    if (support_plane_idx == std::size_t(-1)) {
      support_plane_idx = number_of_support_planes();
      m_support_planes.push_back(new_support_plane);
    }

    intersect_with_bbox(support_plane_idx);

    if (m_sp2input_polygon.size() <= number_of_support_planes())
      m_sp2input_polygon.resize(number_of_support_planes());

    return std::make_pair(support_plane_idx, true);
  }

  template<typename PointRange>
  std::pair<std::size_t, bool> add_support_plane(const PointRange& polygon, const bool is_bbox) {
    const Support_plane new_support_plane(polygon, is_bbox);
    std::size_t support_plane_idx = std::size_t(-1);

    for (std::size_t i = 0; i < number_of_support_planes(); ++i) {
      if (new_support_plane == support_plane(i)) {
        support_plane_idx = i;
        return std::make_pair(support_plane_idx, false);
      }
    }

    if (support_plane_idx == std::size_t(-1)) {
      support_plane_idx = number_of_support_planes();
      m_support_planes.push_back(new_support_plane);
    }

    intersect_with_bbox(support_plane_idx);

    if (m_sp2input_polygon.size() <= number_of_support_planes())
      m_sp2input_polygon.resize(number_of_support_planes());

    return std::make_pair(support_plane_idx, true);
  }

  void intersect_with_bbox(const std::size_t sp_idx) {
    if (is_bbox_support_plane(sp_idx)) return;

    typename Intersection_kernel::FT bbox_center_x = 0, bbox_center_y = 0, bbox_center_z = 0;
    for (std::size_t i = 0; i < 8; i++) {
      IkPoint_3 tmp = point_3(IVertex(i));
      bbox_center_x += tmp.x();
      bbox_center_y += tmp.y();
      bbox_center_z += tmp.z();
    }

    IkPoint_3 bbox_center(bbox_center_x * 0.125, bbox_center_y * 0.125, bbox_center_z * 0.125);

    // Intersect current plane with all bbox iedges.
    IkPoint_3 point;
    Point_3 p1;
    const auto& sp = support_plane(sp_idx);
    const auto& plane = sp.exact_plane();

    using IEdge_vec = std::vector<IEdge>;
    using IPair = std::pair<IVertex, IEdge_vec>;
    using Pair = std::pair<IkPoint_3, IPair>;

    std::vector<Pair> polygon;
    polygon.reserve(3);
    const auto all_iedges = m_intersection_graph.edges();
    for (const auto iedge : all_iedges) {
      const auto segment = segment_3(iedge);

      if (!intersection(plane, segment, point))
        continue;

      const auto isource = source(iedge);
      const auto itarget = target(iedge);
      const bool is_isource = KSP::internal::distance(point, point_3(isource)) == 0;
      const bool is_itarget = KSP::internal::distance(point, point_3(itarget)) == 0;

      std::vector<IEdge> iedges;
      if (is_isource) {
        CGAL_assertion(!is_itarget);
        const auto inc_edges = m_intersection_graph.incident_edges(isource);
        CGAL_assertion(iedges.size() == 0);
        iedges.reserve(inc_edges.size());
        std::copy(inc_edges.begin(), inc_edges.end(), std::back_inserter(iedges));
        CGAL_assertion(iedges.size() == inc_edges.size());
        polygon.push_back(std::make_pair(point, std::make_pair(isource, iedges)));
      }

      if (is_itarget) {
        CGAL_assertion(!is_isource);
        const auto inc_edges = m_intersection_graph.incident_edges(itarget);
        CGAL_assertion(iedges.size() == 0);
        iedges.reserve(inc_edges.size());
        std::copy(inc_edges.begin(), inc_edges.end(), std::back_inserter(iedges));
        CGAL_assertion(iedges.size() == inc_edges.size());
        polygon.push_back(std::make_pair(point, std::make_pair(itarget, iedges)));
      }

      if (!is_isource && !is_itarget) {
        CGAL_assertion(iedges.size() == 0);
        iedges.push_back(iedge);
        polygon.push_back(std::make_pair(point, std::make_pair(null_ivertex(), iedges)));
      }
    }

    // Sort the points to get an oriented polygon.
    typename Intersection_kernel::FT x = 0, y = 0, z = 0, f = 1.0 / polygon.size();
    for (const auto& p : polygon) {
      x += p.first.x();
      y += p.first.y();
      z += p.first.z();
    }
    x *= f;
    y *= f;
    z *= f;
    IkPoint_2 mid = sp.to_2d(IkPoint_3(x, y, z));

    std::sort(polygon.begin(), polygon.end(),
    [&](const Pair& a, const Pair& b) {
      const auto a2 = sp.to_2d(a.first);
      const auto b2 = sp.to_2d(b.first);
      const IkSegment_2 sega(mid, a2);
      const IkSegment_2 segb(mid, b2);
      return (IkDirection_2(sega) < IkDirection_2(segb));
    });

    remove_equal_points(polygon, 0);

    CGAL_assertion(is_valid_polygon(sp_idx, polygon));

    // Find common planes.
    std::vector<IVertex> vertices;
    std::vector<std::size_t> common_bbox_planes_idx;
    std::map<std::size_t, std::size_t> map_lines_idx; // Maps edges between vertices to line_idx

    const std::size_t n = polygon.size();
    common_bbox_planes_idx.reserve(n);
    vertices.reserve(n);

    std::vector< std::set<std::size_t> > all_iplanes;
    all_iplanes.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
      const auto& item = polygon[i].second;
      const auto& iedges = item.second;

      std::set<std::size_t> iplanes;
      for (const auto& iedge : iedges) {
        const auto& planes = m_intersection_graph.intersected_planes(iedge);
        iplanes.insert(planes.begin(), planes.end());
      }

      CGAL_assertion(iplanes.size() >= 2);
      all_iplanes.push_back(iplanes);
    }
    CGAL_assertion(all_iplanes.size() == n);

    for (std::size_t i = 0; i < n; ++i) {
      const std::size_t ip = (i + 1) % n;
      const auto& item = polygon[i].second;

      const auto& iplanes0 = all_iplanes[i];
      const auto& iplanes1 = all_iplanes[ip];

      std::size_t common_bbox_plane_idx = std::size_t(-1);
      bool dump = false;
      const std::function<void(const std::size_t& idx)> lambda =
        [&](const std::size_t& idx) {
        if (idx < 6) {

          if (common_bbox_plane_idx != std::size_t(-1))
            dump = true;
            common_bbox_plane_idx = idx;
          }
        };

      std::set_intersection(
        iplanes0.begin(), iplanes0.end(),
        iplanes1.begin(), iplanes1.end(),
        boost::make_function_output_iterator(lambda)
      );

      if (dump) {
        From_exact from_exact;

        std::ofstream vout("bug.polylines.txt");
        vout.precision(20);
        vout << "2 ";
        vout << " " << from_exact(polygon[i].first);
        vout << " " << from_exact(polygon[ip].first);
        vout << std::endl;
        vout.close();
      }

      CGAL_assertion(common_bbox_plane_idx != std::size_t(-1));
      common_bbox_planes_idx.push_back(common_bbox_plane_idx);

      const auto pair = map_lines_idx.insert(
        std::make_pair(common_bbox_plane_idx, std::size_t(-1)));
      const bool is_inserted = pair.second;
      if (is_inserted) {
        typename Intersection_kernel::Line_3 line;
        bool intersect = intersection(plane, m_support_planes[common_bbox_plane_idx].exact_plane(), line);
        CGAL_assertion(intersect);
        pair.first->second = m_intersection_graph.add_line(line);
      }

      if (item.first != null_ivertex()) {
        vertices.push_back(item.first);
      } else {
        CGAL_assertion(item.first == null_ivertex());
        vertices.push_back(
          m_intersection_graph.add_vertex(polygon[i].first).first);
      }
    }
    CGAL_assertion(common_bbox_planes_idx.size() == n);
    CGAL_assertion(vertices.size() == n);

    // Insert, split iedges.
    for (std::size_t i = 0; i < n; ++i) {
      const std::size_t ip = (i + 1) % n;
      const auto& item = polygon[i].second;

      const std::size_t common_bbox_plane_idx = common_bbox_planes_idx[i];
      const auto new_iedge = m_intersection_graph.add_edge(vertices[i], vertices[ip], sp_idx).first;

      m_intersection_graph.intersected_planes(new_iedge).insert(common_bbox_plane_idx);
      CGAL_assertion(map_lines_idx.find(common_bbox_plane_idx) != map_lines_idx.end());
      m_intersection_graph.set_line(new_iedge, map_lines_idx.at(common_bbox_plane_idx));
      support_plane(sp_idx).unique_iedges().insert(new_iedge);
      support_plane(common_bbox_plane_idx).unique_iedges().insert(new_iedge);

      // No further treatment necessary for exact intersections at vertices of edges.

      if (item.first == null_ivertex()) { // edge case, split
        const auto& iplanes = all_iplanes[i];
        CGAL_assertion(iplanes.size() >= 2);
        CGAL_assertion(item.second.size() == 1);
        const auto& iedge = item.second[0];
        for (const std::size_t plane_idx : iplanes) {
          support_plane(plane_idx).unique_iedges().erase(iedge);
        }
        const auto edges = m_intersection_graph.split_edge(iedge, vertices[i]);

        const auto& iplanes0 = m_intersection_graph.intersected_planes(edges.first);
        for (const std::size_t plane_idx : iplanes0) {
          support_plane(plane_idx).unique_iedges().insert(edges.first);
        }

        const auto& iplanes1 = m_intersection_graph.intersected_planes(edges.second);
        for (const std::size_t plane_idx : iplanes1) {
          support_plane(plane_idx).unique_iedges().insert(edges.second);
        }
      }
    }
  }

  template<typename PointRange>
  void add_bbox_polygon(const PointRange& polygon) {
    bool is_added = true;
    std::size_t support_plane_idx = std::size_t(-1);
    std::tie(support_plane_idx, is_added) = add_support_plane(polygon, true);
    CGAL_assertion(is_added);
    CGAL_assertion(support_plane_idx != std::size_t(-1));

    std::array<IVertex, 4> ivertices;
    std::array<Point_2, 4> points;
    for (std::size_t i = 0; i < 4; ++i) {
      points[i] = from_exact(support_plane(support_plane_idx).to_2d(polygon[i]));
      ivertices[i] = m_intersection_graph.add_vertex(polygon[i]).first;
    }

    const auto vertices =
      support_plane(support_plane_idx).add_bbox_polygon(points, ivertices);

    for (std::size_t i = 0; i < 4; ++i) {
      const auto pair = m_intersection_graph.add_edge(ivertices[i], ivertices[(i + 1) % 4], support_plane_idx);
      const auto& iedge = pair.first;
      const bool is_inserted = pair.second;
      if (is_inserted) {
        typename Intersection_kernel::Line_3 line(polygon[i], polygon[(i + 1) % 4]);
        m_intersection_graph.set_line(iedge, m_intersection_graph.add_line(line));
      }

      support_plane(support_plane_idx).set_iedge(vertices[i], vertices[(i + 1) % 4], iedge);
      support_plane(support_plane_idx).unique_iedges().insert(iedge);
    }
  }

  void add_input_polygon( const std::size_t support_plane_idx, const std::vector<std::size_t>& input_indices, const std::vector<Point_2>& polygon) {
    std::vector< std::pair<Point_2, bool> > points;
    points.reserve(polygon.size());
    for (const auto& point : polygon) {
      points.push_back(std::make_pair(point, true));
    }
    CGAL_assertion(points.size() == polygon.size());

    preprocess(points);
    sort_points_by_direction(points);
    support_plane(support_plane_idx).
      add_input_polygon(points, input_indices);
    for (const std::size_t input_index : input_indices) {
      m_input_polygon_map[input_index] = support_plane_idx;
      m_sp2input_polygon[support_plane_idx].insert(input_index);
    }
  }

  template<typename Pair>
  void preprocess(std::vector<Pair>& points, const FT min_dist = 0, const FT min_angle = FT(10)) const {
    remove_equal_points(points, min_dist);

    remove_collinear_points(points, min_angle);
  }

  template<typename Pair>
  void remove_equal_points(std::vector<Pair>& points, const FT min_dist) const {
    std::vector<Pair> polygon;
    const std::size_t n = points.size();
    for (std::size_t i = 0; i < n; ++i) {
      const auto& first = points[i];
      polygon.push_back(first);

      while (true) {
        const auto& p = points[i].first;
        const std::size_t ip = (i + 1) % n;
        const auto& q = points[ip].first;
        const FT distance = from_exact(KSP::internal::distance(p, q));
        const bool is_small = (distance <= min_dist);
        if (ip == 0 && is_small) break;
        if (is_small) {
          CGAL_assertion(ip != 0);
          i = ip; continue;
        }
        CGAL_assertion(!is_small);
        break;
      };
    }
    CGAL_assertion(polygon.size() >= 3);
    points = polygon;
  }

  template<typename Pair>
  void remove_collinear_points(std::vector<Pair>& points, const FT min_angle) const {
    std::vector<Pair> polygon;
    const std::size_t n = points.size();
    for (std::size_t i = 0; i < n; ++i) {
      const std::size_t im = (i + n - 1) % n;
      const std::size_t ip = (i + 1) % n;

      const auto& p = points[im].first;
      const auto& q = points[i].first;
      const auto& r = points[ip].first;

      Vector_2 vec1(q, r);
      Vector_2 vec2(q, p);
      vec1 = KSP::internal::normalize(vec1);
      vec2 = KSP::internal::normalize(vec2);

      const Direction_2 dir1(vec1);
      const Direction_2 dir2(vec2);
      const FT angle = KSP::internal::angle_2(dir1, dir2);

      if (angle > min_angle) polygon.push_back(points[i]);
    }
    if (polygon.size() >= 3) points = polygon;
    else remove_collinear_points(points, min_angle / FT(2));
  }

  template<typename Pair>
  void sort_points_by_direction(std::vector<Pair>& points) const {
    FT x = FT(0), y = FT(0);
    for (const auto& pair : points) {
      const auto& point = pair.first;
      x += point.x();
      y += point.y();
    }
    x /= static_cast<FT>(points.size());
    y /= static_cast<FT>(points.size());
    const Point_2 mid(x, y);

    std::sort(points.begin(), points.end(),
    [&](const Pair& a, const Pair& b) -> bool {
      const Segment_2 sega(mid, a.first);
      const Segment_2 segb(mid, b.first);
      return ( Direction_2(sega) < Direction_2(segb) );
    });
  }

  /*******************************
  **        PSimplices          **
  ********************************/

  const PVertices pvertices(const std::size_t support_plane_idx) const {
    return PVertices(
      boost::make_transform_iterator(mesh(support_plane_idx).vertices().begin(), Make_PSimplex<PVertex>(support_plane_idx)),
      boost::make_transform_iterator(mesh(support_plane_idx).vertices().end(), Make_PSimplex<PVertex>(support_plane_idx)));
  }

  const PEdges pedges(const std::size_t support_plane_idx) const {
    return PEdges(
      boost::make_transform_iterator(mesh(support_plane_idx).edges().begin(), Make_PSimplex<PEdge>(support_plane_idx)),
      boost::make_transform_iterator(mesh(support_plane_idx).edges().end(), Make_PSimplex<PEdge>(support_plane_idx)));
  }

  const PFaces pfaces(const std::size_t support_plane_idx) const {
    return PFaces(
      boost::make_transform_iterator(mesh(support_plane_idx).faces().begin(), Make_PSimplex<PFace>(support_plane_idx)),
      boost::make_transform_iterator(mesh(support_plane_idx).faces().end(), Make_PSimplex<PFace>(support_plane_idx)));
  }

  // Get prev and next pvertices of the free pvertex.
  const PVertex prev(const PVertex& pvertex) const {
    return PVertex(pvertex.first, support_plane(pvertex).prev(pvertex.second));
  }
  const PVertex next(const PVertex& pvertex) const {
    return PVertex(pvertex.first, support_plane(pvertex).next(pvertex.second));
  }

  void get_and_sort_all_connected_iedges(const std::size_t sp_idx, const IVertex& ivertex, std::vector< std::pair<IEdge, Direction_2> >& iedges) const {
    auto inc_iedges = incident_iedges(ivertex);
    const std::function<void(const IEdge& inc_iedge)> lambda =
      [&](const IEdge& inc_iedge) {
        const auto iplanes = intersected_planes(inc_iedge);
        if (iplanes.find(sp_idx) == iplanes.end()) {
          return;
        }
        const Direction_2 direction(
          point_2(sp_idx, opposite(inc_iedge, ivertex)) -
          point_2(sp_idx, ivertex));
        iedges.push_back(std::make_pair(inc_iedge, direction));
      };

    std::copy(
      inc_iedges.begin(), inc_iedges.end(),
      boost::make_function_output_iterator(lambda)
    );

    std::sort(iedges.begin(), iedges.end(),
      [&](const std::pair<IEdge, Direction_2>& a,
          const std::pair<IEdge, Direction_2>& b) -> bool {
        return a.second < b.second;
      }
    );
    CGAL_assertion(iedges.size() > 0);
  }

  const IFace add_iface(std::size_t support_plane) {
    return m_intersection_graph.add_face(support_plane);;
  }

  const PFace add_iface_to_mesh(std::size_t support_plane, IFace f_idx) {
    typename Intersection_graph::Face_property &f = m_intersection_graph.face(f_idx);

    std::vector<Vertex_index> vertices;
    vertices.reserve(f.vertices.size());
    Support_plane& sp = m_support_planes[support_plane];

    for (auto v : f.vertices) {
      auto& m = sp.ivertex2pvertex();
      std::pair<IVertex, Vertex_index> p(v, Vertex_index());
      auto pair = m.insert(p);
      if (pair.second) {
        pair.first->second = sp.mesh().add_vertex(point_2(support_plane, v));
      }
      sp.set_ivertex(pair.first->second, v);
      vertices.push_back(pair.first->second);
    }

    Face_index fi = sp.mesh().add_face(vertices);
    if (fi == Support_plane::Mesh::null_face()) {
      std::cout << "ERROR: invalid face created!" << std::endl;
      for (std::size_t i = 0; i < f.vertices.size(); i++) {
        std::cout << "2 " << point_3(f.vertices[i]) << " " << point_3(f.vertices[(i + 1) % f.vertices.size()]) << std::endl;
      }
      std::cout << "ERROR: end of invalid face" << std::endl;
    }

    // Link pedges to iedges
    auto h = sp.mesh().halfedge(fi);
    auto first = h;
    do {
      Edge_index e = sp.mesh().edge(h);
      PVertex t = PVertex(support_plane, sp.mesh().target(h));
      PVertex s = PVertex(support_plane, sp.mesh().source(h));
      IVertex it = ivertex(t);
      IVertex is = ivertex(s);
      sp.set_iedge(e, m_intersection_graph.edge(is, it));
      h = sp.mesh().next(h);
    } while (h != first);

    f.part_of_partition = true;

    return PFace(support_plane, fi);
  }

  void clear_pfaces(const std::size_t support_plane_idx) {
    support_plane(support_plane_idx).clear_pfaces();
  }

  void clear_polygon_faces(const std::size_t support_plane_idx) {
    Mesh& m = mesh(support_plane_idx);
    for (const auto& fi : m.faces()) {
      m.remove_face(fi);
    }
    for (const auto& ei : m.edges()) {
      m.remove_edge(ei);
    }
    for (const auto& vi : m.vertices()) {
      m.set_halfedge(vi, Halfedge_index());
    }
  }

  PVertex opposite(const PEdge& pedge, const PVertex& pvertex) const {
    if (mesh(pedge).target(mesh(pedge).halfedge(pedge.second)) == pvertex.second) {
      return PVertex(pedge.first, mesh(pedge).source(mesh(pedge).halfedge(pedge.second)));
    }
    CGAL_assertion(mesh(pedge).source(mesh(pedge).halfedge(pedge.second)) == pvertex.second);
    return PVertex(pedge.first, mesh(pedge).target(mesh(pedge).halfedge(pedge.second)));
  }

  Point_3 centroid_of_pface(const PFace& pface) const {
    const std::function<Point_3(PVertex)> unary_f =
      [&](const PVertex& pvertex) -> Point_3 {
      return point_3(pvertex);
    };
    const std::vector<Point_3> polygon(
      boost::make_transform_iterator(pvertices_of_pface(pface).begin(), unary_f),
      boost::make_transform_iterator(pvertices_of_pface(pface).end(), unary_f));
    CGAL_assertion(polygon.size() >= 3);
    return CGAL::centroid(polygon.begin(), polygon.end());
  }

  PEdges_of_pface pedges_of_pface(const PFace& pface) const {
    const auto pedges = PEdges_of_pface(
      boost::make_transform_iterator(halfedges_around_face(halfedge(pface.second, mesh(pface)), mesh(pface)).begin(), Halfedge_to_pedge(pface.first, mesh(pface))),
      boost::make_transform_iterator(halfedges_around_face(halfedge(pface.second, mesh(pface)), mesh(pface)).end(), Halfedge_to_pedge(pface.first, mesh(pface))));
    CGAL_assertion(pedges.size() >= 3);
    return pedges;
  }

  PEdges_around_pvertex pedges_around_pvertex(const PVertex& pvertex) const {
    const auto pedges = PEdges_around_pvertex(
      boost::make_transform_iterator(halfedges_around_target(halfedge(pvertex.second, mesh(pvertex)), mesh(pvertex)).begin(), Halfedge_to_pedge(pvertex.first, mesh(pvertex))),
      boost::make_transform_iterator(halfedges_around_target(halfedge(pvertex.second, mesh(pvertex)), mesh(pvertex)).end(), Halfedge_to_pedge(pvertex.first, mesh(pvertex))));
    CGAL_assertion(pedges.size() >= 2);
    return pedges;
  }

  PVertices_of_pface pvertices_of_pface(const PFace& pface) const {
    const auto pvertices = PVertices_of_pface(
      boost::make_transform_iterator(halfedges_around_face(halfedge(pface.second, mesh(pface)), mesh(pface)).begin(), Halfedge_to_pvertex(pface.first, mesh(pface))),
      boost::make_transform_iterator(halfedges_around_face(halfedge(pface.second, mesh(pface)), mesh(pface)).end(), Halfedge_to_pvertex(pface.first, mesh(pface))));
    CGAL_assertion(pvertices.size() >= 3);
    return pvertices;
  }

  std::vector<Volume_cell> incident_volumes(const PFace& query_pface) const {
    std::vector<Volume_cell> nvolumes;
    for (const auto& volume : m_volumes) {
      for (const auto& pface : volume.pfaces) {
        if (pface == query_pface) nvolumes.push_back(volume);
      }
    }
    return nvolumes;
  }

  void incident_faces(const IEdge& query_iedge, std::vector<PFace>& nfaces) const {
    nfaces.clear();
    for (const auto plane_idx : intersected_planes(query_iedge)) {
      for (const auto pedge : pedges(plane_idx)) {
        if (iedge(pedge) == query_iedge) {
          const auto& m = mesh(plane_idx);
          const auto he = m.halfedge(pedge.second);
          const auto op = m.opposite(he);
          const auto face1 = m.face(he);
          const auto face2 = m.face(op);
          if (face1 != Support_plane::Mesh::null_face()) {
            nfaces.push_back(PFace(plane_idx, face1));
          }
          if (face2 != Support_plane::Mesh::null_face()) {
            nfaces.push_back(PFace(plane_idx, face2));
          }
        }
      }
    }
  }

  const std::vector<std::size_t>& input(const PFace& pface) const{ return support_plane(pface).input(pface.second); }
  std::vector<std::size_t>& input(const PFace& pface) { return support_plane(pface).input(pface.second); }

  const int& k(const std::size_t support_plane_idx) const { return support_plane(support_plane_idx).k(); }
  int& k(const std::size_t support_plane_idx) { return support_plane(support_plane_idx).k(); }

  const Vector_2& direction(const PVertex& pvertex) const { return support_plane(pvertex).direction(pvertex.second); }
  Vector_2& direction(const PVertex& pvertex) { return support_plane(pvertex).direction(pvertex.second); }

  const FT speed(const PVertex& pvertex) { return support_plane(pvertex).speed(pvertex.second); }

  /*******************************
  **          ISimplices        **
  ********************************/

  static IVertex null_ivertex() { return Intersection_graph::null_ivertex(); }
  static IEdge null_iedge() { return Intersection_graph::null_iedge(); }

  decltype(auto) ivertices() const { return m_intersection_graph.vertices(); }
  decltype(auto) iedges() const { return m_intersection_graph.edges(); }

  std::size_t nb_intersection_lines() const { return m_intersection_graph.nb_lines(); }
  std::size_t line_idx(const IEdge& iedge) const { return m_intersection_graph.line(iedge); }
  std::size_t line_idx(const PVertex& pvertex) const { return line_idx(iedge(pvertex)); }

  const IVertex add_ivertex(const IkPoint_3& point, const std::set<std::size_t>& support_planes_idx) {
    std::vector<std::size_t> vec_planes;
    std::copy(
      support_planes_idx.begin(),
      support_planes_idx.end(),
      std::back_inserter(vec_planes));
    const auto pair = m_intersection_graph.add_vertex(point, vec_planes);
    const auto ivertex = pair.first;
    return ivertex;
  }

  void add_iedge(const std::set<std::size_t>& support_planes_idx, std::vector<IVertex>& vertices) {
    const auto source = m_intersection_graph.point_3(vertices.front());
    std::sort(vertices.begin(), vertices.end(),
      [&](const IVertex& a, const IVertex& b) -> bool {
        const auto ap = m_intersection_graph.point_3(a);
        const auto bp = m_intersection_graph.point_3(b);
        const auto sq_dist_a = CGAL::squared_distance(source, ap);
        const auto sq_dist_b = CGAL::squared_distance(source, bp);
        return (sq_dist_a < sq_dist_b);
      }
    );

    typename Intersection_kernel::Line_3 line;
    auto it = support_planes_idx.begin();
    bool intersect = intersection(m_support_planes[*it++].exact_plane(), m_support_planes[*it++].exact_plane(), line);
    CGAL_assertion(intersect);

    std::size_t line_idx = m_intersection_graph.add_line(line);
    for (std::size_t i = 0; i < vertices.size() - 1; ++i) {

      CGAL_assertion(!is_zero_length_iedge(vertices[i], vertices[i + 1]));
      const auto pair = m_intersection_graph.add_edge(
        vertices[i], vertices[i + 1], support_planes_idx);
      const auto iedge = pair.first;
      CGAL_assertion(pair.second);// is inserted?
      m_intersection_graph.set_line(iedge, line_idx);

      for (const auto support_plane_idx : support_planes_idx) {
        support_plane(support_plane_idx).unique_iedges().insert(iedge);
      }
    }
  }

  const IVertex source(const IEdge& edge) const { return m_intersection_graph.source(edge); }
  const IVertex target(const IEdge& edge) const { return m_intersection_graph.target(edge); }

  const IVertex opposite(const IEdge& edge, const IVertex& ivertex) const {
    const auto out = source(edge);
    if (out == ivertex) {
      return target(edge);
    }
    CGAL_assertion(target(edge) == ivertex);
    return out;
  }

  decltype(auto) incident_iedges(const IVertex& ivertex) const {
    return m_intersection_graph.incident_edges(ivertex);
  }

  const std::vector<IEdge>& iedges(const std::size_t support_plane_idx) const {
    return support_plane(support_plane_idx).iedges();
  }

  std::vector<IEdge>& iedges(const std::size_t support_plane_idx) {
    return support_plane(support_plane_idx).iedges();
  }

  const std::vector<Bbox_2>& ibboxes(const std::size_t support_plane_idx) const {
    return support_plane(support_plane_idx).ibboxes();
  }

  std::vector<Bbox_2>& ibboxes(const std::size_t support_plane_idx) {
    return support_plane(support_plane_idx).ibboxes();
  }

  const std::set<std::size_t>& intersected_planes(const IEdge& iedge) const {
    return m_intersection_graph.intersected_planes(iedge);
  }

  const std::set<std::size_t> intersected_planes(const IVertex& ivertex, const bool keep_bbox = true) const {
    std::set<std::size_t> out;
    for (const auto &incident_iedge : incident_iedges(ivertex)) {
      for (const auto &support_plane_idx : intersected_planes(incident_iedge)) {
        if (!keep_bbox && support_plane_idx < 6) {
          continue;
        }
        out.insert(support_plane_idx);
      }
    }
    return out;
  }

  bool is_zero_length_iedge(const IVertex& a, const IVertex& b) const {
    const auto& p = m_intersection_graph.point_3(a);
    const auto& q = m_intersection_graph.point_3(b);
    return KSP::internal::distance(p, q) == 0;
  }

  bool is_iedge(const IVertex& source, const IVertex& target) const {
    return m_intersection_graph.is_edge(source, target);
  }

  bool is_bbox_iedge(const IEdge& edge) const {
    for (const auto support_plane_idx : m_intersection_graph.intersected_planes(edge)) {
      if (support_plane_idx < 6) {
        return true;
      }
    }
    return false;
  }

  /*******************************
  **          STRINGS           **
  ********************************/

  inline const std::string str(const IEdge& iedge) const {
    std::ostringstream oss; oss << "IEdge" << iedge; return oss.str();
  }

  inline const std::string lstr(const PEdge& pedge) const {
    return "PEdge(" + std::to_string(pedge.first) + ":e" + std::to_string(pedge.second)
      + ")[v" + std::to_string(source(pedge).second) + "->v" + std::to_string(target(pedge).second) + "]";
  }

  /*******************************
  **        CONNECTIVITY        **
  ********************************/

  bool has_ivertex(const PVertex& pvertex) const { return support_plane(pvertex).has_ivertex(pvertex.second); }
  IVertex ivertex(const PVertex& pvertex) const { return support_plane(pvertex).ivertex(pvertex.second); }

  bool has_iedge(const PVertex& pvertex) const { return support_plane(pvertex).has_iedge(pvertex.second); }
  IEdge iedge(const PVertex& pvertex) const { return support_plane(pvertex).iedge(pvertex.second); }

  bool has_iedge(const PEdge& pedge) const { return support_plane(pedge).has_iedge(pedge.second); }
  IEdge iedge(const PEdge& pedge) const { return support_plane(pedge).iedge(pedge.second); }

  /*******************************
  **        CONVERSIONS         **
  ********************************/

  IkPoint_2 to_2d(const std::size_t support_plane_idx, const IVertex& ivertex) const {
    return support_plane(support_plane_idx).to_2d(point_3(ivertex));
  }

  Segment_2 to_2d(const std::size_t support_plane_idx, const Segment_3& segment_3) const {
    return support_plane(support_plane_idx).to_2d(segment_3);
  }

  IkSegment_2 to_2d(const std::size_t support_plane_idx, const IkSegment_3& segment_3) const {
    return support_plane(support_plane_idx).to_2d(segment_3);
  }

  Point_2 to_2d(const std::size_t support_plane_idx, const Point_3& point_3) const {
    return support_plane(support_plane_idx).to_2d(point_3);
  }

  IkPoint_2 to_2d(const std::size_t support_plane_idx, const IkPoint_3& point_3) const {
    return support_plane(support_plane_idx).to_2d(point_3);
  }

  Point_2 point_2(const PVertex& pvertex) const {
    return support_plane(pvertex).point_2(pvertex.second);
  }

  Point_2 point_2(const std::size_t support_plane_idx, const IVertex& ivertex) const {
    return support_plane(support_plane_idx).to_2d(from_exact(point_3(ivertex)));
  }

  IkSegment_2 segment_2(const std::size_t support_plane_idx, const IEdge& iedge) const {
    return support_plane(support_plane_idx).to_2d(segment_3(iedge));
  }

  Point_3 to_3d(const std::size_t support_plane_idx, const Point_2& point_2) const {
    return support_plane(support_plane_idx).to_3d(point_2);
  }

  IkPoint_3 to_3d(const std::size_t support_plane_idx, const IkPoint_2& point_2) const {
    return support_plane(support_plane_idx).to_3d(point_2);
  }

  Point_3 point_3(const PVertex& pvertex) const {
    return support_plane(pvertex).point_3(pvertex.second);
  }

  IkPoint_3 point_3(const IVertex& vertex) const {
    return m_intersection_graph.point_3(vertex);
  }

  Segment_3 segment_3(const PEdge& pedge) const {
    return support_plane(pedge).segment_3(pedge.second);
  }

  IkSegment_3 segment_3(const IEdge& edge) const {
    return m_intersection_graph.segment_3(edge);
  }

  /*******************************
  **    CHECKING PROPERTIES     **
  ********************************/

  template<typename Pair>
  bool is_valid_polygon(const std::size_t sp_idx, const std::vector<Pair>& points) const {
    std::vector< std::pair<IkPoint_2, bool> > polygon;
    polygon.reserve(points.size());
    for (const auto& pair : points) {
      const auto& p = pair.first;
      const auto q = to_2d(sp_idx, p);
      polygon.push_back(std::make_pair(q, true));
    }
    CGAL_assertion(polygon.size() == points.size());

    // const bool is_simple = support_plane(sp_idx).is_simple_polygon(polygon);
    // const bool is_convex = support_plane(sp_idx).is_convex_polygon(polygon);
    const bool is_valid = support_plane(sp_idx).is_valid_polygon(polygon);

    if (!is_valid) {
      for (const auto& pair : polygon) {
        std::cout << to_3d(sp_idx, pair.first) << std::endl;
      }
    }

    // CGAL_assertion_msg(is_simple, "ERROR: POLYGON IS NOT SIMPLE!");
    // CGAL_assertion_msg(is_convex, "ERROR: POLYGON IS NOT CONVEX!");
    CGAL_assertion_msg(is_valid, "ERROR: POLYGON IS NOT VALID!");
    return is_valid;
  }

  bool check_bbox() const {
    for (std::size_t i = 0; i < 6; ++i) {
      const auto pfaces = this->pfaces(i);
      for (const auto pface : pfaces) {
        for (const auto pvertex : pvertices_of_pface(pface)) {
          if (!has_ivertex(pvertex)) {
            std::cout << "debug pvertex: " << std::endl;
            CGAL_assertion_msg(has_ivertex(pvertex), "ERROR: BBOX VERTEX IS MISSING AN IVERTEX!");
            return false;
          }
        }
        for (const auto pedge : pedges_of_pface(pface)) {
          if (!has_iedge(pedge)) {
            std::cout << "debug pedge: " << std::endl;
            CGAL_assertion_msg(has_iedge(pedge), "ERROR: BBOX EDGE IS MISSING AN IEDGE!");
            return false;
          }
        }
      }
    }
    return true;
  }

  bool check_interior() const {
    for (std::size_t i = 6; i < number_of_support_planes(); ++i) {
      const auto pfaces = this->pfaces(i);
      for (const auto pface : pfaces) {
        for (const auto pvertex : pvertices_of_pface(pface)) {
          if (!has_ivertex(pvertex)) {
            std::cout << "debug pvertex " << std::endl;
            CGAL_assertion_msg(has_ivertex(pvertex), "ERROR: INTERIOR VERTEX IS MISSING AN IVERTEX!");
            return false;
          }
        }
        for (const auto pedge : pedges_of_pface(pface)) {
          if (!has_iedge(pedge)) {
            std::cout << "debug pedge " << std::endl;
            CGAL_assertion_msg(has_iedge(pedge), "ERROR: INTERIOR EDGE IS MISSING AN IEDGE!");
            return false;
          }
        }
      }
    }
    return true;
  }

  bool check_vertices() const {
    bool success = true;

    for (const auto vertex : m_intersection_graph.vertices()) {
      const auto nedges = m_intersection_graph.incident_edges(vertex);
      if (nedges.size() <= 2) {
        std::cout << "ERROR: CURRENT NUMBER OF EDGES = " << nedges.size() << std::endl;
        CGAL_assertion_msg(nedges.size() > 2,
        "ERROR: VERTEX MUST HAVE AT LEAST 3 NEIGHBORS!");
        success = false;
      }
    }
    return success;
  }

  bool check_edges() const {
    bool success = true;

    std::vector<PFace> nfaces;
    for (const auto edge : m_intersection_graph.edges()) {
      incident_faces(edge, nfaces);
      if (nfaces.size() == 1) {
        std::cout << "ERROR: CURRENT NUMBER OF FACES = " << nfaces.size() << std::endl;
        std::cout << edge << std::endl;
        std::cout << "PFace(" << nfaces[0].first << " " << nfaces[0].second << ")" << std::endl;
        CGAL_assertion_msg(nfaces.size() != 1,
        "ERROR: EDGE MUST HAVE 0 OR AT LEAST 2 NEIGHBORS!");
        success = false;
      }
    }
    return success;
  }

  bool check_faces() const {
    bool success = true;

    for (std::size_t i = 0; i < number_of_support_planes(); ++i) {
      const auto pfaces = this->pfaces(i);
      for (const auto pface : pfaces) {
        const auto nvolumes = incident_volumes(pface);
        if (nvolumes.size() == 0 || nvolumes.size() > 2) {
          std::cout << "ERROR: CURRENT NUMBER OF VOLUMES = " << nvolumes.size() << std::endl;
          CGAL_assertion_msg(nvolumes.size() == 1 || nvolumes.size() == 2,
            "ERROR: FACE MUST HAVE 1 OR 2 NEIGHBORS!");
          success = false;
        }
      }
    }

    return success;
  }

  bool is_volume_degenerate(const std::vector<PFace>& pfaces) const {
    for (const auto& pface : pfaces) {
      const auto pedges = pedges_of_pface(pface);
      const std::size_t n = pedges.size();

      std::size_t count = 0;
      for (const auto pedge : pedges) {
        CGAL_assertion(has_iedge(pedge));
        const auto iedge = this->iedge(pedge);
        const std::size_t num_found = find_adjacent_pfaces(pface, iedge, pfaces);
        if (num_found == 1) ++count;
      }
      if (count != n) {
        std::cout << "- current number of neighbors " << count << " != " << n << std::endl;
        dump_info(*this, pface, *pedges.begin(), pfaces, "");
        return true;
      }
    }
    return false;
  }
};

#endif //DOXYGEN_RUNNING

} // namespace internal
} // namespace KSP_3
} // namespace CGAL

#endif // CGAL_KSP_3_DATA_STRUCTURE_H
