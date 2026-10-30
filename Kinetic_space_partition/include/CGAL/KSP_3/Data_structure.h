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
  using Vertex = typename Support_plane::Vertex;

  using FT = typename Kernel::FT;
  using IkFT = typename Intersection_kernel::FT;
  using Point_2 = typename Kernel::Point_2;
  using IkPoint_2 = typename Intersection_kernel::Point_2;
  using Point_3 = typename Kernel::Point_3;
  using IkPoint_3 = typename Intersection_kernel::Point_3;
  using Segment_2 = typename Kernel::Segment_2;
  using IkSegment_2 = typename Intersection_kernel::Segment_2;
  using Segment_3 = typename Kernel::Segment_3;
  using IkSegment_3 = typename Intersection_kernel::Segment_3;
  using Vector_2 = typename Kernel::Vector_2;
  using IkVector_2 = typename Intersection_kernel::Vector_2;
  using Direction_2 = typename Kernel::Direction_2;
  using IkDirection_2 = typename Intersection_kernel::Direction_2;
  using Triangle_2  = typename Kernel::Triangle_2;
  using Line_2      = typename Kernel::Line_2;
  using IkLine_2    = typename Intersection_kernel::Line_2;
  using Plane_3     = typename Kernel::Plane_3;

  using Polygon_2  = CGAL::Polygon_2<Kernel>;

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
    Make_PSimplex(const std::size_t sp_idx)
      : support_plane_idx(sp_idx)
    {}

    const result_type operator()(const argument_type& arg) const {
      return result_type(support_plane_idx, arg);
    }
  };

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

    Halfedge_to_pvertex(const std::size_t sp_idx, const Mesh& m)
      : support_plane_idx(sp_idx), mesh(m)
    {}

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

    Halfedge_to_pedge(const std::size_t sp_idx, const Mesh& m)
      : support_plane_idx(sp_idx), mesh(m)
    {}

    const result_type operator()(const argument_type& arg) const {
      return result_type(support_plane_idx, mesh.edge(arg));
    }
  };


  using PEdge_of_pface_iterator = boost::transform_iterator<Halfedge_to_pedge, CGAL::Halfedge_around_face_iterator<Mesh> >;
  using PEdges_of_pface = CGAL::Iterator_range<PEdge_of_pface_iterator>;


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

    void add_pface(const PFace& pface, const int neighbor) {
      pfaces.push_back(pface);
      neighbors.push_back(neighbor);
    }
  };

private:
  std::vector<Support_plane> m_support_planes;
  Intersection_graph m_intersection_graph;

  std::vector<Vertex> m_kinetic_vertices;

  To_exact to_exact;
  From_exact from_exact;


  std::string m_prefix;
  std::string m_log_prefix; // Prefix for each line written to std::cout, e.g., the index of the partition.

  std::vector<Volume_cell> m_volumes;
  std::vector<Point_3> m_vertices;
  std::vector<typename Intersection_kernel::Point_3> m_exact_vertices;
  std::vector<int> m_ivertex2vertex; // Used to map ivertices to m_vertices which only contain vertices of the finalized kinetic partition.
  std::vector<std::pair<int, int> > m_face2volumes;

  std::map<PFace, std::size_t> m_face2index;
  std::vector<std::vector<std::size_t> > m_face2vertices;
  std::vector<bool> m_part_of_initial;
  std::vector<std::size_t> m_face_input; // Index of the input polygon of each face, std::size_t(-1) if none.
  std::map<PFace, std::pair<int, int> > m_pface_neighbors;
  std::vector<std::size_t> m_face2sp;
  std::vector<std::set<std::size_t> > m_sp2input_polygon;
  std::map<std::size_t, std::size_t> m_input_polygon_map; // Maps index of input polygon onto support plane indices.

public:
  Data_structure(const std::string &prefix, const std::string &log_prefix = std::string())
    : to_exact(), from_exact(), m_prefix(prefix), m_log_prefix(log_prefix)
  {}

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
  }

  // Everything that is modified by the propagation and the finalization:
  // the support planes, the kinetic vertices and parts of the intersection graph.
  struct Initial_state {
    std::vector<Support_plane> support_planes;
    std::vector<Vertex> kinetic_vertices;
    typename Intersection_graph::Propagation_state igraph;
  };

  // Saves a copy of the initialized state. Together with restore_initial_state(), it allows
  // to propagate the same initialized data structure several times, e.g., with different k.
  Initial_state save_initial_state() const {
    Initial_state state;
    state.support_planes.reserve(m_support_planes.size());
    for (const Support_plane& sp : m_support_planes)
      state.support_planes.push_back(sp.deep_copy());

    state.kinetic_vertices = m_kinetic_vertices;
    state.igraph = m_intersection_graph.save_propagation_state();

    return state;
  }

  // Restores the state saved by save_initial_state() and discards the results of the finalization.
  void restore_initial_state(Initial_state&& state) {
    CGAL_assertion(state.support_planes.size() == m_support_planes.size());
    m_support_planes = std::move(state.support_planes);
    // The support planes refer to m_kinetic_vertices, so the vector itself has to stay in place.
    m_kinetic_vertices = std::move(state.kinetic_vertices);
    m_intersection_graph.restore_propagation_state(std::move(state.igraph));

    std::vector<Volume_cell>().swap(m_volumes);
    std::vector<Point_3>().swap(m_vertices);
    std::vector<typename Intersection_kernel::Point_3>().swap(m_exact_vertices);
    std::vector<int>().swap(m_ivertex2vertex);
    std::vector<std::pair<int, int> >().swap(m_face2volumes);
    std::map<PFace, std::size_t>().swap(m_face2index);
    std::vector<std::vector<std::size_t> >().swap(m_face2vertices);
    std::vector<bool>().swap(m_part_of_initial);
    std::vector<std::size_t>().swap(m_face_input);
    std::map<PFace, std::pair<int, int> >().swap(m_pface_neighbors);
    std::vector<std::size_t>().swap(m_face2sp);
  }

  /*******************************
  **           ACCESS           **
  ********************************/

  const std::map<std::size_t, std::size_t>& input_polygon_map() const { return m_input_polygon_map; }

  std::vector<Vertex>& kinetic_vertices() {
    return m_kinetic_vertices;
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

  std::vector<std::size_t>& face_to_input_polygon() { return m_face_input; }
  const std::vector<std::size_t>& face_to_input_polygon() const { return m_face_input; }

  std::vector<std::size_t>& face_to_support_plane() { return m_face2sp; }


  const Intersection_graph& igraph() const { return m_intersection_graph; }
  Intersection_graph& igraph() { return m_intersection_graph; }

  const std::string& prefix() const { return m_prefix; }
  const std::string& log_prefix() const { return m_log_prefix; }

  void resize(const std::size_t number_of_items) {
    m_support_planes.resize(number_of_items);
  }

  std::vector<Volume_cell>& volumes() { return m_volumes; }
  const std::vector<Volume_cell>& volumes() const { return m_volumes; }

  const std::vector<std::size_t>& volume(std::size_t volume_index) {
    return m_volumes[volume_index].faces;
  }

  const std::vector<std::size_t> &face(std::size_t face_index) const { return m_face2vertices[face_index]; }
  const Point_3& vertex(std::size_t vertex_index) const { return m_vertices[vertex_index]; }

  /*******************************
  **      SUPPORT PLANES        **
  ********************************/

  const Support_plane& support_plane(const std::size_t idx) const { return m_support_planes[idx]; }

  Support_plane& support_plane(const std::size_t idx) { return m_support_planes[idx]; }

  const Mesh& mesh(const std::size_t support_plane_idx) const { return support_plane(support_plane_idx).mesh(); }

  Mesh& mesh(const std::size_t support_plane_idx) { return support_plane(support_plane_idx).mesh(); }

  std::size_t number_of_support_planes() const {
    return m_support_planes.size();
  }

  bool is_bbox_support_plane(const std::size_t support_plane_idx) const {
    return (support_plane_idx < 6);
  }

  template<typename PointRange>
  std::pair<std::size_t, bool> add_support_plane(const PointRange& polygon, const bool is_bbox, const typename Intersection_kernel::Plane_3& plane) {
    const Support_plane new_support_plane(polygon, is_bbox, plane, m_kinetic_vertices);
    const std::size_t uninitialized = static_cast<std::size_t>(-1);
    std::size_t support_plane_idx = uninitialized;

    for (std::size_t i = 0; i < number_of_support_planes(); ++i) {
      if (new_support_plane == support_plane(i)) {
        support_plane_idx = i;
        return std::make_pair(support_plane_idx, false);
      }
    }

    if (support_plane_idx == uninitialized) {
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
    const Support_plane new_support_plane(polygon, is_bbox, m_kinetic_vertices);
    const std::size_t uninitialized = static_cast<std::size_t>(- 1);
    std::size_t support_plane_idx = uninitialized;

    for (std::size_t i = 0; i < number_of_support_planes(); ++i) {
      if (new_support_plane == support_plane(i)) {
        support_plane_idx = i;
        return std::make_pair(support_plane_idx, false);
      }
    }

    if (support_plane_idx == uninitialized) {
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
    Support_plane &sp = support_plane(sp_idx);

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

    is_valid_polygon(sp_idx, polygon);

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

    // If the plane contains an existing iedge, e.g., an edge of the bbox, consecutive vertices of the polygon are
    // connected by this iedge. It is reused instead of creating a new iedge on the same line.
    std::vector<IEdge> existing(n, null_iedge());
    for (std::size_t i = 0; i < n; ++i) {
      const std::size_t ip = (i + 1) % n;
      if (polygon[i].second.first != null_ivertex() && polygon[ip].second.first != null_ivertex()
        && m_intersection_graph.is_edge(polygon[i].second.first, polygon[ip].second.first))
        existing[i] = m_intersection_graph.edge(polygon[i].second.first, polygon[ip].second.first);
    }

    for (std::size_t i = 0; i < n; ++i) {
      const std::size_t ip = (i + 1) % n;
      const auto& item = polygon[i].second;

      if (existing[i] != null_iedge()) {
        common_bbox_planes_idx.push_back(static_cast<std::size_t>(-1));
        vertices.push_back(item.first);
        continue;
      }

      const auto& iplanes0 = all_iplanes[i];
      const auto& iplanes1 = all_iplanes[ip];

      const std::size_t uninitialized = static_cast<std::size_t>(-1);
      std::size_t common_bbox_plane_idx = uninitialized;
      bool dump = false;
      const std::function<void(const std::size_t& idx)> lambda =
        [&](const std::size_t& idx) {
        if (idx < 6) {

          if (common_bbox_plane_idx != uninitialized)
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

      CGAL_assertion(common_bbox_plane_idx != uninitialized);
      common_bbox_planes_idx.push_back(common_bbox_plane_idx);

      const auto pair = map_lines_idx.insert(
        std::make_pair(common_bbox_plane_idx, uninitialized));
      const bool is_inserted = pair.second;
      if (is_inserted) {
        typename Intersection_kernel::Line_3 line;
        CGAL_assertion_code(bool intersect =)  intersection(plane, m_support_planes[common_bbox_plane_idx].exact_plane(), line);
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

      if (existing[i] != null_iedge()) {
        const IEdge iedge = existing[i];
        m_intersection_graph.intersected_planes(iedge).insert(sp_idx);
        sp.iedges().insert(iedge);
        const std::size_t line_idx = m_intersection_graph.line(iedge);
        const IkLine_2 line(sp.to_2d(point_3(vertices[i])), sp.to_2d(point_3(vertices[ip])));
        sp.data().bbox_lines[line_idx] = line;
        sp.data().lines[line_idx] = line;
        continue;
      }

      const std::size_t common_bbox_plane_idx = common_bbox_planes_idx[i];
      const auto new_iedge = m_intersection_graph.add_edge(vertices[i], vertices[ip], sp_idx).first;

      m_intersection_graph.intersected_planes(new_iedge).insert(common_bbox_plane_idx);
      CGAL_assertion(map_lines_idx.find(common_bbox_plane_idx) != map_lines_idx.end());
      std::size_t line_idx = map_lines_idx.at(common_bbox_plane_idx);
      m_intersection_graph.set_line(new_iedge, line_idx);
      sp.iedges().insert(new_iedge);
      support_plane(common_bbox_plane_idx).iedges().insert(new_iedge);
      auto it = sp.data().bbox_lines.find(line_idx);
      if (it == sp.data().bbox_lines.end()) {
        const typename Intersection_kernel::Line_3 &l3 = m_intersection_graph.line(line_idx);
        const typename Intersection_kernel::Point_2 p1 = sp.exact_plane().to_2d(l3.point());
        const typename Intersection_kernel::Point_2 p2 = sp.exact_plane().to_2d(l3.point() + l3.to_vector());
        IkLine_2 line = IkLine_2(p1, p2);
        sp.data().bbox_lines[line_idx] = line;
        sp.data().lines[line_idx] = line;
      }

      // No further treatment necessary for exact intersections at vertices of edges.

      if (item.first == null_ivertex()) { // edge case, split
        const auto& iplanes = all_iplanes[i];
        CGAL_assertion(iplanes.size() >= 2);
        CGAL_assertion(item.second.size() == 1);
        const auto& iedge = item.second[0];
        for (const std::size_t plane_idx : iplanes) {
          support_plane(plane_idx).iedges().erase(iedge);
        }
        const auto edges = m_intersection_graph.split_edge(iedge, vertices[i]);

        const auto& iplanes0 = m_intersection_graph.intersected_planes(edges.first);
        for (const std::size_t plane_idx : iplanes0) {
          support_plane(plane_idx).iedges().insert(edges.first);
        }

        const auto& iplanes1 = m_intersection_graph.intersected_planes(edges.second);
        for (const std::size_t plane_idx : iplanes1) {
          support_plane(plane_idx).iedges().insert(edges.second);
        }
      }
    }
  }

  template<typename PointRange>
  void add_bbox_polygon(const PointRange& polygon) {
    bool is_added = true;
    const std::size_t uninitialized = static_cast<std::size_t>(-1);
    std::size_t support_plane_idx = uninitialized;
    std::tie(support_plane_idx, is_added) = add_support_plane(polygon, true);
    CGAL_assertion(is_added);
    CGAL_assertion(support_plane_idx != uninitialized);

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
      support_plane(support_plane_idx).iedges().insert(iedge);
    }
  }

  void add_input_polygon( const std::size_t support_plane_idx, const std::vector<std::size_t>& input_indices, const std::vector<typename Intersection_kernel::Point_2>& polygon) {
    std::vector< std::pair<typename Intersection_kernel::Point_2, bool> > points;
    points.reserve(polygon.size());
    for (const auto& point : polygon) {
      points.push_back(std::make_pair(point, true));
    }
    CGAL_assertion(points.size() == polygon.size());

    // The polygon is the exact convex hull of the clipped input polygon, i.e., it is strictly convex and counterclockwise.
    // Thus, no inexact preprocessing is needed. The vertices are sorted exactly around the centroid by the support plane.
    CGAL_assertion(points.size() >= 3);
    support_plane(support_plane_idx).add_input_polygon(points, input_indices.front());
    for (const std::size_t input_index : input_indices) {
      m_input_polygon_map[input_index] = support_plane_idx;
      m_sp2input_polygon[support_plane_idx].insert(input_index);
    }
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

  // Provides the iedges incident to ivertex within the support plane sp_idx in counterclockwise order.
  // Only the cyclic order is defined, not the first iedge.
  // The directions are first sorted by their approximated angles. Directions whose angles cannot be separated by the
  // approximation are then sorted with exact predicates. Thus, the order is correct for arbitrarily short iedges and
  // small angles, while exact computations are only required for nearly parallel iedges.
  void get_and_sort_all_connected_iedges(const std::size_t sp_idx, const IVertex& ivertex, std::vector< std::pair<IEdge, IkDirection_2> >& iedges) const {
    using IkVector_2 = typename Intersection_kernel::Vector_2;
    const Support_plane& sp = support_plane(sp_idx);
    const typename Intersection_kernel::Point_2 center = sp.to_2d(point_3(ivertex));

    struct Entry {
      IEdge edge;
      IkVector_2 v;
      double angle, error;
    };

    std::vector<Entry> entries;
    bool approximation_valid = true;
    for (const IEdge& inc_iedge : incident_iedges(ivertex)) {
      const auto& iplanes = intersected_planes(inc_iedge);
      if (iplanes.find(sp_idx) == iplanes.end())
        continue;

      const IkVector_2 v = sp.to_2d(point_3(opposite(inc_iedge, ivertex))) - center;
      const std::pair<double, double> x = CGAL::to_interval(v.x()), y = CGAL::to_interval(v.y());
      const double mx = (x.first + x.second) / 2, my = (y.first + y.second) / 2;
      const double len = std::sqrt(mx * mx + my * my);
      double angle = std::atan2(my, mx);
      if (angle < 0)
        angle += 2 * CGAL_PI;

      // Bound of the angular error of the approximation, including the rounding of atan2.
      const double width = (std::max)(x.second - x.first, y.second - y.first);
      double error = 4 * width / len + 1e-12;
      if (!(len > 4 * width) || !(error < 0.01)) {
        approximation_valid = false;
        error = 0;
      }
      entries.push_back({ inc_iedge, v, angle, error });
    }
    CGAL_assertion(entries.size() > 0);

    // Exact comparison of two directions within less than a half-turn: a precedes b if b lies counterclockwise of a.
    auto ccw = [](const Entry& a, const Entry& b) {
      return CGAL::orientation(a.v, b.v) == CGAL::LEFT_TURN;
    };

    if (approximation_valid && entries.size() > 1) {
      std::sort(entries.begin(), entries.end(), [](const Entry& a, const Entry& b) { return a.angle < b.angle; });

      // Start after the largest gap, so that no group of nearly parallel directions is split at the start.
      std::size_t start = 0;
      double max_gap = -1;
      for (std::size_t i = 0; i < entries.size(); i++) {
        const std::size_t n = (i + 1) % entries.size();
        double gap = entries[n].angle - entries[i].angle;
        if (n == 0)
          gap += 2 * CGAL_PI;
        if (gap > max_gap) {
          max_gap = gap;
          start = n;
        }
      }
      std::rotate(entries.begin(), entries.begin() + start, entries.end());

      // Groups of consecutive directions that cannot be separated by the approximation are sorted exactly.
      // The errors are small, so each group spans much less than a half-turn and the exact comparison is a strict weak ordering.
      const double cut = entries[0].angle;
      auto rel = [&](const Entry& e) { double a = e.angle - cut; return (a < 0) ? a + 2 * CGAL_PI : a; };
      if (max_gap <= 2 * 0.01)
        approximation_valid = false;
      else {
        std::size_t first = 0;
        for (std::size_t i = 1; i <= entries.size(); i++) {
          if (i == entries.size() || rel(entries[i]) - rel(entries[i - 1]) > entries[i].error + entries[i - 1].error) {
            if (rel(entries[i - 1]) - rel(entries[first]) > CGAL_PI / 2) {
              approximation_valid = false;
              break;
            }
            if (i - first > 1)
              std::sort(entries.begin() + first, entries.begin() + i, ccw);
            first = i;
          }
        }
      }
    }

    if (!approximation_valid) {
      // Directions with an angle in [0, pi) come first, those in [pi, 2 pi) second, sorted exactly within each half-plane.
      auto half = [](const IkVector_2& v) {
        return (CGAL::is_positive(v.y()) || (CGAL::is_zero(v.y()) && CGAL::is_positive(v.x()))) ? 0 : 1;
      };
      std::sort(entries.begin(), entries.end(), [&](const Entry& a, const Entry& b) {
        const int ha = half(a.v), hb = half(b.v);
        if (ha != hb)
          return ha < hb;
        return ccw(a, b);
      });
    }

    iedges.clear();
    for (const Entry& e : entries)
      iedges.push_back(std::make_pair(e.edge, IkDirection_2(e.v)));
  }

  void init_border(IFace f_idx, bool initial = false) {
    typename Intersection_graph::Face_property& f = m_intersection_graph.face(f_idx);
    Support_plane& sp = m_support_planes[f.support_plane];

    if (!f.part_of_partition) {
      const PFace &pface = add_iface_to_mesh(f.support_plane, f_idx);
      if (initial)
        sp.set_initial(pface.second);
    }

    std::list<IEdge>& border = sp.data().borders[f_idx];
    CGAL_assertion(border.empty());

    std::copy(f.edges.begin(), f.edges.end(), std::back_inserter(border));
  }

  IVertex get_target_ivertex(Vertex& v, IEdge edge) {
    Support_plane& sp = m_support_planes[v.sp_idx];
    typename Intersection_kernel::Vector_3 v3d = sp.exact_plane().base1() * v.v[0] + sp.exact_plane().base2() * v.v[1];
    typename Intersection_kernel::Vector_3 ed = m_intersection_graph.point_3(m_intersection_graph.target(edge)) - m_intersection_graph.point_3(m_intersection_graph.source(edge));

    CGAL_assertion_code({
    if (ed.x() != 0 && ed.y() != 0) {
      IkFT x = v3d.x() / ed.x();
      IkFT y = v3d.y() / ed.y();
      CGAL_assertion(x == y);
    }
    if (ed.x() != 0 && ed.z() != 0) {
      IkFT x = v3d.x() / ed.x();
      IkFT z = v3d.z() / ed.z();
      CGAL_assertion(x == z);
    }
    if (ed.y() != 0 && ed.z() != 0) {
      IkFT y = v3d.y() / ed.y();
      IkFT z = v3d.z() / ed.z();
      CGAL_assertion(y == z);
    }
      });

    if (CGAL::is_positive(v3d * ed))
      return m_intersection_graph.target(edge); // bounds checking is not needed as can_cross() should check for that
    else
      return m_intersection_graph.source(edge); // bounds checking is not needed as can_cross() should check for that
  }

  IVertex get_target_ivertex(Vertex& v, IVertex on, std::size_t line_idx) {
    // Find ivertex entry in line
    const std::vector<std::pair<IkFT, IVertex> >& vol = m_intersection_graph.vertices_on_line(line_idx);
    auto it = std::find_if(vol.begin(), vol.end(), [&](const std::pair<IkFT, IVertex>& v) {
      return v.second == on;
      });
    CGAL_assertion(it != vol.end());

    // Check direction for choosing edge
    Support_plane& sp = m_support_planes[v.sp_idx];
    typename Intersection_kernel::Vector_3 v3d = sp.exact_plane().base1() * v.v[0] + sp.exact_plane().base2() * v.v[1];
    if (CGAL::is_positive(v3d * m_intersection_graph.line(line_idx).to_vector()))
      return std::next(it)->second; // bounds checking is not needed as can_cross() should check for that
    else
      return std::prev(it)->second; // bounds checking is not needed as can_cross() should check for that
  }

  IkFT get_u_on_line(IVertex on, std::size_t line_idx) {
    const std::vector<std::pair<IkFT, IVertex> >& vol = m_intersection_graph.vertices_on_line(line_idx);
    auto it = std::find_if(vol.begin(), vol.end(), [&](const std::pair<IkFT, IVertex>& v) {
      return v.second == on;
      });
    CGAL_assertion(it != vol.end());

    return it->first;
  }

  void get_iedge(Vertex& v, std::size_t line_idx, IVertex &s, IVertex &t) {
    typename Intersection_kernel::Line_3 line = m_intersection_graph.line(line_idx);
    Support_plane& sp = m_support_planes[v.sp_idx];
    IkFT u = line.to_vector() * (sp.to_3d(v.p0) - line.point());

    const std::vector<std::pair<IkFT, IVertex> >& vol = m_intersection_graph.vertices_on_line(line_idx);

    if (u < vol[0].first)
      u = vol[0].first;
    else if (vol.back().first < u)
      u = vol.back().first;

    int low = 0, high = static_cast<int>(vol.size() - 1);
    while (low < high - 1) {
      int mid = (low + high) >> 1;
      if (u == vol[mid].first) {
        low = high = mid;
        break;
      }
      if (u < vol[mid].first)
        high = mid;
      else
        low = mid;
      CGAL_assertion(low != high);
    }

    typename Intersection_kernel::Vector_3 v3d = sp.exact_plane().base1() * v.v[0] + sp.exact_plane().base2() * v.v[1];
    if (CGAL::is_positive(v3d * line.to_vector())) {
      if (low == high) {
        CGAL_assertion(low + 1 < signed(vol.size()));
        s = vol[low].second;
        t = vol[low + 1].second;
      }
      else {
        s = vol[low].second;
        t = vol[high].second;
      }
    }
    else {
      if (low == high) {
        CGAL_assertion(0 < low);
        s = vol[low].second;
        t = vol[low - 1].second;
      }
      else {
        s = vol[high].second;
        t = vol[low].second;
      }
    }
  }

  // Does the polygon have zero area at t = 0? (A piece beyond a line touched by the polygon at t = 0.)
  bool is_zero_area(const std::list<std::size_t>& poly) const {
    auto it = poly.begin();
    const IkPoint_2& p = m_kinetic_vertices[*it].p0;
    while (it != poly.end() && m_kinetic_vertices[*it].p0 == p)
      ++it;
    if (it == poly.end())
      return true;
    const IkPoint_2& q = m_kinetic_vertices[*it].p0;
    for (; it != poly.end(); ++it)
      if (!CGAL::collinear(p, q, m_kinetic_vertices[*it].p0))
        return false;
    return true;
  }

  bool split_polygon(std::list<std::size_t>& poly, std::size_t line_idx, std::list<std::size_t>& pos, std::pair<bool, bool>* zero_area = nullptr) {
    return split_polygon(m_kinetic_vertices, poly, line_idx, pos, zero_area);
  }

  // Side of a point moving with p + t * v with respect to a line at t = 0+.
  static CGAL::Oriented_side side_0plus(const IkLine_2& l, const IkPoint_2& p, const IkVector_2& v) {
    CGAL::Oriented_side s = l.oriented_side(p);
    if (s != CGAL::ON_ORIENTED_BOUNDARY)
      return s;
    return CGAL::Oriented_side(CGAL::sign(l.a() * v.x() + l.b() * v.y()));
  }

  // Splits the convex polygon poly along the line line_idx. The side of each vertex is decided at t = 0+.
  // Vertices on the positive side are moved into pos, the others remain in poly.
  // A vertex staying on the line (moving along it or not moving) belongs to both parts and is duplicated.
  // For each edge crossing the line, a pair of vertices (one for each part) is created at the intersection.
  // If the crossing edge lies on another line, the new vertices are fixed at the intersection of both lines.
  // Returns true if the line splits the polygon, otherwise poly is not changed and pos is empty.
  // If zero_area is given, it reports whether the remaining (negative) part and the positive part have no vertex
  // strictly on their side at t = 0, i.e., whether they are zero-area pieces beyond a line touched by the polygon.
  bool split_polygon(std::vector<Vertex> &vts, std::list<std::size_t> &poly, std::size_t line_idx, std::list<std::size_t> &pos, std::pair<bool, bool>* zero_area = nullptr) {
    Support_plane &sp = m_support_planes[vts[poly.front()].sp_idx];
    const IkLine_2 &l = sp.data().lines[line_idx];

    std::vector<std::size_t> idx(poly.begin(), poly.end());
    const std::size_t n = idx.size();
    std::vector<CGAL::Oriented_side> sign(n);
    bool has_pos = false, has_neg = false, strict_pos = false, strict_neg = false;
    for (std::size_t i = 0; i < n; i++) {
      const CGAL::Oriented_side s0 = l.oriented_side(vts[idx[i]].p0);
      sign[i] = (s0 != CGAL::ON_ORIENTED_BOUNDARY) ? s0 : CGAL::Oriented_side(CGAL::sign(l.a() * vts[idx[i]].v.x() + l.b() * vts[idx[i]].v.y()));
      has_pos |= (sign[i] == CGAL::ON_POSITIVE_SIDE);
      has_neg |= (sign[i] == CGAL::ON_NEGATIVE_SIDE);
      strict_pos |= (s0 == CGAL::ON_POSITIVE_SIDE);
      strict_neg |= (s0 == CGAL::ON_NEGATIVE_SIDE);
    }
    if (zero_area != nullptr)
      *zero_area = std::make_pair(!strict_neg, !strict_pos);

    pos.clear();
    if (!has_pos || !has_neg) {
      // Not split. A vertex moving along the line stays on it and is constrained to it.
      for (std::size_t i = 0; i < n; i++)
        if (sign[i] == CGAL::ON_ORIENTED_BOUNDARY && vts[idx[i]].moving)
          vts[idx[i]].constraints.insert(line_idx);
      return false;
    }

    auto make_twins = [&](std::size_t first, std::size_t second) {
      vts[first].other = second;
      vts[second].other = first;
    };

    // Create the vertices at the crossings. The pair for the crossing from the positive to the negative side is
    // created first (keeps the vertex order of the previous implementation).
    std::vector<std::pair<std::size_t, std::size_t> > crossing(n, std::make_pair(std::size_t(-1), std::size_t(-1)));
    for (int pass = 0; pass < 2; pass++)
      for (std::size_t i = 0; i < n; i++) {
        const std::size_t j = (i + 1) % n;
        if (sign[i] != (pass == 0 ? CGAL::ON_POSITIVE_SIDE : CGAL::ON_NEGATIVE_SIDE) || sign[j] != CGAL::opposite(sign[i]))
          continue;

        const std::size_t a = idx[i], b = idx[j];
        std::set<std::size_t> common;
        std::set_intersection(vts[a].constraints.begin(), vts[a].constraints.end(),
          vts[b].constraints.begin(), vts[b].constraints.end(), std::inserter(common, common.begin()));

        IkPoint_2 intersection = vts[a].p0;
        if (vts[a].p0 != vts[b].p0) {
          auto res = CGAL::intersection(l, IkLine_2(vts[a].p0, vts[b].p0));
          CGAL_assertion_code(bool success =) CGAL::assign(intersection, res);
          CGAL_assertion(success);
        }

        Vertex x(vts[a].sp_idx, intersection, IkVector_2(0, 0));
        x.constraints = common;
        x.constraints.insert(line_idx);
        if (common.empty()) {
          CGAL_assertion(vts[a].moving && vts[b].moving);
          x.v = Support_plane::calculate_edge_speed(vts[a], vts[b], intersection, l);
        }
        else
          x.moving = false;

        vts.push_back(x);
        vts.push_back(x);
        crossing[i] = std::make_pair(vts.size() - 2, vts.size() - 1);
        make_twins(vts.size() - 2, vts.size() - 1);
      }

    std::list<std::size_t> neg;
    for (std::size_t i = 0; i < n; i++) {
      const std::size_t a = idx[i];
      if (sign[i] == CGAL::ON_POSITIVE_SIDE)
        pos.push_back(a);
      else if (sign[i] == CGAL::ON_NEGATIVE_SIDE)
        neg.push_back(a);
      else {
        // Vertex staying on the line.
        vts[a].constraints.insert(line_idx);
        Vertex copy = vts[a];
        vts.push_back(copy);
        pos.push_back(a);
        neg.push_back(vts.size() - 1);
        make_twins(a, vts.size() - 1);
      }

      if (crossing[i].first != std::size_t(-1)) {
        pos.push_back(crossing[i].first);
        neg.push_back(crossing[i].second);
      }
    }

    // Let the positive part start with the first vertex of the positive run (keeps the order of the previous implementation).
    for (std::size_t i = 0; i < n; i++)
      if (sign[i] == CGAL::ON_POSITIVE_SIDE && sign[(i + n - 1) % n] != CGAL::ON_POSITIVE_SIDE) {
        auto it = std::find(pos.begin(), pos.end(), idx[i]);
        pos.splice(pos.end(), pos, pos.begin(), it);
        break;
      }

    poly = std::move(neg);
    return true;
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
      KSP::internal::Log_line(m_log_prefix) << "ERROR: invalid face created!";
      for (std::size_t i = 0; i < f.vertices.size(); i++) {
        KSP::internal::Log_line(m_log_prefix) << "2 " << point_3(f.vertices[i]) << " " << point_3(f.vertices[(i + 1) % f.vertices.size()]);
      }
      KSP::internal::Log_line(m_log_prefix) << "ERROR: end of invalid face";
      exit(1);
    }

    sp.add_face(f_idx, fi);

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

  // Face containing the point p + t * v at t = 0+. The point is strictly inside one face at t = 0+, unless it moves along a line.
  IFace locate_iface_0plus(std::size_t sp_idx, const IkPoint_2& p, const IkVector_2& v) const {
    const Support_plane &sp = m_support_planes[sp_idx];
    for (auto& f : sp.ifaces()) {
      const typename Intersection_graph::Face_property& fp = m_intersection_graph.face(f);
      bool inside = true;
      // pts are oriented ccw
      for (std::size_t i = 0; i < fp.pts.size() && inside; i++) {
        const IkPoint_2& a = fp.pts[i];
        const IkPoint_2& b = fp.pts[(i + 1) % fp.pts.size()];
        CGAL::Orientation o = CGAL::orientation(a, b, p);
        if (o == CGAL::COLLINEAR)
          o = CGAL::Orientation(CGAL::sign((b.x() - a.x()) * v.y() - (b.y() - a.y()) * v.x()));
        inside = (o == CGAL::LEFT_TURN);
      }
      if (inside)
        return f;
    }

    return IFace(-1);
  }

  void clear_pfaces(const std::size_t support_plane_idx) {
    support_plane(support_plane_idx).clear_pfaces();
  }

  PVertex opposite(const PEdge& pedge, const PVertex& pvertex) const {
    if (mesh(pedge.first).target(mesh(pedge.first).halfedge(pedge.second)) == pvertex.second) {
      return PVertex(pedge.first, mesh(pedge.first).source(mesh(pedge.first).halfedge(pedge.second)));
    }
    CGAL_assertion(mesh(pedge).source(mesh(pedge.first).halfedge(pedge.second)) == pvertex.second);
    return PVertex(pedge.first, mesh(pedge.first).target(mesh(pedge.first).halfedge(pedge.second)));
  }

  PEdges_of_pface pedges_of_pface(const PFace& pface) const {
    const auto pedges = PEdges_of_pface(
      boost::make_transform_iterator(halfedges_around_face(halfedge(pface.second, mesh(pface.first)), mesh(pface.first)).begin(), Halfedge_to_pedge(pface.first, mesh(pface.first))),
      boost::make_transform_iterator(halfedges_around_face(halfedge(pface.second, mesh(pface.first)), mesh(pface.first)).end(), Halfedge_to_pedge(pface.first, mesh(pface.first))));
    CGAL_assertion(pedges.size() >= 3);
    return pedges;
  }

  PVertices_of_pface pvertices_of_pface(const PFace& pface) const {
    const auto pvertices = PVertices_of_pface(
      boost::make_transform_iterator(halfedges_around_face(halfedge(pface.second, mesh(pface.first)), mesh(pface.first)).begin(), Halfedge_to_pvertex(pface.first, mesh(pface.first))),
      boost::make_transform_iterator(halfedges_around_face(halfedge(pface.second, mesh(pface.first)), mesh(pface.first)).end(), Halfedge_to_pvertex(pface.first, mesh(pface.first))));
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

  /*******************************
  **          ISimplices        **
  ********************************/

  static IVertex null_ivertex() { return Intersection_graph::null_ivertex(); }
  static IEdge null_iedge() { return Intersection_graph::null_iedge(); }

  decltype(auto) ivertices() const { return m_intersection_graph.vertices(); }
  decltype(auto) iedges() const { return m_intersection_graph.edges(); }

  std::size_t line_idx(IEdge iedge) const { return m_intersection_graph.line(iedge); }
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
    typename Intersection_kernel::Line_3 line;
    auto it = support_planes_idx.begin();
    CGAL_assertion_code(bool intersect =) intersection(m_support_planes[*it++].exact_plane(), m_support_planes[*it++].exact_plane(), line);
    CGAL_assertion(intersect);

    std::vector<std::pair<IkFT, IVertex>> sorted_vertices;
    for (const auto& vertex : vertices) {
      const auto p = m_intersection_graph.point_3(vertex);
      IkFT u = line.to_vector() * (p - line.point());
      sorted_vertices.push_back(std::make_pair(u, vertex));
    }

    std::sort(sorted_vertices.begin(), sorted_vertices.end(),
      [&](const std::pair<IkFT, IVertex>& a, const std::pair<IkFT, IVertex>& b) -> bool {
        return (a.first < b.first);
      }
    );

    // An intersection with another line may coincide with a vertex on the bbox.
    sorted_vertices.erase(std::unique(sorted_vertices.begin(), sorted_vertices.end(),
      [](const std::pair<IkFT, IVertex>& a, const std::pair<IkFT, IVertex>& b) { return a.second == b.second; }), sorted_vertices.end());

    // The line may already exist, e.g., if it lies in a side of the bbox and has been created by intersecting the support
    // planes with the bbox. Then, the existing line is extended by the support planes.
    std::size_t line_idx = std::size_t(-1);
    // Such a line passes through the first vertex, which lies on the boundary of the bbox.
    {
      const IVertex first = sorted_vertices.front().second;
      const IVertex last = sorted_vertices.back().second;
      for (const auto& e : m_intersection_graph.incident_edges(first)) {
        const IVertex other = m_intersection_graph.target(e);
        if (other == first)
          continue;
        if (CGAL::collinear(m_intersection_graph.point_3(first), m_intersection_graph.point_3(other), m_intersection_graph.point_3(last))) {
          line_idx = m_intersection_graph.line(e);
          break;
        }
      }
    }

    if (line_idx != std::size_t(-1)) {
      add_to_existing_line(line_idx, support_planes_idx, sorted_vertices);
      return;
    }

    line_idx = m_intersection_graph.add_line(line);

    for (std::size_t sp_idx : support_planes_idx) {
      Support_plane& sp = support_plane(sp_idx);
      sp.data().lines[line_idx] = typename Intersection_kernel::Line_2(sp.to_2d(m_intersection_graph.point_3(sorted_vertices[0].second)), sp.to_2d(m_intersection_graph.point_3(sorted_vertices.back().second)));
    }

    for (std::size_t i = 0; i < sorted_vertices.size() - 1; ++i) {
      const auto pair = m_intersection_graph.add_edge(
        sorted_vertices[i].second, sorted_vertices[i + 1].second, support_planes_idx);
      const auto iedge = pair.first;
      CGAL_assertion(pair.second);// is inserted?
      m_intersection_graph.set_line(iedge, line_idx);

      for (const auto support_plane_idx : support_planes_idx) {
        support_plane(support_plane_idx).iedges().insert(iedge);
      }
    }

    m_intersection_graph.set_vertices_on_line(line_idx, std::move(sorted_vertices));
  }

  // Adds the support planes to an existing line, e.g., a line in a side of the bbox, and inserts the given vertices into it.
  // The existing edges are split at the new vertices.
  void add_to_existing_line(std::size_t line_idx, const std::set<std::size_t>& support_planes_idx, const std::vector<std::pair<IkFT, IVertex> >& new_vertices) {
    const typename Intersection_kernel::Line_3& line = m_intersection_graph.line(line_idx);
    auto u = [&](IVertex v) { return line.to_vector() * (m_intersection_graph.point_3(v) - line.point()); };

    std::vector<std::pair<IkFT, IVertex> > old_vertices(m_intersection_graph.vertices_on_line(line_idx).begin(), m_intersection_graph.vertices_on_line(line_idx).end());
    std::set<IVertex> all;
    for (const auto& p : old_vertices)
      all.insert(p.second);
    for (const auto& p : new_vertices)
      all.insert(p.second);

    std::vector<std::pair<IkFT, IVertex> > sorted;
    for (IVertex v : all)
      sorted.push_back(std::make_pair(u(v), v));
    std::sort(sorted.begin(), sorted.end(), [](const auto& a, const auto& b) { return a.first < b.first; });

    // Split the existing edges at the vertices lying in their interior.
    for (std::size_t i = 0; i + 1 < old_vertices.size(); i++) {
      const IVertex a = old_vertices[i].second;
      const IVertex b = old_vertices[i + 1].second;
      if (!m_intersection_graph.is_edge(a, b))
        continue;
      const IkFT ua = u(a), ub = u(b);
      std::vector<IVertex> interior;
      for (const auto& p : sorted)
        if ((ua < p.first && p.first < ub) || (ub < p.first && p.first < ua))
          interior.push_back(p.second);
      if (interior.empty())
        continue;
      if (ub < ua)
        std::reverse(interior.begin(), interior.end());

      IEdge edge = m_intersection_graph.edge(a, b);
      const std::set<std::size_t> planes = m_intersection_graph.intersected_planes(edge);
      for (std::size_t plane_idx : planes)
        support_plane(plane_idx).iedges().erase(edge);
      IVertex last = a;
      for (IVertex v : interior) {
        m_intersection_graph.split_edge(edge, v);
        CGAL_assertion(m_intersection_graph.is_edge(last, v) && m_intersection_graph.is_edge(v, b));
        for (std::size_t plane_idx : planes)
          support_plane(plane_idx).iedges().insert(m_intersection_graph.edge(last, v));
        edge = m_intersection_graph.edge(v, b);
        last = v;
      }
      for (std::size_t plane_idx : planes)
        support_plane(plane_idx).iedges().insert(edge);
    }

    // Add missing edges and the support planes to all edges between the new vertices.
    for (std::size_t i = 0; i + 1 < sorted.size(); i++) {
      const IVertex a = sorted[i].second, b = sorted[i + 1].second;
      // Only the part of the line between the vertices of the new line belongs to the support planes.
      if (sorted[i].first < u(new_vertices.front().second) && sorted[i].first < u(new_vertices.back().second))
        continue;
      if (u(new_vertices.front().second) < sorted[i + 1].first && u(new_vertices.back().second) < sorted[i + 1].first)
        continue;
      auto pair = m_intersection_graph.add_edge(a, b, support_planes_idx);
      if (pair.second)
        m_intersection_graph.set_line(pair.first, line_idx);
      CGAL_assertion(m_intersection_graph.line(pair.first) == line_idx);
      for (std::size_t plane_idx : m_intersection_graph.intersected_planes(pair.first))
        support_plane(plane_idx).iedges().insert(pair.first);
    }

    for (std::size_t sp_idx : support_planes_idx) {
      Support_plane& sp = support_plane(sp_idx);
      if (sp.data().lines.find(line_idx) == sp.data().lines.end())
        sp.data().lines[line_idx] = typename Intersection_kernel::Line_2(sp.to_2d(m_intersection_graph.point_3(sorted.front().second)), sp.to_2d(m_intersection_graph.point_3(sorted.back().second)));
    }

    m_intersection_graph.set_vertices_on_line(line_idx, std::move(sorted));
  }

  template<typename IVertexRange>
  void set_vertices_on_line(IVertexRange &vertices, std::size_t line_idx) {
    std::vector<std::pair<IkFT, IVertex>> sorted_vertices;
    sorted_vertices.reserve(vertices.size());

    const typename Intersection_kernel::Line_3 &line = m_intersection_graph.line(line_idx);
    for (const auto& vertex : vertices) {
      const auto p = m_intersection_graph.point_3(vertex);
      IkFT u = line.to_vector() * (p - line.point());
      sorted_vertices.push_back(std::make_pair(u, vertex));
    }

    std::sort(sorted_vertices.begin(), sorted_vertices.end(),
      [&](const std::pair<IkFT, IVertex>& a, const std::pair<IkFT, IVertex>& b) -> bool {
        return (a.first < b.first);
      }
    );

    m_intersection_graph.set_vertices_on_line(line_idx, std::move(sorted_vertices));
  }

  const IVertex source(IEdge edge) const { return m_intersection_graph.source(edge); }
  const IVertex target(IEdge edge) const { return m_intersection_graph.target(edge); }

  const IVertex opposite(IEdge edge, IVertex ivertex) const {
    const auto out = source(edge);
    if (out == ivertex) {
      return target(edge);
    }
    CGAL_assertion(target(edge) == ivertex);
    return out;
  }

  decltype(auto) incident_iedges(IVertex ivertex) const {
    return m_intersection_graph.incident_edges(ivertex);
  }

  const typename Support_plane::IEdge_set& iedges(const std::size_t support_plane_idx) const {
    return support_plane(support_plane_idx).iedges();
  }

  typename Support_plane::IEdge_set& iedges(const std::size_t support_plane_idx) {
    return support_plane(support_plane_idx).iedges();
  }

  const std::set<std::size_t>& intersected_planes(IEdge iedge) const {
    return m_intersection_graph.intersected_planes(iedge);
  }

  const std::set<std::size_t> intersected_planes(IVertex ivertex, const bool keep_bbox = true) const {
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

  /*******************************
  **          STRINGS           **
  ********************************/

  inline const std::string str(IEdge iedge) const {
    std::ostringstream oss; oss << "IEdge" << iedge; return oss.str();
  }

  /*******************************
  **        CONNECTIVITY        **
  ********************************/

  bool has_ivertex(const PVertex& pvertex) const { return support_plane(pvertex.first).has_ivertex(pvertex.second); }
  IVertex ivertex(const PVertex& pvertex) const { return support_plane(pvertex.first).ivertex(pvertex.second); }

  bool has_iedge(const PVertex& pvertex) const { return support_plane(pvertex.first).has_iedge(pvertex.second); }
  IEdge iedge(const PVertex& pvertex) const { return support_plane(pvertex.first).iedge(pvertex.second); }

  bool has_iedge(const PEdge& pedge) const { return support_plane(pedge.first).has_iedge(pedge.second); }
  IEdge iedge(const PEdge& pedge) const { return support_plane(pedge.first).iedge(pedge.second); }

  /*******************************
  **        CONVERSIONS         **
  ********************************/

  IkPoint_2 to_2d(const std::size_t support_plane_idx, IVertex ivertex) const {
    return support_plane(support_plane_idx).to_2d(point_3(ivertex));
  }

  Segment_2 to_2d(const std::size_t support_plane_idx, const Segment_3& segment_3) const {
    return support_plane(support_plane_idx).to_2d(segment_3);
  }

  template <class IK = Intersection_kernel>
  auto to_2d(const std::size_t support_plane_idx, const IkSegment_3& segment_3) const
      -> std::enable_if_t<!std::is_same_v<Kernel, IK>, IkSegment_2> {
    return support_plane(support_plane_idx).to_2d(segment_3);
  }

  Point_2 to_2d(const std::size_t support_plane_idx, const Point_3& point_3) const {
    return support_plane(support_plane_idx).to_2d(point_3);
  }

  template <class IK = Intersection_kernel>
  auto to_2d(const std::size_t support_plane_idx, const IkPoint_3& point_3) const
      -> std::enable_if_t<!std::is_same_v<Kernel, IK>, IkPoint_2> {
    return support_plane(support_plane_idx).to_2d(point_3);
  }

  Point_2 point_2(const PVertex& pvertex) const {
    return support_plane(pvertex).point_2(pvertex.second);
  }

  Point_2 point_2(const std::size_t support_plane_idx, IVertex ivertex) const {
    return support_plane(support_plane_idx).to_2d(from_exact(point_3(ivertex)));
  }

  Point_3 to_3d(const std::size_t support_plane_idx, const Point_2& point_2) const {
    return support_plane(support_plane_idx).to_3d(point_2);
  }

  Point_3 point_3(const PVertex& pvertex) const {
    return support_plane(pvertex.first).point_3(pvertex.second);
  }

  IkPoint_3 point_3(const IVertex& vertex) const {
    return m_intersection_graph.point_3(vertex);
  }

  Segment_3 segment_3(const PEdge& pedge) const {
    return support_plane(pedge.first).segment_3(pedge.second);
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
      const auto q = m_support_planes[sp_idx].to_2d(p);
      polygon.push_back(std::make_pair(q, true));
    }
    CGAL_assertion(polygon.size() == points.size());

    const bool is_valid = support_plane(sp_idx).is_valid_polygon(polygon);

    if (!is_valid) {
      for (const auto& pair : polygon) {
        KSP::internal::Log_line(m_log_prefix) << m_support_planes[sp_idx].to_3d(pair.first);
      }
    }

    CGAL_assertion_msg(is_valid, "ERROR: POLYGON IS NOT VALID!");
    return is_valid;
  }

  bool check_bbox() const {
    for (std::size_t i = 0; i < 6; ++i) {
      const auto pfaces = this->pfaces(i);
      for (const auto pface : pfaces) {
        for (const auto pvertex : pvertices_of_pface(pface)) {
          if (!has_ivertex(pvertex)) {
            KSP::internal::Log_line(m_log_prefix) << "debug pvertex: ";
            CGAL_assertion_msg(has_ivertex(pvertex), "ERROR: BBOX VERTEX IS MISSING AN IVERTEX!");
            return false;
          }
        }
        for (const auto pedge : pedges_of_pface(pface)) {
          if (!has_iedge(pedge)) {
            KSP::internal::Log_line(m_log_prefix) << "debug pedge: ";
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
            KSP::internal::Log_line(m_log_prefix) << "debug pvertex ";
            CGAL_assertion_msg(has_ivertex(pvertex), "ERROR: INTERIOR VERTEX IS MISSING AN IVERTEX!");
            return false;
          }
        }
        for (const auto pedge : pedges_of_pface(pface)) {
          if (!has_iedge(pedge)) {
            KSP::internal::Log_line(m_log_prefix) << "debug pedge ";
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
        KSP::internal::Log_line(m_log_prefix) << "ERROR: CURRENT NUMBER OF EDGES = " << nedges.size();
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
        KSP::internal::Log_line(m_log_prefix) << "ERROR: CURRENT NUMBER OF FACES = " << nfaces.size();
        KSP::internal::Log_line(m_log_prefix) << edge;
        KSP::internal::Log_line(m_log_prefix) << "PFace(" << nfaces[0].first << " " << nfaces[0].second << ")";
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
          KSP::internal::Log_line(m_log_prefix) << "ERROR: CURRENT NUMBER OF VOLUMES = " << nvolumes.size();
          CGAL_assertion_msg(nvolumes.size() == 1 || nvolumes.size() == 2,
            "ERROR: FACE MUST HAVE 1 OR 2 NEIGHBORS!");
          success = false;
        }
      }
    }

    return success;
  }

};

#endif //DOXYGEN_RUNNING

} // namespace internal
} // namespace KSP_3
} // namespace CGAL

#endif // CGAL_KSP_3_DATA_STRUCTURE_H
