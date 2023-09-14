// Copyright (c) 2019 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Dmitry Anisimov
//

#ifndef CGAL_BARYCENTRIC_DELAUNAY_DOMAIN_2_H
#define CGAL_BARYCENTRIC_DELAUNAY_DOMAIN_2_H

#include <CGAL/license/Barycentric_coordinates_2.h>

// CGAL includes.
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>

// Internal includes.
#include <CGAL/Barycentric_coordinates_2/internal/utils_2.h>

namespace CGAL {
namespace Barycentric_coordinates {

  /*!
    \ingroup PkgBarycentricCoordinates2RefHarmonic

    \brief 2D Delaunay domain restricted to a simple polygon.

    This class implements a discretized domain restricted to a simple polygon.
    The interior part of the input polygon is triangulated and refined with respect
    to the user-specified shape size parameter. The final triangulation is a constrained
    Delaunay triangulation, where the constraints are the polygon edges.

    Internally, the package \ref PkgMesh2 is used. See it for more details.

    \tparam VertexRange
    a model of `ConstRange` whose iterator type is `RandomAccessIterator`

    \tparam GeomTraits
    a model of `BarycentricTraits_2`

    \tparam PointMap
    a model of `ReadablePropertyMap` whose key type is `VertexRange::value_type` and
    value type is `Point_2`. The default is `CGAL::Identity_property_map`.

    \cgalModels{DiscretizedDomain_2}
  */
  template<
  typename VertexRange,
  typename GeomTraits,
  typename PointMap = CGAL::Identity_property_map<typename GeomTraits::Point_2> >
  class Delaunay_domain_2 {

  public:

    /// \name Types
    /// @{

    /// \cond SKIP_IN_MANUAL
    using Vertex_range = VertexRange;
    using Geom_traits = GeomTraits;
    using Point_map = PointMap;

    using Construct_centroid_2 = typename GeomTraits::Construct_centroid_2;

    struct VI {
      bool is_on_boundary = false;
      std::size_t index = std::size_t(-1);
      std::vector<std::size_t> neighbors;
    };

    using FB  = CGAL::Delaunay_mesh_face_base_2<GeomTraits>;
    using VB  = CGAL::Triangulation_vertex_base_with_info_2<VI, GeomTraits>;
    using TDS = CGAL::Triangulation_data_structure_2<VB, FB>;
    using TAG = CGAL::Exact_predicates_tag;
    using CDT = CGAL::Constrained_Delaunay_triangulation_2<GeomTraits, TDS, TAG>;

    using Vertex_handle = typename CDT::Vertex_handle;

    using Criteria = CGAL::Delaunay_mesh_size_criteria_2<CDT>;
    using Mesher   = CGAL::Delaunay_mesher_2<CDT, Criteria>;
    /// \endcond

    /// Number type.
    typedef typename GeomTraits::FT FT;

    /// Point type.
    typedef typename GeomTraits::Point_2 Point_2;

    /// @}

    /// \name Initialization
    /// @{

    /*!
      \brief initializes all internal data structures.

      \param polygon
      an instance of `VertexRange` with the vertices of a simple polygon

      \param traits
      a traits class with geometric objects, predicates, and constructions;
      the default initialization is provided

      \param point_map
      an instance of `PointMap` that maps a vertex from `polygon` to `Point_2`;
      the default initialization is provided

      \pre polygon.size() >= 3
      \pre polygon is simple
    */
    Delaunay_domain_2(
      const VertexRange& polygon,
      const GeomTraits traits = GeomTraits(),
      const PointMap point_map = PointMap()) :
    m_polygon(polygon),
    m_traits(traits),
    m_point_map(point_map),
    m_construct_centroid_2(m_traits.construct_centroid_2_object()) {

      CGAL_precondition(
        polygon.size() >= 3);
      CGAL_precondition(
        internal::is_simple_2(polygon, traits, point_map));
      clear();
    }

    /*!
      \brief creates a refined Delaunay triangulation restricted to the input polygon.

      After the construction is completed, the first n vertices are the polygon vertices,
      the next m vertices are the vertices generated along the polygon boundary,
      the last k vertices are the vertices generated in the interior part of the polygon.

      \tparam PointRange
      a model of `ConstRange` whose value type is `GeomTraits::Point_2`

      \param max_edge_length
      an upper bound on the length of the longest edge;
      see `Delaunay_mesh_size_criteria_2` for more details

      \param seeds
      contains seed points indicating, which parts of the input polygon
      should be partitioned and subdivided
    */
    template<typename PointRange>
    void create(
      const FT max_edge_length, const PointRange& seeds) {

      create_triangulation();
      refine_triangulation(
        max_edge_length, seeds);
      check_boundaries();
      create_neighbors();
    }

    /// @}

    /// \name Access
    /// @{

    /*!
      \brief computes barycenters of all generated triangles.

      \tparam OutIterator
      a model of `OutputIterator` that accepts points of type `Point_2`

      \param b_begin
      the beginning of the destination range with the computed barycenters

      \return an output iterator to the element in the destination range,
      one past the last barycenter stored
    */
    template<typename OutIterator>
    OutIterator barycenters(OutIterator b_begin) const {

      const std::size_t num_faces = get_number_of_faces();
      if (num_faces == 0) return b_begin;

      for (auto fh = m_cdt.finite_faces_begin();
      fh != m_cdt.finite_faces_end(); ++fh) {
        if (!fh->is_in_domain()) continue;

        const Point_2 b = m_construct_centroid_2(
        fh->vertex(0)->point(),
        fh->vertex(1)->point(),
        fh->vertex(2)->point());
        *(b_begin++) = b;
      }
      return b_begin;
    }

    /*!
      \brief returns the number of triangulation vertices.

      This function implements `DiscretizedDomain_2::number_of_vertices()`.
    */
    std::size_t number_of_vertices() const {

      CGAL_assertion(
        m_vhs.size() == m_cdt.number_of_vertices());
      return m_vhs.size();
    }

    /*!
      \brief returns a const reference to the triangulation vertex with
      the index `query_index`.

      This function implements `DiscretizedDomain_2::vertex()`.

      \param query_index
      an index of the requested vertex

      \pre query_index >= 0 && query_index < number_of_vertices()
    */
    const Point_2& vertex(
      const std::size_t query_index) const {

      CGAL_precondition(query_index < number_of_vertices());
      return m_vhs[query_index]->point();
    }

    /*!
      \brief controls if the triangulation vertex with the index `query_index`
      is on the polygon boundary.

      This function implements `DiscretizedDomain_2::is_on_boundary()`.

      \param query_index
      an index of the query vertex

      \pre query_index >= 0 && query_index < number_of_vertices()
    */
    bool is_on_boundary(
      const std::size_t query_index) const {

      CGAL_precondition(query_index < number_of_vertices());
      return m_vhs[query_index]->info().is_on_boundary;
    }

    /*!
      \brief returns the one-ring neighborhood of the triangulation vertex
      with the index `query_index`.

      This function implements `DiscretizedDomain_2::operator()()`.

      \param query_index
      an index of the query vertex

      \param neighbors
      stores indices of the vertices, which from the one-ring neighborhood
      of the query vertex

      \pre query_index >= 0 && query_index < number_of_vertices()
    */
    void operator()(
      const std::size_t query_index,
      std::vector<std::size_t>& neighbors) const {

      CGAL_precondition(query_index < number_of_vertices());
      const auto vh = m_vhs[query_index];
      neighbors = vh->info().neighbors;
    }

    /*!
      \brief locates a triangle that contains a given query point.

      If `triangle` is empty, the query point does not belong to the domain.

      This function implements `DiscretizedDomain_2::locate()`.

      \param query
      a query point

      \param triangle
      stores indices of the vertices, which form a triangle, that contains
      the query point
    */
    void locate(
      const Point_2& query,
      std::vector<std::size_t>& triangle) const {

      triangle.clear();
      const auto fh = m_cdt.locate(query);
      if (fh->is_in_domain()) {
        for (int i = 0; i < 3; ++i) {
          triangle.push_back(fh->vertex(i)->info().index);
        }
      }
    }

    /// @}

    /// \name Memory Management
    /// @{

    /*!
      \brief clears all internal data structures.
    */
    void clear() {
      m_vhs.clear();
      m_cdt.clear();
    }

    /*!
      \brief releases all memory that is used internally.
    */
    void release_memory() {
      clear();
      m_vhs.shrink_to_fit();
    }

    /// @}

    /// \cond SKIP_IN_MANUAL
    void export_points_2(
      const std::vector<Point_2>& points,
      const std::string file_path) const {

      std::stringstream out;
      out.precision(20);
      const std::size_t num_vertices = points.size();
      add_ply_header_points(num_vertices, out);

      for (const auto& point : points) {
        out << point << " 0 0 0 0" << std::endl;
      }
      save(out, file_path + ".ply");
    }

    template<typename Point_3>
    void export_points_3(
      const std::vector<Point_3>& points,
      const std::string file_path) const {

      std::stringstream out;
      out.precision(20);
      const std::size_t num_vertices = points.size();
      add_ply_header_points(num_vertices, out);

      for (const auto& point : points) {
        out << point << " 0 0 0" << std::endl;
      }
      save(out, file_path + ".ply");
    }

    void export_triangulation(
      const std::string file_path) const {

      std::stringstream out;
      out.precision(20);
      const std::size_t num_faces = get_number_of_faces();
      add_ply_header_mesh(num_faces * 3, num_faces, out);

      for (auto fh = m_cdt.finite_faces_begin();
      fh != m_cdt.finite_faces_end(); ++fh) {
        if (!fh->is_in_domain()) continue;

        const auto& p0 = fh->vertex(0)->point();
        const auto& p1 = fh->vertex(1)->point();
        const auto& p2 = fh->vertex(2)->point();

        out << p0.x() << " " << p0.y() << " 0" << std::endl;
        out << p1.x() << " " << p1.y() << " 0" << std::endl;
        out << p2.x() << " " << p2.y() << " 0" << std::endl;
      }

      std::size_t i = 0;
      for (auto fh = m_cdt.finite_faces_begin();
      fh != m_cdt.finite_faces_end(); ++fh) {
        if (!fh->is_in_domain()) continue;

        out << 3 << " "
        << i + 0 << " " << i + 1 << " " << i + 2 << " "
        << "0 0 0" << std::endl;
        i += 3;
      }
      save(out, file_path + ".ply");
    }
    /// \endcond

  private:

    // Fields.
    const VertexRange& m_polygon;
    const GeomTraits m_traits;
    const PointMap m_point_map;

    const Construct_centroid_2 m_construct_centroid_2;

    CDT m_cdt;
    std::vector<Vertex_handle> m_vhs;

    void create_triangulation() {

      const std::size_t n = m_polygon.size();
      m_cdt.clear(); m_vhs.clear();
      m_vhs.reserve(n);

      for (std::size_t i = 0; i < n; ++i) {
        const auto& p = get(m_point_map, *(m_polygon.begin() + i));
        m_vhs.push_back(m_cdt.insert(p));
      }

      CGAL_assertion(m_vhs.size() == n);
      CGAL_assertion(m_cdt.number_of_vertices() == n);

      for (std::size_t i = 0; i < n; ++i) {
        const std::size_t ip = (i + 1) % n;
        if (m_vhs[i] != m_vhs[ip]) {
          m_cdt.insert_constraint(m_vhs[i], m_vhs[ip]);
        }
      }
    }

    template<typename PointRange>
    void refine_triangulation(
      const FT max_edge_length, const PointRange& seeds) {

      // 0.125 is the default shape bound that corresponds to a bound of 20.6 degrees.
      Mesher mesher(m_cdt);
      mesher.set_seeds(seeds.begin(), seeds.end(), true);
      mesher.set_criteria(Criteria(0.125, max_edge_length));
      mesher.refine_mesh();

      m_vhs.clear();
      m_vhs.reserve(m_cdt.number_of_vertices());

      std::size_t count = 0;
      for (auto vh = m_cdt.finite_vertices_begin();
      vh != m_cdt.finite_vertices_end(); ++vh, ++count) {
        vh->info().index = count;
        m_vhs.push_back(vh);
      }
      CGAL_assertion(m_vhs.size() == m_cdt.number_of_vertices());
    }

    void check_boundaries() {

      for (std::size_t i = 0; i < m_vhs.size(); ++i) {
        const auto vh = m_vhs[i];
        vh->info().is_on_boundary = false;

        auto face = m_cdt.incident_faces(vh);
        CGAL_assertion(!face.is_empty());
        const auto end = face;

        do {
          if (!face->is_in_domain()) {
            vh->info().is_on_boundary = true; break;
          } ++face;
        } while (face != end);
      }
    }

    void create_neighbors() {

      for (std::size_t i = 0; i < m_vhs.size(); ++i) {
        const auto vh = m_vhs[i];

        auto edge = m_cdt.incident_edges(vh);
        CGAL_assertion(!edge.is_empty());
        auto end = edge;

        // Move the pointer to the most left boundary face.
        if (vh->info().is_on_boundary) {
          do {
            const auto fh = edge->first;
            if (!fh->is_in_domain()) {
              end = edge; break;
            } ++edge;
          } while (edge != end);
          do {
            const auto fh = edge->first;
            if (fh->is_in_domain()) {
              end = edge; break;
            } ++edge;
          } while (edge != end);
        }

        // We then traverse edges in the counterclockwise order.
        vh->info().neighbors.clear();
        do {
          const auto fh = edge->first;
          const auto nh = fh->vertex(edge->second);

          if (fh->is_in_domain()) {
            vh->info().neighbors.push_back(nh->info().index);
          } else {
            vh->info().neighbors.push_back(nh->info().index);
            break;
          } ++edge;
        } while (edge != end);
        CGAL_assertion(
          vh->info().neighbors.size() > 0);
      }
    }

    std::size_t get_number_of_faces() const {

      std::size_t num_faces = 0;
      for (auto fh = m_cdt.finite_faces_begin();
      fh != m_cdt.finite_faces_end(); ++fh) {
        if (!fh->is_in_domain()) continue;
        ++num_faces;
      }
      return num_faces;
    }

    void add_ply_header_points(
      const std::size_t num_vertices,
      std::stringstream& out) const {

      out <<
      "ply" << std::endl <<
      "format ascii 1.0" << std::endl <<
      "element vertex " << num_vertices << std::endl <<
      "property double x" << std::endl <<
      "property double y" << std::endl <<
      "property double z" << std::endl <<
      "property uchar red"   << std::endl <<
      "property uchar green" << std::endl <<
      "property uchar blue"  << std::endl <<
      "end_header" << std::endl;
    }

    void add_ply_header_mesh(
      const std::size_t num_vertices,
      const std::size_t num_faces,
      std::stringstream& out) const {

      out <<
      "ply" << std::endl <<
      "format ascii 1.0" << std::endl <<
      "element vertex " << num_vertices << std::endl <<
      "property double x" << std::endl <<
      "property double y" << std::endl <<
      "property double z" << std::endl <<
      "element face " << num_faces << std::endl <<
      "property list uchar int vertex_indices" << std::endl <<
      "property uchar red"   << std::endl <<
      "property uchar green" << std::endl <<
      "property uchar blue"  << std::endl <<
      "end_header" << std::endl;
    }

    void save(
      const std::stringstream& out,
      const std::string file_path) const {

      std::ofstream file(
        file_path.c_str(), std::ios_base::out);
      if (!file) {
        std::cerr << std::endl <<
          "ERROR: error saving file " << file_path
        << "!" << std::endl << std::endl;
      }
      file << out.str();
      file.close();
    }
  };

} // namespace Barycentric_coordinates
} // namespace CGAL

#endif // CGAL_BARYCENTRIC_DELAUNAY_DOMAIN_2_H
