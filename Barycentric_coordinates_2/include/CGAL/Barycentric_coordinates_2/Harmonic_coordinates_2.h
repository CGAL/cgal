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
// Author(s)     : Dmitry Anisimov, Kai Hormann
//

#ifndef CGAL_BARYCENTRIC_HARMONIC_COORDINATES_2_H
#define CGAL_BARYCENTRIC_HARMONIC_COORDINATES_2_H

#include <CGAL/license/Barycentric_coordinates_2.h>

// Eigen includes.
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

// Internal includes.
#include <CGAL/Weights/cotangent_weights.h>
#include <CGAL/Barycentric_coordinates_2/internal/utils_2.h>

// [1] Reference: "P. Joshi, M. Meyer, T. DeRose, B. Green, and T. Sanocki.
// Harmonic coordinates for character articulation.
// ACM Transactions on Graphics, 26(3):71:1-9, 2007.".

namespace CGAL {
namespace Barycentric_coordinates {

  /*!
    \ingroup PkgBarycentricCoordinates2RefHarmonic

    \brief 2D harmonic coordinates.

    This class implements 2D harmonic coordinate functions ( \cite cgal:bc:fhk-gcbcocp-06,
    \cite cgal:bc:jmdgs-hcfca-07 ), which can be evaluated at any point inside a
    simple polygon.

    Harmonic coordinates are well-defined and non-negative in the closure
    of any simple polygon, however they cannot be computed analytically and hence
    they are approximated. The classical way to approximate these coordinates is
    by discretizing over the space of piecewise linear functions with respect to
    a partition of the polygon's interior domain, e.g. a triangulation.

    Once computed at the vertices of the discretized domain, the coordinate functions
    can be evaluated at any point inside a polygon by locating the finite element that
    contains a query point and linearly interpolating within this element. See more details
    in the user manual \ref compute_hm_coord "here".

    \tparam VertexRange
    a model of `ConstRange` whose iterator type is `RandomAccessIterator`

    \tparam DiscretizedDomain
    a model of `DiscretizedDomain_2`. For the moment, we only support domains
    whose partition's finite elements are triangles.

    \tparam GeomTraits
    a model of `BarycentricTraits_2`

    \tparam PointMap
    a model of `ReadablePropertyMap` whose key type is `VertexRange::value_type` and
    value type is `Point_2`. The default is `CGAL::Identity_property_map`.
  */
  template<
  typename VertexRange,
  typename DiscretizedDomain,
  typename GeomTraits,
  typename PointMap = CGAL::Identity_property_map<typename GeomTraits::Point_2> >
  class Harmonic_coordinates_2 {

  public:

    /// \name Types
    /// @{

    /// \cond SKIP_IN_MANUAL
    using Vertex_range = VertexRange;
    using Discretized_domain = DiscretizedDomain;
    using Geom_traits = GeomTraits;
    using Point_map = PointMap;
    /// \endcond

    /// Number type.
    typedef typename GeomTraits::FT FT;

    /// Point type.
    typedef typename GeomTraits::Point_2 Point_2;

    /// \cond SKIP_IN_MANUAL
    using Construct_vector_2 = typename GeomTraits::Construct_vector_2;

    using VectorFT  = Eigen::Matrix<FT, Eigen::Dynamic, Eigen::Dynamic>;
    using MatrixFT  = Eigen::SparseMatrix<FT>;
    using TripletFT = Eigen::Triplet<FT>;

    #if defined(HMC_SparseLU)
    using Solver = Eigen::SparseLU<MatrixFT>;
    #elif defined(HMC_SimplicialLLT)
    using Solver = Eigen::SimplicialLLT<MatrixFT>;
    #else
    using Solver = Eigen::SimplicialLDLT<MatrixFT>;
    #endif
    /// \endcond

    /// @}

    /// \name Initialization
    /// @{

    /*!
      \brief initializes all internal data structures.

      This class implements the behavior of harmonic coordinates
      for 2D query points.

      \param polygon
      an instance of `VertexRange` with the vertices of a simple polygon

      \param domain
      an instance of `DiscretizedDomain` with a partition of the interior part of a simple polygon

      \param traits
      a traits class with geometric objects, predicates, and constructions;
      the default initialization is provided

      \param point_map
      an instance of `PointMap` that maps a vertex from `polygon` to `Point_2`;
      the default initialization is provided

      \pre polygon.size() >= 3
      \pre polygon is simple
    */
    Harmonic_coordinates_2(
      const VertexRange& polygon,
      const DiscretizedDomain& domain,
      const GeomTraits traits = GeomTraits(),
      const PointMap point_map = PointMap()) :
    m_polygon(polygon),
    m_domain(domain),
    m_traits(traits),
    m_point_map(point_map),
    m_construct_vector_2(m_traits.construct_vector_2_object()) {

      CGAL_precondition(
        polygon.size() >= 3);
      CGAL_precondition(
        internal::is_simple_2(polygon, traits, point_map));
      clear();
    }

    /// @}

    /// \name Access
    /// @{

    /*!
      \brief evaluates 2D harmonic coordinates.

      This function fills `c_begin` with harmonic coordinates evaluated at the `query`
      point with respect to the vertices of the input polygon. Evaluation is performed
      by locating the finite element in the input domain that contains `query` and then
      linearly interpolating harmonic coordinates within this element.

      If `query` does not belong to the input domain or the located element has more than
      3 vertices, all coordinates are set to zero.

      The number of returned coordinates equals to the number of polygon vertices.

      After the coordinates \f$b_i\f$ with \f$i = 1\dots n\f$ are computed, where
      \f$n\f$ is the number of polygon vertices, the query point \f$q\f$ can be obtained
      as \f$q = \sum_{i = 1}^{n}b_ip_i\f$, where \f$p_i\f$ are the polygon vertices.

      \tparam OutIterator
      a model of `OutputIterator` that accepts values of type `FT`

      \param query
      a query point

      \param c_begin
      the beginning of the destination range with the computed coordinates

      \return an output iterator to the element in the destination range,
      one past the last coordinate stored
    */
    template<typename OutIterator>
    OutIterator operator()(
      const Point_2& query, OutIterator c_begin) {

      CGAL_precondition(
        m_setup_is_called &&
        m_factorize_is_called &&
        m_solve_is_called);
      if (!(
        m_setup_is_called &&
        m_factorize_is_called &&
        m_solve_is_called)) return c_begin;

      const std::size_t n = m_polygon.size();

      m_element.clear();
      m_domain.locate(query, m_element);
      if (m_element.size() != 3) {
        internal::get_default(n, c_begin);
        return c_begin;
      }

      const std::size_t i0 = m_element[0];
      const std::size_t i1 = m_element[1];
      const std::size_t i2 = m_element[2];

      CGAL_assertion(i0 < m_domain.number_of_vertices());
      CGAL_assertion(i1 < m_domain.number_of_vertices());
      CGAL_assertion(i2 < m_domain.number_of_vertices());

      const auto& p0 = m_domain.vertex(i0);
      const auto& p1 = m_domain.vertex(i1);
      const auto& p2 = m_domain.vertex(i2);

      m_coordinates.clear();
      internal::planar_coordinates_2(
        p0, p1, p2, query, std::back_inserter(m_coordinates), m_traits);
      CGAL_assertion(m_coordinates.size() == 3);

      const auto& b = m_coordinates;
      CGAL_assertion(b[0] >= FT(0) && b[0] <= FT(1));
      CGAL_assertion(b[1] >= FT(0) && b[1] <= FT(1));
      CGAL_assertion(b[2] >= FT(0) && b[2] <= FT(1));

      CGAL_assertion(m_boundary.size() > 0);
      CGAL_assertion(m_interior.size() > 0);

      FT hm0 = FT(0), hm1 = FT(0), hm2 = FT(0);
      for (std::size_t k = 0; k < n; ++k) {
        if (m_domain.is_on_boundary(i0)) hm0 = m_boundary(m_indices[i0], k);
        else hm0 = m_interior(m_indices[i0], k);
        if (m_domain.is_on_boundary(i1)) hm1 = m_boundary(m_indices[i1], k);
        else hm1 = m_interior(m_indices[i1], k);
        if (m_domain.is_on_boundary(i2)) hm2 = m_boundary(m_indices[i2], k);
        else hm2 = m_interior(m_indices[i2], k);

        CGAL_assertion(hm0 >= FT(0) && hm0 <= FT(1));
        CGAL_assertion(hm1 >= FT(0) && hm1 <= FT(1));
        CGAL_assertion(hm2 >= FT(0) && hm2 <= FT(1));
        *(c_begin++) = hm0 * b[0] + hm1 * b[1] + hm2 * b[2];
      }
      return c_begin;
    }

    /*!
      \brief returns 2D harmonic coordinates at one domain vertex.

      This function fills `c_begin` with harmonic coordinates computed at the
      vertex of the input domain with the index `query_index`.

      The number of returned coordinates equals to the number of polygon vertices.

      After the coordinates \f$b_i\f$ with \f$i = 1\dots n\f$ are computed, where
      \f$n\f$ is the number of polygon vertices, the partition vertex \f$q\f$ with
      the index `query_index` can be obtained as \f$q = \sum_{i = 1}^{n}b_ip_i\f$,
      where \f$p_i\f$ are the polygon vertices.

      \tparam OutIterator
      a model of `OutputIterator` that accepts values of type `FT`

      \param query_index
      a domain's vertex index

      \param c_begin
      the beginning of the destination range with the computed coordinates

      \return an output iterator to the element in the destination range,
      one past the last coordinate stored

      \pre query_index >= 0 && query_index < domain.number_of_vertices()
    */
    template<typename OutIterator>
    OutIterator operator()(
      const std::size_t query_index, OutIterator c_begin) {

      CGAL_precondition(
        m_setup_is_called &&
        m_factorize_is_called &&
        m_solve_is_called);
      if (!(
        m_setup_is_called &&
        m_factorize_is_called &&
        m_solve_is_called)) return c_begin;

      CGAL_precondition(
        query_index < m_domain.number_of_vertices());
      CGAL_assertion(m_boundary.size() > 0);
      CGAL_assertion(m_interior.size() > 0);

      // Save harmonic coordinates.
      const std::size_t n = m_polygon.size();
      if (m_domain.is_on_boundary(query_index)) {
        for (std::size_t k = 0; k < n; ++k) {
          *(c_begin++) = m_boundary(m_indices[query_index], k);
        }
      } else {
        for (std::size_t k = 0; k < n; ++k) {
          *(c_begin++) = m_interior(m_indices[query_index], k);
        }
      }
      return c_begin;
    }

    /*!
      \brief returns 2D harmonic coordinates at all domain vertices.

      This function fills `c_begin` with harmonic coordinates computed at
      the vertices of the input domain.

      The number of returned coordinates equals to the number of input domain vertices.

      \tparam OutIterator
      a model of `OutputIterator` that accepts values of type `std::vector<FT>`

      \param c_begin
      the beginning of the destination range with the computed coordinates

      \return an output iterator to the element in the destination range,
      one past the last coordinate set stored
    */
    template<typename OutIterator>
    OutIterator operator()(OutIterator c_begin) {

      std::vector<FT> coordinates;
      coordinates.reserve(m_polygon.size());
      for (std::size_t k = 0; k < m_domain.number_of_vertices(); ++k) {
        coordinates.clear();
        operator()(k, std::back_inserter(coordinates));
        *(c_begin++) = coordinates;
      }
      return c_begin;
    }

    /// @}

    /// \name Computation
    /// @{

    /*!
      \brief computes 2D harmonic coordinates at the vertices of the input domain.
    */
    void compute() {
      setup(); // compute harmonic data
      factorize(); // factorize the matrix A
      solve(); // solve the linear system Ax = b
    }

    /// @}

    /// \name Step by Step Computation
    /// @{

    /*!
      \brief computes all necessary harmonic data.

      This function fills in the left side matrix A and the right side vector b
      of the linear system Ax = b. The matrix A is a sparse symmetric positive
      definite matrix. The solution vector x of this system gives the harmonic
      coordinates at the vertices of the input domain.
    */
    void setup() {

      if (m_setup_is_called) return;
      const std::size_t n = m_polygon.size();

      // Create an index map. It splits interior and boundary vertices.
      const auto pair = create_indices(m_indices);

      // Initialize all containers.
      const std::size_t numB = pair.first;
      m_boundary = VectorFT(numB, n); // boundary

      const std::size_t numI = pair.second;
      m_interior = VectorFT::Zero(numI, n); // interior

      m_A = MatrixFT(numI, numI); // a sparse matrix
      m_b = VectorFT::Zero(numI, n); // boundary conditions

      // Compute harmonic coordinates.
      set_boundary_vector(
        numB, m_indices, m_boundary);
      set_harmonic_data(
        numI, m_indices, m_boundary, m_A, m_b);
      m_setup_is_called = true;
    }

    /*!
      \brief factorizes the matrix A.

      \pre setup() is called
    */
    void factorize() {
      if (m_factorize_is_called) return;
      CGAL_precondition(m_setup_is_called);
      m_solver_ptr = std::make_shared<Solver>();
      m_solver_ptr->compute(m_A);
      m_factorize_is_called = true;
    }

    /*!
      \brief solves the linear system Ax = b.

      \pre factorize() is called
    */
    void solve() {
      if (m_solve_is_called) return;
      CGAL_precondition(m_factorize_is_called);
      m_interior = m_solver_ptr->solve(m_b);
      m_A.resize(0, 0); m_b.resize(0, 0);
      m_solver_ptr = nullptr;
      m_solve_is_called = true;
    }

    /// @}

    /// \name Memory Management
    /// @{

    /*!
      \brief clears all internal data structures.
    */
    void clear() {
      m_indices.clear();
      m_element.clear();
      m_coordinates.clear();
      m_boundary.resize(0, 0);
      m_interior.resize(0, 0);
      m_setup_is_called = false;
      m_factorize_is_called = false;
      m_solve_is_called = false;
    }

    /*!
      \brief releases all memory that is used internally.
    */
    void release_memory() {
      clear();
      m_indices.shrink_to_fit();
      m_element.shrink_to_fit();
      m_coordinates.shrink_to_fit();
    }

    /// @}

  private:

    // Fields.
    const VertexRange& m_polygon;
    const DiscretizedDomain& m_domain;
    const GeomTraits m_traits;
    const PointMap m_point_map;

    const Construct_vector_2 m_construct_vector_2;

    // Indices of the finite element.
    std::vector<std::size_t> m_element;

    // Barycentric coordinates of the query point.
    std::vector<FT> m_coordinates;

    // Harmonic coordinates are stored separately
    // for interior and boundary vertices.
    VectorFT m_interior, m_boundary;

    // Splits boundary and interior vertices.
    std::vector<std::size_t> m_indices;

    // Tags.
    bool m_setup_is_called;
    bool m_factorize_is_called;
    bool m_solve_is_called;

    // Temporary data.
    MatrixFT m_A;
    VectorFT m_b;

    // Sparse solver;
    std::shared_ptr<Solver> m_solver_ptr;

    std::pair<std::size_t, std::size_t> create_indices(
      std::vector<std::size_t>& indices) const {

      const std::size_t N = m_domain.number_of_vertices();

      // Create an index map.
      indices.clear();
      indices.reserve(N);
      std::size_t numB = 0, numI = 0;
      for (std::size_t i = 0; i < N; ++i) {
        if (m_domain.is_on_boundary(i)) {
          indices.push_back(numB); ++numB;
        } else {
          indices.push_back(numI); ++numI;
        }
      }

      CGAL_assertion(indices.size() == N);
      return std::make_pair(numB, numI);
    }

    void set_boundary_vector(
      const std::size_t /* numB */, // this is not used here but the number of elements is numB!
      const std::vector<std::size_t>& indices,
      VectorFT& boundary) const {

      const std::size_t n = m_polygon.size();
      const std::size_t N = m_domain.number_of_vertices();

      // Initialize temporary containers.
      // Can I remove this lambda?
      std::vector<FT> lambda;
      lambda.reserve(n);

      // Traverse boundary vertices of the domain.
      for (std::size_t i = 0; i < N; ++i) {
        if (m_domain.is_on_boundary(i)) {
          const auto& query = m_domain.vertex(i);

          // Find index of the polygon edge that contains the boundary vertex.
          const auto edge_is_found =
            internal::get_edge_index_approximate(
              m_polygon, query, m_traits, m_point_map);
          CGAL_assertion(static_cast<bool>(edge_is_found));
          const auto location = (*edge_is_found).first;
          const auto index = (*edge_is_found).second;

          // Compute barycentric boundary coordinates.
          lambda.clear();
          internal::boundary_coordinates_2(
            m_polygon, query, location, index,
            std::back_inserter(lambda), m_traits, m_point_map);

          // Set boundary vector.
          for (std::size_t k = 0; k < n; ++k) {
            boundary(indices[i], k) = lambda[k];
          }
        }
      }
    }

    void set_harmonic_data(
      const std::size_t numI,
      const std::vector<std::size_t>& indices,
      const VectorFT& boundary,
      MatrixFT& A,
      VectorFT& b) const {

      const std::size_t n = m_polygon.size();
      const std::size_t N = m_domain.number_of_vertices();

      // Initialize temporary containers.
      // 7 is an average number of neighbors.
      std::vector<TripletFT> triplet_list;
      triplet_list.reserve(numI * 7);

      std::vector<std::size_t> neighbors;
      neighbors.reserve(7);

      // Traverse interior vertices of the domain.
      for (std::size_t i = 0; i < N; ++i) {
        if (!m_domain.is_on_boundary(i)) {
          const auto& query = m_domain.vertex(i);

          // Find one-ring neighborhood of the interior vertex.
          neighbors.clear();
          m_domain(i, neighbors);
          const std::size_t nn = neighbors.size(); // nn is about 7
          CGAL_assertion(nn > 0);

          // Set the right side vector b of the system
          // and fill in a triplet list.
          FT W = FT(0);
          for (std::size_t j = 0; j < nn; ++j) {
            const std::size_t jm  = (j + nn - 1) % nn;
            const std::size_t jp  = (j + 1) % nn;

            const std::size_t idx = neighbors[j];
            const auto& p0 = m_domain.vertex(neighbors[jm]);
            const auto& p1 = m_domain.vertex(neighbors[j]);
            const auto& p2 = m_domain.vertex(neighbors[jp]);
            const FT w = -Weights::cotangent_weight(
              p0, p1, p2, query, m_traits) / FT(2);
            W -= w;

            if (m_domain.is_on_boundary(idx)) {
              for (std::size_t k = 0; k < n; ++k) {
                b(indices[i], k) -= boundary(indices[idx], k) * w;
              }
            } else {
              triplet_list.push_back(
                TripletFT(int(indices[i]), int(indices[idx]), w));
            }
          }
          triplet_list.push_back(
            TripletFT(int(indices[i]), int(indices[i]), W));
        }
      }

      // Set the sparse matrix A. The left side of the system.
      A.setFromTriplets(
        triplet_list.begin(), triplet_list.end());
      A.makeCompressed();
    }
  };

} // namespace Barycentric_coordinates
} // namespace CGAL

#endif // CGAL_BARYCENTRIC_HARMONIC_COORDINATES_2_H
