// Copyright (c) 2014 GeometryFactory SARL (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Yin Xu, Andreas Fabri and Ilker O. Yaz

#ifndef CGAL_ARAP_WEIGHTS_H
#define CGAL_ARAP_WEIGHTS_H

#include "cotangent_weights.h"

namespace CGAL {

/// \ingroup PkgWeightsRefARAPWeights
///@brief Deformation algorithm type
enum Deformation_algorithm_tag
{
  ORIGINAL_ARAP,   /**< use original as-rigid-as possible algorithm */
  SPOKES_AND_RIMS, /**< use spokes and rims version of as-rigid-as possible algorithm */
  SRE_ARAP         /**< use smooth rotation enhanced As-rigid-as-possible */
};

namespace internal {

template<typename TriangleMesh, typename VertexPointMap,
  Deformation_algorithm_tag deformation_algorithm_tag>
struct Types_selectors;

template<typename TriangleMesh, typename VertexPointMap>
struct Types_selectors<TriangleMesh, VertexPointMap, CGAL::SPOKES_AND_RIMS>
{
  typedef CGAL::Weights::Single_cotangent_weight<TriangleMesh> Weight_calculator;

  struct ARAP_visitor
  {
    void init(const TriangleMesh, VertexPointMap) {}

    void rotation_matrix_pre(typename boost::graph_traits<TriangleMesh>::vertex_descriptor,
      const TriangleMesh&) {}

    template <class Square_matrix_3>
    void update_covariance_matrix(Square_matrix_3&,
      const Square_matrix_3&) {}

    void set_sre_arap_alpha(double) {}
  };
};

template<class TriangleMesh, typename VertexPointMap>
struct Types_selectors<TriangleMesh, VertexPointMap, CGAL::ORIGINAL_ARAP>
{
  typedef CGAL::Weights::Cotangent_weight<TriangleMesh> Weight_calculator;

  typedef typename Types_selectors<TriangleMesh, VertexPointMap, CGAL::SPOKES_AND_RIMS>::ARAP_visitor ARAP_visitor;
};

template<class TriangleMesh, typename VertexPointMap>
struct Types_selectors<TriangleMesh, VertexPointMap, CGAL::SRE_ARAP>
{
  typedef CGAL::Weights::Cotangent_weight<TriangleMesh> Weight_calculator;

  class ARAP_visitor
  {
    double m_nb_edges_incident;
    double m_area;
    double m_alpha;

  public:
    ARAP_visitor() : m_alpha(0.02) {}

    void init(const TriangleMesh triangle_mesh, const VertexPointMap& vpmap)
    {
      // calculate area
      m_area = 0;
      typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;
      for (face_descriptor f : faces(triangle_mesh))
      {
        typename boost::graph_traits<TriangleMesh>::halfedge_descriptor
          h = halfedge(f, triangle_mesh);
        m_area += std::sqrt(CGAL::squared_area(
          get(vpmap, source(h, triangle_mesh)),
          get(vpmap, target(h, triangle_mesh)),
          get(vpmap, target(next(h, triangle_mesh), triangle_mesh))));
      }
    }

    void rotation_matrix_pre(
      typename boost::graph_traits<TriangleMesh>::vertex_descriptor vi,
      const TriangleMesh& hg)
    {
      typename boost::graph_traits<TriangleMesh>::in_edge_iterator e, e_end;
      std::tie(e, e_end) = in_edges(vi, hg);
      m_nb_edges_incident = (double)std::distance(e, e_end);
    }

    template <class Square_matrix_3>
    void update_covariance_matrix(
      Square_matrix_3& cov,
      const Square_matrix_3& rot_mtr)
    {
      // add neighbor rotation
      cov += m_alpha * m_area * rot_mtr.transpose() / m_nb_edges_incident;
    }

    void set_sre_arap_alpha(double a) { m_alpha = a; }
  };
};

}

}

#endif
