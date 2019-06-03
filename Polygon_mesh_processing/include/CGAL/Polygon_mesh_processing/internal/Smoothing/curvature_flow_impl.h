// Copyright (c) 2018 GeometryFactory (France).
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
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Konstantinos Katrioplas (konst.katrioplas@gmail.com)

#ifndef CGAL_POLYGON_MESH_PROCESSING_INTERNAL_CURVATURE_FLOW_IMPL_H
#define CGAL_POLYGON_MESH_PROCESSING_INTERNAL_CURVATURE_FLOW_IMPL_H

#include <CGAL/license/Polygon_mesh_processing/meshing_hole_filling.h>

#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/Weights.h>
#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>
#include <CGAL/Polygon_mesh_processing/internal/Smoothing/smoothing_helpers.h>
#include <CGAL/utility.h>

#if defined(CGAL_EIGEN3_ENABLED)
#include <CGAL/Eigen_solver_traits.h>
#endif

#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>

#include <algorithm>
#include <iostream>
#include <utility>
#include <vector>

namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {

template<typename TriangleMesh,
         typename VertexPointMap,
         typename VertexConstraintMap,
         typename SparseLinearSolver,
         typename GeomTraits>
class Shape_smoother
{
  typedef typename GeomTraits::FT                                                 NT;
  typedef typename GeomTraits::Point_3                                            Point;
  typedef typename boost::property_traits<VertexPointMap>::reference              Point_ref;
  typedef CGAL::Triple<int, int, double>                                          Triplet;

  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor           vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor         halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor             edge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor             face_descriptor;
  typedef typename boost::property_map<TriangleMesh, boost::vertex_index_t>::type IndexMap;

  // linear system
  typedef typename SparseLinearSolver::Matrix                                     Eigen_matrix;
  typedef typename SparseLinearSolver::Vector                                     Eigen_vector;

public:
  Shape_smoother(TriangleMesh& mesh,
                 VertexPointMap& vpmap,
                 VertexConstraintMap& vcmap)
    :
      mesh_(mesh),
      vpmap_(vpmap),
      vcmap_(vcmap),
      vimap_(get(boost::vertex_index, mesh_)),
      weight_calculator_(mesh, vpmap)
  {}

  template<typename FaceRange>
  void init_smoothing(const FaceRange& face_range)
  {
    set_face_range(face_range);
    nb_vert_ = vrange_.size();
  }

  void setup_system(Eigen_matrix& A,
                    Eigen_vector& bx, Eigen_vector& by, Eigen_vector& bz,
                    std::vector<Triplet>& stiffness_elements,
                    const double& time)
  {
    compute_coefficient_matrix(A, stiffness_elements, time);
    compute_rhs(bx, by, bz);
  }

  void solve_system(const Eigen_matrix& A,
                    Eigen_vector& Xx, Eigen_vector& Xy, Eigen_vector& Xz,
                    const Eigen_vector& bx, const Eigen_vector& by, const Eigen_vector& bz,
                    SparseLinearSolver& solver)
  {
    NT D;
    // calls compute once to factorize with the preconditioner
    if(!solver.factor(A, D))
    {
      std::cerr << "Could not factorize linear system with preconditioner." << std::endl;
      return;
    }

    if(!solver.linear_solver(bx, Xx) ||
       !solver.linear_solver(by, Xy) ||
       !solver.linear_solver(bz, Xz) )
    {
      std::cerr << "Could not solve linear system." << std::endl;
      return;
    }
  }

  void calculate_stiffness_matrix_elements(std::vector<Triplet>& stiffness_elements)
  {
    CGAL_assertion(stiffness_elements.empty());
    stiffness_elements.reserve(8 * nb_vert_); //estimation

    boost::unordered_map<int, NT> diag_coeff;
    for(face_descriptor f : frange_)
    {
      for(halfedge_descriptor hi : halfedges_around_face(halfedge(f, mesh_), mesh_))
      {
        halfedge_descriptor hi_opp = opposite(hi, mesh_);
        if(!is_border(hi_opp, mesh_) && hi > hi_opp) continue;
        vertex_descriptor v_source = source(hi, mesh_);
        vertex_descriptor v_target = target(hi, mesh_);

        if(!is_constrained(v_source) && !is_constrained(v_target))
        {
          int i_source = static_cast<int>(vimap_[v_source]);
          int i_target = static_cast<int>(vimap_[v_target]);
          NT Lij = weight_calculator_(hi);
          if (!is_border(hi_opp, mesh_))
            Lij+= weight_calculator_(hi_opp);
          stiffness_elements.push_back(Triplet(i_source, i_target, Lij));
          stiffness_elements.push_back(Triplet(i_target, i_source, Lij));
          diag_coeff.insert(std::make_pair(i_source,0)).first->second -= Lij;
          diag_coeff.insert(std::make_pair(i_target,0)).first->second -= Lij;
        }
      }
    }

    for(typename boost::unordered_map<int, NT>::iterator p = diag_coeff.begin();
        p != diag_coeff.end(); ++p)
    {
      stiffness_elements.push_back(Triplet(p->first, p->first, p->second));
    }
  }

  void update_mesh(Eigen_vector& Xx, Eigen_vector& Xy, Eigen_vector& Xz)
  {
    for(vertex_descriptor v : vrange_)
    {
      std::size_t index = get(vimap_, v);
      NT x_new = Xx[index];
      NT y_new = Xy[index];
      NT z_new = Xz[index];
      put(vpmap_, v, Point(x_new, y_new, z_new));
    }
  }

private:

  // compute linear system
  // -----------------------------------
  void compute_coefficient_matrix(Eigen_matrix& A,
                                  std::vector<Triplet>& stiffness_elements,
                                  const double& time)
  {
    CGAL_assertion(A.row_dimension() != 0);
    CGAL_assertion(A.column_dimension() != 0);
    CGAL_assertion(std::size_t(A.row_dimension()) == nb_vert_);
    CGAL_assertion(std::size_t(A.column_dimension()) == nb_vert_);

    calculate_D_diagonal(diagonal_);
    CGAL_assertion(diagonal_.size() == nb_vert_);

    // fill A = D - time * L
    for(const Triplet& t : stiffness_elements)
    {
      if (t.get<0>() != t.get<1>())
        A.set_coef(t.get<0>(), t.get<1>(), -time * t.get<2>(), true);
      else
      {
        A.set_coef(t.get<0>(), t.get<1>(), diagonal_[t.get<0>()] -time * t.get<2>(), true);
      }
    }

    apply_constraints(A);
    // we do not call A.assemble_matrix here
    // Eigen's compute during factorization does the building correctly,
    // and without assemble_matrix the reference A can be used in the next iterations.
  }

  void calculate_D_diagonal(std::vector<double>& diagonal)
  {
    // calculates all elements on the diagonal of D
    // and marks the entries that correspond to constrained vertices

    diagonal.assign(nb_vert_, 0);
    constrained_flags_.assign(nb_vert_, false);

    for(face_descriptor f : frange_)
    {
      double area = face_area(f, mesh_);
      for(vertex_descriptor v : vertices_around_face(halfedge(f, mesh_), mesh_))
      {
        std::size_t idx = vimap_[v];
        if(!is_constrained(v))
        {
          diagonal[idx] += area / 6.0; // 2*area/12
        }
        else
        {
          constrained_flags_[idx] = true;
          diagonal[idx] = 1.0; // needed at multiplication D*b when computing rhs
        }
      }
    }
  }

  void compute_rhs(Eigen_vector& bx, Eigen_vector& by, Eigen_vector& bz)
  {
    CGAL_assertion(diagonal_.size() == nb_vert_);
    CGAL_assertion(std::size_t(bx.size()) == nb_vert_);
    CGAL_assertion(std::size_t(by.size()) == nb_vert_);
    CGAL_assertion(std::size_t(bz.size()) == nb_vert_);

    for(vertex_descriptor vi : vrange_)
    {
      std::size_t index = vimap_[vi];
      Point_ref p = get(vpmap_, vi);
      bx.set(index, p.x());
      by.set(index, p.y());
      bz.set(index, p.z());
    }

    // multiply D * b
    for(std::size_t i = 0; i < nb_vert_; ++i)
    {
      bx[i] *= diagonal_[i];
      by[i] *= diagonal_[i];
      bz[i] *= diagonal_[i];
    }
  }

  // handling constrains
  // -----------------------------------
  void apply_constraints(Eigen_matrix& A)
  {
    CGAL_assertion(constrained_flags_.size() == nb_vert_);
    for(std::size_t i = 0; i < nb_vert_; ++i)
    {
      if(constrained_flags_[i] == true)
        A.set_coef(i, i, 1.0, true);
    }
  }

  bool is_constrained(const vertex_descriptor& v)
  {
    return get(vcmap_, v);
  }

  template<typename FaceRange>
  void set_face_range(const FaceRange& face_range)
  {
    frange_.assign(face_range.begin(), face_range.end());
    vrange_.reserve(3 * face_range.size());
    for(face_descriptor f : face_range)
    {
      for(vertex_descriptor v : vertices_around_face(halfedge(f, mesh_), mesh_))
        vrange_.push_back(v);
    }
    // get rid of duplicate vertices
    std::sort(vrange_.begin(), vrange_.end());
    vrange_.erase(std::unique(vrange_.begin(), vrange_.end()), vrange_.end());
  }

  private:
  // data members
  // ------------
  std::size_t nb_vert_;
  TriangleMesh& mesh_;
  VertexPointMap& vpmap_;
  VertexConstraintMap& vcmap_;
  IndexMap vimap_;

  // linear system data
  std::vector<double> diagonal_; // index of vector -> index of vimap_
  std::vector<bool> constrained_flags_;

  std::vector<face_descriptor> frange_;
  std::vector<vertex_descriptor> vrange_;
  Edge_cotangent_weight<TriangleMesh, VertexPointMap> weight_calculator_;
};

} // internal
} // PMP
} // CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_INTERNAL_CURVATURE_FLOW_NEW_IMPL_H
