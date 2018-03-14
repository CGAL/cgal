// Copyright (c) 2017 GeometryFactory (France).
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
// Author(s)     : Konstantinos Katrioplas (konst.katrioplas@gmail.com)


#ifndef CURVATURE_FLOW_NEW_IMPL_H
#define CURVATURE_FLOW_NEW_IMPL_H

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>
#include <CGAL/Polygon_mesh_processing/Weights.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>
#include <CGAL/barycenter.h>
#include <CGAL/utility.h>

#if defined(CGAL_EIGEN3_ENABLED)
#include <CGAL/Eigen_solver_traits.h>
#include <Eigen/Sparse>
#endif

#include <fstream>
#include <unordered_set>
#include <unordered_map>

#include <CGAL/Polygon_mesh_processing/internal/Smoothing/constraints_map.h>

namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {

template<typename PolygonMesh,
         typename VertexPointMap,
         typename VertexConstraintMap,
         typename SparseLinearSolver,
         typename GeomTraits>
class Shape_smoother_new{

private:

  typedef typename GeomTraits::FT NT;
  typedef typename GeomTraits::Point_3 Point;
  typedef typename GeomTraits::Triangle_3 Triangle;
  typedef CGAL::Triple<int, int, double> Triplet;

  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::edge_descriptor edge_descriptor;
  typedef typename boost::property_map<PolygonMesh, boost::vertex_index_t>::type IndexMap;


  // linear system
  typedef typename SparseLinearSolver::Matrix Eigen_matrix;
  typedef typename SparseLinearSolver::Vector Eigen_vector;


public:
  Shape_smoother_new(PolygonMesh& mesh,
                     VertexPointMap& vpmap,
                     VertexConstraintMap& vcmap) :
      mesh_(mesh),
      vpmap_(vpmap),
      vcmap_(vcmap),
      weight_calculator_(mesh, vpmap),
      inc_areas_calculator_(mesh),
      nb_vert_(static_cast<int>(vertices(mesh).size())) {}

  template<typename FaceRange>
  void init_smoothing(const FaceRange& face_range)
  {
    check_face_range(face_range);
  }

  void setup_system(Eigen_matrix& A,
                    Eigen_vector& bx, Eigen_vector& by, Eigen_vector& bz,
                    const double& time)
  {
    compute_coefficient_matrix(A, time);
    compute_rhs(bx, by, bz);
  }

  void solve_system(const Eigen_matrix& A,
                    Eigen_vector& Xx, Eigen_vector& Xy, Eigen_vector& Xz,
                    const Eigen_vector& bx, const Eigen_vector& by, const Eigen_vector& bz)
  {
    SparseLinearSolver solver;
    NT Dx, Dy, Dz;

    if(!solver.linear_solver(A, bx, Xx, Dx) ||
       !solver.linear_solver(A, by, Xy, Dy) ||
       !solver.linear_solver(A, bz, Xz, Dz) )
    {
      std::cerr << "Could not solve linear system." << std::endl;
      return;
    }

  }

  void calculate_stiffness_matrix_elements()
  {
    CGAL_assertion(stiffness_elements_.empty());
    stiffness_elements_.reserve(8 * nb_vert_); //estimation


    BOOST_FOREACH(face_descriptor f, frange_)
    {
      BOOST_FOREACH(halfedge_descriptor hi, halfedges_around_face(halfedge(f, mesh_), mesh_))
      {
        vertex_descriptor v_source = source(hi, mesh_);
        vertex_descriptor v_target = target(hi, mesh_);

        if(!is_constrained(v_source) && !is_constrained(v_target))
        {
          int i_source = vimap_[v_source];
          int i_target = vimap_[v_target];
          NT Lij = weight_calculator_(hi);
          stiffness_elements_.push_back(Triplet(i_source, i_target, Lij));
          stiffness_elements_.push_back(Triplet(i_target, i_source, Lij));
          stiffness_elements_.push_back(Triplet(i_source, i_source, -Lij));
          stiffness_elements_.push_back(Triplet(i_target, i_target, -Lij));
        }
      }
    }


    /*

    boost::unordered_map<std::size_t, NT> diag_coeff;
    BOOST_FOREACH(face_descriptor f, frange_)
    {
      BOOST_FOREACH(halfedge_descriptor hi, halfedges_around_face(halfedge(f, mesh_), mesh_))
      {
        halfedge_descriptor hi_opp = opposite(hi, mesh_);
        if(!is_border(hi_opp, mesh_) && hi > hi_opp) continue;
        vertex_descriptor v_source = source(hi, mesh_);
        vertex_descriptor v_target = target(hi, mesh_);

        if(!is_constrained(v_source) && !is_constrained(v_target))
        {
          auto i_source = vimap_[v_source];
          auto i_target = vimap_[v_target];
          NT Lij = weight_calculator_(hi);
          if (!is_border(hi_opp, mesh_))
            Lij+= weight_calculator_(hi_opp);
          stiffness_elements_.push_back(Triplet(i_source, i_target, Lij));
          stiffness_elements_.push_back(Triplet(i_target, i_source, Lij));
          diag_coeff.insert(std::make_pair(i_source,0)).first->second -= Lij;
          diag_coeff.insert(std::make_pair(i_target,0)).first->second -= Lij;
        }
      }
    }
    for (auto p : diag_coeff)
      stiffness_elements_.push_back(Triplet(p.first, p.first, p.second));

    */

  }

  void update_mesh(Eigen_vector& Xx, Eigen_vector& Xy, Eigen_vector& Xz)
  {
    update_map(Xx, Xy, Xz); // TODO: REMOVE the function
  }

private:

  void compute_coefficient_matrix(Eigen_matrix& A, const double& time)
  {
    CGAL_assertion(A.row_dimension() != 0);
    CGAL_assertion(A.column_dimension() != 0);
    CGAL_assertion(A.row_dimension() == nb_vert_);
    CGAL_assertion(A.column_dimension() == nb_vert_);

    calculate_D_diagonal(diagonal_);
    CGAL_assertion(diagonal_.size() == nb_vert_);


    // fill A = D - time * L



    for(int idx = 0; idx < diagonal_.size(); ++idx)
    {
      A.set_coef(idx, idx, diagonal_[idx], true); // A is empty
    }

    for(Triplet t : stiffness_elements_)
    {
      A.add_coef(t.get<0>(), t.get<1>(), -time * t.get<2>());
    }


    /*

    std::vector<bool> diag_set(diagonal_.size(), false);

    for(Triplet t : stiffness_elements_)
    {
      if (t.get<0>() != t.get<1>())
        A.set_coef(t.get<0>(), t.get<1>(), -time * t.get<2>(), true);
      else
      {
        A.set_coef(t.get<0>(), t.get<1>(), diagonal_[t.get<0>()] -time * t.get<2>(), true);
        diag_set[t.get<0>()]=true;
      }
    }

    for(int idx = 0; idx < diagonal_.size(); ++idx)
    {
      if (!diag_set[idx])
      {
        std::cout << "diag_NOT_set\n"  ;
        A.set_coef(idx, idx, diagonal_[idx], true); // A is empty
      }
    }
    */



    A.assemble_matrix(); // does setFromTriplets and some compression



    /*
    std::ofstream out("data/traits_A.dat");
    for(auto j = 0 ; j < A.column_dimension(); ++j)
    {
      for(auto i = 0; i < A.row_dimension(); ++i)
      {

        out << A.get_coef(i, j) << " ";
      }
      out << std::endl;
    }
    */



  }

  void calculate_D_diagonal(std::vector<double>& diagonal)
  {
    diagonal.clear();
    //diagonal.resize(nb_vert_);
    diagonal.assign(nb_vert_, 0);
    for(face_descriptor f : frange_)
    {
      double area = face_area(f, mesh_);
      for(vertex_descriptor v : vertices_around_face(halfedge(f, mesh_), mesh_))
      {
        int idx = vimap_[v];
        if(!is_constrained(v))
        {
          diagonal[idx] += 2.0 * area;
        }
        else
        {
          diagonal[idx] = 1.0;
        }
      }
    }

    for(int i = 0; i < diagonal.size(); ++i)
    {
      diagonal[i] /= 12.0; // constraints (= 1) as well ?
    }

  }

  void compute_rhs(Eigen_vector& bx, Eigen_vector& by, Eigen_vector& bz)
  {
    CGAL_assertion(diagonal_.size() == nb_vert_);
    CGAL_assertion(bx.size() == nb_vert_);
    CGAL_assertion(by.size() == nb_vert_);
    CGAL_assertion(bz.size() == nb_vert_);


    for(vertex_descriptor vi : vrange_)
    {
      int index = vimap_[vi];
      Point p = get(vpmap_, vi);
      bx.set(index, p.x());
      by.set(index, p.y());
      bz.set(index, p.z());
    }

    for(int i = 0; i < nb_vert_; ++i)
    {
      bx[i] *= diagonal_[i];
      by[i] *= diagonal_[i];
      bz[i] *= diagonal_[i];
    }


    /*

    std::ofstream outbx("data/outb-traits-x");
    for(int i = 0 ; i < nb_vert_; ++i)
    {
      outbx << bx[i] << std::endl;
    }
    std::ofstream outby("data/outb-traits-y");
    for(int i = 0 ; i < by.dimension(); ++i)
    {
      outby << by[i] << std::endl;
    }
    std::ofstream outbz("data/outb-traits-z");
    for(int i = 0 ; i < bz.dimension(); ++i)
    {
      outbz << bz[i] << std::endl;
    }
    */



  }

  // update mesh
  // -----------------------------------
  Triangle triangle(face_descriptor f) const
  {
    halfedge_descriptor h = halfedge(f, mesh_);
    vertex_descriptor v1  = target(h, mesh_);
    vertex_descriptor v2  = target(next(h, mesh_), mesh_);
    vertex_descriptor v3  = target(next(next(h, mesh_), mesh_), mesh_);
    return Triangle(get(vpmap_, v1), get(vpmap_, v2), get(vpmap_, v3));
  }

  void update_map(Eigen_vector& Xx, Eigen_vector& Xy, Eigen_vector& Xz)
  {
    for (vertex_descriptor v : vertices(mesh_))
    {
      int index = get(vimap_, v);
      NT x_new = Xx[index];
      NT y_new = Xy[index];
      NT z_new = Xz[index];
      put(vpmap_, v, Point(x_new, y_new, z_new));
    }
  }

  // handling constrains
  // -----------------------------------
  bool is_constrained(const vertex_descriptor& v)
  {
    return get(vcmap_, v);
  }

  template<typename FaceRange>
  void check_face_range(const FaceRange& face_range)
  {
    BOOST_FOREACH(face_descriptor f, face_range)
    {
      frange_.insert(f);
      BOOST_FOREACH(vertex_descriptor v, vertices_around_face(halfedge(f, mesh_), mesh_))
        vrange_.insert(v);
    }
  }

  private:
  // data members
  // ------------
  PolygonMesh& mesh_;
  VertexPointMap& vpmap_;
  std::size_t nb_vert_;

  // linear system data
  std::vector<Triplet> stiffness_elements_;
  std::vector<double> diagonal_; // index of vector -> index of vimap_

  std::set<face_descriptor> frange_;
  std::set<vertex_descriptor> vrange_;
  IndexMap vimap_ = get(boost::vertex_index, mesh_);
  VertexConstraintMap vcmap_;
  Edge_cotangent_weight<PolygonMesh, VertexPointMap> weight_calculator_;
  Incident_area<PolygonMesh> inc_areas_calculator_;

};


} // internal
} // PMP
} // CGAL















#endif // CURVATURE_FLOW_NEW_IMPL_H
