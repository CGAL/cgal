// Copyright (c) 2015 GeometryFactory (France).
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
// Author(s)     : Jane Tournois

#ifndef CGAL_POLYGON_MESH_PROCESSING_FAIR_POLYHEDRON_3_H
#define CGAL_POLYGON_MESH_PROCESSING_FAIR_POLYHEDRON_3_H

#include <CGAL/license/Polygon_mesh_processing/meshing_hole_filling.h>


#include <map>
#include <set>
#include <CGAL/assertions.h>
#include <CGAL/Polygon_mesh_processing/Weights.h>
#ifdef CGAL_PMP_FAIR_DEBUG
#include <CGAL/Timer.h>
#endif
#include <iterator>

namespace CGAL {

namespace Polygon_mesh_processing {

namespace internal {

// [On Linear Variational Surface Deformation Methods-2008]
template<class PolygonMesh,
         class SparseLinearSolver,
         class WeightCalculator,
         class VertexPointMap>
class Fair_Polyhedron_3 {
  // typedefs
  typedef typename VertexPointMap::value_type Point_3;
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
  typedef  Halfedge_around_target_circulator<PolygonMesh>  Halfedge_around_vertex_circulator;

  typedef SparseLinearSolver Sparse_linear_solver;
  typedef typename Sparse_linear_solver::Matrix Solver_matrix;
  typedef typename Sparse_linear_solver::Vector Solver_vector;

// members
  PolygonMesh& pmesh;
  WeightCalculator weight_calculator;
  VertexPointMap ppmap;

public:
  Fair_Polyhedron_3(PolygonMesh& pmesh
      , VertexPointMap vpmap
      , WeightCalculator weight_calculator)
    : pmesh(pmesh)
    , weight_calculator(weight_calculator)
    , ppmap(vpmap)
  { }
  
private:
  double sum_weight(vertex_descriptor v) {
  double weight = 0;
  Halfedge_around_vertex_circulator circ(halfedge(v,pmesh),pmesh), done(circ);
  do {
    weight += weight_calculator.w_ij(*circ);
    } while(++circ != done);
    return weight;
  }

  // recursively computes a row (use depth parameter to compute L, L^2, L^3)
  // Equation 6 in [On Linear Variational Surface Deformation Methods]
  void compute_row(
    vertex_descriptor v,
    int row_id,                            // which row to insert in [ frees stay left-hand side ]
    Solver_matrix& matrix, 
    double& x, double& y, double& z,               // constants transfered to right-hand side
    double multiplier,
    const std::map<vertex_descriptor, std::size_t>& vertex_id_map,
    unsigned int depth)
  {
    if(depth == 0) {
      typename std::map<vertex_descriptor, std::size_t>::const_iterator vertex_id_it = vertex_id_map.find(v);
      if(vertex_id_it != vertex_id_map.end()) {
        int col_id = static_cast<int>(vertex_id_it->second);
        matrix.add_coef(row_id, col_id, multiplier);
      }
      else { 
        typename boost::property_traits<VertexPointMap>::reference p = get(ppmap, v);
        x += multiplier * - to_double(p.x());
        y += multiplier * - to_double(p.y());
        z += multiplier * - to_double(p.z());
      }
    }
    else {
      double w_i = weight_calculator.w_i(v);

      Halfedge_around_vertex_circulator circ(halfedge(v,pmesh),pmesh), done(circ);
      do {
        double w_i_w_ij = w_i * weight_calculator.w_ij(*circ) ;

        vertex_descriptor nv = target(opposite(*circ,pmesh),pmesh);
        compute_row(nv, row_id, matrix, x, y, z, -w_i_w_ij*multiplier, vertex_id_map, depth-1);
      } while(++circ != done);

      double w_i_w_ij_sum = w_i * sum_weight(v);
      compute_row(v, row_id, matrix, x, y, z, w_i_w_ij_sum*multiplier, vertex_id_map, depth-1);
    }
  }

public:
  template<class VertexRange>
  bool fair(const VertexRange& vertices
    , SparseLinearSolver solver
    , unsigned int fc)
  {
    int depth = static_cast<int>(fc) + 1;
    if(depth < 0 || depth > 3) {
      CGAL_warning(!"Continuity should be between 0 and 2 inclusively!");
      return false; 
    }

    std::set<vertex_descriptor> interior_vertices(boost::begin(vertices),
                                                  boost::end(vertices));
    if(interior_vertices.empty()) { return true; }

    #ifdef CGAL_PMP_FAIR_DEBUG
    CGAL::Timer timer; timer.start();
    #endif
    const std::size_t nb_vertices = interior_vertices.size();
    Solver_vector X(nb_vertices), Bx(nb_vertices);
    Solver_vector Y(nb_vertices), By(nb_vertices);
    Solver_vector Z(nb_vertices), Bz(nb_vertices);

    std::map<vertex_descriptor, std::size_t> vertex_id_map;
    std::size_t id = 0;
    BOOST_FOREACH(vertex_descriptor vd, interior_vertices)
    {
      if( !vertex_id_map.insert(std::make_pair(vd, id)).second ) {
        CGAL_warning(!"Duplicate vertex is found!");
        return false;
      }
      ++id;
    }

    Solver_matrix A(nb_vertices);

    BOOST_FOREACH(vertex_descriptor vd, interior_vertices)
    {
      int v_id = static_cast<int>(vertex_id_map[vd]);
      compute_row(vd, v_id, A, Bx[v_id], By[v_id], Bz[v_id], 1, vertex_id_map, depth);
    }
    #ifdef CGAL_PMP_FAIR_DEBUG
    std:cerr << "**Timer** System construction: " << timer.time() << std::endl; timer.reset();
    #endif

    // factorize
    double D;
    bool prefactor_ok = solver.factor(A, D);
    if(!prefactor_ok) {
      CGAL_warning(!"pre_factor failed!");
      return false;
    }
    #ifdef CGAL_PMP_FAIR_DEBUG
    std::cerr << "**Timer** System factorization: " << timer.time() << std::endl; timer.reset();
    #endif

    // solve
    bool is_all_solved = solver.linear_solver(Bx, X) && solver.linear_solver(By, Y) && solver.linear_solver(Bz, Z);
    if(!is_all_solved) {
      CGAL_warning(!"linear_solver failed!"); 
      return false; 
    }
    #ifdef CGAL_PMP_FAIR_DEBUG
    std::cerr << "**Timer** System solver: " << timer.time() << std::endl; timer.reset();
    #endif

    
    /* This relative error is to large for cases that the results are not good */ 
    /*
    double rel_err_x = (A.eigen_object()*X - Bx).norm() / Bx.norm();
    double rel_err_y = (A.eigen_object()*Y - By).norm() / By.norm();
    double rel_err_z = (A.eigen_object()*Z - Bz).norm() / Bz.norm();
    CGAL_TRACE_STREAM << "rel error: " << rel_err_x 
                                << " " << rel_err_y
                                << " " << rel_err_z << std::endl;
                                */

    // update 
    id = 0;
    BOOST_FOREACH(vertex_descriptor vd, interior_vertices)
    {
      put(ppmap, vd, Point_3(X[id], Y[id], Z[id]));
      ++id;
    }
    return true;
  }
};

}//namespace internal

}//namespace Polygon_mesh_processing

}//namespace CGAL
#endif //CGAL_POLYGON_MESH_PROCESSING_FAIR_POLYHEDRON_3_H
