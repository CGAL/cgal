// Copyright (c) 2016  GeometryFactory (France).
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
// Author(s)     :

#ifndef CGAL_SURFACE_MESH_PARAMETERIZATION_ORBITAL_TUTTE_PARAMETERIZER_3_H
#define CGAL_SURFACE_MESH_PARAMETERIZATION_ORBITAL_TUTTE_PARAMETERIZER_3_H

#include <CGAL/Surface_mesh_parameterization/internal/angles.h>
#include <CGAL/Surface_mesh_parameterization/internal/kernel_traits.h>
#include <CGAL/Surface_mesh_parameterization/internal/orbital_cone_helper.h>
#include <CGAL/Surface_mesh_parameterization/IO/File_off.h>

#include <CGAL/Surface_mesh_parameterization/Error_code.h>

#include <CGAL/Polygon_mesh_processing/Weights.h>

#include <CGAL/assertions.h>
#include <CGAL/circulator.h>
#include <CGAL/Eigen_solver_traits.h>
#include <CGAL/Timer.h>

#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#ifdef CGAL_SMP_USE_SPARSESUITE_SOLVERS
#include <Eigen/SPQRSupport>
#include <Eigen/UmfPackSupport>
#endif

#include <boost/array.hpp>
#include <boost/foreach.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <utility>

/// \file Orbital_Tutte_parameterizer_3.h

// #define CGAL_SMP_OUTPUT_ORBITAL_MATRICES

// @todo checks that cones are different, are on seams, seam is one connected
//       component
// @todo Should the order of cones provided in entry matter ? Map the first cone
//       to [-1, -1] for example ?

namespace CGAL {

namespace Surface_mesh_parameterization {

enum Weight_type
{
  Cotangent = 0,
  Mean_value
};

template
<
  typename SeamMesh,
  typename SparseLinearAlgebraTraits_d
#ifdef CGAL_SMP_USE_SPARSESUITE_SOLVERS
    = Eigen_solver_traits<Eigen::UmfPackLU<Eigen_sparse_matrix<double>::EigenType> >
#else
    = Eigen_solver_traits<Eigen::SparseLU<Eigen_sparse_matrix<double>::EigenType> >
#endif
>
class Orbital_Tutte_parameterizer_3
{
private:
  typedef typename boost::graph_traits<SeamMesh>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::graph_traits<SeamMesh>::halfedge_descriptor  halfedge_descriptor;
  typedef typename boost::graph_traits<SeamMesh>::face_descriptor      face_descriptor;

  typedef typename boost::graph_traits<SeamMesh>::vertex_iterator      vertex_iterator;
  typedef typename boost::graph_traits<SeamMesh>::face_iterator        face_iterator;

  // SparseLinearAlgebraTraits_d subtypes:
  typedef SparseLinearAlgebraTraits_d                               Sparse_LA;
  typedef typename Sparse_LA::Vector                                Vector;
  typedef typename Sparse_LA::Matrix                                Matrix;

  // Kernel subtypes
  typedef typename internal::Kernel_traits<SeamMesh>::Kernel    Kernel;
  typedef typename internal::Kernel_traits<SeamMesh>::PPM       PPM;
  typedef typename Kernel::FT                                       NT;
  typedef typename Kernel::Vector_2                                 Vector_2;
  typedef typename Kernel::Vector_3                                 Vector_3;
  typedef typename Kernel::Point_2                                  Point_2;
  typedef typename Kernel::Point_3                                  Point_3;

  Orbifold_type orb_type;
  Weight_type weight_type;

private:
  /// check input's correctness.
  template<typename ConeMap>
  Error_code check_cones(ConeMap cmap) const
  {
    if(orb_type == Parallelogram) {
      if(cmap.size() != 6) {
        std::cerr << "Using orb_type '" << get_orbifold_type(orb_type)
                  << "' requires 4 vertices marked as cones (thus 6 in the seam mesh)" << std::endl;
        std::cerr << "currently: " << cmap.size() << std::endl;
        return ERROR_WRONG_PARAMETER;
      }
    } else if(cmap.size() != 4){ // orb_type == Square, Diamond, Triangle
      std::cerr << "Using orb_type '" << get_orbifold_type(orb_type)
                << "' requires 3 vertices marked as cones (thus 4 in the seam mesh)" << std::endl;
      std::cerr << "currently: " << cmap.size() << std::endl;
      return ERROR_WRONG_PARAMETER;
    }

    std::cout << "cones and ids" << std::endl;
    typename ConeMap::const_iterator cit = cmap.begin(), cend = cmap.end();
    for(; cit!=cend; ++cit) {
//      std::cout << target(halfedge(cit->first,mesh), mesh.mesh())
//                << " n°: " << get(vimap, cit->first) << std::endl;
    }

    return OK;
  }

  // Linear system
  /// Compute the number of linear constraints in the system.
  int number_of_linear_constraints(const SeamMesh& mesh) const
  {
    if(orb_type == Parallelogram) {
      // number of constraints for orb I, II, III is the number of seam edges
      // and 3 constrained cones.
      return 3 + static_cast<int>(mesh.number_of_seam_edges());
    }
    else { // orb_type == Square, Diamond, Triangle
      // number of constraints for orb I, II, III is the number of seam edges
      // and 2 constrained cones.
      return 2 + static_cast<int>(mesh.number_of_seam_edges());
    }
  }

  /// Adds a positional constraint on a vertex x_ind, so that x_ind * w = rhs.
  void addConstraint(Matrix& A, Vector& B, int& id_r, int id_c, double w, Point_2 rhs) const
  {
    std::cout << "Constraining " << id_c << " to " << rhs << std::endl;

    A.set_coef(2*id_r, 2*id_c, w, true /*new_coef*/);
    A.set_coef(2*id_r + 1, 2*id_c + 1, w, true /*new_coef*/);

    B[2*id_r] = rhs.x();
    B[2*id_r + 1] = rhs.y();

    ++id_r; // current line index in A is increased
  }

  /// Adds constraints so that T * x_sinds = x_tinds, where T is a 2x2
  /// matrix, and the Transformation T is modified to affine from
  /// linear by requiring that T * x_si - x_ti = T * x_s1 - x_t1.
  void addTransConstraints(int s0, int t0, int s, int t,
                           int& id_r,
                           const Eigen::Matrix2d& T,
                           Matrix& A, Vector& B) const
  {
    std::cout << "transconstraints: " << s << " " << t << std::endl;

    // Matlab lines are commented for comparison.
    // Matlab fills together 2*x-1 and 2*x, but C++ fills 2*x and 2*x+1,
    // as everything (including loops!) starts at 0 and not 1.

    // iterate on both rows ot the 2x2 matrix T
    for(int vert_ind=0; vert_ind<2; ++vert_ind) {
      // building up the equations by summing up the terms

      // <T(vert_ind,:), x_si>
      // obj.A(end+1, 2*sinds(ind)+[-1,0]) = T(vert_ind,:);
      A.set_coef(2*id_r + vert_ind, 2*s, T(vert_ind, 0), true /*new_coef*/);
      A.set_coef(2*id_r + vert_ind, 2*s + 1, T(vert_ind, 1), true /*new_coef*/);

      // -<T(vert_ind,:), x_s1>
      // obj.A(end, 2*sinds(1)+[-1,0]) = obj.A(end, 2*sinds(1)+[-1,0]) - T(vert_ind,:);
      A.add_coef(2*id_r + vert_ind, 2*s0, -T(vert_ind, 0));
      A.add_coef(2*id_r + vert_ind, 2*s0 + 1, -T(vert_ind, 1));

      //  - x_ti
      // obj.A(end, 2*tinds(ind)+vert_ind-2) = obj.A(end, 2*tinds(ind)+vert_ind-2)-1;
      A.add_coef(2*id_r + vert_ind, 2*t + vert_ind, -1);

      // + x_t1
      // obj.A(end, 2*tinds(1)+vert_ind-2) = obj.A(end, 2*tinds(1)+vert_ind-2)+1;
      A.add_coef(2*id_r + vert_ind, 2*t0 + vert_ind, 1);

      // left hand side is zero
      // obj.b=[obj.b; 0];
      B[2*id_r + vert_ind] = 0;
    }

    ++id_r; // current line index in A is increased
  }

  void constrain_seam_segment(const std::vector<std::pair<int, int> >& seam_segment,
                              NT ang, int& current_line_id_in_A,
                              Matrix& A, Vector& B) const
  {
    std::cout << "constraining segment of length " << seam_segment.size() << std::endl;

    // check that if there is a common vertex, it is at the beginning
    bool is_reversed = (seam_segment.back().first == seam_segment.back().second);

    if(is_reversed) {
      ang *= -1;
    }

    // the rotation matrix according to the angle 'ang'
    Eigen::Matrix2d R;
    R(0,0) = std::cos(2 * CGAL_PI / ang); R(0,1) = - std::sin(2 * CGAL_PI / ang);
    R(1,0) = std::sin(2 * CGAL_PI / ang); R(1,1) = std::cos(2 * CGAL_PI / ang);

    int s0 = is_reversed ? seam_segment.back().first : seam_segment.front().first;
    int t0 = is_reversed ? seam_segment.back().second : seam_segment.front().second;

    std::cout << "s0/t0: " << s0 << " " << t0 << std::endl;

    typename std::vector<std::pair<int, int> >::const_iterator it = seam_segment.begin(),
                                                               end = seam_segment.end();

    // ignore the first entry of the seam segment (they correspond to a constrained point)
    if(is_reversed)
      --end;
    else
      ++it;

    for(; it!=end; ++it) {
      int s = it->first;
      int t = it->second;
      std::cout << "v1/v2: " << s << " " << t << std::endl;
      CGAL_assertion(s != t);

      // sending s to t (and _not_ t to s !)
      addTransConstraints(t0, s0, t, s, current_line_id_in_A, R, A, B);
    }
  }

  /// Computes the rotational constraint on the border of the mesh.
  /// Cone constraints are also added.
  template<typename ConeMap,
           typename VertexIndexMap>
  void AddRotationalConstraint(const SeamMesh& mesh,
                               const ConeMap& cmap,
                               VertexIndexMap vimap,
                               Matrix& A, Vector& B) const
  {
    // positions of the cones in the plane
    typedef std::vector<Point_2>                      Point_container;
    const Point_container& tcoords =
                  get_cones_parameterized_coordinates<Point_container>(orb_type);

    // angles at the cones
    typedef std::vector<NT>                           Angle_container;
    const Angle_container& angs = get_angles_at_cones<Angle_container>(orb_type);

    // the index of the line in A that we are filling next
    int current_line_id_in_A = 0.;

    // Initialize some variables used in the seam walk
    int start_cone_index = -1; // index of the beginning of the seam
    vertex_descriptor start_cone;
    internal::find_start_cone(cmap, vimap, start_cone, start_cone_index);
    CGAL_postcondition(start_cone != vertex_descriptor() && start_cone_index != -1);

//    std::cout << "initial cone is " << start_cone << std::endl;

    // parameterize the initial cone
    addConstraint(A, B, current_line_id_in_A, start_cone_index,
                  1. /*entry in A*/, tcoords[0]);

    // by property of the seam mesh, the canonical halfedge that points to start_cone
    // is on the seam, and is not on the border
    halfedge_descriptor hd = halfedge(start_cone, mesh);
    CGAL_precondition(mesh.has_on_seam(hd));
    halfedge_descriptor bhd = opposite(hd, mesh);
    CGAL_precondition(is_border(bhd, mesh));

    // points between two cones, and the corresponding points on the opposite side of the seam
    std::vector<std::pair<int, int> > seam_segment;
    std::size_t segment_index = 0; // counting the segments (3 max)

    // Go through the seam, marking rotation and cone constraints
    while(true) { // breaking at the last cone
      // Get the two halfedges on each side of the stream
      halfedge_descriptor hd1 = bhd; // only for clarity

      // the non-border halfedge with same vertices (in the underlying mesh of the seam
      // mesh) as bhd is simply bhd with the 'seam' boolean set to false
      halfedge_descriptor hd2 = bhd;
      hd2.seam = false;

      // Compute the corresponding indices
      vertex_descriptor hd1_source = source(hd1, mesh);
      vertex_descriptor hd2_source = source(hd2, mesh);
      int hd1s_index = get(vimap, hd1_source);
      int hd2s_index = get(vimap, hd2_source);
      std::cout << hd1s_index << " " << hd2s_index << std::endl;

      // If orbifold type IV and it is second cone in flattening, add constraint
      if(orb_type == Parallelogram && cmap.find(hd1_source) != cmap.end()
                                   && segment_index == 1) {
        addConstraint(A, B, current_line_id_in_A, hd1s_index,
                      1. /*entry in A*/, tcoords[1]);
      }

      // Add the pair to the seam segment
      seam_segment.push_back(std::make_pair(hd1s_index, hd2s_index));

      // Check if we have reached a cone
      vertex_descriptor bhd_target = target(bhd, mesh);
      typename ConeMap::const_iterator is_in_map = cmap.find(bhd_target);
      if(is_in_map != cmap.end()) {
        // add the target to finish the seam segment
        vertex_descriptor hd1_target = target(hd1, mesh);
        vertex_descriptor hd2_target = target(hd2, mesh);
        int hd1t_index = get(vimap, hd1_target);
        int hd2t_index = get(vimap, hd2_target);

        seam_segment.push_back(std::make_pair(hd1t_index, hd2t_index));

        CGAL_assertion(segment_index < angs.size());
        NT ang = angs[segment_index];
        constrain_seam_segment(seam_segment, ang, current_line_id_in_A, A, B);

        // Check if we have reached the end of the seam
        if(is_in_map->second == Second_unique_cone) {
          CGAL_assertion(hd1_target == hd2_target);
        // the last cone of the seam is constrained
          addConstraint(A, B, current_line_id_in_A, hd1t_index,
                        1. /*entry in A*/, tcoords.back());
          break;
        }

        std::cout << "-------------------------" << std::endl;
        seam_segment.clear();
        ++segment_index;
      }

      // move to the next halfedge couple (walking on the border of the seam)
      bhd = next(bhd, mesh);
      CGAL_postcondition(mesh.has_on_seam(bhd) && is_border(bhd, mesh));
    }

    std::cout << current_line_id_in_A << " vs " << number_of_linear_constraints(mesh) << std::endl;
    CGAL_postcondition(current_line_id_in_A == number_of_linear_constraints(mesh));
  }

  // MVC computations
  NT compute_w_ij_mvc(const Point_3& pi, const Point_3& pj, const Point_3& pk) const
  {
    //                                                               ->     ->
    // Compute the angle (pj, pi, pk), the angle between the vectors ij and ik
    NT angle = internal::compute_angle_rad<Kernel>(pj, pi, pk);
    NT weight = std::tan(0.5 * angle);

    return weight;
  }

  /// Computes the coefficients of the mean value Laplacian matrix for the edge.
  /// `ij` in the face `ijk`
  void fill_mvc_matrix(const Point_3& pi, int i,
                       const Point_3& pj, int j,
                       const Point_3& pk, int k, Matrix& L) const
  {
//    std::cout << "compute mvc coef at: " << i << " " << j << " " << k << std::endl;

    // For MVC, the entry of L(i,j) is - [ tan(gamma_ij/2) + tan(delta_ij)/2 ] / |ij|
    // where gamma_ij and delta_ij are the angles at i around the edge ij

    // This function computes the angle alpha at i, and add
    // -- A(i,j) += tan(alpha / 2) / |ij|
    // -- A(i,k) += tan(alpha / 2) / |ik|
    // -- A(i,i) -= A(i,j) + A(i,k)

    // The other parts of A(i,j) and A(i,k) will be added when this function
    // is called from the neighboring faces of F_ijk that share the vertex i

    // Compute: - tan(alpha / 2)
    NT w_i_base = 1.0 * compute_w_ij_mvc(pi, pj, pk);

    // @fixme unefficient: lengths are computed (and inversed!) twice per edge

    // Set w_ij in matrix
    Vector_3 edge_ij = pi - pj;
    NT len_ij = CGAL::sqrt(edge_ij * edge_ij);
    CGAL_assertion(len_ij != 0.0); // two points are identical!
    NT w_ij = w_i_base / len_ij;
    L.add_coef(2*i, 2*j, w_ij);
    L.add_coef(2*i +1, 2*j + 1, w_ij);

    // Set w_ik in matrix
    Vector_3 edge_ik = pi - pk;
    NT len_ik = CGAL::sqrt(edge_ik * edge_ik);
    CGAL_assertion(len_ik != 0.0); // two points are identical!
    NT w_ik = w_i_base / len_ik;
    L.add_coef(2*i, 2*k, w_ik);
    L.add_coef(2*i + 1, 2*k + 1, w_ik);

    // Add to w_ii (w_ii = - sum w_ij)
    NT w_ii = - w_ij - w_ik;
    L.add_coef(2*i, 2*i, w_ii);
    L.add_coef(2*i + 1, 2*i + 1, w_ii);
  }

  /// Compute the mean value Laplacian matrix.
  template<typename VertexIndexMap>
  void mean_value_laplacian(const SeamMesh& mesh,
                            VertexIndexMap vimap,
                            Matrix& L) const
  {
    const PPM ppmap = get(vertex_point, mesh);

    BOOST_FOREACH(face_descriptor fd, faces(mesh)) {
      halfedge_descriptor hd = halfedge(fd, mesh);

      vertex_descriptor vd_i = target(hd, mesh);
      vertex_descriptor vd_j = target(next(hd, mesh), mesh);
      vertex_descriptor vd_k = source(hd, mesh);
      const Point_3& pi = get(ppmap, vd_i);
      const Point_3& pj = get(ppmap, vd_j);
      const Point_3& pk = get(ppmap, vd_k);
      int i = get(vimap, vd_i);
      int j = get(vimap, vd_j);
      int k = get(vimap, vd_k);

      fill_mvc_matrix(pi, i, pj, j, pk, k, L);
      fill_mvc_matrix(pj, j, pk, k, pi, i, L);
      fill_mvc_matrix(pk, k, pi, i, pj, j, L);
    }
  }

  template<typename VertexIndexMap>
  void cotangent_laplacien(SeamMesh& mesh,
                           VertexIndexMap vimap,
                           Matrix& L) const
  {
    const PPM ppmap = get(vertex_point, mesh);

    // not exactly sure which cotan weights should be used:
    // 0.5 (cot a + cot b) ? 1/T1 cot a + 1/T2 cot b ? 1/Vor(i) (cot a + cot b?)
    // Comparing to the matlab code, the basic Cotangent_weight gives the same results.
    typedef CGAL::internal::Cotangent_weight<SeamMesh>                      Cotan_weights;
//    typedef CGAL::internal::Cotangent_weight_with_triangle_area<SeamMesh>   Cotan_weights;

    Cotan_weights cotan_weight_calculator(mesh, ppmap);

    BOOST_FOREACH(halfedge_descriptor hd, halfedges(mesh)) {
      vertex_descriptor vi = source(hd, mesh);
      vertex_descriptor vj = target(hd, mesh);
      int i = get(vimap, vi);
      int j = get(vimap, vj);

      if(i > j)
        continue;

      // times 2 because Cotangent_weight returns 1/2 (cot alpha + cot beta)...
      double w_ij = 2 * cotan_weight_calculator(hd);

      // ij
      L.set_coef(2*i, 2*j, w_ij, true /* new coef */);
      L.set_coef(2*i +1, 2*j + 1, w_ij, true /* new coef */);

      // ji
      L.set_coef(2*j, 2*i, w_ij, true /* new coef */);
      L.set_coef(2*j +1, 2*i + 1, w_ij, true /* new coef */);

      // ii
      L.add_coef(2*i, 2*i, - w_ij);
      L.add_coef(2*i + 1, 2*i + 1, - w_ij);

      // jj
      L.add_coef(2*j, 2*j, - w_ij);
      L.add_coef(2*j + 1, 2*j + 1, - w_ij);
    }
  }

  /// Copy the solution into the UV property map.
  template <typename VertexIndexMap, typename VertexUVMap>
  void assign_solution(const SeamMesh& mesh,
                       const Vector& X,
                       VertexUVMap uvmap,
                       const VertexIndexMap vimap) const
  {
    std::cout << "size of X: " << X.size() << std::endl;
    CGAL_assertion(X.size() == static_cast<int>(2 * num_vertices(mesh)));

    BOOST_FOREACH(vertex_descriptor vd, vertices(mesh)) {
      int index = get(vimap, vd);
      NT u = X(2*index);
      NT v = X(2*index + 1);

      put(uvmap, vd, Point_2(u, v));
    }
  }

  /// Solves the linear system.
  template <typename VertexUVMap,
            typename VertexIndexMap>
  Error_code computeFlattening(const SeamMesh& mesh,
                               const Matrix& A, const Vector& B,
                               const Matrix& L,
                               VertexUVMap uvmap, VertexIndexMap vimap) const
  {
    CGAL_precondition(A.column_dimension() == L.column_dimension());
    std::size_t n = L.column_dimension();
    std::size_t big_n = n + A.row_dimension();

    std::cout << "n: " << n << " and big_n: " << big_n << std::endl;

    // Fill the large system M.Xf = Bf:
    // ( L A' ) ( Xf ) = ( B )
    // ( A 0  ) ( Xf ) = ( 0 )
    Matrix M(big_n, big_n);
    Vector Bf(big_n);
    NT D;
    Vector Xf(big_n);

    // full B
    for(int i=0; i<B.size(); ++i) {
      Bf[n+i] = B[i];
    }

    std::cout << "filled Bf" << std::endl;

    // matrix M
    for(int k=0; k<L.eigen_object().outerSize(); ++k) {
      for(typename Eigen::SparseMatrix<double>::InnerIterator
                                            it(L.eigen_object(), k); it; ++it) {
        M.set_coef(it.row(), it.col(), it.value(), true /*new_coef*/);
//        std::cout <<  it.row() << " " << it.col() << " " << it.value() << '\n';
      }
    }

    std::cout << "filled topleft of M" << std::endl;

    for(int k=0; k<A.eigen_object().outerSize(); ++k) {
      for(typename Eigen::SparseMatrix<double>::InnerIterator
                                            it(A.eigen_object(), k); it; ++it) {
        M.set_coef(it.col(), it.row() + n, it.value(), true /*new_coef*/); // A
        M.set_coef(it.row() + n, it.col(), it.value(), true /*new_coef*/); // A'
//        std::cout <<  it.row() << " " << it.col() << " " << it.value() << '\n';
      }
    }

    std::cout << "Filled M and Bf" << std::endl;

#ifdef CGAL_SMP_OUTPUT_ORBITAL_MATRICES
    std::ofstream outM("matrices/M.txt");
    outM.precision(20);
    outM << M.row_dimension() << " " << M.column_dimension() << std::endl;
    for(int k=0; k<M.eigen_object().outerSize(); ++k) {
      for(typename Eigen::SparseMatrix<double>::InnerIterator
                                            it(M.eigen_object(), k); it; ++it) {
        outM <<  it.row() << " " << it.col() << " " << it.value() << '\n';
      }
    }

//    std::ofstream outfM("matrices/fullM.txt");
//    outfM << M.eigen_object() << std::endl;

    std::ofstream outBf("matrices/Bf.txt");
    outBf.precision(20);
    outBf << Bf.size() << std::endl;
    outBf << Bf << std::endl;

#ifdef CGAL_SMP_OUTPUT_ORBITAL_MATRICES_FOR_MATLAB
    std::ofstream mat_outM("matrices/matrixM.dat");
    mat_outM.precision(20);
    for(int k=0; k<M.eigen_object().outerSize(); ++k) {
      for(typename Eigen::SparseMatrix<double>::InnerIterator
                                            it(M.eigen_object(), k); it; ++it) {
        mat_outM <<  it.row()+1 << " " << it.col()+1 << " " << it.value() << '\n';
      }
    }

    std::ofstream matoutBf("matrices/vectorB.dat");
    matoutBf.precision(20);
    matoutBf << Bf << std::endl;
#endif

#endif

    CGAL::Timer task_timer;
    task_timer.start();

    std::cout << "Solving..." << std::endl;
    SparseLinearAlgebraTraits_d solver;
    if(!solver.linear_solver(M, Bf, Xf, D)) {
      std::cout << "Could not solve linear system" << std::endl;
      return ERROR_CANNOT_SOLVE_LINEAR_SYSTEM;
    }
    CGAL_assertion(D == 1.0);
    std::cout << "Solved the linear system in " << task_timer.time() << " seconds" << std::endl;

    Vector X(n);
    for(std::size_t i=0; i<n; ++i) {
      X[i] = Xf[i];
    }

#ifdef CGAL_SMP_OUTPUT_ORBITAL_MATRICES
    std::ofstream outf("matrices/X.txt");
    for(std::size_t i=0; i<n; ++i) {
      outf << X[i] << " ";
      outf << X[++i] << std::endl;
    }
#endif

    assign_solution(mesh, X, uvmap, vimap);

    return OK;
  }

public:
  /// Flattens the mesh to one of the orbifolds. In the end, the
  /// position of each vertex is stored in the property map `uvmap`.
  ///
  /// \param mesh a model of the `FaceGraph` concept
  /// \param bhd a halfedge on the border of the seam mesh
  /// \param cmap
  /// \param uvmap an instanciation of the class `VertexUVmap`.
  /// \param vimap an instanciation of the class `VertexIndexMap`.
  /// \param orb the type of orbifold mapping
  /// \param convexBoundary - omitted or false for free boundary, of one of
  ///        the 4 sphere orbifolds as specified in the constructor; true
  ///        for the classic Tutte embedding with fixed boundary into a disk;
  ///        'square' for "classic" Tutte on a square with prefixed boundary
  ///        map; 'freesquare' for the square disk orbifold; 'freetri' for
  ///        the triangle disk orbifold.
  /// \param wt weight type (cotan weights or MVC weights).
  ///
  /// \pre cones and seams must be valid.
  template<typename ConeMap,
           typename VertexIndexMap,
           typename VertexUVMap>
  Error_code parameterize(SeamMesh& mesh,
                          halfedge_descriptor bhd,
                          ConeMap cmap,
                          VertexUVMap uvmap,
                          VertexIndexMap vimap) const
  {
    std::cout << "Flattening" << std::endl;
    Error_code status;

    status = check_cones(cmap);
    if(status != OK) {
      return status;
    }

    // %%%%%%%%%%%%%%%%%%%%%%%
    //   Boundary conditions
    // %%%%%%%%%%%%%%%%%%%%%%%
    int lcn = number_of_linear_constraints(mesh);
    int nbVertices = static_cast<int>(num_vertices(mesh));
    Matrix A(2 * lcn, 2 * nbVertices); // change me 2 * n_of_constraints
    Vector B(2 * lcn);

    // add rotational constraints
    AddRotationalConstraint(mesh, cmap, vimap, A, B);

#ifdef CGAL_SMP_OUTPUT_ORBITAL_MATRICES
    std::cout << "A and B are filled" << std::endl;
    std::ofstream outA("matrices/A.txt"), outB("matrices/B.txt");

    for(int k=0; k<A.eigen_object().outerSize(); ++k) {
      for(Eigen::SparseMatrix<double>::InnerIterator it(A.eigen_object(), k); it; ++it) {
        outA << "(" << it.row() << ", " << it.col() << ") " << it.value() << std::endl;
      }
    }
    outB << B << std::endl;
#endif

    // %%%%%%%%%%%%%%%%%%%%
    //  Energy (Laplacian)
    // %%%%%%%%%%%%%%%%%%%%

    Matrix L(2 * nbVertices, 2 * nbVertices);
    if(weight_type == Cotangent)
      cotangent_laplacien(mesh, vimap, L);
    else // weight_type == Mean_value
      mean_value_laplacian(mesh, vimap, L);

#ifdef CGAL_SMP_OUTPUT_ORBITAL_MATRICES
    std::ofstream outL("matrices/L.txt");
    outL.precision(20);
    outL << L.row_dimension() << " " << L.column_dimension() << std::endl;
    for(int k=0; k<L.eigen_object().outerSize(); ++k) {
      for(typename Eigen::SparseMatrix<double>::InnerIterator
                                            it(L.eigen_object(), k); it; ++it) {
        outL <<  it.row() << " " << it.col() << " " << it.value() << '\n';
      }
    }
#endif

    // compute the flattening by solving the boundary conditions
    // while satisfying the convex combination property with L
    std::cout << "Solving the system ";
#ifdef CGAL_SMP_USE_SPARSESUITE_SOLVERS
    std::cout << "with a sparse linear solver from Sparsesuite." << std::endl;
#else
    std::cout << "with a sparse linear solver from Eigen." << std::endl;
#endif

    status = computeFlattening(mesh, A, B, L, uvmap, vimap);
    if(status != OK)
      return status;

    std::ofstream out("orbital_result.off");
    IO::output_uvmap_to_off(mesh, bhd, uvmap, out);

    return OK;
  }

/// Constructor
public:
  Orbital_Tutte_parameterizer_3(Orbifold_type orb_type = Square,
                                Weight_type weight_type = Cotangent)
    :
      orb_type(orb_type),
      weight_type(weight_type)
  { }
};

} // namespace Surface_mesh_parameterization

} // namespace CGAL

#endif // CGAL_SURFACE_MESH_PARAMETERIZATION_ORBITAL_TUTTE_PARAMETERIZER_3_H
