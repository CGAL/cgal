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

enum Orbifold_type
{
  Square = 0,
  Diamond,
  Triangle,
  Parallelogram
};

enum Weight_type
{
  Cotan,
  Mvc
};

template
<
  typename TriangleMesh,
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
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor  halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor      face_descriptor;

  typedef typename boost::graph_traits<TriangleMesh>::vertex_iterator      vertex_iterator;
  typedef typename boost::graph_traits<TriangleMesh>::face_iterator        face_iterator;

  // SparseLinearAlgebraTraits_d subtypes:
  typedef SparseLinearAlgebraTraits_d                               Sparse_LA;
  typedef typename Sparse_LA::Vector                                Vector;
  typedef typename Sparse_LA::Matrix                                Matrix;

  // Kernel subtypes
  typedef typename internal::Kernel_traits<TriangleMesh>::Kernel    Kernel;
  typedef typename internal::Kernel_traits<TriangleMesh>::PPM       PPM;
  typedef typename Kernel::FT                                       NT;
  typedef typename Kernel::Vector_2                                 Vector_2;
  typedef typename Kernel::Vector_3                                 Vector_3;
  typedef typename Kernel::Point_2                                  Point_2;
  typedef typename Kernel::Point_3                                  Point_3;

  Orbifold_type orb_type;

public:
  // Orbifold type functions
  std::vector<Point_2> get_cones_parameterized_coordinates() const
  {
    std::vector<Point_2> tcoords;
    if(orb_type == Square) {
      tcoords.push_back(Point_2(-1, -1));
      tcoords.push_back(Point_2(1, 1));
    } else if(orb_type == Parallelogram) {
      tcoords.push_back(Point_2(0, -0.5));
      tcoords.push_back(Point_2(1, -0.5));
      tcoords.push_back(Point_2(0, 0.5));
    } else { // if(orb_type == Diamond || orb_type == Triangle)
      tcoords.push_back(Point_2(-1, 1));
      tcoords.push_back(Point_2(-1, -1));
    }
    return tcoords;
  }

  /// Angles are minus as we go around the seam border in a counterclockwise manner
  std::vector<NT> get_angles_at_cones() const
  {
    std::vector<NT> angs;
    if(orb_type == Square) {
      angs.push_back(-4.);
      angs.push_back(-4.);
    } else if(orb_type == Diamond) {
      angs.push_back(-3.);
      angs.push_back(-3.);
    } else if(orb_type == Triangle) {
      angs.push_back(-6.);
      angs.push_back(-2.);
    } else { // if(orb_type == Parallelogram)
      angs.push_back(-2);
      angs.push_back(-1);
      angs.push_back(-2);
    }
    return angs;
  }

  // Linear system
  template<typename ConeMap,
           typename VertexIndexMap>
  void find_start_cone(const ConeMap& cmap,
                       VertexIndexMap vimap,
                       vertex_descriptor& cone,
                       int& cone_index) const
  {
    typename ConeMap::const_iterator cmit = cmap.begin(), cend = cmap.end();
    for(; cmit!=cend; ++cmit) {
      if(cmit->second != Unique_cone)
        continue;

      cone = cmit->first;
      cone_index = get(vimap, cone);
      return;
    }
  }

  /// Compute the number of linear constraints in the system.
  int number_of_linear_constraints(const TriangleMesh& mesh) const
  {
    // number of constraints for orb I, II, III is the number of edges
    // on the seam and 2 constrained cones.
    if(orb_type == Parallelogram)
      return 3 + static_cast<int>(mesh.number_of_seam_edges());
    else // orb_type == Square, Diamond, Triangle
      return 2 + static_cast<int>(mesh.number_of_seam_edges());
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

      addTransConstraints(s0, t0, s, t, current_line_id_in_A, R, A, B);
    }
  }

  /// Computes the rotational constraint on the border of the mesh.
  /// Cone constraints are also added.
  template<typename ConeMap,
           typename VertexIndexMap>
  void AddRotationalConstraint(const TriangleMesh& mesh,
                               const ConeMap& cmap,
                               VertexIndexMap vimap,
                               Matrix& A, Vector& B) const
  {
    // positions of the cones in the plane
    const std::vector<Point_2>& tcoords = get_cones_parameterized_coordinates();

    // angles at the cones
    const std::vector<NT>& angs = get_angles_at_cones();

    // the index of the line in A that we are filling next
    int current_line_id_in_A = 0.;

    // Initialize some variables used in the seam walk
    int start_cone_index = -1; // index of the beginning of the seam
    vertex_descriptor start_cone;
    find_start_cone(cmap, vimap, start_cone, start_cone_index);
    CGAL_postcondition(start_cone != vertex_descriptor() && start_cone_index != -1);

    // parameterize the initial cone
    addConstraint(A, B, current_line_id_in_A, start_cone_index,
                  1. /*entry in A*/, tcoords[0]);

    // by property of the seam mesh, the canonical halfedge that points to `hd`
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

        // the last cone of the seam is constrained
        if(is_in_map->second == Unique_cone) { // reached the end of the seam
          CGAL_assertion(hd1_target == hd2_target);
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
  void mean_value_laplacian(const TriangleMesh& mesh,
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

  /// Copy the solution into the UV property map.
  template <typename VertexIndexMap, typename VertexUVMap>
  void assign_solution(const TriangleMesh& mesh,
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
  Error_code computeFlattening(const TriangleMesh& mesh,
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
    std::ofstream outM("M.txt");
    outM.precision(20);
    outM << M.row_dimension() << " " << M.column_dimension() << std::endl;
    for(int k=0; k<M.eigen_object().outerSize(); ++k) {
      for(typename Eigen::SparseMatrix<double>::InnerIterator
                                            it(M.eigen_object(), k); it; ++it) {
        outM <<  it.row() << " " << it.col() << " " << it.value() << '\n';
      }
    }

    std::ofstream outBf("Bf.txt");
    outBf.precision(20);
    outBf << Bf.size() << std::endl;
    outBf << Bf << std::endl;
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
    std::ofstream outf("solution.txt");
    outf << X << std::endl;
#endif

    assign_solution(mesh, X, uvmap, vimap);

    return OK;
  }

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
  Error_code parameterize(const TriangleMesh& mesh,
                          halfedge_descriptor bhd,
                          ConeMap cmap,
                          VertexUVMap uvmap,
                          VertexIndexMap vimap) const
  {
    if(orb_type == Parallelogram)
      CGAL_precondition(cmap.size() == 6);
    else // orb_type == Square, Diamond, Triangle
      CGAL_precondition(cmap.size() == 4);

    std::cout << "Flattening" << std::endl;
    Error_code status;

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
    std::ofstream outA("A.txt"), outB("B.txt");

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
    mean_value_laplacian(mesh, vimap, L);

#ifdef CGAL_SMP_OUTPUT_ORBITAL_MATRICES
    std::ofstream outL("L.txt");
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
  Orbital_Tutte_parameterizer_3(Orbifold_type orb_type) : orb_type(orb_type) { }
};

} // namespace Surface_mesh_parameterization

} // namespace CGAL

#endif // CGAL_SURFACE_MESH_PARAMETERIZATION_ORBITAL_TUTTE_PARAMETERIZER_3_H
