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
#include <CGAL/Surface_mesh_parameterization/IO/File_off.h>

#include <CGAL/Surface_mesh_parameterization/Error_code.h>

#include <CGAL/circulator.h>
#include <CGAL/Eigen_solver_traits.h>
#include <CGAL/Timer.h>

#include <Eigen/Dense>
#include <Eigen/SPQRSupport>

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

// @todo checks that cones are different, are on seams, seam is one connected
//       component

namespace CGAL {

namespace Surface_mesh_parameterization {

enum Cone_type
{
  Unique_cone,
  Duplicated_cone
};

enum Orbifold_type
{
  Square,
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
    // ~1s (for bear.off)
    = Eigen_solver_traits<Eigen::SPQR<Eigen_sparse_symmetric_matrix<double>::EigenType> >

    // ~15s
//    = Eigen_solver_traits<Eigen::BiCGSTAB<Eigen_sparse_symmetric_matrix<double>::EigenType> >

    // ~20s
//    = Eigen_solver_traits< >

    // ~35s
//    = Eigen_solver_traits<Eigen::ConjugateGradient<Eigen_sparse_symmetric_matrix<double>::EigenType> >

    // Can't use SuperLU / UmfPackLU, and SparseQR / SparseLU are way too slow
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

public:
  // Linear system
  /// adds a positional constraint on a vertex x_ind, so that x_ind*w=rhs
  void addConstraint(Matrix& A, Vector& B, int ind, double w, Vector_2 rhs) const
  {
    std::cout << "Constraining " << ind << std::endl;

    A.set_coef(2*ind, 2*ind, w, true /*new_coeff*/);
    A.set_coef(2*ind + 1, 2*ind + 1, w, true /*new_coeff*/);

    B[2*ind] = rhs[0];
    B[2*ind + 1] = rhs[1];
  }

  /// adds constraints so that T * x_sinds = x_tinds, where T is a 2x2
  /// matrix, and the Transformation T is modified to affine from
  /// linear by requiring that T * x_si - x_ti = T * x_s1 - x_t1
  void addTransConstraints(int s0, int s, int t,
                           const Eigen::Matrix2d& T,
                           Matrix& A, Vector& B) const
  {
    // Matlab lines are commented for comparison.
    // Matlab fills together 2*x-1 and 2*x, but C++ fills 2*x and 2*x+1,
    // as everything (including loops!) starts at 0 and not 1.

    int t0 = s0; // for clarity

    // iterate on both rows ot the 2x2 matrix T
    for(int vert_ind=0; vert_ind<2; ++vert_ind) {
      // building up the equations by summing up the terms

      // <T(vert_ind,:), x_si>
      // obj.A(end+1, 2*sinds(ind)+[-1,0]) = T(vert_ind,:);
      A.set_coef(2*s + vert_ind, 2*s, T(vert_ind, 0), true /*new_coeff*/);
      A.set_coef(2*s + vert_ind, 2*s + 1, T(vert_ind, 1), true /*new_coeff*/);

      // -<T(vert_ind,:), x_s1>
      // obj.A(end, 2*sinds(1)+[-1,0]) = obj.A(end, 2*sinds(1)+[-1,0]) - T(vert_ind,:);
      A.add_coef(2*s + vert_ind, 2*s0, - T(vert_ind, 0));
      A.add_coef(2*s + vert_ind, 2*s0 + 1, -T(vert_ind, 1));

      //  - x_ti
      // obj.A(end, 2*tinds(ind)+vert_ind-2) = obj.A(end, 2*tinds(ind)+vert_ind-2)-1;
      A.add_coef(2*s + vert_ind, 2*t + vert_ind, -1);

      // + x_t1
      // obj.A(end, 2*tinds(1)+vert_ind-2) = obj.A(end, 2*tinds(1)+vert_ind-2)+1;
      A.add_coef(2*s + vert_ind, 2*t0 + vert_ind, 1);

      // left hand side is zero
      // obj.b=[obj.b; 0];
      B[2*s + vert_ind] = 0;
    }
  }

  /// Compute the rotational constraint on the border of the mesh.
  /// Cone constraints are also added.
  template<typename ConeMap,
           typename VertexIndexMap>
  void AddRotationalConstraint(const TriangleMesh& mesh,
                               const ConeMap& cmap,
                               VertexIndexMap vimap,
                               Matrix& A, Vector& B) const
  {
    // positions of the cones in the plane TMP
    std::vector<Vector_2> tcoords(2);
    tcoords[0] = Vector_2(-1, -1);
    tcoords[1] = Vector_2(1, 1);

    // rotations of the seams TMP
    // angles at the cones
    std::vector<NT> angs(2);
    angs[0] = 4; // singularity is [4; 4] for orb1
    angs[1] = 4; // singularity is [4; 4] for orb1

    // Go through the seam
    boost::unordered_set<int> constrained_cones;
    for(int i=0; i<2; ++i) { // TMP
      // Find a non-duplicated cone that has not yet been constrained
      vertex_descriptor start_cone;
      int start_cone_index = -1;

      BOOST_FOREACH(vertex_descriptor vd, vertices(mesh)) {
        typename ConeMap::const_iterator it = cmap.find(vd);
        if(it != cmap.end() && it->second == Unique_cone) {
          start_cone_index = get(vimap, vd);
          if(constrained_cones.find(start_cone_index) == constrained_cones.end() ) {
            start_cone = vd;
            constrained_cones.insert(start_cone_index);
            break;
          }
        }
      }
      CGAL_postcondition(start_cone != vertex_descriptor() && start_cone_index != -1);

      std::cout << "starting from " << start_cone_index << std::endl;

      // Mark 'start_cone' as constraint
      addConstraint(A, B, start_cone_index, 1. /*entry in A*/, tcoords[i]);

      // Go through the seam, marking rotation and cone constraints
      halfedge_descriptor hd = halfedge(start_cone, mesh);
      CGAL_precondition(mesh.has_on_seam(hd));
      halfedge_descriptor bhd = opposite(hd, mesh);
      CGAL_precondition(is_border(bhd, mesh));

      // the rotation angle
      double ang = angs[i];

      // the rotation matrix according to the angle 'ang'
      Eigen::Matrix2d R;
      R(0,0) = std::cos(2 * CGAL_PI / ang); R(0,1) = - std::sin(2 * CGAL_PI / ang);
      R(1,0) = std::sin(2 * CGAL_PI / ang); R(1,1) = std::cos(2 * CGAL_PI / ang);

      // move through the seam from start_cone_index till we reach a duplicated cone
      while(true) {
        // Get the two halfedges on each side of the stream
        halfedge_descriptor hd1 = bhd; // only for clarity

        // the non-border halfedge with same vertices (in the underlying mesh of the seam
        // mesh) as bhd is simply bhd with the 'seam' boolean set to false
        halfedge_descriptor hd2 = bhd;
        hd2.seam = false;

        // Add the rotational constraint between the two halfedges on the seam
        vertex_descriptor hd1_target = target(hd1, mesh);
        vertex_descriptor hd2_target = target(hd2, mesh);
        int hd1t_index = get(vimap, hd1_target);
        int hd2t_index = get(vimap, hd2_target);

        std::cout << "hd1/hd2: " << hd1t_index << " " << hd2t_index << std::endl;

        addTransConstraints(start_cone_index, hd1t_index, hd2t_index, R, A, B);

        // Check if we have reached the duplicated cone
        vertex_descriptor bhd_target = target(bhd, mesh); // also known as 'hd1_target'
        typename ConeMap::const_iterator is_in_map = cmap.find(bhd_target);
        if(is_in_map != cmap.end()) {
          // starting from a unique cone, the next cone must be the duplicated cone
          CGAL_assertion(is_in_map->second == Duplicated_cone);
          break;
        }

        // move to the next halfedge couple (walking on the border of the seam)
        bhd = next(bhd, mesh);
        CGAL_postcondition(mesh.has_on_seam(bhd) && is_border(bhd, mesh));
      }
    }
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

  /// Compute the coefficients of the mean value Laplacian matrix for the edge
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
    NT w_i_base = -1.0 * compute_w_ij_mvc(pi, pj, pk);

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
    CGAL_assertion(X.size() == 2*num_vertices(mesh) );

    BOOST_FOREACH(vertex_descriptor vd, vertices(mesh)) {
      int index = get(vimap, vd);
      NT u = X(2*index);
      NT v = X(2*index + 1);

      put(uvmap, vd, Point_2(u, v));
    }
  }

  /// Solve the linear system.
  template <typename VertexUVMap,
            typename VertexIndexMap>
  Error_code computeFlattening(const TriangleMesh& mesh,
                               const Matrix& A, const Vector& B,
                               const Matrix& L,
                               VertexUVMap uvmap, VertexIndexMap vimap) const
  {
    std::size_t n = B.size();

    // Fill the large system M.Xf = Bf:
    // ( L A' ) ( Xf ) = ( B )
    // ( A 0  ) ( Xf ) = ( 0 )
    Matrix M(2*n, 2*n);
    Vector Bf(2*n);
    NT D;
    Vector Xf(2*n);

    // full B
    for(std::size_t i=0; i<n; ++i) {
      Bf[n+i] = B[i];
    }

    std::cout << "filled Bf" << std::endl;

    // matrix M
    for(int k=0; k<L.eigen_object().outerSize(); ++k) {
      for(typename Eigen::SparseMatrix<double>::InnerIterator
                                            it(L.eigen_object(), k); it; ++it) {
        M.set_coef(it.row(), it.col(), it.value(), true /*new_coeff*/);
//        std::cout <<  it.row() << " " << it.col() << " " << it.value() << '\n';
      }
    }

    std::cout << "filled topleft of M" << std::endl;

    for(int k=0; k<A.eigen_object().outerSize(); ++k) {
      for(typename Eigen::SparseMatrix<double>::InnerIterator
                                            it(A.eigen_object(), k); it; ++it) {
        M.set_coef(it.col(), it.row() + n, it.value(), true /*new_coeff*/); // A
        M.set_coef(it.row() + n, it.col(), it.value(), true /*new_coeff*/); // A'
//        std::cout <<  it.row() << " " << it.col() << " " << it.value() << '\n';
      }
    }

    std::cout << "Filled M and Bf" << std::endl;

#ifdef SMP_OUTPUT_ORBITAL_MATRICES
    std::ofstream outM("linears_M.txt");
    outM << M.eigen_object() << std::endl;
#endif

    CGAL::Timer task_timer;
    task_timer.start();

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

#ifdef SMP_OUTPUT_ORBITAL_MATRICES
    std::ofstream outf("solution.txt");
    outf << X << std::endl;
#endif

    assign_solution(mesh, X, uvmap, vimap);

    return OK;
  }

  /// Flatten the mesh to one of the orbifolds. In the end, the
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
    std::cout << "Flattening" << std::endl;
    Error_code status;

    // %%%%%%%%%%%%%%%%%%%%%%%
    //   Boundary conditions
    // %%%%%%%%%%%%%%%%%%%%%%%
    std::size_t nbVertices = num_vertices(mesh);
    Matrix A(2 * nbVertices, 2 * nbVertices);
    Vector B(2 * nbVertices);

    // add rotational constraints
    AddRotationalConstraint(mesh, cmap, vimap, A, B);

#ifdef SMP_OUTPUT_ORBITAL_MATRICES
    std::cout << "A and B are filled" << std::endl;
    std::ofstream outA("A.txt"), outB("B.txt");

    for(int k=0; k<A.eigen_object().outerSize(); ++k) {
      for(Eigen::SparseMatrix<double>::InnerIterator it(A.eigen_object(), k); it; ++it) {
        outA << "(" << it.row() << ", " << it.col() << ") " << it.value() << std::endl;
      }
    }
    outA << std::endl << A.eigen_object() << std::endl << std::endl;

    outB << B << std::endl;
#endif

    // %%%%%%%%%%%%%%%%%%%%
    //  Energy (Laplacian)
    // %%%%%%%%%%%%%%%%%%%%

    Matrix L(2*nbVertices, 2*nbVertices);
    mean_value_laplacian(mesh, vimap, L);

#ifdef SMP_OUTPUT_ORBITAL_MATRICES
    std::ofstream outL("MVCc.txt");
    outL << L.eigen_object() << std::endl;
#endif

    // compute the flattening by solving the boundary conditions
    // while satisfying the convex combination property with L
    std::cout << "solving system" << std::endl;
    status = computeFlattening(mesh, A, B, L, uvmap, vimap);
    if(status != OK)
      return status;

    std::ofstream out("orbital_result.off");
    IO::output_uvmap_to_off(mesh, bhd, uvmap, out);

    return OK;
  }
};

} // namespace Surface_mesh_parameterization

} // namespace CGAL

#endif // CGAL_SURFACE_MESH_PARAMETERIZATION_ORBITAL_TUTTE_PARAMETERIZER_3_H
