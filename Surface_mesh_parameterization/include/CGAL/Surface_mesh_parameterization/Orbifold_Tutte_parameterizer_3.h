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
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_SURFACE_MESH_PARAMETERIZATION_ORBIFOLD_TUTTE_PARAMETERIZER_3_H
#define CGAL_SURFACE_MESH_PARAMETERIZATION_ORBIFOLD_TUTTE_PARAMETERIZER_3_H

#include <CGAL/license/Surface_mesh_parameterization.h>

#include <CGAL/Surface_mesh_parameterization/internal/angles.h>
#include <CGAL/Surface_mesh_parameterization/internal/kernel_traits.h>
#include <CGAL/Surface_mesh_parameterization/internal/orbifold_cone_helper.h>
#include <CGAL/Surface_mesh_parameterization/IO/File_off.h>

#include <CGAL/Surface_mesh_parameterization/Error_code.h>

#include <CGAL/Polygon_mesh_processing/Weights.h>

#include <CGAL/assertions.h>
#include <CGAL/circulator.h>
#include <CGAL/Default.h>
#include <CGAL/Timer.h>

#if defined(CGAL_EIGEN3_ENABLED)
#include <CGAL/Eigen_solver_traits.h>
#ifdef CGAL_SMP_USE_SPARSESUITE_SOLVERS
#include <Eigen/UmfPackSupport>
#endif
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

/// \file Orbifold_Tutte_parameterizer_3.h

// @todo checks that cones are different, are on seams, seam is one connected
//       component
// @todo Should the order of cones provided in entry matter ? Map the first cone
//       to [-1, -1] for example ?

namespace CGAL {

namespace Surface_mesh_parameterization {

/// \ingroup PkgSurfaceParameterizationEnums
///
/// Weight type used in the parameterization computation.
///
/// MVC weights are guaranteed to generate positive edge weights, and the parameterization
/// is guaranteed to be injective.
///
/// In case the cotangent weights are used, the orbifold-Tutte embedding globally
/// minimizes the Dirichlet energy and approximates conformal mappings.
enum Weight_type
{
  Cotangent = 0,
  Mean_value
};

/// \ingroup  PkgSurfaceParameterizationMethods
///
/// The class `Orbifold_Tutte_parameterizer_3` implements <em>Orbifold Tutte Planar
/// Embeddings</em> \cgalCite{aigerman2015orbifold}.
///
/// This is a borderless parameterization. A one-to-one mapping is guaranteed.
///
/// The main function of the class `Orbifold_Tutte_parameterizer_3` is `parameterize()`,
/// to which the user provides a `Seam_mesh` with marked edges (the seams)
/// and a set of vertices of the mesh (the cones). The choice of cones influences
/// the resulting parameterization, but not the choice of the seam path between these cones.
///
/// The example \ref Surface_mesh_parameterization/orbifold.cpp "orbifold.cpp"
/// shows how to select cones on the input mesh and automatically construct
/// the seams and the cones on the `Seam_mesh`.
///
/// \cgalModels `Parameterizer_3`
///
/// \tparam SeamMesh must be a `Seam_mesh`, with underlying mesh any model of `FaceListGraph` and `HalfedgeListGraph`.
///
/// \tparam SolverTraits_ must be a model of `SparseLinearAlgebraTraits_d`.<br>
///         <b>%Default:</b> If \ref thirdpartyEigen "Eigen" 3.1 (or greater) is available
///         and `CGAL_EIGEN3_ENABLED` is defined, then an overload of `Eigen_solver_traits`
///         is provided as default parameter:
/// \code
///   CGAL::Eigen_solver_traits<
///           Eigen::SparseLU<Eigen_sparse_matrix<double>::EigenType> >
/// \endcode
///         Moreover, if SparseSuite solvers are available, which is greatly preferable for speed,
///         then the default parameter is:
/// \code
///   CGAL::Eigen_solver_traits<
///           Eigen::UmfPackLU<Eigen_sparse_matrix<double>::EigenType> >
/// \endcode
///
/// \sa `CGAL::Surface_mesh_parameterization::ARAP_parameterizer_3<TriangleMesh, BorderParameterizer, SolverTraits>`
/// \sa `CGAL::Surface_mesh_parameterization::Barycentric_mapping_parameterizer_3<TriangleMesh, BorderParameterizer, SolverTraits>`
/// \sa `CGAL::Surface_mesh_parameterization::Discrete_authalic_parameterizer_3<TriangleMesh, BorderParameterizer, SolverTraits>`
/// \sa `CGAL::Surface_mesh_parameterization::Discrete_conformal_map_parameterizer_3<TriangleMesh, BorderParameterizer, SolverTraits>`
/// \sa `CGAL::Surface_mesh_parameterization::LSCM_parameterizer_3<TriangleMesh, BorderParameterizer>`
/// \sa `CGAL::Surface_mesh_parameterization::Mean_value_coordinates_parameterizer_3<TriangleMesh, BorderParameterizer, SolverTraits>`
///
template < typename SeamMesh,
           typename SolverTraits_ = Default>
class Orbifold_Tutte_parameterizer_3
{
public:
#ifndef DOXYGEN_RUNNING
  typedef typename Default::Get<
    SolverTraits_,
  #if defined(CGAL_EIGEN3_ENABLED)
    #ifdef CGAL_SMP_USE_SPARSESUITE_SOLVERS
      CGAL::Eigen_solver_traits<
        Eigen::UmfPackLU<Eigen_sparse_matrix<double>::EigenType> >
    #else
      CGAL::Eigen_solver_traits<
        Eigen::SparseLU<Eigen_sparse_matrix<double>::EigenType> >
    #endif
  #else
    #pragma message("Error: You must either provide 'SolverTraits_' or link CGAL with the Eigen library")
    SolverTraits_ // no parameter provided, and Eigen is not enabled: so don't compile!
  #endif
  >::type                                                     Solver_traits;
#else
  typedef SolverTraits_                                       Solver_traits;
#endif

private:
  typedef typename boost::graph_traits<SeamMesh>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::graph_traits<SeamMesh>::halfedge_descriptor  halfedge_descriptor;
  typedef typename boost::graph_traits<SeamMesh>::face_descriptor      face_descriptor;

  typedef typename boost::graph_traits<SeamMesh>::vertex_iterator      vertex_iterator;
  typedef typename boost::graph_traits<SeamMesh>::face_iterator        face_iterator;

  // Solver traits subtypes:
  typedef typename Solver_traits::Vector                               Vector;
  typedef typename Solver_traits::Matrix                               Matrix;

  // Kernel subtypes
  typedef typename internal::Kernel_traits<SeamMesh>::Kernel        Kernel;
  typedef typename internal::Kernel_traits<SeamMesh>::PPM           PPM;
  typedef typename Kernel::FT                                       NT;
  typedef typename Kernel::Vector_2                                 Vector_2;
  typedef typename Kernel::Vector_3                                 Vector_3;
  typedef typename Kernel::Point_2                                  Point_2;
  typedef typename Kernel::Point_3                                  Point_3;

  const Orbifold_type orb_type;
  const Weight_type weight_type;

private:
  // Check input's correctness.
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

    return OK;
  }

  // Compute the number of linear constraints in the system.
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

  // Adds a positional constraint on a vertex x_ind, so that x_ind * w = rhs.
  void addConstraint(Matrix& M, Vector& B, int& id_r, int id_c, double w, Point_2 rhs) const
  {
    M.set_coef(2*id_r, 2*id_c, w, true /*new_coef*/);
    M.set_coef(2*id_r + 1, 2*id_c + 1, w, true /*new_coef*/);

    // Since we are filling the big system M.Xf = B with
    // ( L A' ) ( Xf ) = ( C )
    // ( A 0  ) ( Xf ) = ( 0 )
    // we already add the transposed here:
    M.set_coef(2*id_c, 2*id_r, w, true /*new_coef*/);
    M.set_coef(2*id_c + 1, 2*id_r + 1, w, true /*new_coef*/);

    B[2*id_r] = rhs.x();
    B[2*id_r + 1] = rhs.y();

    ++id_r; // current line index in A is increased
  }

  // Adds constraints so that T * x_sinds = x_tinds, where T is a 2x2
  // matrix, and the Transformation T is modified to affine from
  // linear by requiring that T * x_si - x_ti = T * x_s1 - x_t1.
  void addTransConstraints(int s0, int t0, int s, int t,
                           int& id_r,
                           const std::vector<double>& T,
                           Matrix& M, Vector& B) const
  {
    // Everything is duplicated since we are filling the big system M.Xf = B with
    // ( L A' ) ( Xf ) = ( C )
    // ( A 0  ) ( Xf ) = ( 0 )

    // Iterate on both rows ot the 2x2 matrix T
    for(int vert_ind=0; vert_ind<2; ++vert_ind) {
      // building up the equations by summing up the terms

      // Matlab lines are commented for comparison.
      // Matlab fills together 2*x-1 and 2*x, but C++ fills 2*x and 2*x+1,
      // as everything (including loops!) starts at 0 and not 1.

      // <T(vert_ind,:), x_si>
      // obj.A(end+1, 2*sinds(ind)+[-1,0]) = T(vert_ind,:);
      M.set_coef(2*id_r + vert_ind, 2*s, T[2 * vert_ind], true /*new_coef*/);
      M.set_coef(2*id_r + vert_ind, 2*s + 1, T[2 * vert_ind + 1], true /*new_coef*/);

      M.set_coef(2*s, 2*id_r + vert_ind, T[2 * vert_ind], true /*new_coef*/);
      M.set_coef(2*s + 1, 2*id_r + vert_ind, T[2 * vert_ind + 1], true /*new_coef*/);

      // -<T(vert_ind,:), x_s1>
      // obj.A(end, 2*sinds(1)+[-1,0]) = obj.A(end, 2*sinds(1)+[-1,0]) - T(vert_ind,:);
      M.add_coef(2*id_r + vert_ind, 2*s0, - T[2 * vert_ind]);
      M.add_coef(2*id_r + vert_ind, 2*s0 + 1, - T[2 * vert_ind + 1]);

      M.add_coef(2*s0, 2*id_r + vert_ind, - T[2 * vert_ind]);
      M.add_coef(2*s0 + 1, 2*id_r + vert_ind, - T[2 * vert_ind + 1]);

      //  - x_ti
      // obj.A(end, 2*tinds(ind)+vert_ind-2) = obj.A(end, 2*tinds(ind)+vert_ind-2)-1;
      M.add_coef(2*id_r + vert_ind, 2*t + vert_ind, -1);

      M.add_coef(2*t + vert_ind, 2*id_r + vert_ind, -1);

      // + x_t1
      // obj.A(end, 2*tinds(1)+vert_ind-2) = obj.A(end, 2*tinds(1)+vert_ind-2)+1;
      M.add_coef(2*id_r + vert_ind, 2*t0 + vert_ind, 1);

      M.add_coef(2*t0 + vert_ind, 2*id_r + vert_ind, 1);

      // left hand side is zero
      // obj.b=[obj.b; 0];
      B[2*id_r + vert_ind] = 0;
    }

    ++id_r; // current line index in M is increased
  }

  // Add the constraints from a seam segment to the linear system.
  void constrain_seam_segment(const std::vector<std::pair<int, int> >& seam_segment,
                              NT ang, int& current_line_id_in_M,
                              Matrix& M, Vector& B) const
  {
    // check that if there is a common vertex, it is at the beginning
    const bool is_reversed = (seam_segment.back().first == seam_segment.back().second);

    if(is_reversed) {
      ang *= -1;
    }

    // The rotation matrix according to the angle 'ang'. Put in a vector
    // because we need to access it later and Matrix does not provide read access...
    std::vector<double> R(4);
    R[0] = std::cos(2 * CGAL_PI / ang);
    R[1] = - std::sin(2 * CGAL_PI / ang);
    R[2] = std::sin(2 * CGAL_PI / ang);
    R[3] = std::cos(2 * CGAL_PI / ang);

    const int s0 = is_reversed ? seam_segment.back().first : seam_segment.front().first;
    const int t0 = is_reversed ? seam_segment.back().second : seam_segment.front().second;

    typename std::vector<std::pair<int, int> >::const_iterator it = seam_segment.begin(),
                                                               end = seam_segment.end();

    // ignore the first entry of the seam segment (they correspond to a constrained point)
    if(is_reversed)
      --end;
    else
      ++it;

    for(; it!=end; ++it) {
      const int s = it->first;
      const int t = it->second;
      CGAL_assertion(s != t);

      // sending s to t (and _not_ t to s !)
      addTransConstraints(t0, s0, t, s, current_line_id_in_M, R, M, B);
    }
  }

  // Computes the rotational constraint on the border of the mesh.
  // Cone constraints are also added.
  template<typename ConeMap,
           typename VertexIndexMap>
  void AddRotationalConstraint(const SeamMesh& mesh,
                               const ConeMap& cmap,
                               VertexIndexMap vimap,
                               Matrix& M, Vector& B) const
  {
    // positions of the cones in the plane
    typedef std::vector<Point_2>                      Point_container;
    const Point_container& tcoords =
                  get_cones_parameterized_coordinates<Point_container>(orb_type);

    // angles at the cones
    typedef std::vector<NT>                           Angle_container;
    const Angle_container& angs = get_angles_at_cones<Angle_container>(orb_type);

    // The index of the line in M that we are filling next.

    // Since we are filling the big system M.Xf = B with
    // ( L A' ) ( Xf ) = ( C )
    // ( A 0  ) ( Xf ) = ( 0 )
    // we do not start at 0, but at the first line below the matrix L
    // (note that this should thus be 2*num_vertices, but in the filling functions
    // we use 2*line_number to fill two at the time...)
    int current_line_id_in_M = static_cast<int>(num_vertices(mesh));
    CGAL_postcondition_code(int initial_line_id = current_line_id_in_M;)

    // Initialize some variables used in the seam walk
    int start_cone_index = -1; // index of the beginning of the seam
    vertex_descriptor start_cone;
    internal::find_start_cone(cmap, vimap, start_cone, start_cone_index);
    CGAL_postcondition(start_cone != vertex_descriptor() && start_cone_index != -1);

    // parameterize the initial cone
    addConstraint(M, B, current_line_id_in_M, start_cone_index,
                  1. /*entry in M*/, tcoords[0]);

    // by property of the seam mesh, the canonical halfedge that points to start_cone
    // is on the seam, and is not on the border
    const halfedge_descriptor hd = halfedge(start_cone, mesh);
    CGAL_precondition(mesh.has_on_seam(hd));
    halfedge_descriptor bhd = opposite(hd, mesh);
    CGAL_precondition(is_border(bhd, mesh));

    // points between two cones, and the corresponding points on the opposite side of the seam
    std::vector<std::pair<int, int> > seam_segment;
    std::size_t segment_index = 0; // counting the segments (3 max)

    // Go through the seam, marking rotation and cone constraints
    while(true) { // breaking at the last cone
      // Get the two halfedges on each side of the stream
      const halfedge_descriptor hd1 = bhd; // only for clarity

      // the non-border halfedge with same vertices (in the underlying mesh of the seam
      // mesh) as bhd is simply bhd with the 'seam' boolean set to false
      const halfedge_descriptor hd2(bhd, false /* not on seam*/);

      // Compute the corresponding indices
      const vertex_descriptor hd1_source = source(hd1, mesh);
      const vertex_descriptor hd2_source = source(hd2, mesh);
      const int hd1s_index = get(vimap, hd1_source);
      const int hd2s_index = get(vimap, hd2_source);

      // If orbifold type IV and it is second cone in flattening, add constraint
      if(orb_type == Parallelogram && cmap.find(hd1_source) != cmap.end()
                                   && segment_index == 1) {
        addConstraint(M, B, current_line_id_in_M, hd1s_index,
                      1. /*entry in M*/, tcoords[1]);
      }

      // Add the pair to the seam segment
      seam_segment.push_back(std::make_pair(hd1s_index, hd2s_index));

      // Check if we have reached a cone
      const vertex_descriptor bhd_target = target(bhd, mesh);
      typename ConeMap::const_iterator is_in_map = cmap.find(bhd_target);
      if(is_in_map != cmap.end()) {
        // add the target to finish the seam segment
        const vertex_descriptor hd1_target = target(hd1, mesh);
        const vertex_descriptor hd2_target = target(hd2, mesh);
        const int hd1t_index = get(vimap, hd1_target);
        const int hd2t_index = get(vimap, hd2_target);

        seam_segment.push_back(std::make_pair(hd1t_index, hd2t_index));

        CGAL_assertion(segment_index < angs.size());
        NT ang = angs[segment_index];
        constrain_seam_segment(seam_segment, ang, current_line_id_in_M, M, B);

        // Check if we have reached the end of the seam
        if(is_in_map->second == Second_unique_cone) {
          CGAL_assertion(hd1_target == hd2_target);
        // the last cone of the seam is constrained
          addConstraint(M, B, current_line_id_in_M, hd1t_index,
                        1. /*entry in M*/, tcoords.back());
          break;
        }

        seam_segment.clear();
        ++segment_index;
      }

      // move to the next halfedge couple (walking on the border of the seam)
      bhd = next(bhd, mesh);
      CGAL_postcondition(mesh.has_on_seam(bhd) && is_border(bhd, mesh));
    }

    CGAL_postcondition(current_line_id_in_M - initial_line_id == number_of_linear_constraints(mesh));
  }

  // MVC computations
  NT compute_w_ij_mvc(const Point_3& pi, const Point_3& pj, const Point_3& pk) const
  {
    //                                                               ->     ->
    // Compute the angle (pj, pi, pk), the angle between the vectors ij and ik
    const NT angle = internal::compute_angle_rad<Kernel>(pj, pi, pk);
    const NT weight = std::tan(0.5 * angle);

    return weight;
  }

  // Computes the coefficients of the mean value Laplacian matrix for the edge.
  // `ij` in the face `ijk`
  void fill_mvc_matrix(const Point_3& pi, int i,
                       const Point_3& pj, int j,
                       const Point_3& pk, int k, Matrix& M) const
  {
    // For MVC, the entry of M(i,j) is - [ tan(gamma_ij/2) + tan(delta_ij)/2 ] / |ij|
    // where gamma_ij and delta_ij are the angles at i around the edge ij

    // This function computes the angle alpha at i, and add
    // -- M(i,j) += tan(alpha / 2) / |ij|
    // -- M(i,k) += tan(alpha / 2) / |ik|
    // -- M(i,i) -= M(i,j) + M(i,k)

    // The other parts of M(i,j) and M(i,k) will be added when this function
    // is called from the neighboring faces of F_ijk that share the vertex i

    // Compute: - tan(alpha / 2)
    const NT w_i_base = 1.0 * compute_w_ij_mvc(pi, pj, pk);

    // @fixme unefficient: lengths are computed (and inversed!) twice per edge

    // Set w_ij in matrix
    const Vector_3 edge_ij = pi - pj;
    const NT len_ij = CGAL::sqrt(edge_ij * edge_ij);
    CGAL_assertion(len_ij != 0.0); // two points are identical!
    const NT w_ij = w_i_base / len_ij;
    M.add_coef(2*i, 2*j, w_ij);
    M.add_coef(2*i +1, 2*j + 1, w_ij);

    // Set w_ik in matrix
    Vector_3 edge_ik = pi - pk;
    const NT len_ik = CGAL::sqrt(edge_ik * edge_ik);
    CGAL_assertion(len_ik != 0.0); // two points are identical!
    const NT w_ik = w_i_base / len_ik;
    M.add_coef(2*i, 2*k, w_ik);
    M.add_coef(2*i + 1, 2*k + 1, w_ik);

    // Add to w_ii (w_ii = - sum w_ij)
    const NT w_ii = - w_ij - w_ik;
    M.add_coef(2*i, 2*i, w_ii);
    M.add_coef(2*i + 1, 2*i + 1, w_ii);
  }

  // Compute the mean value Laplacian matrix.
  template<typename VertexIndexMap>
  void mean_value_laplacian(const SeamMesh& mesh,
                            VertexIndexMap vimap,
                            Matrix& M) const
  {
    const PPM ppmap = get(vertex_point, mesh);

    BOOST_FOREACH(face_descriptor fd, faces(mesh)) {
      const halfedge_descriptor hd = halfedge(fd, mesh);

      const vertex_descriptor vd_i = target(hd, mesh);
      const vertex_descriptor vd_j = target(next(hd, mesh), mesh);
      const vertex_descriptor vd_k = source(hd, mesh);
      const Point_3& pi = get(ppmap, vd_i);
      const Point_3& pj = get(ppmap, vd_j);
      const Point_3& pk = get(ppmap, vd_k);
      const int i = get(vimap, vd_i);
      const int j = get(vimap, vd_j);
      const int k = get(vimap, vd_k);

      fill_mvc_matrix(pi, i, pj, j, pk, k, M);
      fill_mvc_matrix(pj, j, pk, k, pi, i, M);
      fill_mvc_matrix(pk, k, pi, i, pj, j, M);
    }
  }

  // Compute the system weights using a Cotangent Laplacian.
  template<typename VertexIndexMap>
  void cotangent_laplacien(SeamMesh& mesh,
                           VertexIndexMap vimap,
                           Matrix& M) const
  {
    const PPM ppmap = get(vertex_point, mesh);

    // not exactly sure which cotan weights should be used:
    // 0.5 (cot a + cot b) ? 1/T1 cot a + 1/T2 cot b ? 1/Vor(i) (cot a + cot b?)
    // Comparing to the matlab code, the basic Cotangent_weight gives the same results.
    typedef CGAL::internal::Cotangent_weight<SeamMesh>                      Cotan_weights;
//    typedef CGAL::internal::Cotangent_weight_with_triangle_area<SeamMesh>   Cotan_weights;

    Cotan_weights cotan_weight_calculator(mesh, ppmap);

    BOOST_FOREACH(halfedge_descriptor hd, halfedges(mesh)) {
      const vertex_descriptor vi = source(hd, mesh);
      const vertex_descriptor vj = target(hd, mesh);
      const int i = get(vimap, vi);
      const int j = get(vimap, vj);

      if(i > j)
        continue;

      // times 2 because Cotangent_weight returns 1/2 (cot alpha + cot beta)...
      const NT w_ij = 2 * cotan_weight_calculator(hd);

      // ij
      M.set_coef(2*i, 2*j, w_ij, true /* new coef */);
      M.set_coef(2*i +1, 2*j + 1, w_ij, true /* new coef */);

      // ji
      M.set_coef(2*j, 2*i, w_ij, true /* new coef */);
      M.set_coef(2*j +1, 2*i + 1, w_ij, true /* new coef */);

      // ii
      M.add_coef(2*i, 2*i, - w_ij);
      M.add_coef(2*i + 1, 2*i + 1, - w_ij);

      // jj
      M.add_coef(2*j, 2*j, - w_ij);
      M.add_coef(2*j + 1, 2*j + 1, - w_ij);
    }
  }

  // Copy the solution into the UV property map.
  template <typename VertexIndexMap, typename VertexUVMap>
  void assign_solution(const SeamMesh& mesh,
                       const Vector& X,
                       VertexUVMap uvmap,
                       const VertexIndexMap vimap) const
  {
    CGAL_assertion(X.dimension() == static_cast<int>(2 * num_vertices(mesh)));

    BOOST_FOREACH(vertex_descriptor vd, vertices(mesh)) {
      const int index = get(vimap, vd);
      const NT u = X(2*index);
      const NT v = X(2*index + 1);

      put(uvmap, vd, Point_2(u, v));
    }
  }

  // Solves the linear system.
  template <typename VertexUVMap,
            typename VertexIndexMap>
  Error_code computeFlattening(const SeamMesh& mesh,
                               const Matrix& M, const Vector& B,
                               VertexUVMap uvmap, VertexIndexMap vimap) const
  {
    CGAL_precondition(M.row_dimension() == M.column_dimension());
    CGAL_precondition(M.row_dimension() == B.dimension());

    const int big_n = M.row_dimension();
    const std::size_t n = 2 * num_vertices(mesh);

    NT D;
    Vector Xf(big_n);

    CGAL::Timer task_timer;
    task_timer.start();

    Solver_traits solver;
    if(!solver.linear_solver(M, B, Xf, D)) {
      std::cerr << "Could not solve linear system" << std::endl;
      return ERROR_CANNOT_SOLVE_LINEAR_SYSTEM;
    }
    CGAL_assertion(D == 1.0);

    Vector X(n);
    for(std::size_t i=0; i<n; ++i) {
      X[i] = Xf[i];
    }

#ifdef CGAL_SMP_OUTPUT_ORBIFOLD_MATRICES
    std::ofstream outf("matrices/X.txt");
    for(std::size_t i=0; i<n; ++i) {
      outf << X[i] << " ";
      outf << X[++i] << '\n';
    }
#endif

    assign_solution(mesh, X, uvmap, vimap);
    return OK;
  }

public:
  /// Compute a one-to-one mapping from a triangular 3D surface mesh
  /// to a piece of the 2D space.
  /// The mapping is piecewise linear (linear in each triangle).
  /// The result is the (u,v) pair image of each vertex of the 3D surface.
  ///
  /// \tparam VertexUVmap must be a model of `ReadWritePropertyMap` with
  ///         `boost::graph_traits<Seam_mesh>::%vertex_descriptor` as key type and
  ///         %Point_2 (type deduced from `Seam_mesh` using `Kernel_traits`)
  ///         as value type.
  /// \tparam VertexIndexMap must be a model of `ReadablePropertyMap` with
  ///         `boost::graph_traits<Seam_mesh>::%vertex_descriptor` as key type and
  ///         a unique integer as value type.
  /// \tparam VertexParameterizedMap must be a model of `ReadWritePropertyMap` with
  ///         `boost::graph_traits<Seam_mesh>::%vertex_descriptor` as key type and
  ///         a Boolean as value type.
  ///
  /// \param mesh a `Seam_mesh` parameterized by any model of a `FaceGraph`
  /// \param bhd a halfedge on the border of the seam mesh
  /// \param cmap a mapping of the `vertex_descriptor`s of `mesh` that are cones
  ///             to their respective `Cone_type`.
  /// \param uvmap an instanciation of the class `VertexUVmap`.
  /// \param vimap an instanciation of the class `VertexIndexMap`.
  ///
  /// \pre `mesh` must be a triangular mesh.
  /// \pre The underlying mesh of `mesh` is a topological ball.
  /// \pre The vertices must be indexed (vimap must be initialized).
  /// \pre The cones are vertices of `mesh` and their number is adapted to
  ///      the orbifold type (4 for types I, II or III and 6 for type IV).
  /// \pre The seam edges form a set of segments that contains the different cones,
  ///      and starts and ends at two different cones.
  /// \pre The seam edges form a set of segments that is homotopic to a
  ///      line. Specifically, the paths between cones must not self intersect
  ///      or intersect other paths (see Figure below).
  ///
  /// \cgalFigureBegin{Surface_mesh_parameterizationfigorbifold, orbifold_path.svg}
  /// Bad (left) and good (right) seam paths. The seam edges are shown in dark red.
  /// Cones are marked in yellow and blue.
  /// \cgalFigureEnd
  template<typename ConeMap,
           typename VertexIndexMap,
           typename VertexUVMap>
  Error_code parameterize(SeamMesh& mesh,
                          halfedge_descriptor bhd,
                          ConeMap cmap,
                          VertexUVMap uvmap,
                          VertexIndexMap vimap) const
  {
    Error_code status;

    status = check_cones(cmap);
    if(status != OK) {
      return status;
    }

    const int lcn = number_of_linear_constraints(mesh);
    const int nbVertices = static_cast<int>(num_vertices(mesh));

    // To avoid concatenating matrices, we will write everything directly
    // in the big system that we will (eventually solve). It is:

    //  M.Xf = B with
    // ( L A' ) ( Xf ) = ( C )
    // ( A 0  ) ( Xf ) = ( 0 )

    // where:
    // -- L is a (2 * nbVertices, 2 * nbVertices) Laplacian matrix
    // -- A is a (2 * lcn,        2 * nbVertices) constraint matrix
    // -- C is a (2 * lcn) right hand side vector

    // Thus M is a (2 * (nbVertices + lcn), 2 * (nbVertices + lcn)) matrix
    const int total_size = 2 * (lcn + nbVertices);
    Matrix M(total_size, total_size);
    Vector B(total_size);

    // %%%%%%%%%%%%%%%%%%%%%%%
    //   Boundary conditions
    // %%%%%%%%%%%%%%%%%%%%%%%

    // add rotational constraints
    AddRotationalConstraint(mesh, cmap, vimap, M, B);

    // %%%%%%%%%%%%%%%%%%%%
    //  Energy (Laplacian)
    // %%%%%%%%%%%%%%%%%%%%
    if(weight_type == Cotangent)
      cotangent_laplacien(mesh, vimap, M);
    else // weight_type == Mean_value
      mean_value_laplacian(mesh, vimap, M);

#ifdef CGAL_SMP_OUTPUT_ORBIFOLD_MATRICES
    std::ofstream outM("matrices/M.txt");
    outM.precision(20);
    outM << total_size << " " << total_size << std::endl;

  #ifdef CGAL_EIGEN3_ENABLED
    for(int k=0; k<M.eigen_object().outerSize(); ++k) {
      for(typename Eigen::SparseMatrix<double>::InnerIterator
                                            it(M.eigen_object(), k); it; ++it) {
        outM <<  it.row() << " " << it.col() << " " << it.value() << '\n';
      }
    }
  #else // !CGAL_EIGEN3_ENABLED
    for(int k=0; k<M.row_dimension(); k++)
      for(int l=0; l<M.column_dimension(); l++)
        outM << k << " " << l << " " << M.get_coef(k, l) << '\n';
  #endif

    std::ofstream outB("matrices/B.txt");
    outB.precision(20);
    outB << total_size << std::endl;
    outB << B << std::endl;
#endif // CGAL_SMP_OUTPUT_ORBIFOLD_MATRICES

    // compute the flattening by solving the boundary conditions
    // while satisfying the convex combination property with L
    status = computeFlattening(mesh, M, B, uvmap, vimap);
    if(status != OK)
      return status;

    std::ofstream out("orbifold_result.off");
    IO::output_uvmap_to_off(mesh, bhd, uvmap, out);

    return OK;
  }

public:
  /// Constructor.
  Orbifold_Tutte_parameterizer_3(const Orbifold_type orb_type = Square,
                                 const Weight_type weight_type = Cotangent)
    :
      orb_type(orb_type),
      weight_type(weight_type)
  { }
};

} // namespace Surface_mesh_parameterization

} // namespace CGAL

#endif // CGAL_SURFACE_MESH_PARAMETERIZATION_ORBIFOLD_TUTTE_PARAMETERIZER_3_H
