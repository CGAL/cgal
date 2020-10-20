// Copyright (c) 2016  GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_SURFACE_MESH_PARAMETERIZATION_ARAP_PARAMETERIZER_3_H
#define CGAL_SURFACE_MESH_PARAMETERIZATION_ARAP_PARAMETERIZER_3_H

#include <CGAL/disable_warnings.h>

#include <CGAL/license/Surface_mesh_parameterization.h>

#include <CGAL/Surface_mesh_parameterization/internal/Bool_property_map.h>
#include <CGAL/Surface_mesh_parameterization/internal/Containers_filler.h>
#include <CGAL/Surface_mesh_parameterization/internal/kernel_traits.h>
#include <CGAL/Surface_mesh_parameterization/internal/validity.h>
#include <CGAL/Surface_mesh_parameterization/IO/File_off.h>

#include <CGAL/Surface_mesh_parameterization/Mean_value_coordinates_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/LSCM_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/MVC_post_processor_3.h>
#include <CGAL/Surface_mesh_parameterization/Two_vertices_parameterizer_3.h>

#include <CGAL/Surface_mesh_parameterization/parameterize.h>

#include <CGAL/Polygon_mesh_processing/connected_components.h>

#if defined(CGAL_EIGEN3_ENABLED)
#include <CGAL/Eigen_solver_traits.h>
#ifdef CGAL_SMP_USE_SPARSESUITE_SOLVERS
#include <Eigen/UmfPackSupport>
#endif
#endif

#include <CGAL/assertions.h>
#include <CGAL/basic.h>
#include <CGAL/circulator.h>
#include <CGAL/Default.h>
#include <CGAL/number_utils.h>

// Below are two macros that can be used to improve the accuracy of optimal Lt
// matrices.
// Note that at least one of these macros should be defined. If:
//
// - CGAL_SMP_SOLVE_CUBIC_EQUATION is defined: a cubic equation is solved instead of the
//   complete bivariate system. Although less accurate, it is usually sufficient.
// - CGAL_SMP_SOLVE_CUBIC_EQUATION and CGAL_SMP_SOLVE_EQUATIONS_WITH_GMP are defined:
//   the same cubic is solved but using GMP and CGAL's algebraic kernel.
// - CGAL_SMP_SOLVE_EQUATIONS_WITH_GMP is defined: a bivariate system is solved,
//   using GMP and CGAL's algebraic kernel.
//
// Using CGAL_SMP_SOLVE_EQUATIONS_WITH_GMP requires GMP, MPFI, and linking CGAL
// with Core and MPFI. This can be simply done in 'CMakeLists.txt' by using:
// 'find_package(CGAL QUIET COMPONENTS Core MPFI)'

// -----------------------------------------------------------------------------
#define CGAL_SMP_SOLVE_CUBIC_EQUATION
//#define CGAL_SMP_SOLVE_EQUATIONS_WITH_GMP
// -----------------------------------------------------------------------------

#if !defined(CGAL_SMP_SOLVE_CUBIC_EQUATION) && !defined(CGAL_SMP_SOLVE_EQUATIONS_WITH_GMP)
#error "Either 'CGAL_SMP_SOLVE_CUBIC_EQUATION' or 'CGAL_SMP_SOLVE_EQUATIONS_WITH_GMP' must be defined."
#endif

#if defined(CGAL_SMP_SOLVE_EQUATIONS_WITH_GMP)
#if !defined(CGAL_USE_GMP) || !defined(CGAL_USE_MPFI)
#error "'CGAL_SMP_SOLVE_EQUATIONS_WITH_GMP' cannot be defined if GMP or MPFI is not present."
#endif
#endif

#ifdef CGAL_SMP_SOLVE_EQUATIONS_WITH_GMP
#include <CGAL/GMP_arithmetic_kernel.h>
#include <CGAL/Algebraic_kernel_d_2.h>
#elif defined(CGAL_SMP_SOLVE_CUBIC_EQUATION)
#include <CGAL/Kernel/Conic_misc.h> // used to solve conic equations
#endif

#include <boost/function_output_iterator.hpp>
#include <boost/functional/hash.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/unordered_set.hpp>

#include <iostream>
#include <fstream>
#include <limits>
#include <set>
#include <string>
#include <utility>
#include <vector>

/// \file ARAP_parameterizer_3.h

// @todo Handle the case cot = 0 with a local parameterization aligned with the axes
//       (this produces C2=0 which is problematic to compute a & b)
// @todo Add distortion measures
// @todo Parallelize the local phase?

namespace CGAL {

namespace Surface_mesh_parameterization {

// ------------------------------------------------------------------------------------
// Declaration
// ------------------------------------------------------------------------------------

/// \ingroup PkgSurfaceMeshParameterizationMethods
///
/// The class `ARAP_parameterizer_3` implements the
/// *Local/Global Approach to Mesh Parameterization* \cgalCite{liu2008local}.
///
/// This parameterization allows the user to prioritize angle preservation,
/// shape preservation, or a balance of both.
/// A parameter \f$\lambda\f$ controls whether the priority is given to angle
/// or to shape preservation: when \f$\lambda=0\f$, the parameterization is
/// as-similar-as-possible (ASAP) and is equivalent to the (conforming) LSCM parameterization.
/// As \f$\lambda\f$ grows, the shape preservation becomes more and more important,
/// yielding, when \f$\lambda\f$ goes to infinity, a parameterization that is as-rigid-as-possible (ARAP).
///
/// This is a free border parameterization. There is no need to map the border of the surface
/// onto a convex polygon.
/// When \f$\lambda=0\f$, only two pinned vertices are needed to ensure a unique solution.
/// When \f$\lambda\f$ is non-null, the border does not need to be parameterized and
/// a random vertex is pinned.
///
/// If flips are present in the initial parameterization, a post-processing step
/// is applied using `CGAL::Surface_mesh_parameterization::MVC_post_processor_3<TriangleMesh_, SolverTraits_>`
/// to attempt to obtain a valid final embedding.
///
/// A one-to-one mapping is *not* guaranteed.
///
/// \cgalModels `Parameterizer_3`
///
/// \tparam TriangleMesh_ must be a model of `FaceGraph`.
///
/// \tparam BorderParameterizer_ is a Strategy to parameterize the surface border
///         and must be a model of `Parameterizer_3`.<br>
///         <b>%Default:</b>
/// \code
///   Two_vertices_parameterizer_3<TriangleMesh_>
/// \endcode
///
/// \tparam SolverTraits_ must be a model of `SparseLinearAlgebraTraits_d`.<br>
///         Note that the system is *not* symmetric.<br>
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
/// \sa `CGAL::Surface_mesh_parameterization::Fixed_border_parameterizer_3<TriangleMesh, BorderParameterizer, SolverTraits>`
/// \sa `CGAL::Iterative_authalic_parameterizer_3<TriangleMesh, BorderParameterizer, SolverTraits>`
///
template < class TriangleMesh_,
           class BorderParameterizer_ = Default,
           class SolverTraits_ = Default>
class ARAP_parameterizer_3
{
public:
#ifndef DOXYGEN_RUNNING
  typedef typename Default::Get<
    BorderParameterizer_,
    Two_vertices_parameterizer_3<TriangleMesh_> >::type       Border_parameterizer;

  #if !defined(CGAL_EIGEN3_ENABLED)
  CGAL_static_assertion_msg(!(boost::is_same<SolverTraits_, Default>::value),
                            "Error: You must either provide 'SolverTraits_' or link CGAL with the Eigen library");
  #endif

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
    SolverTraits_ // no parameter provided, and Eigen is not enabled: so don't compile!
  #endif
  >::type                                                     Solver_traits;
#else
  /// Border parameterizer type
  typedef Border_parameterizer_                               Border_parameterizer;

  /// Solver traits type
  typedef SolverTraits_                                       Solver_traits;

  /// Number type, deduced from the internal vertex point map of `Triangle_mesh`
  typedef unspecified_type                                    NT;
#endif

  /// Triangle mesh type
  typedef TriangleMesh_                                       Triangle_mesh;

  typedef TriangleMesh_                                       TriangleMesh;

  /// Mesh halfedge type
  typedef typename boost::graph_traits<Triangle_mesh>::halfedge_descriptor  halfedge_descriptor;

// Private types
private:
  typedef typename boost::graph_traits<Triangle_mesh>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::graph_traits<Triangle_mesh>::face_descriptor      face_descriptor;
  typedef typename boost::graph_traits<Triangle_mesh>::face_iterator        face_iterator;
  typedef typename boost::graph_traits<Triangle_mesh>::vertex_iterator      vertex_iterator;

  typedef CGAL::Halfedge_around_target_circulator<Triangle_mesh>    halfedge_around_target_circulator;
  typedef CGAL::Halfedge_around_face_circulator<Triangle_mesh>      halfedge_around_face_circulator;

  typedef boost::unordered_set<vertex_descriptor>                   Vertex_set;
  typedef std::vector<face_descriptor>                              Faces_vector;

  // Traits subtypes:
  typedef typename internal::Kernel_traits<Triangle_mesh>::Kernel   Kernel;
  typedef typename Kernel::FT                                       NT;
  typedef typename Kernel::Point_2                                  Point_2;
  typedef typename Kernel::Point_3                                  Point_3;
  typedef typename Kernel::Vector_2                                 Vector_2;
  typedef typename Kernel::Vector_3                                 Vector_3;

  // Solver traits subtypes:
  typedef typename Solver_traits::Vector                            Vector;
  typedef typename Solver_traits::Matrix                            Matrix;

  // Memory maps
    // Each triangle is associated a linear transformation matrix
  typedef std::pair<NT, NT>                                         Lt_matrix;
  typedef CGAL::Unique_hash_map<face_descriptor,
                                Lt_matrix,
                                boost::hash<face_descriptor> >      Lt_hash_map;
  typedef boost::associative_property_map<Lt_hash_map>              Lt_map;

    // Each angle (uniquely determined by the opposite half edge) has a cotangent
  typedef CGAL::Unique_hash_map<halfedge_descriptor, NT,
                                boost::hash<halfedge_descriptor> >  Cot_hm;
  typedef boost::associative_property_map<Cot_hm>                   Cot_map;

    // Each face has a local 2D isometric parameterization
  typedef std::pair<int, int>                                       Local_indices;
  typedef CGAL::Unique_hash_map<halfedge_descriptor,
                                Local_indices,
                                boost::hash<halfedge_descriptor> >  Lp_hm;
  typedef boost::associative_property_map<Lp_hm>                    Lp_map;
  typedef std::vector<Point_2>                                      Local_points;

// Private fields
private:
  // %Object that maps (at least two) border vertices onto a 2D space
  Border_parameterizer m_borderParameterizer;

  // Traits object to solve a sparse linear system
  Solver_traits m_linearAlgebra;

  // Controlling parameters
  const NT m_lambda;
  const NT m_lambda_tolerance; // used to determine when we switch to full ARAP

  // Energy minimization parameters
  const unsigned int m_iterations;
  const NT m_tolerance;

// Private accessors
private:
  // Get the object that maps the surface's border onto a 2D space.
  Border_parameterizer& get_border_parameterizer() { return m_borderParameterizer; }

  // Get the sparse linear algebra (traits object to access the linear system).
  Solver_traits& get_linear_algebra_traits() { return m_linearAlgebra; }

// Private utilities
private:
  // Print the parameterized mesh.
  template <typename VertexUVMap,
            typename VertexIndexMap>
  void output_uvmap(const std::string filename,
                    const Triangle_mesh& mesh,
                    const Vertex_set& vertices,
                    const Faces_vector& faces,
                    const VertexUVMap uvmap,
                    const VertexIndexMap vimap) const
  {
    std::ofstream out(filename.c_str());
    IO::output_uvmap_to_off(mesh, vertices, faces, uvmap, vimap, out);
  }

  // Print the parameterized mesh.
  template <typename VertexUVMap,
            typename VertexIndexMap>
  void output_uvmap(const std::string filename,
                    const unsigned int iter,
                    const Triangle_mesh& mesh,
                    const Vertex_set& vertices,
                    const Faces_vector& faces,
                    const VertexUVMap uvmap,
                    const VertexIndexMap vimap) const
  {
    std::ostringstream out_ss;
    out_ss << filename << iter << ".off" << std::ends;
    std::ofstream out(out_ss.str().c_str());
    IO::output_uvmap_to_off(mesh, vertices, faces, uvmap, vimap, out);
  }

  // Copy the data from two vectors to the UVmap.
  template <typename VertexUVMap,
            typename VertexIndexMap,
            typename VertexParameterizedMap>
  void assign_solution(const Vector& Xu,
                       const Vector& Xv,
                       const Vertex_set& vertices,
                       VertexUVMap uvmap,
                       const VertexIndexMap vimap,
                       const VertexParameterizedMap vpmap)
  {
    for(vertex_descriptor vd : vertices) {
      // The solver might not have managed to exactly constrain the vertex that was marked
      // as constrained; simply don't update its position.
      if(get(vpmap, vd))
        continue;

      int index = get(vimap, vd);
      NT u = Xu(index);
      NT v = Xv(index);
      put(uvmap, vd, Point_2(u, v));
    }
  }

// Private operations
private:
  // Store the vertices and faces of the mesh in memory.
  Error_code initialize_containers(const Triangle_mesh& mesh,
                                   halfedge_descriptor bhd,
                                   Vertex_set& vertices,
                                   Faces_vector& faces) const
  {
    internal::Containers_filler<Triangle_mesh> fc(mesh, vertices, &faces);
    Polygon_mesh_processing::connected_component(
                                      face(opposite(bhd, mesh), mesh),
                                      mesh,
                                      boost::make_function_output_iterator(fc));

    if (vertices.empty() || faces.empty())
      return ERROR_EMPTY_MESH;

    return OK;
  }

  // Initialize the UV values with a first parameterization of the input.
  template <typename VertexUVMap,
            typename VertexIndexMap>
  Error_code compute_initial_uv_map(Triangle_mesh& mesh,
                                    halfedge_descriptor bhd,
                                    VertexUVMap uvmap,
                                    VertexIndexMap vimap) const
  {
    Error_code status;

    unsigned int number_of_borders =
                         CGAL::Polygon_mesh_processing::number_of_borders(mesh);
    if(number_of_borders == 0) {
      status = ERROR_BORDER_TOO_SHORT;
      return status;
    }

    // temporary vpmap since we do not need it in the future
    boost::unordered_set<vertex_descriptor> vs;
    internal::Bool_property_map<boost::unordered_set<vertex_descriptor> > vpmap(vs);

    // According to the paper, MVC is better for single border and LSCM is better
    // when there are multiple borders
    if(number_of_borders == 1) {
      typedef Mean_value_coordinates_parameterizer_3<Triangle_mesh>     MVC_parameterizer;
      MVC_parameterizer mvc_parameterizer;
      status = mvc_parameterizer.parameterize(mesh, bhd, uvmap, vimap, vpmap);
    } else {
      typedef LSCM_parameterizer_3<Triangle_mesh, Border_parameterizer> LSCM_parameterizer;
      LSCM_parameterizer lscm_parameterizer;
      status = lscm_parameterizer.parameterize(mesh, bhd, uvmap, vimap, vpmap);
    }

    return status;
  }

  // Parameterize the border. The number of fixed vertices depends on the value
  // of the parameter lambda.
  template <typename VertexUVMap,
            typename VertexIndexMap,
            typename VertexParameterizedMap>
  Error_code parameterize_border(const Triangle_mesh& mesh,
                                 const Vertex_set& vertices,
                                 halfedge_descriptor bhd,
                                 VertexIndexMap vimap,
                                 VertexParameterizedMap vpmap)
  {
    Error_code status = OK;
    CGAL_precondition(!vertices.empty());

    if(m_lambda != 0.) {
      // Fix a random vertex, the value in uvmap is already set
      vertex_descriptor vd = *(vertices.begin());
      put(vpmap, vd, true);
    } else {
      // This fixes two vertices that are far in the original geometry

      // A local uvmap (that is then discarded) is passed to avoid changing
      // the values of the real uvmap. Since we can't get VertexUVMap::C,
      // we build a map with the same key and value types
      typedef boost::unordered_map<typename VertexUVMap::key_type,
                                   typename VertexUVMap::value_type> Useless_map;
      typedef boost::associative_property_map<Useless_map>           Useless_pmap;

      Useless_map useless_map;
      Useless_pmap useless_uvpmap(useless_map);
      status = get_border_parameterizer().parameterize(mesh, bhd, useless_uvpmap,
                                                       vimap, vpmap);
    }

    return status;
  }

  // Compute the cotangent of the angle between the vectors ij and ik.
  void compute_cotangent_angle(const Triangle_mesh& mesh,
                               halfedge_descriptor hd,
                               vertex_descriptor vi,
                               vertex_descriptor vj,
                               vertex_descriptor vk,
                               Cot_map ctmap) const
  {
    typedef typename boost::property_map<Triangle_mesh,
                                         boost::vertex_point_t>::const_type PPmap;
    const PPmap ppmap = get(vertex_point, mesh);

    const Point_3& position_vi = get(ppmap, vi);
    const Point_3& position_vj = get(ppmap, vj);
    const Point_3& position_vk = get(ppmap, vk);

    NT cot = internal::cotangent<Kernel>(position_vi, position_vj, position_vk);
    put(ctmap, hd, cot);
  }

  // Fill the map 'ctmap' with the cotangents of the angles of the faces of 'mesh'.
  Error_code compute_cotangent_angles(const Triangle_mesh& mesh,
                                      const Faces_vector& faces,
                                      Cot_map ctmap) const
  {
    for(face_descriptor fd : faces) {
      halfedge_descriptor hd = halfedge(fd, mesh), hdb = hd;

      vertex_descriptor vi = target(hd, mesh);
      hd = next(hd, mesh);
      vertex_descriptor vj = target(hd, mesh);
      hd = next(hd, mesh);
      vertex_descriptor vk = target(hd, mesh);
      hd = next(hd, mesh);

      if(hd != hdb) { // make sure that it is a triangular face
        return ERROR_NON_TRIANGULAR_MESH;
      }

      compute_cotangent_angle(mesh, hd, vk, vj, vi , ctmap); // angle at v_j
      compute_cotangent_angle(mesh, next(hd, mesh), vi, vk, vj , ctmap); // angle at v_k
      compute_cotangent_angle(mesh, prev(hd, mesh), vj, vi, vk , ctmap); // angle at v_i
    }

    return OK;
  }

  // computes `w_ij`, the `(i,j)`-coefficient of matrix `A` for `j` neighbor vertex of `i`.
  NT compute_w_ij(const Triangle_mesh& mesh,
                  halfedge_descriptor hd,
                  const Cot_map ctmap) const
  {
    // Note that BGL halfedge and vertex circulators move clockwise!!

    // coefficient corresponding to the angle at vk if vk is the vertex before vj
    // while circulating around vi
    NT c_k = get(ctmap, opposite(hd, mesh));

    // coefficient corresponding to the angle at vl if vl is the vertex after vj
    // while circulating around vi
    NT c_l = get(ctmap, hd);

    NT weight = c_k + c_l;
    return weight;
  }

  // Compute the line i of matrix A
  // - call compute_w_ij() to compute the A coefficient w_ij for each neighbor v_j.
  // - compute w_ii = - sum of w_ijs.
  //
  // \pre Vertices must be indexed.
  // \pre Vertex i musn't be already parameterized.
  // \pre Line i of A must contain only zeros.
  template <typename VertexIndexMap>
  Error_code fill_linear_system_matrix(Matrix& A,
                                       const Triangle_mesh& mesh,
                                       vertex_descriptor vertex,
                                       const Cot_map ct_map,
                                       VertexIndexMap vimap) const
  {
    int i = get(vimap, vertex);

    // Circulate over the vertices around 'vertex' to compute w_ii and w_ijs
    NT w_ii = 0;
    int vertexIndex = 0;

    halfedge_around_target_circulator hc(vertex, mesh), end = hc;
    CGAL_For_all(hc, end) {
      halfedge_descriptor hd = *hc;
      CGAL_assertion(target(hd, mesh) == vertex);

      NT w_ij = -1.0 * compute_w_ij(mesh, hd, ct_map);

      // w_ii = - sum of w_ijs
      w_ii -= w_ij;

      // Get j index
      vertex_descriptor v_j = source(hd, mesh);
      int j = get(vimap, v_j);

      // Set w_ij in matrix
      A.set_coef(i, j, w_ij, true /*new*/);
      vertexIndex++;
    }

    if (vertexIndex < 2)
      return ERROR_NON_TRIANGULAR_MESH;

    // Set w_ii in matrix
    A.set_coef(i, i, w_ii, true /*new*/);
    return OK;
  }

  // Initialize the (constant) matrix A in the linear system "A*X = B",
  // after (at least two) border vertices parameterization.
  template <typename VertexIndexMap,
            typename VertexParameterizedMap>
  Error_code initialize_matrix_A(const Triangle_mesh& mesh,
                                 const Vertex_set& vertices,
                                 const Cot_map ctmap,
                                 VertexIndexMap vimap,
                                 VertexParameterizedMap vpmap,
                                 Matrix& A) const
  {
    Error_code status = OK;

    // compute A
    unsigned int count = 0;
    for(vertex_descriptor vd : vertices) {
      if(!get(vpmap, vd)) { // not yet parameterized
        // Compute the line i of the matrix A
        status = fill_linear_system_matrix(A, mesh, vd, ctmap, vimap);
        if(status != OK)
          return status;
      } else { // fixed vertices
        int index = get(vimap, vd);
        A.set_coef(index, index, 1, true /*new*/);
        ++count;
      }
    }
    return status;
  }

  // Solves the cubic equation a3 x^3 + a2 x^2 + a1 x + a0 = 0.
  std::size_t solve_cubic_equation(const NT a3, const NT a2, const NT a1, const NT a0,
                                   std::vector<NT>& roots) const
  {
    CGAL_precondition(roots.empty());
    NT r1 = 0, r2 = 0, r3 = 0; // roots of the cubic equation

    int root_n = CGAL::solve_cubic(a3, a2, a1, a0, r1, r2, r3);
    CGAL_postcondition(root_n > 0);

    roots.push_back(r1);
    if(root_n > 1)
      roots.push_back(r2);
    if(root_n > 2)
      roots.push_back(r3);

    return roots.size();
  }

#if defined(CGAL_SMP_SOLVE_CUBIC_EQUATION) && defined(CGAL_SMP_SOLVE_EQUATIONS_WITH_GMP)
  // Solves the equation a3 x^3 + a2 x^2 + a1 x + a0 = 0, using CGAL's algebraic kernel.
  std::size_t solve_cubic_equation_with_AK(const NT a3, const NT a2,
                                           const NT a1, const NT a0,
                                           std::vector<NT>& roots) const
  {
    CGAL_precondition(roots.empty());

    typedef CGAL::Gmpq                                        GMP_NT;
    typedef CGAL::Algebraic_kernel_d_1<GMP_NT>                Algebraic_kernel_d_1;
    typedef typename Algebraic_kernel_d_1::Polynomial_1       Polynomial_1;
    typedef typename Algebraic_kernel_d_1::Algebraic_real_1   Algebraic_real_1;
    typedef typename Algebraic_kernel_d_1::Coefficient        Coefficient;
    typedef typename Algebraic_kernel_d_1::Multiplicity_type  Multiplicity_type;

    typedef CGAL::Polynomial_traits_d<Polynomial_1>           Polynomial_traits_1;

    typedef Algebraic_kernel_d_1::Solve_1                     Solve_1;

    Algebraic_kernel_d_1 ak_1;
    const Solve_1 solve_1 = ak_1.solve_1_object();

    GMP_NT a3q(a3);
    GMP_NT a2q(a2);
    GMP_NT a1q(a1);
    GMP_NT a0q(a0);

    typename Polynomial_traits_1::Construct_polynomial construct_polynomial_1;
    std::pair<CGAL::Exponent_vector, Coefficient> coeffs_x[1]
      = {std::make_pair(CGAL::Exponent_vector(1,0), Coefficient(1))};
    Polynomial_1 x=construct_polynomial_1(coeffs_x, coeffs_x+1);
    Polynomial_1 pol = a3q*CGAL::ipower(x,3) + a2q*CGAL::ipower(x,2) + a1q*x + a0q;

    std::vector<std::pair<Algebraic_real_1, Multiplicity_type> > ak_roots;
    solve_1(pol, std::back_inserter(ak_roots));

    for(std::size_t i=0; i<ak_roots.size(); ++i)
      roots.push_back(CGAL::to_double(ak_roots[i].first));

    return roots.size();
  }
#endif // defined(CGAL_SMP_SOLVE_CUBIC_EQUATION) && defined(CGAL_SMP_SOLVE_EQUATIONS_WITH_GMP)

#if !defined(CGAL_SMP_SOLVE_CUBIC_EQUATION) && defined(CGAL_SMP_SOLVE_EQUATIONS_WITH_GMP)
  // Solve the bivariate system
  // { C1 * a + 2 * lambda * a ( a^2 + b^2 - 1 ) = C2
  // { C1 * b + 2 * lambda * b ( a^2 + b^2 - 1 ) = C3
  // using CGAL's algebraic kernel.
  std::size_t solve_bivariate_system(const NT C1, const NT C2, const NT C3,
                                     std::vector<NT>& a_roots, std::vector<NT>& b_roots) const
  {
    typedef CGAL::Gmpq                                        GMP_NT;
    typedef CGAL::Algebraic_kernel_d_2<GMP_NT>                Algebraic_kernel_d_2;
    typedef typename Algebraic_kernel_d_2::Polynomial_2       Polynomial_2;
    typedef typename Algebraic_kernel_d_2::Algebraic_real_2   Algebraic_real_2;
    typedef typename Algebraic_kernel_d_2::Multiplicity_type  Multiplicity_type;
    typedef typename Algebraic_kernel_d_2::Coefficient        Coefficient;

    typedef CGAL::Polynomial_traits_d<Polynomial_2>           Polynomial_traits_2;

    typedef Algebraic_kernel_d_2::Solve_2             Solve_2;
    typedef Algebraic_kernel_d_2::Is_coprime_2        Is_coprime_2;
    typedef Algebraic_kernel_d_2::Make_coprime_2      Make_coprime_2;

    Algebraic_kernel_d_2 ak_2;
    const Is_coprime_2 is_coprime_2 = ak_2.is_coprime_2_object();
    const Solve_2 solve_2 = ak_2.solve_2_object();
    const Make_coprime_2 make_coprime_2 = ak_2.make_coprime_2_object();

    typename Polynomial_traits_2::Construct_polynomial construct_polynomial_2;
    std::pair<CGAL::Exponent_vector, Coefficient> coeffs_x[1]
      = {std::make_pair(CGAL::Exponent_vector(1,0), Coefficient(1))};
    Polynomial_2 x=construct_polynomial_2(coeffs_x,coeffs_x+1);
    std::pair<CGAL::Exponent_vector, Coefficient> coeffs_y[1]
      = {std::make_pair(CGAL::Exponent_vector(0,1), Coefficient(1))};
    Polynomial_2 y=construct_polynomial_2(coeffs_y,coeffs_y+1);

    GMP_NT C1q(C1);
    GMP_NT C2q(C2);
    GMP_NT C3q(C3);

    Polynomial_2 pol1 = C1q * x + 2 * m_lambda * x * ( x*x + y*y - 1) - C2q;
    Polynomial_2 pol2 = C1q * y + 2 * m_lambda * y * ( x*x + y*y - 1) - C3q;

    CGAL_precondition_code(typedef Algebraic_kernel_d_2::Is_square_free_2    Is_square_free_2;)
    CGAL_precondition_code(const Is_square_free_2 is_square_free_2 = ak_2.is_square_free_2_object();)
    CGAL_precondition(is_square_free_2(pol1));
    CGAL_precondition(is_square_free_2(pol2));
    if(!is_coprime_2(pol1, pol2)) {
      CGAL_assertion(false); // @todo handle that case

      Polynomial_2 g, q1, q2;
      make_coprime_2(pol1, pol2, g, q1, q2);
      pol1 = q1;
      pol2 = q2;
    }

    std::vector<std::pair<Algebraic_real_2, Multiplicity_type> > roots;
    solve_2(pol1, pol2, std::back_inserter(roots));

    for(std::size_t i=0; i<roots.size(); i++)
    {
      a_roots.push_back(roots[i].first.to_double().first);
      b_roots.push_back(roots[i].first.to_double().second);
    }

    return a_roots.size();
  }
#endif // !defined(CGAL_SMP_SOLVE_CUBIC_EQUATION) && defined(CGAL_SMP_SOLVE_EQUATIONS_WITH_GMP)

  // Compute the root that gives the lowest face energy.
  template <typename VertexUVMap>
  std::size_t compute_root_with_lowest_energy(const Triangle_mesh& mesh,
                                              face_descriptor fd,
                                              const Cot_map ctmap,
                                              const Local_points& lp,
                                              const Lp_map lpmap,
                                              const VertexUVMap uvmap,
                                              const NT C2_denom, const NT C3,
                                              const std::vector<NT>& roots) const
  {
    CGAL_precondition(!roots.empty());
    NT E_min = (std::numeric_limits<NT>::max)();
    std::size_t index_arg = 0;
    for(std::size_t i=0; i<roots.size(); ++i)
    {
      const NT a = roots[i];
      const NT b = C3 * C2_denom * a;
      NT Ef = compute_current_face_energy(mesh, fd, ctmap, lp, lpmap, uvmap, a, b);
      if(Ef < E_min) {
        E_min = Ef;
        index_arg = i;
      }
    }
    return index_arg;
  }

  // Compute the root that gives the lowest face energy.
  template <typename VertexUVMap>
  std::size_t compute_root_with_lowest_energy(const Triangle_mesh& mesh,
                                              face_descriptor fd,
                                              const Cot_map ctmap,
                                              const Local_points& lp,
                                              const Lp_map lpmap,
                                              const VertexUVMap uvmap,
                                              const std::vector<NT>& a_roots,
                                              const std::vector<NT>& b_roots) const
  {
    CGAL_precondition(!a_roots.empty() && a_roots.size() == b_roots.size());
    NT E_min = (std::numeric_limits<NT>::max)();
    std::size_t index_arg = -1;
    for(std::size_t i=0; i<a_roots.size(); ++i)
    {
      NT Ef = compute_current_face_energy(mesh, fd, ctmap, lp, lpmap, uvmap,
                                          a_roots[i], b_roots[i]);
      if(Ef < E_min) {
        E_min = Ef;
        index_arg = i;
      }
    }
    return index_arg;
  }

  // Compute the optimal values of the linear transformation matrices Lt.
  template <typename VertexUVMap>
  Error_code compute_optimal_Lt_matrices(const Triangle_mesh& mesh,
                                         const Faces_vector& faces,
                                         const Cot_map ctmap,
                                         const Local_points& lp,
                                         const Lp_map lpmap,
                                         const VertexUVMap uvmap,
                                         Lt_map ltmap) const
  {
    Error_code status = OK;

    for(face_descriptor fd : faces) {
      // Compute the coefficients C1, C2, C3
      NT C1 = 0., C2 = 0., C3 = 0.;

      halfedge_around_face_circulator hc(halfedge(fd, mesh), mesh), end(hc);
      CGAL_For_all(hc, end) {
        halfedge_descriptor hd = *hc;
        NT c = get(ctmap, hd);

        // UV positions
        const Point_2& uvpi = get(uvmap, source(hd, mesh));
        const Point_2& uvpj = get(uvmap, target(hd, mesh));
        NT diff_x = uvpi.x() - uvpj.x();
        NT diff_y = uvpi.y() - uvpj.y();
//        CGAL_warning(diff_x == 0. && diff_y == 0.);

        // local positions (in the isometric 2D param)
        const Local_indices& li = get(lpmap, hd);
        const Point_2& ppi = lp[ li.first ];
        const Point_2& ppj = lp[ li.second ];
        NT p_diff_x = ppi.x() - ppj.x();
        NT p_diff_y = ppi.y() - ppj.y();
        CGAL_precondition(p_diff_x != 0. || p_diff_y != 0.);

        C1 += c * ( p_diff_x*p_diff_x + p_diff_y*p_diff_y );
        C2 += c * ( diff_x*p_diff_x + diff_y*p_diff_y );
        C3 += c * ( diff_x*p_diff_y - diff_y*p_diff_x );
      }

      // Compute a and b
      NT a = 0., b = 0.;

      if(m_lambda == 0.) { // ASAP
        CGAL_precondition(C1 != 0.);
        a = C2 / C1;
        b = C3 / C1;
      }
      else if( std::abs(C1) < m_lambda_tolerance * m_lambda &&
               std::abs(C2) < m_lambda_tolerance * m_lambda ) { // ARAP
        // If lambda is large compared to C1 and C2, the cubic equation that
        // determines a and b can be simplified to a simple quadric equation

        CGAL_precondition(C2*C2 + C3*C3 != 0.);
        NT denom = 1. / CGAL::sqrt(C2*C2 + C3*C3);
        a = C2 * denom;
        b = C3 * denom;
      }
      else { // general case
#ifdef CGAL_SMP_SOLVE_CUBIC_EQUATION
        CGAL_precondition(C2 != 0.);
        NT C2_denom = 1. / C2;
        NT a3_coeff = 2. * m_lambda * (C2 * C2 + C3 * C3) * C2_denom * C2_denom;

        std::vector<NT> roots;
#ifdef CGAL_SMP_SOLVE_EQUATIONS_WITH_GMP
        solve_cubic_equation_with_AK(a3_coeff, 0., (C1 - 2. * m_lambda), -C2, roots);
#else // !CGAL_SMP_SOLVE_EQUATIONS_WITH_GMP
        solve_cubic_equation(a3_coeff, 0., (C1 - 2. * m_lambda), -C2, roots);
#endif
        std::size_t ind = compute_root_with_lowest_energy(mesh, fd,
                                                          ctmap, lp, lpmap, uvmap,
                                                          C2_denom, C3, roots);

        a = roots[ind];
        b = C3 * C2_denom * a;
#else // !CGAL_SMP_SOLVE_CUBIC_EQUATION, solve the bivariate system
        std::vector<NT> a_roots;
        std::vector<NT> b_roots;
        solve_bivariate_system(C1, C2, C3, a_roots, b_roots);

        std::size_t ind = compute_root_with_lowest_energy(mesh, fd,
                                                          ctmap, lp, lpmap, uvmap,
                                                          a_roots, b_roots);
        a = a_roots[ind];
        b = b_roots[ind];
#endif
      }

      // Update the map faces --> optimal Lt matrices
      Lt_matrix ltm = std::make_pair(a, b);
      put(ltmap, fd, ltm);
    }

    return status;
  }

  // Computes the coordinates of the vertices p0, p1, p2
  // in a local 2D orthonormal basis of the triangle's plane.
  void project_triangle(const Point_3& p0, const Point_3& p1, const Point_3& p2,  // in
                        Point_2& z0, Point_2& z1, Point_2& z2) const              // out
  {
    Vector_3 X = p1 - p0;
    NT X_norm = CGAL::sqrt(X * X);
    if(X_norm != 0.0)
      X = X / X_norm;

    Vector_3 Z = CGAL::cross_product(X, p2 - p0);
    NT Z_norm = CGAL::sqrt(Z * Z);
    if(Z_norm != 0.0)
      Z = Z / Z_norm;

    Vector_3 Y = CGAL::cross_product(Z, X);

    NT x0 = 0;
    NT y0 = 0;
    NT x1 = CGAL::sqrt( (p1 - p0)*(p1 - p0) );
    NT y1 = 0;
    NT x2 = (p2 - p0) * X;
    NT y2 = (p2 - p0) * Y;

    z0 = Point_2(x0, y0);
    z1 = Point_2(x1, y1);
    z2 = Point_2(x2, y2);
  }

  // Compute the local parameterization (2D) of a face and store it in memory.
  void project_face(const Triangle_mesh& mesh,
                    vertex_descriptor vi,
                    vertex_descriptor vj,
                    vertex_descriptor vk,
                    Local_points& lp) const
  {
    typedef typename boost::property_map<Triangle_mesh,
                                         boost::vertex_point_t>::const_type PPmap;
    const PPmap ppmap = get(vertex_point, mesh);

    const Point_3& position_vi = get(ppmap, vi);
    const Point_3& position_vj = get(ppmap, vj);
    const Point_3& position_vk = get(ppmap, vk);

    Point_2 pvi, pvj, pvk;
    project_triangle(position_vi, position_vj, position_vk,
                     pvi, pvj, pvk);

    // Add the newly computed 2D points to the vector of local positions
    lp.push_back(pvi);
    lp.push_back(pvj);
    lp.push_back(pvk);
  }

  // Utility for fill_linear_system_rhs():
  // Compute the local isometric parameterization (2D) of the faces of the mesh.
  Error_code compute_local_parameterization(const Triangle_mesh& mesh,
                                            const Faces_vector& faces,
                                            Local_points& lp,
                                            Lp_map lpmap) const
  {
    int global_index = 0;

    for(face_descriptor fd : faces) {
      halfedge_descriptor hd = halfedge(fd, mesh), hdb = hd;

      vertex_descriptor vi = target(hd, mesh); // hd is k -- > i
      put(lpmap, hd, std::make_pair(global_index + 2, global_index));

      hd = next(hd, mesh);
      vertex_descriptor vj = target(hd, mesh); // hd is i -- > j
      put(lpmap, hd, std::make_pair(global_index, global_index + 1));

      hd = next(hd, mesh);
      vertex_descriptor vk = target(hd, mesh); // hd is j -- > k
      put(lpmap, hd, std::make_pair(global_index + 1, global_index + 2));

      hd = next(hd, mesh);
      if(hd != hdb) { // to make sure that it is a triangular face
        return ERROR_NON_TRIANGULAR_MESH;
      }

      project_face(mesh, vi, vj, vk, lp);
      global_index += 3;
    }

    if(lp.size() != 3 * faces.size())
      return ERROR_NON_TRIANGULAR_MESH;

    return OK;
  }

  // Compute the coefficient b_ij = (i,j) of the right hand side vector B,
  // for j neighbor vertex of i.
  void compute_b_ij(const Triangle_mesh& mesh,
                    halfedge_descriptor hd,
                    const Cot_map ctmap,
                    const Local_points& lp,
                    const Lp_map lpmap,
                    const Lt_map ltmap,
                    NT& x, NT& y) const // x for Bu and y for Bv
  {
    CGAL_precondition(x == 0.0 && y == 0.0);

    // Note that :
    // -- Circulators move in a clockwise manner
    // -- The halfedge hd points to vi

    // -- Face IJK with vk before vj while circulating around vi
    halfedge_descriptor hd_opp = opposite(hd, mesh); // hd_opp points to vj
    face_descriptor fd_k = face(hd_opp, mesh);

    // Get the matrix L_t corresponding to the face ijk
    const Lt_matrix& ltm_k = get(ltmap, fd_k);
    NT a_k = ltm_k.first;
    NT b_k = ltm_k.second;

    // Get the local parameterization in the triangle ijk
    const Local_indices& loc_indices_k = get(lpmap, hd_opp);
    const Point_2& pvi_k = lp[ loc_indices_k.first ];
    const Point_2& pvj_k = lp[ loc_indices_k.second ];

    // x_i - x_j in the local parameterization of the triangle ijk
    NT diff_k_x = pvi_k.x() - pvj_k.x();
    NT diff_k_y = pvi_k.y() - pvj_k.y();

    // get the cotangent angle at vk
    NT ct_k = get(ctmap, hd_opp);

    // cot_k * Lt_k * (xi - xj)_k
    x = ct_k * (  a_k * diff_k_x + b_k * diff_k_y );
    y = ct_k * ( -b_k * diff_k_x + a_k * diff_k_y );

    // -- Face IJL with vl after vj while circulating around vi
    face_descriptor fd_l = face(hd, mesh);

    // Get the matrix L_t corresponding to the face ijl
    const Lt_matrix& ltm_l = get(ltmap, fd_l);
    NT a_l = ltm_l.first;
    NT b_l = ltm_l.second;

    // Get the local parameterization in the triangle ijl
    const Local_indices& loc_indices_l = get(lpmap, hd);
    const Point_2& pvi_l = lp[ loc_indices_l.second ]; // because hd points to vi
    const Point_2& pvj_l = lp[ loc_indices_l.first ];

    // x_i - x_j in the local parameterization of the triangle ijl
    NT diff_l_x = pvi_l.x() - pvj_l.x();
    NT diff_l_y = pvi_l.y() - pvj_l.y();

    // get the cotangent angle at vl
    NT ct_l = get(ctmap, hd);

    // cot_l * Lt_l * (xi - xj)_l
    x += ct_l * (  a_l * diff_l_x + b_l * diff_l_y );
    y += ct_l * ( -b_l * diff_l_x + a_l * diff_l_y );
  }

  // Compute the line i of right hand side vectors Bu and Bv
  // - call compute_b_ij() for each neighbor v_j to compute the B coefficient b_i
  //
  // \pre Vertices must be indexed.
  // \pre Vertex i musn't be already parameterized.
  // \pre Lines i of Bu and Bv must be zero.
  template <typename VertexIndexMap>
  Error_code fill_linear_system_rhs(const Triangle_mesh& mesh,
                                    vertex_descriptor vertex,
                                    const Cot_map ctmap,
                                    const Local_points& lp,
                                    const Lp_map lpmap,
                                    const Lt_map ltmap,
                                    VertexIndexMap vimap,
                                    Vector& Bu, Vector& Bv) const
  {
    int i = get(vimap, vertex);

    // Circulate over vertices around 'vertex' to compute w_ii and w_ijs
    NT bu_i = 0, bv_i = 0;
    int vertexIndex = 0;

    halfedge_around_target_circulator hc(vertex, mesh), end = hc;
    CGAL_For_all(hc, end) {
      halfedge_descriptor hd = *hc;
      CGAL_assertion(target(hd, mesh) == vertex);

      NT x = 0., y = 0.;
      compute_b_ij(mesh, hd, ctmap, lp, lpmap, ltmap, x, y);

      bu_i += x;
      bv_i += y;
      vertexIndex++;
    }

    // Set the entries for the right hand side vectors
    Bu.set(i, bu_i);
    Bv.set(i, bv_i);

    if (vertexIndex < 2)
      return ERROR_NON_TRIANGULAR_MESH;

    return OK;
  }

  // Compute the entries of the right hand side of the ARAP linear system.
  //
  // \pre Vertices must be indexed.
  template <typename VertexUVMap,
            typename VertexIndexMap,
            typename VertexParameterizedMap>
  Error_code compute_rhs(const Triangle_mesh& mesh,
                         const Vertex_set& vertices,
                         const Cot_map ctmap,
                         const Local_points& lp,
                         const Lp_map lpmap,
                         const Lt_map ltmap,
                         VertexUVMap uvmap,
                         VertexIndexMap vimap,
                         VertexParameterizedMap vpmap,
                         Vector& Bu, Vector& Bv) const
  {
    // Initialize the right hand side B in the linear system "A*X = B"
    Error_code status = OK;

    unsigned int count = 0;
    for(vertex_descriptor vd : vertices) {
      if(!get(vpmap, vd)) { // not yet parameterized
        // Compute the lines i of the vectors Bu and Bv
        status = fill_linear_system_rhs(mesh, vd, ctmap, lp, lpmap,
                                                  ltmap, vimap, Bu, Bv);
        if(status != OK)
          return status;
      } else { // fixed vertices
        int index = get(vimap, vd);
        const Point_2& uv = get(uvmap, vd);
        Bu.set(index, uv.x());
        Bv.set(index, uv.y());
        ++count;
      }
    }
    return status;
  }

  // Compute the right hand side and solve the linear system to obtain the
  // new UV coordinates.
  template <typename VertexUVMap,
            typename VertexIndexMap,
            typename VertexParameterizedMap>
  Error_code update_solution(const Triangle_mesh& mesh,
                             const Vertex_set& vertices,
                             const Cot_map ctmap,
                             const Local_points& lp,
                             const Lp_map lpmap,
                             const Lt_map ltmap,
                             VertexUVMap uvmap,
                             VertexIndexMap vimap,
                             VertexParameterizedMap vpmap,
                             const Matrix& A)
  {
    Error_code status = OK;

    // Create two sparse linear systems "A*Xu = Bu" and "A*Xv = Bv"
    int nbVertices = static_cast<int>(vertices.size());
    Vector Xu(nbVertices), Xv(nbVertices), Bu(nbVertices), Bv(nbVertices);

    // Fill the right hand side vectors
    status = compute_rhs(mesh, vertices, ctmap, lp, lpmap, ltmap,
                                         uvmap, vimap, vpmap, Bu, Bv);
    if (status != OK)
      return status;

    // Solve "A*Xu = Bu". On success, the solution is (1/Du) * Xu.
    // Solve "A*Xv = Bv". On success, the solution is (1/Dv) * Xv.
    NT Du, Dv;
    if(!get_linear_algebra_traits().linear_solver(A, Bu, Xu, Du) ||
       !get_linear_algebra_traits().linear_solver(A, Bv, Xv, Dv)) {
      std::cerr << "Could not solve linear system" << std::endl;
      status = ERROR_CANNOT_SOLVE_LINEAR_SYSTEM;
      return status;
    }

    // WARNING: this package does not support homogeneous coordinates!
    CGAL_assertion(Du == 1.0);
    CGAL_assertion(Dv == 1.0);

    CGAL_postcondition_code
    (
      // make sure that the constrained vertices have not been moved
      for(vertex_descriptor vd : vertices) {
        if(get(vpmap, vd)) {
          int index = get(vimap, vd);
          CGAL_warning(std::abs(Xu[index] - Bu[index] ) < 1e-7);
          CGAL_warning(std::abs(Xv[index] - Bv[index] ) < 1e-7);
        }
      }
    )

    assign_solution(Xu, Xv, vertices, uvmap, vimap, vpmap);
    return status;
  }


  // Compute the current energy of a face, given a linear transformation matrix.
  template <typename VertexUVMap>
  NT compute_current_face_energy(const Triangle_mesh& mesh,
                                 face_descriptor fd,
                                 const Cot_map ctmap,
                                 const Local_points& lp,
                                 const Lp_map lpmap,
                                 const VertexUVMap uvmap,
                                 const NT a, const NT b) const
  {
    NT Ef = 0.;

    halfedge_around_face_circulator hc(halfedge(fd, mesh), mesh), end(hc);
    CGAL_For_all(hc, end) {
      halfedge_descriptor hd = *hc;
      NT cot = get(ctmap, hd);
      NT nabla_x = 0., nabla_y = 0.;

      // UV positions
      const Point_2& pi = get(uvmap, source(hd, mesh));
      const Point_2& pj = get(uvmap, target(hd, mesh));
      NT diff_x = pi.x() - pj.x();
      NT diff_y = pi.y() - pj.y();

      // local positions (in the 2D param)
      const Local_indices& li = get(lpmap, hd);
      const Point_2& ppi = lp[ li.first ];
      const Point_2& ppj = lp[ li.second ];
      NT p_diff_x = ppi.x() - ppj.x();
      NT p_diff_y = ppi.y() - ppj.y();

      nabla_x = diff_x - (  a * p_diff_x + b * p_diff_y );
      nabla_y = diff_y - ( -b * p_diff_x + a * p_diff_y );

      NT sq_nabla_norm = nabla_x * nabla_x + nabla_y * nabla_y;
      Ef += cot * sq_nabla_norm;
    }

    NT s = a * a + b * b - 1;
    Ef += m_lambda * s * s;
    return Ef;
  }

  // Compute the current energy of a face.
  template <typename VertexUVMap>
  NT compute_current_face_energy(const Triangle_mesh& mesh,
                                 face_descriptor fd,
                                 const Cot_map ctmap,
                                 const Local_points& lp,
                                 const Lp_map lpmap,
                                 const Lt_map ltmap,
                                 const VertexUVMap uvmap) const
  {
    const Lt_matrix& ltm = get(ltmap, fd); // the (current) optimal linear transformation
    NT a = ltm.first;
    NT b = ltm.second;

    return compute_current_face_energy(mesh, fd, ctmap, lp, lpmap, uvmap, a, b);
  }

  // Compute the current energy of the parameterization.
  template <typename VertexUVMap>
  NT compute_current_energy(const Triangle_mesh& mesh,
                            const Faces_vector& faces,
                            const Cot_map ctmap,
                            const Local_points& lp,
                            const Lp_map lpmap,
                            const Lt_map ltmap,
                            const VertexUVMap uvmap) const
  {
    NT E = 0.;

    for(face_descriptor fd : faces) {
      NT Ef = compute_current_face_energy(mesh, fd, ctmap, lp, lpmap,
                                          ltmap, uvmap);
      E += Ef;
    }

    E *= 0.5;
    return E;
  }

// Post processing functions
  // Use the convex virtual boundary algorithm of Karni et al.[2005] to fix
  // the (hopefully few) flips in the result.
  template <typename VertexUVMap,
            typename VertexIndexMap>
  Error_code post_process(const Triangle_mesh& mesh,
                          const Vertex_set& vertices,
                          const Faces_vector& faces,
                          halfedge_descriptor bhd,
                          VertexUVMap uvmap,
                          const VertexIndexMap vimap) const
  {
    typedef MVC_post_processor_3<Triangle_mesh>     Post_processor;

    Post_processor p;
    Error_code status = p.parameterize(mesh, vertices, faces, bhd, uvmap, vimap);
    if(status != OK)
      return status;

#ifdef CGAL_SMP_ARAP_DEBUG
    output_uvmap("ARAP_final_post_processing.off", mesh, vertices, faces, uvmap, vimap);
#endif

    return OK;
  }

// Public operations
public:
  /// returns whether the 3D -> 2D mapping is one-to-one.
  template <typename FaceRange, typename VertexUVMap>
  bool is_one_to_one_mapping(const Triangle_mesh& mesh,
                             const FaceRange& faces,
                             const VertexUVMap uvmap) const
  {
    return internal::is_one_to_one_mapping(mesh, faces, uvmap);
  }

  /// computes a mapping from a triangular 3D surface mesh to a piece of the 2D space.
  /// The mapping is piecewise linear (linear in each triangle).
  /// The result is the (u,v) pair image of each vertex of the 3D surface.
  ///
  /// \tparam VertexUVmap must be a model of `ReadWritePropertyMap` with
  ///         `boost::graph_traits<Triangle_mesh>::%vertex_descriptor` as key type and
  ///         %Point_2 (type deduced from `Triangle_mesh` using `Kernel_traits`)
  ///         as value type.
  /// \tparam VertexIndexMap must be a model of `ReadablePropertyMap` with
  ///         `boost::graph_traits<Triangle_mesh>::%vertex_descriptor` as key type and
  ///         a unique integer as value type.
  /// \tparam VertexParameterizedMap must be a model of `ReadWritePropertyMap` with
  ///         `boost::graph_traits<Triangle_mesh>::%vertex_descriptor` as key type and
  ///         a Boolean as value type.
  ///
  /// \param mesh a triangulated surface.
  /// \param bhd a halfedge descriptor on the boundary of `mesh`.
  /// \param uvmap an instanciation of the class `VertexUVmap`.
  /// \param vimap an instanciation of the class `VertexIndexMap`.
  /// \param vpmap an instanciation of the class `VertexParameterizedMap`.
  ///
  /// \pre `mesh` must be a triangular mesh.
  /// \pre The vertices must be indexed (vimap must be initialized).
  ///
  template <typename VertexUVMap,
            typename VertexIndexMap,
            typename VertexParameterizedMap>
  Error_code parameterize(Triangle_mesh& mesh,
                          halfedge_descriptor bhd,
                          VertexUVMap uvmap,
                          VertexIndexMap vimap,
                          VertexParameterizedMap vpmap)
  {
    CGAL_precondition(is_valid_polygon_mesh(mesh));
    CGAL_precondition(is_triangle_mesh(mesh));
    CGAL_precondition(bhd != boost::graph_traits<Triangle_mesh>::null_halfedge() && is_border(bhd, mesh));

    Error_code status = OK;

    // vertices and faces containers
    Faces_vector faces;
    Vertex_set vertices;
    status = initialize_containers(mesh, bhd, vertices, faces);
    if(status != OK)
      return status;

    // linear transformation matrices L_t
    // Only need to store 2 indices since the complete matrix is {{a,b},{-b,a}}
    Lt_hash_map lt_hm;
    Lt_map ltmap(lt_hm); // will be filled in 'compute_optimal_Lt_matrices()'

    // Compute the initial parameterization of the mesh
    status = compute_initial_uv_map(mesh, bhd, uvmap, vimap);

#ifdef CGAL_SMP_ARAP_DEBUG
    output_uvmap("ARAP_initial_param.off", mesh, vertices, faces, uvmap, vimap);
#endif

    if(status != OK)
      return status;

    // Handle the boundary condition depending on lambda
    status = parameterize_border<VertexUVMap>(mesh, vertices, bhd, vimap, vpmap);
    if(status != OK)
      return status;

    // Compute all cotangent angles
    Cot_hm cthm;
    Cot_map ctmap(cthm);
    status = compute_cotangent_angles(mesh, faces, ctmap);
    if(status != OK)
      return status;

    // Compute all local 2D parameterizations
    Lp_hm lphm;
    Lp_map lpmap(lphm);
    Local_points lp;
    status = compute_local_parameterization(mesh, faces, lp, lpmap);
    if(status != OK)
      return status;

    // The matrix A is constant and can be initialized outside of the loop
    int nbVertices = static_cast<int>(vertices.size());
    Matrix A(nbVertices, nbVertices); // the constant matrix using in the linear system A*X = B
    status = initialize_matrix_A(mesh, vertices, ctmap, vimap, vpmap, A);
    if(status != OK)
      return status;

    NT energy_this = compute_current_energy(mesh, faces, ctmap, lp, lpmap,
                                            ltmap, uvmap);
    NT energy_last;

#ifdef CGAL_PARAMETERIZATION_ARAP_VERBOSE
    std::cout << "Initial energy: " << energy_this << std::endl;
    std::cout << m_iterations << " max iterations" << std::endl;
#endif

    // main loop
    unsigned int ite = 1;
    for(;;)
    {
      compute_optimal_Lt_matrices(mesh, faces, ctmap, lp, lpmap, uvmap, ltmap);
      status = update_solution(mesh, vertices, ctmap, lp, lpmap, ltmap,
                                               uvmap, vimap, vpmap, A);

      // Output the current parameterization
#ifdef CGAL_SMP_ARAP_DEBUG
      output_uvmap("ARAP_iteration_", ite, mesh, vertices, faces, uvmap, vimap);
#endif

      if(status != OK)
        return status;

      // energy based termination
      if(m_tolerance > 0. && ite <= m_iterations) { // if tolerance <= 0, don't compute energy
        energy_last = energy_this;
        energy_this = compute_current_energy(mesh, faces, ctmap, lp, lpmap, ltmap, uvmap);

#ifdef CGAL_PARAMETERIZATION_ARAP_VERBOSE
        std::cout << "Energy at iteration " << ite << " : " << energy_this << std::endl;
#endif

        if(energy_this < 0) {
          // numerical issues can make it so you may get an energy of -1e-17,
          // but it shouldn't be too wrong
          CGAL_assertion(energy_this >= - std::numeric_limits<NT>::epsilon());
          break;
        }

        double energy_diff = std::abs((energy_last - energy_this) / energy_this);
        if(energy_diff < m_tolerance) {
          break;
        }
      }

      if(ite >= m_iterations)
        break;
      else
        ++ite;
    }

#ifdef CGAL_PARAMETERIZATION_ARAP_VERBOSE
    std::cout << "Minimization process ended after: " << ite << " iterations. " << std::endl;
#endif

#ifdef CGAL_SMP_ARAP_DEBUG
    output_uvmap("ARAP_final_pre_processing.off", mesh, vertices, faces, uvmap, vimap);
#endif

    if(!is_one_to_one_mapping(mesh, faces, uvmap)) {
     // Use post processing to handle flipped elements
#ifdef CGAL_PARAMETERIZATION_ARAP_VERBOSE
      std::cout << "Parameterization is not valid; calling post processor" << std::endl;
#endif
      status = post_process(mesh, vertices, faces, bhd, uvmap, vimap);
    }

    return status;
  }

public:
  /// Constructor taking only the parameter \f$\lambda\f$.
  ARAP_parameterizer_3(NT lambda)
    :
      m_borderParameterizer(Border_parameterizer()),
      m_linearAlgebra(Solver_traits()),
      m_lambda(lambda),
      m_lambda_tolerance(1e-10),
      m_iterations(50),
      m_tolerance(1e-6)
  { }

  /// %Default Constructor.
  ///
  /// \param border_param %Object that maps the surface's border to the 2D space.
  /// \param sparse_la %Traits object to access a sparse linear system.
  /// \param lambda Parameter to give importance to shape or angle preservation.
  /// \param iterations Maximal number of iterations in the energy minimization process.
  /// \param tolerance Minimal energy difference between two iterations for the minimization
  ///        process to continue.
  ///
  ARAP_parameterizer_3(Border_parameterizer border_param = Border_parameterizer(),
                       Solver_traits sparse_la = Solver_traits(),
                       NT lambda = 1000.,
                       unsigned int iterations = 50,
                       NT tolerance = 1e-6)
    :
      m_borderParameterizer(border_param),
      m_linearAlgebra(sparse_la),
      m_lambda(lambda),
      m_lambda_tolerance(1e-10),
      m_iterations(iterations),
      m_tolerance(tolerance)
  { }

  // Default copy constructor and operator=() are fine
};

} // namespace Surface_mesh_parameterization

} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_SURFACE_MESH_PARAMETERIZATION_ARAP_PARAMETERIZER_3_H
