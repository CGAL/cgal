// Copyright (c) 2005  INRIA (France).
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
//
// Author(s)     : Laurent Saboret, Pierre Alliez, Bruno Levy

#ifndef CGAL_SURFACE_MESH_PARAMETERIZATION_FIXED_BORDER_PARAMETERIZER_3_H
#define CGAL_SURFACE_MESH_PARAMETERIZATION_FIXED_BORDER_PARAMETERIZER_3_H

#include <CGAL/license/Surface_mesh_parameterization.h>

#include <CGAL/Surface_mesh_parameterization/internal/angles.h>
#include <CGAL/Surface_mesh_parameterization/internal/Containers_filler.h>
#include <CGAL/Surface_mesh_parameterization/internal/kernel_traits.h>

#include <CGAL/Surface_mesh_parameterization/Error_code.h>
#include <CGAL/Surface_mesh_parameterization/Circular_border_parameterizer_3.h>

#include <CGAL/Polygon_mesh_processing/connected_components.h>

#include <CGAL/boost/graph/properties.h>
#include <CGAL/circulator.h>
#include <CGAL/Default.h>

#if defined(CGAL_EIGEN3_ENABLED)
#include <CGAL/Eigen_solver_traits.h>
#endif

#include <boost/foreach.hpp>
#include <boost/function_output_iterator.hpp>
#include <boost/unordered_set.hpp>

/// \file Fixed_border_parameterizer_3.h

namespace CGAL {

namespace Surface_mesh_parameterization {

// ------------------------------------------------------------------------------------
// Declaration
// ------------------------------------------------------------------------------------

/// \ingroup PkgSurfaceParameterizationMethods
///
/// The class `Fixed_border_parameterizer_3`
/// is the base class of fixed border parameterization methods (Tutte, Floater, ...).
///
/// A one-to-one mapping is guaranteed if the border of the surface is mapped onto a convex polygon.
///
/// This class is a pure virtual class and thus cannot be instantiated.
/// Nevertheless, it implements most of the parameterization algorithm `parameterize()`.
/// Subclasses are *Strategies* \cgalCite{cgal:ghjv-dpero-95} that modify the behavior of this algorithm:
/// - They provide the template parameters `BorderParameterizer_` and `SolverTraits_`.
/// - They implement `compute_w_ij()` to compute w_ij = (i, j), coefficient of matrix A
///   for j neighbor vertex of i.
///
// @todo `Fixed_border_parameterizer_3` should remove border vertices
// from the linear systems in order to have a symmetric positive definite
// matrix for Tutte Barycentric Mapping and Discrete Conformal Map algorithms.
///
/// \cgalModels `Parameterizer_3`
///
/// \tparam TriangleMesh_ must be a model of `FaceGraph`.
///
/// \tparam BorderParameterizer_ is a Strategy to parameterize the surface border
///         and must be a model of `Parameterizer_3`.<br>
///         <b>%Default:</b>
/// \code
///   Circular_border_arc_length_parameterizer_3<TriangleMesh_>
/// \endcode
///
/// \tparam SolverTraits_ must be a model of `SparseLinearAlgebraTraits_d`.<br>
///         Note that the system is *not* symmetric because `Fixed_border_parameterizer_3`
///         does not remove border vertices from the system.<br>
///         <b>%Default:</b> If \ref thirdpartyEigen "Eigen" 3.1 (or greater) is available
///         and `CGAL_EIGEN3_ENABLED` is defined, then an overload of `Eigen_solver_traits`
///         is provided as default parameter:
/// \code
///   CGAL::Eigen_solver_traits<
///           Eigen::BiCGSTAB<Eigen_sparse_matrix<double>::EigenType,
///                           Eigen::IncompleteLUT< double > > >
/// \endcode
///
/// \sa `CGAL::Surface_mesh_parameterization::ARAP_parameterizer_3<TriangleMesh, BorderParameterizer, SolverTraits>`
/// \sa `CGAL::Surface_mesh_parameterization::Barycentric_mapping_parameterizer_3<TriangleMesh, BorderParameterizer, SolverTraits>`
/// \sa `CGAL::Surface_mesh_parameterization::Discrete_authalic_parameterizer_3<TriangleMesh, BorderParameterizer, SolverTraits>`
/// \sa `CGAL::Surface_mesh_parameterization::Discrete_conformal_map_parameterizer_3<TriangleMesh, BorderParameterizer, SolverTraits>`
/// \sa `CGAL::Surface_mesh_parameterization::Mean_value_coordinates_parameterizer_3<TriangleMesh, BorderParameterizer, SolverTraits>`
///
template < typename TriangleMesh_,
           class BorderParameterizer_ = Default,
           class SolverTraits_ = Default>
class Fixed_border_parameterizer_3
{
public:
#ifndef DOXYGEN_RUNNING
  typedef typename Default::Get<
    BorderParameterizer_,
    Circular_border_arc_length_parameterizer_3<TriangleMesh_> >::type  Border_parameterizer;

  typedef typename Default::Get<
    SolverTraits_,
  #if defined(CGAL_EIGEN3_ENABLED)
    CGAL::Eigen_solver_traits<
      Eigen::BiCGSTAB<Eigen_sparse_matrix<double>::EigenType,
                      Eigen::IncompleteLUT<double> > >
  #else
    #pragma message("Error: You must either provide 'SolverTraits_' or link CGAL with the Eigen library");
    SolverTraits_ // no parameter provided, and Eigen is not enabled: so don't compile!
  #endif
  >::type                                                     Solver_traits;
#else
  typedef Border_parameterizer_                               Border_parameterizer;
  typedef SolverTraits_                                       Solver_traits;
#endif

  typedef TriangleMesh_                                       TriangleMesh;

// Private types
private:
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;

  typedef CGAL::Vertex_around_target_circulator<TriangleMesh> vertex_around_target_circulator;
  typedef CGAL::Vertex_around_face_circulator<TriangleMesh> vertex_around_face_circulator;

// Protected types
protected:
  // Traits subtypes:
  typedef typename internal::Kernel_traits<TriangleMesh>::Kernel    Kernel;
  typedef typename internal::Kernel_traits<TriangleMesh>::PPM       PPM;
  typedef typename Kernel::FT                                       NT;
  typedef typename Kernel::Point_2                                  Point_2;
  typedef typename Kernel::Point_3                                  Point_3;
  typedef typename Kernel::Vector_3                                 Vector_3;

  // Solver traits subtypes:
  typedef typename Solver_traits::Vector                            Vector;
  typedef typename Solver_traits::Matrix                            Matrix;

// Public operations
public:
  /// Constructor
  Fixed_border_parameterizer_3(Border_parameterizer border_param = Border_parameterizer(),
                               ///< %Object that maps the surface's border to 2D space
                               Solver_traits sparse_la = Solver_traits())
                               ///< Traits object to access a sparse linear system
    : m_borderParameterizer(border_param), m_linearAlgebra(sparse_la)
  { }

  /// Destructor of base class should be virtual.
  virtual ~Fixed_border_parameterizer_3() { }

  // Default copy constructor and operator =() are fine

  /// Compute a one-to-one mapping from a triangular 3D surface mesh
  /// to a piece of the 2D space.
  /// The mapping is piecewise linear (linear in each triangle).
  /// The result is the (u,v) pair image of each vertex of the 3D surface.
  ///
  /// \tparam VertexUVmap must be a model of `ReadWritePropertyMap` with
  ///         `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key type and
  ///         %Point_2 (type deduced from `TriangleMesh` using `Kernel_traits`)
  ///         as value type.
  /// \tparam VertexIndexMap must be a model of `ReadablePropertyMap` with
  ///         `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key type and
  ///         a unique integer as value type.
  /// \tparam VertexParameterizedMap must be a model of `ReadWritePropertyMap` with
  ///         `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key type and
  ///         a Boolean as value type.
  ///
  /// \param mesh a triangulated surface.
  /// \param bhd a halfedge descriptor on the boundary of `mesh`.
  /// \param uvmap an instanciation of the class `VertexUVmap`.
  /// \param vimap an instanciation of the class `VertexIndexMap`.
  /// \param vpmap an instanciation of the class `VertexParameterizedMap`.
  ///
  /// \pre `mesh` must be a triangular mesh.
  /// \pre The mesh border must be mapped onto a convex polygon.
  /// \pre The vertices must be indexed (`vimap` must be initialized)
  ///
  template <typename VertexUVmap,
            typename VertexIndexMap,
            typename VertexParameterizedMap>
  Error_code parameterize(TriangleMesh& mesh,
                          halfedge_descriptor bhd,
                          VertexUVmap uvmap,
                          VertexIndexMap vimap,
                          VertexParameterizedMap vpmap)
  {
    Error_code status = OK;

    typedef boost::unordered_set<vertex_descriptor> Vertex_set;
    Vertex_set vertices;

    internal::Containers_filler<TriangleMesh> fc(mesh, vertices);
    Polygon_mesh_processing::connected_component(
                                        face(opposite(bhd, mesh), mesh),
                                        mesh,
                                        boost::make_function_output_iterator(fc));

    // Count vertices
    int nbVertices = static_cast<int>(vertices.size());

    if (nbVertices == 0)
      return ERROR_EMPTY_MESH;

    // Compute (u,v) for border vertices and mark them as "parameterized"
    status = get_border_parameterizer().parameterize(mesh, bhd, uvmap, vimap, vpmap);

    if (status != OK)
      return status;

    // Create two sparse linear systems "A*Xu = Bu" and "A*Xv = Bv" (one line/column per vertex)
    Matrix A(nbVertices, nbVertices);
    Vector Xu(nbVertices), Xv(nbVertices), Bu(nbVertices), Bv(nbVertices);

    // Initialize A, Xu, Xv, Bu and Bv after border parameterization
    // Fill the border vertices' lines in both linear systems:
    // "u = constant" and "v = constant"
    //
    // @todo Fixed_border_parameterizer_3 should remove border vertices
    // from the linear systems in order to have a symmetric positive definite
    // matrix for Tutte Barycentric Mapping and Discrete Conformal Map algorithms.
    initialize_system_from_mesh_border(A, Bu, Bv, mesh, bhd, uvmap, vimap);

    // Fill the matrix for the inner vertices v_i: compute A's coefficient
    // w_ij for each neighbor j; then w_ii = - sum of w_ijs
    boost::unordered_set<vertex_descriptor> main_border;

    BOOST_FOREACH(vertex_descriptor v, vertices_around_face(bhd,mesh)){
      main_border.insert(v);
    }

    int count = 0;
    BOOST_FOREACH(vertex_descriptor v, vertices){
      // inner vertices only
      if(main_border.find(v) == main_border.end()){
        // Compute the line i of matrix A for i inner vertex
        status = setup_inner_vertex_relations(A, Bu, Bv, mesh, v, vimap);
        if(status != OK)
          return status;
      } else {
        count++;
      }
    }

    // Solve "A*Xu = Bu". On success, solution is (1/Du) * Xu.
    // Solve "A*Xv = Bv". On success, solution is (1/Dv) * Xv.
    NT Du, Dv;
    if(!get_linear_algebra_traits().linear_solver(A, Bu, Xu, Du) ||
       !get_linear_algebra_traits().linear_solver(A, Bv, Xv, Dv))
    {
      status = ERROR_CANNOT_SOLVE_LINEAR_SYSTEM;
    }

    if(status != OK)
        return status;

    // WARNING: this package does not support homogeneous coordinates!
    CGAL_assertion(Du == 1.0);
    CGAL_assertion(Dv == 1.0);

    // Copy Xu and Xv coordinates into the (u,v) pair of each vertex
    BOOST_FOREACH(vertex_descriptor v, vertices)
    {
      // inner vertices only
      if(main_border.find(v) == main_border.end()){
        int index = get(vimap,v);
        put(uvmap,v,Point_2(Xu[index],Xv[index]));
        put(vpmap,v,true);
      }
    }

    // Check postconditions
    // AF status = check_parameterize_postconditions(amesh, A, Bu, Bv);

    if(status != OK)
      return status;

    return status;
  }

// Protected operations
protected:
  /// Initialize A, Bu and Bv after border parameterization.
  /// Fill the border vertices' lines in both linear systems:
  /// "u = constant" and "v = constant".
  ///
  /// \tparam VertexUVmap must be a model of `ReadWritePropertyMap` with
  ///         `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key type and
  ///         %Point_2 (type deduced from `TriangleMesh` using `Kernel_traits`)
  ///         as value type.
  /// \tparam VertexIndexMap must be a model of `ReadablePropertyMap` with
  ///         `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key type and
  ///         a unique integer as value type.
  ///
  /// \param A the matrix in both linear system
  /// \param Bu the right hand side vector in the linear system of x coordinates
  /// \param Bv the right hand side vector in the linear system of y coordinates
  /// \param mesh a triangulated surface.
  /// \param bhd a halfedge descriptor on the boundary of `mesh`.
  /// \param uvmap an instanciation of the class `VertexUVmap`.
  /// \param vimap an instanciation of the class `VertexIndexMap`.
  ///
  /// \pre Vertices must be indexed (`vimap` must be initialized).
  /// \pre A, Bu and Bv must be allocated.
  /// \pre Border vertices must be parameterized.
  template <typename VertexUVmap, typename VertexIndexMap>
  void initialize_system_from_mesh_border(Matrix& A, Vector& Bu, Vector& Bv,
                                          const TriangleMesh& mesh,
                                          halfedge_descriptor bhd,
                                          VertexUVmap uvmap,
                                          VertexIndexMap vimap) const
  {
    BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_face(bhd, mesh)){
      // Get vertex index in sparse linear system
      int index = get(vimap, target(hd, mesh));
      // Write a diagonal coefficient of A
      A.set_coef(index, index, 1, true /*new*/);
      // get the halfedge uv
      // Write constant in Bu and Bv
      const Point_2& uv = get(uvmap, target(hd, mesh));
      Bu[index] = uv.x();
      Bv[index] = uv.y();
    }
  }

  /// Compute w_ij, coefficient of matrix A for j neighbor vertex of i.
  /// Implementation note: Subclasses must at least implement compute_w_ij().
  ///
  /// \param mesh a triangulated surface.
  /// \param main_vertex_v_i the vertex of `mesh` with index `i`
  /// \param neighbor_vertex_v_j the vertex of `mesh` with index `j`
  virtual NT compute_w_ij(const TriangleMesh& mesh,
                          vertex_descriptor main_vertex_v_i,
                          vertex_around_target_circulator neighbor_vertex_v_j) const
  = 0;

  /// Compute the line i of matrix A for i inner vertex:
  /// - call compute_w_ij() to compute the A coefficient w_ij for each neighbor v_j.
  /// - compute w_ii = - sum of w_ijs.
  ///
  /// \pre Vertices must be indexed.
  /// \pre Vertex i musn't be already parameterized.
  /// \pre Line i of A must contain only zeros.
  // TODO: check if this must be virtual
  // virtual
  template <typename VertexIndexMap>
  Error_code setup_inner_vertex_relations(Matrix& A,
                                          Vector&,
                                          Vector&,
                                          const TriangleMesh& mesh,
                                          vertex_descriptor vertex,
                                          VertexIndexMap vimap) const
  {
    int i = get(vimap,vertex);

    // circulate over vertices around 'vertex' to compute w_ii and w_ijs
    NT w_ii = 0;
    int vertexIndex = 0;

    vertex_around_target_circulator v_j(halfedge(vertex, mesh), mesh), end = v_j;
    CGAL_For_all(v_j, end){
      // Call to virtual method to do the actual coefficient computation
      NT w_ij = -1.0 * compute_w_ij(mesh, vertex, v_j);
      // w_ii = - sum of w_ijs
      w_ii -= w_ij;

      // Get j index
      int j = get(vimap, *v_j);

      // Set w_ij in matrix
      A.set_coef(i,j, w_ij, true /*new*/);
      vertexIndex++;
    }

    if (vertexIndex < 2)
      return ERROR_NON_TRIANGULAR_MESH;

    // Set w_ii in matrix
    A.set_coef(i,i, w_ii, true /*new*/);
    return OK;
  }

// Protected accessors
protected:
  /// Get the object that maps the surface's border onto a 2D space.
  Border_parameterizer& get_border_parameterizer() { return m_borderParameterizer; }

  /// Get the sparse linear algebra (traits object to access the linear system).
  Solver_traits& get_linear_algebra_traits() { return m_linearAlgebra; }

// Fields
private:
  // Object that maps the surface's border onto a 2D space.
  Border_parameterizer m_borderParameterizer;

  // Traits object to solve a sparse linear system
  Solver_traits m_linearAlgebra;
};

} // namespace Surface_mesh_parameterization

} // namespace CGAL

#endif // CGAL_SURFACE_MESH_PARAMETERIZATION_FIXED_BORDER_PARAMETERIZER_3_H
