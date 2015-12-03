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
// $URL$
// $Id$
//
//
// Author(s)     : Laurent Saboret, Pierre Alliez, Bruno Levy


#ifndef CGAL_FIXED_BORDER_PARAMETERIZER_3_H
#define CGAL_FIXED_BORDER_PARAMETERIZER_3_H

#include <CGAL/license/Surface_mesh_parameterization.h>


#include <CGAL/circulator.h>
#include <CGAL/Timer.h>
#include <CGAL/Eigen_solver_traits.h>


#include <CGAL/Parameterizer_traits_3.h>
#include <CGAL/Circular_border_parameterizer_3.h>
#include <CGAL/surface_mesh_parameterization_assertions.h>
#include <boost/foreach.hpp>
#include <iostream>
#include <set>

/// \file Fixed_border_parameterizer_3.h

namespace CGAL {


// ------------------------------------------------------------------------------------
// Declaration
// ------------------------------------------------------------------------------------


/// \ingroup  PkgSurfaceParameterizationMethods
///
/// The class `Fixed_border_parameterizer_3`
/// is the base class of fixed border parameterization methods (Tutte, Floater, ...).
///
/// One-to-one mapping is guaranteed if surface's border is mapped onto a convex polygon.
///
/// This class is a pure virtual class, thus cannot be instantiated.
/// Anyway, it implements most of the parameterization algorithm `parameterize()`.
/// Subclasses are Strategies that modify the behavior of this algorithm:
/// - They provide `BorderParameterizer_3` and `SparseLinearAlgebraTraits_d` template
///   parameters.
/// - They implement `compute_w_ij()` to compute w_ij = (i, j) coefficient of matrix A
///   for j neighbor vertex of i.
/// - They may implement an optimized version of `is_one_to_one_mapping()`.
///
// @todo `Fixed_border_parameterizer_3` should remove border vertices
/// from the linear systems in order to have a symmetric positive definite
/// matrix for Tutte Barycentric Mapping and Discrete Conformal Map algorithms.
///
/// \cgalModels `ParameterizerTraits_3`
///
///
///
/// \sa `CGAL::Parameterizer_traits_3<ParameterizationMesh_3>`
/// \sa `CGAL::Barycentric_mapping_parameterizer_3<ParameterizationMesh_3, BorderParameterizer_3, SparseLinearAlgebraTraits_d>`
/// \sa `CGAL::Discrete_authalic_parameterizer_3<ParameterizationMesh_3, BorderParameterizer_3, SparseLinearAlgebraTraits_d>`
/// \sa `CGAL::Discrete_conformal_map_parameterizer_3<ParameterizationMesh_3, BorderParameterizer_3, SparseLinearAlgebraTraits_d>`
/// \sa `CGAL::LSCM_parameterizer_3<ParameterizationMesh_3, BorderParameterizer_3, SparseLinearAlgebraTraits_d>`
/// \sa `CGAL::Mean_value_coordinates_parameterizer_3<ParameterizationMesh_3, BorderParameterizer_3, SparseLinearAlgebraTraits_d>`

template
<
    class ParameterizationMesh_3,       ///< 3D surface mesh
    class BorderParameterizer_3         ///< Strategy to parameterize the surface border
                = Circular_border_arc_length_parameterizer_3<ParameterizationMesh_3>,
    class SparseLinearAlgebraTraits_d   ///< Traits class to solve a sparse linear system
                = Eigen_solver_traits<Eigen::BiCGSTAB<Eigen_sparse_matrix<double>::EigenType, Eigen::IncompleteLUT< double > > > >

class Fixed_border_parameterizer_3
    : public Parameterizer_traits_3<ParameterizationMesh_3>
{
  // Private types
private:
  typedef Parameterizer_traits_3<ParameterizationMesh_3> Base;

  // Public types
public:
  // We have to repeat the types exported by superclass
  /// @cond SKIP_IN_MANUAL
  typedef typename Base::Error_code       Error_code;
  typedef ParameterizationMesh_3          TriangleMesh;
  /// @endcond
  
  /// Export BorderParameterizer_3 template parameter.
  typedef BorderParameterizer_3           Border_param;
  /// Export SparseLinearAlgebraTraits_d template parameter.
  typedef SparseLinearAlgebraTraits_d     Sparse_LA;
  
  // Private types
private:
  
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;
  
  typedef CGAL::Vertex_around_target_circulator<TriangleMesh> vertex_around_target_circulator;
  typedef CGAL::Vertex_around_face_circulator<TriangleMesh> vertex_around_face_circulator;
  typedef CGAL::Halfedge_around_target_circulator<TriangleMesh> halfedge_around_target_circulator;
  
  // Mesh_Adaptor_3 subtypes:
  typedef typename Base::NT            NT;
  typedef typename Base::Point_2       Point_2;
  typedef typename Base::Point_3       Point_3;
  typedef typename Base::Vector_3      Vector_3;
  
  
  // SparseLinearAlgebraTraits_d subtypes:
  typedef typename Sparse_LA::Vector      Vector;
  typedef typename Sparse_LA::Matrix      Matrix;
  
protected:
  // Using statements needed for derived class
  using Base::compute_angle_rad;
  using Base::cotangent;
  
  // Public operations
public:
  /// Constructor
  Fixed_border_parameterizer_3(Border_param border_param = Border_param(),
                               ///< Object that maps the surface's border to 2D space
                               Sparse_LA sparse_la = Sparse_LA())
    ///< Traits object to access a sparse linear system
    : m_borderParameterizer(border_param), m_linearAlgebra(sparse_la)
  {}
  
  // Default copy constructor and operator =() are fine
  
  /// Compute a one-to-one mapping from a triangular 3D surface mesh
  /// to a piece of the 2D space.
  /// The mapping is linear by pieces (linear in each triangle).
  /// The result is the (u,v) pair image of each vertex of the 3D surface.
  ///
  /// \pre `mesh` must be a surface with one connected component.
  /// \pre `mesh` must be a triangular mesh.
  /// \pre The mesh border must be mapped onto a convex polygon.
  
  template <typename VertexUVmap, typename VertexIndexMap, typename VertexParameterizedMap>
  Error_code  parameterize(TriangleMesh& mesh,
                           halfedge_descriptor,
                           VertexUVmap,
                           VertexIndexMap vimap,
                           VertexParameterizedMap vpm);
  
  // Protected operations
protected:
  /// Check parameterize() preconditions:
  /// - `mesh` must be a surface with one connected component.
  /// - `mesh` must be a triangular mesh.
  /// - The mesh border must be mapped onto a convex polygon.
  /// virtual Error_code  check_parameterize_preconditions(TriangleMesh& mesh);
  
  /// Initialize A, Bu and Bv after border parameterization.
  /// Fill the border vertices' lines in both linear systems:
  /// "u = constant" and "v = constant".
  ///
  /// \pre Vertices must be indexed.
  /// \pre A, Bu and Bv must be allocated.
  /// \pre Border vertices must be parameterized.
  template <typename VertexUVmap, typename VertexIndexMap >
  void  initialize_system_from_mesh_border (Matrix& A, Vector& Bu, Vector& Bv,
                                            const TriangleMesh& tmesh,
                                            halfedge_descriptor bhd,
                                            VertexUVmap uvmap,
                                            VertexIndexMap vimap)
  {
    // AF: loop over border halfedges
    BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_face(bhd,tmesh)){
      // AF: get the halfedge-as-vertex index
      // Get vertex index in sparse linear system
      int index = get(vimap, target(opposite(next(hd,tmesh),tmesh),tmesh));
      // Write a diagonal coefficient of A
      A.set_coef(index, index, 1, true /*new*/);
      //std::cerr << index << "  " << index << " 1" << std::endl;
      // get the halfedge uv
      // Write constant in Bu and Bv
      Point_2 uv = get(uvmap, target(opposite(next(hd,tmesh),tmesh),tmesh));
      Bu[index] = uv.x();
      Bv[index] = uv.y();
    }
  }

  /// Compute w_ij = (i, j) coefficient of matrix A for j neighbor vertex of i.
  /// Implementation note: Subclasses must at least implement compute_w_ij().
  virtual NT compute_w_ij(const TriangleMesh& mesh,
                          vertex_descriptor main_vertex_v_i,
                          halfedge_around_target_circulator neighbor_vertex_v_j)
  = 0;

  /// Compute the line i of matrix A for i inner vertex:/// - call compute_w_ij() to compute the A coefficient w_ij for each neighbor v_j.
  /// - compute w_ii = - sum of w_ijs.
  ///
  /// \pre Vertices must be indexed.
  /// \pre Vertex i musn't be already parameterized.
  /// \pre Line i of A must contain only zeros.
  // TODO: check if this must be virtual 
  // virtual 
  template <typename VertexIndexMap>
  Error_code setup_inner_vertex_relations(Matrix& A,
                                          Vector& Bu,
                                          Vector& Bv,
                                          const TriangleMesh& mesh,
                                          vertex_descriptor vertex,
                                          VertexIndexMap vimap)
  {
    int i = get(vimap,vertex);
  
    // circulate over vertices around 'vertex' to compute w_ii and w_ijs
    // use halfedge_around_target to get the right "vertex" if it is on a seam
    NT w_ii = 0;
    int vertexIndex = 0;
  
    halfedge_around_target_circulator v_j(halfedge(vertex,mesh), mesh), end = v_j;
    CGAL_For_all(v_j, end){
      // Call to virtual method to do the actual coefficient computation
      NT w_ij = -1.0 * compute_w_ij(mesh, vertex, v_j);
    
      // w_ii = - sum of w_ijs
      w_ii -= w_ij;
    
      // Get j index
      int j = get(vimap, target(opposite(*v_j,mesh),mesh));
    
      // Set w_ij in matrix
      A.set_coef(i,j, w_ij, true /*new*/);
      //std::cout << i << " " << j << " "  << w_ij << std::endl;
      vertexIndex++;
    }
    if (vertexIndex < 2)
      return Base::ERROR_NON_TRIANGULAR_MESH;
  
    // Set w_ii in matrix
    A.set_coef(i,i, w_ii, true /*new*/);
    // std::cout << i << " " << i << " "  << w_ii << std::endl;
    return Base::OK;
  }
  


#if 0
  /// Check parameterize() postconditions:
  /// - 3D -> 2D mapping is one-to-one.
  virtual Error_code check_parameterize_postconditions(const TriangleMesh& mesh,
                                                       const Matrix& A,
                                                       const Vector& Bu,
                                                       const Vector& Bv);

  /// Check if 3D -> 2D mapping is one-to-one.
  /// The default implementation checks each normal.
  virtual bool  is_one_to_one_mapping(const TriangleMesh& mesh,
                                      const Matrix& A,
                                      const Vector& Bu,
                                      const Vector& Bv);
#endif 

  // Protected accessors
protected:
  /// Get the object that maps the surface's border onto a 2D space.
  Border_param&   get_border_parameterizer()    { return m_borderParameterizer; }

  /// Get the sparse linear algebra (traits object to access the linear system).
  Sparse_LA&      get_linear_algebra_traits() { return m_linearAlgebra; }

  // Fields
private:
  /// Object that maps the surface's border onto a 2D space.
  Border_param    m_borderParameterizer;

  /// Traits object to solve a sparse linear system
  Sparse_LA       m_linearAlgebra;
};


// ------------------------------------------------------------------------------------
// Implementation
// ------------------------------------------------------------------------------------

// Compute a one-to-one mapping from a triangular 3D surface mesh
// to a piece of the 2D space.
// The mapping is linear by pieces (linear in each triangle).
// The result is the (u,v) pair image of each vertex of the 3D surface.
//
// Preconditions:
// - `mesh` must be a surface with one connected component.
// - `mesh` must be a triangular mesh.
// - The mesh border must be mapped onto a convex polygon.
template <class TriangleMesh, class Border_param, class Sparse_LA> 
template <typename VertexUVmap, typename VertexIndexMap, typename VertexParameterizedMap>
typename Fixed_border_parameterizer_3<TriangleMesh, Border_param, Sparse_LA>::Error_code
Fixed_border_parameterizer_3<TriangleMesh, Border_param, Sparse_LA>::
parameterize(TriangleMesh& mesh, halfedge_descriptor bhd, VertexUVmap uvmap, VertexIndexMap vimap,
                           VertexParameterizedMap vpm)
{

#ifdef DEBUG_TRACE
    // Create timer for traces
    CGAL::Timer timer;
    timer.start();
#endif

    // Check preconditions
    Error_code status = Base::OK; // AF check_parameterize_preconditions(amesh);
#ifdef DEBUG_TRACE
    std::cerr << "  parameterization preconditions: " << timer.time() << " seconds." << std::endl;
    timer.reset();
#endif
    if (status != Base::OK)
        return status;

    // Count vertices
    int nbVertices= num_vertices(mesh);

    // Compute (u,v) for border vertices
    // and mark them as "parameterized"
    status = get_border_parameterizer().parameterize_border(mesh,bhd,uvmap,vpm);
#ifdef DEBUG_TRACE
    std::cerr << "  border vertices parameterization: " << timer.time() << " seconds." << std::endl;
    timer.reset();
#endif
    if (status != Base::OK)
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
    initialize_system_from_mesh_border (A, Bu, Bv, mesh, bhd, uvmap, vimap);

    // AF: no change, as this are only concerns inner vertices
    // Fill the matrix for the inner vertices v_i: compute A's coefficient
    // w_ij for each neighbor j; then w_ii = - sum of w_ijs
    std::set<vertex_descriptor> main_border;

    BOOST_FOREACH(vertex_descriptor v, vertices_around_face(bhd,mesh)){
      main_border.insert(v);
    }

    BOOST_FOREACH(vertex_descriptor v, vertices(mesh)){
      // inner vertices only
      if( main_border.find(v) == main_border.end() )
        {
            // Compute the line i of matrix A for i inner vertex
            status = setup_inner_vertex_relations(A, Bu, Bv,
                                                  mesh,
                                                  v,
                                                  vimap);
            if (status != Base::OK)
                return status;
        }
    }
#ifdef DEBUG_TRACE
    std::cerr << "  matrix filling (" << nbVertices << " x " << nbVertices << "): "
              << timer.time() << " seconds." << std::endl;
    timer.reset();
#endif

    // Solve "A*Xu = Bu". On success, solution is (1/Du) * Xu.
    // Solve "A*Xv = Bv". On success, solution is (1/Dv) * Xv.
    NT Du, Dv;
    if ( !get_linear_algebra_traits().linear_solver(A, Bu, Xu, Du) ||
         !get_linear_algebra_traits().linear_solver(A, Bv, Xv, Dv) )
    {
        status = Base::ERROR_CANNOT_SOLVE_LINEAR_SYSTEM;
    }
#ifdef DEBUG_TRACE
    std::cerr << "  solving two linear systems: "
              << timer.time() << " seconds." << std::endl;
    timer.reset();
#endif
    if (status != Base::OK)
        return status;

    // WARNING: this package does not support homogeneous coordinates!
    CGAL_surface_mesh_parameterization_assertion(Du == 1.0);
    CGAL_surface_mesh_parameterization_assertion(Dv == 1.0);

    // Copy Xu and Xv coordinates into the (u,v) pair of each vertex

    BOOST_FOREACH(vertex_descriptor v, vertices(mesh))
    {
      // inner vertices only
      if( main_border.find(v) == main_border.end() )
        {
          int index = get(vimap,v);
          put(uvmap,v,Point_2(Xu[index],Xv[index]));
          put(vpm,v,true);
          /*
          BOOST_FOREACH(halfedge_descriptor hd, halfedges_around_target(v,mesh)){
            int index = get(vimap,hd);
            put(uvmap,hd,Point_2(Xu[index],Xv[index]));
          }
          */
        }
    }
#ifdef DEBUG_TRACE
    std::cerr << "  copy computed UVs to mesh :"
              << timer.time() << " seconds." << std::endl;
    timer.reset();
#endif

    // Check postconditions
    // AF status = check_parameterize_postconditions(amesh, A, Bu, Bv);
#ifdef DEBUG_TRACE
    std::cerr << "  parameterization postconditions: " << timer.time() << " seconds." << std::endl;
#endif
    if (status != Base::OK)
        return status;

    return status;
}

#if 0
// Check parameterize() preconditions:
// - `mesh` must be a surface with one connected component.
// - `mesh` must be a triangular mesh.
// - The mesh border must be mapped onto a convex polygon.
template<class TriangleMesh, class Border_param, class Sparse_LA>
inline
typename Fixed_border_parameterizer_3<TriangleMesh, Border_param, Sparse_LA>::Error_code
Fixed_border_parameterizer_3<TriangleMesh, Border_param, Sparse_LA>::
check_parameterize_preconditions(TriangleMesh& amesh)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_iterator vertex_iterator;

    const TriangleMesh& mesh = amesh.get_adapted_mesh();

    Error_code status = Base::OK;	    // returned value

    // Helper class to compute genus or count borders, vertices, ...
    //typedef Parameterization_mesh_feature_extractor<TriangleMesh>
    //                                        Mesh_feature_extractor;
    //Mesh_feature_extractor feature_extractor(amesh);

    // Check that mesh is not empty
    vertex_iterator b, e;
    boost::tie(b,e) = vertices(mesh);
    if (b == e)
        status = Base::ERROR_EMPTY_MESH;
    if (status != Base::OK)
        return status;

    // The whole surface parameterization package is restricted to triangular meshes
    status = is_triangle_mesh(mesh) ? Base::OK
                                    : Base::ERROR_NON_TRIANGULAR_MESH;
    if (status != Base::OK)
        return status;

    // The whole package is restricted to surfaces: genus = 0,
    // one connected component and at least one border
    int genus = feature_extractor.get_genus();
    int nb_borders = feature_extractor.get_nb_borders();
    int nb_components = feature_extractor.get_nb_connex_components();
    status = (genus == 0 && nb_borders >= 1 && nb_components == 1)
           ? Base::OK
           : Base::ERROR_NO_TOPOLOGICAL_DISC;
    if (status != Base::OK)
        return status;

    // One-to-one mapping is guaranteed if all w_ij coefficients are > 0 (for j vertex neighbor of i)
    // and if the surface border is mapped onto a 2D convex polygon
    status = get_border_parameterizer().is_border_convex()
           ? Base::OK
           : Base::ERROR_NON_CONVEX_BORDER;
    if (status != Base::OK)
        return status;

    return status;
}
#endif

#if 0
// Initialize A, Bu and Bv after border parameterization.
// Fill the border vertices' lines in both linear systems: "u = constant" and "v = constant".
//
// Preconditions:
// - Vertices must be indexed.
// - A, Bu and Bv must be allocated.
// - Border vertices must be parameterized.
template<class TriangleMesh, class Border_param, class Sparse_LA>
inline
template<class VertexUVmap>
void Fixed_border_parameterizer_3<TriangleMesh, Border_param, Sparse_LA>::
initialize_system_from_mesh_border (Matrix& A, Vector& Bu, Vector& Bv,
                                    const TriangleMesh& tmesh,
                                    halfedge_descriptor bhd,
                                    VertexUVmap uvmap)
#endif

#if 0
// Compute the line i of matrix A for i inner vertex:
// - call compute_w_ij() to compute the A coefficient w_ij for each neighbor v_j.
// - compute w_ii = - sum of w_ijs.
//
// Preconditions:
// - Vertices must be indexed.
// - Vertex i must not be already parameterized.
// - Line i of A must contain only zeros.
template<class TriangleMesh, class Border_param, class Sparse_LA>
inline
template <typename VertexIndexMap>
typename Fixed_border_parameterizer_3<TriangleMesh, Border_param, Sparse_LA>::Error_code
Fixed_border_parameterizer_3<TriangleMesh, Border_param, Sparse_LA>::
setup_inner_vertex_relations(Matrix& A,
                             Vector& ,
                             Vector& ,
                             const TriangleMesh& mesh,
                             vertex_descriptor vertex,
                             VertexIndexMap vimap)

// Copy Xu and Xv coordinates into the (u,v) pair of each surface vertex.
// AF: this only concerns vertices NOT on the border
template<class TriangleMesh, class Border_param, class Sparse_LA>
#endif


#if 0
// Check parameterize() postconditions:
// - 3D -> 2D mapping is one-to-one.
template<class TriangleMesh, class Border_param, class Sparse_LA>
inline
typename Fixed_border_parameterizer_3<TriangleMesh, Border_param, Sparse_LA>::Error_code
Fixed_border_parameterizer_3<TriangleMesh, Border_param, Sparse_LA>::
check_parameterize_postconditions(const TriangleMesh& mesh,
                                  const Matrix& A,
                                  const Vector& Bu,
                                  const Vector& Bv)
{
    Error_code status = Base::OK;

    // Check if 3D -> 2D mapping is one-to-one
    status = is_one_to_one_mapping(mesh, A, Bu, Bv)
           ? Base::OK
           : Base::ERROR_NO_1_TO_1_MAPPING;
    if (status != Base::OK)
        return status;

    return status;
}

// Check if 3D -> 2D mapping is one-to-one.
// The default implementation checks each normal.
template<class TriangleMesh, class Border_param, class Sparse_LA>
inline
bool Fixed_border_parameterizer_3<TriangleMesh, Border_param, Sparse_LA>::
is_one_to_one_mapping(const TriangleMesh& mesh,
                      const Matrix& ,
                      const Vector& ,
                      const Vector& )
{
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;
  
    Vector_3 first_triangle_normal = NULL_VECTOR; // initialize to avoid warning

    BOOST_FOREACH(face_descriptor f, faces(mesh))
    {
        // Get 3 vertices of the facet
        vertex_descriptor v0, v1, v2;
        int vertexIndex = 0;
        vertex_around_face_circulator cir(halfedge(f,mesh),mesh), first(cir), end(cir);
        CGAL_For_all(cir, end)
        {
            if (vertexIndex == 0)
                v0 = *cir;
            else if (vertexIndex == 1)
                v1 = *cir;
            else if (vertexIndex == 2)
                v2 = *cir;

            vertexIndex++;
        }
        CGAL_surface_mesh_parameterization_assertion(vertexIndex >= 3);

        // Get the 3 vertices position IN 2D
        Point_2 p0 = amesh.get_vertex_uv(v0) ;
        Point_2 p1 = amesh.get_vertex_uv(v1) ;
        Point_2 p2 = amesh.get_vertex_uv(v2) ;

        // Compute the facet normal
        Point_3 p0_3D(p0.x(), p0.y(), 0);
        Point_3 p1_3D(p1.x(), p1.y(), 0);
        Point_3 p2_3D(p2.x(), p2.y(), 0);
        Vector_3 v01_3D = p1_3D - p0_3D;
        Vector_3 v02_3D = p2_3D - p0_3D;
        Vector_3 normal = CGAL::cross_product(v01_3D, v02_3D);

        // Check that all normals are oriented the same way
        // => no 2D triangle is flipped
        if (cir == first)
        {
            first_triangle_normal = normal;
        }
        else
        {
            if (first_triangle_normal * normal < 0)
                return false;
        }
    }

    return true;            // OK if we reach this point
}
#endif 

} //namespace CGAL

#endif //CGAL_FIXED_BORDER_PARAMETERIZER_3_H
