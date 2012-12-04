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


#ifndef CGAL_LSCM_PARAMETERIZER_3_H
#define CGAL_LSCM_PARAMETERIZER_3_H

#include <CGAL/circulator.h>
#include <CGAL/Timer.h>
#include <CGAL/OpenNL/linear_solver.h>

#include <CGAL/Parameterizer_traits_3.h>
#include <CGAL/Two_vertices_parameterizer_3.h>
#include <CGAL/surface_mesh_parameterization_assertions.h>

#include <iostream>

/// \file LSCM_parameterizer_3.h

namespace CGAL {


// ------------------------------------------------------------------------------------
// Declaration
// ------------------------------------------------------------------------------------

/// \ingroup  PkgSurfaceParameterizationMethods
///
/// The class LSCM_parameterizer_3 implements the
/// *Least Squares Conformal Maps (LSCM)* parameterization  \cite cgal:lprm-lscm-02.
///
/// This is a conformal parameterization, i.e. it attempts to preserve angles.
///
/// This is a free border parameterization. No need to map the surface's border
/// onto a convex polygon (only two pinned vertices are needed to ensure a
/// unique solution), but one-to-one mapping is *not* guaranteed.
///
/// \cgalModels `ParameterizerTraits_3`
///
///
/// \sa `CGAL::Parameterizer_traits_3<ParameterizationMesh_3>`
/// \sa `CGAL::Fixed_border_parameterizer_3<ParameterizationMesh_3, BorderParameterizer_3, SparseLinearAlgebraTraits_d>`
/// \sa `CGAL::Barycentric_mapping_parameterizer_3<ParameterizationMesh_3, BorderParameterizer_3, SparseLinearAlgebraTraits_d>`
/// \sa `CGAL::Discrete_authalic_parameterizer_3<ParameterizationMesh_3, BorderParameterizer_3, SparseLinearAlgebraTraits_d>`
/// \sa `CGAL::Discrete_conformal_map_parameterizer_3<ParameterizationMesh_3, BorderParameterizer_3, SparseLinearAlgebraTraits_d>`
/// \sa `CGAL::Mean_value_coordinates_parameterizer_3<ParameterizationMesh_3, BorderParameterizer_3, SparseLinearAlgebraTraits_d>`

template
<
    class ParameterizationMesh_3,     ///< 3D surface mesh.
    class BorderParameterizer_3
                = Two_vertices_parameterizer_3<ParameterizationMesh_3>,
                                      ///< Strategy to parameterize the surface border.
                                      ///< The minimum is to parameterize two vertices.
    class SparseLinearAlgebraTraits_d
                = OpenNL::SymmetricLinearSolverTraits<typename ParameterizationMesh_3::NT>
                                      ///< Traits class to solve a sparse linear system.
                                      ///< We may use a symmetric definite positive solver because LSCM
                                      ///< solves the system in the least squares sense.
>
class LSCM_parameterizer_3
    : public Parameterizer_traits_3<ParameterizationMesh_3>
{
// Private types
private:
    // Superclass
    typedef Parameterizer_traits_3<ParameterizationMesh_3>
                                            Base;

// Public types
public:
    // We have to repeat the types exported by superclass
    /// @cond SKIP_IN_MANUAL
    typedef typename Base::Error_code       Error_code;
    typedef ParameterizationMesh_3          Adaptor;
    /// @endcond

    /// Export BorderParameterizer_3 template parameter.
    typedef BorderParameterizer_3           Border_param;
    /// Export SparseLinearAlgebraTraits_d template parameter.
    typedef SparseLinearAlgebraTraits_d     Sparse_LA;

// Private types
private:
    // Mesh_Adaptor_3 subtypes:
    typedef typename Adaptor::NT            NT;
    typedef typename Adaptor::Point_2       Point_2;
    typedef typename Adaptor::Point_3       Point_3;
    typedef typename Adaptor::Vector_2      Vector_2;
    typedef typename Adaptor::Vector_3      Vector_3;
    typedef typename Adaptor::Facet         Facet;
    typedef typename Adaptor::Facet_handle  Facet_handle;
    typedef typename Adaptor::Facet_const_handle
                                            Facet_const_handle;
    typedef typename Adaptor::Facet_iterator Facet_iterator;
    typedef typename Adaptor::Facet_const_iterator
                                            Facet_const_iterator;
    typedef typename Adaptor::Vertex        Vertex;
    typedef typename Adaptor::Vertex_handle Vertex_handle;
    typedef typename Adaptor::Vertex_const_handle
                                            Vertex_const_handle;
    typedef typename Adaptor::Vertex_iterator Vertex_iterator;
    typedef typename Adaptor::Vertex_const_iterator
                                            Vertex_const_iterator;
    typedef typename Adaptor::Border_vertex_iterator
                                            Border_vertex_iterator;
    typedef typename Adaptor::Border_vertex_const_iterator
                                            Border_vertex_const_iterator;
    typedef typename Adaptor::Vertex_around_facet_circulator
                                            Vertex_around_facet_circulator;
    typedef typename Adaptor::Vertex_around_facet_const_circulator
                                            Vertex_around_facet_const_circulator;
    typedef typename Adaptor::Vertex_around_vertex_circulator
                                            Vertex_around_vertex_circulator;
    typedef typename Adaptor::Vertex_around_vertex_const_circulator
                                            Vertex_around_vertex_const_circulator;

    // SparseLinearAlgebraTraits_d subtypes:
    typedef typename Sparse_LA::Vector      Vector;
    typedef typename Sparse_LA::Matrix      Matrix;

    typedef typename OpenNL::LinearSolver<Sparse_LA>
                                            LeastSquaresSolver ;

// Public operations
public:
    /// Constructor
    LSCM_parameterizer_3(Border_param border_param = Border_param(),
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
    virtual Error_code  parameterize(Adaptor& mesh);

// Private operations
private:
    /// Check parameterize() preconditions:
    /// - `mesh` must be a surface with one connected component.
    /// - `mesh` must be a triangular mesh.
    virtual Error_code  check_parameterize_preconditions(Adaptor& mesh);

    /// Initialize "A*X = B" linear system after
    /// (at least two) border vertices are parameterized.
    ///
    /// \pre Vertices must be indexed.
    /// \pre X and B must be allocated and empty.
    /// \pre At least 2 border vertices must be parameterized.
    void initialize_system_from_mesh_border(LeastSquaresSolver& solver,
                                            const Adaptor& mesh) ;

    /// Utility for setup_triangle_relations():
    /// Computes the coordinates of the vertices of a triangle
    /// in a local 2D orthonormal basis of the triangle's plane.
    void project_triangle(const Point_3& p0, const Point_3& p1, const Point_3& p2,  // in
                          Point_2& z0, Point_2& z1, Point_2& z2);                   // out

    /// Create two lines in the linear system per triangle (one for u, one for v).
    ///
    /// \pre vertices must be indexed.
    Error_code setup_triangle_relations(LeastSquaresSolver& solver,
                                        const Adaptor& mesh,
                                        Facet_const_handle facet) ;

    /// Copy X coordinates into the (u,v) pair of each vertex
    void set_mesh_uv_from_system(Adaptor& mesh,
                                 const LeastSquaresSolver& solver) ;

    /// Check parameterize() postconditions:
    /// - 3D -> 2D mapping is one-to-one.
    virtual Error_code check_parameterize_postconditions(const Adaptor& mesh,
                                                         const LeastSquaresSolver& solver);

    /// Check if 3D -> 2D mapping is one-to-one
    bool  is_one_to_one_mapping(const Adaptor& mesh,
                                 const LeastSquaresSolver& solver);

// Private accessors
private:
    /// Get the object that maps the surface's border onto a 2D space.
    Border_param&   get_border_parameterizer()    { return m_borderParameterizer; }

    /// Get the sparse linear algebra (traits object to access the linear system).
    Sparse_LA&      get_linear_algebra_traits() { return m_linearAlgebra; }

// Fields
private:
    /// Object that maps (at least two) border vertices onto a 2D space
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
//
// Implementation note: Outline of the algorithm:
// 1) Find an initial solution by projecting on a plane.
// 2) Lock two vertices of the mesh.
// 3) Copy the initial u,v coordinates to OpenNL.
// 3) Construct the LSCM equation with OpenNL.
// 4) Solve the equation with OpenNL.
// 5) Copy OpenNL solution to the u,v coordinates.
template<class Adaptor, class Border_param, class Sparse_LA>
inline
typename LSCM_parameterizer_3<Adaptor, Border_param, Sparse_LA>::Error_code
LSCM_parameterizer_3<Adaptor, Border_param, Sparse_LA>::
parameterize(Adaptor& mesh)
{
#ifdef DEBUG_TRACE
    // Create timer for traces
    CGAL::Timer timer;
    timer.start();
#endif

    // Check preconditions
    Error_code status = check_parameterize_preconditions(mesh);
#ifdef DEBUG_TRACE
    std::cerr << "  parameterization preconditions: " << timer.time() << " seconds." << std::endl;
    timer.reset();
#endif
    if (status != Base::OK)
        return status;

    // Count vertices
    int nbVertices = mesh.count_mesh_vertices();

    // Index vertices from 0 to nbVertices-1
    mesh.index_mesh_vertices();

    // Mark all vertices as *not* "parameterized"
    Vertex_iterator vertexIt;
    for (vertexIt = mesh.mesh_vertices_begin();
        vertexIt != mesh.mesh_vertices_end();
        vertexIt++)
    {
        mesh.set_vertex_parameterized(vertexIt, false);
    }

    // Compute (u,v) for (at least two) border vertices
    // and mark them as "parameterized"
    status = get_border_parameterizer().parameterize_border(mesh);
#ifdef DEBUG_TRACE
    std::cerr << "  border vertices parameterization: " << timer.time() << " seconds." << std::endl;
    timer.reset();
#endif
    if (status != Base::OK)
        return status;

    // Create sparse linear system "A*X = B" of size 2*nbVertices x 2*nbVertices
    // (in fact, we need only 2 lines per triangle x 1 column per vertex)
    LeastSquaresSolver solver(2*nbVertices);
    solver.set_least_squares(true) ;

    // Initialize the "A*X = B" linear system after
    // (at least two) border vertices parameterization
    initialize_system_from_mesh_border(solver, mesh);

    // Fill the matrix for the other vertices
    solver.begin_system() ;
    for (Facet_iterator facetIt = mesh.mesh_facets_begin();
         facetIt != mesh.mesh_facets_end();
         facetIt++)
    {
        // Create two lines in the linear system per triangle (one for u, one for v)
        status = setup_triangle_relations(solver, mesh, facetIt);
            if (status != Base::OK)
            return status;
    }
    solver.end_system() ;
#ifdef DEBUG_TRACE
    std::cerr << "  matrix filling (" << 2*mesh.count_mesh_facets() << " x " << nbVertices << "): "
              << timer.time() << " seconds." << std::endl;
    timer.reset();
#endif

    // Solve the "A*X = B" linear system in the least squares sense
    if ( ! solver.solve() )
        status = Base::ERROR_CANNOT_SOLVE_LINEAR_SYSTEM;
#ifdef DEBUG_TRACE
    std::cerr << "  solving linear system: "
              << timer.time() << " seconds." << std::endl;
    timer.reset();
#endif
    if (status != Base::OK)
        return status;

    // Copy X coordinates into the (u,v) pair of each vertex
    set_mesh_uv_from_system(mesh, solver);
#ifdef DEBUG_TRACE
    std::cerr << "  copy computed UVs to mesh :"
              << timer.time() << " seconds." << std::endl;
    timer.reset();
#endif

    // Check postconditions
    status = check_parameterize_postconditions(mesh, solver);
#ifdef DEBUG_TRACE
    std::cerr << "  parameterization postconditions: " << timer.time() << " seconds." << std::endl;
#endif
    if (status != Base::OK)
        return status;

    return status;
}


// Check parameterize() preconditions:
// - `mesh` must be a surface with one connected component
// - `mesh` must be a triangular mesh
template<class Adaptor, class Border_param, class Sparse_LA>
inline
typename LSCM_parameterizer_3<Adaptor, Border_param, Sparse_LA>::Error_code
LSCM_parameterizer_3<Adaptor, Border_param, Sparse_LA>::
check_parameterize_preconditions(Adaptor& mesh)
{
    Error_code status = Base::OK;	    // returned value

    // Helper class to compute genus or count borders, vertices, ...
    typedef Parameterization_mesh_feature_extractor<Adaptor>
                                            Mesh_feature_extractor;
    Mesh_feature_extractor feature_extractor(mesh);

    // Check that mesh is not empty
    if (mesh.mesh_vertices_begin() == mesh.mesh_vertices_end())
        status = Base::ERROR_EMPTY_MESH;
    if (status != Base::OK)
        return status;

    // The whole surface parameterization package is restricted to triangular meshes
    status = mesh.is_mesh_triangular() ? Base::OK
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

    return status;
}

// Initialize "A*X = B" linear system after
// (at least two) border vertices are parameterized
//
// Preconditions:
// - Vertices must be indexed
// - X and B must be allocated and empty
// - At least 2 border vertices must be parameterized
template<class Adaptor, class Border_param, class Sparse_LA>
inline
void LSCM_parameterizer_3<Adaptor, Border_param, Sparse_LA>::
initialize_system_from_mesh_border(LeastSquaresSolver& solver,
                                   const Adaptor& mesh)
{
    for (Vertex_const_iterator it = mesh.mesh_vertices_begin();
        it != mesh.mesh_vertices_end();
        it++)
    {
        // Get vertex index in sparse linear system
        int index = mesh.get_vertex_index(it);

        // Get vertex (u,v) (meaningless if vertex is not parameterized)
        Point_2 uv = mesh.get_vertex_uv(it);

        // Write (u,v) in X (meaningless if vertex is not parameterized)
        // Note  : 2*index     --> u
        //         2*index + 1 --> v
        solver.variable(2*index    ).set_value(uv.x()) ;
        solver.variable(2*index + 1).set_value(uv.y()) ;

        // Copy (u,v) in B if vertex is parameterized
        if (mesh.is_vertex_parameterized(it)) {
            solver.variable(2*index    ).lock() ;
            solver.variable(2*index + 1).lock() ;
        }
    }
}

// Utility for setup_triangle_relations():
// Computes the coordinates of the vertices of a triangle
// in a local 2D orthonormal basis of the triangle's plane.
template<class Adaptor, class Border_param, class Sparse_LA>
inline
void
LSCM_parameterizer_3<Adaptor, Border_param, Sparse_LA>::
project_triangle(const Point_3& p0, const Point_3& p1, const Point_3& p2,   // in
                 Point_2& z0, Point_2& z1, Point_2& z2)                     // out
{
    Vector_3 X = p1 - p0 ;
    NT X_norm = std::sqrt(X*X);
    if (X_norm != 0.0)
        X = X / X_norm;

    Vector_3 Z = CGAL::cross_product(X, p2 - p0) ;
    NT Z_norm = std::sqrt(Z*Z);
    if (Z_norm != 0.0)
        Z = Z / Z_norm;

    Vector_3 Y = CGAL::cross_product(Z, X) ;

    const Point_3& O = p0 ;

    NT x0 = 0 ;
    NT y0 = 0 ;
    NT x1 = std::sqrt( (p1 - O)*(p1 - O) ) ;
    NT y1 = 0 ;
    NT x2 = (p2 - O) * X ;
    NT y2 = (p2 - O) * Y ;

    z0 = Point_2(x0,y0) ;
    z1 = Point_2(x1,y1) ;
    z2 = Point_2(x2,y2) ;
}


// Create two lines in the linear system per triangle (one for u, one for v)
//
// Precondition: vertices must be indexed
//
// Implementation note: LSCM equation is:
//       (Z1 - Z0)(U2 - U0) = (Z2 - Z0)(U1 - U0)
// where Uk = uk + i.v_k is the complex number corresponding to (u,v) coords
//       Zk = xk + i.yk is the complex number corresponding to local (x,y) coords
// cool: no divide with this expression; makes it more numerically stable
// in presence of degenerate triangles
template<class Adaptor, class Border_param, class Sparse_LA>
inline
typename LSCM_parameterizer_3<Adaptor, Border_param, Sparse_LA>::Error_code
LSCM_parameterizer_3<Adaptor, Border_param, Sparse_LA>::
setup_triangle_relations(LeastSquaresSolver& solver,
                         const Adaptor& mesh,
                         Facet_const_handle facet)
{
    // Get the 3 vertices of the triangle
    Vertex_const_handle v0, v1, v2;
    int vertexIndex = 0;
    Vertex_around_facet_const_circulator cir = mesh.facet_vertices_begin(facet),
                                         end = cir;
    CGAL_For_all(cir, end)
    {
        if (vertexIndex == 0)
            v0 = cir;
        else if (vertexIndex == 1)
            v1 = cir;
        else if (vertexIndex == 2)
            v2 = cir;

        vertexIndex++;
    }
    if (vertexIndex != 3)
        return Base::ERROR_NON_TRIANGULAR_MESH;

    // Get the vertices index
    int id0 = mesh.get_vertex_index(v0) ;
    int id1 = mesh.get_vertex_index(v1) ;
    int id2 = mesh.get_vertex_index(v2) ;

    // Get the vertices position
    const Point_3& p0 = mesh.get_vertex_position(v0) ;
    const Point_3& p1 = mesh.get_vertex_position(v1) ;
    const Point_3& p2 = mesh.get_vertex_position(v2) ;

    // Computes the coordinates of the vertices of a triangle
    // in a local 2D orthonormal basis of the triangle's plane.
    Point_2 z0,z1,z2 ;
    project_triangle(p0,p1,p2,  //in
                     z0,z1,z2); // out
    Vector_2 z01 = z1 - z0 ;
    Vector_2 z02 = z2 - z0 ;
    NT a = z01.x() ;
    NT b = z01.y() ;
    NT c = z02.x() ;
    NT d = z02.y() ;
    CGAL_surface_mesh_parameterization_assertion(b == 0.0) ;

    // Create two lines in the linear system per triangle (one for u, one for v)
    // LSCM equation is:
    //       (Z1 - Z0)(U2 - U0) = (Z2 - Z0)(U1 - U0)
    // where Uk = uk + i.v_k is the complex number corresponding to (u,v) coords
    //       Zk = xk + i.yk is the complex number corresponding to local (x,y) coords
    //
    // Note  : 2*index     --> u
    //         2*index + 1 --> v
    int u0_id = 2*id0     ;
    int v0_id = 2*id0 + 1 ;
    int u1_id = 2*id1     ;
    int v1_id = 2*id1 + 1 ;
    int u2_id = 2*id2     ;
    int v2_id = 2*id2 + 1 ;
    //
    // Real part
    // Note: b = 0
    solver.begin_row() ;
    solver.add_coefficient(u0_id, -a+c)  ;
    solver.add_coefficient(v0_id,  b-d)  ;
    solver.add_coefficient(u1_id,   -c)  ;
    solver.add_coefficient(v1_id,    d)  ;
    solver.add_coefficient(u2_id,    a) ;
    solver.end_row() ;
    //
    // Imaginary part
    // Note: b = 0
    solver.begin_row() ;
    solver.add_coefficient(u0_id, -b+d) ;
    solver.add_coefficient(v0_id, -a+c) ;
    solver.add_coefficient(u1_id,   -d) ;
    solver.add_coefficient(v1_id,   -c) ;
    solver.add_coefficient(v2_id,    a) ;
    solver.end_row() ;

    return Base::OK;
}

// Copy X coordinates into the (u,v) pair of each vertex
template<class Adaptor, class Border_param, class Sparse_LA>
inline
void LSCM_parameterizer_3<Adaptor, Border_param, Sparse_LA>::
set_mesh_uv_from_system(Adaptor& mesh,
                        const LeastSquaresSolver& solver)
{
    Vertex_iterator vertexIt;
    for (vertexIt = mesh.mesh_vertices_begin();
         vertexIt != mesh.mesh_vertices_end();
         vertexIt++)
    {
        int index = mesh.get_vertex_index(vertexIt);

        // Note  : 2*index     --> u
        //         2*index + 1 --> v
        NT u = solver.variable(2*index    ).value() ;
        NT v = solver.variable(2*index + 1).value() ;

        // Fill vertex (u,v) and mark it as "parameterized"
        mesh.set_vertex_uv(vertexIt, Point_2(u,v));
        mesh.set_vertex_parameterized(vertexIt, true);
    }
}

// Check parameterize() postconditions:
// - 3D -> 2D mapping is one-to-one.
template<class Adaptor, class Border_param, class Sparse_LA>
inline
typename LSCM_parameterizer_3<Adaptor, Border_param, Sparse_LA>::Error_code
LSCM_parameterizer_3<Adaptor, Border_param, Sparse_LA>::
check_parameterize_postconditions(const Adaptor& mesh,
                                  const LeastSquaresSolver& solver)
{
    Error_code status = Base::OK;

    // Check if 3D -> 2D mapping is one-to-one
    status = is_one_to_one_mapping(mesh, solver)
           ? Base::OK
           : Base::ERROR_NO_1_TO_1_MAPPING;
    if (status != Base::OK)
        return status;

    return status;
}

// Check if 3D -> 2D mapping is one-to-one.
template<class Adaptor, class Border_param, class Sparse_LA>
inline
bool LSCM_parameterizer_3<Adaptor, Border_param, Sparse_LA>::
is_one_to_one_mapping(const Adaptor& mesh,
                      const LeastSquaresSolver& )
{
    Vector_3    first_triangle_normal;

    for (Facet_const_iterator facetIt = mesh.mesh_facets_begin();
         facetIt != mesh.mesh_facets_end();
         facetIt++)
    {
        // Get 3 vertices of the facet
        Vertex_const_handle v0, v1, v2;
        int vertexIndex = 0;
        Vertex_around_facet_const_circulator cir = mesh.facet_vertices_begin(facetIt),
                                             end = cir;
        CGAL_For_all(cir, end)
        {
            if (vertexIndex == 0)
                v0 = cir;
            else if (vertexIndex == 1)
                v1 = cir;
            else if (vertexIndex == 2)
                v2 = cir;

            vertexIndex++;
        }
        CGAL_surface_mesh_parameterization_assertion(vertexIndex >= 3);

        // Get the 3 vertices position IN 2D
        Point_2 p0 = mesh.get_vertex_uv(v0) ;
        Point_2 p1 = mesh.get_vertex_uv(v1) ;
        Point_2 p2 = mesh.get_vertex_uv(v2) ;

        // Compute the facet normal
        Point_3 p0_3D(p0.x(), p0.y(), 0);
        Point_3 p1_3D(p1.x(), p1.y(), 0);
        Point_3 p2_3D(p2.x(), p2.y(), 0);
        Vector_3 v01_3D = p1_3D - p0_3D;
        Vector_3 v02_3D = p2_3D - p0_3D;
        Vector_3 normal = CGAL::cross_product(v01_3D, v02_3D);

        // Check that all normals are oriented the same way
        // => no 2D triangle is flipped
        if (cir == mesh.facet_vertices_begin(facetIt))
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


} //namespace CGAL

#endif //CGAL_LSCM_PARAMETERIZER_3_H
