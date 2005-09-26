// Copyright (c) 2005  INRIA (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$
// $Name$
//
// Author(s)     : Laurent Saboret, Pierre Alliez


#ifndef CGAL_LSCM_PARAMETIZER_3_H
#define CGAL_LSCM_PARAMETIZER_3_H

#include <CGAL/circulator.h>
#include <OpenNL/linear_solver.h>

#include <CGAL/Parametizer_traits_3.h>
#include <CGAL/Two_vertices_parametizer_3.h>
#include <CGAL/parameterization_assertions.h>

CGAL_BEGIN_NAMESPACE


//
// Declaration
//

/// Class LSCM_parametizer_3
/// implements Least Square Conformal Maps parameterization (Levy et al).
/// No need to map the surface's border onto a convex polygon
/// but 1 to 1 mapping is NOT guaranteed.
/// This is a conformal parameterization, i.e. it attempts to preserve angles.
///
/// Concept: Model of the ParametizerTraits_3 concept.
///
/// Design pattern:
/// ParametizerTraits_3 models are Strategies (see [GOF95]): they implement
/// a strategy of surface parameterization for models of MeshAdaptor_3.

template
<
    class MeshAdaptor_3,              ///< 3D surface mesh
    class BorderParametizer_3         ///< Strategy to parameterize the surface border
                = Two_vertices_parametizer_3<MeshAdaptor_3>,
                                      ///< Class to parameterize 2 border vertices
    class SparseLinearAlgebraTraits_d ///< Traits class to solve a sparse linear system
                = OpenNL::SymmetricLinearSolverTraits<typename MeshAdaptor_3::NT>
                                      ///< Symmetric solver for solving a sparse linear
                                      ///< system in the least squares sense
>
class LSCM_parametizer_3
    : public Parametizer_traits_3<MeshAdaptor_3>
{
// Private types
private:

    /// Superclass
    typedef Parametizer_traits_3<MeshAdaptor_3>
                                            Base;

// Public types
public:
    // Export Mesh_Adaptor_3, BorderParametizer_3
    // and SparseLinearAlgebraTraits_d types
    typedef MeshAdaptor_3                   Adaptor;
    typedef typename Base::Error_code       Error_code;
    typedef typename Adaptor::NT            NT;
    typedef typename Adaptor::Facet_handle  Facet_handle;
    typedef typename Adaptor::Facet_const_handle
                                            Facet_const_handle;
    typedef typename Adaptor::Vertex_handle Vertex_handle;
    typedef typename Adaptor::Vertex_const_handle
                                            Vertex_const_handle;
    typedef typename Adaptor::Point_3       Point_3;
    typedef typename Adaptor::Point_2       Point_2;
    typedef typename Adaptor::Vector_3      Vector_3;
    typedef typename Adaptor::Vector_2      Vector_2;
    typedef typename Adaptor::Facet_iterator Facet_iterator;
    typedef typename Adaptor::Facet_const_iterator
                                            Facet_const_iterator;
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
    typedef BorderParametizer_3             Border_param;
    typedef SparseLinearAlgebraTraits_d     Sparse_LA;
    typedef typename Sparse_LA::Vector      Vector;
    typedef typename Sparse_LA::Matrix      Matrix;

// Public operations
public:
    /// Constructor
    LSCM_parametizer_3(Border_param border_param = Border_param(),
                                            ///< Object that maps the surface's border to 2D space
                       Sparse_LA sparse_la = Sparse_LA())
                                            ///< Traits object to access a sparse linear system
        : m_borderParametizer(border_param), m_linearAlgebra(sparse_la)
    {}

    // Default copy constructor and operator =() are fine

    /// Compute a 1 to 1 mapping from a triangular 3D surface 'mesh'
    /// to a piece of the 2D space.
    /// The mapping is linear by pieces (linear in each triangle).
    /// The result is the (u,v) pair image of each vertex of the 3D surface.
    ///
    /// Preconditions:
    /// - 'mesh' must be a surface with 1 connected component.
    /// - 'mesh' must be a triangular mesh.
    virtual Error_code  parameterize(Adaptor* mesh);

// Private types
private:
    typedef typename OpenNL::LinearSolver<Sparse_LA>
                                            LeastSquaresSolver ;

// Private operations
private:
    /// Check parameterize() preconditions:
    /// - 'mesh' must be a surface with 1 connected component.
    /// - 'mesh' must be a triangular mesh.
    virtual Error_code  check_parameterize_preconditions(Adaptor* mesh);

    /// Initialize "A*X = B" linear system after
    /// (at least 2) border vertices are parameterized.
    ///
    /// Preconditions:
    /// - vertices must be indexed.
    /// - X and B must be allocated and empty.
    /// - (at least 2) border vertices must be parameterized.
    void initialize_system_from_mesh_border(LeastSquaresSolver* solver,
                                            const Adaptor& mesh) ;

    /// Utility for setup_triangle_relations():
    /// Computes the coordinates of the vertices of a triangle
    /// in a local 2D orthonormal basis of the triangle's plane.
    void project_triangle(const Point_3& p0, const Point_3& p1, const Point_3& p2,
                          Point_2* z0, Point_2* z1, Point_2* z2) ;

    /// Create 2 lines in the linear system per triangle (1 for u, 1 for v).
    ///
    /// Preconditions:
    /// - vertices must be indexed.
    Error_code setup_triangle_relations(LeastSquaresSolver* solver,
                                        const Adaptor& mesh,
                                        Facet_const_handle facet) ;

    /// Copy X coordinates into the (u,v) pair of each vertex
    void set_mesh_uv_from_system(Adaptor* mesh,
                                 const LeastSquaresSolver& solver) ;

    /// Check parameterize() postconditions:
    /// - "A*X = B" system is solvable (in the least squares sense)
    ///    with a good conditioning.
    /// - 3D -> 2D mapping is 1 to 1.
    virtual Error_code check_parameterize_postconditions(const Adaptor& mesh,
                                                         const LeastSquaresSolver& solver);

    /// Check if 3D -> 2D mapping is 1 to 1
    bool  is_one_to_one_mapping(const Adaptor& mesh,
                                 const LeastSquaresSolver& solver);

// Fields
private:
    /// Object that maps (at least 2) border vertices onto a 2D space
    Border_param    m_borderParametizer;

    /// Traits object to solve a sparse linear system
    Sparse_LA       m_linearAlgebra;
};


//
// Implementation
//

/// Compute a 1 to 1 mapping from a triangular 3D surface 'mesh'
/// to a piece of the 2D space.
/// The mapping is linear by pieces (linear in each triangle).
/// The result is the (u,v) pair image of each vertex of the 3D surface.
///
/// Preconditions:
/// - 'mesh' must be a surface with 1 connected component.
/// - 'mesh' must be a triangular mesh.
///
/// Implementation note: Outline of the algorithm:
/// 1) Find an initial solution by projecting on a plane.
/// 2) Lock two vertices of the mesh.
/// 3) Copy the initial u,v coordinates to OpenNL.
/// 3) Construct the LSCM equation with OpenNL.
/// 4) Solve the equation with OpenNL.
/// 5) Copy OpenNL solution to the u,v coordinates.
template<class Adaptor, class Border_param, class Sparse_LA>
inline
typename Parametizer_traits_3<Adaptor>::Error_code
LSCM_parametizer_3<Adaptor, Border_param, Sparse_LA>::
parameterize(Adaptor* mesh)
{
    CGAL_parameterization_assertion(mesh != NULL);

    // Check preconditions
    Error_code status = check_parameterize_preconditions(mesh);
    if (status != Base::OK)
        return status;

    // Count vertices
    int nbVertices = mesh->count_mesh_vertices();

    // Index vertices from 0 to nbVertices-1
    mesh->index_mesh_vertices();

    // Create sparse linear system "A*X = B" of size 2*nbVertices x 2*nbVertices
    // (in fact, we need only 2 lines per triangle x 1 column per vertex)
    LeastSquaresSolver solver(2*nbVertices);
    solver.set_least_squares(true) ;

    // Mark all vertices as NOT "parameterized"
    for (Vertex_iterator vertexIt = mesh->mesh_vertices_begin();
        vertexIt != mesh->mesh_vertices_end();
        vertexIt++)
    {
        mesh->set_vertex_parameterized(vertexIt, false);
    }

    // Compute (u,v) for (at least 2) border vertices
    // and mark them as "parameterized"
    status = m_borderParametizer.parameterize_border(mesh);
    if (status != Base::OK)
        return status;

    // Initialize the "A*X = B" linear system after
    // (at least 2) border vertices parameterization
    initialize_system_from_mesh_border(&solver, *mesh);

    // Fill the matrix for the other vertices
    fprintf(stderr,"  matrix filling (%d x %d)...\n",
                      int(2*mesh->count_mesh_facets()),
                      nbVertices);
    solver.begin_system() ;
    for (Facet_iterator facetIt = mesh->mesh_facets_begin();
         facetIt != mesh->mesh_facets_end();
         facetIt++)
    {
        // Create 2 lines in the linear system per triangle (1 for u, 1 for v)
        status = setup_triangle_relations(&solver, *mesh, facetIt);
            if (status != Base::OK)
            return status;
    }
    solver.end_system() ;
    fprintf(stderr,"    matrix filling ok\n");

    // Solve the "A*X = B" linear system in the least squares sense
    std::cerr << "  solver..." << std::endl;
    if ( ! solver.solve() )
    {
        std::cerr << "  error ERROR_CANNOT_SOLVE_LINEAR_SYSTEM!" << std::endl;
        return Base::ERROR_CANNOT_SOLVE_LINEAR_SYSTEM;
    }
    std::cerr << "    solver ok" << std::endl;

    // Copy X coordinates into the (u,v) pair of each vertex
    set_mesh_uv_from_system(mesh, solver);

    // Check postconditions
    status = check_parameterize_postconditions(*mesh, solver);
    if (status != Base::OK)
        return status;

    return status;
}


/// Check parameterize() preconditions:
/// - 'mesh' must be a surface with 1 connected component
/// - 'mesh' must be a triangular mesh
template<class Adaptor, class Border_param, class Sparse_LA>
inline
typename Parametizer_traits_3<Adaptor>::Error_code
LSCM_parametizer_3<Adaptor, Border_param, Sparse_LA>::
check_parameterize_preconditions(Adaptor* mesh)
{
    Error_code status = Base::OK;			// returned value

    // Helper class to compute genus or count boundaries, vertices, ...
    typedef Mesh_adaptor_feature_extractor<Adaptor>
                                            Mesh_feature_extractor;
    Mesh_feature_extractor feature_extractor(mesh);

    // Allways check that mesh is not empty
    if (mesh->mesh_vertices_begin() == mesh->mesh_vertices_end())
        status = Base::ERROR_EMPTY_MESH;
    if (status != Base::OK) {
        std::cerr << "  error ERROR_EMPTY_MESH!" << std::endl;
        return status;
    }

    // The whole surface parameterization package is restricted to triangular meshes
    CGAL_parameterization_expensive_precondition_code(                         \
        status = mesh->is_mesh_triangular() ? Base::OK                         \
                                            : Base::ERROR_NON_TRIANGULAR_MESH; \
    );
    if (status != Base::OK) {
        std::cerr << "  error ERROR_NON_TRIANGULAR_MESH!" << std::endl;
        return status;
    }

    // The whole package is restricted to surfaces: genus = 0,
    // 1 connected component and at least 1 boundary
    CGAL_parameterization_expensive_precondition_code(                      \
        int genus = feature_extractor.get_genus();                          \
        int nb_boundaries = feature_extractor.get_nb_boundaries();          \
        int nb_components = feature_extractor.get_nb_connex_components();   \
        status = (genus == 0 && nb_boundaries >= 1 && nb_components == 1)   \
               ? Base::OK                                                   \
               : Base::ERROR_NO_SURFACE_MESH;                               \
    );
    if (status != Base::OK) {
        std::cerr << "  error ERROR_NO_SURFACE_MESH!" << std::endl;
        return status;
    }

    return status;
}

/// Initialize "A*X = B" linear system after
/// (at least 2) border vertices are parameterized
///
/// Preconditions:
/// - vertices must be indexed
/// - X and B must be allocated and empty
/// - (at least 2) border vertices must be parameterized
template<class Adaptor, class Border_param, class Sparse_LA>
inline
void LSCM_parametizer_3<Adaptor, Border_param, Sparse_LA>::
initialize_system_from_mesh_border(LeastSquaresSolver* solver,
                                   const Adaptor& mesh)
{
    CGAL_parameterization_assertion(solver != NULL);
    CGAL_parameterization_assertion(solver != NULL);

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
        solver->variable(2*index    ).set_value(uv.x()) ;
        solver->variable(2*index + 1).set_value(uv.y()) ;

        // Copy (u,v) in B if vertex is parameterized
        if (mesh.is_vertex_parameterized(it)) {
            solver->variable(2*index    ).lock() ;
            solver->variable(2*index + 1).lock() ;
        }
    }
}

/// Utility for setup_triangle_relations():
/// Computes the coordinates of the vertices of a triangle
/// in a local 2D orthonormal basis of the triangle's plane.
template<class Adaptor, class Border_param, class Sparse_LA>
inline
void
LSCM_parametizer_3<Adaptor, Border_param, Sparse_LA>::
project_triangle(const Point_3& p0, const Point_3& p1, const Point_3& p2,
                 Point_2* z0, Point_2* z1, Point_2* z2)
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

    *z0 = Point_2(x0,y0) ;
    *z1 = Point_2(x1,y1) ;
    *z2 = Point_2(x2,y2) ;
}


/// Create 2 lines in the linear system per triangle (1 for u, 1 for v)
///
/// Preconditions:
/// - vertices must be indexed
///
/// Implementation note: LSCM equation is:
///       (Z1 - Z0)(U2 - U0) = (Z2 - Z0)(U1 - U0)
/// where Uk = uk + i.vk is the complex number corresponding to (u,v) coords
///       Zk = xk + i.yk is the complex number corresponding to local (x,y) coords
/// cool: no divide with this expression; makes it more numerically stable
/// in presence of degenerate triangles
template<class Adaptor, class Border_param, class Sparse_LA>
inline
typename Parametizer_traits_3<Adaptor>::Error_code
LSCM_parametizer_3<Adaptor, Border_param, Sparse_LA>::
setup_triangle_relations(LeastSquaresSolver* solver,
                         const Adaptor& mesh,
                         Facet_const_handle facet)
{
    CGAL_parameterization_assertion(solver != NULL);

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
    {
        std::cerr << "  error ERROR_NON_TRIANGULAR_MESH!" << std::endl;
        return Base::ERROR_NON_TRIANGULAR_MESH;
    }

    // Get the vertices index
    int id0 = mesh.get_vertex_index(v0) ;
    int id1 = mesh.get_vertex_index(v1) ;
    int id2 = mesh.get_vertex_index(v2) ;

#ifdef DEBUG_TRACE
    fprintf(stderr,"    Fill line for triangle (#%d, #%d, #%d): \n", id0, id1, id2);
#endif

    // Get the vertices position
    const Point_3& p0 = mesh.get_vertex_position(v0) ;
    const Point_3& p1 = mesh.get_vertex_position(v1) ;
    const Point_3& p2 = mesh.get_vertex_position(v2) ;

    // Computes the coordinates of the vertices of a triangle
    // in a local 2D orthonormal basis of the triangle's plane.
    Point_2 z0,z1,z2 ;
    project_triangle(p0,p1,p2, &z0,&z1,&z2) ;
    Vector_2 z01 = z1 - z0 ;
    Vector_2 z02 = z2 - z0 ;
    NT a = z01.x() ;
    NT b = z01.y() ;
    NT c = z02.x() ;
    NT d = z02.y() ;
    CGAL_parameterization_assertion(b == 0.0) ;

#ifdef DEBUG_TRACE
    fprintf(stderr,"      a=%f b=%f c=%f d=%f\n", (float)a, (float)b, (float)c, (float)d);
#endif

    // Create 2 lines in the linear system per triangle (1 for u, 1 for v)
    // LSCM equation is:
    //       (Z1 - Z0)(U2 - U0) = (Z2 - Z0)(U1 - U0)
    // where Uk = uk + i.vk is the complex number corresponding to (u,v) coords
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
    solver->begin_row() ;
    solver->add_coefficient(u0_id, -a+c)  ;
    solver->add_coefficient(v0_id,  b-d)  ;
    solver->add_coefficient(u1_id,   -c)  ;
    solver->add_coefficient(v1_id,    d)  ;
    solver->add_coefficient(u2_id,    a) ;
    solver->end_row() ;
    //
    // Imaginary part
    // Note: b = 0
    solver->begin_row() ;
    solver->add_coefficient(u0_id, -b+d) ;
    solver->add_coefficient(v0_id, -a+c) ;
    solver->add_coefficient(u1_id,   -d) ;
    solver->add_coefficient(v1_id,   -c) ;
    solver->add_coefficient(v2_id,    a) ;
    solver->end_row() ;

#ifdef DEBUG_TRACE
    fprintf(stderr,"      ok\n");
#endif

    return Base::OK;
}

/// Copy X coordinates into the (u,v) pair of each vertex
template<class Adaptor, class Border_param, class Sparse_LA>
inline
void LSCM_parametizer_3<Adaptor, Border_param, Sparse_LA>::
set_mesh_uv_from_system(Adaptor* mesh,
                        const LeastSquaresSolver& solver)
{
    fprintf(stderr,"  copy computed UVs to mesh\n");
    Vertex_iterator vertexIt;
    for (vertexIt = mesh->mesh_vertices_begin();
         vertexIt != mesh->mesh_vertices_end();
         vertexIt++)
    {
        int index = mesh->get_vertex_index(vertexIt);

        // Note  : 2*index     --> u
        //         2*index + 1 --> v
        NT u = solver.variable(2*index    ).value() ;
        NT v = solver.variable(2*index + 1).value() ;

        // Fill vertex (u,v) and mark it as "parameterized"
        mesh->set_vertex_uv(vertexIt, Point_2(u,v));
        mesh->set_vertex_parameterized(vertexIt, true);
    }
    fprintf(stderr,"    ok\n");
}

/// Check parameterize() postconditions:
/// - "A*X = B" system is solvable (in the least squares sense)
///   with a good conditioning
/// - 3D -> 2D mapping is 1 to 1
template<class Adaptor, class Border_param, class Sparse_LA>
inline
typename Parametizer_traits_3<Adaptor>::Error_code
LSCM_parametizer_3<Adaptor, Border_param, Sparse_LA>::
check_parameterize_postconditions(const Adaptor& mesh,
                                  const LeastSquaresSolver& solver)
{
    Error_code status = Base::OK;

    /* LS 02/2005: commented out this section because OpenNL::LinearSolver
     *             does not provide a is_solvable() method
     *
    // Check if "A*Xu = Bu" and "A*Xv = Bv" systems
    // are solvable with a good conditioning
    CGAL_parameterization_expensive_postcondition_code(         \
        status = get_linear_algebra_traits().is_solvable(A, Bu) \
               ? Base::OK                                       \
               : Base::ERROR_BAD_MATRIX_CONDITIONING;           \
    );
    if (status != Base::OK) {
        std::cerr << "  error ERROR_BAD_MATRIX_CONDITIONING!" << std::endl;
        return status;
    }
    CGAL_parameterization_expensive_postcondition_code(         \
        status = get_linear_algebra_traits().is_solvable(A, Bv) \
               ? Base::OK                                       \
               : Base::ERROR_BAD_MATRIX_CONDITIONING;           \
    );
    if (status != Base::OK) {
        std::cerr << "  error ERROR_BAD_MATRIX_CONDITIONING!" << std::endl;
        return status;
    }
     */

    // Check if 3D -> 2D mapping is 1 to 1
    CGAL_parameterization_postcondition_code( 		\
        status = is_one_to_one_mapping(mesh, solver) 	\
               ? Base::OK                               \
               : Base::ERROR_NO_1_TO_1_MAPPING;         \
    );
    if (status != Base::OK) {
        std::cerr << "  error ERROR_NO_1_TO_1_MAPPING!" << std::endl;
        return status;
    }

    return status;
}

/// Check if 3D -> 2D mapping is 1 to 1
/// The default implementation checks each normal
template<class Adaptor, class Border_param, class Sparse_LA>
inline
bool LSCM_parametizer_3<Adaptor, Border_param, Sparse_LA>::
is_one_to_one_mapping(const Adaptor& mesh,
                      const LeastSquaresSolver& solver)
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
        CGAL_parameterization_assertion(vertexIndex >= 3);

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


CGAL_END_NAMESPACE

#endif //CGAL_LSCM_PARAMETIZER_3_H

