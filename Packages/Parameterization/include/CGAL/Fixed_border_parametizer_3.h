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


#ifndef CGAL_FIXED_BORDER_PARAMETIZER_3_H
#define CGAL_FIXED_BORDER_PARAMETIZER_3_H

#include <CGAL/circulator.h>
#include <OpenNL/linear_solver.h>

#include <CGAL/Parametizer_3.h>
#include <CGAL/Circular_border_parametizer_3.h>
#include <CGAL/Mesh_adaptor_feature_extractor.h>
#include <CGAL/parameterization_assertions.h>

CGAL_BEGIN_NAMESPACE


// ------------------------------------------------------------------------------------
// Declaration
// ------------------------------------------------------------------------------------

// Class Fixed_border_parametizer_3
// Model of the Parametizer_3 concept.
//
// Base class of fixed border parameterization methods (Tutte, Floater, ...)
// 1 to 1 mapping is guaranteed if surface's border is mapped onto a convex polygon.
//
// Implementation notes: - Subclasses must at least implement compute_wij()
//                       - The current implementation does not remove border vertices
//                         from the linear systems => A cannot be symmetric.
template
<
    class MeshAdaptor_3,              // 3D surface mesh
    class BorderParametizer_3         // Strategy to parameterize the surface border
                = Circular_border_arc_length_parametizer_3<MeshAdaptor_3>,
    class SparseLinearAlgebraTraits_d // Traits class to solve a sparse linear system
                = OpenNL::DefaultLinearSolverTraits<typename MeshAdaptor_3::NT>
>
class Fixed_border_parametizer_3
    : public Parametizer_3<MeshAdaptor_3>
{
// Private types
private:

    // Superclass
    typedef Parametizer_3<MeshAdaptor_3>    Base;

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
    // Constructor
    // @param border_param  Object that maps the surface's border to 2D space
    // @param sparse_la     Traits object to access a sparse linear system
    Fixed_border_parametizer_3(Border_param border_param = Border_param(),
                               Sparse_LA sparse_la = Sparse_LA())
        : m_borderParametizer(border_param), m_linearAlgebra(sparse_la)
    {}

    // Default copy constructor and operator =() are fine

    // Compute a 1 to 1 mapping from a triangular 3D surface 'mesh'
    // to a piece of the 2D space.
    // The mapping is linear by pieces (linear in each triangle).
    // The result is the (u,v) pair image of each vertex of the 3D surface.
    //
    // Preconditions:
    // * 'mesh' must be a surface with 1 connected component and no hole
    // * 'mesh' must be a triangular mesh
    // * the mesh border must be mapped onto a convex polygon
    virtual Error_code  parameterize(Adaptor* mesh);

// Protected operations
protected:
    // Check parameterize() preconditions:
    // * 'mesh' must be a surface with 1 connected component and no hole
    // * 'mesh' must be a triangular mesh
    // * the mesh border must be mapped onto a convex polygon
    virtual Error_code  check_parameterize_preconditions(Adaptor* mesh);

    // Initialize A, Bu and Bv after boundary parameterization
    // Fill the border vertices' lines in both linear systems:
    // "u = constant" and "v = constant"
    //
    // Preconditions:
    // * vertices must be indexed
    // * A, Bu and Bv must be allocated
    // * border vertices must be parameterized
    void  initialize_system_from_mesh_border (Matrix* A, Vector* Bu, Vector* Bv,
                                              const Adaptor& mesh);

    // compute wij = (i,j) coefficient of matrix A for j neighbor vertex of i
    // Implementation note: Subclasses must at least implement compute_wij()
    virtual NT  compute_wij(const Adaptor& mesh,
                            Vertex_const_handle main_vertex_Vi,
                            Vertex_around_vertex_const_circulator neighbor_vertex_Vj)
    = 0;

    // Compute the line i of matrix A for i inner vertex
    // * call compute_wij() to compute the A coefficient Wij for each neighbor Vj
    // * compute Wii = - sum of Wij
    //
    // Preconditions:
    // * vertices must be indexed
    // * vertex i musn't be already parameterized
    // * line i of A must contain only zeros
    virtual Error_code setup_inner_vertex_relations(Matrix* A,
                                                    Vector* Bu,
                                                    Vector* Bv,
                                                    const Adaptor& mesh,
                                                    Vertex_const_handle vertex);

    // Copy Xu and Xv coordinates into the (u,v) pair of each surface vertex
    void  set_mesh_uv_from_system (Adaptor* mesh,
                                   const Vector& Xu, const Vector& Xv);

    // Check parameterize() postconditions:
    // * "A*Xu = Bu" and "A*Xv = Bv" systems are solvable with a good conditioning
    // * 3D -> 2D mapping is 1 to 1
    virtual Error_code check_parameterize_postconditions(const Adaptor& mesh,
                                                         const Matrix& A,
                                                         const Vector& Bu,
                                                         const Vector& Bv);

    // Check if 3D -> 2D mapping is 1 to 1
    //
    // The default implementation checks each normal
    virtual bool  is_one_to_one_mapping(const Adaptor& mesh,
                                        const Matrix& A,
                                        const Vector& Bu,
                                        const Vector& Bv);

// Protected accessors
protected:
    // Get the object that maps the surface's border onto a 2D space
    Border_param&   get_border_parametizer()    { return m_borderParametizer; }

    // Get the sparse linear algebra (traits object to access the linear system)
    Sparse_LA&      get_linear_algebra_traits() { return m_linearAlgebra; }

// Fields
private:
    // Object that maps the surface's border onto a 2D space
    Border_param    m_borderParametizer;

    // Traits object to access the linear system
    Sparse_LA       m_linearAlgebra;
};


// ------------------------------------------------------------------------------------
// Implementation
// ------------------------------------------------------------------------------------

// Compute a 1 to 1 mapping from a triangular 3D surface 'mesh'
// to a piece of the 2D space.
// The mapping is linear by pieces (linear in each triangle).
// The result is the (u,v) pair image of each vertex of the 3D surface.
//
// Preconditions:
// * 'mesh' must be a surface with 1 connected component and no hole
// * 'mesh' must be a triangular mesh
// * the mesh border must be mapped onto a convex polygon
template<class Adaptor, class Border_param, class Sparse_LA>
inline
typename Parametizer_3<Adaptor>::Error_code
Fixed_border_parametizer_3<Adaptor, Border_param, Sparse_LA>::
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

    // Create 2 sparse linear systems "A*Xu = Bu" and "A*Xv = Bv" (1 line/column per vertex)
    Matrix A(nbVertices, nbVertices);
    Vector Xu(nbVertices), Xv(nbVertices), Bu(nbVertices), Bv(nbVertices);

    // Mark all vertices as NOT "parameterized"
    Vertex_iterator vertexIt;
    for (vertexIt = mesh->mesh_vertices_begin();
        vertexIt != mesh->mesh_vertices_end();
        vertexIt++)
    {
        mesh->set_vertex_parameterized(vertexIt, false);
    }

    // compute (u,v) for border vertices and mark them as "parameterized"
    status = get_border_parametizer().parameterize_border(mesh);
    if (status != Base::OK)
        return status;

    // Initialize A, Xu, Xv, Bu and Bv after boundary parameterization
    // Fill the border vertices' lines in both linear systems:
    // "u = constant" and "v = constant"
    //
    // Implementation note: the current implementation does not remove
    // border vertices from the linear systems => A cannot be symmetric
    initialize_system_from_mesh_border (&A, &Bu, &Bv, *mesh);

    // Fill the matrix for the inner vertices Vi: compute A's coefficient
    // Wij for each neighbor j; then Wii = - sum of Wij
    fprintf(stderr,"  matrix filling (%d x %d)...\n",nbVertices,nbVertices);
    for (vertexIt = mesh->mesh_vertices_begin();
         vertexIt != mesh->mesh_vertices_end();
         vertexIt++)
    {
        CGAL_parameterization_assertion(mesh->is_vertex_on_main_border(vertexIt)
                                     == mesh->is_vertex_parameterized(vertexIt));

        // inner vertices only
        if( ! mesh->is_vertex_on_main_border(vertexIt) )
        {
            // Compute the line i of matrix A for i inner vertex
            status = setup_inner_vertex_relations(&A, &Bu, &Bv,
                                                  *mesh,
                                                  vertexIt);
            if (status != Base::OK)
                return status;
        }
    }
    fprintf(stderr,"    matrix filling ok\n");

    // Solve "A*Xu = Bu". On success, solution is (1/Du) * Xu.
    // Solve "A*Xv = Bv". On success, solution is (1/Dv) * Xv.
    std::cerr << "  solver..." << std::endl;
    NT Du, Dv;
    if ( !get_linear_algebra_traits().linear_solver(A, Bu, Xu, Du) ||
         !get_linear_algebra_traits().linear_solver(A, Bv, Xv, Dv) )
    {
        std::cerr << "    error ERROR_CANNOT_SOLVE_LINEAR_SYSTEM!" << std::endl;
        return Base::ERROR_CANNOT_SOLVE_LINEAR_SYSTEM;
    }
    // WARNING: this package does not support homogeneous coordinates!
    CGAL_parameterization_assertion(Du == 1.0);
    CGAL_parameterization_assertion(Dv == 1.0);
    std::cerr << "    solver ok" << std::endl;

    // Copy Xu and Xv coordinates into the (u,v) pair of each vertex
    set_mesh_uv_from_system (mesh, Xu, Xv);

    // Check postconditions
    status = check_parameterize_postconditions(*mesh, A, Bu, Bv);
    if (status != Base::OK)
        return status;

    return status;
}


// Check parameterize() preconditions:
// * 'mesh' must be a surface with 1 connected component and no hole
// * 'mesh' must be a triangular mesh
// * the mesh border must be mapped onto a convex polygon
template<class Adaptor, class Border_param, class Sparse_LA>
inline
typename Parametizer_3<Adaptor>::Error_code
Fixed_border_parametizer_3<Adaptor, Border_param, Sparse_LA>::
check_parameterize_preconditions(Adaptor* mesh)
{
    Error_code status = Base::OK;			// returned value

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

    // 1 to 1 mapping is guaranteed if all Wij coefficients are > 0 (for j vertex neighbor of i)
    // and if the surface boundary is mapped onto a 2D convex polygon
    CGAL_parameterization_expensive_precondition_code(          \
        status = get_border_parametizer().is_border_convex()    \
               ? Base::OK                                       \
               : Base::ERROR_INVALID_BOUNDARY;                  \
    );
    if (status != Base::OK) {
        std::cerr << "  error ERROR_INVALID_BOUNDARY!" << std::endl;
        return status;
    }

    return status;
}

// Initialize A, Bu and Bv after boundary parameterization
// Fill the border vertices' lines in both linear systems: "u = constant" and "v = constant"
//
// Preconditions:
// * vertices must be indexed
// * A, Bu and Bv must be allocated
// * border vertices must be parameterized
template<class Adaptor, class Border_param, class Sparse_LA>
inline
void Fixed_border_parametizer_3<Adaptor, Border_param, Sparse_LA>::
initialize_system_from_mesh_border (Matrix* A, Vector* Bu, Vector* Bv,
                                    const Adaptor& mesh)
{
    CGAL_parameterization_assertion(A != NULL);
    CGAL_parameterization_assertion(Bu != NULL);
    CGAL_parameterization_assertion(Bv != NULL);

    for (Border_vertex_const_iterator it = mesh.mesh_main_border_vertices_begin();
         it != mesh.mesh_main_border_vertices_end();
         it++)
    {
        CGAL_parameterization_assertion(mesh.is_vertex_parameterized(it));

        // Get vertex index in sparse linear system
        int index = mesh.get_vertex_index(it);

        // Write 1 as diagonal coefficient of A
        A->set_coef(index, index, 1);

        // Write constant in Bu and Bv
        Point_2 uv = mesh.get_vertex_uv(it);
        (*Bu)[index] = uv.x();
        (*Bv)[index] = uv.y();
    }
}

// Compute the line i of matrix A for i inner vertex
// * call compute_wij() to compute the A coefficient Wij for each neighbor Vj
// * compute Wii = - sum of Wij
//
// Preconditions:
// * vertices must be indexed
// * vertex i musn't be already parameterized
// * line i of A must contain only zeros
template<class Adaptor, class Border_param, class Sparse_LA>
inline
typename Parametizer_3<Adaptor>::Error_code
Fixed_border_parametizer_3<Adaptor, Border_param, Sparse_LA>::
setup_inner_vertex_relations(Matrix* A,
                             Vector* Bu,
                             Vector* Bv,
                             const Adaptor& mesh,
                             Vertex_const_handle vertex)
{
    CGAL_parameterization_assertion( ! mesh.is_vertex_on_main_border(vertex) );
    CGAL_parameterization_assertion( ! mesh.is_vertex_parameterized(vertex) );

    int i = mesh.get_vertex_index(vertex);

#ifdef DEBUG_TRACE
    fprintf(stderr,"    Fill line #%d: \n", i);
#endif

    // circulate over vertices around 'vertex' to compute Wii and Wijs
    NT Wii = 0;
    int vertexIndex = 0;
    Vertex_around_vertex_const_circulator vj = mesh.vertices_around_vertex_begin(vertex);
    Vertex_around_vertex_const_circulator end = vj;
    CGAL_For_all(vj, end)
    {
        // Call to virtual method to do the actual coefficient computation
        NT Wij = -1.0 * compute_wij(mesh, vertex, vj);

        // Wii = - sum of Wij
        Wii -= Wij;

        // Get j index
        int j = mesh.get_vertex_index(vj);

        // Set Wij in matrix
        A->set_coef(i,j, Wij);

#ifdef DEBUG_TRACE
        fprintf(stderr,"      neighbor #%d -> wij = %5.2f\n", j, (float)Wij);
#endif

        vertexIndex++;
    }
    if (vertexIndex < 2)
    {
        std::cerr << "  error ERROR_NON_TRIANGULAR_MESH!" << std::endl;
        return Base::ERROR_NON_TRIANGULAR_MESH;
    }

#ifdef DEBUG_TRACE
    fprintf(stderr,"      ok\n");
#endif

    // Set Wii in matrix
    A->set_coef(i,i, Wii);

    return Base::OK;
}

// Copy Xu and Xv coordinates into the (u,v) pair of each surface vertex
template<class Adaptor, class Border_param, class Sparse_LA>
inline
void Fixed_border_parametizer_3<Adaptor, Border_param, Sparse_LA>::
set_mesh_uv_from_system(Adaptor* mesh,
                        const Vector& Xu, const Vector& Xv)
{
    fprintf(stderr,"  copy computed UVs to mesh\n");
    Vertex_iterator vertexIt;
    for (vertexIt = mesh->mesh_vertices_begin();
        vertexIt != mesh->mesh_vertices_end();
        vertexIt++)
    {
        int index = mesh->get_vertex_index(vertexIt);

        NT u = Xu[index];
        NT v = Xv[index];

        // Fill vertex (u,v) and mark it as "parameterized"
        mesh->set_vertex_uv(vertexIt, Point_2(u,v));
        mesh->set_vertex_parameterized(vertexIt, true);
    }
    fprintf(stderr,"    ok\n");
}

// Check parameterize() postconditions:
// * "A*Xu = Bu" and "A*Xv = Bv" systems are solvable with a good conditioning
// * 3D -> 2D mapping is 1 to 1
template<class Adaptor, class Border_param, class Sparse_LA>
inline
typename Parametizer_3<Adaptor>::Error_code
Fixed_border_parametizer_3<Adaptor, Border_param, Sparse_LA>::
check_parameterize_postconditions(const Adaptor& mesh,
                                  const Matrix& A,
                                  const Vector& Bu,
                                  const Vector& Bv)
{
    Error_code status = Base::OK;

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

    // Check if 3D -> 2D mapping is 1 to 1
    CGAL_parameterization_expensive_postcondition_code( \
        status = is_one_to_one_mapping(mesh, A, Bu, Bv) \
               ? Base::OK                               \
               : Base::ERROR_NO_1_TO_1_MAPPING;         \
    );
    if (status != Base::OK) {
        std::cerr << "  error ERROR_NO_1_TO_1_MAPPING!" << std::endl;
        //CGAL_parameterization_postcondition(false);
        return status;
    }

    return status;
}

// Check if 3D -> 2D mapping is 1 to 1
// The default implementation checks each normal
template<class Adaptor, class Border_param, class Sparse_LA>
inline
bool Fixed_border_parametizer_3<Adaptor, Border_param, Sparse_LA>::
is_one_to_one_mapping(const Adaptor& mesh,
                      const Matrix& A,
                      const Vector& Bu,
                      const Vector& Bv)
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

#endif //CGAL_FIXED_BORDER_PARAMETIZER_3_H

