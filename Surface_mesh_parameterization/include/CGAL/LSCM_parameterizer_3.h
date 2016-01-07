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

#include <CGAL/license/Surface_mesh_parameterization.h>


#include <CGAL/circulator.h>
#include <CGAL/Timer.h>
#include <CGAL/OpenNL/linear_solver.h>

#ifdef CGAL_EIGEN3_ENABLED
#include <CGAL/Eigen_solver_traits.h>
#endif

#include <CGAL/Parameterizer_traits_3.h>
#include <CGAL/Two_vertices_parameterizer_3.h>
#include <CGAL/surface_mesh_parameterization_assertions.h>
#include <CGAL/Parameterization_mesh_feature_extractor.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>

#include <iostream>

/// \file LSCM_parameterizer_3.h

namespace CGAL {


// ------------------------------------------------------------------------------------
// Declaration
// ------------------------------------------------------------------------------------

/// \ingroup  PkgSurfaceParameterizationMethods
///
/// The class `LSCM_parameterizer_3` implements the
/// *Least Squares Conformal Maps (LSCM)* parameterization  \cgalCite{cgal:lprm-lscm-02}.
///
/// This is a conformal parameterization, i.e. it attempts to preserve angles.
///
/// This is a free border parameterization. No need to map the border of the surface
/// onto a convex polygon (only two pinned vertices are needed to ensure a
/// unique solution), but one-to-one mapping is *not* guaranteed.
///
/// Note that his parametrization method has no template parameter for changing the sparse solver.
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
#if defined(CGAL_EIGEN3_ENABLED) || defined(DOXYGEN_RUNNING)
  = Eigen_solver_traits<Eigen::SimplicialLDLT<Eigen_sparse_symmetric_matrix<double>::EigenType> >
#else
  = OpenNL::SymmetricLinearSolverTraits<typename ParameterizationMesh_3::NT>
#endif
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
    typedef ParameterizationMesh_3          TriangleMesh;
    /// @endcond

    /// Export BorderParameterizer_3 template parameter.
    typedef BorderParameterizer_3           Border_param;
  
// Private types
private:  
    typedef SparseLinearAlgebraTraits_d     Sparse_LA;

// Private types
private:

  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_iterator face_iterator;
  typedef typename boost::graph_traits<TriangleMesh>::vertex_iterator vertex_iterator;
 
  typedef CGAL::Vertex_around_target_circulator<TriangleMesh> vertex_around_target_circulator;
  typedef CGAL::Vertex_around_face_circulator<TriangleMesh> vertex_around_face_circulator;

    // Mesh_Adaptor_3 subtypes:
  typedef Parameterizer_traits_3<TriangleMesh> Traits;
    typedef typename Traits::NT            NT;
    typedef typename Traits::Point_2       Point_2;
    typedef typename Traits::Point_3       Point_3;
    typedef typename Traits::Vector_2      Vector_2;
    typedef typename Traits::Vector_3      Vector_3;
   
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
template <typename VertexUVmap, typename VertexIndexMap, typename VertexParameterizedMap>
Error_code  parameterize(TriangleMesh& tmesh,
                         halfedge_descriptor bhd,
                         VertexUVmap uvmap,
                         VertexIndexMap vimap,
                         VertexParameterizedMap vpm)
{
#ifdef DEBUG_TRACE
    // Create timer for traces
    CGAL::Timer timer;
    timer.start();
#endif

    // Check preconditions
    Error_code status = check_parameterize_preconditions(tmesh);
#ifdef DEBUG_TRACE
    std::cerr << "  parameterization preconditions: " << timer.time() << " seconds." << std::endl;
    timer.reset();
#endif
    if (status != Base::OK)
        return status;

    // Count vertices
    int nbVertices = num_vertices(tmesh);

    // Index vertices from 0 to nbVertices-1
    // TODO mesh.index_mesh_vertices();

    // Compute (u,v) for (at least two) border vertices
    // and mark them as "parameterized"
    status = get_border_parameterizer().parameterize_border(tmesh,bhd,uvmap,vpm);
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
    initialize_system_from_mesh_border(solver, tmesh, uvmap, vimap,vpm);

    // Fill the matrix for the other vertices
    solver.begin_system() ;
    std::vector<face_descriptor> ccfaces;
    CGAL::Polygon_mesh_processing::connected_component(face(opposite(bhd,tmesh),tmesh),
                                                       tmesh,
                                                       std::back_inserter(ccfaces));
    BOOST_FOREACH(face_descriptor fd, ccfaces) // faces(tmesh)
    {
        // Create two lines in the linear system per triangle (one for u, one for v)
      status = setup_triangle_relations(solver, tmesh, fd,vimap);
            if (status != Base::OK)
            return status;
    }
    solver.end_system() ;


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
    //set_mesh_uv_from_system(tmesh, solver, uvmap);

    BOOST_FOREACH(vertex_descriptor vd, vertices(tmesh))
    {
      int index = get(vimap,vd);
      NT u = solver.variable(2*index    ).value() ;
      NT v = solver.variable(2*index + 1).value() ;
      std::cerr << u << " " << v << std::endl;
      put(uvmap, vd, Point_2(u,v));
    }

#ifdef DEBUG_TRACE
    std::cerr << "  copy computed UVs to mesh :"
              << timer.time() << " seconds." << std::endl;
    timer.reset();
#endif

    // Check postconditions
    status = check_parameterize_postconditions(tmesh, solver);
#ifdef DEBUG_TRACE
    std::cerr << "  parameterization postconditions: " << timer.time() << " seconds." << std::endl;
#endif
    if (status != Base::OK)
        return status;

    return status;
}


// Private operations
private:
    /// Check parameterize() preconditions:
    /// - `mesh` must be a surface with one connected component.
    /// - `mesh` must be a triangular mesh.
    virtual Error_code  check_parameterize_preconditions(TriangleMesh& mesh);

    /// Initialize "A*X = B" linear system after
    /// (at least two) border vertices are parameterized.
    ///
    /// \pre Vertices must be indexed.
    /// \pre X and B must be allocated and empty.
    /// \pre At least 2 border vertices must be parameterized.
  template <typename UVmap, typename VertexIndexMap,  typename VertexParameterizedMap>
    void initialize_system_from_mesh_border(LeastSquaresSolver& solver,
                                            const TriangleMesh& tmesh,
                                            UVmap uvmap,
                                            VertexIndexMap vimap,
                                            VertexParameterizedMap vpm)
{ 
    BOOST_FOREACH(vertex_descriptor v, vertices(tmesh)){
        // Get vertex index in sparse linear system
      int index = get(vimap, v);

        // Get vertex (u,v) (meaningless if vertex is not parameterized)
      Point_2 uv = get(uvmap, v);
      // TODO: it is meaningless but must it be called for non-border vertices??
        // Write (u,v) in X (meaningless if vertex is not parameterized)
        // Note  : 2*index     --> u
        //         2*index + 1 --> v
        solver.variable(2*index    ).set_value(uv.x()) ;
        solver.variable(2*index + 1).set_value(uv.y()) ;

        // Copy (u,v) in B if vertex is parameterized
        if (get(vpm,v)) {
          std::cerr << "vertex is parameterized"<< std::endl;
          solver.variable(2*index    ).lock() ;
          solver.variable(2*index + 1).lock() ;
        }
    }
 }
    


    /// Utility for setup_triangle_relations():
    /// Computes the coordinates of the vertices of a triangle
    /// in a local 2D orthonormal basis of the triangle's plane.
    void project_triangle(const Point_3& p0, const Point_3& p1, const Point_3& p2,  // in
                          Point_2& z0, Point_2& z1, Point_2& z2);                   // out

    /// Create two lines in the linear system per triangle (one for u, one for v).
    ///
    /// \pre vertices must be indexed.
  template <typename HalfedgeAsVertexIndexMap >
    Error_code setup_triangle_relations(LeastSquaresSolver& solver,
                                        const TriangleMesh& mesh,
                                        face_descriptor facet,
                                        HalfedgeAsVertexIndexMap) ;

    /// Copy X coordinates into the (u,v) pair of each vertex
  template <typename HalfedgeUVmap>
    void set_mesh_uv_from_system(TriangleMesh& mesh,
                                 const LeastSquaresSolver& solver,
                                 HalfedgeUVmap uvmap) ;

    /// Check parameterize() postconditions:
    /// - 3D -> 2D mapping is one-to-one.
    Error_code check_parameterize_postconditions(const TriangleMesh& mesh,
                                                         const LeastSquaresSolver& solver);

    /// Check if 3D -> 2D mapping is one-to-one
    bool  is_one_to_one_mapping(const TriangleMesh& mesh,
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
/*
/template<class TriangleMesh, class Border_param, class Sparse_LA>
inline
typename LSCM_parameterizer_3<TriangleMesh, Border_param, Sparse_LA>::Error_code
LSCM_parameterizer_3<TriangleMesh, Border_param, Sparse_LA>::
parameterize(TriangleMesh& tmesh, halfedge_descriptor bhd, HalfedgeUVmap uvmap, HalfedgeAsVertexIndexMap hvimap)
*/

// Check parameterize() preconditions:
// - `mesh` must be a surface with one connected component
// - `mesh` must be a triangular mesh
template<class TriangleMesh, class Border_param, class Sparse_LA>
inline
typename LSCM_parameterizer_3<TriangleMesh, Border_param, Sparse_LA>::Error_code
LSCM_parameterizer_3<TriangleMesh, Border_param, Sparse_LA>::
check_parameterize_preconditions(TriangleMesh& tmesh)
{
    Error_code status = Base::OK;	    // returned value
#if 0
    // Helper class to compute genus or count borders, vertices, ...
    typedef Parameterization_mesh_feature_extractor<TriangleMesh>
                                            Mesh_feature_extractor;
    Mesh_feature_extractor feature_extractor(mesh);

    // Check that mesh is not empty
    vertex_iterator b,e;
    boost::tie(b,e) = vertices(tmesh);
    if (b == e)
        status = Base::ERROR_EMPTY_MESH;
    if (status != Base::OK)
        return status;

    // The whole surface parameterization package is restricted to triangular meshes
    status = is_triangle_mesh(tmesh) ? Base::OK
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
#endif 
    return status;
}


#if 0
// Initialize "A*X = B" linear system after
// (at least two) border vertices are parameterized
//
// Preconditions:
// - Vertices must be indexed
// - X and B must be allocated and empty
// - At least 2 border vertices must be parameterized
template<class TriangleMesh, class Border_param, class Sparse_LA>
inline
void LSCM_parameterizer_3<TriangleMesh, Border_param, Sparse_LA>::
initialize_system_from_mesh_border(LeastSquaresSolver& solver,
                                   const TriangleMesh& tmesh)
#endif

// Utility for setup_triangle_relations():
// Computes the coordinates of the vertices of a triangle
// in a local 2D orthonormal basis of the triangle's plane.
template<class TriangleMesh, class Border_param, class Sparse_LA>
inline
void
LSCM_parameterizer_3<TriangleMesh, Border_param, Sparse_LA>::
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
template<class TriangleMesh, class Border_param, class Sparse_LA>
  template <typename HalfedgeAsVertexIndexMap >
inline
typename LSCM_parameterizer_3<TriangleMesh, Border_param, Sparse_LA>::Error_code
LSCM_parameterizer_3<TriangleMesh, Border_param, Sparse_LA>::
setup_triangle_relations(LeastSquaresSolver& solver,
                         const TriangleMesh& tmesh,
                         face_descriptor facet,
                         HalfedgeAsVertexIndexMap hvimap)
{
    typedef typename boost::property_map<TriangleMesh, boost::vertex_point_t>::const_type PPmap;
    PPmap ppmap = get(vertex_point, tmesh);
    // Get the 3 vertices of the triangle
    vertex_descriptor v0, v1, v2;
    halfedge_descriptor h0 = halfedge(facet,tmesh);
    v0 = target(h0,tmesh);
    halfedge_descriptor h1 = next(h0,tmesh);
    v1 = target(h1,tmesh);
    halfedge_descriptor h2 = next(h1,tmesh);
    v2 = target(h2,tmesh);

    // Get the vertices index
    int id0 = get(hvimap,h0) ;
    int id1 = get(hvimap,h1) ;
    int id2 = get(hvimap,h2) ;

    // Get the vertices position
    const Point_3& p0 = get(ppmap,v0) ;
    const Point_3& p1 = get(ppmap,v1) ;
    const Point_3& p2 = get(ppmap,v2) ;

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

#if 0
// Copy X coordinates into the (u,v) pair of each vertex
template<class TriangleMesh, class Border_param, class Sparse_LA>
template <typename HalfedgeUVmap>
inline
void LSCM_parameterizer_3<TriangleMesh, Border_param, Sparse_LA>::
set_mesh_uv_from_system(TriangleMesh& tmesh,
                        const LeastSquaresSolver& solver,
                        HalfedgeUVmap uvmap)
{
    BOOST_FOREACH(halfedge_descriptor vd, vertices(tmesh))
    {
        int index = mesh.get_vertex_index(vd);

        // Note  : 2*index     --> u
        //         2*index + 1 --> v
        NT u = solver.variable(2*index    ).value() ;
        NT v = solver.variable(2*index + 1).value() ;

        // Fill vertex (u,v) and mark it as "parameterized"
        mesh.set_vertex_uv(vd, Point_2(u,v));
        mesh.set_vertex_parameterized(vd, true);
    }
}
#endif

// Check parameterize() postconditions:
// - 3D -> 2D mapping is one-to-one.
template<class TriangleMesh, class Border_param, class Sparse_LA>
inline
typename LSCM_parameterizer_3<TriangleMesh, Border_param, Sparse_LA>::Error_code
LSCM_parameterizer_3<TriangleMesh, Border_param, Sparse_LA>::
check_parameterize_postconditions(const TriangleMesh& mesh,
                                  const LeastSquaresSolver& solver)
{
    Error_code status = Base::OK;
#if 0
    // Check if 3D -> 2D mapping is one-to-one
    status = is_one_to_one_mapping(mesh, solver)
           ? Base::OK
           : Base::ERROR_NO_1_TO_1_MAPPING;
    if (status != Base::OK)
        return status;
#endif
    return status;
}

// Check if 3D -> 2D mapping is one-to-one.
template<class TriangleMesh, class Border_param, class Sparse_LA>
inline
bool LSCM_parameterizer_3<TriangleMesh, Border_param, Sparse_LA>::
is_one_to_one_mapping(const TriangleMesh& mesh,
                      const LeastSquaresSolver& )
{
    const TriangleMesh& tmesh = mesh.get_adapted_mesh();
    Vector_3    first_triangle_normal(0., 0., 0.);

    BOOST_FOREACH(face_descriptor fd, faces(tmesh))
    {
        // Get 3 vertices of the facet
        vertex_descriptor v0, v1, v2;
        int vertexIndex = 0;
        vertex_around_face_circulator cir(halfedge(fd,tmesh),tmesh), first(cir), end(cir);
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


} //namespace CGAL

#endif //CGAL_LSCM_PARAMETERIZER_3_H
