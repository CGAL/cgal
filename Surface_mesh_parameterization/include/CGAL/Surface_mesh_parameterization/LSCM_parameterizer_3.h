// Copyright (c) 2005  INRIA (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Laurent Saboret, Pierre Alliez, Bruno Levy

#ifndef CGAL_SURFACE_MESH_PARAMETERIZATION_LSCM_PARAMETERIZER_3_H
#define CGAL_SURFACE_MESH_PARAMETERIZATION_LSCM_PARAMETERIZER_3_H

#include <CGAL/license/Surface_mesh_parameterization.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Surface_mesh_parameterization/internal/Containers_filler.h>
#include <CGAL/Surface_mesh_parameterization/internal/kernel_traits.h>
#include <CGAL/Surface_mesh_parameterization/internal/validity.h>

#include <CGAL/Surface_mesh_parameterization/Two_vertices_parameterizer_3.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>

#include <CGAL/circulator.h>
#include <CGAL/Default.h>

#ifdef CGAL_EIGEN3_ENABLED
#include <CGAL/Eigen_solver_traits.h>
#endif

#include <boost/iterator/function_output_iterator.hpp>

#include <unordered_set>
#include <vector>

/// \file LSCM_parameterizer_3.h

namespace CGAL {

namespace Surface_mesh_parameterization {

// ------------------------------------------------------------------------------------
// Declaration
// ------------------------------------------------------------------------------------

/// \ingroup  PkgSurfaceMeshParameterizationMethods
///
/// The class `LSCM_parameterizer_3` implements the
/// *Least Squares Conformal Maps (LSCM)* parameterization \cgalCite{cgal:lprm-lscm-02}.
///
/// This is a conformal parameterization, i.e. it attempts to preserve angles.
///
/// This is a free border parameterization. There is no need to map the border
/// of the surface onto a convex polygon (only two pinned vertices are needed
/// to ensure a unique solution), but a one-to-one mapping is *not* guaranteed.
///
/// \cgalModels{Parameterizer_3}
///
/// \tparam TriangleMesh_ must be a model of `FaceGraph`.
///
/// \tparam BorderParameterizer_ is a strategy to parameterize the surface border
///         and must be a model of `Parameterizer_3`.<br>
///         <b>%Default:</b>
/// \code
///   Two_vertices_parameterizer_3<TriangleMesh_>
/// \endcode
///
/// \tparam SolverTraits_ must be a model of `SparseLinearAlgebraTraits_d`.<br>
///         Note: We may use a symmetric definite positive solver because LSCM
///         solves the system in the least squares sense.<br>
///         <b>%Default:</b> If \ref thirdpartyEigen "Eigen" 3.1 (or greater) is available
///         and `CGAL_EIGEN3_ENABLED` is defined, then an overload of `Eigen_solver_traits`
///         is provided as default parameter:
/// \code
///   CGAL::Eigen_solver_traits<Eigen::SimplicialLDLT< Eigen::SparseMatrix<double> > >
/// \endcode
///
///
/// \sa `CGAL::Surface_mesh_parameterization::Two_vertices_parameterizer_3<TriangleMesh, BorderParameterizer, SolverTraits>`
///
template < class TriangleMesh_,
           class BorderParameterizer_ = Default,
           class SolverTraits_ = Default>
class LSCM_parameterizer_3
{
public:
#ifndef DOXYGEN_RUNNING
  typedef typename Default::Get<
    BorderParameterizer_,
    Two_vertices_parameterizer_3<TriangleMesh_> >::type  Border_parameterizer;

#if !defined(CGAL_EIGEN3_ENABLED)
  static_assert(!(std::is_same<SolverTraits_, Default>::value),
                            "Error: You must either provide 'SolverTraits_' or link CGAL with the Eigen library");
  #endif

  typedef typename Default::Get<
    SolverTraits_,
  #if defined(CGAL_EIGEN3_ENABLED)
    CGAL::Eigen_solver_traits<
            Eigen::SimplicialLDLT<Eigen_sparse_symmetric_matrix<double>::EigenType> >
  #else
    SolverTraits_ // no parameter provided, and Eigen is not enabled: so don't compile!
  #endif
  >::type                                                     Solver_traits;
#else
  /// The border parameterizer
  typedef Border_parameterizer_                               Border_parameterizer;

  /// Solver traits type
  typedef SolverTraits_                                       Solver_traits;
#endif

  /// Triangle mesh type
  typedef TriangleMesh_                                       Triangle_mesh;

  typedef TriangleMesh_                                       TriangleMesh;

  /// Mesh halfedge type
  typedef typename boost::graph_traits<Triangle_mesh>::halfedge_descriptor halfedge_descriptor;

// Private types
private:
  typedef typename boost::graph_traits<Triangle_mesh>::vertex_descriptor   vertex_descriptor;
  typedef typename boost::graph_traits<Triangle_mesh>::face_descriptor     face_descriptor;

  // Traits subtypes:
  typedef typename internal::Kernel_traits<Triangle_mesh>::PPM      PPM;
  typedef typename internal::Kernel_traits<Triangle_mesh>::Kernel   Kernel;
  typedef typename Kernel::FT                                       NT;
  typedef typename Kernel::Point_2                                  Point_2;
  typedef typename Kernel::Point_3                                  Point_3;
  typedef typename Kernel::Vector_2                                 Vector_2;
  typedef typename Kernel::Vector_3                                 Vector_3;

  // Solver traits subtypes:
  typedef typename Solver_traits::Vector                            Vector;
  typedef typename Solver_traits::Matrix                            Matrix;

// Fields
private:
  // %Object that maps (at least two) border vertices onto a 2D space
  Border_parameterizer m_borderParameterizer;

  // Traits object to solve a sparse linear system
  Solver_traits m_linearAlgebra;

// Public operations
public:
  /// Constructor
  LSCM_parameterizer_3(Border_parameterizer border_param = Border_parameterizer(),
                       ///< %Object that maps the surface's border to 2D space
                       Solver_traits sparse_la = Solver_traits())
                       ///< Traits object to access a sparse linear system
    : m_borderParameterizer(border_param), m_linearAlgebra(sparse_la)
  { }

  // Default copy constructor and operator =() are fine

// Private accessors
private:
  // Get the object that maps the surface's border onto a 2D space.
  Border_parameterizer& get_border_parameterizer() { return m_borderParameterizer; }

  // Get the sparse linear algebra (traits object to access the linear system).
  Solver_traits& get_linear_algebra_traits() { return m_linearAlgebra; }

public:
  /// returns whether the 3D -> 2D mapping is one-to-one.
  template <typename VertexUVMap>
  bool is_one_to_one_mapping(const Triangle_mesh& mesh,
                             halfedge_descriptor bhd,
                             const VertexUVMap uvmap) const
  {
    return internal::is_one_to_one_mapping(mesh, bhd, uvmap);
  }

  /// computes a one-to-one mapping from a triangular 3D surface mesh
  /// to a piece of the 2D space.
  /// The mapping is piecewise linear (linear in each triangle).
  /// The result is the `(u,v)` pair image of each vertex of the 3D surface.
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
  /// \param uvmap an instantiation of the class `VertexUVmap`.
  /// \param vimap an instantiation of the class `VertexIndexMap`.
  /// \param vpmap an instantiation of the class `VertexParameterizedMap`.
  ///
  /// \pre `mesh` must be a triangular mesh.
  /// \pre The vertices must be indexed (`vimap` must be initialized).
  ///
  template <typename VertexUVmap, typename VertexIndexMap, typename VertexParameterizedMap>
  Error_code parameterize(Triangle_mesh& mesh,
                          halfedge_descriptor bhd,
                          VertexUVmap uvmap,
                          VertexIndexMap vimap,
                          VertexParameterizedMap vpmap)
  {
    CGAL_precondition(is_valid_polygon_mesh(mesh));
    CGAL_precondition(is_triangle_mesh(mesh));
    CGAL_precondition(bhd != boost::graph_traits<Triangle_mesh>::null_halfedge() && is_border(bhd, mesh));

    // Fill containers
    std::unordered_set<vertex_descriptor> ccvertices;
    std::vector<face_descriptor> ccfaces;

    internal::Containers_filler<Triangle_mesh> fc(mesh, ccvertices, &ccfaces);
    Polygon_mesh_processing::connected_component(
                                      face(opposite(bhd, mesh), mesh),
                                      mesh,
                                      boost::make_function_output_iterator(fc));

    // Count vertices
    int nbVertices = static_cast<int>(ccvertices.size());

    if (ccvertices.empty() || ccfaces.empty())
      return ERROR_EMPTY_MESH;

    // Compute (u,v) for (at least two) border vertices
    // and mark them as "parameterized"
    Error_code status = get_border_parameterizer().parameterize(mesh, bhd, uvmap, vimap, vpmap);

    if(status != OK)
      return status;

    // Vertices that are already parameterized are not part of the system
    std::vector<int> fvi(nbVertices); // vimap to free vertex index
    int index = 0;
    for(vertex_descriptor vd : ccvertices) {
      if(get(vpmap,vd)) {
        fvi[get(vimap, vd)] = -1;
      } else {
        fvi[get(vimap, vd)] = index++;
      }
    }

    // Create sparse linear system "A*X = B"
    int nbFreeVertices = index;
    int nbVariables = 2 * nbFreeVertices;
    Matrix A(nbVariables, nbVariables); // the constant matrix used in the linear system A*X = B
    Vector X(nbVariables), B(nbVariables);

    for(face_descriptor fd : ccfaces) {
      // Create two lines in the linear system per triangle (one for u, one for v)
      status = setup_triangle_relations(mesh, fd, uvmap, vimap, vpmap, fvi, A, B);
      if (status != OK)
        return status;
    }

    // Solve the "A*X = B" linear system in the least squares sense
    double D;
    if(!get_linear_algebra_traits().linear_solver(A, B, X, D)) {
      std::cerr << "Could not solve linear system" << std::endl;
      status = ERROR_CANNOT_SOLVE_LINEAR_SYSTEM;
      return status;
    }

    // Copy X coordinates into the (u,v) pair of each vertex
    for(vertex_descriptor vd : ccvertices) {
      if(get(vpmap,vd))
        continue;

      int index = get(vimap,vd);
      NT u = X[2 * fvi[index]];
      NT v = X[2 * fvi[index] + 1];
      put(uvmap, vd, Point_2(u,v));
    }

    return status;
  }

// Private operations
private:
  // Utility for setup_triangle_relations():
  // Computes the coordinates of the vertices of a triangle
  // in a local 2D orthonormal basis of the triangle's plane.
  void project_triangle(const Point_3& p0, const Point_3& p1, const Point_3& p2, // in
                        Point_2& z0, Point_2& z1, Point_2& z2) const             // out
  {
    Vector_3 X = p1 - p0;
    NT X_norm = std::sqrt(X * X);
    if(X_norm != 0.0)
    X = X / X_norm;

    Vector_3 Z = CGAL::cross_product(X, p2 - p0);
    NT Z_norm = std::sqrt(Z * Z);
    if(Z_norm != 0.0)
    Z = Z / Z_norm;

    Vector_3 Y = CGAL::cross_product(Z, X);

    const Point_3& O = p0;

    NT x0 = 0;
    NT y0 = 0;
    NT x1 = std::sqrt( (p1 - O)*(p1 - O) );
    NT y1 = 0;
    NT x2 = (p2 - O) * X;
    NT y2 = (p2 - O) * Y;

    z0 = Point_2(x0, y0);
    z1 = Point_2(x1, y1);
    z2 = Point_2(x2, y2);
  }

  // Create two lines in the linear system per triangle (one for u, one for v).
  //
  // \pre vertices of `mesh` must be indexed.
  //
  // Implementation note: LSCM equation is:
  //       (Z1 - Z0)(U2 - U0) = (Z2 - Z0)(U1 - U0)
  // where Uk = uk + i.v_k is the complex number corresponding to (u,v) coords
  //       Zk = xk + i.yk is the complex number corresponding to local (x,y) coords
  // cool: no divide with this expression; makes it more numerically stable
  // in presence of degenerate triangles
  template <typename UVmap,
            typename VertexIndexMap,
            typename VertexParameterizedMap,
            typename FreeVertexIndices>
  Error_code setup_triangle_relations(const Triangle_mesh& mesh,
                                      face_descriptor facet,
                                      UVmap uvmap,
                                      VertexIndexMap vimap,
                                      VertexParameterizedMap vpmap,
                                      const FreeVertexIndices& fvi,
                                      Matrix& A,
                                      Vector& B) const
  {
    const PPM ppmap = get(vertex_point, mesh);

    // Get the 3 vertices of the triangle
    vertex_descriptor v0, v1, v2;
    halfedge_descriptor h0 = halfedge(facet, mesh);
    v0 = target(h0, mesh);
    halfedge_descriptor h1 = next(h0, mesh);
    v1 = target(h1, mesh);
    halfedge_descriptor h2 = next(h1, mesh);
    v2 = target(h2, mesh);

    // Get the vertices position
    const Point_3& p0 = get(ppmap, v0);
    const Point_3& p1 = get(ppmap, v1);
    const Point_3& p2 = get(ppmap, v2);

    // Computes the coordinates of the vertices of a triangle
    // in a local 2D orthonormal basis of the triangle's plane.
    Point_2 z0, z1, z2;
    project_triangle(p0, p1, p2, //in
                     z0, z1, z2); // out
    Vector_2 z01 = z1 - z0;
    Vector_2 z02 = z2 - z0;
    NT a = z01.x();
    NT b = z01.y();
    NT c = z02.x();
    NT d = z02.y();
    CGAL_assertion(b == 0.0);

    // Create two lines in the linear system per triangle (one for u, one for v)
    // LSCM equation is:
    //       (Z1 - Z0)(U2 - U0) = (Z2 - Z0)(U1 - U0)
    // where Uk = uk + i.v_k is the complex number corresponding to (u,v) coords
    //       Zk = xk + i.yk is the complex number corresponding to local (x,y) coords
    //
    // Note  : 2*index     --> u
    //         2*index + 1 --> v

    std::vector<double> vals;
    std::vector<unsigned int> ids;
    std::vector<double> fvals;
    std::vector<double> pvals;

    auto add_coefficient = [&](const vertex_descriptor v, const double a, const bool is_u)
    {
      if(get(vpmap, v)) {
        const Point_2& uv = get(uvmap, v);
        fvals.push_back(a);
        pvals.push_back(is_u ? uv.x() : uv.y());
      } else {
        const int index = 2 * fvi[get(vimap, v)] + !is_u;
        CGAL_assertion(index >= 0);
        vals.push_back(a);
        ids.push_back(index);
      }
    };

    auto add_u_coefficient = [&](const vertex_descriptor v, const double a) {
      return add_coefficient(v, a, true /*is u*/);
    };
    auto add_v_coefficient = [&](const vertex_descriptor v, const double a) {
      return add_coefficient(v, a, false /*is not u*/);
    };

    auto accumulate = [&]()
    {
      const std::size_t nf = vals.size();
      const std::size_t nl = fvals.size();

      for(std::size_t i=0; i<nf; ++i)
        for(std::size_t j=0; j<nf; ++j)
          A.add_coef(ids[i], ids[j], vals[i] * vals[j]);

      double s = 0.;
      for(std::size_t i=0; i<nl; ++i)
        s += fvals[i] * pvals[i];

      for(std::size_t i=0; i<nf; ++i)
        B[ids[i]] -= vals[i] * s;

      vals.clear();
      ids.clear();
      fvals.clear();
      pvals.clear();
    };

    // Real part
    // Note: b = 0
    add_u_coefficient(v0, -a+c);
    add_v_coefficient(v0,  b-d);
    add_u_coefficient(v1, -c);
    add_v_coefficient(v1,  d);
    add_u_coefficient(v2,  a);
    accumulate();

    // Imaginary part
    // Note: b = 0
    add_u_coefficient(v0, -b+d);
    add_v_coefficient(v0, -a+c);
    add_u_coefficient(v1, -d);
    add_v_coefficient(v1, -c);
    add_v_coefficient(v2,  a);
    accumulate();

    return OK;
  }
};

} // namespace Surface_mesh_parameterization

} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_SURFACE_MESH_PARAMETERIZATION_LSCM_PARAMETERIZER_3_H
