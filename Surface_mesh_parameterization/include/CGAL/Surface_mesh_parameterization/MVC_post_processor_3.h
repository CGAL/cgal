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

#ifndef CGAL_SURFACE_MESH_PARAMETERIZATION_MVC_POST_PROCESSOR_3_H
#define CGAL_SURFACE_MESH_PARAMETERIZATION_MVC_POST_PROCESSOR_3_H

#include <CGAL/license/Surface_mesh_parameterization.h>

#include <CGAL/Surface_mesh_parameterization/internal/Bool_property_map.h>
#include <CGAL/Surface_mesh_parameterization/internal/Containers_filler.h>
#include <CGAL/Surface_mesh_parameterization/internal/kernel_traits.h>

#include <CGAL/Surface_mesh_parameterization/Two_vertices_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/parameterize.h>

#include <CGAL/Constrained_triangulation_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>

#ifdef CGAL_EIGEN3_ENABLED
#include <CGAL/Eigen_solver_traits.h>
#endif

#include <CGAL/Default.h>

#include <boost/unordered_set.hpp>

#include <vector>
#include <fstream>
#include <iostream>

/// \file MVC_post_processor_3.h

// @todo Determine the proper name of this file
// @todo Handle non-simple boundary

namespace CGAL {

namespace Surface_mesh_parameterization {

// ------------------------------------------------------------------------------------
// Declaration
// ------------------------------------------------------------------------------------

///
/// The class `MVC_post_processor_3` implements
/// the *Free boundary linear Parameterization* algorithm \cgalCite{kami2005free}.
///
/// This parameterizer provides a post processing step to other parameterizers
/// that do not necessarily return a valid embedding. It is based on
/// the convexification of the initial (2D) parameterization and the resolution
/// of a linear system with coefficients based on Mean Value Coordinates.
///
/// \cgalModels `Parameterizer_3`
///
/// \tparam TriangleMesh_ must be a model of `FaceGraph`.
///
/// \tparam SolverTraits_ must be a model of `SparseLinearAlgebraTraits_d`.<br>
///         <b>%Default:</b> If \ref thirdpartyEigen "Eigen" 3.1 (or greater) is available
///         and `CGAL_EIGEN3_ENABLED` is defined, then an overload of `Eigen_solver_traits`
///         is provided as default parameter:
/// \code
///   CGAL::Eigen_solver_traits<
///           Eigen::BiCGSTAB<Eigen_sparse_matrix<double>::EigenType,
///                           Eigen::IncompleteLUT< double > > >
/// \endcode
///
/// \sa `CGAL::Surface_mesh_parameterization::ARAP_parameterizer_3<TriangleMesh, BorderParameterizer_, Solver_traits>`
///
template < class TriangleMesh_,
           class SolverTraits_ = Default>
class MVC_post_processor_3
{
public:
#ifndef DOXYGEN_RUNNING
  typedef typename Default::Get<
    SolverTraits_,
  #if defined(CGAL_EIGEN3_ENABLED)
    CGAL::Eigen_solver_traits<
      Eigen::BiCGSTAB<Eigen_sparse_matrix<double>::EigenType,
                      Eigen::IncompleteLUT<double> > >
  #else
    SolverTraits_ // no parameter provided, and Eigen is not enabled: so don't compile!
  #endif
  >::type                                                     Solver_traits;
#else
  /// Solver traits type
  typedef SolverTraits_                                       Solver_traits;
#endif

  /// Triangle mesh type
  typedef TriangleMesh_                                       Triangle_mesh;

  typedef TriangleMesh_                                       TriangleMesh;

// Private types
private:
  // This class
  typedef MVC_post_processor_3<Triangle_mesh, Solver_traits>  Self;

// Private types
private:
  typedef typename boost::graph_traits<Triangle_mesh>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::graph_traits<Triangle_mesh>::halfedge_descriptor  halfedge_descriptor;
  typedef typename boost::graph_traits<Triangle_mesh>::face_descriptor      face_descriptor;
  typedef typename boost::graph_traits<Triangle_mesh>::face_iterator        face_iterator;
  typedef typename boost::graph_traits<Triangle_mesh>::vertex_iterator      vertex_iterator;

  typedef boost::unordered_set<vertex_descriptor>       Vertex_set;
  typedef std::vector<face_descriptor>                  Faces_vector;

  // Traits subtypes:
  typedef typename internal::Kernel_traits<Triangle_mesh>::Kernel   Kernel;
  typedef typename Kernel::FT                                       NT;
  typedef typename Kernel::Point_2                                  Point_2;
  typedef typename Kernel::Vector_2                                 Vector_2;
  typedef typename Kernel::Segment_2                                Segment_2;

  // Solver traits subtypes:
  typedef typename Solver_traits::Vector                            Vector;
  typedef typename Solver_traits::Matrix                            Matrix;

  // Types used for the convexification of the mesh
    // Each triangulation vertex is associated its corresponding vertex_descriptor
  typedef CGAL::Triangulation_vertex_base_with_info_2<vertex_descriptor, Kernel>  Vb;
    // Each triangulation face is associated a color (inside/outside information)
  typedef CGAL::Triangulation_face_base_with_info_2<int, Kernel>                  Fb;
  typedef CGAL::Constrained_triangulation_face_base_2<Kernel, Fb>                 Cfb;
  typedef CGAL::Triangulation_data_structure_2<Vb, Cfb>                           TDS;
  typedef CGAL::No_constraint_intersection_requiring_constructions_tag            Itag;

    // Can choose either a triangulation or a Delaunay triangulation
  typedef CGAL::Constrained_triangulation_2<Kernel, TDS, Itag>                    CT;
//  typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel, TDS, Itag>           CT;

// Private fields
private:
  // Traits object to solve a sparse linear system
  Solver_traits m_linearAlgebra;

// Private accessors
private:
  // Get the sparse linear algebra (traits object to access the linear system).
  Solver_traits& get_linear_algebra_traits() { return m_linearAlgebra; }

// Private utility
  // Print the exterior faces of the constrained triangulation.
  template <typename CT>
  void output_ct_exterior_faces(const CT& ct) const
  {
    std::ofstream out("constrained_triangulation_exterior.txt");
    typename CT::Finite_faces_iterator fit = ct.finite_faces_begin(),
                                       fend = ct.finite_faces_end();
    for(; fit!=fend; ++fit) {
      if(fit->info() != -1) // only output exterior (finite) faces
        continue;

      out << "4 ";
      for(std::size_t i=0; i<4; ++i) {
        out << fit->vertex(i%3)->point() << " 0 ";
      }
      out << '\n';
    }
    out << std::endl;
  }

  // Copy the data from two vectors to the UVmap.
  template <typename VertexUVMap,
            typename VertexIndexMap>
  void assign_solution(const Vector& Xu,
                       const Vector& Xv,
                       const Vertex_set& vertices,
                       VertexUVMap uvmap,
                       const VertexIndexMap vimap) const
  {
    for(vertex_descriptor vd : vertices) {
      int index = get(vimap, vd);
      NT u = Xu(index);
      NT v = Xv(index);
      put(uvmap, vd, Point_2(u, v));
    }
  }

// Private operations
private:
  // Store the vertices and faces of the mesh in memory.
  void initialize_containers(const Triangle_mesh& mesh,
                             halfedge_descriptor bhd,
                             Vertex_set& vertices,
                             Faces_vector& faces) const
  {
    internal::Containers_filler<Triangle_mesh> fc(mesh, vertices, &faces);
    CGAL::Polygon_mesh_processing::connected_component(
                                      face(opposite(bhd, mesh), mesh),
                                      mesh,
                                      boost::make_function_output_iterator(fc));
  }

  // Checks whether the polygon's border is simple.
  template <typename VertexUVMap>
  bool is_polygon_simple(const Triangle_mesh& mesh,
                         halfedge_descriptor bhd,
                         const VertexUVMap uvmap) const
  {
    // @fixme unefficient: use sweep line algorithms instead of brute force

    for(halfedge_descriptor hd_1 : halfedges_around_face(bhd, mesh)) {
      for(halfedge_descriptor hd_2 : halfedges_around_face(bhd, mesh)) {
        if(hd_1 == hd_2 || // equality
           next(hd_1, mesh) == hd_2 || next(hd_2, mesh) == hd_1) // adjacency
          continue;

        if(CGAL::do_intersect(Segment_2(get(uvmap, source(hd_1, mesh)),
                                        get(uvmap, target(hd_1, mesh))),
                              Segment_2(get(uvmap, source(hd_2, mesh)),
                                        get(uvmap, target(hd_2, mesh))))) {
#ifdef CGAL_SMP_ARAP_DEBUG
          std::ofstream out("non-simple.txt"); // polygon lines
          out << "2 " << get(uvmap, source(hd_1, mesh)) << " 0 "
                      << get(uvmap, target(hd_1, mesh)) << " 0" << std::endl;
          out << "2 " << get(uvmap, source(hd_2, mesh)) << " 0 "
                      << get(uvmap, target(hd_2, mesh)) << " 0" << std::endl;
#endif
          return false;
        }
      }
    }
    return true;
  }

  // Spread the inside / outside coloring from a Face to its neighbors
  // depending on whether the common edge is constrained.
  template <typename CT>
  void spread(CT& ct,
              const typename CT::Face_handle fh) const
  {
    typedef typename CT::Edge           Edge;
    typedef typename CT::Face_handle    Face_handle;

    int fh_color = fh->info();
    CGAL_precondition(fh_color != 0); // fh must be colored

    // look at the three neighbors for potential further spreading
    for(int i=0; i<3; ++i)
    {
      // ignore infinite faces and faces that already have a color
      Face_handle neigh_fh = fh->neighbor(i);

      if(ct.is_infinite(neigh_fh)) // infinite
        continue;

      if(neigh_fh->info() != 0) // already colored
        continue;

      bool is_common_edge_constrained = ct.is_constrained(Edge(fh, i));

      // Only constrain the exterior faces of the ct; since we have started
      // from an exterior face, we must not go over any constrained edge
      if(is_common_edge_constrained)
        continue;
      else
        neigh_fh->info() = fh_color;

      spread(ct, neigh_fh);
    }
  }

  // Triangulate the convex hull of the border of the parameterization.
  template <typename CT,
            typename VertexUVMap>
  Error_code triangulate_convex_hull(const Triangle_mesh& mesh,
                                     halfedge_descriptor bhd,
                                     const VertexUVMap uvmap,
                                     CT& ct) const
  {
    // Build the constrained triangulation

    // Since the border is closed and we are interested in triangles that are outside
    // of the border, we actually only need to insert points on the border
    for(halfedge_descriptor hd : halfedges_around_face(bhd, mesh)) {
      vertex_descriptor s = source(hd, mesh);
      const Point_2& sp = get(uvmap, s);

      typename CT::Vertex_handle vh = ct.insert(sp);
      vh->info() = s;
    }

    // Insert constraints (the border)
    for(halfedge_descriptor hd : halfedges_around_face(bhd, mesh)) {
      vertex_descriptor s = source(hd, mesh), t = target(hd, mesh);
      const Point_2& sp = get(uvmap, s), tp = get(uvmap, t);

      ct.insert_constraint(sp, tp);
    }

#ifdef CGAL_SMP_ARAP_DEBUG
    std::ofstream out("constrained_triangulation.cgal");
    out << ct;
#endif

    return OK;
  }

  // Color the (finite) faces of the constrained triangulation with an outside (-1) tag
  template <typename CT>
  Error_code color_faces(CT& ct) const
  {
    typedef typename CT::Edge                               Edge;
    typedef typename CT::Face_handle                        Face_handle;

    // Initialize all colors to 0
    typename CT::Finite_faces_iterator fit = ct.finite_faces_begin(),
                                       fend = ct.finite_faces_end();
    for(; fit!=fend; ++fit)
      fit->info() = 0;

    // start from infinite faces and check if the neighboring finite face is
    // inside or outside 'mesh'. If it is outside, mark it and spread to other
    // neighboring exterior faces
    typename CT::Vertex_handle infinite_vertex = ct.infinite_vertex();
    typename CT::Face_circulator fc = ct.incident_faces(infinite_vertex), done(fc);
    do {
      CGAL_precondition(ct.is_infinite(fc));
      int index_of_inf_vertex = fc->index(infinite_vertex);

      Face_handle mirror_face = fc->neighbor(index_of_inf_vertex);
      if(mirror_face->info() != 0)
        continue;

      bool is_edge_constrained = ct.is_constrained(Edge(fc, index_of_inf_vertex));
      if(!is_edge_constrained) {
        // edge is not constrained so the face is part of the convex hull but
        // geometrically outside of 'mesh'
        mirror_face->info() = -1; // outside
        spread<CT>(ct, mirror_face);
      }
    } while(++fc != done);

#ifdef CGAL_SMP_ARAP_DEBUG
    // Output the exterior faces of the constrained triangulation
    output_ct_exterior_faces(ct);
#endif

    return OK;
  }

  //                                                      -> ->
  // Return angle (in radians) of of (P,Q,R) corner (i.e. QP,QR angle).
  double compute_angle_rad(const Point_2& P,
                           const Point_2& Q,
                           const Point_2& R) const
  {
    Vector_2 u = P - Q;
    Vector_2 v = R - Q;

    double angle = std::atan2(v.y(), v.x()) - std::atan2(u.y(), u.x());
    if(angle < 0)
      angle += 2 * CGAL_PI;

    return angle;
  }

  // Fix vertices that are on the convex hull.
  template <typename CT,
            typename VertexParameterizedMap>
  void fix_convex_hull_border(const CT& ct,
                              VertexParameterizedMap vpmap) const
  {
    // All the vertices incident to the infinite vertex are on the convex hull
    typename CT::Vertex_circulator vc = ct.incident_vertices(ct.infinite_vertex()),
                                   vend = vc;
    do{
      vertex_descriptor vd = vc->info();
      put(vpmap, vd, true);
    } while (++vc != vend);
  }

  NT compute_w_ij_mvc(const Point_2& pi, const Point_2& pj, const Point_2& pk) const
  {
    //                                                               ->     ->
    // Compute the angle (pj, pi, pk), the angle between the vectors ij and ik
    NT angle = compute_angle_rad(pj, pi, pk);

    // For flipped triangles, the connectivity is inversed and thus the angle
    // computed by the previous function is not the one we need. Instead,
    // we need the explementary angle.
    if(angle > CGAL_PI) { // flipped triangle
      angle = 2 * CGAL_PI - angle;
    }
    NT weight = std::tan(0.5 * angle);

    return weight;
  }

  void fill_linear_system_matrix_mvc_from_points(const Point_2& pi, int i,
                                                 const Point_2& pj, int j,
                                                 const Point_2& pk, int k,
                                                 Matrix& A) const
  {
    // For MVC, the entry of A(i,j) is - [ tan(gamma_ij/2) + tan(delta_ij)/2 ] / |ij|
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
    Vector_2 edge_ij = pi - pj;
    double len_ij = std::sqrt(edge_ij * edge_ij);
    CGAL_assertion(len_ij != 0.0); // two points are identical!
    NT w_ij = w_i_base / len_ij;
    A.add_coef(i, j, w_ij);

    // Set w_ik in matrix
    Vector_2 edge_ik = pi - pk;
    double len_ik = std::sqrt(edge_ik * edge_ik);
    CGAL_assertion(len_ik != 0.0); // two points are identical!
    NT w_ik = w_i_base / len_ik;
    A.add_coef(i, k, w_ik);

    // Add to w_ii (w_ii = - sum w_ij)
    NT w_ii = - w_ij - w_ik;
    A.add_coef(i, i, w_ii);
  }

  // Partially fill the matrix coefficients A(i,i), A(i,j) and A(i,k)
  // Precondition: i, j, and k are ordered counterclockwise
  template <typename CT,
            typename VertexUVMap,
            typename VertexIndexMap,
            typename VertexParameterizedMap>
  void fill_linear_system_matrix_mvc_from_ct_vertices(typename CT::Vertex_handle vh_i,
                                                      typename CT::Vertex_handle vh_j,
                                                      typename CT::Vertex_handle vh_k,
                                                      const VertexUVMap uvmap,
                                                      const VertexIndexMap vimap,
                                                      const VertexParameterizedMap vpmap,
                                                      Matrix& A) const
  {
    vertex_descriptor vd_i = vh_i->info();
    vertex_descriptor vd_j = vh_j->info();
    vertex_descriptor vd_k = vh_k->info();

    // Coordinates are grabbed from the uvmap
    const Point_2& pi = get(uvmap, vd_i);
    const Point_2& pj = get(uvmap, vd_j);
    const Point_2& pk = get(uvmap, vd_k);
    int i = get(vimap, vd_i);
    int j = get(vimap, vd_j);
    int k = get(vimap, vd_k);

    // if vh_i is fixed, there is nothing to do: A(i,i)=1 and A(i,j)=0 for j!=i
    if(get(vpmap, vd_i))
    {
      // @fixme unefficient: A(i,i) is written as many times as i has neighbors
      A.set_coef(i, i, 1);
      return;
    }

    // vh_i is NOT fixed
    fill_linear_system_matrix_mvc_from_points(pi, i, pj, j, pk, k, A);
  }

  // Add the corresponding coefficients in A for all the edges of the face fh
  template <typename CT,
            typename VertexUVMap,
            typename VertexIndexMap,
            typename VertexParameterizedMap>
  void fill_linear_system_matrix_mvc_from_ct_face(const CT& CGAL_precondition_code(ct),
                                                  typename CT::Finite_faces_iterator fh,
                                                  const VertexUVMap uvmap,
                                                  const VertexIndexMap vimap,
                                                  const VertexParameterizedMap vpmap,
                                                  Matrix& A) const
  {
    CGAL_precondition(!ct.is_infinite(fh));
    typedef typename CT::Vertex_handle                    Vertex_handle;

    // Doing it explicitely rather than a loop for clarity
    Vertex_handle vh0 = fh->vertex(0);
    Vertex_handle vh1 = fh->vertex(1);
    Vertex_handle vh2 = fh->vertex(2);

    fill_linear_system_matrix_mvc_from_ct_vertices<CT>(vh0, vh1, vh2,
                                                       uvmap, vimap, vpmap, A);
    fill_linear_system_matrix_mvc_from_ct_vertices<CT>(vh1, vh2, vh0,
                                                       uvmap, vimap, vpmap, A);
    fill_linear_system_matrix_mvc_from_ct_vertices<CT>(vh2, vh0, vh1,
                                                       uvmap, vimap, vpmap, A);
  }

  template <typename VertexUVMap,
            typename VertexIndexMap,
            typename VertexParameterizedMap>
  void fill_linear_system_matrix_mvc_from_mesh_halfedge(const Triangle_mesh& mesh,
                                                        halfedge_descriptor hd,
                                                        const VertexUVMap uvmap,
                                                        const VertexIndexMap vimap,
                                                        const VertexParameterizedMap vpmap,
                                                        Matrix& A) const
  {
    vertex_descriptor vd_i = target(hd, mesh);
    vertex_descriptor vd_j = target(next(hd, mesh), mesh);
    vertex_descriptor vd_k = source(hd, mesh);
    const Point_2& pi = get(uvmap, vd_i);
    const Point_2& pj = get(uvmap, vd_j);
    const Point_2& pk = get(uvmap, vd_k);
    int i = get(vimap, vd_i);
    int j = get(vimap, vd_j);
    int k = get(vimap, vd_k);

    // if vh_i is fixed, there is nothing to do: A(i,i)=1 and A(i,j)=0 for j!=i
    if(get(vpmap, vd_i))
    {
      // @fixme unefficient A(i,i) is written as many times as i has neighbors
      A.set_coef(i, i, 1);
      return;
    }

    // Below, vh_i is NOT fixed
    fill_linear_system_matrix_mvc_from_points(pi, i, pj, j, pk, k, A);
  }

  // Fill the matrix A in an MVC linear system with the face 'fd' of 'mesh'.
  template <typename VertexUVMap,
            typename VertexIndexMap,
            typename VertexParameterizedMap>
  void fill_linear_system_matrix_mvc_from_mesh_face(const Triangle_mesh& mesh,
                                                    face_descriptor fd,
                                                    const VertexUVMap uvmap,
                                                    const VertexIndexMap vimap,
                                                    const VertexParameterizedMap vpmap,
                                                    Matrix& A) const
  {
    halfedge_descriptor hd = halfedge(fd, mesh);

    fill_linear_system_matrix_mvc_from_mesh_halfedge(mesh, hd, uvmap, vimap, vpmap, A);
    fill_linear_system_matrix_mvc_from_mesh_halfedge(mesh, next(hd, mesh),
                                                     uvmap, vimap, vpmap, A);
    fill_linear_system_matrix_mvc_from_mesh_halfedge(mesh, prev(hd, mesh),
                                                     uvmap, vimap, vpmap, A);
  }

  // Compute the matrix A in the MVC linear system using the exterior faces
  // of the constrained triangulation 'ct' and the graph 'mesh'.
  template <typename CT,
            typename VertexUVMap,
            typename VertexIndexMap,
            typename VertexParameterizedMap>
  Error_code compute_mvc_matrix(const CT& ct,
                                const Triangle_mesh& mesh,
                                const Faces_vector& faces,
                                const VertexUVMap uvmap,
                                const VertexIndexMap vimap,
                                const VertexParameterizedMap vpmap,
                                Matrix& A) const
  {
    Error_code status = OK;

    // The constrained triangulation has only "real" faces outside of the border
    // of mesh.

    // Thus, coefficients will come from 'ct' for faces that are in the convex hull
    // but not in 'mesh' nor in the holes of 'mesh'.
    // The coefficients will come from 'mesh' for faces that are in the convex hull
    // but not in 'ct'.

    // Loop over the "outside" faces of ct
    typename CT::Finite_faces_iterator fit = ct.finite_faces_begin(),
                                       fend = ct.finite_faces_end();
    for(; fit!=fend; ++fit)
    {
      if(fit->info() != -1) // not an outside face
        continue;

      fill_linear_system_matrix_mvc_from_ct_face(ct, fit, uvmap, vimap, vpmap, A);
    }

    // Loop over the faces of 'mesh'
    for(face_descriptor fd : faces) {
      fill_linear_system_matrix_mvc_from_mesh_face(mesh, fd, uvmap, vimap, vpmap, A);
    }

    return status;
  }

  // Compute the right hand side of a MVC linear system.
  template <typename VertexUVMap,
            typename VertexIndexMap,
            typename VertexParameterizedMap>
  void compute_mvc_rhs(const Vertex_set& vertices,
                       const VertexUVMap uvmap,
                       const VertexIndexMap vimap,
                       const VertexParameterizedMap vpmap,
                       Vector& Bu, Vector& Bv) const
  {
    for(vertex_descriptor vd : vertices) {
      int index = get(vimap, vd);
      const Point_2& uv = get(uvmap, vd);
      if(!get(vpmap, vd)) { // not yet parameterized
        Bu[index] = 0.; // might not be needed
        Bv[index] = 0.;
      } else { // fixed vertices
        Bu[index] = uv.x();
        Bv[index] = uv.y();
      }
    }
  }

  // Solve the two linear systems A*Xu=Bu and A*Xv=Bv.
  Error_code solve_mvc(const Matrix& A,
                       const Vector& Bu, const Vector& Bv,
                       Vector& Xu, Vector& Xv)
  {
    Error_code status = OK;

    NT Du, Dv;
    if(!get_linear_algebra_traits().linear_solver(A, Bu, Xu, Du) ||
       !get_linear_algebra_traits().linear_solver(A, Bv, Xv, Dv)) {
      status = ERROR_CANNOT_SOLVE_LINEAR_SYSTEM;
    }

    return status;
  }

  // Color the faces with inside/outside information and fix the border.
  template <typename CT, typename VertexParameterizedMap>
  Error_code prepare_CT_for_parameterization(CT& ct,
                                             VertexParameterizedMap vpmap) const
  {
    Error_code status = OK;

    // Gather the finite faces of the CT that are outside the (main) border of 'mesh'
    status = color_faces(ct);
    if(status != OK)
      return status;

    // Gather the vertices that are on the border of the convex hull and will be fixed
    fix_convex_hull_border(ct, vpmap);

    return status;
  }

  // Run an MVC parameterization on the (2D) ARAP UV map and the convexified mesh.
  template <typename CT,
            typename VertexUVMap,
            typename VertexIndexMap,
            typename VertexParameterizedMap>
  Error_code parameterize_convex_hull_with_MVC(const Triangle_mesh& mesh,
                                               const Vertex_set& vertices,
                                               const Faces_vector& faces,
                                               const CT& ct,
                                               VertexUVMap uvmap,
                                               const VertexIndexMap vimap,
                                               const VertexParameterizedMap vpmap)
  {
    Error_code status = OK;

    // Matrices and vectors needed in the resolution
    int nbVertices = static_cast<int>(vertices.size());
    Matrix A(nbVertices, nbVertices);
    Vector Xu(nbVertices), Xv(nbVertices), Bu(nbVertices), Bv(nbVertices);

    // Compute the (constant) matrix A
    compute_mvc_matrix(ct, mesh, faces, uvmap, vimap, vpmap, A);

    // Compute the right hand side of the linear system
    compute_mvc_rhs(vertices, uvmap, vimap, vpmap, Bu, Bv);

    // Solve the linear system
    status = solve_mvc(A, Bu, Bv, Xu, Xv);
    if(status != OK)
      return status;

    CGAL_postcondition_code
    (
      // make sure that the constrained vertices have not been moved
      for(vertex_descriptor vd : vertices) {
        if(get(vpmap, vd)) {
          int index = get(vimap, vd);
          CGAL_postcondition(std::abs(Xu[index] - Bu[index] ) < 1e-10);
          CGAL_postcondition(std::abs(Xv[index] - Bv[index] ) < 1e-10);
        }
      }
    )

    // Assign the UV values
    assign_solution(Xu, Xv, vertices, uvmap, vimap);

    return OK;
  }

public:
  template <typename VertexUVMap,
            typename VertexIndexMap>
  Error_code parameterize(const Triangle_mesh& mesh,
                          const Vertex_set& vertices,
                          const Faces_vector& faces,
                          halfedge_descriptor bhd,
                          VertexUVMap uvmap,
                          const VertexIndexMap vimap)
  {
    // Check if the border forms a simple polygon.
    // Note that there can be self-intersections in other borders,
    // but it is irrelevant and potential holes do not matter.
    const bool is_param_border_simple = is_polygon_simple(mesh, bhd, uvmap);

    // Not sure how to handle non-simple yet @fixme
    if(!is_param_border_simple) {
      std::cerr << "Border is not simple!" << std::endl;
      return ERROR_NON_CONVEX_BORDER;
    }

    // Compute the convex hull of the border of 'mesh'
    CT ct;
    triangulate_convex_hull<CT>(mesh, bhd, uvmap, ct);

    // Prepare the constrained triangulation: collect exterior faces (faces in
    // the convex hull but not -- geometrically -- in 'mesh').
    boost::unordered_set<vertex_descriptor> vs;
    internal::Bool_property_map<boost::unordered_set<vertex_descriptor> > vpmap(vs);
    prepare_CT_for_parameterization(ct, vpmap);

    // Run the MVC
    parameterize_convex_hull_with_MVC(mesh, vertices, faces, ct, uvmap, vimap, vpmap);

    return OK;
  }

  /// computes a one-to-one mapping from a triangular 2D surface mesh
  /// that is not necessarily embedded to a piece of the 2D space.
  ///
  /// \tparam VertexUVmap must be a model of `ReadWritePropertyMap` with
  ///         `boost::graph_traits<Triangle_mesh>::%vertex_descriptor` as key type and
  ///         %Point_2 (type deduced from `Triangle_mesh` using `Kernel_traits`)
  ///         as value type.
  /// \tparam VertexIndexMap must be a model of `ReadablePropertyMap` with
  ///         `boost::graph_traits<Triangle_mesh>::%vertex_descriptor` as key type and
  ///         a unique integer as value type.
  ///
  /// \param mesh a triangulated surface.
  /// \param bhd a halfedge descriptor on the boundary of `mesh`.
  /// \param uvmap an instanciation of the class `VertexUVmap`.
  /// \param vimap an instanciation of the class `VertexIndexMap`.
  ///
  template <typename VertexUVMap,
            typename VertexIndexMap>
  Error_code parameterize(const Triangle_mesh& mesh,
                          halfedge_descriptor bhd,
                          VertexUVMap uvmap,
                          const VertexIndexMap vimap)
  {
    Vertex_set vertices;
    Faces_vector faces;
    initialize_containers(mesh, bhd, vertices, faces);

    return parameterize(mesh, vertices, faces, bhd, uvmap, vimap);
  }

public:
  /// Constructor
  ///
  /// \param sparse_la %Traits object to access a sparse linear system.
  ///
  MVC_post_processor_3(Solver_traits sparse_la = Solver_traits())
    :
      m_linearAlgebra(sparse_la)
  { }
};

} // namespace Surface_mesh_parameterization

} // namespace CGAL

#endif // CGAL_SURFACE_MESH_PARAMETERIZATION_MVC_POST_PROCESSOR_3_H
