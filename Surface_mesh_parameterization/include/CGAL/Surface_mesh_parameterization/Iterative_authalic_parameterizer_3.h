// Copyright (c) 2020 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

// This file has been adopted from CGAL, to integrate the below mentioned paper
//
// Paper         : Learning to Reconstruct Symmetric Shapes using Planar Parameterization of 3D Surface
// Author(s)     : Hardik Jain, Manuel Wöllhaf, Olaf Hellwich
// Conference    : IEEE International Conference on Computer Vision Workshops (ICCVW) 2019
//

#ifndef CGAL_SURFACE_MESH_PARAMETERIZATION_ITERATIVE_AUTHALIC_PARAMETERIZER_3_H
#define CGAL_SURFACE_MESH_PARAMETERIZATION_ITERATIVE_AUTHALIC_PARAMETERIZER_3_H

#include <CGAL/license/Surface_mesh_parameterization.h>

#include <CGAL/Surface_mesh_parameterization/internal/angles.h>
#include <CGAL/Surface_mesh_parameterization/internal/Bool_property_map.h>
#include <CGAL/Surface_mesh_parameterization/internal/Containers_filler.h>
#include <CGAL/Surface_mesh_parameterization/internal/kernel_traits.h>
#include <CGAL/Surface_mesh_parameterization/IO/File_off.h>
#include <CGAL/Surface_mesh_parameterization/Error_code.h>
#include <CGAL/Surface_mesh_parameterization/Circular_border_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Mean_value_coordinates_parameterizer_3.h>

#include <CGAL/Default.h>
#include <CGAL/Dynamic_property_map.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/Weights.h>
#include <CGAL/number_type_config.h>

#if defined(CGAL_EIGEN3_ENABLED)
#include <CGAL/Eigen_solver_traits.h>
#endif

#include <boost/function_output_iterator.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <unordered_set>

#define DEBUG_L0 1 // @fixme

/// \file Iterative_authalic_parameterizer_3.h

namespace CGAL {
namespace Surface_mesh_parameterization {

/// \ingroup  PkgSurfaceMeshParameterizationMethods
///
/// The class `Iterative_authalic_parameterizer_3` implements the *Iterative Parameterization* algorithm,
/// as described by Jain et al. \cgalCite{cgal:j-lrsspp-19}.
///
/// This parameterization is a fixed border parameterization and is part of the authalic
/// parameterization family, meaning that it aims to minimize area distortion
/// between the input surface mesh and the parameterized output.
/// More precisely, the approach used by this parameterizer is to iteratively redistribute
/// the \f$ L_2\f$ stretch - as defined by Sander et al. \cgalCite{cgal:ssgh-tmpm-01} - over the mesh.
///
/// \tparam TriangleMesh_ must be a model of `FaceGraph`.
/// \tparam BorderParameterizer_ is a Strategy to parameterize the surface border
///         and must be a model of `Parameterizer_3`.<br>
///         <b>%Default:</b>
/// \code
///   Circular_border_arc_length_parameterizer_3<TriangleMesh_>
/// \endcode
///
/// \tparam SolverTraits_ must be a model of `SparseLinearAlgebraTraits_d`.<br>
///         Note that the system is *not* symmetric because border vertices are not removed from the system.<br>
///         <b>%Default:</b> If \ref thirdpartyEigen "Eigen" 3.1 (or greater) is available
///         and `CGAL_EIGEN3_ENABLED` is defined, then an overload of `Eigen_solver_traits`
///         is provided as default parameter:
/// \code
///   CGAL::Eigen_solver_traits<
///           Eigen::BiCGSTAB<Eigen_sparse_matrix<double>::EigenType,
///                           Eigen::IncompleteLUT< double > > >
/// \endcode
///
/// \sa `CGAL::Surface_mesh_parameterization::Discrete_authalic_parameterizer_3<TriangleMesh, BorderParameterizer, SolverTraits>`
/// \sa `CGAL::Surface_mesh_parameterization::ARAP_parameterizer_3<TriangleMesh, BorderParameterizer, SolverTraits>`
///
template <class TriangleMesh_,
          class BorderParameterizer_ = Default,
          class SolverTraits_ = Default>
class Iterative_authalic_parameterizer_3
{
public:
#ifndef DOXYGEN_RUNNING
  typedef typename Default::Get<BorderParameterizer_,
                                Circular_border_arc_length_parameterizer_3<TriangleMesh_> >::type  Border_parameterizer;

  typedef typename Default::Get<SolverTraits_,
#if defined(CGAL_EIGEN3_ENABLED)
                                CGAL::Eigen_solver_traits<
                                  Eigen::BiCGSTAB<Eigen_sparse_matrix<double>::EigenType,
                                                  Eigen::IncompleteLUT<double> > >
#else
                                SolverTraits_ // no parameter provided, and Eigen is not enabled: so don't compile!
#endif
  >::type                                                                        Solver_traits;
#else
  /// Border parameterizer type
  typedef Border_parameterizer_                                                  Border_parameterizer;

  /// Solver traits type
  typedef SolverTraits_                                                          Solver_traits;
#endif

  /// Triangle mesh type
  typedef TriangleMesh_                                                          Triangle_mesh;

  typedef TriangleMesh_                                                          TriangleMesh;

  // Private types
private:
  typedef typename boost::graph_traits<Triangle_mesh>::vertex_descriptor         vertex_descriptor;
  typedef typename boost::graph_traits<Triangle_mesh>::halfedge_descriptor       halfedge_descriptor;
  typedef typename boost::graph_traits<Triangle_mesh>::face_descriptor           face_descriptor;

  typedef CGAL::Vertex_around_target_circulator<Triangle_mesh>                   vertex_around_target_circulator;
  typedef CGAL::Face_around_target_circulator<Triangle_mesh>                     face_around_target_circulator;

  typedef typename internal::Kernel_traits<Triangle_mesh>::Kernel                Kernel;

  typedef typename Kernel::FT                                                    NT;
  typedef typename Kernel::Point_2                                               Point_2;
  typedef typename Kernel::Point_3                                               Point_3;
  typedef typename Kernel::Vector_3                                              Vector_3;

  typedef typename internal::Kernel_traits<Triangle_mesh>::PPM                   PPM;
  typedef typename boost::property_traits<PPM>::reference                        PPM_ref;

  typedef std::unordered_set<vertex_descriptor>                                  Vertex_set;

  // Solver traits subtypes:
  typedef typename Solver_traits::Vector                                         Vector;
  typedef typename Solver_traits::Matrix                                         Matrix;

  typedef CGAL::dynamic_vertex_property_t<double>                                Vertex_double_tag;
  typedef typename boost::property_map<Triangle_mesh, Vertex_double_tag>::type   Vertex_Double_map;
  typedef CGAL::dynamic_vertex_property_t<bool>                                  Vertex_bool_tag;
  typedef typename boost::property_map<Triangle_mesh, Vertex_bool_tag>::type     Vertex_bool_map;
  typedef CGAL::dynamic_vertex_property_t<int>                                   Vertex_int_tag;
  typedef typename boost::property_map<Triangle_mesh, Vertex_int_tag>::type      Vertex_int_map;
  typedef CGAL::dynamic_vertex_property_t<Point_2>                               Vertex_point2_tag;
  typedef typename boost::property_map<Triangle_mesh, Vertex_point2_tag>::type   Vertex_point2_map;
  typedef CGAL::dynamic_face_property_t<NT>                                      Face_NT_tag;
  typedef typename boost::property_map<Triangle_mesh, Face_NT_tag>::type         Face_NT_map;

  // Fields
private:
  // Object that maps the surface's border onto a 2D space.
  Border_parameterizer m_border_parameterizer;

  // Traits object to solve a sparse linear system
  Solver_traits m_linear_algebra;

  // Object that keeps the last best UV Map
  Vertex_point2_map m_last_best_uv_map;

  // Counter to keep track of failure of Linear solver
  int m_linear_solver_failures;

  Vertex_Double_map m_vertex_L2_map;
  Face_NT_map m_face_L2_map;
  Face_NT_map m_face_areas;

  // Protected accessors
protected:
  /// Get the object that maps the surface's border onto a 2D space.
  Border_parameterizer& get_border_parameterizer() { return m_border_parameterizer; }

  /// Get the sparse linear algebra (traits object to access the linear system).
  Solver_traits& get_linear_algebra_traits() { return m_linear_algebra; }

  // Public operations
public:
  /// Constructor
  ///
  /// \param border_parameterizer %Object that maps the surface's border to 2D space
  /// \param sparse_la Traits object to access a sparse linear system
  ///
  Iterative_authalic_parameterizer_3(Border_parameterizer border_parameterizer = Border_parameterizer(),
                                     Solver_traits sparse_la = Solver_traits())
    :
      m_border_parameterizer(border_parameterizer),
      m_linear_algebra(sparse_la)
  { }

  // Disable copy constructor and operator =() due to property maps
  Iterative_authalic_parameterizer_3(const Iterative_authalic_parameterizer_3&) = delete;
  Iterative_authalic_parameterizer_3& operator=(const Iterative_authalic_parameterizer_3&) = delete;

  // Distortion functions
public:
  // Measure L2 stretch
  template <typename FaceRange, typename VertexUVmap>
  NT compute_area_distortion(const FaceRange& face_range,
                             const NT A_3D,
                             Triangle_mesh& tmesh,
                             const VertexUVmap uvmap)
  {
    Face_NT_map area_2DMap = get(Face_NT_tag(), tmesh);

    std::vector<NT> area_dist;
    NT A_2D = 0;

    for(face_descriptor f : face_range)
    {
      // get area in parameterised mesh
      const halfedge_descriptor h = halfedge(f, tmesh);
      const NT a_2D = abs(CGAL::area(get(uvmap, source(h, tmesh)),
                                     get(uvmap, target(h, tmesh)),
                                     get(uvmap, target(next(h, tmesh), tmesh))));
      put(area_2DMap, f, a_2D);
      A_2D += a_2D;
    }

    for(face_descriptor f : face_range)
    {
      const NT a_3D = get(m_face_areas, f);
      const NT a_2D = get(area_2DMap, f);

      area_dist.push_back(abs(a_3D/A_3D - a_2D/A_2D));
    }

    return std::accumulate(area_dist.begin(), area_dist.end(), NT(0));
  }

  // IO Helpers
public:
  template <typename VertexIndexMap>
  void print_matrix(const Vertex_set& vertices,
                    const VertexIndexMap vimap,
                    const Matrix& A,
                    const std::string name)
  {
    std::cout << "Matrix " << name << "(" << A.row_dimension() << "x" << A.column_dimension() << ")" << std::endl;

    Matrix A1(A.row_dimension(), A.column_dimension());
    int r=0, c=0;
    for(vertex_descriptor v1 : vertices)
    {
      int i = get(vimap, v1);
      for(vertex_descriptor v2 : vertices)
      {
        int j = get(vimap, v2);
        A1.set_coef(r, c, A.get_coef(i, j));
        ++c;
      }

      ++r;
      c = 0;
    }

    for(int r=0; r<A.row_dimension(); ++r)
    {
      for(int c=0; c<A.column_dimension(); ++c)
        std::cout << std::setw(10) << A1.get_coef(r, c) << "\t" << std::flush;
      std::cout << std::endl;
    }
  }

  template <typename VertexIndexMap>
  void print_vector(const Vertex_set& vertices,
                    const VertexIndexMap vimap,
                    const Vector& A,
                    const std::string name)
  {
    std::cout << "Vector " << name << "(" << A.size() << ")" << std::endl;
    Vector A1(A.size());
    int r = 0;
    for(vertex_descriptor v1 : vertices)
    {
      int i = get(vimap,v1);
      A1.set(r, A(i));
      ++r;
    }

    for(int r=0; r<A.size(); ++r)
      std::cout << A1(r) << std::endl;
  }

  void print(Matrix& A, Vector& Xu, Vector& Bu)
  {
    std::cout << "Matrix " << "(" << A.row_dimension() << "x" << A.column_dimension() << ")" << std::endl;
    for(int r=0; r<A.row_dimension(); ++r)
    {
      for(int c=0; c<A.column_dimension(); ++c)
        std::cout << std::setw(10) << A.get_coef(r, c) << "\t" << std::flush;
      std::cout << "\t\t"  << Xu(r) << "\t\t" << Bu(r) << std::endl;
    }
  }

  // Computation helpers
protected:
  // `operator=(onst Matrix& other)` isn't part of the concept...
  template <typename VertexIndexMap>
  void copy_sparse_matrix(const Matrix& src,
                          Matrix& dest,
                          const Triangle_mesh& tmesh,
                          const Vertex_set& vertices,
                          const VertexIndexMap vimap)
  {
    CGAL_precondition(src.row_dimension() == dest.row_dimension());
    CGAL_precondition(src.column_dimension() == dest.column_dimension());

    for(vertex_descriptor vertex : vertices)
    {
      const int i = get(vimap, vertex);
      vertex_around_target_circulator v_j(halfedge(vertex, tmesh), tmesh), end = v_j;
      CGAL_For_all(v_j, end)
      {
        const int j = get(vimap, *v_j);
        dest.set_coef(i, j, src.get_coef(i, j), false);
      }
    }
  }

private:
  double compute_vertex_L2(const Triangle_mesh& tmesh,
                           const vertex_descriptor v) const
  {
    NT phi = 0, local_area = 0;

    for(face_descriptor f : CGAL::faces_around_target(halfedge(v, tmesh), tmesh))
    {
      if(f == boost::graph_traits<Triangle_mesh>::null_face())
        continue;

      phi += CGAL::square(get(m_face_L2_map, f)) * get(m_face_areas, f);
      local_area += get(m_face_areas, f);
    }

    return sqrt(phi / local_area);
  }

  void compute_vertices_L2(const Vertex_set& vertices,
                           Triangle_mesh& tmesh)
  {
    m_vertex_L2_map = get(Vertex_double_tag(), tmesh);
    for(vertex_descriptor v : vertices)
    {
      put(m_vertex_L2_map, v, compute_vertex_L2(tmesh, v));
//      std::cout << "Vertex L2: " << v << " = " << compute_vertex_L2(tmesh, v) << std::endl;
    }
  }

  NT get_A(const std::array<const Point_2*, 3>& uv_points) const
  {
    NT A = (((uv_points[1]->x() - uv_points[0]->x()) * (uv_points[2]->y() - uv_points[0]->y()))
          - ((uv_points[2]->x() - uv_points[0]->x()) * (uv_points[1]->y() - uv_points[0]->y()))) / NT(2);

    CGAL_warning(A != NT(0)); // means degenerate face in the param space
    if(A == NT(0))
      return NT(1);

    return A;
  }

  Point_3 get_Ss(const std::array<const Point_3*, 3>& mesh_points,
                 const std::array<const Point_2*, 3>& uv_points,
                 const NT den) const
  {
    const NT dt0 = uv_points[1]->y() - uv_points[2]->y();
    const NT dt1 = uv_points[2]->y() - uv_points[0]->y();
    const NT dt2 = uv_points[0]->y() - uv_points[1]->y();
    Point_3 Ss(den * (mesh_points[0]->x()*dt0 + mesh_points[1]->x()*dt1 + mesh_points[2]->x()*dt2),
               den * (mesh_points[0]->y()*dt0 + mesh_points[1]->y()*dt1 + mesh_points[2]->y()*dt2),
               den * (mesh_points[0]->z()*dt0 + mesh_points[1]->z()*dt1 + mesh_points[2]->z()*dt2));
    return Ss;
  }

  Point_3 get_St(const std::array<const Point_3*, 3>& mesh_points,
                 const std::array<const Point_2*, 3>& uv_points,
                 const NT den) const
  {
    const NT ds0 = uv_points[2]->x() - uv_points[1]->x();
    const NT ds1 = uv_points[0]->x() - uv_points[2]->x();
    const NT ds2 = uv_points[1]->x() - uv_points[0]->x();
    Point_3 St(den * (mesh_points[0]->x()*ds0 + mesh_points[1]->x()*ds1 +mesh_points[2]->x()*ds2),
               den * (mesh_points[0]->y()*ds0 + mesh_points[1]->y()*ds1 +mesh_points[2]->y()*ds2),
               den * (mesh_points[0]->z()*ds0 + mesh_points[1]->z()*ds1 +mesh_points[2]->z()*ds2));
    return St;
  }

  NT compute_inner_product(const Point_3& pointA, const Point_3& pointB) const
  {
    return ((pointA.x())*(pointB.x()) + (pointA.y())*(pointB.y()) + (pointA.z())*(pointB.z()));
  }

  template <typename VertexUVMap>
  NT compute_face_L2(const face_descriptor f,
                     const Triangle_mesh& tmesh,
                     const VertexUVMap uvmap,
                     const PPM ppmap) const
  {
    std::array<const Point_2*, 3> uv_points;
    std::array<const Point_3*, 3> mesh_points;

    halfedge_descriptor h = halfedge(f, tmesh);
    for(std::size_t i=0; i<3; ++i)
    {
      vertex_descriptor v = target(h, tmesh);

      // just to be safe in case of weird VPM returning temporaries
      typename boost::property_traits<VertexUVMap>::reference uvp = get(uvmap, v);
      PPM_ref p = get(ppmap, v);

      uv_points[i] = &uvp;
      mesh_points[i] = &p;

      h = next(h, tmesh);
    }

    // Formula from Sander et al. 'Texture Mapping Progressive Meshes'
    const NT A = get_A(uv_points);
    const NT den = NT(1) / (NT(2) * A);
    const Point_3 Ss = get_Ss(mesh_points, uv_points, den);
    const Point_3 St = get_St(mesh_points, uv_points, den);

    const NT a = compute_inner_product(Ss, Ss);
    const NT c = compute_inner_product(St, St);

    return sqrt((a+c) / NT(2));
  }

  template <typename FaceRange, typename VertexUVMap>
  void compute_faces_L2(const FaceRange& face_range,
                        Triangle_mesh& tmesh,
                        const VertexUVMap uvmap)
  {
    m_face_L2_map = get(Face_NT_tag(), tmesh);
    const PPM ppmap = get(vertex_point, tmesh);

    for(face_descriptor f : face_range)
    {
      put(m_face_L2_map, f, compute_face_L2(f, tmesh, uvmap, ppmap));
//      std::cout << "Face L2: " << f << " = " << compute_face_L2(f, tmesh, uvmap, ppmap) << std::endl;
    }
  }

  template <typename FaceRange>
  NT initialize_faces_areas(const FaceRange& face_range,
                            Triangle_mesh& tmesh)
  {
    m_face_areas = get(Face_NT_tag(), tmesh);
    NT total_area = 0;

    for(face_descriptor f : face_range)
    {
      const NT f_area = Polygon_mesh_processing::face_area(f, tmesh);
      put(m_face_areas, f, f_area);
      total_area += f_area;
    }

    return total_area;
  }

  struct Neighbor_list
  {
    vertex_descriptor vertex;
    double angle;
    double length;
    Vector_3 vector;
    Point_2 uv;
    double weight;
  };

  NT determinant(Point_2& v0, Point_2& v1) const
  {
    return (v0.x() * v1.y() - v1.x() * v0.y());
  }

  void get_bary_coords(const Point_2& uv0, const Point_2& uv1, const Point_2& uv2,
                       NT& tau0, NT& tau1, NT& tau2) const
  {
    const NT det0 = determinant(uv1, uv2);
    const NT det1 = determinant(uv2, uv0);
    const NT det2 = determinant(uv0, uv1);
    const NT det3 = CGAL::determinant(Vector_2<Kernel>(uv1.x()-uv0.x(), uv1.y()-uv0.y()),
                                      Vector_2<Kernel>(uv2.x()-uv0.x(), uv2.y()-uv0.y()));
    CGAL_assertion(det3 > NT(0));
    if(det3 <= NT(0))
      det3 = NT(1);

    tau0 = det0 / det3;
    tau1 = det1 / det3;
    tau2 = det2 / det3;
  }

  double angle(Vector_3& v0, Vector_3& v1) const // @fixme use helper's?
  {
    return std::acos(v0*v1 / (CGAL::sqrt(v0*v0) * CGAL::sqrt(v1*v1)));
  }

  template <typename VertexIndexMap>
  Error_code setup_inner_vertex_relations(Matrix& A,
                                          Matrix& A_prev,
                                          Vector&,
                                          Vector&,
                                          const Triangle_mesh& tmesh,
                                          vertex_descriptor v,
                                          VertexIndexMap vimap)
  {
    const PPM ppmap = get(vertex_point, tmesh);

    // circulate over vertices around 'vertex' to compute w_ii and w_ijs
    std::vector<Neighbor_list> neighbor_list;
    int neighborsCounter = 0;
    double theta_sum = 0.;

    // create Neighbor_list vector with vertex and vector
    vertex_around_target_circulator v_j(halfedge(v, tmesh), tmesh), end_v_j = v_j;
    CGAL_For_all(v_j, end_v_j)
    {
      Neighbor_list NL;
      NL.vertex = *v_j;
      NL.vector = Vector_3(get(ppmap, v), tmesh.point(*v_j));
      NL.length = sqrt(NL.vector.squared_length());
      neighbor_list.push_back(NL);
      ++neighborsCounter;
    }

    if(neighborsCounter < 2)
      return ERROR_NON_TRIANGULAR_MESH;

    if(neighborsCounter == 2 && is_border(v, tmesh))
    {
      std::cout << "Encountered inner border with valency-2 vertex (" << v << "), "
                << "initializing with Tutte weights, this can affect optimization" << std::endl;

      // Tutte weights
      for(int k=0; k<neighborsCounter; ++k)
        neighbor_list[k].weight = 1.0;
    }
    else
    {
      for(int n=0; n<neighborsCounter; ++n)
      {
        int n_prev = (n==0 ? neighborsCounter-1 : n-1);
        double theta = angle(neighbor_list[n].vector, neighbor_list[n_prev].vector);
        neighbor_list[n].angle = theta;
        theta_sum += theta;
      }

      // Normalise the angle
      double factor = 2. / theta_sum;
      factor *= CGAL_PI;
      for(int n=0; n<neighborsCounter; ++n)
        neighbor_list[n].angle *= factor;

      neighbor_list[0].angle = 0.;
      for(int n=1; n<neighborsCounter; ++n)
        neighbor_list[n].angle += neighbor_list[n-1].angle;

      for(int n=0; n<neighborsCounter; ++n)
        neighbor_list[n].uv = Point_2(neighbor_list[n].length * cos(neighbor_list[n].angle),
                                      neighbor_list[n].length * sin(neighbor_list[n].angle));

      for(int j=0; j<neighborsCounter; ++j)
      {
        // Given the j-th neighbour of node i, find the two neighbours by intersecting the
        // line through nodes i and j with all segments of the polygon made by the neighbours.
        // Take the two neighbours on either side. Only one segment intersects this line.
        for(int k=0; k<neighborsCounter; ++k)
        {
          int kk = (k == neighborsCounter-1 ? 0 : k+1);
          if(k == j || kk == j)
            continue;

          NT cross1 = determinant(neighbor_list[j].uv, neighbor_list[k].uv);
          NT cross2 = determinant(neighbor_list[j].uv, neighbor_list[kk].uv);

          if(cross1 * cross2 <= NT(0))
          {
            NT tau0, tau1, tau2;
            get_bary_coords(neighbor_list[j].uv, neighbor_list[k].uv, neighbor_list[kk].uv, tau0, tau1, tau2);
            neighbor_list[j].weight += tau0;
            neighbor_list[k].weight += tau1;
            neighbor_list[kk].weight += tau2;
            break;
          }
        }
      }
    }

    // Scale the weights so that they sum to 1.
    double sum = 0.;
    for(int j=0; j<neighborsCounter; ++j)
      sum += neighbor_list[j].weight;

    double ratio = 1.0 / sum;
    for(int j=0; j<neighborsCounter; ++j)
      neighbor_list[j].weight *= ratio;

    // assign these weights to the edge pairs now
    NT w_ii = 0;
    const int i = get(vimap, v);
    for(int n=0; n<neighborsCounter; ++n)
    {
      NT w_ij = -1.0 * neighbor_list[n].weight;
      if(w_ij > 0)
        w_ij *= -1.0;
      w_ii -= w_ij;

      // Get j index
      const int j = get(vimap, neighbor_list[n].vertex);

      // Set w_ij in matrix
      A.set_coef(i, j, w_ij, true /*new*/);
      A_prev.set_coef(i, j, w_ij, true);
    }

    // Set w_ii in matrix
    A.set_coef(i,i, w_ii, true /*new*/);
    return OK;
  }

  /// computes `w_ij`, coefficient of matrix `A` for `j` neighbor vertex of `i`.
  ///
  /// \param mesh a triangulated surface.
  /// \param main_vertex_v_i the vertex of `mesh` with index `i`
  /// \param neighbor_vertex_v_j the vertex of `mesh` with index `j`
  NT compute_w_ij(const Triangle_mesh& tmesh,
                  vertex_descriptor main_vertex_v_i,
                  Vertex_around_target_circulator<Triangle_mesh> neighbor_vertex_v_j) const
  {
    const PPM ppmap = get(vertex_point, tmesh);

    const PPM_ref position_v_i = get(ppmap, main_vertex_v_i);
    const PPM_ref position_v_j = get(ppmap, *neighbor_vertex_v_j);

    // Compute the square norm of v_j -> v_i vector
    Vector_3 edge = position_v_i - position_v_j;
    NT square_len = edge*edge;

    // Compute cotangent of (v_k,v_j,v_i) corner (i.e. cotan of v_j corner)
    // if v_k is the vertex before v_j when circulating around v_i
    vertex_around_target_circulator previous_vertex_v_k = neighbor_vertex_v_j;
    --previous_vertex_v_k;
    const PPM_ref position_v_k = get(ppmap, *previous_vertex_v_k);
//    NT cotg_psi_ij = internal::cotangent<Kernel>(position_v_k, position_v_j, position_v_i);
    NT cotg_beta_ij = internal::cotangent<Kernel>(position_v_i, position_v_k, position_v_j);

    // Compute cotangent of (v_i,v_j,v_l) corner (i.e. cotan of v_j corner)
    // if v_l is the vertex after v_j when circulating around v_i
    vertex_around_target_circulator next_vertex_v_l = neighbor_vertex_v_j;
    ++next_vertex_v_l;

    const Point_3 position_v_l = get(ppmap, *next_vertex_v_l);
//    NT cotg_theta_ij = internal::cotangent<Kernel>(position_v_i, position_v_j, position_v_l);
    NT cotg_alpha_ij = internal::cotangent<Kernel>(position_v_j, position_v_l, position_v_i);

    NT weight = 0;
    CGAL_assertion(square_len > NT(0)); // two points are identical!
    if(square_len != NT(0))
      weight = cotg_beta_ij + cotg_alpha_ij;

    return weight;
  }

  template <typename VertexIndexMap>
  Error_code setup_inner_vertex_relations_cotangent(Matrix& A,
                                                    Vector&,
                                                    Vector&,
                                                    const Triangle_mesh& tmesh,
                                                    vertex_descriptor v,
                                                    VertexIndexMap& vimap)
  {
    const int i = get(vimap, v);

    // circulate over vertices around 'vertex' to compute w_ii and w_ijs
    NT w_ii = 0;

    bool use_uniform_weights = false;
    const int neighborsCounter = degree(v, tmesh);
    if(neighborsCounter < 2)
    {
      return ERROR_NON_TRIANGULAR_MESH;
    }
    else if(neighborsCounter == 2 && is_border(v, tmesh))
    {
      use_uniform_weights = true;
      std::cerr << "Encountered inner border with valency-2 vertex (" << get(vimap, v) << ") "
                << "initializing with uniform weights, this can affect optimization" << std::endl;
    }

    vertex_around_target_circulator vj(halfedge(v, tmesh), tmesh), end = vj;
    CGAL_For_all(vj, end)
    {
      NT w_ij = NT(-1);
      if(!use_uniform_weights)
        w_ij *= compute_w_ij(tmesh, v, vj);

      // w_ii = - sum of w_ijs
      w_ii -= w_ij;

      // Get j index
      const int j = get(vimap, *vj);
//      std::cout << "W[" << i << ", " << j << "] = " << w_ij << std::endl;

      // Set w_ij in matrix
      A.set_coef(i, j, w_ij, true /*new*/);
    }

    // Set w_ii in matrix
    A.set_coef(i, i, w_ii, true /*new*/);

    return OK;
  }

  // Initialize the UV values with a first parameterization of the input.
  template <typename VertexIndexMap>
  Error_code setup_inner_vertex_relations_MVC(Matrix& A,
                                              Vector&,
                                              Vector&,
                                              Triangle_mesh& tmesh,
                                              vertex_descriptor v,
                                              VertexIndexMap& vimap) const
  {
    auto vpm = get_const_property_map(CGAL::vertex_point, tmesh);
    CGAL::internal::Mean_value_weight<Triangle_mesh, decltype(vpm)> compute_mvc(tmesh, vpm);

    const int i = get(vimap, v);

    // circulate over vertices around 'vertex' to compute w_ii and w_ijs
    NT w_ii = 0;
    int vertexIndex = 0;

    for(halfedge_descriptor h : CGAL::halfedges_around_target(v, tmesh))
    {
      NT w_ij = NT(-1) * compute_mvc(h);
      // w_ii = - sum of w_ijs
      w_ii -= w_ij;

      // Get j index
      int j = get(vimap, source(h, tmesh));

      // Set w_ij in matrix
      A.set_coef(i, j, w_ij, true /*new*/);
      vertexIndex++;
    }

    if (vertexIndex < 2)
      return ERROR_NON_TRIANGULAR_MESH;

    // Set w_ii in matrix
    A.set_coef(i,i, w_ii, true /*new*/);
    return OK;
  }

  template <typename VertexIndexMap>
  Error_code setup_iter_inner_vertex_relations(Matrix& A,
                                               Matrix& A_prev,
                                               Vector&,
                                               Vector&,
                                               const Triangle_mesh& tmesh,
                                               vertex_descriptor v,
                                               VertexIndexMap vimap,
                                               NT gamma)
  {
    const int i = get(vimap, v);

    // circulate over vertices around 'v' to compute w_ii and w_ij's
    NT w_ii = 0;
    int vi = 0;

    vertex_around_target_circulator v_j(halfedge(v, tmesh), tmesh), end = v_j;
    CGAL_For_all(v_j, end)
    {
      // Get j index
      const int j = get(vimap, *v_j);

      // NT w_ij = A_prev.get_coef(i, j) / compute_sig_ij(v, *v_j) / gamma;
      NT w_ij = A_prev.get_coef(i, j) / pow(compute_sig_ij(v, *v_j, NT(1)), gamma);

      // w_ii = - sum of w_ij's
      w_ii -= w_ij;

      // Set w_ij in matrix
      A.set_coef(i, j, w_ij, false);

      ++vi;
    }

    if(vi < 2)
      return ERROR_NON_TRIANGULAR_MESH;

    // Set w_ii in matrix
    A.set_coef(i,i, w_ii, true /*new*/);
    return OK;
  }

  NT compute_sig_ij(const vertex_descriptor v_i,
                        const vertex_descriptor v_j,
                        const NT gamma)
  {
    const NT sig = (pow(get(m_vertex_L2_map, v_i), gamma) + pow(get(m_vertex_L2_map, v_j), gamma)) / 2.;
    CGAL_assertion(sig > NT(0));
    return sig;
  }

  template <typename VertexUVMap>
  void store_new_best_uv_map(const Vertex_set& vertices,
                             VertexUVMap uvmap)
  {
    for(vertex_descriptor v : vertices)
      put(m_last_best_uv_map, v, get(uvmap, v));
  }

public:
  /// initializes `A`, `Bu`, and `Bv` after border parameterization.
  /// Fill the border vertices' lines in both linear systems:
  /// "u = constant" and "v = constant".
  ///
  /// \tparam VertexUVmap must be a model of `ReadWritePropertyMap` with
  ///         `boost::graph_traits<Triangle_mesh>::%vertex_descriptor` as key type and
  ///         %Point_2 (type deduced from `Triangle_mesh` using `Kernel_traits`)
  ///         as value type.
  /// \tparam VertexIndexMap must be a model of `ReadablePropertyMap` with
  ///         `boost::graph_traits<Triangle_mesh>::%vertex_descriptor` as key type and
  ///         a unique integer as value type.
  ///
  /// \param A the matrix in both linear system
  /// \param Bu the right hand side vector in the linear system of x coordinates
  /// \param Bv the right hand side vector in the linear system of y coordinates
  /// \param tmesh a triangulated surface
  /// \param bhd a halfedge descriptor on the boundary of `mesh`
  /// \param uvmap an instanciation of the class `VertexUVmap`
  /// \param vimap an instanciation of the class `VertexIndexMap`
  ///
  /// \pre Vertices must be indexed (`vimap` must be initialized).
  /// \pre `A`, `Bu`, and `Bv` must be allocated.
  /// \pre Border vertices must be parameterized.
  template <typename VertexUVmap, typename VertexIndexMap>
  void initialize_system_from_mesh_border(Matrix& A, Vector& Bu, Vector& Bv,
                                          const Triangle_mesh& tmesh,
                                          halfedge_descriptor bhd,
                                          VertexUVmap uvmap,
                                          VertexIndexMap vimap) const
  {
    for(halfedge_descriptor hd : halfedges_around_face(bhd, tmesh))
    {
      // Get vertex index in sparse linear system
      int index = get(vimap, target(hd, tmesh));

      // Write a diagonal coefficient of A
      A.set_coef(index, index, 1, true /*new*/);

      // get the halfedge uv
      // Write constant in Bu and Bv
      const Point_2& uv = get(uvmap, target(hd, tmesh));
      Bu[index] = uv.x();
      Bv[index] = uv.y();
    }
  }

  template <typename VertexUVmap,
            typename VertexIndexMap,
            typename VertexParameterizedMap>
  Error_code parameterize(Triangle_mesh& tmesh,
                          halfedge_descriptor bhd,
                          VertexUVmap uvmap,
                          VertexIndexMap vimap,
                          VertexParameterizedMap vpmap,
                          unsigned int& iterations,
                          double& error)
  {
    CGAL_precondition(is_valid_polygon_mesh(tmesh));
    CGAL_precondition(is_triangle_mesh(tmesh));
    CGAL_precondition(is_border(bhd, tmesh));

    Error_code status = OK;

    Vertex_set cc_vertices;
    std::vector<face_descriptor> cc_faces;
    cc_faces.reserve(num_faces(tmesh));

    internal::Containers_filler<Triangle_mesh, Vertex_set> fc(tmesh, cc_vertices, &cc_faces);
    Polygon_mesh_processing::connected_component(face(opposite(bhd, tmesh), tmesh), tmesh,
                                                 boost::make_function_output_iterator(fc));

    std::size_t nv = cc_vertices.size();
    if(nv == 0)
      return ERROR_EMPTY_MESH;

    // Compute (u,v) for border vertices and mark them as "parameterized"
    status = get_border_parameterizer().parameterize(tmesh, bhd, uvmap, vimap, vpmap);
    if (status != OK)
      return status;

    // Create two sparse linear systems "A*Xu = Bu" and "A*Xv = Bv" (one line/column per vertex)
    Matrix A(nv, nv);
    Matrix A_prev(nv, nv);
    Vector Xu(nv), Xv(nv), Bu(nv), Bv(nv);
    std::vector<double> err(iterations);

    // Initialize A, Xu, Xv, Bu and Bv after border parameterization
    // Fill the border vertices' lines in both linear systems:
    // "u = constant" and "v = constant"
    initialize_system_from_mesh_border(A, Bu, Bv, tmesh, bhd, uvmap, vimap);

    // Fill the matrix for the inner vertices v_i: compute A's coefficient
    // w_ij for each neighbor j; then w_ii = - sum of w_ijs
    std::unordered_set<vertex_descriptor> main_border;
    for(vertex_descriptor v : vertices_around_face(bhd, tmesh))
      main_border.insert(v); // @todo use marks?

    m_last_best_uv_map = get(Vertex_point2_tag(), tmesh);

    NT area_3D = initialize_faces_areas(cc_faces, tmesh);

    if(DEBUG_L0)
      std::cout << std::endl;

    unsigned int last_best_i = 0;
    NT gamma = 1; // @todo what value should that be
    bool is_changed = false;

    // iterate it with the new weights
    unsigned int i = 0;
    while(i < iterations)
    {
      if(DEBUG_L0)
        std::cout << "Iteration " << i << ", gamma = " << gamma << std::flush;

      // update weights for inner vertices
      for(vertex_descriptor v : cc_vertices)
      {
        // inner vertices only
        if(main_border.count(v) == 0)
        {
          // Compute the line i of matrix A for i inner vertex
          if(i == 0)
          {
#if 0
            status = setup_inner_vertex_relations(A, A_prev, Bu, Bv, tmesh, v, vimap);
#elif 1
            status = setup_inner_vertex_relations_cotangent(A, Bu, Bv, tmesh, v, vimap);
#else
            status = setup_inner_vertex_relations_MVC(A, Bu, Bv, tmesh, v, vimap);
#endif

            if(status != OK)
              return status;
          }
          else
          {
            status = setup_iter_inner_vertex_relations(A, A_prev, Bu, Bv, tmesh, v, vimap, gamma);
            if(status != OK)
              return status;
          }
        }
      }

      // solve linear equations
      // Solve "A*Xu = Bu". On success, solution is (1/Du) * Xu.
      // Solve "A*Xv = Bv". On success, solution is (1/Dv) * Xv.
      NT Du = 0, Dv = 0;
      if(!get_linear_algebra_traits().linear_solver(A, Bu, Xu, Du) ||
         !get_linear_algebra_traits().linear_solver(A, Bv, Xv, Dv))
      {
        if(DEBUG_L0)
          std::cout << " Linear solver failure #" << m_linear_solver_failures << std::endl;

        status = ERROR_CANNOT_SOLVE_LINEAR_SYSTEM;
      }
      else
      {
        m_linear_solver_failures = 0;
      }

      if(status != OK)
      {
        if(m_linear_solver_failures < 4)
        {
          // modify the weights and re-try the linear solver
          ++m_linear_solver_failures;
          gamma /= 2;
          continue;
        }
        else
        {
          status = OK;
          break;
        }
      }

      // WARNING: this package does not support homogeneous coordinates!
      CGAL_postcondition(Du == NT(1));
      CGAL_postcondition(Dv == NT(1));

//      for(std::size_t i=0; i<nv; ++i)
//        std::cout << "Sol[" << i << "] = " << Xu[i] << " " << Xv[i] << std::endl;

      // Copy A to A_prev, it is a computationally inefficient task but neccesary
      copy_sparse_matrix(A, A_prev, tmesh, cc_vertices, vimap);

      // Copy Xu and Xv coordinates into the (u,v) pair of each vertex
      for(vertex_descriptor v : cc_vertices)
      {
        // inner vertices only
        if(main_border.count(v) == 0)
        {
          int index = get(vimap,v);
          put(uvmap, v, Point_2(Xu[index], Xv[index]));
          put(vpmap, v, true);
        }
      }

      if(DEBUG_L0)
      {
        std::ofstream out("last_solve.off");
        out.precision(17);
        IO::output_uvmap_to_off(tmesh, bhd, uvmap, out);
        out.close();
      }

      compute_faces_L2(cc_faces, tmesh, uvmap);
      compute_vertices_L2(cc_vertices, tmesh);

      err[i] = compute_area_distortion(cc_faces, area_3D, tmesh, uvmap);

      if(DEBUG_L0)
        std::cout << " err " << err[i] << std::flush;

      if(err[i] <= err[last_best_i])
      {
        store_new_best_uv_map(cc_vertices, uvmap);

        last_best_i = i;
        is_changed = false;

        if(DEBUG_L0)
          std::cout << " *****" << std::flush;
      }
      else if(err[i] > 100) // @fixme is that reasonnable
      {
        break;
      }
      else
      {
        if(!is_changed)
        {
          gamma /= 2;
          is_changed = true;
        }
      }

      ++i;
      std::cout << std::endl;
    }

    // Check postconditions
    if(status != OK)
      return status;

    if(i == 0 && i != iterations)
    {
      // means that the computation terminated for the first iteration may be because system was unsolvable
      return ERROR_CANNOT_SOLVE_LINEAR_SYSTEM;
    }

    for(vertex_descriptor v : cc_vertices)
      put(uvmap, v, get(m_last_best_uv_map, v));

    iterations = last_best_i;
    error = err[last_best_i];

    return status;
  }

  template <typename VertexUVmap,
            typename VertexIndexMap,
            typename VertexParameterizedMap>
  Error_code parameterize(Triangle_mesh& tmesh,
                          halfedge_descriptor bhd,
                          VertexUVmap uvmap,
                          VertexIndexMap vimap,
                          VertexParameterizedMap vpmap,
                          unsigned int& iterations)
  {
    double unused_error;
    return parameterize(tmesh, bhd, uvmap, vimap, vpmap, iterations, unused_error);
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
  /// \param tmesh a triangulated surface
  /// \param bhd a halfedge descriptor on the boundary of `mesh`
  /// \param uvmap an instanciation of the class `VertexUVmap`
  /// \param vimap an instanciation of the class `VertexIndexMap`
  /// \param vpmap an instanciation of the class `VertexParameterizedMap`
  /// \param iterations an integer number of iterations to run the parameterization
  ///
  /// \pre `tmesh` must be a triangular mesh.
  /// \pre The mesh border must be mapped onto a convex polygon.
  /// \pre The vertices must be indexed (`vimap` must be initialized).
  ///
  template <typename VertexUVmap,
            typename VertexIndexMap,
            typename VertexParameterizedMap>
  Error_code parameterize(Triangle_mesh& tmesh,
                          halfedge_descriptor bhd,
                          VertexUVmap uvmap,
                          VertexIndexMap vimap,
                          VertexParameterizedMap vpmap,
                          // the '&' below is important, otherwise the function just calls itself
                          const unsigned int& iterations = 15)
  {
    unsigned int iter = iterations; // need a non-const ref
    return parameterize(tmesh, bhd, uvmap, vimap, vpmap, iter);
  }

  template <typename VertexUVmap>
  Error_code parameterize(Triangle_mesh& tmesh,
                          halfedge_descriptor bhd,
                          VertexUVmap uvmap,
                          const unsigned int iterations = 15)
  {
    Vertex_int_map vimap = get(Vertex_int_tag(), tmesh);
    internal::fill_index_map_of_cc(bhd, tmesh, vimap);

    Vertex_bool_map vpmap = get(Vertex_bool_tag(), tmesh);

    return parameterize(tmesh, bhd, uvmap, vimap, vpmap, iterations);
  }

  template <typename VertexUVmap>
  Error_code parameterize(Triangle_mesh& tmesh,
                          VertexUVmap uvmap,
                          const unsigned int iterations = 15)
  {
    const halfedge_descriptor bhd = Polygon_mesh_processing::longest_border(tmesh).first;

    return parameterize(tmesh, bhd, uvmap, iterations);
  }
};

} // namespace Surface_mesh_parameterization
} // namespace CGAL

#endif // CGAL_SURFACE_MESH_PARAMETERIZATION_ITERATIVE_AUTHALIC_PARAMETERIZER_3_H
