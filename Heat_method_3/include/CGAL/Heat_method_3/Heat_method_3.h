// Copyright (c) 2018  Carnegie Mellon University (USA), GeometryFactory (France)
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
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s) : Christina Vaz, Keenan Crane, Andreas Fabri

#ifndef CGAL_HEAT_METHOD_3_HEAT_METHOD_3_H
#define CGAL_HEAT_METHOD_3_HEAT_METHOD_3_H

#include <CGAL/license/Heat_method_3.h>
#include <CGAL/Heat_method_3/Intrinsic_Delaunay_triangulation_3.h>
#include <CGAL/Heat_method_3/internal/V2V.h>
#include <CGAL/disable_warnings.h>

#include <CGAL/property_map.h>
#include <CGAL/double.h>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/Dynamic_property_map.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/number_utils.h>
#ifdef CGAL_EIGEN3_ENABLED
#include <Eigen/Cholesky>
#include <CGAL/Eigen_matrix.h>
#include <CGAL/Eigen_vector.h>
#include <CGAL/Eigen_solver_traits.h>
#endif

#include <boost/foreach.hpp>

#include <vector>
#include <set>

#ifdef CGAL_TESTSUITE
struct Heat_method_3_private_tests;
#endif

namespace CGAL {
namespace Heat_method_3 {


/**
 * \ingroup PkgHeatMethod
 * 
 * Class `Heat_method_3` is an implementation of the Heat Method by Crane, et al, an algorithm that computes geodesic distance.
 * \tparam TriangleMesh a triangulated surface mesh, model of `FaceGraph` and `HalfedgeListGraph`
 * \tparam Traits a model of HeatMethodTraits_3
 * \tparam LA a model of `SparseLinearAlgebraWithFactorTraits_d`.

 * \tparam VertexPointMap a model of `ReadablePropertyMap` with
 *        `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key and
 *        `Traits::Point_3` as value type.
 *        The default is `typename boost::property_map<TriangleMesh, vertex_point_t>::%type`.
 *
 */
template <typename TriangleMesh,
          typename Traits,
          typename VertexDistanceMap,
#ifdef CGAL_EIGEN3_ENABLED
          typename LA = Eigen_solver_traits<Eigen::SimplicialLDLT<typename Eigen_sparse_matrix<double>::EigenType > >,
#else
          typename LA,
#endif
          typename VertexPointMap = typename boost::property_map< TriangleMesh, vertex_point_t>::const_type>
class Heat_method_3
{
#ifdef CGAL_TESTSUITE
  friend Heat_method_3_private_tests;
#endif
  /// Polygon_mesh typedefs
  typedef typename boost::graph_traits<TriangleMesh>               graph_traits;
  typedef typename graph_traits::vertex_descriptor            vertex_descriptor;
  typedef typename graph_traits::edge_descriptor                edge_descriptor;
  typedef typename graph_traits::halfedge_descriptor        halfedge_descriptor;
  typedef typename graph_traits::face_descriptor                face_descriptor;
  typedef typename std::set<vertex_descriptor>::iterator        vertex_iterator;
  /// Geometric typedefs
  typedef typename Traits::Point_3                                      Point_3;
  typedef typename Traits::FT                                                FT;
  typedef typename Traits::Vector_3                                    Vector_3;
  typedef typename Traits::Point_2                                      Point_2;

  // Property map typedefs
  typedef typename boost::property_traits<VertexPointMap>::reference VertexPointMap_reference;

  typedef typename LA::Matrix Matrix;
  typedef typename LA::Vector Vector;
  typedef typename LA::Index Index;

  // The Vertex_id_map is a property map where you can associate an index to a vertex_descriptor
  typedef typename boost::graph_traits<TriangleMesh>::vertices_size_type vertices_size_type;

  typedef CGAL::dynamic_vertex_property_t<Index> Vertex_property_tag;
  typedef typename boost::property_map<TriangleMesh, Vertex_property_tag >::const_type Vertex_id_map;
  Vertex_id_map vertex_id_map;

  typedef CGAL::dynamic_face_property_t<Index> Face_property_tag;
  typedef typename boost::property_map<TriangleMesh, Face_property_tag >::const_type Face_id_map;
  Face_id_map face_id_map;

public:

  /*!
    \brief Constructor
  */
  Heat_method_3(const TriangleMesh& tm, VertexDistanceMap vdm)
    : vertex_id_map(get(Vertex_property_tag(),tm)), face_id_map(get(Face_property_tag(),tm)), v2v(tm), tm(tm), vdm(vdm), vpm(get(vertex_point,tm))
  {
    build();
  }

  
  /*!
    \brief Constructor
  */
  Heat_method_3(const TriangleMesh& tm, VertexDistanceMap vdm, VertexPointMap vpm)
    : v2v(tm), tm(tm), vdm(vdm), vpm(vpm)
  {
    build();
  }

  
  /**
   * returns the triangle mesh the algorithm is running on.
   */
  const TriangleMesh& triangle_mesh() const{
    return tm;
  }


private:
  
  const VertexDistanceMap&
  vertex_distance_map() const
  {
    return vdm;
  }

  
  const Matrix&
  mass_matrix() const
  {
    return m_mass_matrix;
  }

  
  const Matrix&
  cotan_matrix() const
  {
    return m_cotan_matrix;
  }

  
  const VertexPointMap&
  vertex_point_map() const
  {
    return vpm;
  }

  
  const Vertex_id_map&
  get_vertex_id_map() const
  {
    return vertex_id_map;
  }

public:

  /**
   * adds `vd` to the source set, returning `false` if `vd` is already in the set.
   */
  template <typename VD>
  bool
  add_source(VD vd)
  {
    source_change_flag = true;
    return sources.insert(v2v(vd)).second;
  }

  
  /**
   * removes vd` from the source set, returning 'true' if `vd` was in the set.
   */
  template <typename VD>
  bool
  remove_source(VD vd)
  {
    source_change_flag = true;
    return (sources.erase(v2v(vd)) == 1);
  }


  /**
   * clears the current source set.
   */
  void
  clear_sources()
  {
    source_change_flag = true;
    sources.clear();
    return;
  }


  /**
   * get distance from the current source set to a vertex `vd`.
   */
  double
  distance(vertex_descriptor vd) const
  {
    return get(vdm,vd);
  }


  /**
   * returns an iterator to the first vertex in the source set.
   */
  vertex_iterator
  sources_begin()
  {
    return sources.begin();
  }

  
  /**
   * returns past-the-end iterator of the source set.
   */
  vertex_iterator
  sources_end()
  {
    return sources.end();
  }

private:
  
  double
  summation_of_edges() const
  {
    double edge_sum = 0;
    CGAL::Vertex_around_face_iterator<TriangleMesh> vbegin, vend, vmiddle;
    BOOST_FOREACH(face_descriptor f, faces(tm)) {
      halfedge_descriptor hd = halfedge(f,tm);

      vertex_descriptor neighbor_one = target(hd,tm);
      halfedge_descriptor hd2 = next(hd,tm);
      vertex_descriptor neighbor_two = target(hd2,tm);
      halfedge_descriptor hd3 = next(hd2,tm);
      vertex_descriptor current = target(hd3,tm);

      VertexPointMap_reference pi = get(vpm,current);
      VertexPointMap_reference pj = get(vpm, neighbor_one);
      VertexPointMap_reference pk = get(vpm, neighbor_two);
      edge_sum += (CGAL::is_border(opposite(hd,tm),tm)?1.0:0.5) * std::sqrt(CGAL::squared_distance(pi,pj));
      edge_sum += (CGAL::is_border(opposite(hd2,tm),tm)?1.0:0.5) * std::sqrt(CGAL::squared_distance(pj,pk)) ;
      edge_sum += (CGAL::is_border(opposite(hd3,tm),tm)?1.0:0.5) * std::sqrt(CGAL::squared_distance(pk,pi)) ;
    }
    return edge_sum;
  }

  
  double
  time_step() const
  {
    return m_time_step;
  }

  
  void
  update_kronecker_delta()
  {
    //currently just working with a single vertex in source set, add the first one for now
    Index i;
    Matrix K(num_vertices(tm), 1);
    if(sources.empty()) {
      i = 0;
      K.set_coef(i,0, 1, true);
      source_index.insert(0);
    } else {
      vertex_descriptor current = *(sources.begin());
      i = get(vertex_id_map, current);
      vertex_iterator vd =sources.begin();
      while(vd!=sources.end()){
        i = get(vertex_id_map, *(vd));
        K.set_coef(i,0, 1, true);
        vd = ++vd;
      }
    }
    kronecker.swap(K);
  }


  const Matrix&
  kronecker_delta() const
  {
    return kronecker;
  }

  
  void solve_cotan_laplace()
  {
    Matrix A, A0;
    A0 = m_time_step * m_cotan_matrix;
    A = m_mass_matrix + A0;

    double d=0;
    
    if(! la.factor(A,d)) {
      // decomposition failed
      CGAL_error_msg("Eigen Decomposition in cotan failed");
    }
         
    if(! la.linear_solver(kronecker, solved_u)) {
      // solving failed
      CGAL_error_msg("Eigen Solving in cotan failed");
    }
  }

  
  void
  compute_unit_gradient()
  {
    if(X.empty()){
      X.resize(num_faces(tm));
    }
    CGAL::Vertex_around_face_iterator<TriangleMesh> vbegin, vend, vmiddle;
    BOOST_FOREACH(face_descriptor f, faces(tm)) {
      boost::tie(vbegin, vend) = vertices_around_face(halfedge(f,tm),tm);
      vertex_descriptor current = *(vbegin);
      vertex_descriptor neighbor_one = *(++vbegin);
      vertex_descriptor neighbor_two = *(++vbegin);
      VertexPointMap_reference p_i = get(vpm,current);
      VertexPointMap_reference p_j = get(vpm, neighbor_one);
      VertexPointMap_reference p_k = get(vpm, neighbor_two);
      Index face_i = get(face_id_map, f);

      Index i = get(vertex_id_map, current);
      Index j = get(vertex_id_map, neighbor_one);
      Index k = get(vertex_id_map, neighbor_two);

      //get area of face_i
      //get outward unit normal
      //cross that with eij, ejk, eki
      //so (Ncross eij) *uk and so on
      //sum all of those then multiply by 1./(2a)
      Vector_3 cross = CGAL::cross_product((p_j-p_i), (p_k-p_i));
      double N_cross = (CGAL::sqrt(cross*cross));
      Vector_3 unit_cross = cross/N_cross;
      double area_face = N_cross * (1./2);
      double u_i = CGAL::abs(solved_u(i));
      double u_j = CGAL::abs(solved_u(j));
      double u_k = CGAL::abs(solved_u(k));
      double r_Mag = 1./(std::max)((std::max)(u_i, u_j),u_k);
      /* normalize heat values so that they have roughly unit magnitude */
      if(!std::isinf(r_Mag)) {
        u_i = u_i * r_Mag;
        u_j = u_j * r_Mag;
        u_k = u_k * r_Mag;
      }
      Vector_3 edge_sums = u_k * CGAL::cross_product(unit_cross,(p_j-p_i));
      edge_sums = edge_sums + u_i * (CGAL::cross_product(unit_cross, (p_k-p_j)));
      edge_sums = edge_sums + u_j * CGAL::cross_product(unit_cross, (p_i-p_k));
      edge_sums = edge_sums * (1./area_face);
      double e_magnitude = CGAL::sqrt(edge_sums*edge_sums);
      X[face_i] = edge_sums*(1./e_magnitude);
    }
  }

  
  double
  dot_eigen_vector(const Vector_3& a, const Vector_3& b) const
  {
    return (a.x()*b.x() + a.y()*b.y() + a.z()*b.z());
  }


  void
  compute_divergence()
  {
    Matrix indexD(dimension,1);
    CGAL::Vertex_around_face_iterator<TriangleMesh> vbegin, vend, vmiddle;
    BOOST_FOREACH(face_descriptor f, faces(tm)) {
      boost::tie(vbegin, vend) = vertices_around_face(halfedge(f,tm),tm);
      vertex_descriptor current = *(vbegin);
      vertex_descriptor neighbor_one = *(++vbegin);
      vertex_descriptor neighbor_two = *(++vbegin);
      Index i = get(vertex_id_map, current);
      Index j = get(vertex_id_map, neighbor_one);
      Index k = get(vertex_id_map, neighbor_two);
      VertexPointMap_reference p_i = get(vpm,current);
      VertexPointMap_reference p_j = get(vpm, neighbor_one);
      VertexPointMap_reference p_k = get(vpm, neighbor_two);
      Index face_i = get(face_id_map, f);

      Vector_3 cross = CGAL::cross_product((p_j-p_i), (p_k-p_i));
      double norm_cross = (CGAL::sqrt(cross*cross));
      double dot = (p_j-p_i)*(p_k-p_i);
      double cotan_i = dot/norm_cross;

      cross = CGAL::cross_product((p_i-p_j), (p_k-p_j));
      dot = to_double((p_i-p_j)*(p_k-p_j));
      double cotan_j = dot/norm_cross;

      cross = CGAL::cross_product((p_i-p_k), (p_j-p_k));
      dot = to_double((p_i-p_k)*(p_j-p_k));
      double cotan_k = dot/norm_cross;

      const Vector_3& a  = X[face_i];
      double i_entry = cotan_k*(dot_eigen_vector(a,(p_j-p_i))) + cotan_j*(dot_eigen_vector(a,(p_k-p_i)));
      double j_entry = cotan_i*(dot_eigen_vector(a,(p_k-p_j))) + cotan_k*(dot_eigen_vector(a,(p_i-p_j)));
      double k_entry = cotan_j*(dot_eigen_vector(a,(p_i-p_k))) + cotan_i*(dot_eigen_vector(a,(p_j-p_k)));

      indexD.add_coef(i, 0, (1./2)*i_entry);
      indexD.add_coef(j, 0, (1./2)*j_entry);
      indexD.add_coef(k, 0, (1./2)*k_entry);
    }
    indexD.swap(m_index_divergence);
  }

    
  void
  value_at_source_set(const Vector& phi)
  {
    Vector source_set_val(dimension);
    if(sources.empty()) {
      for(int k = 0; k<dimension; k++) {
        source_set_val(k,0) = phi.coeff(0,0);
      }
      source_set_val = phi-source_set_val;
    } else {
      for(int i = 0; i<dimension; i++) {
        double min_val = (std::numeric_limits<double>::max)();
        vertex_iterator current;
        Index current_Index;
        current = sources.begin();
        //go through the distances to the sources and leave the minimum distance;
        for(int j = 0; j<sources.size(); j++) {
          current_Index = get(vertex_id_map, *current);
          double new_d = CGAL::abs(-phi.coeff(current_Index,0)+phi.coeff(i,0));
          if(phi.coeff(current_Index,0)==phi.coeff(i,0)) {
            min_val = 0.;
          }
          if(new_d < min_val) {
            min_val = new_d;
          }  
          current = ++current;
        }
        source_set_val(i,0) = min_val;
      }
    }
    solved_phi.swap(source_set_val);
  }


  void
  solve_phi()
  {
    Vector phi;
    double d = 1;
    if(! la.factor(m_cotan_matrix, d)) {
      // decomposition failed
      CGAL_error_msg("Eigen Decomposition in phi failed");
    }
      
    if(! la.linear_solver(m_index_divergence, phi)) {
      // solving failed
      CGAL_error_msg("Eigen Solving in phi failed");
    }  
    value_at_source_set(phi);
  }

  
  // this function returns a (number of vertices)x1 vector where
  // the ith index has the distance from the first vertex to the ith vertex
  const Vector&
  distances() const
  {
    return solved_phi;
  }

public:

  /**
   *  Updates the distance property map after changes in the source set. 
   **/
  void update()
  {
    double d=0;
    if(source_change_flag) {
      //don't need to recompute Mass matrix, cotan matrix or timestep reflect that in this function
      update_kronecker_delta();
      solve_cotan_laplace();
      compute_unit_gradient();
      compute_divergence();
      solve_phi();
      source_change_flag = false;
      std::cout<<"sources changed, recompute\n";
    }
    BOOST_FOREACH(vertex_descriptor vd, vertices(tm)){
      Index i_d = get(vertex_id_map, vd);
      d = solved_phi(i_d,0);
      put(vdm,vd,d);
    }
  }

private:

  void
  build()
  {
    source_change_flag = false;

    CGAL_precondition(is_triangle_mesh(tm));
    Index i = 0;
    BOOST_FOREACH(vertex_descriptor vd, vertices(tm)){
      put(vertex_id_map, vd, i++);
    }
    Index face_i = 0;
    BOOST_FOREACH(face_descriptor fd, faces(tm)){
      put(face_id_map, fd, face_i++);
    }
    dimension = static_cast<int>(num_vertices(tm));

    m_mass_matrix.eigen_object().resize(dimension, dimension);
    m_cotan_matrix.eigen_object().resize(dimension, dimension);
    CGAL::Vertex_around_face_iterator<TriangleMesh> vbegin, vend, vmiddle;

    BOOST_FOREACH(face_descriptor f, faces(tm)) {
      boost::tie(vbegin, vend) = vertices_around_face(halfedge(f,tm),tm);
      vertex_descriptor current = *(vbegin);
      vertex_descriptor neighbor_one = *(++vbegin);
      vertex_descriptor neighbor_two = *(++vbegin);
      Index i = get(vertex_id_map, current);
      Index j = get(vertex_id_map, neighbor_one);
      Index k = get(vertex_id_map, neighbor_two);
      Point_3 pi, pj, pk;

      VertexPointMap_reference p_i = get(vpm,current);
      VertexPointMap_reference p_j = get(vpm, neighbor_one);
      VertexPointMap_reference p_k = get(vpm, neighbor_two);
      pi = p_i;
      pj = p_j;
      pk = p_k;

      Vector_3 cross = CGAL::cross_product((pj-pi), (pk-pi));
      double dot = (pj-pi)*(pk-pi);

      double norm_cross = (CGAL::sqrt(cross*cross));

      double cotan_i = dot/norm_cross;
      m_cotan_matrix.add_coef(j,k ,-(1./2)*cotan_i);
      m_cotan_matrix.add_coef(k,j,-(1./2)* cotan_i);
      m_cotan_matrix.add_coef(j,j,(1./2)*cotan_i);
      m_cotan_matrix.add_coef(k,k,(1./2)* cotan_i);

      cross = CGAL::cross_product((pi-pj), (pk-pj));
      dot = to_double((pi-pj)*(pk-pj));
      double cotan_j = dot/norm_cross;
      m_cotan_matrix.add_coef(i,k ,-(1./2)*cotan_j);
      m_cotan_matrix.add_coef(k,i,-(1./2)* cotan_j);
      m_cotan_matrix.add_coef(i,i,(1./2)* cotan_j);
      m_cotan_matrix.add_coef(k,k,(1./2)* cotan_j);

      cross = CGAL::cross_product((pi-pk), (pj-pk));
      dot = to_double((pi-pk)*(pj-pk));
      double cotan_k = dot/norm_cross;
      m_cotan_matrix.add_coef(i,j,-(1./2)*cotan_k);
      m_cotan_matrix.add_coef(j,i,-(1./2)* cotan_k);
      m_cotan_matrix.add_coef(i,i,(1./2)* cotan_k);
      m_cotan_matrix.add_coef(j,j,(1./2)* cotan_k);

      //double area_face = CGAL::Polygon_mesh_processing::face_area(f,tm);
      //cross is 2*area
      m_mass_matrix.add_coef(i,i, (1./6.)*norm_cross);
      m_mass_matrix.add_coef(j,j, (1./6.)*norm_cross);
      m_mass_matrix.add_coef(k,k, (1./6.)*norm_cross);
      m_cotan_matrix.add_coef(i,i, 1e-8);
      m_cotan_matrix.add_coef(j,j, 1e-8);
      m_cotan_matrix.add_coef(k,k, 1e-8);
    }

    m_time_step = 1./(num_edges(tm));
    m_time_step = m_time_step*summation_of_edges();
    m_time_step = m_time_step*m_time_step;
    update_kronecker_delta();
    solve_cotan_laplace();
    compute_unit_gradient();
    compute_divergence();
    solve_phi();
  }

  LA la;
  int dimension;
  V2V<TriangleMesh> v2v;
  const TriangleMesh& tm;
  VertexDistanceMap vdm;
  VertexPointMap vpm;
  std::set<vertex_descriptor> sources;
  double m_time_step;
  Matrix kronecker;
  Matrix m_mass_matrix, m_cotan_matrix;
  Vector solved_u;
  std::vector<Vector_3> X;
  Matrix m_index_divergence;
  Vector solved_phi;
  std::set<Index> source_index;
  bool source_change_flag;
};


  /*! \addtogroup PkgHeatMethod
   *
   * @{
   */

/// \sa CGAL::Heat_method_3::Heat_method_3
/// computes for each vertex  of the triangle mesh `tm` the geodesic distance to a given source vertex. 
template <typename TriangleMesh, typename VertexDistanceMap>
void
geodesic_distances_3(const TriangleMesh& tm,
                     VertexDistanceMap vdm,
                     typename boost::graph_traits<TriangleMesh>::vertex_descriptor source)
{
  typedef typename boost::property_map<TriangleMesh, vertex_point_t>::type PPM;
  typedef typename boost::property_traits<PPM>::value_type Point_3;
  typedef typename CGAL::Kernel_traits<Point_3>::Kernel Kernel;
  typedef CGAL::Heat_method_3::Heat_method_3<TriangleMesh,Kernel,VertexDistanceMap> Heat_method;
  
  Heat_method hm(tm,vdm);
  hm.add_source(source);
  hm.update();
}


/// \sa CGAL::Heat_method_3::Heat_method_3
/// computes for each vertex of the triangle mesh `tm` the geodesic distance to a given source vertex. This version computes better results when `tm` has triangles with small angles. 
template <typename TriangleMesh, typename VertexDistanceMap>
void
geodesic_distances_with_intrinsic_Delaunay_triangulation_3(const TriangleMesh& tm,
                                                           VertexDistanceMap vdm,
                                                           typename boost::graph_traits<TriangleMesh>::vertex_descriptor source)
{
  typedef typename boost::property_map<TriangleMesh, vertex_point_t>::type PPM;
  typedef typename boost::property_traits<PPM>::value_type Point_3;
  typedef typename CGAL::Kernel_traits<Point_3>::Kernel Kernel;
  typedef CGAL::Heat_method_3::Intrinsic_Delaunay_triangulation_3<TriangleMesh,Kernel, VertexDistanceMap> Idt;
  typedef CGAL::Heat_method_3::Heat_method_3<Idt,Kernel,typename Idt::Vertex_distance_map> Heat_method;
  
  Idt idt(tm, vdm);
  Heat_method hm(idt,idt.vertex_distance_map());
  hm.add_source(source);
  hm.update();
}
  

 /*! @} */

} // namespace Heat_method_3
} // namespace CGAL

#include <CGAL/enable_warnings.h>
#endif // CGAL_HEAT_METHOD_3_HEAT_METHOD_3_H
