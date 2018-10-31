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

#ifndef CGAL_HEAT_METHOD_3_SURFACE_MESH_GEODESIC_DISTANCES_3_H
#define CGAL_HEAT_METHOD_3_SURFACE_MESH_GEODESIC_DISTANCES_3_H

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
 * A tag class used to specify that the heat method is applied to the input mesh.
 */
struct Direct
{};

/**
 * A tag class used to specify that the heat method applies the intrinsic Delaunay triangulation to the input mesh.
 */
struct Intrinsic_Delaunay
{};

namespace internal {
template <typename TriangleMesh,
          typename Traits,
          typename LA,
          typename VertexPointMap>
class Surface_mesh_geodesic_distances_3
{
protected:
#ifdef CGAL_TESTSUITE
  friend Surface_mesh_geodesic_distances_3_private_tests;
#endif
  /// Polygon_mesh typedefs
  typedef typename boost::graph_traits<TriangleMesh>               graph_traits;
  typedef typename graph_traits::vertex_descriptor            vertex_descriptor;
  typedef typename graph_traits::edge_descriptor                edge_descriptor;
  typedef typename graph_traits::halfedge_descriptor        halfedge_descriptor;
  typedef typename graph_traits::face_descriptor                face_descriptor;
  typedef typename std::set<vertex_descriptor>                     Vertex_range;
  /// Geometric typedefs
  typedef typename Traits::Point_3                                      Point_3;
  typedef typename Traits::FT                                                FT;
  typedef typename Traits::Vector_3                                    Vector_3;

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
  Surface_mesh_geodesic_distances_3(const TriangleMesh& tm)
    : vertex_id_map(get(Vertex_property_tag(),tm)), face_id_map(get(Face_property_tag(),tm)), v2v(tm), tm(tm), vpm(get(vertex_point,tm))
  {
    build();
  }


  /*!
    \brief Constructor
  */
  Surface_mesh_geodesic_distances_3(const TriangleMesh& tm, VertexPointMap vpm)
    : v2v(tm), tm(tm), vpm(vpm)
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
    return m_sources.insert(v2v(vd)).second;
  }

  /**
   * adds vertices in the `vrange` to the source set'.
   */
  template <typename VertexRange>
  void
  add_sources(VertexRange vrange)
  {
    typedef typename std::iterator_traits<typename VertexRange::const_iterator>::value_type value_type;
    source_change_flag = true;
    BOOST_FOREACH(value_type vd, vrange){
      m_sources.insert(v2v(vd));
    }
  }


  /**
   * removes `vd` from the source set, returning `true` if `vd` was in the set.
   */
  template <typename VD>
  bool
  remove_source(VD vd)
  {
    source_change_flag = true;
    return (m_sources.erase(v2v(vd)) == 1);
  }


  /**
   * clears the current source set.
   */
  void
  clear_sources()
  {
    source_change_flag = true;
    m_sources.clear();
    return;
  }

  /**
   * returns the set of  source vertices.
   */

  const Vertex_range&
  sources() const
  {
    return m_sources;
  }

private:

  double
  summation_of_edges() const
  {
    typename Traits::Compute_squared_distance_3 squared_distance = Traits().compute_squared_distance_3_object();
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
      edge_sum += (CGAL::is_border(opposite(hd,tm),tm)?1.0:0.5) * std::sqrt(squared_distance(pi,pj));
      edge_sum += (CGAL::is_border(opposite(hd2,tm),tm)?1.0:0.5) * std::sqrt(squared_distance(pj,pk)) ;
      edge_sum += (CGAL::is_border(opposite(hd3,tm),tm)?1.0:0.5) * std::sqrt(squared_distance(pk,pi)) ;
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
    Matrix K(static_cast<int>(num_vertices(tm)), 1);
    if(sources().empty()) {
      i = 0;
      K.set_coef(i,0, 1, true);
      source_index.insert(0);
    } else {
      vertex_descriptor current = *(sources().begin());
      i = get(vertex_id_map, current);
      BOOST_FOREACH(vertex_descriptor vd, sources()){
        i = get(vertex_id_map, vd);
        K.set_coef(i,0, 1, true);
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
    typename Traits::Construct_vector_3 construct_vector = Traits().construct_vector_3_object();
    typename Traits::Construct_sum_of_vectors_3 sum = Traits().construct_sum_of_vectors_3_object();
    typename Traits::Compute_scalar_product_3 scalar_product = Traits().compute_scalar_product_3_object();
    typename Traits::Construct_cross_product_vector_3 cross_product = Traits().construct_cross_product_vector_3_object();
    typename Traits::Construct_scaled_vector_3 scale = Traits().construct_scaled_vector_3_object();
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

      Vector_3 v_ij = construct_vector(p_i,p_j);
      Vector_3 v_ik = construct_vector(p_i,p_k);

      Vector_3 cross = cross_product(v_ij, v_ik);
      double N_cross = (CGAL::sqrt(scalar_product(cross,cross)));
      Vector_3 unit_cross = scale(cross, 1./N_cross);
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
      Vector_3 edge_sums = scale(cross_product(unit_cross,v_ij), u_k);
      edge_sums = sum(edge_sums, scale(cross_product(unit_cross, construct_vector(p_j,p_k)), u_i));
      edge_sums = sum(edge_sums, scale(cross_product(unit_cross, construct_vector(p_k,p_i)), u_j));
      edge_sums = scale(edge_sums, (1./area_face));
      double e_magnitude = CGAL::sqrt(scalar_product(edge_sums,edge_sums));
      X[face_i] = scale(edge_sums,(1./e_magnitude));
    }
  }


  void
  compute_divergence()
  {
    typename Traits::Construct_cross_product_vector_3 cross_product = Traits().construct_cross_product_vector_3_object();
    typename Traits::Compute_scalar_product_3 scalar_product = Traits().compute_scalar_product_3_object();
    typename Traits::Construct_vector_3 construct_vector = Traits().construct_vector_3_object();
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

      Vector_3 v_ij = construct_vector(p_i,p_j);
      Vector_3 v_ik = construct_vector(p_i,p_k);
      Vector_3 cross = cross_product(v_ij, v_ik);
      double norm_cross = CGAL::sqrt(scalar_product(cross,cross));
      double dot = scalar_product(v_ij, v_ik);
      double cotan_i = dot/norm_cross;

      Vector_3 v_ji = construct_vector(p_j, p_i);
      Vector_3 v_jk = construct_vector(p_j, p_k);

      cross = cross_product(v_ji, v_jk);
      dot = to_double(scalar_product(v_ji, v_jk));
      double cotan_j = dot/norm_cross;

      Vector_3 v_ki = construct_vector(p_k,p_i);
      Vector_3 v_kj = construct_vector(p_k,p_j);

      cross = cross_product(v_ki, v_kj);
      dot = to_double(scalar_product(v_ki,v_kj));
      double cotan_k = dot/norm_cross;

      const Vector_3& a  = X[face_i];
      double i_entry = (scalar_product(a,v_ij) * cotan_k) + (scalar_product(a,v_ik) *  cotan_j);
      double j_entry = (scalar_product(a,v_jk) * cotan_i) + (scalar_product(a,v_ji) * cotan_k);
      double k_entry = (scalar_product(a,v_ki) * cotan_j) + (scalar_product(a,v_kj) * cotan_i);

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
    if(sources().empty()) {
      for(int k = 0; k<dimension; k++) {
        source_set_val(k,0) = phi.coeff(0,0);
      }
      source_set_val = phi-source_set_val;
    } else {
      for(int i = 0; i<dimension; i++) {
        double min_val = (std::numeric_limits<double>::max)();
        Index vd_index;
        //go through the distances to the sources and leave the minimum distance;
        BOOST_FOREACH(vertex_descriptor vd, sources()){
          vd_index = get(vertex_id_map, vd);
          double new_d = CGAL::abs(-phi.coeff(vd_index,0)+phi.coeff(i,0));
          if(phi.coeff(vd_index,0)==phi.coeff(i,0)) {
            min_val = 0.;
          }
          if(new_d < min_val) {
            min_val = new_d;
          }
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
  template<class VertexDistanceMap>
  void estimate_geodesic_distances(VertexDistanceMap vdm)
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
    typename Traits::Construct_cross_product_vector_3 cross_product = Traits().construct_cross_product_vector_3_object();
    typename Traits::Compute_scalar_product_3 scalar_product = Traits().compute_scalar_product_3_object();
    typename Traits::Construct_vector_3 construct_vector = Traits().construct_vector_3_object();

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
    {
      Matrix M(dimension,dimension);
      m_mass_matrix.swap(M);
    }
    {
      Matrix M(dimension,dimension);
      m_cotan_matrix.swap(M);
    }
    //m_mass_matrix.eigen_object().resize(dimension, dimension);
    //m_cotan_matrix.eigen_object().resize(dimension, dimension);
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

      Vector_3 v_ij = construct_vector(p_i,p_j);
      Vector_3 v_ik = construct_vector(p_i,p_k);

      Vector_3 cross = cross_product(v_ij, v_ik);
      double dot = scalar_product(v_ij,v_ik);

      double norm_cross = (CGAL::sqrt(scalar_product(cross,cross)));

      double cotan_i = dot/norm_cross;
      m_cotan_matrix.add_coef(j,k ,-(1./2)*cotan_i);
      m_cotan_matrix.add_coef(k,j,-(1./2)* cotan_i);
      m_cotan_matrix.add_coef(j,j,(1./2)*cotan_i);
      m_cotan_matrix.add_coef(k,k,(1./2)* cotan_i);

      Vector_3 v_ji = construct_vector(p_j,p_i);
      Vector_3 v_jk = construct_vector(p_j,p_k);

      cross = cross_product(v_ji, v_jk);
      dot = to_double(scalar_product(v_ji, v_jk));
      double cotan_j = dot/norm_cross;
      m_cotan_matrix.add_coef(i,k ,-(1./2)*cotan_j);
      m_cotan_matrix.add_coef(k,i,-(1./2)* cotan_j);
      m_cotan_matrix.add_coef(i,i,(1./2)* cotan_j);
      m_cotan_matrix.add_coef(k,k,(1./2)* cotan_j);

      Vector_3 v_ki = construct_vector(p_k,p_i);
      Vector_3 v_kj = construct_vector(p_k,p_j);
      cross = cross_product(v_ki, v_kj);
      dot = to_double(scalar_product(v_ki,v_kj));
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
    m_cotan_matrix.assemble_matrix();
    m_mass_matrix.assemble_matrix();

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
  VertexPointMap vpm;
  std::set<vertex_descriptor> m_sources;
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

template <typename TriangleMesh,
          typename Traits,
          typename Mode,
          typename LA,
          typename VertexPointMap>
struct Base_helper
  : public Surface_mesh_geodesic_distances_3<TriangleMesh, Traits, LA, VertexPointMap>
{
  typedef Surface_mesh_geodesic_distances_3<TriangleMesh, Traits, LA, VertexPointMap> type;

  Base_helper(const TriangleMesh& tm, VertexPointMap vpm)
    : type(tm, vpm)
  {}

  Base_helper(const TriangleMesh& tm)
    : type(tm)
  {}

  type& base()
  {
    return static_cast<type&>(*this);
  }

  const type& base() const
  {
    return static_cast<const type&>(*this);
  }

  template <class VertexDistanceMap>
  void estimate_geodesic_distances(VertexDistanceMap vdm)
  {
    base().estimate_geodesic_distances(vdm);
  }
};

template<class TriangleMesh,
         class Traits,
         class VertexPointMap>
struct Idt_storage
{
  Intrinsic_Delaunay_triangulation_3<TriangleMesh, Traits> m_idt;

  Idt_storage(const TriangleMesh& tm, VertexPointMap vpm)
    : m_idt(tm, vpm)
  {}

  Idt_storage(const TriangleMesh& tm)
    : m_idt(tm)
  {}
};

template <typename TriangleMesh,
          typename Traits,
          typename LA,
          typename VertexPointMap>
struct Base_helper<TriangleMesh, Traits, Intrinsic_Delaunay, LA, VertexPointMap>
  : public Idt_storage<TriangleMesh, Traits, VertexPointMap>
  , public Surface_mesh_geodesic_distances_3<Intrinsic_Delaunay_triangulation_3<TriangleMesh, Traits>,
                                             Traits,
                                             LA,
                                             typename Intrinsic_Delaunay_triangulation_3<TriangleMesh, Traits>::Vertex_point_map>
{
  typedef CGAL::Heat_method_3::Intrinsic_Delaunay_triangulation_3<TriangleMesh, Traits> Idt;
  typedef Idt_storage<TriangleMesh, Traits, VertexPointMap> Idt_wrapper;

  typedef Surface_mesh_geodesic_distances_3<Idt, Traits, LA, typename Idt::Vertex_point_map> type;

  Base_helper(const TriangleMesh& tm, VertexPointMap vpm)
    : Idt_wrapper(tm, vpm)
    , type(this->m_idt)
  {}

  Base_helper(const TriangleMesh& tm)
    :  Idt_wrapper(tm)
    , type(this->m_idt)
  {}

  type& base()
  {
    return static_cast<type&>(*this);
  }

  const type& base() const
  {
    return static_cast<const type&>(*this);
  }

  template <class VertexDistanceMap>
  void estimate_geodesic_distances(VertexDistanceMap vdm)
  {
    base().estimate_geodesic_distances(this->m_idt.vertex_distance_map(vdm));
  }
};

} // namespace internal


/**
 * \ingroup PkgHeatMethod
 *
 * Class `Surface_mesh_geodesic_distances_3` computes estimated geodesic distances for a set of source vertices where sources can be added and removed.
 * The class performs a preprocessing step that does only depend on the mesh, so that the distance computation takes less
 * time after changes of the set of sources.
 *
 * \tparam TriangleMesh a triangulated surface mesh, model of `FaceListGraph` and `HalfedgeListGraph`
 * \tparam Mode must be `Intrinsic_Delaunay` to indicate that an intrinsic Delaunay triangulation is internally constructed
 *                              or `Direct`to indicate that the input mesh should be used as is.
 *                              If `Intrinsic_Delaunay` the type `TriangleMesh` must have an internal property for `vertex_point`
 *                              and its value type must be the same as the value type of `VertexPointMap`.
 * \tparam VertexPointMap a model of `ReadablePropertyMap` with
 *        `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key and
 *        `Traits::Point_3` as value type.
 *        The default is `typename boost::property_map<TriangleMesh, vertex_point_t>::%type`.
 * \tparam LA a model of `SparseLinearAlgebraWithFactorTraits_d`.
 * \tparam Traits a model of `HeatMethodTraits_3`
 *
 */
template <typename TriangleMesh,
          typename Mode = Direct,
          typename VertexPointMap = typename boost::property_map< TriangleMesh, vertex_point_t>::const_type,
#ifdef CGAL_EIGEN3_ENABLED
          typename LA = Eigen_solver_traits<Eigen::SimplicialLDLT<typename Eigen_sparse_matrix<double>::EigenType > >,
#else
          typename LA = Default,
#endif
          typename Traits = typename Kernel_traits< typename boost::property_traits<VertexPointMap>::value_type>::Kernel >
class Surface_mesh_geodesic_distances_3
#ifndef DOXYGEN_RUNNING
  : public internal::Base_helper<TriangleMesh, Traits, Mode, LA, VertexPointMap>
#endif
{
  typedef internal::Base_helper<TriangleMesh, Traits, Mode, LA, VertexPointMap> Base_helper;

  const typename Base_helper::type& base() const
  {
    return Base_helper::base();
  }

  typename Base_helper::type& base()
  {
    return Base_helper::base();
  }

public:

  /// Vertex descriptor type
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;
  #ifndef DOXYGEN_RUNNING
  /// Source vertex iterator type
  typedef typename Base_helper::type::Vertex_range Vertex_range;
  #else
  /// a model of `ConstRange` with an iterator that is model of `ForwardIterator` with `vertex_descriptor` as value type.
  typedef unspecified_type Vertex_range;
  #endif

  /*!
    \brief Constructor
  */
  Surface_mesh_geodesic_distances_3(const TriangleMesh& tm)
    : Base_helper(tm)
  {}

  /*!
    \brief Constructor
  */
  Surface_mesh_geodesic_distances_3(const TriangleMesh& tm, VertexPointMap vpm)
    : Base_helper(tm, vpm)
  {}

  /**
   * returns the triangle mesh the algorithm is running on.
   */
  const TriangleMesh& triangle_mesh() const
  {
    return base().triangle_mesh();
  }

  /**
   * adds `vd` to the source set, returning `false` if `vd` is already in the set.
   */
  bool
  add_source(vertex_descriptor vd)
  {
    return base().add_source(vd);
  }

  /**
   * adds the range of vertices to the source set.
   */
  template <typename VertexRange>
  void
  add_sources(const Vertex_range& vrange)
  {
    base().add_sources(vrange);
  }

  /**
   * removes `vd` from the source set, returning `true` if `vd` was in the set.
   */
  bool
  remove_source(vertex_descriptor vd)
  {
    return base().remove_source(vd);
  }

  /**
   * clears the current source set.
   */
  void
  clear_sources()
  {
    base().clear_sources();
  }

  /**
   * get estimated distance from the current source set to a vertex `vd`.
   */
  double
  estimate_geodesic_distance(vertex_descriptor vd) const
  {
    return base().estimate_geodesic_distance(vd);
  }

  /**
   * returns the source set.
   */
  const Vertex_range&
  sources() const
  {
    return base().sources();
  }


  /**
   * fills the distance property map with the estimated geodesic distance of each vertex to the closest source vertex.
   * \tparam VertexDistanceMap a property map model of `WritablePropertyMap`
   * with `vertex_descriptor` as key type and `double` as value type.
   * \param vdm the vertex distance map to be filled
   **/
  template <class VertexDistanceMap>
  void estimate_geodesic_distances(VertexDistanceMap vdm)
  {
    Base_helper::estimate_geodesic_distances(vdm);
  }
};


/// \ingroup PkgHeatMethod
/// computes for each vertex of the triangle mesh `tm` the estimated geodesic distance to a given source vertex.
/// \tparam TriangleMesh a triangulated surface mesh, model of `FaceListGraph` and `HalfedgeListGraph`.
///         It must have an internal vertex point property map with the value type being a 3D point from a cgal Kernel model
/// \tparam VertexDistanceMap a property map model of `WritablePropertyMap`
///         with `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key type and `double` as value type.
/// \tparam Mode either the tag `Direct` or `Intrinsic_Delaunay`, which determines if the geodesic distance
///              is computed directly on the mesh, or if the intrinsic Delaunay triangulation is applied first.
///              The default is `Direct`.
///
/// \sa CGAL::Heat_method_3::Surface_mesh_geodesic_distances_3
template <typename TriangleMesh, typename VertexDistanceMap, typename Mode>
void
estimate_geodesic_distances(const TriangleMesh& tm,
                            VertexDistanceMap vdm,
                            typename boost::graph_traits<TriangleMesh>::vertex_descriptor source,
                            Mode)
{
  CGAL::Heat_method_3::Surface_mesh_geodesic_distances_3<TriangleMesh, Mode> hm(tm);
  hm.add_source(source);
  hm.estimate_geodesic_distances(vdm);
}


#ifndef DOXYGEN_RUNNING
template <typename TriangleMesh, typename VertexDistanceMap>
void
estimate_geodesic_distances(const TriangleMesh& tm,
                            VertexDistanceMap vdm,
                            typename boost::graph_traits<TriangleMesh>::vertex_descriptor source)
{
  CGAL::Heat_method_3::Surface_mesh_geodesic_distances_3<TriangleMesh, Direct> hm(tm);
  hm.add_source(source);
  hm.estimate_geodesic_distances(vdm);
}
#endif


/// \ingroup PkgHeatMethod
/// computes for each vertex  of the triangle mesh `tm` the estimated geodesic distance to a given set of source vertices.
/// \tparam TriangleMesh a triangulated surface mesh, model of `FaceListGraph` and `HalfedgeListGraph`
///         It must have an internal vertex point property map with the value type being a 3D point from a cgal Kernel model
/// \tparam VertexDistanceMap a property map model of `WritablePropertyMap`
///         with `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key type and `double` as value type.
/// \tparam VertexRange a range with value type `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
/// \tparam Mode either the tag `Direct` or `Intrinsic_Delaunay`, which determines if the geodesic distance
///              is computed directly on the mesh, or if the intrinsic Delaunay triangulation is applied first.
///              The default is `Direct`.
///
/// \sa CGAL::Heat_method_3::Surface_mesh_geodesic_distances_3
template <typename TriangleMesh, typename VertexDistanceMap, typename VertexRange, typename Mode>
void
estimate_geodesic_distances(const TriangleMesh& tm,
                            VertexDistanceMap vdm,
                            const VertexRange sources,
                            Mode)
{
  CGAL::Heat_method_3::Surface_mesh_geodesic_distances_3<TriangleMesh, Mode> hm(tm);
  hm.add_sources(sources);
  hm.estimate_geodesic_distances(vdm);
}


} // namespace Heat_method_3
} // namespace CGAL

#include <CGAL/enable_warnings.h>
#endif // CGAL_HEAT_METHOD_3_SURFACE_MESH_GEODESIC_DISTANCES_3_H
