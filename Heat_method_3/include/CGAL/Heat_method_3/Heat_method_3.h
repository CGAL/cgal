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

#include <CGAL/disable_warnings.h>
#include <set>

#include <CGAL/property_map.h>
#include <CGAL/double.h>
#include <Eigen/Cholesky>
#include <Eigen/Sparse>

#include <boost/foreach.hpp>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/Dynamic_property_map.h>
#include <vector>
#include <CGAL/squared_distance_3.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/number_utils.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/Heat_method_3/Intrinsic_Delaunay_Triangulation_3.h>



namespace CGAL {
namespace Heat_method_3 {

  // This class will later go into another file
  // It encapsulates what we use from Eigen so that one potentially can use another LA library
  struct Heat_method_Eigen_traits_3 {
    typedef Eigen::SparseMatrix<double> SparseMatrix;
    typedef Eigen::Triplet<double> T;
    typedef int Index;
  };


  /**
   * Class `Heat_method_3` is a ...
   * \tparam TriangleMesh a triangulated surface mesh, model of `FaceGraph` and `HalfedgeListGraph`
   * \tparam Traits a model of HeatMethodTraits_3
   * \tparam VertexPointMap a model of `ReadablePropertyMap` with
   *        `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key and
   *        `Traits::Point_3` as value type.
   *        The default is `typename boost::property_map< TriangleMesh, vertex_point_t>::%type`.
   *
   */
  template <typename TriangleMesh,
            typename Traits,
            typename VertexDistanceMap,
            typename VertexPointMap = typename boost::property_map< TriangleMesh, vertex_point_t>::const_type,
            typename FaceIndexMap = typename boost::property_map< TriangleMesh, face_index_t>::const_type,
            typename LA = Heat_method_Eigen_traits_3>
  class Heat_method_3
  {
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

    typedef typename LA::SparseMatrix Matrix;
    typedef typename LA::Index Index;
    typedef typename LA::T triplet;

    // The Vertex_id_map is a property map where you can associate an index to a vertex_descriptor
    typedef typename boost::graph_traits<TriangleMesh>::vertices_size_type vertices_size_type;

    typedef CGAL::dynamic_vertex_property_t<Index> Vertex_property_tag;
    typedef typename boost::property_map<TriangleMesh, Vertex_property_tag >::type Vertex_id_map;
    Vertex_id_map vertex_id_map;

    typedef CGAL::dynamic_face_property_t<Index> Face_property_tag;
    typedef typename boost::property_map<TriangleMesh, Face_property_tag >::type Face_id_map;
    Face_id_map face_id_map;

    //types for intrinsic delaunay triangulation
    typedef CGAL::dynamic_halfedge_property_t<Point_2> Halfedge_coordinate_tag;
    typedef typename boost::property_map<TriangleMesh, Halfedge_coordinate_tag >::type Halfedge_coordinate_map;
    typedef CGAL::dynamic_vertex_property_t<double> Vertex_distance_tag;
    typedef typename boost::property_map<TriangleMesh, Vertex_distance_tag >::type Vertex_distance_map;

    // AF    typedef CGAL::Intrinsic_Delaunay_Triangulation_3::Intrinsic_Delaunay_Triangulation_3<TriangleMesh,Traits, Halfedge_coordinate_map> IDT;

  public:

    Heat_method_3(const TriangleMesh& tm, VertexDistanceMap vdm)
      : vertex_id_map(get(Vertex_property_tag(),const_cast<TriangleMesh&>(tm))), face_id_map(get(Face_property_tag(),tm)), tm(tm), vdm(vdm), vpm(get(vertex_point,tm)), fpm(get(boost::face_index,tm))
    {
      build();
    }

    Heat_method_3(const TriangleMesh& tm, VertexDistanceMap vdm, VertexPointMap vpm, FaceIndexMap fpm)
      : tm(tm), vdm(vdm), vpm(vpm), fpm(fpm)
    {
      build();
    }

    /**
     * add `vd` to the source set, returning `false` if `vd` is already in the set.
     */

     const VertexDistanceMap& get_vertex_distance_map() const
     {
       return vdm;
     }

     const Matrix& mass_matrix() const
     {
       return m_mass_matrix;
     }

     const Matrix& cotan_matrix() const
     {
       return m_cotan_matrix;
     }

     const VertexPointMap& vertex_point_map() const
     {
        return vpm;
     }

     const Vertex_id_map& get_vertex_id_map() const
     {
       return vertex_id_map;
     }

    bool add_source(vertex_descriptor vd)
    {
      source_change_flag = true;
      return sources.insert(vd).second;
    }

    /**
     * remove 'vd' from the source set, returning 'true' if 'vd' was in the set
     */
    bool remove_source(vertex_descriptor vd)
    {
      source_change_flag = true;
      return (sources.erase(vd) == 1);
    }

    /**
     *return current source set
    */
    const std::set<vertex_descriptor>& get_sources() const
    {
      return sources;
    }

    /**
     * clear the current source set
     */
    void clear_sources()
    {
      source_change_flag = true;
      sources.clear();
      return;
    }

    /**
     * return vertex_descriptor to first vertex in the source set
     */
    vertex_iterator sources_begin()
    {
      return sources.begin();
    }
    /**
     * return vertex_descriptor to last vertex in the source set
     */
    vertex_iterator sources_end()
    {
      return sources.end();
    }

    /**
     * get distance from the current source set to a vertex `vd`.
     */
    double distance(vertex_descriptor vd) const
    {
      return 0;
    }

    double summation_of_edges() const
    {
      double edge_sum = 0;
      BOOST_FOREACH(edge_descriptor ed, edges(tm))
      {
        edge_sum += Polygon_mesh_processing::edge_length(halfedge(ed,tm), tm, CGAL::parameters::vertex_point_map(vpm));
      }
      return edge_sum;
    }

    double time_step() const
    {
      return m_time_step;
    }

    Matrix kronecker_delta(const std::set<vertex_descriptor>& sources)
    {
      //currently just working with a single vertex in source set, add the first one for now
      Index i;
      Matrix K(num_vertices(tm), 1);
      if(sources.empty())
      {
        i = 0;
        K.insert(i,0) = 1;
        source_index.insert(0);
      }
      else
      {
        vertex_descriptor current = *(sources.begin());
        i = get(vertex_id_map, current);
        vertex_iterator vd =sources.begin();
        while(vd!=sources.end())
        {

          i = get(vertex_id_map, *(vd));
          K.insert(i,0) = 1;
          vd = ++vd;
        }
      }
      return K;
    }


    const Matrix& kronecker_delta() const
    {
      return kronecker;
    }

    Eigen::VectorXd solve_cotan_laplace(const Matrix& M, const Matrix& c, const Matrix& x, double a_time_step, int dimension) const
    {
      Eigen::VectorXd u;
      Matrix A = (M+ a_time_step*c);
      Eigen::SimplicialLLT<Matrix> solver;
      solver.compute(A);
      if(solver.info()!=Eigen::Success) {
        // decomposition failed
        CGAL_error_msg("Eigen Decomposition in cotan failed");
      }
      u = solver.solve(x);
      if(solver.info()!=Eigen::Success) {
        // solving failed
        CGAL_error_msg("Eigen Solving in cotan failed");
      }
      return u;
    }

    Eigen::MatrixXd compute_unit_gradient(const Eigen::VectorXd& u) const
    {
      Eigen::MatrixXd X(num_faces(tm), 3);
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
        //get area of face_i
        //get outward unit normal
        //cross that with eij, ejk, eki
        //so (Ncross eij) *uk and so on
        //sum all of those then multiply by 1./(2a)
        Vector_3 cross = CGAL::cross_product((p_j-p_i), (p_k-p_i));
        double N_cross = (CGAL::sqrt(cross*cross));
        Vector_3 unit_cross = cross/N_cross;
        double area_face = N_cross * (1./2);
        Vector_3 edge_sums = u(k) * CGAL::cross_product(unit_cross,(p_j-p_i));
        edge_sums = edge_sums + u(i) * (CGAL::cross_product(unit_cross, (p_k-p_j)));
        edge_sums = edge_sums + u(j) * CGAL::cross_product(unit_cross, (p_i-p_k));
        edge_sums = edge_sums * (1./area_face);
        double e_magnitude = CGAL::sqrt(edge_sums*edge_sums);
        Vector_3 unit_grad = edge_sums*(1./e_magnitude);
        X(face_i, 0) = unit_grad.x();
        X(face_i, 1) = unit_grad.y();
        X(face_i, 2) = unit_grad.z();
      }
      return X;
    }

    double dot_eigen_vector(const Eigen::Vector3d& a, const Vector_3& b) const
    {
      return (a(0)*b.x() + a(1)*b.y() + a(2)*b.z());
    }


    Matrix compute_divergence(const Eigen::MatrixXd& X, int rows) const
    {
      Matrix indexD;
      std::vector<triplet> d_matrix_entries;
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

        Eigen::VectorXd a = X.row(face_i);
        double i_entry = cotan_k*(dot_eigen_vector(a,(p_j-p_i))) + cotan_j*(dot_eigen_vector(a,(p_k-p_i)));
        double j_entry = cotan_i*(dot_eigen_vector(a,(p_k-p_j))) + cotan_k*(dot_eigen_vector(a,(p_i-p_j)));
        double k_entry = cotan_j*(dot_eigen_vector(a,(p_i-p_k))) + cotan_i*(dot_eigen_vector(a,(p_j-p_k)));

        d_matrix_entries.push_back(triplet(i,0, (1./2)*i_entry));
        d_matrix_entries.push_back(triplet(j,0, (1./2)*j_entry));
        d_matrix_entries.push_back(triplet(k,0, (1./2)*k_entry));
      }
      indexD.resize(rows,1);
      indexD.setFromTriplets(d_matrix_entries.begin(), d_matrix_entries.end());
      return indexD;
    }

    Eigen::VectorXd value_at_source_set(const Eigen::VectorXd& phi, int dimension) const
    {
      Eigen::VectorXd source_set_val(dimension,1);
      if(sources.empty())
      {
        for(int k = 0; k<dimension; k++)
        {
          source_set_val(k,0) = phi.coeff(0,0);
        }
        source_set_val = phi-source_set_val;
      }
      else
      {
        for(int i = 0; i<dimension; i++)
        {
          double min_val = INT_MAX;
          vertex_iterator current;
          Index current_Index;
          current = sources.begin();
          //go through the distances to the sources and leave the minimum distance;
          for(int j = 0; j<sources.size(); j++)
          {
            current_Index = get(vertex_id_map, *current);
            double new_d = fabs(-phi.coeff(current_Index,0)+phi.coeff(i,0));
            if(new_d < min_val)
            {
              min_val = new_d;
            }
            current = ++current;
          }
          source_set_val(i,0) = min_val;
        }
      }
      return source_set_val;
    }


    Eigen::VectorXd solve_phi(Matrix c, Matrix divergence, int dimension) const
    {

      Eigen::VectorXd phi;
      Eigen::SimplicialLDLT<Matrix> solver;
      solver.compute(c);
      if(solver.info()!=Eigen::Success) {
        // decomposition failed
        CGAL_error_msg("Eigen Decomposition in phi failed");
      }
      phi = solver.solve(divergence);
      if(solver.info()!=Eigen::Success) {
        // solving failed
        CGAL_error_msg("Eigen Solving in phi failed");
      }
      return value_at_source_set(phi, dimension);
    }

    //currently, this function will return a (number of vertices)x1 vector where
    //the ith index has the distance from the first vertex to the ith vertex
    const Eigen::VectorXd& distances() const
    {
      return solved_phi;
    }


    void update()
    {
      double d=0;
      if(source_change_flag)
      {
        //don't need to recompute Mass matrix, cotan matrix or timestep reflect that in this function
        sources = get_sources();
        kronecker = kronecker_delta(sources);
        solved_u = solve_cotan_laplace(m_mass_matrix, m_cotan_matrix, kronecker, m_time_step, dimension);
        X = compute_unit_gradient(solved_u);
        index_divergence = compute_divergence(X, dimension);
        solved_phi = solve_phi(m_cotan_matrix, index_divergence, dimension);
        source_change_flag = false;
        std::cout<<"sources changed, recompute\n";
      }
      BOOST_FOREACH(vertex_descriptor vd, vertices(tm)){
        Index i_d = get(vertex_id_map, vd);
        d = solved_phi(i_d,0);
        std::cout<<d<<"\n";
        put(vdm,vd,d);
      }
    }

  private:

    void build()
    {
      source_change_flag = false;
      /*
      TriangleMesh idt_copy;
      if(idf)
      {
        CGAL::copy_face_graph(tm, idt_copy);
        std::cout<<"copied graph\n";
        halfedge_coord_map = get(Halfedge_coordinate_tag(), idt_copy);
        AF        IDT im(idt_copy, halfedge_coord_map);
      }
      */

      CGAL_precondition(is_triangle_mesh(tm));
      //vertex_id_map = get(Vertex_property_tag(),const_cast<TriangleMesh&>(tm));
      Index i = 0;
      BOOST_FOREACH(vertex_descriptor vd, vertices(tm)){
        put(vertex_id_map, vd, i++);
      }
      //face_id_map = get(Face_property_tag(), const_cast<TriangleMesh&>(tm));
      Index face_i = 0;
      BOOST_FOREACH(face_descriptor fd, faces(tm)){
        put(face_id_map, fd, face_i++);
      }
      int m = static_cast<int>(num_vertices(tm));
      dimension = m;
      //mass matrix entries
      std::vector<triplet> A_matrix_entries;
      //cotan matrix entries
      std::vector<triplet> c_matrix_entries;
      CGAL::Vertex_around_face_iterator<TriangleMesh> vbegin, vend, vmiddle;
      //Go through each face on the mesh
      //if(!idf)
      {
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
          c_matrix_entries.push_back(triplet(j,k ,-(1./2)*cotan_i));
          c_matrix_entries.push_back(triplet(k,j,-(1./2)* cotan_i));
          c_matrix_entries.push_back(triplet(j,j,(1./2)*cotan_i));
          c_matrix_entries.push_back(triplet(k,k,(1./2)* cotan_i));

          cross = CGAL::cross_product((pi-pj), (pk-pj));
          dot = to_double((pi-pj)*(pk-pj));
          double cotan_j = dot/norm_cross;
          c_matrix_entries.push_back(triplet(i,k ,-(1./2)*cotan_j));
          c_matrix_entries.push_back(triplet(k,i,-(1./2)* cotan_j));
          c_matrix_entries.push_back(triplet(i,i,(1./2)* cotan_j));
          c_matrix_entries.push_back(triplet(k,k,(1./2)* cotan_j));

          cross = CGAL::cross_product((pi-pk), (pj-pk));
          dot = to_double((pi-pk)*(pj-pk));
          double cotan_k = dot/norm_cross;
          c_matrix_entries.push_back(triplet(i,j,-(1./2)*cotan_k));
          c_matrix_entries.push_back(triplet(j,i,-(1./2)* cotan_k));
          c_matrix_entries.push_back(triplet(i,i,(1./2)* cotan_k));
          c_matrix_entries.push_back(triplet(j,j,(1./2)* cotan_k));

          //double area_face = CGAL::Polygon_mesh_processing::face_area(f,tm);
          //cross is 2*area
          A_matrix_entries.push_back(triplet(i,i, (1./6.)*norm_cross));
          A_matrix_entries.push_back(triplet(j,j, (1./6.)*norm_cross));
          A_matrix_entries.push_back(triplet(k,k, (1./6.)*norm_cross));
        }
      }
      /*
      else
      {
        BOOST_FOREACH(face_descriptor f, faces(idt_copy)) {
          boost::tie(vbegin, vend) = vertices_around_face(halfedge(f,tm),tm);
          vertex_descriptor current = *(vbegin);
          vertex_descriptor neighbor_one = *(++vbegin);
          vertex_descriptor neighbor_two = *(++(vbegin));
          Index i = get(vertex_id_map, current);
          Index j = get(vertex_id_map, neighbor_one);
          Index k = get(vertex_id_map, neighbor_two);
          halfedge_descriptor first_h = next(halfedge(f, idt_copy), idt_copy);
          halfedge_descriptor second_h = next(first_h, idt_copy);
          halfedge_descriptor third_h = next(second_h, idt_copy);

          //add a check to make sure halfedge->face is the face we are looking at
          Point_2 p_i = get(halfedge_coord_map, first_h);
          Point_2 p_j = get(halfedge_coord_map, second_h);
          Point_2 p_k = get(halfedge_coord_map, third_h);
          Point_3 pi(p_i.x(), p_i.y(),0);
          Point_3 pj(p_j.x(), p_j.y(),0);
          Point_3 pk(p_k.x(), p_k.y(),0);

          Vector_3 cross = CGAL::cross_product((pj-pi), (pk-pi));
          double dot = (pj-pi)*(pk-pi);
          double norm_cross = (CGAL::sqrt(cross*cross));

          double cotan_i = dot/norm_cross;
          c_matrix_entries.push_back(triplet(j,k ,-(1./2)*cotan_i));
          c_matrix_entries.push_back(triplet(k,j,-(1./2)* cotan_i));
          c_matrix_entries.push_back(triplet(j,j,(1./2)*cotan_i));
          c_matrix_entries.push_back(triplet(k,k,(1./2)* cotan_i));

          cross = CGAL::cross_product((pi-pj), (pk-pj));
          dot = to_double((pi-pj)*(pk-pj));
          double cotan_j = dot/norm_cross;
          c_matrix_entries.push_back(triplet(i,k ,-(1./2)*cotan_j));
          c_matrix_entries.push_back(triplet(k,i,-(1./2)* cotan_j));
          c_matrix_entries.push_back(triplet(i,i,(1./2)* cotan_j));
          c_matrix_entries.push_back(triplet(k,k,(1./2)* cotan_j));

          cross = CGAL::cross_product((pi-pk), (pj-pk));
          dot = to_double((pi-pk)*(pj-pk));
          double cotan_k = dot/norm_cross;
          c_matrix_entries.push_back(triplet(i,j,-(1./2)*cotan_k));
          c_matrix_entries.push_back(triplet(j,i,-(1./2)* cotan_k));
          c_matrix_entries.push_back(triplet(i,i,(1./2)* cotan_k));
          c_matrix_entries.push_back(triplet(j,j,(1./2)* cotan_k));
          //double area_face = CGAL::Polygon_mesh_processing::face_area(f,tm);
          //cross is 2*area
          A_matrix_entries.push_back(triplet(i,i, (1./6.)*norm_cross));
          A_matrix_entries.push_back(triplet(j,j, (1./6.)*norm_cross));
          A_matrix_entries.push_back(triplet(k,k, (1./6.)*norm_cross));
        
   } 
      */
      m_mass_matrix.resize(m,m);
      m_mass_matrix.setFromTriplets(A_matrix_entries.begin(), A_matrix_entries.end());
      m_cotan_matrix.resize(m,m);
      m_cotan_matrix.setFromTriplets(c_matrix_entries.begin(), c_matrix_entries.end());

      m_time_step = 1./(num_edges(tm));
      m_time_step = m_time_step*summation_of_edges();

      sources = get_sources();
      kronecker = kronecker_delta(sources);
      solved_u = solve_cotan_laplace(m_mass_matrix, m_cotan_matrix, kronecker, m_time_step, m);
      X = compute_unit_gradient(solved_u);
      index_divergence = compute_divergence(X, m);
      solved_phi = solve_phi(m_cotan_matrix, index_divergence, m);
    }
    
    int dimension;
    const TriangleMesh& tm;
    VertexDistanceMap vdm;
    VertexPointMap vpm;
    FaceIndexMap fpm;
    std::set<vertex_descriptor> sources;
    double m_time_step;
    Matrix kronecker;
    Matrix m_mass_matrix, m_cotan_matrix;
    Eigen::VectorXd solved_u;
    Eigen::MatrixXd X;
    Matrix index_divergence;
    Eigen::VectorXd solved_phi;
    std::set<Index> source_index;
    bool source_change_flag;
    // AF    Halfedge_coordinate_map halfedge_coord_map;
  };

} // namespace Heat_method_3
} // namespace CGAL
#include <CGAL/enable_warnings.h>
#endif // CGAL_HEAT_METHOD_3_HEAT_METHOD_3_H
