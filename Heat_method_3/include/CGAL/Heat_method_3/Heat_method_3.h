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

#include <Eigen/Cholesky>
#include <Eigen/Sparse>

#include <boost/foreach.hpp>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Dynamic_property_map.h>
#include <vector>
#include <CGAL/Vector_3.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/number_utils.h>



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
            typename VertexPointMap = typename boost::property_map< TriangleMesh, vertex_point_t>::const_type,
            typename FacePointMap = typename boost::property_map< TriangleMesh, face_index_t>::type,
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

    typedef typename Simple_cartesian<double>::Vector_3                    vector;


    typedef typename LA::SparseMatrix Matrix;
    typedef typename LA::Index Index;
    typedef typename LA::T triplet;

    // The Vertex_id_map is a property map where you can associate an index to a vertex_descriptor
    typedef typename boost::graph_traits<TriangleMesh>::vertices_size_type vertices_size_type;

    typedef CGAL::dynamic_vertex_property_t<Index> Vertex_property_tag;
    typedef typename boost::property_map<TriangleMesh, Vertex_property_tag >::type Vertex_id_map;
    Vertex_id_map vertex_id_map;


    typedef typename boost::graph_traits<TriangleMesh>::faces_size_type face_size_type;

    typedef CGAL::dynamic_face_property_t<Index> Face_property_tag;
    typedef typename boost::property_map<TriangleMesh, Face_property_tag >::type Face_id_map;
    Face_id_map face_id_map;

  public:

    Heat_method_3(const TriangleMesh& tm)
      : tm(tm), vpm(get(vertex_point,tm))
    {
      build();
    }

    Heat_method_3(const TriangleMesh& tm, VertexPointMap vpm, FacePointMap fpm)
      : tm(tm), vpm(vpm), fpm(fpm)
    {
      build();
    }

    /**
     * add `vd` to the source set, returning `false` if `vd` is already in the set.
     */

     const Matrix& get_mass_matrix()
     {
       return mass_matrix;
     }

     const Matrix& get_cotan_matrix()
     {
       return cotan_matrix;
     }

    bool add_source(vertex_descriptor vd)
    {
      return sources.insert(vd).second;
    }

    /**
     * remove 'vd' from the source set, returning 'true' if 'vd' was in the set
     */
    bool remove_source(vertex_descriptor vd)
    {
      return (sources.erase(vd) == 1);

      // The following code is not only longer but has a bug.
      // sources.end() does not return an iterator to the last element but is a past-the-end iterator
      // that is you cannot dereference it.
      // sources.find(vd) returns the past-the-end iterator of vd cannot be found
      // Otherwise it returns the iterator it with *it == vd
      // So the test would be if(sources.find() != sources,end()
      if((sources.find(vd)) != (sources.end()) || vd == *(sources.end()))
      {
        sources.erase(vd);
        return true;
      }
      else
      {
        return false;
      }
    }

    /**
     *return current source set
    */
    const std::set<vertex_descriptor>& get_sources()
    {
      return sources;
    }

    /**
     * clear the current source set
     */
    void clear_sources()
    {
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
    double distance(vertex_descriptor vd)
    {
      return 0;
    }

    double summation_of_edges()
    {
      double edge_sum = 0;
      BOOST_FOREACH(edge_descriptor ed, edges(tm))
      {
        edge_sum += Polygon_mesh_processing::edge_length(halfedge(ed,tm), tm);
      }
      return edge_sum;
    }

    double get_time_step()
    {
      return time_step;
    }

    Matrix kronecker_delta(std::set<vertex_descriptor> sources)
    {
      //currently just working with a single vertex in source set, add the first one for now
      Index i;
      if(sources.empty())
      {
        add_source(*(vertices(tm).begin()));
        i = 0;
      }
      else
      {
        i = get(vertex_id_map, *(sources.begin()));
      }
      Matrix K(num_vertices(tm), 1);
      K.insert(i,0) = 1;
      return K;
    }


    const Matrix& get_kronecker_delta()
    {
      return kronecker;
    }



    Eigen::VectorXd solve_cotan_laplace(Matrix M, Matrix c, Matrix x, double time_step, int dimension)
    {
      Eigen::VectorXd u;
      Matrix A = (M+ time_step*c);
      Eigen::SimplicialLLT<Matrix> solver;
      solver.compute(A);
      if(solver.info()!=Eigen::Success) {
        // decomposition failed
        CGAL_error_msg("Eigen Decomposition failed");
        CGAL_error();
      }
      u = solver.solve(x);
      if(solver.info()!=Eigen::Success) {
        // solving failed
        CGAL_error_msg("Eigen Solving failed");
        CGAL_error();
      }
      // solve for another right hand side:
      return u;
    }

  private:

    void build()
    {
      CGAL_precondition(is_triangle_mesh(tm));

      vertex_id_map = get(Vertex_property_tag(),const_cast<TriangleMesh&>(tm));
      Index i = 0;
      BOOST_FOREACH(vertex_descriptor vd, vertices(tm)){
        put(vertex_id_map, vd, i++);
      }
      face_id_map = get(Face_property_tag(), const_cast<TriangleMesh&>(tm));
      Index face_i = 0;
      BOOST_FOREACH(face_descriptor fd, faces(tm)){
        put(face_id_map, fd, face_i++);
      }

      int m = static_cast<int>(num_vertices(tm));

      //mass matrix entries
      std::vector<triplet> A_matrix_entries;
      //cotan matrix entries
      std::vector<triplet> c_matrix_entries;
      CGAL::Vertex_around_face_iterator<TriangleMesh> vbegin, vend, vmiddle;
      //Go through each face on the mesh
      BOOST_FOREACH(face_descriptor f, faces(tm)) {
        //Prior assumption that it is a triangle mesh
        boost::tie(vbegin, vend) = vertices_around_face(halfedge(f,tm),tm);
        vertex_descriptor current = *(vbegin);
        vertex_descriptor neighbor_one = *(++vbegin);
        vertex_descriptor neighbor_two = *(++vbegin);
        Index i = get(vertex_id_map, current);
        Index j = get(vertex_id_map, neighbor_one);
        Index k = get(vertex_id_map, neighbor_two);
        Point_3 p_i = get(vpm,current);
        Point_3 p_j = get(vpm, neighbor_one);
        Point_3 p_k = get(vpm, neighbor_two);

        //If the passed in mesh is not a triangle mesh, the algorithm breaks here
        vector cross = CGAL::cross_product((p_j-p_i), (p_k-p_i));
        double dot = (p_j-p_i)*(p_k-p_i);

        double norm_cross = (CGAL::sqrt(cross*cross));

        double cotan_i = dot/norm_cross;
        c_matrix_entries.push_back(triplet(j,k ,-(1./2)*cotan_i));
        c_matrix_entries.push_back(triplet(k,j,-(1./2)* cotan_i));
        c_matrix_entries.push_back(triplet(j,j,(1./2)*cotan_i));
        c_matrix_entries.push_back(triplet(k,k,(1./2)* cotan_i));

        cross = CGAL::cross_product((p_i-p_j), (p_k-p_j));
        dot = to_double((p_i-p_j)*(p_k-p_j));
        double cotan_j = dot/norm_cross;
        c_matrix_entries.push_back(triplet(i,k ,-(1./2)*cotan_j));
        c_matrix_entries.push_back(triplet(k,i,-(1./2)* cotan_j));
        c_matrix_entries.push_back(triplet(i,i,(1./2)* cotan_j));
        c_matrix_entries.push_back(triplet(k,k,(1./2)* cotan_j));

        cross = CGAL::cross_product((p_i-p_k), (p_j-p_k));
        dot = to_double((p_i-p_k)*(p_j-p_k));
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
      mass_matrix.resize(m,m);
      mass_matrix.setFromTriplets(A_matrix_entries.begin(), A_matrix_entries.end());
      cotan_matrix.resize(m,m);
      cotan_matrix.setFromTriplets(c_matrix_entries.begin(), c_matrix_entries.end());

      time_step = 1./(num_edges(tm));
      time_step = time_step*summation_of_edges();

      // AF: This segfaults as sources is empty
      sources = get_sources();
      kronecker = kronecker_delta(sources);
      solved_u = solve_cotan_laplace(mass_matrix, cotan_matrix, kronecker, time_step, m);


    }




    const TriangleMesh& tm;
    VertexPointMap vpm;
    FacePointMap fpm;
    std::set<vertex_descriptor> sources;
    double time_step;
    Matrix kronecker;
    Eigen::VectorXd solved_u;
    Matrix mass_matrix, cotan_matrix;
  };

} // namespace Heat_method_3
} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif CGAL_HEAT_METHOD_3_HEAT_METHOD_3_H
