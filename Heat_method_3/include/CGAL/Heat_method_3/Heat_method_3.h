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
            typename VertexPointMap = typename boost::property_map< TriangleMesh, vertex_point_t>::type,
            typename LA = Heat_method_Eigen_traits_3>
  class Heat_method_3
  {
    /// Polygon_mesh typedefs
    typedef typename boost::graph_traits<TriangleMesh>               graph_traits;
    typedef typename graph_traits::vertex_descriptor            vertex_descriptor;
    typedef typename graph_traits::edge_descriptor                edge_descriptor;
    typedef typename graph_traits::halfedge_descriptor        halfedge_descriptor;
    typedef typename graph_traits::face_descriptor                face_descriptor;
    typedef typename graph_traits::vertex_iterator                vertex_iterator;
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
  public:

    Heat_method_3(const TriangleMesh& tm)
      : tm(tm), vpm(get(vertex_point,tm))
    {
      build();
    }

    Heat_method_3(const TriangleMesh& tm, VertexPointMap vpm)
      : tm(tm), vpm(vpm)
    {
      build();
    }

    /**
     * add `vd` to the source set, returning `false` if `vd` is already in the set.
     */
    bool add_source(vertex_descriptor vd)
    {
      return sources.insert(vd).second;
    }

    /**
     * remove 'vd' from the source set, returning 'true' if 'vd' was in the set
     */
    bool remove_source(vertex_descriptor vd)
    {
      if(sources.find(vd))
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
    vertex_descriptor get_sources()
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
    vertex_descriptor sources_end()
    {
      return sources.end();
    }

    /**
     * get distance from the current source set to a vertex ` vd`.
     */
    double distance(vertex_descriptor vd)
    {
      return 0;
    }

  private:

    void build()
    {
      vertex_id_map = get(Vertex_property_tag(),const_cast<TriangleMesh&>(tm));
      Index i = 0;
      BOOST_FOREACH(vertex_descriptor vd, vertices(tm)){
        put(vertex_id_map, vd, i++);
      }

      int m = num_vertices(tm);
      Matrix c(m,m);
      Matrix A(m,m);
      std::vector<triplet> A_matrix_entries;
      std::vector<triplet> c_matrix_entries;
      CGAL::Vertex_around_face_iterator<TriangleMesh> vbegin, vend;
      BOOST_FOREACH(face_descriptor f, faces(tm)) {
        boost::tie(vbegin, vend) = vertices_around_face(halfedge(f,tm),tm);
        vertex_descriptor current = *(vbegin);
        vertex_descriptor neighbor_one = *(vbegin++);
        vertex_descriptor neighbor_two = *(vend);
        Index i = get(vertex_id_map, current);
        Index j = get(vertex_id_map, neighbor_one);
        Index k = get(vertex_id_map, neighbor_two);
        Point_3 p_i = get(vpm,current);
        Point_3 p_j = get(vpm, neighbor_one);
        Point_3 p_k = get(vpm, neighbor_two);
        //If the passed in mesh is not a triangle mesh, the algorithm breaks here
        vector cross = CGAL::cross_product((p_j-p_i), (p_k-p_i));
        double dot = to_double((p_j-p_i)*(p_k-p_i));
        double squared_cross = to_double(CGAL::approximate_sqrt(cross*cross));
        double cotan_i = dot/squared_cross;
        c_matrix_entries.push_back(triplet(j,k ,-.5*cotan_i));
        c_matrix_entries.push_back(triplet(k,j,-.5* cotan_i));
        c_matrix_entries.push_back(triplet(j,j,.5* cotan_i));
        c_matrix_entries.push_back(triplet(k,k,.5* cotan_i));

        cross = CGAL::cross_product((p_i-p_j), (p_k-p_j));
        dot = to_double((p_i-p_j)*(p_k-p_j));
        squared_cross = to_double(CGAL::approximate_sqrt(cross*cross));
        double cotan_j = dot/squared_cross;
        c_matrix_entries.push_back(triplet(i,k ,-.5*cotan_j));
        c_matrix_entries.push_back(triplet(k,i,-.5* cotan_j));
        c_matrix_entries.push_back(triplet(i,i,.5* cotan_j));
        c_matrix_entries.push_back(triplet(k,k,.5* cotan_j));

        cross = CGAL::cross_product((p_i-p_k), (p_j-p_k));
        dot = to_double((p_i-p_k)*(p_j-p_k));
        squared_cross = to_double(CGAL::approximate_sqrt(cross*cross));
        double cotan_k = dot/squared_cross;
        c_matrix_entries.push_back(triplet(i,j,-.5*cotan_k));
        c_matrix_entries.push_back(triplet(j,i,-.5* cotan_k));
        c_matrix_entries.push_back(triplet(i,i,.5* cotan_k));
        c_matrix_entries.push_back(triplet(j,j,.5* cotan_k));

        double area_face = CGAL::Polygon_mesh_processing::face_area(f,tm);
        A_matrix_entries.push_back(triplet(i,i, (1./3.)*area_face));
        A_matrix_entries.push_back(triplet(j,j, (1./3.)*area_face));
        A_matrix_entries.push_back(triplet(k,k, (1./3.)*area_face));
      }

      A.setFromTriplets(A_matrix_entries.begin(), A_matrix_entries.end());
      c.setFromTriplets(c_matrix_entries.begin(), c_matrix_entries.end());




    }

    const TriangleMesh& tm;
    VertexPointMap vpm;
    std::set<vertex_descriptor> sources;

    Matrix m;
  };

} // namespace Heat_method_3
} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif CGAL_HEAT_METHOD_3_HEAT_METHOD_3_H
