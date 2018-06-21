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


#ifndef CGAL_INTRINSIC_DELAUNAY_TRIANGULATION_3_H
#define CGAL_INTRINSIC_DELAUNAY_TRIANGULATION_3_H

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
#include <CGAL/number_utils.h>
#include <stack>


namespace CGAL {
namespace Intrinsic_Delaunay_Triangulation_3 {

    struct Intrinsic_Delaunay_Triangulation_Eigen_traits_3 {
      typedef Eigen::SparseMatrix<double> SparseMatrix;
      typedef Eigen::Triplet<double> T;
      typedef int Index;
    };


    /**
     * Class `Intrinsic_Delaunay_Triangulation_3` is a ...
     * \tparam TriangleMesh a triangulated surface mesh, model of `FaceGraph` and `HalfedgeListGraph`
     * \tparam Traits a model of IntrinsicDelaunayTriangulation_3
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
               typename EdgeIndexMap = typename boost::property_map< TriangleMesh, boost::edge_index_t>::const_type,
               typename LA = Intrinsic_Delaunay_Triangulation_Eigen_traits_3>
     class Intrinsic_Delaunay_Triangulation_3
     {
       typedef typename boost::graph_traits<TriangleMesh>               graph_traits;
       typedef typename graph_traits::vertex_descriptor            vertex_descriptor;
       typedef typename graph_traits::edge_descriptor                edge_descriptor;
       typedef typename graph_traits::halfedge_descriptor        halfedge_descriptor;
       typedef typename graph_traits::face_descriptor                face_descriptor;
       typedef typename std::set<vertex_descriptor>::iterator        vertex_iterator;
       typedef typename std::set<edge_descriptor>::iterator            edge_iterator;
       /// Geometric typedefs
       typedef typename Traits::Point_3                                      Point_3;
       typedef typename Traits::FT                                                FT;
       typedef typename Traits::Vector_3                                    Vector_3;

       typedef int Index;
       typedef typename boost::property_traits<VertexPointMap>::reference VertexPointMap_reference;

       typedef typename boost::graph_traits<TriangleMesh>::vertices_size_type vertices_size_type;

       typedef CGAL::dynamic_vertex_property_t<Index> Vertex_property_tag;
       typedef typename boost::property_map<TriangleMesh, Vertex_property_tag >::type Vertex_id_map;
       Vertex_id_map vertex_id_map;


       typedef typename boost::graph_traits<TriangleMesh>::faces_size_type face_size_type;

       typedef CGAL::dynamic_face_property_t<Index> Face_property_tag;
       typedef typename boost::property_map<TriangleMesh, Face_property_tag >::type Face_id_map;
       Face_id_map face_id_map;


       typedef CGAL::dynamic_edge_property_t<Index> Edge_property_tag;
       typedef typename boost::property_map<TriangleMesh, Edge_property_tag >::type Edge_id_map;
       Edge_id_map edge_id_map;

       std::stack<edge_iterator, std::list<edge_iterator> > stack;

     public:

      Intrinsic_Delaunay_Triangulation_3(TriangleMesh tm, VertexDistanceMap vdm)
        : tm(tm), vdm(vdm), vpm(get(vertex_point,tm))
      {
        build();
      }



       Intrinsic_Delaunay_Triangulation_3(TriangleMesh tm, VertexDistanceMap vdm, VertexPointMap vpm, FaceIndexMap fpm, EdgeIndexMap epm)
         : tm(tm), vdm(vdm), vpm(vpm), fpm(fpm), epm(epm)
       {
         build();
       }


       //return true if edge is locally delaunay (opposing angles are less than pi)
       bool is_edge_locally_delaunay(edge_descriptor ed)
       {
         //two ways of doing this: taking angles directly (not good with virtual edges)
         //OR: taking edge length and using law of cosines
         //the second way also finds cotan weights which can be used to populate c


         return true;
       }

       //law of cosines
       double new_edge_length(double side_A, double side_B, double angle)
       {
         return CGAL::sqrt(side_A*side_A + side_B*side_B - 2*side_A*side_B*std::cos(angle));
       }

       //Heron's formula
       double face_area(double a, double b, double c)
       {
         double S = (a+b+c)/2;
         return CGAL::sqrt(S*(S-a)*(S-b)*(S-c));
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
         edge_id_map = get(Edge_property_tag(), const_cast<TriangleMesh&>(tm));
         Index edge_i = 0;
         BOOST_FOREACH(edge_descriptor ed, edges(tm)){
           put(edge_id_map, ed, edge_i++);
         }


       }
       TriangleMesh tm;
       VertexPointMap vpm;
       FaceIndexMap fpm;
       EdgeIndexMap epm;
       VertexDistanceMap vdm;





     };


} // namespace Intrinsic_Delaunay_Triangulation_3
} // namespace CGAL
#include <CGAL/enable_warnings.h>
#endif CGAL_HEAT_METHOD_3_HEAT_METHOD_3_H
