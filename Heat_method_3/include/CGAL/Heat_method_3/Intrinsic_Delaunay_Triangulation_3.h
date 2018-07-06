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
#include <CGAL/Surface_mesh.h>
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
#include <CGAL/boost/graph/helpers.h>
#include <stack>
#include <boost/graph/filtered_graph.hpp>
#include <fstream>
#include <array>
#include <math.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_parameterization/IO/File_off.h>
#include <CGAL/Surface_mesh_parameterization/Circular_border_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Discrete_authalic_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Error_code.h>
#include <CGAL/Surface_mesh_parameterization/parameterize.h>
#include <CGAL/Polygon_mesh_processing/measure.h>




namespace CGAL {
namespace Intrinsic_Delaunay_Triangulation_3 {

    struct Intrinsic_Delaunay_Triangulation_Eigen_traits_3 {
      typedef Eigen::SparseMatrix<double> SparseMatrix;
      typedef Eigen::Triplet<double> T;
      typedef int Index;
    };


  template <typename IDT>
  struct IDT_vertex_point_property_map;
  
  
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
               typename HalfedgeCoordinateMap,
               typename VertexPointMap = typename boost::property_map< TriangleMesh, vertex_point_t>::const_type,
               typename FaceIndexMap = typename boost::property_map< TriangleMesh, face_index_t>::const_type,
               typename EdgeIndexMap = typename boost::property_map< TriangleMesh, boost::edge_index_t>::const_type,
               typename LA = Intrinsic_Delaunay_Triangulation_Eigen_traits_3>
     class Intrinsic_Delaunay_Triangulation_3
     {
       typedef Intrinsic_Delaunay_Triangulation_3<TriangleMesh,Traits,HalfedgeCoordinateMap,VertexPointMap,FaceIndexMap,EdgeIndexMap,LA> Self;
       
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

       typedef typename Traits::Point_2                                     Point_2;

       typedef typename LA::SparseMatrix Matrix;
       typedef typename LA::Index Index;

       typedef typename boost::property_traits<VertexPointMap>::reference VertexPointMap_reference;

       typedef typename boost::graph_traits<TriangleMesh>::vertices_size_type vertices_size_type;

       typedef CGAL::dynamic_vertex_property_t<Index> Vertex_property_tag;
       typedef typename boost::property_map<TriangleMesh, Vertex_property_tag >::type Vertex_id_map;
       Vertex_id_map vertex_id_map;

       typedef CGAL::dynamic_face_property_t<Index> Face_property_tag;
       typedef typename boost::property_map<TriangleMesh, Face_property_tag >::type Face_id_map;
       Face_id_map face_id_map;

       typedef CGAL::dynamic_edge_property_t<Index> Edge_property_tag;
       typedef typename boost::property_map<TriangleMesh, Edge_property_tag >::type Edge_id_map;
       Edge_id_map edge_id_map;
       typedef typename std::stack<edge_descriptor, std::list<edge_descriptor> > edge_stack;

       friend struct IDT_vertex_point_property_map<Self>;

       
     public: // for the BGL functions below. They should maybe become friend?
       struct Vertex_descriptor {
         halfedge_descriptor hd;

         Vertex_descriptor(const halfedge_descriptor& hd)
           : hd(hd)
         {}
       };

     public:
       Intrinsic_Delaunay_Triangulation_3(TriangleMesh& tm, HalfedgeCoordinateMap hcm)
        : tm(tm), hcm(hcm)
        {
          build();
        }



       Intrinsic_Delaunay_Triangulation_3(TriangleMesh& tm, HalfedgeCoordinateMap hcm, FaceIndexMap fpm, EdgeIndexMap epm)
         : tm(tm), hcm(hcm), fpm(fpm), epm(epm)
       {
         build();
       }


       typedef TriangleMesh Triangle_mesh;
       
       const Triangle_mesh& triangle_mesh() const
       {
         return tm;
       }

       Triangle_mesh& triangle_mesh()
       {
         return tm;
       }

       const HalfedgeCoordinateMap& hcmap() const
       {
         return hcm;
       }
       
       double get_cotan_weight(edge_descriptor ed)
       {
          double cotan_weight = 0;
          halfedge_descriptor hd = halfedge(ed, tm);
          halfedge_descriptor hd2 = next(hd,tm);
          halfedge_descriptor hd3 = next(hd2,tm);
          Index a_i = get(edge_id_map, ed);
          Index b_i = get(edge_id_map, edge(hd2,tm));
          Index c_i = get(edge_id_map, edge(hd3,tm));
          double a = edge_lengths(a_i,0) + 0.0;
          double b = edge_lengths(b_i,0) + 0.0;
          double c = edge_lengths(c_i,0) + 0.0;

          double tan2 = CGAL::sqrt(CGAL::abs(((a-b+c)*(a+b-c))/((a+b+c)*(-a+b+c))));
          cotan_weight+=(1-(tan2*tan2))/(2*tan2);

          hd = opposite(hd,tm);
          hd2 =next(hd,tm);
          hd3 = next(hd2,tm);
          b_i = get(edge_id_map, edge(hd2,tm));
          c_i = get(edge_id_map, edge(hd3,tm));
          b = edge_lengths(b_i,0) + 0.0;
          c = edge_lengths(c_i,0) + 0.0;
          tan2 = CGAL::sqrt(CGAL::abs(((a-b+c)*(a+b-c))/((a+b+c)*(-a+b+c))));
          cotan_weight+=(1-(tan2*tan2))/(2*tan2);
          return cotan_weight;
        }

         //return true if edge is locally delaunay (opposing angles are less than pi)
        bool is_edge_locally_delaunay(edge_descriptor ed)
        {
           //two ways of doing this: taking angles directly (not good with virtual edges)
           //OR: taking edge length and using law of cosines
           //the second way checks cotan weights
           if(get_cotan_weight(ed)>=0)
           {
             return true;
           }
          else
          {
            return false;
          }
        }

        void change_edge_length(Index i, edge_descriptor ed)
        {
          halfedge_descriptor hd = halfedge(ed,tm);
          halfedge_descriptor hd2 = next(hd,tm);
          halfedge_descriptor hd3 = next(hd2,tm);
          Index b_i = get(edge_id_map, edge(hd2,tm));
          Index c_i = get(edge_id_map, edge(hd3,tm));
          double a = edge_lengths(i,0);
          double b1 = edge_lengths(b_i,0);
          double c1 = edge_lengths(c_i,0);
          double tan2a = CGAL::sqrt(CGAL::abs(((c1-a+b1)*(-b1+a+c1))/((a+b1+c1)*(b1+a-c1))));
          hd = opposite(hd,tm);
          hd2 =next(hd,tm);
          hd3 = next(hd2,tm);
          b_i = get(edge_id_map, edge(hd2,tm));
          c_i = get(edge_id_map, edge(hd3,tm));
          double b2 = edge_lengths(b_i,0);
          double c2 = edge_lengths(c_i,0);
          double tan2d = CGAL::sqrt(CGAL::abs(((-a+b2+c2)*(a+b2-c2))/((a+b2+c2)*(a-b2+c2))));
          double tan2ad = (tan2a + tan2d)/(1-tan2a*tan2d);
          double cosad = (1-tan2ad*tan2ad)/(1+tan2ad*tan2ad);
          double new_length = CGAL::sqrt( CGAL::abs(b1*b1 + c2*c2 - 2*b1*c2*cosad));
          edge_lengths(i,0) = new_length;
        }


         //Heron's formula
        double face_area(double a, double b, double c)
        {
          double S = (a+b+c)/2;
          return CGAL::sqrt(S*(S-a)*(S-b)*(S-c));
        }

        void loop_over_edges(edge_stack stack, Eigen::VectorXd marked_edges)
        {
          int a = 0;
          while(!stack.empty())
          {
            edge_descriptor ed = stack.top();
            stack.pop();

            Index edge_i = get(edge_id_map,ed);

            marked_edges(edge_i,0)=0;
            //if the edge itself is not locally delaunay, go back
            if(!(is_edge_locally_delaunay(ed)))
            {
              if(!(is_border(ed,tm)))
              {
               a++;
               change_edge_length(edge_i,ed);
               halfedge_descriptor hd = (halfedge(ed, tm));
               CGAL::Euler::flip_edge(hd, tm);
               edge_descriptor next_edge= edge(next(hd,tm),tm);
               Index next_edge_i =  get(edge_id_map, next_edge);

               //if edge was already checked, go back and check again
               //for the 4 surrounding edges, since local 'geometry' changed,
               if(!(marked_edges(next_edge_i,0)))
               {
                 stack.push(next_edge);
                 marked_edges(next_edge_i,0) = 1;
               }
               next_edge = edge(prev(hd,tm),tm);
               next_edge_i = get(edge_id_map,next_edge);
               if(!(marked_edges(next_edge_i,0)))
               {
                 stack.push(next_edge);
                 marked_edges(next_edge_i,0) = 1;
               }
               next_edge = edge(next(opposite(hd,tm),tm),tm);
               next_edge_i = get(edge_id_map,next_edge);
               if(!(marked_edges(next_edge_i,0)))
               {
                 stack.push(next_edge);
                 marked_edges(next_edge_i,0) = 1;
               }
               next_edge = edge(prev(opposite(hd,tm),tm),tm);
               next_edge_i = get(edge_id_map,next_edge);
               if(!(marked_edges(next_edge_i,0)))
               {
                 stack.push(next_edge);
                 marked_edges(next_edge_i,0) = 1;
               }
             }
               //then go back to top of the stack
            }
          }
          std::cout<<"this many edges were flipped: "<< a <<"\n";
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
           edge_stack stack;
           edge_lengths.resize(number_of_edges, 1);
           mark_edges.resize(number_of_edges, 1);
           edge_id_map = get(Edge_property_tag(), const_cast<TriangleMesh&>(tm));
           Index edge_i = 0;
           BOOST_FOREACH(edge_descriptor ed, edges(tm)){
             mark_edges(edge_i,0)=1;
             edge_lengths(edge_i,0) = Polygon_mesh_processing::edge_length(halfedge(ed,tm),tm);
             put(edge_id_map, ed, edge_i++);
             stack.push(ed);
           }
           loop_over_edges(stack, mark_edges);
           //now that edges are calculated, go through and for each face, calculate the vertex positions around it

           BOOST_FOREACH(face_descriptor f, faces(tm))
           {
             CGAL::Vertex_around_face_iterator<TriangleMesh> vbegin, vend, vmiddle;
             Index face_i = get(face_id_map, f);

             boost::tie(vbegin, vend) = vertices_around_face(halfedge(f,tm),tm);
             halfedge_descriptor hd = halfedge(f,tm);
             if(face(hd,tm) != f)
             {
               hd = opposite(hd,tm);
             }
             hd = next(hd,tm);
             //each 'local' set of coordinates will have 0,0 at the first vertex/halfedge
             Point_2 p11(0,0);
             put(hcm, hd,p11);
             edge_descriptor ed1 = edge(hd, tm);
             hd = next(hd,tm);
             //the second local coordinate will be edge_length(first edge),0
             Point_2 p21(edge_lengths(get(edge_id_map,ed1),0),0);
             put(hcm,hd,p21);


             //use basic trigonometry to compute third coordinate
             edge_descriptor ed2 = edge(hd, tm);
             hd = next(hd,tm);
             edge_descriptor ed3 = edge(hd, tm);
             Index e1 = get(edge_id_map, ed1);
             Index e2 = get(edge_id_map, ed2);
             Index e3 = get(edge_id_map, ed3);
             double e1_len = edge_lengths(e1,0);
             double e2_len = edge_lengths(e2,0);
             double e3_len = edge_lengths(e3,0);
             double angle_a = -(e2_len*e2_len) + e3_len*e3_len + e1_len*e1_len;
             angle_a = angle_a/(2*e3_len*e1_len);
             Point_2 p31(e3_len*std::cos(angle_a), e3_len*std::sin(angle_a));
             put(hcm,hd,p31);

           }
         }
         //todo:: determine which can be const
         TriangleMesh& tm;
         VertexPointMap vpm;
         FaceIndexMap fpm;
         EdgeIndexMap epm;
         HalfedgeCoordinateMap hcm;
         int number_of_edges = num_edges(tm);
         int number_of_faces = num_faces(tm);
         Eigen::VectorXd edge_lengths;
         Eigen::VectorXd mark_edges;
      };

} // namespace Intrinsic_Delaunay_Triangulation_3
} // namespace CGAL

namespace boost {
  
template <typename TM,
          typename T,
          typename HCM,
          typename VPM,
          typename FIM,
          typename EIM,
           typename LA>
struct graph_traits<CGAL::Intrinsic_Delaunay_Triangulation_3::Intrinsic_Delaunay_Triangulation_3<TM,T,HCM,VPM,FIM,EIM,LA> > {

  typedef CGAL::Intrinsic_Delaunay_Triangulation_3::Intrinsic_Delaunay_Triangulation_3<TM,T,HCM,VPM,FIM,EIM,LA> Mesh;
  typedef typename Mesh::Vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<TM>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<TM>::halfedge_iterator halfedge_iterator;
  typedef typename boost::graph_traits<TM>::edge_descriptor edge_descriptor;
  typedef typename boost::graph_traits<TM>::edge_iterator edge_iterator;
  typedef typename boost::graph_traits<TM>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<TM>::face_iterator face_iterator;

  typedef typename boost::graph_traits<TM>::vertices_size_type vertices_size_type;
};
  
} // namespace boost


namespace CGAL {
namespace Intrinsic_Delaunay_Triangulation_3 {
    
template <typename TM,
          typename T,
          typename HCM,
          typename VPM,
          typename FIM,
          typename EIM,
          typename LA>
typename boost::graph_traits<TM>::vertices_size_type
num_vertices(const Intrinsic_Delaunay_Triangulation_3<TM,T,HCM,VPM,FIM,EIM,LA>& idt)
{
  return num_vertices(idt.triangle_mesh());
}
  
template <typename TM,
          typename T,
          typename HCM,
          typename VPM,
          typename FIM,
          typename EIM,
          typename LA>
typename std::pair<typename boost::graph_traits<TM>::vertex_iterator,
                   typename boost::graph_traits<TM>::vertex_iterator>
vertices(const Intrinsic_Delaunay_Triangulation_3<TM,T,HCM,VPM,FIM,EIM,LA>& idt)
 {
   return vertices(idt.triangle_mesh());
 }
  
template <typename TM,
          typename T,
          typename HCM,
          typename VPM,
          typename FIM,
          typename EIM,
          typename LA>
typename std::pair<typename boost::graph_traits<TM>::halfedge_iterator,
                   typename boost::graph_traits<TM>::halfedge_iterator>
halfedges(const Intrinsic_Delaunay_Triangulation_3<TM,T,HCM,VPM,FIM,EIM,LA>& idt)
 {
   return halfedges(idt.triangle_mesh());
 }

  
template <typename TM,
          typename T,
          typename HCM,
          typename VPM,
          typename FIM,
          typename EIM,
          typename LA>
typename std::pair<typename boost::graph_traits<TM>::edge_iterator,
                   typename boost::graph_traits<TM>::edge_iterator>
edges(const Intrinsic_Delaunay_Triangulation_3<TM,T,HCM,VPM,FIM,EIM,LA>& idt)
 {
   return edges(idt.triangle_mesh());
 }

  
template <typename TM,
          typename T,
          typename HCM,
          typename VPM,
          typename FIM,
          typename EIM,
          typename LA>
typename boost::graph_traits<Intrinsic_Delaunay_Triangulation_3<TM,T,HCM,VPM,FIM,EIM,LA> >::vertex_descriptor
vertex(typename boost::graph_traits<Intrinsic_Delaunay_Triangulation_3<TM,T,HCM,VPM,FIM,EIM,LA> >::halfedge_descriptor hd,
       const Intrinsic_Delaunay_Triangulation_3<TM,T,HCM,VPM,FIM,EIM,LA>& )
{
  return boost::graph_traits<Intrinsic_Delaunay_Triangulation_3<TM,T,HCM,VPM,FIM,EIM,LA> >::vertex_descriptor(hd);
}

  
template <typename TM,
          typename T,
          typename HCM,
          typename VPM,
          typename FIM,
          typename EIM,
          typename LA>
typename boost::graph_traits<Intrinsic_Delaunay_Triangulation_3<TM,T,HCM,VPM,FIM,EIM,LA> >::halfedge_descriptor
halfedge(typename boost::graph_traits<Intrinsic_Delaunay_Triangulation_3<TM,T,HCM,VPM,FIM,EIM,LA> >::face_descriptor fd,
       const Intrinsic_Delaunay_Triangulation_3<TM,T,HCM,VPM,FIM,EIM,LA>& idt)
{
  return halfedge(fd, idt.triangle_mesh());
}

  
template <typename TM,
          typename T,
          typename HCM,
          typename VPM,
          typename FIM,
          typename EIM,
          typename LA>
typename boost::graph_traits<Intrinsic_Delaunay_Triangulation_3<TM,T,HCM,VPM,FIM,EIM,LA> >::halfedge_descriptor
next(typename boost::graph_traits<Intrinsic_Delaunay_Triangulation_3<TM,T,HCM,VPM,FIM,EIM,LA> >::halfedge_descriptor hd,
       const Intrinsic_Delaunay_Triangulation_3<TM,T,HCM,VPM,FIM,EIM,LA>& idt)
{
  return next(hd, idt.triangle_mesh());
}
  
  
  template <typename IDT>
  struct IDT_vertex_point_property_map {
    const IDT& idt;
    typedef typename IDT::Triangle_mesh TM;
    typedef typename boost::graph_traits<IDT>::vertex_descriptor key_type;
    typedef typename IDT::Point_3 value_type;
    typedef value_type reference;
    
    IDT_vertex_point_property_map(const IDT& idt)
      : idt(idt)
      {}

    friend value_type get(const IDT_vertex_point_property_map<IDT>& pm,
                          key_type vd)
    {
      typename IDT::Point_2 p = get(pm.idt.hcmap(), vd.hd);
      return value_type(p.x(), p.y(), 0);
    }
  };

  
  template <typename IDT>
  struct IDT_vertex_index_property_map {
    const IDT& idt;
    typedef typename IDT::Triangle_mesh TM;
    typedef typename boost::graph_traits<IDT>::vertex_descriptor key_type;
    typedef typename boost::graph_traits<IDT>::vertices_size_type value_type;
    typedef value_type reference;
    
    IDT_vertex_index_property_map(const IDT& idt)
      : idt(idt)
      {}
    
    friend value_type get(const IDT_vertex_index_property_map<IDT>& pm,
                          key_type vd)
    {
      typename TM::vertex_descriptor tm_vd = target(vd.hd,pm.idt.triangle_mesh());
      
      return 0;// AF todo: return the index of tm_vd
    }
  };

} // namespace Intrinsic_Delaunay_Triangulation_3
} // namespace CGAL

namespace boost {
  
template <typename TM,
          typename T,
          typename HCM,
          typename VPM,
          typename FIM,
          typename EIM,
          typename LA>
 struct property_map<CGAL::Intrinsic_Delaunay_Triangulation_3::Intrinsic_Delaunay_Triangulation_3<TM,T,HCM,VPM,FIM,EIM,LA>,
                     CGAL::vertex_point_t > {
   typedef CGAL::Intrinsic_Delaunay_Triangulation_3::Intrinsic_Delaunay_Triangulation_3<TM,T,HCM,VPM,FIM,EIM,LA> IDT;
   typedef CGAL::Intrinsic_Delaunay_Triangulation_3::IDT_vertex_point_property_map<IDT> type;
   typedef type const_type;
 };

template <typename TM,
          typename T,
          typename HCM,
          typename VPM,
          typename FIM,
          typename EIM,
          typename LA>
CGAL::Intrinsic_Delaunay_Triangulation_3::IDT_vertex_point_property_map<CGAL::Intrinsic_Delaunay_Triangulation_3::Intrinsic_Delaunay_Triangulation_3<TM,T,HCM,VPM,FIM,EIM,LA> >
get(CGAL::vertex_point_t,
    const CGAL::Intrinsic_Delaunay_Triangulation_3::Intrinsic_Delaunay_Triangulation_3<TM,T,HCM,VPM,FIM,EIM,LA>& idt)
{
  return CGAL::Intrinsic_Delaunay_Triangulation_3::IDT_vertex_point_property_map<CGAL::Intrinsic_Delaunay_Triangulation_3::Intrinsic_Delaunay_Triangulation_3<TM,T,HCM,VPM,FIM,EIM,LA> >(idt);
}

  
template <typename TM,
          typename T,
          typename HCM,
          typename VPM,
          typename FIM,
          typename EIM,
          typename LA>
CGAL::Intrinsic_Delaunay_Triangulation_3::IDT_vertex_index_property_map<CGAL::Intrinsic_Delaunay_Triangulation_3::Intrinsic_Delaunay_Triangulation_3<TM,T,HCM,VPM,FIM,EIM,LA> >
get(boost::vertex_index_t,
    const CGAL::Intrinsic_Delaunay_Triangulation_3::Intrinsic_Delaunay_Triangulation_3<TM,T,HCM,VPM,FIM,EIM,LA>& idt)
{
  return CGAL::Intrinsic_Delaunay_Triangulation_3::IDT_vertex_index_property_map<CGAL::Intrinsic_Delaunay_Triangulation_3::Intrinsic_Delaunay_Triangulation_3<TM,T,HCM,VPM,FIM,EIM,LA> >(idt);
}

  
template <typename TM,
          typename T,
          typename HCM,
          typename VPM,
          typename FIM,
          typename EIM,
          typename LA,
          typename dT>
 struct property_map<CGAL::Intrinsic_Delaunay_Triangulation_3::Intrinsic_Delaunay_Triangulation_3<TM,T,HCM,VPM,FIM,EIM,LA>,
                     CGAL::dynamic_vertex_property_t<dT> > {
  typedef CGAL::Intrinsic_Delaunay_Triangulation_3::Intrinsic_Delaunay_Triangulation_3<TM,T,HCM,VPM,FIM,EIM,LA> IDT;
  typedef typename property_map<TM, CGAL::dynamic_vertex_property_t<dT> >::type type;
  typedef typename property_map<TM, CGAL::dynamic_vertex_property_t<dT> >::const_type const_type;
 };

} // namespace boost


namespace CGAL {
  
template <typename TM,
          typename T,
          typename HCM,
          typename VPM,
          typename FIM,
          typename EIM,
          typename LA,
          typename dT>
typename boost::property_map<TM, CGAL::dynamic_vertex_property_t<dT> >::const_type
get(CGAL::dynamic_vertex_property_t<dT> dvp,
    const CGAL::Intrinsic_Delaunay_Triangulation_3::Intrinsic_Delaunay_Triangulation_3<TM,T,HCM,VPM,FIM,EIM,LA>& idt)
{
  return get(dvp,idt.triangle_mesh());
}

template <typename TM,
          typename T,
          typename HCM,
          typename VPM,
          typename FIM,
          typename EIM,
          typename LA,
          typename dT>
typename boost::property_map<TM, CGAL::dynamic_vertex_property_t<dT> >::type
get(CGAL::dynamic_vertex_property_t<dT> dvp,
    CGAL::Intrinsic_Delaunay_Triangulation_3::Intrinsic_Delaunay_Triangulation_3<TM,T,HCM,VPM,FIM,EIM,LA>& idt)
{
  return get(dvp,idt.triangle_mesh());
}
  

} // namespace CGAL

#include <CGAL/enable_warnings.h>
#endif // CGAL_HEAT_METHOD_3_HEAT_METHOD_3_H
