// Copyright (c) 2016 GeometryFactory (France).
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
//
//
// Author(s)     : Sebastien Loriot

#ifndef CGAL_POLYGON_MESH_PROCESSING_INTERNAL_COREFINEMENT_INTERSECTION_NODES_H
#define CGAL_POLYGON_MESH_PROCESSING_INTERNAL_COREFINEMENT_INTERSECTION_NODES_H

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <boost/type_traits/is_floating_point.hpp>
#include <CGAL/Kernel_traits.h>

namespace CGAL {
namespace Corefinement {

// A class responsible for storing the intersection nodes of the intersection
// polylines. Different specializations are available depending whether
// predicates on constructions are needed.
template <class TriangleMesh,
          class VertexPointMap,
          bool Predicates_on_constructions_needed,
          bool Has_exact_constructions=
          !boost::is_floating_point<
            typename Kernel_traits<
              typename boost::property_traits<VertexPointMap>::value_type
            >::Kernel::FT
           >::value >
class Intersection_nodes;

//Store only the double version of the intersection points.
template <class TriangleMesh,
          class VertexPointMap>
class Intersection_nodes<TriangleMesh,VertexPointMap,false,false>
{
//typedefs
  typedef typename boost::property_traits<VertexPointMap>::value_type   Point_3;
  typedef typename Kernel_traits<Point_3>::Kernel                  Input_kernel;
  typedef std::vector <Point_3>                                    Nodes_vector;
  typedef CGAL::Exact_predicates_exact_constructions_kernel        Exact_kernel;
  typedef CGAL::Cartesian_converter<Exact_kernel,Input_kernel>  Exact_to_double;
  typedef boost::graph_traits<TriangleMesh>                                  GT;
  typedef typename GT::halfedge_descriptor                  halfedge_descriptor;
  typedef typename GT::face_descriptor                          face_descriptor;

//members
  Nodes_vector nodes;
  Exact_kernel ek;
  Exact_to_double exact_to_double;

  typename Exact_kernel::Point_3
  to_exact(const typename Input_kernel::Point_3& p) const
  {
    return typename Exact_kernel::Point_3(p.x(), p.y(), p.z());
  }

public:
  const TriangleMesh &tm1, &tm2;
  VertexPointMap vpm1, vpm2;

  Intersection_nodes(const TriangleMesh& tm1_,
                     const TriangleMesh& tm2_,
                     const VertexPointMap& vpm1_,
                     const VertexPointMap& vpm2_)
  : tm1(tm1_)
  , tm2(tm2_)
  , vpm1(vpm1_)
  , vpm2(vpm2_)
  {}

  const Point_3& operator[](std::size_t i) const {
    return nodes[i];
  }

  size_t size() const {return nodes.size();}

  void add_new_node(const Point_3& p)
  {
    nodes.push_back(p);
  }

  void add_new_node(const typename Exact_kernel::Point_3& p)
  {
    nodes.push_back(  exact_to_double(p) );
  }

  //add a new node in the final graph.
  //it is the intersection of the triangle with the segment
  void add_new_node(halfedge_descriptor h_a,
                    face_descriptor f_b,
                    const TriangleMesh& tm_a,
                    const TriangleMesh& tm_b,
                    const VertexPointMap vpm_a,
                    const VertexPointMap& vpm_b)
  {
    halfedge_descriptor h_b = halfedge(f_b, tm_b);
    add_new_node(
      typename Exact_kernel::Construct_plane_line_intersection_point_3()(
        to_exact( get(vpm_b, source(h_b,tm_b)) ),
        to_exact( get(vpm_b, target(h_b,tm_b)) ),
        to_exact( get(vpm_b, target(next(h_b,tm_b),tm_b)) ),
        to_exact( get(vpm_a, source(h_a,tm_a)) ),
        to_exact( get(vpm_a, target(h_a,tm_a)) ) ) );
  }

  void add_new_node(halfedge_descriptor edge_1, face_descriptor face_2)
  {
    add_new_node(edge_1, face_2, tm1, tm2, vpm1, vpm2);
  }

}; // end specialization
     // Intersection_nodes<Polyhedron,Kernel,No_predicates_on_constructions,false>

//second specializations: store an exact copy of the points so that we can answer exactly predicates
//FYI, it used to have two specializations (one in the case the polyhedron
//can be edited and on if it cannot) building exact representation on demand.
//In the former case, we were using facet and halfedge while in the latter
//triple of vertex_handle and pair of vertex_handle
template <class TriangleMesh, class VertexPointMap>
class Intersection_nodes<TriangleMesh,VertexPointMap,true,false>
{
//typedefs
public:
  typedef CGAL::Simple_cartesian<CGAL::Interval_nt<false> >     Interval_kernel;
  typedef CGAL::Exact_predicates_exact_constructions_kernel        Exact_kernel;
private:
  typedef typename boost::property_traits<VertexPointMap>::value_type   Point_3;
  typedef typename Kernel_traits<Point_3>::Kernel                  Input_kernel;

  typedef Cartesian_converter<Interval_kernel,Input_kernel>  Interval_to_double;
  typedef Cartesian_converter<Input_kernel,Interval_kernel>  Double_to_interval;
  typedef Cartesian_converter<Exact_kernel,Interval_kernel>   Exact_to_interval;
  typedef Cartesian_converter<Input_kernel,Exact_kernel>        Double_to_exact;

  typedef boost::graph_traits<TriangleMesh>                                  GT;
  typedef typename GT::halfedge_descriptor                  halfedge_descriptor;
  typedef typename GT::face_descriptor                          face_descriptor;

  typedef std::vector <Interval_kernel::Point_3>                 Interval_nodes;
  typedef std::vector <Exact_kernel::Point_3>                       Exact_nodes;
//members
  Interval_nodes  inodes;
  Exact_nodes enodes;

  Interval_to_double  interval_to_double;
  Exact_to_interval   exact_to_interval;
  Double_to_interval  double_to_interval;
  Double_to_exact double_to_exact;
  Exact_kernel        ek;

public:
  const TriangleMesh &tm1, &tm2;
  VertexPointMap vpm1, vpm2;
  typedef CGAL::Interval_nt<false>::Protector                         Protector;

  Intersection_nodes(const TriangleMesh& tm1_,
                     const TriangleMesh& tm2_,
                     const VertexPointMap& vpm1_,
                     const VertexPointMap& vpm2_)
  : tm1(tm1_)
  , tm2(tm2_)
  , vpm1(vpm1_)
  , vpm2(vpm2_)
  {}


  Point_3
  operator[](std::size_t i) const
  {
    return interval_to_double(inodes[i]);
  }

  const Interval_kernel::Point_3&
  interval_node(std::size_t i) const
  {
    return inodes[i];
  }

  Interval_kernel::Point_3
  to_interval(const Point_3& p) const
  {
    return double_to_interval(p);
  }

  const Exact_kernel::Point_3
  exact_node(std::size_t i) const
  {
    return enodes[i];
  }

  Exact_kernel::Point_3
  to_exact(const Point_3& p) const
  {
    return typename Exact_kernel::Point_3(p.x(), p.y(), p.z());
  }

  size_t size() const {return enodes.size();}

  void add_new_node(const Exact_kernel::Point_3& p)
  {
    const Interval_kernel::Point_3& p_approx=p.approx();
    const double precision =
      Lazy_exact_nt<typename Exact_kernel::FT>::get_relative_precision_of_to_double();
    if ( !has_smaller_relative_precision(p_approx.x(),precision) ||
         !has_smaller_relative_precision(p_approx.y(),precision) ||
         !has_smaller_relative_precision(p_approx.z(),precision) )
    {
      p.exact();
    }
    enodes.push_back(p);
    inodes.push_back( exact_to_interval(p) );
  }

  //add a new node in the final graph.
  //it is the intersection of the triangle with the segment
  void add_new_node(halfedge_descriptor h_a,
                    face_descriptor f_b,
                    const TriangleMesh& tm_a,
                    const TriangleMesh& tm_b,
                    const VertexPointMap vpm_a,
                    const VertexPointMap& vpm_b)
  {
    halfedge_descriptor h_b = halfedge(f_b, tm_b);
    add_new_node(
      typename Exact_kernel::Construct_plane_line_intersection_point_3()(
        to_exact( get(vpm_b, source(h_b,tm_b)) ),
        to_exact( get(vpm_b, target(h_b,tm_b)) ),
        to_exact( get(vpm_b, target(next(h_b,tm_b),tm_b)) ),
        to_exact( get(vpm_a, source(h_a,tm_a)) ),
        to_exact( get(vpm_a, target(h_a,tm_a)) ) ) );
  }

  void add_new_node(halfedge_descriptor edge_1, face_descriptor face_2)
  {
    add_new_node(edge_1, face_2, tm1, tm2, vpm1, vpm2);
  }

  //the point is an input
  void add_new_node(const Point_3& p){
    enodes.push_back(to_exact(p));
    inodes.push_back( double_to_interval(p) );
  }
}; // end specialization
     // Intersection_nodes<Polyhedron,Kernel,Predicates_on_constructions,false>


//Third specialization: The kernel already has exact constructions.
template <class TriangleMesh,class VertexPointMap,bool Predicates_on_constructions_needed>
class Intersection_nodes<TriangleMesh,VertexPointMap,Predicates_on_constructions_needed,true>
{
//typedefs
  typedef typename boost::property_traits<VertexPointMap>::value_type   Point_3;
  typedef typename Kernel_traits<Point_3>::Kernel                  Input_kernel;
  typedef std::vector <Point_3>                                    Nodes_vector;

  typedef boost::graph_traits<TriangleMesh>                                  GT;
  typedef typename GT::halfedge_descriptor                  halfedge_descriptor;
  typedef typename GT::face_descriptor                          face_descriptor;
//members
  Nodes_vector nodes;
  Input_kernel k;
public:
  typedef Input_kernel                                          Interval_kernel;
  typedef Input_kernel                                             Exact_kernel;
  typedef void*                                                       Protector;

  const TriangleMesh &tm1, &tm2;
  VertexPointMap vpm1, vpm2;

  Intersection_nodes(const TriangleMesh& tm1_,
                     const TriangleMesh& tm2_,
                     const VertexPointMap& vpm1_,
                     const VertexPointMap& vpm2_)
  : tm1(tm1_)
  , tm2(tm2_)
  , vpm1(vpm1_)
  , vpm2(vpm2_)
  {}

  const Point_3& operator[](std::size_t i) const
  {
    return nodes[i];
  }

  size_t size() const {return nodes.size();}
  const Point_3& exact_node(std::size_t i) const {return nodes[i];}
  const Point_3& interval_node(std::size_t i) const {return nodes[i];}

  //add a new node in the final graph.
  //it is the intersection of the triangle with the segment
  void add_new_node(halfedge_descriptor h_a,
                    face_descriptor f_b,
                    const TriangleMesh& tm_a,
                    const TriangleMesh& tm_b,
                    const VertexPointMap vpm_a,
                    const VertexPointMap& vpm_b)
  {
    halfedge_descriptor h_b=halfedge(f_b,tm_b);

    add_new_node(
      typename Exact_kernel::Construct_plane_line_intersection_point_3()(
        get(vpm_b, source(h_b,tm_b)),
        get(vpm_b, target(h_b,tm_b)),
        get(vpm_b, target(next(h_b,tm_b),tm_b)),
        get(vpm_a, source(h_a,tm_a)),
        get(vpm_a, target(h_a,tm_a)) ) );
  }

  void add_new_node(halfedge_descriptor edge_1, face_descriptor face_2)
  {
    add_new_node(edge_1, face_2, tm1, tm2, vpm1, vpm2);
  }


  void add_new_node(const Point_3& p)
  {
    nodes.push_back(p);
  }

  const Point_3& to_interval(const Point_3& p) const { return p; }
  const Point_3& to_exact(const Point_3& p) const { return p; }

}; // end specialization


} } // end of namespace CGAL::Corefinement

#endif // CGAL_POLYGON_MESH_PROCESSING_INTERNAL_COREFINEMENT_INTERSECTION_NODES_H
