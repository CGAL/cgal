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

template <class TriangleMesh,
          class VertexPointMap,
          class Exact_kernel>
typename Exact_kernel::Point_3
compute_triangle_segment_intersection_point(
  typename boost::graph_traits<TriangleMesh>::halfedge_descriptor hd,
  typename boost::graph_traits<TriangleMesh>::face_descriptor fd,
  const Exact_kernel& ek,
  const TriangleMesh& tm_h,
  const TriangleMesh& tm_f,
  const VertexPointMap& vpm_h,
  const VertexPointMap& vpm_f)
{
  typedef typename Kernel_traits<
   typename boost::property_traits<VertexPointMap>::value_type
  >::Kernel Input_kernel;
  typedef boost::graph_traits<TriangleMesh> GT;
  typedef typename GT::vertex_descriptor vertex_descriptor;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;
  typedef typename Exact_kernel::Line_3 Line_3;
  typedef typename Exact_kernel::Plane_3 Plane_3;

  CGAL::Cartesian_converter<Input_kernel,Exact_kernel> to_exact;
  halfedge_descriptor hdf=halfedge(fd,tm_f);
  vertex_descriptor vf1=source(hdf,tm_f),
                    vf2=target(hdf,tm_f),
                    vf3=target(next(hdf,tm_f),tm_f),
                    vh1=source(hd,tm_h),
                    vh2=target(hd,tm_h);

  Plane_3 plane(to_exact( get(vpm_f, vf1) ),
                to_exact( get(vpm_f, vf2) ),
                to_exact( get(vpm_f, vf3) ) );

  Line_3 line( to_exact( get(vpm_h, vh1) ), to_exact( get(vpm_h, vh2) ) );

  typename cpp11::result_of<typename Exact_kernel::Intersect_3(Plane_3,Line_3)>::type
    res = ek.intersect_3_object()(plane,line);
  CGAL_assertion(res!=boost::none);
  const typename Exact_kernel::Point_3* e_pt =
    boost::get<typename Exact_kernel::Point_3>(&(*res));
  CGAL_assertion(e_pt!=NULL);
  return *e_pt;
}

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
  // typedef CGAL::Interval_nt<true>::Protector                          Protector;
  typedef boost::graph_traits<TriangleMesh>                                  GT;
  typedef typename GT::halfedge_descriptor                  halfedge_descriptor;
  typedef typename GT::face_descriptor                          face_descriptor;

//members
  Nodes_vector nodes;
  Exact_kernel ek;
  Exact_to_double exact_to_double;
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

  const Point_3& operator[](int i) const {
    return nodes[i];
  }

  // const Point_3& exact_node(int i) const {return nodes[i];}
  // const Point_3& interval_node(int i) const {return nodes[i];}
  // const Point_3& to_exact(const typename Kernel::Point_3& p) const {return p;}
  // const Point_3& to_interval(const typename Kernel::Point_3& p) const {return p;}

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
  void add_new_node(halfedge_descriptor edge_1, face_descriptor face_2)
  {
    add_new_node(
      compute_triangle_segment_intersection_point(edge_1, face_2, ek, tm1, tm2, vpm1, vpm2)
    );
  }

  void add_new_node(halfedge_descriptor h_a,
                    face_descriptor f_b,
                    const TriangleMesh& tm_a,
                    const TriangleMesh& tm_b,
                    const VertexPointMap vpm_a,
                    const VertexPointMap& vpm_b)
  {
      add_new_node(
        compute_triangle_segment_intersection_point(h_a, f_b, ek, tm_a, tm_b, vpm_a, vpm_b)
      );
  }

}; // end specialization
     // Intersection_nodes<Polyhedron,Kernel,No_predicates_on_constructions,false>

#if 0
//second specializations: store an exact copy of the points so that we can answer exactly predicates
//FYI, it used to have two specializations (one in the case the polyhedron
//can be edited and on if it cannot) building exact representation on demand.
//In the former case, we were using facet and halfedge while in the latter
//triple of vertex_handle and pair of vertex_handle
template <class Polyhedron, class PolyhedronPointPMap, class Kernel>
class Intersection_nodes<Polyhedron,PolyhedronPointPMap,Kernel,Predicates_on_constructions,false>
{
//typedefs
public:
  typedef CGAL::Simple_cartesian<CGAL::Interval_nt<false> >  Ikernel;
  typedef CGAL::Exact_predicates_exact_constructions_kernel  Exact_kernel;
private:
  typedef CGAL::Cartesian_converter<Ikernel,Kernel>          Interval_to_double;
  typedef CGAL::Cartesian_converter<Kernel,Ikernel>          Double_to_interval;
  typedef CGAL::Cartesian_converter<Exact_kernel,Ikernel>    Exact_to_interval;
  typedef CGAL::Cartesian_converter<Kernel,Exact_kernel>     Double_to_exact;

  typedef typename Polyhedron::Vertex_const_handle           Vertex_handle;
  typedef typename Polyhedron::Halfedge_const_handle         Halfedge_handle;
  typedef typename Polyhedron::Facet_const_handle            Facet_handle;

  typedef std::vector <Ikernel::Point_3>      Interval_nodes;
  typedef std::vector <Exact_kernel::Point_3> Exact_nodes;


//members
  Interval_nodes  inodes;
  Exact_nodes enodes;

  Interval_to_double  interval_to_double;
  Exact_to_interval   exact_to_interval;
  Double_to_interval  double_to_interval;
  Double_to_exact double_to_exact;
  Exact_kernel        ek;
  PolyhedronPointPMap ppmap;

public:

  Intersection_nodes(PolyhedronPointPMap ppmap):
    ppmap(ppmap){}

  typedef CGAL::Interval_nt<false>::Protector                 Protector;

  typename Kernel::Point_3
  operator[](int i) const {
    return interval_to_double(inodes[i]);
  }

  const typename Ikernel::Point_3&
  interval_node(int i) const {
    return inodes[i];
  }

  typename Ikernel::Point_3
  to_interval(const typename Kernel::Point_3& p) const {
    return double_to_interval(p);
  }

  const Exact_kernel::Point_3
  exact_node(int i) const {
    return enodes[i];
  }

  typename Exact_kernel::Point_3
  to_exact(const typename Kernel::Point_3& p) const {
    return double_to_exact(p);
  }


  size_t size() const {return enodes.size();}

  void add_new_node(const Exact_kernel::Point_3& p){
    const Ikernel::Point_3& p_approx=p.approx();
    if ( !has_smaller_relative_precision(p_approx.x(),Lazy_exact_nt<typename Exact_kernel::FT>::get_relative_precision_of_to_double()) ||
         !has_smaller_relative_precision(p_approx.y(),Lazy_exact_nt<typename Exact_kernel::FT>::get_relative_precision_of_to_double()) ||
         !has_smaller_relative_precision(p_approx.z(),Lazy_exact_nt<typename Exact_kernel::FT>::get_relative_precision_of_to_double()) ) p.exact();
    enodes.push_back(p);
    inodes.push_back( exact_to_interval(p) );
  }

  void add_new_node(Halfedge_handle edge,Facet_handle facet)
  {
    add_new_node( compute_triangle_segment_intersection_point<Polyhedron,Kernel>(edge,facet,ek,ppmap) );
  }

  //the point is an input
  void add_new_node(const typename Kernel::Point_3& p){
    enodes.push_back(to_exact(p));
    inodes.push_back( double_to_interval(p) );
  }
}; // end specialization
     // Intersection_nodes<Polyhedron,Kernel,Predicates_on_constructions,false>

//Third specialization: The kernel already has exact constructions.
template <class Polyhedron,class PolyhedronPointPMap,class Kernel,class Node_storage>
class Intersection_nodes<Polyhedron,PolyhedronPointPMap,Kernel,Node_storage,true>
{
//typedefs
  typedef std::vector <typename Kernel::Point_3>             Nodes_vector;
  typedef typename Polyhedron::Halfedge_const_handle         Halfedge_handle;
  typedef typename Polyhedron::Facet_const_handle            Facet_handle;
//members
  Nodes_vector nodes;
  Kernel k;
  PolyhedronPointPMap ppmap;
public:
  typedef Kernel Ikernel;
  typedef Kernel Exact_kernel;
  typedef void* Protector;

  Intersection_nodes(PolyhedronPointPMap ppmap):
    ppmap(ppmap){}

  const typename Kernel::Point_3&
  operator[](int i) const {
    return nodes[i];
  }

  size_t size() const {return nodes.size();}
  const typename Kernel::Point_3& exact_node(int i) const {return nodes[i];}
  const typename Kernel::Point_3& interval_node(int i) const {return nodes[i];}

  //add a new node in the final graph.
  //it is the intersection of the triangle with the segment
  void add_new_node(Halfedge_handle edge,Facet_handle facet)
  {
    nodes.push_back (
      compute_triangle_segment_intersection_point<Polyhedron,Kernel>(edge,facet,k,ppmap)
    );
  }

  void add_new_node(const typename Kernel::Point_3& p)
  {
    nodes.push_back(p);
  }

  const typename Kernel::Point_3& to_interval(const typename Kernel::Point_3& p) const { return p; }
  const typename Kernel::Point_3& to_exact(const typename Kernel::Point_3& p) const { return p; }

}; // end specialization
#endif

} } // end of namespace CGAL::Corefinement

#endif // CGAL_POLYGON_MESH_PROCESSING_INTERNAL_COREFINEMENT_INTERSECTION_NODES_H
