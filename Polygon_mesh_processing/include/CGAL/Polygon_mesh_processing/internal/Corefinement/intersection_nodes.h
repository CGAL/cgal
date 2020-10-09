// Copyright (c) 2016 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sebastien Loriot

#ifndef CGAL_POLYGON_MESH_PROCESSING_INTERNAL_COREFINEMENT_INTERSECTION_NODES_H
#define CGAL_POLYGON_MESH_PROCESSING_INTERNAL_COREFINEMENT_INTERSECTION_NODES_H

#include <CGAL/license/Polygon_mesh_processing/corefinement.h>


#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <boost/type_traits/is_floating_point.hpp>
#include <CGAL/Kernel_traits.h>

namespace CGAL {
namespace Polygon_mesh_processing {
namespace Corefinement {

// A class responsible for storing the intersection nodes of the intersection
// polylines. Different specializations are available depending whether
// predicates on constructions are needed.
template <class TriangleMesh,
          class VertexPointMap1, class VertexPointMap2,
          bool Predicates_on_constructions_needed,
          bool Has_exact_constructions =
          !boost::is_floating_point<
            typename Kernel_traits<
              typename boost::property_traits<VertexPointMap1>::value_type
            >::Kernel::FT
           >::value >
class Intersection_nodes;

//Store only the double version of the intersection points.
template <class TriangleMesh,
          class VertexPointMap1, class VertexPointMap2>
class Intersection_nodes<TriangleMesh, VertexPointMap1, VertexPointMap2, false, false>
{
//typedefs
  typedef typename boost::property_traits<VertexPointMap1>::value_type  Point_3;
  CGAL_static_assertion((std::is_same<typename boost::property_traits<VertexPointMap1>::value_type,
                                      typename boost::property_traits<VertexPointMap2>::value_type>::value));

  typedef typename Kernel_traits<Point_3>::Kernel                  Input_kernel;
  typedef std::vector <Point_3>                                    Nodes_vector;
  typedef CGAL::Exact_predicates_exact_constructions_kernel        Exact_kernel;
  typedef CGAL::Cartesian_converter<Exact_kernel,Input_kernel>  Exact_to_double;
  typedef boost::graph_traits<TriangleMesh>                                  GT;
  typedef typename GT::halfedge_descriptor                  halfedge_descriptor;
  typedef typename GT::face_descriptor                          face_descriptor;
  typedef typename GT::vertex_descriptor                      vertex_descriptor;

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
  const VertexPointMap1& vpm1;
  const VertexPointMap2& vpm2;

  Intersection_nodes(const TriangleMesh& tm1_,
                     const TriangleMesh& tm2_,
                     const VertexPointMap1& vpm1_,
                     const VertexPointMap2& vpm2_)
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
  template <class VPM_A, class VPM_B> // VertexPointMap1 or VertexPointMap2
  void add_new_node(halfedge_descriptor h_a,
                    face_descriptor f_b,
                    const TriangleMesh& tm_a,
                    const TriangleMesh& tm_b,
                    const VPM_A& vpm_a,
                    const VPM_B& vpm_b)
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

  template <class VPM> // VertexPointMap1 or VertexPointMap2
  void call_put(const VPM& vpm, vertex_descriptor vd, std::size_t i, TriangleMesh&)
  {
    put(vpm, vd, nodes[i]);
  }

  void all_nodes_created(){}
  void finalize() {}

}; // end specialization
     // Intersection_nodes<Polyhedron,Kernel,No_predicates_on_constructions,false>

// second specializations: store an exact copy of the points so
// that we can answer exactly predicates
template <class TriangleMesh, class VertexPointMap1, class VertexPointMap2>
class Intersection_nodes<TriangleMesh, VertexPointMap1, VertexPointMap2, true, false>
{
//typedefs
public:
  typedef CGAL::Exact_predicates_exact_constructions_kernel        Exact_kernel;

private:
  typedef typename boost::property_traits<VertexPointMap1>::value_type  Point_3;
  CGAL_static_assertion((std::is_same<typename boost::property_traits<VertexPointMap1>::value_type,
                                      typename boost::property_traits<VertexPointMap2>::value_type>::value));

  typedef typename Kernel_traits<Point_3>::Kernel                  Input_kernel;

  typedef Cartesian_converter<Input_kernel,Exact_kernel>        Double_to_exact;
  typedef Cartesian_converter<Exact_kernel, Input_kernel>       Exact_to_double;

  typedef boost::graph_traits<TriangleMesh>                                  GT;
  typedef typename GT::halfedge_descriptor                  halfedge_descriptor;
  typedef typename GT::face_descriptor                          face_descriptor;
  typedef typename GT::vertex_descriptor                      vertex_descriptor;

  typedef std::vector <Exact_kernel::Point_3>                       Exact_nodes;
//members
  Exact_nodes enodes;

  Double_to_exact double_to_exact;
  Exact_to_double exact_to_double;
  Exact_kernel        ek;
  Exact_kernel::Intersect_3 exact_intersection;
  std::vector<vertex_descriptor> tm1_vertices, tm2_vertices;
  const bool doing_autorefinement;

public:
  const TriangleMesh &tm1, &tm2;
  const VertexPointMap1& vpm1;
  const VertexPointMap2& vpm2;

  Intersection_nodes(const TriangleMesh& tm1_,
                     const TriangleMesh& tm2_,
                     const VertexPointMap1& vpm1_,
                     const VertexPointMap2& vpm2_)
  : doing_autorefinement(&tm1_ == &tm2_)
  , tm1(tm1_)
  , tm2(tm2_)
  , vpm1(vpm1_)
  , vpm2(vpm2_)
  {}

  Point_3
  operator[](std::size_t i) const
  {
    return exact_to_double(enodes[i]);
  }

  const Exact_kernel::Point_3
  exact_node(std::size_t i) const
  {
    return enodes[i];
  }

  Exact_kernel::Point_3
  to_exact(const Point_3& p) const
  {
    return Exact_kernel::Point_3(p.x(), p.y(), p.z());
  }

  size_t size() const {return enodes.size();}

  void add_new_node(const Exact_kernel::Point_3& p)
  {
    const Exact_kernel::Approximate_kernel::Point_3& p_approx=p.approx();
    const double precision =
      Lazy_exact_nt<Exact_kernel::FT>::get_relative_precision_of_to_double();
    if ( !has_smaller_relative_precision(p_approx.x(),precision) ||
         !has_smaller_relative_precision(p_approx.y(),precision) ||
         !has_smaller_relative_precision(p_approx.z(),precision) )
    {
      p.exact();
    }
    enodes.push_back(p);
  }

  //add a new node in the final graph.
  //it is the intersection of the triangle with the segment
  template <class VPM_A, class VPM_B> // VertexPointMap1 or VertexPointMap2
  void add_new_node(halfedge_descriptor h_a,
                    face_descriptor f_b,
                    const TriangleMesh& tm_a,
                    const TriangleMesh& tm_b,
                    const VPM_A vpm_a,
                    const VPM_B vpm_b)
  {
    halfedge_descriptor h_b = halfedge(f_b, tm_b);
    add_new_node(
      Exact_kernel::Construct_plane_line_intersection_point_3()(
        to_exact( get(vpm_b, source(h_b,tm_b)) ),
        to_exact( get(vpm_b, target(h_b,tm_b)) ),
        to_exact( get(vpm_b, target(next(h_b,tm_b),tm_b)) ),
        to_exact( get(vpm_a, source(h_a,tm_a)) ),
        to_exact( get(vpm_a, target(h_a,tm_a)) ) ) );
  }

  // use to resolve intersection of 3 faces in autorefinement only
  template <class VPM>
  void add_new_node(halfedge_descriptor h1,
                    halfedge_descriptor h2,
                    halfedge_descriptor h3,
                    const TriangleMesh& tm,
                    const VPM& vpm)
  {
    // TODO Far from optimal!
    typedef Exact_kernel::Plane_3 Plane_3;
    Plane_3 p1(to_exact( get(vpm, source(h1,tm)) ),
               to_exact( get(vpm, target(h1,tm)) ),
               to_exact( get(vpm, target(next(h1,tm),tm)))),
            p2(to_exact( get(vpm, source(h2,tm)) ),
               to_exact( get(vpm, target(h2,tm)) ),
               to_exact( get(vpm, target(next(h2,tm),tm)))),
            p3(to_exact( get(vpm, source(h3,tm)) ),
               to_exact( get(vpm, target(h3,tm)) ),
               to_exact( get(vpm, target(next(h3,tm),tm))));
    typename cpp11::result_of<
      Exact_kernel::Intersect_3(Plane_3, Plane_3, Plane_3)
    >::type inter_res = exact_intersection(p1, p2, p3);

    CGAL_assertion(inter_res != boost::none);
    const Exact_kernel::Point_3* pt =
      boost::get<Exact_kernel::Point_3>(&(*inter_res));
    CGAL_assertion(pt!=nullptr);
    add_new_node(*pt);
  }

  //the point is an input
  void add_new_node(const Point_3& p){
    enodes.push_back(to_exact(p));
  }

  void all_nodes_created()
  {
    tm1_vertices.resize(enodes.size(), GT::null_vertex());
    tm2_vertices.resize(enodes.size(), GT::null_vertex());
  }

  template <class VPM> // VertexPointMap1 or VertexPointMap2
  void call_put(const VPM& vpm, vertex_descriptor vd, std::size_t i, TriangleMesh& tm)
  {
    put(vpm, vd, exact_to_double(enodes[i]));
    if (&tm1==&tm)
    {
      if (  tm1_vertices[i] == GT::null_vertex() )
      {
        tm1_vertices[i] = vd;
        return;
      }
      if (doing_autorefinement)
        tm2_vertices[i] = vd;
    }
    else
      tm2_vertices[i] = vd;
  }

  void finalize()
  {
    for (std::size_t i=0, e=enodes.size(); i!=e; ++i)
    {
      Point_3 pt = exact_to_double(enodes[i]);
      if ( tm1_vertices[i] != GT::null_vertex() )
        put(vpm1, tm1_vertices[i], pt);
      if ( tm2_vertices[i] != GT::null_vertex() )
        put(vpm2, tm2_vertices[i], pt);
    }
  }
}; // end specialization
     // Intersection_nodes<Polyhedron,Kernel,Predicates_on_constructions,false>


//Third specialization: The kernel already has exact constructions.
template <class TriangleMesh, class VertexPointMap1, class VertexPointMap2,
          bool Predicates_on_constructions_needed>
class Intersection_nodes<TriangleMesh, VertexPointMap1, VertexPointMap2,
                         Predicates_on_constructions_needed, true>
{
//typedefs
  typedef typename boost::property_traits<VertexPointMap1>::value_type  Point_3;
  CGAL_static_assertion((std::is_same<typename boost::property_traits<VertexPointMap1>::value_type,
                                      typename boost::property_traits<VertexPointMap2>::value_type>::value));

  typedef typename Kernel_traits<Point_3>::Kernel                  Input_kernel;
  typedef std::vector <Point_3>                                    Nodes_vector;

  typedef boost::graph_traits<TriangleMesh>                                  GT;
  typedef typename GT::halfedge_descriptor                  halfedge_descriptor;
  typedef typename GT::face_descriptor                          face_descriptor;
  typedef typename GT::vertex_descriptor                      vertex_descriptor;
//members
  Nodes_vector nodes;
  Input_kernel k;
  typename Input_kernel::Intersect_3 intersection;
public:
  typedef Input_kernel                                             Exact_kernel;

  const TriangleMesh &tm1, &tm2;
  const VertexPointMap1& vpm1;
  const VertexPointMap2& vpm2;

  Intersection_nodes(const TriangleMesh& tm1_,
                     const TriangleMesh& tm2_,
                     const VertexPointMap1& vpm1_,
                     const VertexPointMap2& vpm2_)
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

  template <class VPM>
  void add_new_node(halfedge_descriptor h1,
                    halfedge_descriptor h2,
                    halfedge_descriptor h3,
                    const TriangleMesh& tm,
                    const VPM& vpm)
  {
    // TODO Far from optimal!
    typedef typename Exact_kernel::Plane_3 Plane_3;
    Plane_3 p1( get(vpm, source(h1,tm)),
                get(vpm, target(h1,tm)),
                get(vpm, target(next(h1,tm),tm))),
            p2(get(vpm, source(h2,tm)),
               get(vpm, target(h2,tm)),
               get(vpm, target(next(h2,tm),tm))),
            p3(get(vpm, source(h3,tm)),
               get(vpm, target(h3,tm)),
               get(vpm, target(next(h3,tm),tm)));
    typename cpp11::result_of<
      typename Exact_kernel::Intersect_3(Plane_3, Plane_3, Plane_3)
    >::type inter_res = intersection(p1, p2, p3);

    CGAL_assertion(inter_res != boost::none);
    const Point_3* pt =
      boost::get<Point_3>(&(*inter_res));
    CGAL_assertion(pt!=nullptr);
    add_new_node(*pt);
  }

  //add a new node in the final graph.
  //it is the intersection of the triangle with the segment
  template <class VPM_A, class VPM_B> // VertexPointMap1 or VertexPointMap2
  void add_new_node(halfedge_descriptor h_a,
                    face_descriptor f_b,
                    const TriangleMesh& tm_a,
                    const TriangleMesh& tm_b,
                    const VPM_A& vpm_a,
                    const VPM_B& vpm_b)
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

  void add_new_node(const Point_3& p)
  {
    nodes.push_back(p);
  }

  const Point_3& to_exact(const Point_3& p) const { return p; }

  template <class VPM> // VertexPointMap1 or VertexPointMap2
  void call_put(const VPM& vpm, vertex_descriptor vd, std::size_t i, TriangleMesh&)
  {
    put(vpm, vd, nodes[i]);
  }

  void all_nodes_created(){}
  void finalize() {}


}; // end specialization


} } } // CGAL::Polygon_mesh_processing::Corefinement

#endif // CGAL_POLYGON_MESH_PROCESSING_INTERNAL_COREFINEMENT_INTERSECTION_NODES_H
