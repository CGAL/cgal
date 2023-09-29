// Copyright (c) 2014-2019  GeometryFactory (France). All rights reserved.
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
// Author(s)     : Maxime Gimeno
//

#ifndef CGAL_POLYGON_MESH_PROCESSING_WALK_TO_SELECT_H
#define CGAL_POLYGON_MESH_PROCESSING_WALK_TO_SELECT_H

// #include <CGAL/license/Polygon_mesh_processing/repair.h>

#include <CGAL/Polyhedron_simplex_type.h>
#include <CGAL/property_map.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/walk_in_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/intersections.h>
#include <CGAL/boost/graph/property_maps.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/boost/graph/dijkstra_shortest_paths.h>

#include <boost/optional.hpp>

#include <bitset>

#include <unordered_map>

//todo: if no point is between the planes, nothing is selected.
namespace CGAL {

namespace Polygon_mesh_processing {

namespace walker_internal{
template <class PolygonMesh, class VertexPointMap>
class Edge_intersection_with_plane
{
  typedef typename Kernel_traits<
    typename boost::property_traits<VertexPointMap>::value_type
  >::Kernel K;

  typedef typename K::Point_3 Point_3;
  typedef typename K::Vector_3 Vector_3;
  typedef typename K::FT FT;

  typedef boost::variant<bool, FT> Variant;
  typedef boost::optional< std::pair<FT, Variant> > result_type;

  Point_3 m_point;
  Vector_3 m_normal;
  const PolygonMesh* m_pm_ptr;
  VertexPointMap m_vpm;
  Vector_3 m_other_normal;
  mutable FT res_scalar;
  mutable int step;

  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor
  halfedge_descriptor;

public:

  typedef FT Barycentric_NT;

  Edge_intersection_with_plane(const Point_3& pt,
                               const Vector_3& n,
                               const Vector_3& other_n,
                               const PolygonMesh& pm,
                               const VertexPointMap& vpm)
    : m_point(pt)
    , m_normal(n)
    ,  m_pm_ptr(&pm)
    ,  m_vpm(vpm)
    , m_other_normal(other_n)
    , res_scalar(1)
    , step(0)
  {}

  result_type operator()(halfedge_descriptor h) const
  {
    typename K::Oriented_side_3 oriented_side;
    Point_3 b0 = get(m_vpm, source(h, *m_pm_ptr) );
    Point_3 b1 = get(m_vpm, target(h, *m_pm_ptr) );

    Oriented_side s0 = oriented_side(m_point, m_normal, b0);
    Oriented_side s1 = oriented_side(m_point, m_normal, b1);

    //the path is jagged, so we take the biggest projection of 2 steps to
    //get the walking direction
    if(step == 0){
      res_scalar = scalar_prod(h);
    }
    if(step == 1){
      FT new_dir= scalar_prod(h);
      if(std::abs(new_dir) > std::abs(res_scalar))
        res_scalar = new_dir;

    }
    if (s0 == ON_ORIENTED_BOUNDARY){return std::make_pair(0.5, Variant(false));}
    if (s1 == ON_ORIENTED_BOUNDARY){return std::make_pair(0.5, Variant(true)); }
    if(s0==s1) return boost::none;
    FT num = ( (m_point - b1) * m_normal );
    FT denum = ((b0-b1) * m_normal);

    step++;
    // the edge is almost in the plane. This is an arbitrary choice
    if (denum==0)
      return std::make_pair(0.5, Variant( 0.5 ));

    FT alpha =  num / denum;

    if (alpha<=0) std::make_pair(0.5, Variant(false));
    if (alpha>=1) std::make_pair(0.5, Variant(true));
    return std::make_pair(0.5, Variant( alpha ));
  }

  FT scalar_prod(const halfedge_descriptor& input) const
  {
    Point_3 b0 = get(m_vpm, source(input, *m_pm_ptr) );
    Point_3 b1 = get(m_vpm, target(input, *m_pm_ptr) );
    Vector_3 dir(b0, b1);
    return dir*m_other_normal;
  }

  bool wrong_direction() const
  {
    return res_scalar <= 0;
  }

  void reset()
  {
    res_scalar = 1;
    step = 0;
  }

  bool do_break() const
  {
    return step > 1 && wrong_direction();
  }
};

template <class PolygonMesh, class Kernel>
class Intersection_info
{
  // typedefs
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Point_3 Point_3;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;

  // data members
  bool m_is_vertex;
  halfedge_descriptor m_hedge;
  Point_3 m_point;

  // internal function
  template <class VertexPointMap>
  Point_3
  compute_point(const FT& alpha, const PolygonMesh& pm, VertexPointMap vpm)
  {
    CGAL_assertion(!m_is_vertex);
    CGAL_assertion(alpha<1 && alpha >0);
    // get canonical_edge
    halfedge_descriptor hedge = m_hedge < opposite(m_hedge, pm) ?
          m_hedge : opposite(m_hedge, pm);

    const Point_3 b0 = get(vpm, source(hedge, pm));
    const Point_3 b1 = get(vpm, target(hedge, pm));

    return b1+alpha*(b0-b1);
  }

public:

  Intersection_info()
    : m_is_vertex(true)
  {}

  template <class VertexPointMap>
  Intersection_info(      halfedge_descriptor h,
                          const PolygonMesh& pm,
                          VertexPointMap vpm)
    : m_is_vertex(true)
    , m_hedge(h)
    , m_point(get(vpm, target(h, pm)))
  {}

  template <class VertexPointMap>
  Intersection_info(      halfedge_descriptor h,
                          const FT& alpha,
                          const PolygonMesh& pm,
                          VertexPointMap vpm)
    : m_is_vertex(false)
    , m_hedge(h)
    , m_point(compute_point(alpha, pm, vpm))
  {}


  Intersection_info(      halfedge_descriptor h,
                    const Point_3& p)
    : m_is_vertex(false)
    , m_hedge(h)
    , m_point(p)
  {}

  const Point_3& point() const
  {
    return m_point;
  }

  bool is_vertex() const
  {
    return m_is_vertex;
  }

  halfedge_descriptor halfedge() const
  {
    return m_hedge;
  }

};

template <class TriangleMesh, class VertexPointMap>
class Walk_in_polygon_mesh_visitor {
  typedef boost::graph_traits<TriangleMesh> BGT;
  typedef typename BGT::vertex_descriptor vertex_descriptor;
  typedef typename BGT::halfedge_descriptor halfedge_descriptor;
  typedef typename BGT::face_descriptor face_descriptor;

  typedef typename boost::property_traits<VertexPointMap>::value_type Point_3;
  typedef typename Kernel_traits<Point_3>::Kernel K;
  typedef typename K::FT FT;

  typedef Intersection_info<TriangleMesh, K> Inter_info;
  // to collect the target vertices of the polyline during the walk (but the last one)
  std::vector<Inter_info> target_vertices;
  // the set of halfedges showing the walk
  std::vector< halfedge_descriptor > m_walked_halfedges;
  // the set of faces modified
  std::vector< face_descriptor > m_split_faces;
  // the set of faces created
  std::vector< face_descriptor > m_new_faces;
  // the face where the walk ended (if not on an edge)
  face_descriptor m_last_face;
  // The triangle mesh to be refined
  const TriangleMesh* m_tm_ptr;
  //The point property map
  VertexPointMap m_vpm;
  //the distance to reach whiPoint_set5le walking
  FT m_dist;
  //the distance currently reached
  FT cur_dist;
  //the last value of cur_dist (used to find the end_point)
  FT last_dist;
  //the initial center
  Point_3 text_center;

public:

  Walk_in_polygon_mesh_visitor(const TriangleMesh& tm,
                               VertexPointMap vpm,
                               FT dist,
                               const Point_3& center)
    : m_tm_ptr(&tm)
    , m_vpm(vpm)
    , m_dist(dist)
    , cur_dist(0)
    , last_dist(0)
    , text_center(center)
  {}

  // intersected_edge points in the face into which the edge should be added
  void found_edge_intersection(halfedge_descriptor intersected_edge,
                               const FT& barycentric_coord )
  {
    m_walked_halfedges.push_back( intersected_edge );
    target_vertices.push_back( Inter_info(intersected_edge, barycentric_coord, *m_tm_ptr, m_vpm) );
    Point_3 new_p = target_vertices.back().point();

    //update distance
    if(target_vertices.size() == 1)
    {
      cur_dist += CGAL::sqrt(CGAL::squared_distance(text_center,new_p));
    }
    else {
      Point_3 last_point = target_vertices[target_vertices.size() -2].point();
      last_dist = cur_dist;
      cur_dist += CGAL::sqrt(CGAL::squared_distance(last_point,new_p));
    }
  }

  void found_vertex_intersection(halfedge_descriptor h)
  {
    m_walked_halfedges.push_back( h );
    target_vertices.push_back( Inter_info(h, *m_tm_ptr, m_vpm) );
    Point_3 new_p = target_vertices.back().point();
    //update distance
    if(target_vertices.size() == 1)
    {
      cur_dist += CGAL::sqrt(CGAL::squared_distance(text_center,new_p));
    }
    else {
      Point_3 last_point = target_vertices[target_vertices.size() -2].point();
      last_dist = cur_dist;
      cur_dist += CGAL::sqrt(CGAL::squared_distance(last_point,new_p));
    }
  }

  void register_initial_intersection(const Inter_info& i)
  {
    target_vertices.push_back(i);
    Point_3 new_p = target_vertices.back().point();
    //update distance
    CGAL_assertion(target_vertices.size() == 1);
    cur_dist = CGAL::sqrt(CGAL::squared_distance(text_center,new_p));
  }

  void on_walk_end(face_descriptor f)
  {
    m_last_face = f;
  }

  template <class T>
  void set_beta(T)
  {}

  //called when the walk is over but the last segment is on an existing edge
  void on_walk_end_on_edge(halfedge_descriptor) {}

  /// indicates that hp1 is on m_hp2, with the same orientation.
  /// if hp1 is not initialized, it indicates that m_source and m_target
  /// or both on the edge and a geometric predicates need to be used
  /// to get a correct hp1
  void on_input_edge(halfedge_descriptor h)
  {
    m_walked_halfedges.push_back(h);
    target_vertices.push_back(Inter_info());
  }

  const std::vector<halfedge_descriptor>&
  walked_halfedges() const
  {
    return m_walked_halfedges;
  }

  bool do_break()
  {
//    if (cur_dist >= m_dist)
//      std::cout << "stop at " << cur_dist << " (vs. " << m_dist << ")\n";
    return cur_dist >= m_dist;
  }

  Point_3 end_point()
  {
    FT alpha = 0;
    if(is_wrapping())
    {
      alpha = 1;
    }
    else
    {
      FT wanted_dist = m_dist - last_dist;
      FT gotten_dist = cur_dist - last_dist;
      alpha = wanted_dist / gotten_dist;
    }

    CGAL_assertion(alpha<=1 && alpha >=0);
    Point_3 b0;
    if(target_vertices.size() <2)
      b0 = text_center;
    else
      b0 = target_vertices[target_vertices.size() - 2].point();
    Point_3 b1;
    if(target_vertices.size() == 0)
      b1 = text_center;
    else
      b1= target_vertices.back().point();
    return b0+alpha*(b1-b0);
  }

  bool is_wrapping()
  {
    return (cur_dist > 0 && (m_dist > cur_dist));
  }

  template <class OutputIterator>
  OutputIterator get_intersection_points(OutputIterator out) const
  {
    for (const Inter_info& i : target_vertices)
      *out++=i.point();
    return out;
  }
};

} // end of Remove_caps_impl


// 0: negative side of axis_2 and 1: positive side of axis_2
template <class K, class TriangleMesh, class Vpm>
boost::optional< std::array<typename boost::graph_traits<TriangleMesh>::halfedge_descriptor, 2> >
get_sorted_intersected_halfedges(typename boost::graph_traits<TriangleMesh>::face_descriptor start_face,
                                 const TriangleMesh& tm,
                                 const Vpm& vpm,
                                 const typename K::Plane_3& axis_1,
                                 const typename K::Plane_3& axis_2)
{
  typedef typename boost::graph_traits<TriangleMesh> BGT;
  typedef typename BGT::halfedge_descriptor halfedge_descriptor;

  halfedge_descriptor h = halfedge(start_face, tm);
  std::array<typename K::Point_3, 3> triangle_points = {
    get(vpm, source(h, tm)),
    get(vpm, target(h, tm)),
    get(vpm, target(next(h, tm), tm))
  };
  std::array<Oriented_side, 3> oriented_sides;
  for (int i=0; i<3; ++i)
  {
    oriented_sides[i] = axis_1.oriented_side( triangle_points[i] );
    // TODO: handle when the plane goes through a vertex/edge
    if (oriented_sides[i] == ON_ORIENTED_BOUNDARY) return boost::none;
  }

  std::array<halfedge_descriptor, 2> int_edges;
  int k=-1;
  if (oriented_sides[0]!=oriented_sides[1]) int_edges[++k]=h;
  if (oriented_sides[1]!=oriented_sides[2]) int_edges[++k]=next(h, tm);
  if (oriented_sides[2]!=oriented_sides[0]) int_edges[++k]=prev(h, tm);

  if (k!=1) return boost::none;

  // TODO: avoid the constructions
  boost::optional<boost::variant<typename K::Point_3, typename K::Segment_3> >
    inter_res = intersection( typename K::Segment_3(get(vpm, source(int_edges[0], tm)),
                                                    get(vpm, target(int_edges[0], tm))),
                              axis_1 );
  CGAL_assertion( inter_res != boost::none );
  const typename K::Point_3* pt  = boost::get<typename K::Point_3>(&(*inter_res));
  CGAL_assertion( pt!=NULL );

  if (axis_2.oriented_side( *pt )!= ON_NEGATIVE_SIDE)
    std::swap(int_edges[0], int_edges[1]);
  return int_edges;
}

template <class K, class TriangleMesh, class NamedParameters>
std::size_t
two_side_walk_and_intersection_point_collection(const TriangleMesh& tm,
                                                const typename K::Point_3& text_center,
                                                typename boost::graph_traits<TriangleMesh>::face_descriptor start_face,
                                                const typename K::Plane_3& axis_1,
                                                const typename K::Plane_3& axis_2,
                                                double width,
                                                std::vector<typename K::Point_3>& inter_points,
                                                const NamedParameters& np)
{
  typedef boost::graph_traits<TriangleMesh> BGT;
  typedef typename BGT::halfedge_descriptor halfedge_descriptor;
  using parameters::choose_parameter;
  using parameters::get_parameter;

  // Vertex point maps
  typedef typename GetVertexPointMap<TriangleMesh,
      NamedParameters>::const_type Vpm;
  Vpm vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(boost::vertex_point, tm));

  typedef walker_internal::Edge_intersection_with_plane<TriangleMesh, Vpm> Predicate;
  typedef walker_internal::Walk_in_polygon_mesh_visitor<TriangleMesh, Vpm> Visitor;

  boost::optional< std::array<halfedge_descriptor, 2> > opt_int_edges =
    get_sorted_intersected_halfedges<K>(start_face, tm, vpm, axis_1, axis_2);

  if (opt_int_edges == boost::none)
    return 0;

  // 1) walk width/2 to the left
  //    first compute the out point from the center triangle (TODO: make more robust)
  boost::optional<boost::variant<typename K::Point_3, typename K::Segment_3> >
  inter_res = intersection( typename K::Segment_3(get(vpm, source((*opt_int_edges)[0], tm)),
                                                  get(vpm, target((*opt_int_edges)[0], tm))),
                            axis_1 );
  if ( inter_res == boost::none )
  {
    CGAL_warning(!"Intersection 1 not found");
    return 0;
  }
  const typename K::Point_3* pt = boost::get<typename K::Point_3>(&(*inter_res));
  if ( pt == nullptr )
  {
    CGAL_warning(!"Intersection 1 is not a point");
    return 0;
  }

  //    walk only if needed
  if ( CGAL::compare_squared_distance(text_center, *pt, width*width / 4.0) == SMALLER )
  {
    Predicate predicate(text_center, axis_1.orthogonal_vector(),
                                     axis_2.orthogonal_vector(), tm, vpm);
    Visitor visitor = Visitor(tm, vpm, width/2.0, text_center);
    visitor.register_initial_intersection({(*opt_int_edges)[0], *pt});
    walk_in_polygon_mesh(tm,
                         (*opt_int_edges)[0], CGAL::POLYHEDRON_EDGE,
                         BGT::null_halfedge(), CGAL::POLYHEDRON_NONE,
                         predicate, visitor);
    visitor.get_intersection_points(std::back_inserter(inter_points));
    inter_points.back() = visitor.end_point(); // replace the last point by the real endpoint of the walk
    // reverse sequence
    std::reverse(inter_points.begin(), inter_points.end());
  }
  else
    inter_points.push_back(*pt);

  // add the center + save its position
  std::size_t center_pos = inter_points.size();
  if (inter_points.back()!=text_center)
    inter_points.push_back(text_center);
  else
    --center_pos;

  //2)walk width/2 to the right
  //    first compute the out point from the center triangle, other direction (TODO: make more robust)
  inter_res = intersection( typename K::Segment_3(get(vpm, source((*opt_int_edges)[1], tm)),
                                                  get(vpm, target((*opt_int_edges)[1], tm))),
                              axis_1 );
  if ( inter_res == boost::none )
  {
    CGAL_warning(!"Intersection 2 not found");
    return 0;
  }
  pt = boost::get<typename K::Point_3>(&(*inter_res));
  if ( pt == nullptr )
  {
    CGAL_warning(!"Intersection 2 is not a point");
    return 0;
  }

  //    walk only if needed
  if ( CGAL::compare_squared_distance(text_center, *pt, width*width / 4.0) == SMALLER )
  {
    Predicate predicate(text_center, axis_1.orthogonal_vector(),
                                    -axis_2.orthogonal_vector(), tm, vpm);
    Visitor visitor(tm, vpm, width/2.0, text_center);
    visitor.register_initial_intersection({(*opt_int_edges)[1], *pt});

    walk_in_polygon_mesh(tm,
                         (*opt_int_edges)[1], CGAL::POLYHEDRON_EDGE,
                         BGT::null_halfedge(), CGAL::POLYHEDRON_NONE,
                         predicate, visitor);
    visitor.get_intersection_points(std::back_inserter(inter_points));
    inter_points.back() = visitor.end_point(); // replace the last point by the real endpoint of the walk
  }
  else
    inter_points.push_back(*pt);

  return center_pos;
}

template <class TriangleMesh, class ECM>
struct Weight_map{
  const TriangleMesh& tm;
  const ECM& ecm;
  typedef double value_type;
  typedef double reference;
  typedef boost::readable_property_map_tag category;
  typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor key_type;

  Weight_map(const TriangleMesh& tm, const ECM& ecm)
    : tm(tm)
    , ecm(ecm)
  {}

  friend
  double get(const Weight_map& wm, key_type ed)
  {
    if (get(wm.ecm, ed)) return std::numeric_limits<double>::max(); // do not use constrained edges (no walk back);
    return std::sqrt( squared_distance(wm.tm.point(source(ed, wm.tm)), wm.tm.point(target(ed, wm.tm))) );
  }
};

template <class K, class Tree, class OutputIterator, class TriangleMesh, class NamedParameters>
OutputIterator
select_to_walk( const TriangleMesh& tm,
                const Tree& aabb_tree,
                typename K::Point_3 text_center,
                typename K::Plane_3 axis_1,
                typename K::Plane_3 axis_2,
                double width,
                double height,
                OutputIterator out,
                const NamedParameters& np)
{
  typedef boost::graph_traits<TriangleMesh> BGT;
  typedef typename BGT::halfedge_descriptor halfedge_descriptor;
  typedef typename BGT::face_descriptor face_descriptor;
  typedef typename BGT::vertex_descriptor vertex_descriptor;
  typedef typename BGT::edge_descriptor edge_descriptor;
  using boost::choose_param;
  using boost::get_param;

  // Vertex point maps
  typedef typename GetVertexPointMap<TriangleMesh,
      NamedParameters>::const_type Vpm;
  Vpm vpm = choose_param(get_param(np, internal_np::vertex_point),
                         get_const_property_map(boost::vertex_point, tm));

  typedef walker_internal::Edge_intersection_with_plane<TriangleMesh, Vpm> Predicate;
  typedef walker_internal::Walk_in_polygon_mesh_visitor<TriangleMesh, Vpm> Visitor;

  face_descriptor start_face = aabb_tree.closest_point_and_primitive(text_center).second;

  // code use to get the top/bottom, left/right points surrounding the letter
  // on the input mesh
#if 0
  typedef typename K::Plane_3 Plane_3;
  typedef typename K::Point_3 Point_3;

  std::array<halfedge_descriptor, 2> int_edges =
    get_sorted_intersected_halfedges<K>(start_face, tm, vpm, axis_1, axis_2);

  // 1) walk width/2 to the right
  Predicate predicate(text_center, axis_1.orthogonal_vector(),
                                   axis_2.orthogonal_vector(), tm, vpm);
  Visitor visitor = Visitor(tm, vpm, width/2.0, text_center);
  Point_3 new_p;
  walk_in_polygon_mesh(tm,
                       int_edges[1], CGAL::POLYHEDRON_EDGE,
                       BGT::null_halfedge(), CGAL::POLYHEDRON_NONE,
                       predicate, visitor);
  new_p= visitor.end_point();
  std::cout << new_p << "\n";
  Plane_3 r_plane(new_p, -axis_2.orthogonal_vector());
  //2)walk width/2 to the left
  predicate =  Predicate(text_center, axis_1.orthogonal_vector(),
                                     -axis_2.orthogonal_vector(), tm, vpm);
  visitor = Visitor(tm, vpm, width/2.0, text_center);
  //get last point orientation to know where to walk next.
  walk_in_polygon_mesh(tm,
                       int_edges[0], CGAL::POLYHEDRON_EDGE,
                       BGT::null_halfedge(), CGAL::POLYHEDRON_NONE,
                       predicate, visitor);

  new_p= visitor.end_point();
  std::cout << new_p << "\n";
  //bool is_l_wrapping = visitor.is_wrapping();
  Plane_3 l_plane(new_p, axis_2.orthogonal_vector());

  int_edges = get_sorted_intersected_halfedges<K>(start_face, tm, vpm, axis_2, axis_1);
  //3)walk height /2 to the top
  predicate =  Predicate(text_center, axis_2.orthogonal_vector(),
                                      axis_1.orthogonal_vector(),tm, vpm);
  visitor = Visitor(tm, vpm, height/2.0, text_center);
  walk_in_polygon_mesh(tm,
                       int_edges[1], CGAL::POLYHEDRON_EDGE,
                       BGT::null_halfedge(), CGAL::POLYHEDRON_NONE,
                       predicate, visitor);
  new_p= visitor.end_point();
  std::cout << new_p << "\n";
  //bool is_t_wrapping = visitor.is_wrapping();
  Plane_3 t_plane(new_p, -axis_1.orthogonal_vector());
  //4)walk height /2 to the bottom
  predicate = Predicate(text_center, axis_2.orthogonal_vector(),
                                    -axis_1.orthogonal_vector(), tm, vpm);
  visitor = Visitor(tm, vpm, height/2.0, text_center);
  walk_in_polygon_mesh(tm,
                       int_edges[0], CGAL::POLYHEDRON_EDGE,
                       BGT::null_halfedge(), CGAL::POLYHEDRON_NONE,
                       predicate, visitor);
  new_p= visitor.end_point();
  std::cout << new_p << "\n";
  //bool is_b_wrapping = visitor.is_wrapping();
  Plane_3 b_plane(new_p, axis_1.orthogonal_vector());
  CGAL_assertion( l_plane.has_on_positive_side(r_plane.point()) &&
                  t_plane.has_on_positive_side(b_plane.point()) );
#endif

  // walk in the input mesh to get the 4 corners of the bbox of the letter
  // Then using, these corner we extract the next halfedge that would have
  // been intersected if the walk wasn't ended. Using these edges, we are
  // trying to use Dijkstra to close the selection area.
  // The result is not satisfactory on meshes with low triangle resolution
  std::vector<vertex_descriptor> targets(4), sources(4);

  typename K::Vector_3 tmp_axis_1 = axis_1.orthogonal_vector(),
                       tmp_axis_2 = axis_2.orthogonal_vector();
  tmp_axis_1 /= std::sqrt(tmp_axis_1*tmp_axis_1);
  tmp_axis_2 /= std::sqrt(tmp_axis_2*tmp_axis_2);

  typename K::Vector_3 v1 = tmp_axis_1 + tmp_axis_2;
  typename K::Vector_3 v2 = -tmp_axis_1 + tmp_axis_2;

  std::array<halfedge_descriptor, 2> int_edges =
    get_sorted_intersected_halfedges<K>(start_face, tm, vpm, typename K::Plane_3(text_center,v1), typename K::Plane_3(text_center,v2));

  double D =  std::sqrt( ( width * width  + height * height )/ 4. );

  halfedge_descriptor hedge = boost::graph_traits<TriangleMesh>::null_halfedge();

  typedef typename boost::property_map<TriangleMesh, dynamic_edge_property_t<bool> >::const_type ECM;
  ECM ecm = get(dynamic_edge_property_t<bool>(), tm);
  for (edge_descriptor ed : edges(tm))
    put(ecm, ed, false);

  // 1) walk to the corner 0
  Predicate predicate(text_center, v1,
                                   v2, tm, vpm);
  Visitor visitor = Visitor(tm, vpm, D, text_center);
  walk_in_polygon_mesh(tm,
                       int_edges[1], CGAL::POLYHEDRON_EDGE,
                       BGT::null_halfedge(), CGAL::POLYHEDRON_NONE,
                       predicate, visitor);
  hedge = visitor.walked_halfedges().back();
  put(ecm, edge(hedge, tm), true);
  targets[0] = target(hedge, tm);
  sources[3] = source(hedge, tm);
#ifndef NO_DEBUG
  std::cout << "* " << visitor.end_point() << "\n";
  std::cout << "2 " << tm.point(source(hedge, tm)) << " " << tm.point(target(hedge, tm)) << "\n";
#endif
  //2)walk to the corner 2
  predicate =  Predicate(text_center, v1,
                                     -v2, tm, vpm);
  visitor = Visitor(tm, vpm, D, text_center);
  walk_in_polygon_mesh(tm,
                       int_edges[0], CGAL::POLYHEDRON_EDGE,
                       BGT::null_halfedge(), CGAL::POLYHEDRON_NONE,
                       predicate, visitor);
  hedge = visitor.walked_halfedges().back();
  put(ecm, edge(hedge, tm), true);
  targets[2] = target(hedge, tm);
  sources[1] = source(hedge, tm);
#ifndef NO_DEBUG
  std::cout << "* " << visitor.end_point() << "\n";
  std::cout << "2 " << tm.point(source(hedge, tm)) << " " << tm.point(target(hedge, tm)) << "\n";
#endif
  int_edges = get_sorted_intersected_halfedges<K>(start_face, tm, vpm, typename K::Plane_3(text_center,v2), typename K::Plane_3(text_center,v1));
  //3)walk to the corner 1
  predicate =  Predicate(text_center, v2,
                                      v1,tm, vpm);
  visitor = Visitor(tm, vpm, D, text_center);
  walk_in_polygon_mesh(tm,
                       int_edges[1], CGAL::POLYHEDRON_EDGE,
                       BGT::null_halfedge(), CGAL::POLYHEDRON_NONE,
                       predicate, visitor);
  hedge = visitor.walked_halfedges().back();
  put(ecm, edge(hedge, tm), true);
  targets[1] = target(hedge, tm);
  sources[0] = source(hedge, tm);
#ifndef NO_DEBUG
  std::cout << "* " << visitor.end_point() << "\n";
  std::cout << "2 " << tm.point(source(hedge, tm)) << " " << tm.point(target(hedge, tm)) << "\n";
#endif
  //4)walk to the corner 3
  predicate = Predicate(text_center, v2,
                                    -v1, tm, vpm);
  visitor = Visitor(tm, vpm, D, text_center);
  walk_in_polygon_mesh(tm,
                       int_edges[0], CGAL::POLYHEDRON_EDGE,
                       BGT::null_halfedge(), CGAL::POLYHEDRON_NONE,
                       predicate, visitor);
  hedge = visitor.walked_halfedges().back();
  put(ecm, edge(hedge, tm), true);
  targets[3] = target(hedge, tm);
  sources[2] = source(hedge, tm);
#ifndef NO_DEBUG
  std::cout << "* " << visitor.end_point() << "\n";
  std::cout << "2 " << tm.point(source(hedge, tm)) << " " << tm.point(target(hedge, tm)) << "\n";
#endif

  class Dijkstra_end_exception : public std::exception
  {
    const char* what() const throw ()
    {
      return "Dijkstra shortest path: reached the target vertex.";
    }
  };

  class Stop_at_target_Dijkstra_visitor : boost::default_dijkstra_visitor
  {
    vertex_descriptor destination_vd;

  public:
    Stop_at_target_Dijkstra_visitor(vertex_descriptor destination_vd)
      : destination_vd(destination_vd)
    { }

    void initialize_vertex(const vertex_descriptor& /*s*/, const TriangleMesh& /*mesh*/) const { }
    void examine_vertex(const vertex_descriptor& /*s*/, const TriangleMesh& /*mesh*/) const { }
    void examine_edge(const edge_descriptor& /*e*/, const TriangleMesh& /*mesh*/) const { }
    void edge_relaxed(const edge_descriptor& /*e*/, const TriangleMesh& /*mesh*/) const { }
    void discover_vertex(const vertex_descriptor& /*s*/, const TriangleMesh& /*mesh*/) const { }
    void edge_not_relaxed(const edge_descriptor& /*e*/, const TriangleMesh& /*mesh*/) const { }
    void finish_vertex(const vertex_descriptor &vd, const TriangleMesh& /* mesh*/) const
    {
      if(vd == destination_vd)
        throw Dijkstra_end_exception();
    }
  };

  Weight_map<TriangleMesh, ECM> wm(tm, ecm);

// close the loop
  for (std::size_t i=0; i<4; ++i)
  {
    typedef boost::unordered_map<vertex_descriptor, vertex_descriptor>     Pred_umap;
    typedef boost::associative_property_map<Pred_umap>                     Pred_pmap;
    vertex_descriptor src = sources[i], tgt=targets[i];
    if (src==tgt) continue; // TODO: another case to handle is when hedge was the same for 2 consecutive corners
    Pred_umap pred_umap;
    Pred_pmap pred_pmap(pred_umap);
    Stop_at_target_Dijkstra_visitor vis(tgt);

    try{
      boost::dijkstra_shortest_paths(tm, src,
                                     boost::predecessor_map(pred_pmap).visitor(vis).weight_map(wm));
    }
    catch(const Dijkstra_end_exception&)
    {}

    vertex_descriptor prev = tgt;
    do{
      CGAL_assertion( pred_umap.count(prev) == 1 );
      vertex_descriptor vd  = pred_umap[prev];
      halfedge_descriptor hd = halfedge(vd, prev, tm).first;
      CGAL_assertion( hd != boost::graph_traits<TriangleMesh>::null_halfedge() );
      put(ecm, edge(hd, tm), true);
#ifndef NO_DEBUG
      std::cout << "2 " << tm.point(prev) << " " << tm.point(vd) << "\n";
#endif
      prev = vd;
    }
    while(prev != src);
    std::cout << "--\n";
  }

#ifndef NO_DEBUG
  for (edge_descriptor ed : edges(tm))
    if (get(ecm, ed))
      std::cout << "2 " << tm.point(source(ed, tm)) << " " << tm.point(target(ed, tm)) << "\n";
#endif
  connected_component(start_face, tm, out, parameters::edge_is_constrained_map(ecm));

#if 1

#if 0


// select face by expension from the center face, cross an edge if one of its endpoint is inside the selection
  //TODO: use dynamic maps

  std::array<const Plane_3*, 4> planes = { &r_plane, &l_plane, &t_plane, &b_plane };

  std::vector<bool> visited_faces(num_faces(tm),false),
                    visited_vertices(num_vertices(tm), false);
  std::vector < std::bitset<4> >
                    vertices_pos(vertices(tm).size());

  std::vector<halfedge_descriptor> queue;
  visited_faces[start_face]=true;
  *out++=start_face;
  queue.push_back( halfedge(start_face, tm) );
  queue.push_back( next(queue.back(), tm) );
  queue.push_back( next(queue.back(), tm) );

  while(!queue.empty())
  {
    halfedge_descriptor h = opposite(queue.back(), tm);
    queue.pop_back();
    face_descriptor f = face(h, tm);
    if (f==BGT::null_face()) continue;
    if (visited_faces[f]) continue;
    visited_faces[f] = true;

    std::array<std::bitset<4>, 3> f_v_selected;
    for (int i=0; i<3; ++i)
    {
      vertex_descriptor v = source(h, tm);
      if ( !visited_vertices[v] )
      {
        for (int k=0;k<4;++k)
          vertices_pos[v][k] = !planes[k]->has_on_negative_side(get(vpm,v));
        visited_vertices[v] = true;
      }
      f_v_selected[i] = vertices_pos[v];
      h = next(h, tm);
    }
    *out++=f;

    // cross any edge having a point in the plane intersection
    // or if it is intersected by at least two consecutive planes
    std::bitset<4> bs = f_v_selected[1] ^ f_v_selected[2];
    if ( f_v_selected[1].all() || f_v_selected[2].all() |
         (bs[0] && bs[2]) || (bs[0] && bs[3]) || (bs[1] && bs[2]) || (bs[1] && bs[3]) )
    {
      queue.push_back(next(h,tm));
    }
    bs = f_v_selected[0] ^ f_v_selected[2];
    if ( f_v_selected[0].all() || f_v_selected[2].all() ||
         (bs[0] && bs[2]) || (bs[0] && bs[3]) || (bs[1] && bs[2]) || (bs[1] && bs[3]) )
    {
      queue.push_back(prev(h,tm));
    }
  }
#endif
#else
  bool w_planes_are_reversed = l_plane.has_on_negative_side(r_plane.point());
  bool h_planes_are_reversed = t_plane.has_on_negative_side(b_plane.point());

  //if a walk on 1 direction made a whole turn of a cc (e.g. wrapping a sphere ),
  //disable the corresponding planes test, as they are combined
  bool w_wrapping = false;
  bool h_wrapping = false;

  std::vector<face_descriptor> all_selected_faces;
  CGAL::Triangle_from_face_descriptor_map<TriangleMesh, Vpm> triangulator(&tm, vpm) ;
  //first pass : all inside 4 planes
  for(auto fd : faces(tm))
  {
    if (fd==start_face)
    {
      all_selected_faces.push_back(fd);
      continue;
    }
    bool do_select = false;
    //get faces with at least 1 point inside planes
    for(auto vd : CGAL::vertices_around_face(halfedge(fd, tm), tm))
    {
      Point_3 p = get(vpm, vd);
      if(
         (w_wrapping ||
          (w_planes_are_reversed
           ? (r_plane.has_on_negative_side(p) && l_plane.has_on_negative_side(p))
           : (r_plane.has_on_positive_side(p) && l_plane.has_on_positive_side(p)))
          ) && (h_wrapping ||
                (h_planes_are_reversed
                 ? (t_plane.has_on_negative_side(p) && b_plane.has_on_negative_side(p))
                 : (t_plane.has_on_positive_side(p) && b_plane.has_on_positive_side(p)) )))
      {
        do_select = true;
        break;
      }
    }
    //if face is not selected : check if intersecting a plane (to avoid no point between point -> no selection)
    if(!do_select)
    {
      if(!w_wrapping &&
         (CGAL::do_intersect(l_plane, get(triangulator, fd)) ||
          CGAL::do_intersect(r_plane, get(triangulator, fd))))
      {
        for(auto vd : CGAL::vertices_around_face(halfedge(fd, tm), tm)){
          Point_3 p = get(vpm, vd);
          if( h_planes_are_reversed
              ? (t_plane.has_on_negative_side(p) && b_plane.has_on_negative_side(p))
              : (t_plane.has_on_positive_side(p) && b_plane.has_on_positive_side(p)) ){
            do_select = true;
            break;
          }
        }
      }
      if(!h_wrapping &&
         (CGAL::do_intersect(b_plane, get(triangulator, fd)) ||
          CGAL::do_intersect(t_plane, get(triangulator, fd))))
      {
        for(auto vd : CGAL::vertices_around_face(halfedge(fd, tm), tm)){
          Point_3 p = get(vpm, vd);
          if( w_planes_are_reversed
              ? (r_plane.has_on_negative_side(p) && l_plane.has_on_negative_side(p))
              : (r_plane.has_on_positive_side(p) && l_plane.has_on_positive_side(p)) ){
            do_select = true;
            break;
          }
        }
      }
    }
    if(do_select)
      all_selected_faces.push_back(fd);
  }
  //add a "same-CC" constraint:the one that contains start_face

  //create a mesh from the selected faces
  CGAL::Face_filtered_graph<TriangleMesh> ffg(tm, all_selected_faces);
  TriangleMesh subsm;
  face_descriptor subsm_center_face;
  CGAL::copy_face_graph(ffg, subsm,
                        CGAL::Emptyset_iterator(), CGAL::Emptyset_iterator(),
                        boost::make_function_output_iterator(
                        [&subsm_center_face, start_face, tm](
                          const std::pair<face_descriptor, face_descriptor>& p)
                        {
                          if (p.first==start_face){ subsm_center_face=p.second; }
                        }));
  //get fcc map
  std::vector<std::size_t> fccmap(num_faces(tm));  // TODO: use dynamic map
  CGAL::Polygon_mesh_processing::connected_components(ffg,CGAL::make_property_map(fccmap));

  //create ffg from CC containing the start_face and give its faces as output
  CGAL::Face_filtered_graph<CGAL::Face_filtered_graph<TriangleMesh> > selected_cc(ffg,fccmap[start_face], CGAL::make_property_map(fccmap));

  //cut the band
  bool take_r = (width > height);
  //filter out the rest
  for( auto fd : faces(selected_cc))
  {
    bool do_select = true;
    //to avoid non sphere topo
    Plane_3 test_plane(text_center, r_plane.orthogonal_vector());
    if(!take_r)
      test_plane = Plane_3(text_center, t_plane.orthogonal_vector());

    if(do_select)
      *out++ = fd;
  }
#endif
  return out;
}

}
}
#endif // CGAL_POLYGON_MESH_PROCESSING_WALK_TO_SELECT_H
