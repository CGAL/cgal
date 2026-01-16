// Copyright (c) 2015,2016 GeometryFactory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_MESH_3_POLYLINES_TO_PROTECT_H
#define CGAL_MESH_3_POLYLINES_TO_PROTECT_H

#include <CGAL/license/Mesh_3.h>

#include <vector>
#include <map>
#include <utility> // std::swap
#include <algorithm> // std::min

#include <CGAL/number_utils.h>
#include <CGAL/squared_distance_3.h>

#include <CGAL/boost/graph/split_graph_into_polylines.h>
#include <CGAL/Mesh_3/internal/Graph_manipulations.h>
#include <boost/graph/adjacency_list.hpp>

#include <CGAL/Search_traits_3.h>
#include <CGAL/Orthogonal_incremental_neighbor_search.h>

#include <CGAL/Mesh_3/Null_subdomain_index.h>

#include <type_traits>
#include <optional>

namespace CGAL {
namespace Mesh_3 {
namespace internal {

template <typename K, typename NT>
struct Returns_midpoint {
  typedef typename K::Point_3 Point_3;

  Point_3 operator()(const Point_3& a,
                     const Point_3& b,
                     const NT,
                     const NT) const
  {
    typename K::Construct_midpoint_3 midpt
      = K().construct_midpoint_3_object();
    return a < b ? midpt(a, b) : midpt(b, a);
  }
  double operator()(const NT,
                const NT) const
  {
    return 0.5;
  }
};

template <typename K, typename NT>
struct Linear_interpolator {
  typedef typename K::Point_3 Point_3;

  Point_3 interpolate(const Point_3& pa,
                      const Point_3& pb,
                      const NT a,
                      const NT b) const
  {
    return pa + (a / (a - b)) * ( pb - pa);
  }


  Point_3 operator()(const Point_3& pa,
                     const Point_3& pb,
                     const NT a,
                     const NT b) const
  {
    return
      (a < b) ?
      interpolate(pa, pb, a, b) :
      interpolate(pb, pa, b, a);
  }
  double operator()(const NT a,
                    const NT b) const
  {
    return a / ( a - b );
  }
};


class Isoline_equation {
  double a;
  double b;
  double c;
  double d;

public:
  Isoline_equation(double a, double b, double c, double d)
    : a(a)
    , b(b)
    , c(c)
    , d(d)
  {}

  double y(double x) const {
    return ( x*(a-b)-a ) /
      ( c-a+x*(a-b-c+d) );
    // If (a+d) == (b+c), then the curve is a straight line.
  }
  double x(double y) const {
    return ( y*(a-c)-a ) /
      ( b-a+y*(a-b-c+d) );
    // If (a+d) == (b+c), then the curve is a straight line.
  }
}; // end of class Isoline_equation

template <typename G_manip, typename vertex_descriptor, typename Geom_traits>
class Insert_curve_in_graph {
  typedef typename Geom_traits::Point_3  Point;
  typedef typename Geom_traits::Vector_3 Vector;

  G_manip& g_manip;
  const vertex_descriptor null_vertex;
  const Geom_traits& gt;
  const typename Geom_traits::Construct_translated_point_3 translate;
  const typename Geom_traits::Construct_scaled_vector_3 scale;
  const int prec;
  const double max_squared_distance;

public:
  Insert_curve_in_graph(G_manip& g_manip,
                        const Geom_traits& gt,
                        const int prec = 10,
                        const double max_squared_distance = 0)
    : g_manip(g_manip)
    , null_vertex(g_manip.null_vertex())
    , gt(gt)
    , translate(gt.construct_translated_point_3_object())
    , scale(gt.construct_scaled_vector_3_object())
    , prec(prec)
    , max_squared_distance(max_squared_distance)
  {}

  struct Coords {
    Coords(double x, double y) : x(x), y(y) {}
    double x;
    double y;
  };

  void insert_curve(Isoline_equation equation,
                    vertex_descriptor start_v,
                    vertex_descriptor end_v,
                    Coords start,
                    Coords end,
                    Point  p00,
                    Vector vx,
                    Vector vy)
  {
    if(CGAL::abs(start.x - end.x) >= CGAL::abs(start.y - end.y)) {
      insert_curve(equation, &Isoline_equation::y,
                   start_v, end_v,
                   start.x, end.x, p00, vx, vy);
    } else {
      insert_curve(equation, &Isoline_equation::x,
                   start_v, end_v,
                   start.y, end.y, p00, vy, vx);
    }
  }
private:
  void insert_curve(Isoline_equation equation,
                    double (Isoline_equation::* f)(double) const,
                    vertex_descriptor start_v,
                    vertex_descriptor end_v,
                    double start,
                    double end,
                    Point  p00,
                    Vector vx,
                    Vector vy)
  {
#ifdef CGAL_MESH_3_DEBUG_POLYLINES_TO_PROTECT
    std::cerr << "New curve:\n"
              << "  base("
              << p00 << " , "
              << p00+vx << " , "
              << p00+vy << ")\n"
              << "  vectors: "
              << "( " << vx << " ) "
              << " ( " << vy << " )\n";
#endif // CGAL_MESH_3_DEBUG_POLYLINES_TO_PROTECT
    const double step = (end - start)/prec;
    const double stop = end-step/2;
    const bool step_is_positive = (step > 0);
    vertex_descriptor old = start_v;
    vertex_descriptor v_int = null_vertex;
    for(double x = start + step;
        (step_is_positive ? x < stop : x > stop);
        x += step)
    {
      const double y = (equation.*f)(x);
#ifdef CGAL_MESH_3_DEBUG_POLYLINES_TO_PROTECT
      std::cerr << "  (" << x << ", " << y << ") -> ";
#endif // CGAL_MESH_3_DEBUG_POLYLINES_TO_PROTECT
      const Point inter_p =
        translate(translate(p00,
                            scale(vx, x)),
                  scale(vy, y));
#ifdef CGAL_MESH_3_DEBUG_POLYLINES_TO_PROTECT
      std::cerr << "( " << inter_p << ")\n";
#endif // CGAL_MESH_3_DEBUG_POLYLINES_TO_PROTECT
      v_int = g_manip.get_vertex(inter_p, false);
      g_manip.try_add_edge(old, v_int);

      CGAL_assertion_msg(max_squared_distance == 0 ||
                         CGAL::squared_distance(g_manip.g[old].point,
                                                g_manip.g[v_int].point) <
                         max_squared_distance,
                       ([this, old, v_int] {
                         std::stringstream s;
                         s << "Problem at segment ("
                           << this->g_manip.g[old].point << ", "
                           << this->g_manip.g[v_int].point << ")";
                         return s.str();
                       }().c_str()));
      old = v_int;
    }
    if(null_vertex != v_int) {
      // v_int can be null if the curve is degenerated into one point.
      g_manip.try_add_edge(v_int, end_v);
      CGAL_assertion_msg(max_squared_distance == 0 ||
                         CGAL::squared_distance(g_manip.g[end_v].point,
                                                g_manip.g[v_int].point) <
                         max_squared_distance,
                         ([this, end_v, v_int] {
                           std::stringstream s;
                           s << "Problem at segment ("
                             << this->g_manip.g[end_v].point << ", "
                             << this->g_manip.g[v_int].point << ")";
                           return s.str();
                         }().c_str()));
    }
  }
};

template <typename Pixel,
          typename Point,
          typename Domain_type,
          typename Image_word_type>
struct Enriched_pixel {
  Pixel pixel;
  Point point;
  Domain_type domain;
  Image_word_type word;
  bool on_edge_of_the_cube;
  bool on_corner_of_the_cube;
}; // end struct template Enriched_pixel<Pix,P,D,C>

} // end namespace internal

template<typename P, typename G>
struct Polyline_visitor
{
  std::vector<std::vector<P> >& polylines;
  G& graph;
  Polyline_visitor(typename std::vector<std::vector<P> >& lines, G& p_graph)
    : polylines(lines), graph(p_graph)
  {}

  ~Polyline_visitor()
  {//DEBUG
#if CGAL_MESH_3_PROTECTION_DEBUG & 2
    std::ofstream og("polylines_graph.polylines.txt");
    og.precision(17);
    for(const std::vector<P>& poly : polylines)
    {
      og << poly.size() << " ";
      for(const P& p : poly)
        og << p << " ";
      og << std::endl;
    }
#endif // CGAL_MESH_3_PROTECTION_DEBUG & 2
  }

  void start_new_polyline()
  {
    std::vector<P> V;
    polylines.push_back(V);
  }

  void add_node(typename boost::graph_traits<G>::vertex_descriptor vd)
  {
    std::vector<P>& polyline = polylines.back();
    polyline.push_back(graph[vd].point);
  }

  void end_polyline()
  {
    // ignore degenerated polylines
    if(polylines.back().size() < 2)
      polylines.resize(polylines.size() - 1);
  }
};

template <typename Kernel>
struct Angle_tester
{
  const double m_angle_sq_cosine;// squared cosine of `std:max(90, angle_deg)`

  Angle_tester()
    : m_angle_sq_cosine(0)
  {}

  Angle_tester(const double angle_deg)//angle given in degrees for readability
    : m_angle_sq_cosine(CGAL::square(std::cos((std::max)(90.,angle_deg) * CGAL_PI / 180.)))
  {}

  template <typename vertex_descriptor, typename Graph>
  bool operator()(vertex_descriptor& v, const Graph& g) const
  {
    typedef typename boost::graph_traits<Graph>::out_edge_iterator out_edge_iterator;
    if (g[v].force_terminal || out_degree(v, g) != 2)
      return true;
    else
    {
      out_edge_iterator out_edge_it, out_edges_end;
      std::tie(out_edge_it, out_edges_end) = out_edges(v, g);

      vertex_descriptor v1 = target(*out_edge_it++, g);
      vertex_descriptor v2 = target(*out_edge_it++, g);
      CGAL_assertion(out_edge_it == out_edges_end);

      const typename Kernel::Point_3& p = g[v].point;
      const typename Kernel::Point_3& p1 = g[v1].point;
      const typename Kernel::Point_3& p2 = g[v2].point;

      //if angle at v is acute, v must be considered as a terminal vertex
      // to ensure termination
      if (CGAL::angle(p1, p, p2) == CGAL::ACUTE)
        return true;
      else if (m_angle_sq_cosine > 0.)//check angle only if angle is > 90.
      {
        const typename Kernel::Vector_3 e1 = p1 - p;
        const typename Kernel::Vector_3 e2 = p2 - p;

        const auto scalar_product = e1 * e2;
        if (CGAL::is_positive(scalar_product))
          return true;

        const auto sq_scalar_product = CGAL::square(scalar_product);
        if (sq_scalar_product <= m_angle_sq_cosine * (e1 * e1) * (e2 * e2))
          return true;
      }
    }
    return false;
  }
};

}//namespace Mesh_3

// this function is overloaded for when `PolylineInputIterator` is `int`.
template <typename K,
          typename Graph,
          typename PolylineInputIterator>
void snap_graph_vertices(Graph& graph,
                         const double vx, const double vy, const double vz,
                         PolylineInputIterator existing_polylines_begin,
                         PolylineInputIterator existing_polylines_end,
                         K)
{
  const double dist_bound = (std::min)(vx,
                                       (std::min)(vy, vz)) / 256;
  const double sq_dist_bound = dist_bound * dist_bound;

  typedef CGAL::Search_traits_3<K> Tree_traits;
  typedef CGAL::Orthogonal_incremental_neighbor_search<Tree_traits> NN_search;
  typedef typename NN_search::Tree Tree;

  Tree tree;
  // insert the extremities of the polylines in the kd-tree
  for(PolylineInputIterator poly_it = existing_polylines_begin;
      poly_it != existing_polylines_end; ++poly_it)
  {
    if(poly_it->begin() != poly_it->end()) {
      tree.insert(*poly_it->begin());
      if(std::next(poly_it->begin()) != poly_it->end()) {
        tree.insert(*std::prev(poly_it->end()));
      }
    }
  }
  if(tree.size() == 0) return;

  for(typename boost::graph_traits<Graph>::vertex_descriptor v :
                make_range(vertices(graph)))
  {
    const typename K::Point_3 p = graph[v].point;
    NN_search nn(tree, p);
    CGAL_assertion(nn.begin() != nn.end());
    if(squared_distance(nn.begin()->first, p) < sq_dist_bound) {
      graph[v].point = nn.begin()->first;
    }
  }
}

template <typename K,
          typename Graph>
void snap_graph_vertices(Graph&, double, double, double, int, int, K)
{}

template <typename Graph>
struct Less_for_Graph_vertex_descriptors
{
  const Graph& graph;
  Less_for_Graph_vertex_descriptors(const Graph& graph) : graph(graph) {}

  template <typename vertex_descriptor>
  bool operator()(vertex_descriptor v1, vertex_descriptor v2) const {
    return graph[v1].point < graph[v2].point;
  }
}; // end of Less_for_Graph_vertex_descriptors<Graph>

namespace internal {
namespace polylines_to_protect_namespace {
  template <typename Point>
  struct Vertex_info {
    Point point;
    bool force_terminal;
    Vertex_info() : force_terminal(false) {}
    Vertex_info(const Point& p) : point(p), force_terminal(false) {}
  };
} // end namespace polylines_to_protect_namespace
} // end namespace internal

template <typename P,
          typename PolylineInputIterator>
void
polylines_to_protect(std::vector<std::vector<P> >& polylines,
                     PolylineInputIterator existing_polylines_begin,
                     PolylineInputIterator existing_polylines_end,
                     const double& angle = 90.)//when not provided, check only for acute angles
{
  typedef P Point_3;
  typedef typename Kernel_traits<P>::Kernel K;
  using CGAL::internal::polylines_to_protect_namespace::Vertex_info;
  typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS,
                                Vertex_info<Point_3> > Graph;
  typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  typedef typename std::iterator_traits<PolylineInputIterator>::value_type Polyline;

  Graph graph;
  typedef Mesh_3::internal::Returns_midpoint<K, int> Midpoint_fct;
  Mesh_3::internal::Graph_manipulations<Graph,
                                        Point_3,
                                        int,
                                        Midpoint_fct> g_manip(graph);

  for (PolylineInputIterator poly_it = existing_polylines_begin;
       poly_it != existing_polylines_end; ++poly_it)
  {
    const Polyline& polyline = *poly_it;
    if (polyline.size() < 2)
      continue;

    typename Polyline::const_iterator pit = polyline.begin();
    while (std::next(pit) != polyline.end())
    {
      vertex_descriptor v = g_manip.get_vertex(*pit, false);
      vertex_descriptor w = g_manip.get_vertex(*std::next(pit), false);
      g_manip.try_add_edge(v, w);
      ++pit;
    }
  }

  Mesh_3::Polyline_visitor<Point_3, Graph> visitor(polylines, graph);
  Less_for_Graph_vertex_descriptors<Graph> less(graph);
  const Graph& const_graph = graph;
  typedef typename Kernel_traits<P>::Kernel K;
  Mesh_3::Angle_tester<K> angle_tester(angle);
  split_graph_into_polylines(const_graph, visitor,
                             angle_tester, less);
}

} // namespace CGAL

#endif // CGAL_MESH_3_POLYLINES_TO_PROTECT_H
