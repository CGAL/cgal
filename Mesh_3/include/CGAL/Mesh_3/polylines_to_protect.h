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
#include <CGAL/tuple.h>
#include <CGAL/Image_3.h>
#include <CGAL/boost/graph/split_graph_into_polylines.h>
#include <CGAL/internal/Mesh_3/Graph_manipulations.h>
#include <boost/graph/adjacency_list.hpp>
#include <CGAL/Labeled_mesh_domain_3.h> // for CGAL::Null_subdomain_index
#include <boost/utility.hpp> // for boost::prior
#include <boost/optional.hpp>

#include <CGAL/Search_traits_3.h>
#include <CGAL/Orthogonal_incremental_neighbor_search.h>

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
  template <typename vertex_descriptor, typename Graph>
  bool operator()(vertex_descriptor& v, const Graph& g) const
  {
    typedef typename boost::graph_traits<Graph>::out_edge_iterator out_edge_iterator;
    if (g[v].force_terminal || out_degree(v, g) != 2)
      return true;
    else
    {
      out_edge_iterator out_edge_it, out_edges_end;
      boost::tie(out_edge_it, out_edges_end) = out_edges(v, g);

      vertex_descriptor v1 = target(*out_edge_it++, g);
      vertex_descriptor v2 = target(*out_edge_it++, g);
      CGAL_assertion(out_edge_it == out_edges_end);

      const typename Kernel::Point_3& p = g[v].point;
      const typename Kernel::Point_3& p1 = g[v1].point;
      const typename Kernel::Point_3& p2 = g[v2].point;

      if(CGAL::angle(p1, p, p2) == CGAL::ACUTE) {
        // const typename Kernel::Vector_3 e1 = p1 - p;
        // const typename Kernel::Vector_3 e2 = p2 - p;
        // std::cerr << "At point " << p << ": the angle is "
        //           << ( std::acos(e1 * e2
        //                          / CGAL::sqrt(e1*e1)
        //                          / CGAL::sqrt(e2*e2))
        //                * 180 / CGAL_PI ) << std::endl;
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
      if(boost::next(poly_it->begin()) != poly_it->end()) {
        tree.insert(*boost::prior(poly_it->end()));
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
          typename Image_word_type,
          typename Null_subdomain_index,
          typename DomainFunctor,
          typename InterpolationFunctor,
          typename PolylineInputIterator>
void
polylines_to_protect
(const CGAL::Image_3& cgal_image,
 const double vx, const double vy, const double vz,
 std::vector<std::vector<P> >& polylines,
 Image_word_type*,
 Null_subdomain_index null,
 DomainFunctor domain_fct,
 InterpolationFunctor interpolate,
 PolylineInputIterator existing_polylines_begin,
 PolylineInputIterator existing_polylines_end,
 boost::optional<Image_word_type> scalar_interpolation_value = boost::none,
 int prec = 10)
{
  typedef typename DomainFunctor::result_type Domain_type;
  typedef typename Kernel_traits<P>::Kernel K;
  typedef P Point_3;

  using CGAL::internal::polylines_to_protect_namespace::Vertex_info;
  typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS,
                                Vertex_info<Point_3> > Graph;

  typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex_descriptor;
  // typedef typename boost::graph_traits<Graph>::edge_iterator edge_iterator;

  const int xdim = static_cast<int>(cgal_image.xdim());
  const int ydim = static_cast<int>(cgal_image.ydim());
  const int zdim = static_cast<int>(cgal_image.zdim());
  const int image_dims[3] = { xdim, ydim, zdim };

  Graph graph;

  using namespace CGAL::Mesh_3::internal;

  typedef Graph_manipulations<Graph,
                              Point_3,
                              Image_word_type,
                              InterpolationFunctor> G_manip;

  G_manip g_manip(graph, interpolate);

  const float& tx = cgal_image.image()->tx;
  const float& ty = cgal_image.image()->ty;
  const float& tz = cgal_image.image()->tz;
  double max_squared_distance =
    (std::max)( (std::max)(vx, vy), vz );
  max_squared_distance *= max_squared_distance;
  max_squared_distance *= 2;

  typedef Insert_curve_in_graph<G_manip, vertex_descriptor, K> Insert_c_in_g;
  Insert_c_in_g insert_curve_in_graph(g_manip, K(), prec,
                                      max_squared_distance);
  typedef typename Insert_c_in_g::Coords Coords;

  std::size_t
    case4 = 0, // 4 colors
    case211 = 0, case121 = 0, // 3 colors
    case31 = 0, case22 = 0, // 2 colors
    case1 = 0; // 1 color

  typename K::Construct_midpoint_3 midpoint =
    K().construct_midpoint_3_object();
  typename K::Construct_vector_3 vector =
    K().construct_vector_3_object();
  typename K::Construct_translated_point_3 translate =
    K().construct_translated_point_3_object();

  for(int axis = 0; axis < 3; ++axis)
  {
#ifdef CGAL_MESH_3_DEBUG_POLYLINES_TO_PROTECT
    std::cerr << "axis = " << axis << "\n";
#endif // CGAL_MESH_3_DEBUG_POLYLINES_TO_PROTECT
    for(int i = 0; i < xdim; i+= (axis == 0 ? (std::max)(1, xdim-1) : 1 ) )
      for(int j = 0; j < ydim; j+= (axis == 1 ? (std::max)(1, ydim-1) : 1 ) )
        for(int k = 0; k < zdim; k+= (axis == 2 ? (std::max)(1, zdim-1) : 1 ) )
        {

          using std::array;
          using std::tuple;
          using std::get;

          typedef array<int, 3> Pixel;

#ifdef CGAL_MESH_3_DEBUG_POLYLINES_TO_PROTECT
          std::cerr << "Pixel(" << i << ", " << j << ", " << k << ")\n";
#endif // CGAL_MESH_3_DEBUG_POLYLINES_TO_PROTECT
          Pixel pix00 = {{i  , j  , k  }},
            pix10 = pix00, pix01 = pix00, pix11 = pix00;

          const int axis_xx = (axis + 1) % 3;
          const int axis_yy = (axis + 2) % 3;

          ++pix10[axis_xx];
          ++pix11[axis_xx]; ++pix11[axis_yy];
          ++pix01[axis_yy];
          if(pix11[0] >= xdim || pix11[1] >= ydim || pix11[2] >= zdim) {
            // we have gone too far
#ifdef CGAL_MESH_3_DEBUG_POLYLINES_TO_PROTECT
            std::cerr << "  ... continue\n";
#endif // CGAL_MESH_3_DEBUG_POLYLINES_TO_PROTECT
            continue;
          }

          typedef Enriched_pixel<Pixel,
                                 Point_3,
                                 Domain_type,
                                 Image_word_type> Enriched_pixel;

          array<array<Enriched_pixel, 2>, 2> square =
            {{ {{ { pix00, Point_3(), Domain_type(), 0, false },
                  { pix01, Point_3(), Domain_type(), 0, false } }},
               {{ { pix10, Point_3(), Domain_type(), 0, false },
                  { pix11, Point_3(), Domain_type(), 0, false } }} }};

          std::map<Domain_type, int> pixel_values_set;
          for(int ii = 0; ii < 2; ++ii) {
            for(int jj = 0; jj < 2; ++jj) {
              const Pixel& pixel = square[ii][jj].pixel;
              double x = pixel[0] * vx + tx;
              double y = pixel[1] * vy + ty;
              double z = pixel[2] * vz + tz;
              square[ii][jj].on_edge_of_the_cube =
                ( ( ( 0 == pixel[0] || (xdim - 1) == pixel[0] ) ? 1 : 0 )
                  +
                  ( ( 0 == pixel[1] || (ydim - 1) == pixel[1] ) ? 1 : 0 )
                  +
                  ( ( 0 == pixel[2] || (zdim - 1) == pixel[2] ) ? 1 : 0 ) > 1 );
#ifdef CGAL_MESH_3_DEBUG_POLYLINES_TO_PROTECT
              if(square[ii][jj].on_edge_of_the_cube) {
                std::cerr << "  Pixel(" << pixel[0] << ", " << pixel[1] << ", "
                          << pixel[2] << ") is on edge\n";
              }
#endif // CGAL_MESH_3_DEBUG_POLYLINES_TO_PROTECT

              square[ii][jj].point = Point_3(x, y, z);
              square[ii][jj].word =
                CGAL::IMAGEIO::static_evaluate<Image_word_type>
                (cgal_image.image(),
                 pixel[0],
                 pixel[1],
                 pixel[2]);
              square[ii][jj].domain = domain_fct(square[ii][jj].word);
              if(scalar_interpolation_value != boost::none) {
                square[ii][jj].word =
                  Image_word_type(square[ii][jj].word -
                                  (*scalar_interpolation_value));
              }
              ++pixel_values_set[square[ii][jj].domain];
            }
          }

          const Point_3& p00 = square[0][0].point;
          const Point_3& p10 = square[1][0].point;
          const Point_3& p01 = square[0][1].point;
          const Point_3& p11 = square[1][1].point;

          const Image_word_type& v00 = square[0][0].word;
          const Image_word_type& v10 = square[1][0].word;
          const Image_word_type& v01 = square[0][1].word;
          const Image_word_type& v11 = square[1][1].word;

          {
            bool out00 = null(square[0][0].domain);
            bool out10 = null(square[1][0].domain);
            bool out01 = null(square[0][1].domain);
            bool out11 = null(square[1][1].domain);

            bool on_edge00 = square[0][0].on_edge_of_the_cube;
            bool on_edge10 = square[1][0].on_edge_of_the_cube;
            bool on_edge01 = square[0][1].on_edge_of_the_cube;
            bool on_edge11 = square[1][1].on_edge_of_the_cube;

            //
            // Protect the edges of the cube
            //
            if(pix00[axis_xx] == 0 &&
               ! ( out00 && out01 ) )
            {
              g_manip.try_add_edge(g_manip.get_vertex(p00, on_edge00),
                                   g_manip.get_vertex(p01, on_edge01));
            }
            if(pix11[axis_xx] == image_dims[axis_xx]-1 &&
               ! ( out10 && out11 ) )
            {
              g_manip.try_add_edge(g_manip.get_vertex(p10, on_edge10),
                                   g_manip.get_vertex(p11, on_edge11));
            }
            if(pix00[axis_yy] == 0 &&
               ! ( out00 && out10 ) )
            {
              g_manip.try_add_edge(g_manip.get_vertex(p00, on_edge00),
                                   g_manip.get_vertex(p10, on_edge10));
            }
            if(pix11[axis_yy] == image_dims[axis_yy]-1 &&
               ! ( out01 && out11 ) )
            {
              g_manip.try_add_edge(g_manip.get_vertex(p01, on_edge01),
                                   g_manip.get_vertex(p11, on_edge11));
            }
          } // end of scope for outIJ and on_edgeIJ

          //
          // Protect lines inside the square
          //
          switch(pixel_values_set.size()) {
          case 4: {
            CGAL_assertion(square[0][0].domain != square[0][1].domain);
            CGAL_assertion(square[0][0].domain != square[1][0].domain);
            CGAL_assertion(square[0][0].domain != square[1][1].domain);
            CGAL_assertion(square[1][0].domain != square[1][1].domain);
            CGAL_assertion(square[0][1].domain != square[1][1].domain);
            CGAL_assertion(square[0][1].domain != square[1][0].domain);

case_4:
            // case 4 or case 2-2
            //
            // ---------------
            // |      |      |
            // |      |      |
            // |______|______|
            // |      |      |
            // |      |      |
            // |      |      |
            // ---------------
            //
            ++case4;
            vertex_descriptor left   = g_manip.split(square[0][0],
                                                     square[0][1], null);
            vertex_descriptor right  = g_manip.split(square[1][0],
                                                     square[1][1], null);
            vertex_descriptor top    = g_manip.split(square[0][1],
                                                     square[1][1], null);
            vertex_descriptor bottom = g_manip.split(square[0][0],
                                                     square[1][0], null);

            vertex_descriptor vmid = g_manip.split(graph[left].point,
                                                   graph[right].point,
                                                   v00, v10,
                                                   false, false,
                                                   false, false);
            g_manip.try_add_edge(left   , vmid);
            g_manip.try_add_edge(right  , vmid);
            g_manip.try_add_edge(top    , vmid);
            g_manip.try_add_edge(bottom , vmid);
          }
            break;
          case 3: {
            if(square[0][0].domain == square[1][1].domain) {
              // Diagonal, but the wrong one.
              // Vertical swap
              std::swap(square[0][1], square[0][0]);
              std::swap(square[1][1], square[1][0]);
            }
            if(square[0][1].domain == square[1][0].domain) {
              // diagonal case 1-2-1
              //
              // A-------------C
              // |        |    |
              // |         \   |
              // |___       \__|
              // |   \         |
              // |    \        |
              // |     |       |
              // B-------------A
              //
              // two curves: -1 ------------- 0   -1 ------------- 1
              //               |             |      |        |    |
              //               |             |      |         \   |
              //               |___          |      |          \__|
              //               |   \         |      |             |
              //               |    \        |      |             |
              //               |     |       |      |             |
              //              1 ------------- -1   0 ------------- -1
              //
              CGAL_assertion(square[0][1].domain == square[1][0].domain);
              CGAL_assertion(square[1][1].domain != square[0][0].domain);
              CGAL_assertion(square[0][1].domain != square[0][0].domain);
              CGAL_assertion(square[0][1].domain != square[1][1].domain);
case_1_2_1:
              ++case121;
              vertex_descriptor left   = g_manip.split(square[0][0],
                                                       square[0][1], null);
              vertex_descriptor right  = g_manip.split(square[1][0],
                                                       square[1][1], null);
              vertex_descriptor top    = g_manip.split(square[0][1],
                                                       square[1][1], null);
              vertex_descriptor bottom = g_manip.split(square[0][0],
                                                       square[1][0], null);

              Isoline_equation equation =
                (scalar_interpolation_value == boost::none) ?
                Isoline_equation(1, -1, -1, 0) :
                Isoline_equation(v00, v10, v01, v11);
              insert_curve_in_graph.insert_curve(equation,
                                                 left, bottom,
                                                 Coords(0., interpolate(v00, v01) ),
                                                 Coords(interpolate(v00, v10), 0. ),
                                                 p00,
                                                 p10 - p00,
                                                 p01 - p00);
              if(scalar_interpolation_value == boost::none) {
                equation = Isoline_equation(0, -1, -1, 1);
              }
              insert_curve_in_graph.insert_curve(equation,
                                                 top, right,
                                                 Coords(interpolate(v01, v11), 1. ),
                                                 Coords(1., interpolate(v10, v11) ),
                                                 p00,
                                                 p10 - p00,
                                                 p01 - p00);
            } else {
              // case 2-1-1
              if(square[0][0].domain == square[1][0].domain) {
                // Diagonal swap
                std::swap(square[0][1], square[1][0]);
              } else
              if(square[0][1].domain == square[1][1].domain) {
                // The other diagonal swap
                std::swap(square[0][0], square[1][1]);
              } else
              if(square[1][0].domain == square[1][1].domain) {
                // Vertical swap
                std::swap(square[0][0], square[1][0]);
                std::swap(square[0][1], square[1][1]);
              }
              CGAL_assertion(square[0][0].domain == square[0][1].domain);
              CGAL_assertion(square[0][0].domain != square[1][0].domain);
              CGAL_assertion(square[0][0].domain != square[1][1].domain);
              CGAL_assertion(square[1][0].domain != square[1][1].domain);
              ++case211;
              // Note: this case cannot occur with non-segmented scalar
              // images, because it needs three domains.
              //
              // A-------------B
              // |        |    |
              // |         \   |
              // |          \__|
              // |          /  |
              // |         /   |
              // |        |    |
              // A-------------C
              //
              // Two curves portions:
              //
              //  1 ------------- 0     1 ------------- -1
              //   |            /|       |       |     |
              //   |           / |       |        \    |
              //   |          /  |       |         \   |
              //   |         /   |       |          \  |
              //   |        /    |       |           \ |
              //   |       |     |       |            \|
              //  1 ------------- -1    1 ------------- 0
              //
              // ...and a segment.
              Point_3 midleft =  midpoint(p00, p01);
              Point_3 midright = midpoint(p10, p11);
              Point_3 inter = translate(midleft
                              , (2./3) * vector(midleft, midright));
              vertex_descriptor v_inter = g_manip.get_vertex(inter, false);
              vertex_descriptor right  = g_manip.split(square[1][0],
                                                       square[1][1], null);
              vertex_descriptor top    = g_manip.split(square[0][1],
                                                       square[1][1], null);
              vertex_descriptor bottom = g_manip.split(square[0][0],
                                                       square[1][0], null);

              insert_curve_in_graph.insert_curve(Isoline_equation(1, -1, 1, 0),
                                                 bottom, v_inter,
                                                 Coords( 0.5  , 0.  ),
                                                 Coords( 2./3., 0.5 ),
                                                 p00,
                                                 p10 - p00,
                                                 p01 - p00);
              insert_curve_in_graph.insert_curve(Isoline_equation(1, 0, 1, -1),
                                                 top, v_inter,
                                                 Coords( 0.5  , 1.  ),
                                                 Coords( 2./3., 0.5 ),
                                                 p00,
                                                 p10 - p00,
                                                 p01 - p00);
              g_manip.try_add_edge(right,  v_inter);
            } // end case 2-1-1
          } // end `case 3:`
            break;
          case 2: {
            if(pixel_values_set.begin()->second ==
               pixel_values_set.rbegin()->second)
            {
              // Case of two colors with two pixels each.

              if(square[0][0].domain==square[1][0].domain) {
                // case 2-2, diagonal swap
                std::swap(square[0][1], square[1][0]);
                CGAL_assertion(square[0][0].domain==square[0][1].domain);
              }
              if(square[1][0].domain==square[1][1].domain) {
                // case 2-2, vertical swap
                std::swap(square[0][1], square[1][1]);
                std::swap(square[0][0], square[1][0]);
                CGAL_assertion(square[0][0].domain==square[0][1].domain);
              }
              if(square[0][1].domain==square[1][1].domain) {
                // case 2-2, diagonal swap
                std::swap(square[0][0], square[1][1]);
                CGAL_assertion(square[0][0].domain==square[0][1].domain);
              }

              if(square[0][0].domain==square[0][1].domain) {
                // vertical case 2-2
                ++case22;
                // case 2-2
                //
                // A-------------B
                // |      |      |
                // |      |      |
                // |      |      |
                // |      |      |
                // |      |      |
                // |      |      |
                // A-------------B
                //
                CGAL_assertion(square[1][0].domain==square[1][1].domain);
                CGAL_assertion(square[1][0].domain!=square[0][1].domain);
                vertex_descriptor top    = g_manip.split(square[0][1],
                                                         square[1][1], null);
                vertex_descriptor bottom = g_manip.split(square[0][0],
                                                         square[1][0], null);
                if(scalar_interpolation_value == boost::none) {
                  g_manip.try_add_edge(top, bottom);
                } else {
                  insert_curve_in_graph.insert_curve(Isoline_equation(v00, v10,
                                                                      v01, v11),
                                                     top, bottom,
                                                     Coords(interpolate(v01,v11), 1.),
                                                     Coords(interpolate(v00,v10), 0.),
                                                     p00,
                                                     p10 - p00,
                                                     p01 - p00);
                }
              } else {
                // Else diagonal case case 2-2
                // Same as the case with 4 colors
                CGAL_assertion(square[0][0].domain==square[1][1].domain);
                CGAL_assertion(square[1][0].domain==square[0][1].domain);
                CGAL_assertion(square[0][0].domain!=square[0][1].domain);

                if(scalar_interpolation_value != boost::none) {
                  // Compute the squared distance between the two branches of
                  // the hyperbola.
                  const double discrimant = double(v00) * v11 - double(v01) * v10;
                  const double squared_distance =
                    8 * CGAL::abs(discrimant) / CGAL::square(double(v00) - v10 - v01 + v11);
#ifdef CGAL_MESH_3_DEBUG_POLYLINES_TO_PROTECT
                  std::cerr << "squared_distance: " << squared_distance << "\n";
#endif // CGAL_MESH_3_DEBUG_POLYLINES_TO_PROTECT
                  if(CGAL::square(prec) * squared_distance > 1)
                  {
                    // In that case, the case 1-2-1 will be applied
                    if(discrimant > 0.) {
                      // Vertical swap
                      std::swap(square[0][1], square[0][0]);
                      std::swap(square[1][1], square[1][0]);
                    }
                    goto case_1_2_1;
                  }
                } // scalar interpolation

                goto case_4;
              }
            }
            else {
              // case of two colors with one pixel green and three red
              Domain_type value_alone;
              if(pixel_values_set.begin()->second == 1) {
                value_alone = pixel_values_set.begin()->first;
              } else {
                CGAL_assertion(pixel_values_set.begin()->second == 3);
                CGAL_assertion(pixel_values_set.rbegin()->second ==1);
                value_alone = pixel_values_set.rbegin()->first;
              }
              if(square[0][1].domain == value_alone) {
                // central symmetry
                std::swap(square[0][1], square[1][0]);
                std::swap(square[0][0], square[1][1]);
                CGAL_assertion(square[1][0].domain == value_alone);
              }
              if(square[1][1].domain == value_alone) {
                // vertical swap
                std::swap(square[0][0], square[0][1]);
                std::swap(square[1][0], square[1][1]);
                CGAL_assertion(square[1][0].domain == value_alone);
              }
              if(square[0][0].domain == value_alone) {
                // horizontal swap
                std::swap(square[0][1], square[1][1]);
                std::swap(square[0][0], square[1][0]);
                CGAL_assertion(square[1][0].domain == value_alone);
              }
              ++case31;
              //
              // A-------------A    1 ------------- 1
              // |             |     |             |
              // |             |     |             |
              // |         ____|     |         ____|
              // |        /    |     |        /    |
              // |       /     |     |       /     |
              // |      |      |     |      |      |
              // A-------------B    1 ------------- -1
              //
              CGAL_assertion(square[1][0].domain == value_alone);
              CGAL_assertion(square[1][0].domain != square[0][0].domain);
              CGAL_assertion(square[1][1].domain == square[0][0].domain);
              CGAL_assertion(square[0][1].domain == square[0][0].domain);
              vertex_descriptor bottom = g_manip.split(square[0][0],
                                                       square[1][0], null);
              vertex_descriptor right  = g_manip.split(square[1][0],
                                                       square[1][1], null);

              Isoline_equation equation =
                (scalar_interpolation_value == boost::none) ?
                Isoline_equation(1, -1, 1, 1) :
                Isoline_equation(v00, v10, v01, v11);

              insert_curve_in_graph.insert_curve(equation,
                                                 bottom, right,
                                                 Coords(interpolate(v00, v10), 0. ),
                                                 Coords(1., interpolate(v10, v11) ),
                                                 p00,
                                                 p10 - p00,
                                                 p01 - p00);
            }
          }
            break;
          default: // case 1
            ++case1;
            // nothing to do
            break;
          }
        }
  }
  // std::cerr << "case 4:     " << case4 << std::endl;
  // std::cerr << "case 2-1-1: " << case211 << std::endl;
  // std::cerr << "case 1-2-1: " << case121 << std::endl;
  // std::cerr << "case 3-1:   " << case31 << std::endl;
  // std::cerr << "case 2-2:   " << case22 << std::endl;
  // std::cerr << "case 1:     " << case1 << std::endl;

  const std::ptrdiff_t nb_facets =
    case4 + case211 + case121 + case31 + case22 + case1;
  const std::ptrdiff_t expected_nb_facets =
    ((xdim != 1 && ydim != 1 && zdim != 1) ? 2 : 1)
    *
    ((xdim-1)*(ydim-1) + (ydim-1)*(zdim-1) + (xdim-1)*(zdim-1));

  // std::cerr << "nb of facets:           " << nb_facets << std::endl
  //           << " expected nb of facets: " << expected_nb_facets << std::endl;

  CGAL_assertion(nb_facets == expected_nb_facets);
  CGAL_USE(nb_facets); CGAL_USE(expected_nb_facets);

  snap_graph_vertices(graph,
                      vx, vy, vz,
                      existing_polylines_begin, existing_polylines_end,
                      K());

  Mesh_3::Polyline_visitor<Point_3, Graph> visitor(polylines, graph);
  Less_for_Graph_vertex_descriptors<Graph> less(graph);
  const Graph& const_graph = graph;
  split_graph_into_polylines(const_graph, visitor,
                             Mesh_3::Angle_tester<K>(), less);
}

template <typename P,
          typename PolylineInputIterator>
void
polylines_to_protect(std::vector<std::vector<P> >& polylines,
                     PolylineInputIterator existing_polylines_begin,
                     PolylineInputIterator existing_polylines_end)
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
    Polyline polyline = *poly_it;
    if (polyline.size() < 2)
      continue;

    typename Polyline::iterator pit = polyline.begin();
    while (boost::next(pit) != polyline.end())
    {
      vertex_descriptor v = g_manip.get_vertex(*pit, false);
      vertex_descriptor w = g_manip.get_vertex(*boost::next(pit), false);
      g_manip.try_add_edge(v, w);
      ++pit;
    }
  }

  Mesh_3::Polyline_visitor<Point_3, Graph> visitor(polylines, graph);
  Less_for_Graph_vertex_descriptors<Graph> less(graph);
  const Graph& const_graph = graph;
  typedef typename Kernel_traits<P>::Kernel K;
  split_graph_into_polylines(const_graph, visitor,
                             Mesh_3::Angle_tester<K>(), less);
}

template <typename P, typename Image_word_type, typename Null_subdomain_index>
void
polylines_to_protect(const CGAL::Image_3& cgal_image,
                     std::vector<std::vector<P> >& polylines,
                     Image_word_type* word_type,
                     Null_subdomain_index null)
{
  polylines_to_protect<P>
    (cgal_image,
     polylines,
     word_type,
     null,
     0,
     0);
}

template <typename P, typename Image_word_type>
void
polylines_to_protect(const CGAL::Image_3& cgal_image,
                     std::vector<std::vector<P> >& polylines)
{
  polylines_to_protect<P, Image_word_type>(cgal_image, polylines, 0, 0);
}

template <typename P,
          typename Image_word_type,
          typename Null_subdomain_index,
          typename PolylineInputIterator>
void
polylines_to_protect(const CGAL::Image_3& cgal_image,
                     std::vector<std::vector<P> >& polylines,
                     Image_word_type* word_type,
                     Null_subdomain_index null,
                     PolylineInputIterator existing_polylines_begin,
                     PolylineInputIterator existing_polylines_end)
{
  typedef typename Kernel_traits<P>::Kernel K;
  polylines_to_protect<P>
    (cgal_image,
     cgal_image.vx(), cgal_image.vy(),cgal_image.vz(),
     polylines,
     word_type,
     null,
     CGAL::Identity<Image_word_type>(),
     Mesh_3::internal::Returns_midpoint<K, Image_word_type>(),
     existing_polylines_begin,
     existing_polylines_end);
}

template <typename P,
          typename Image_word_type,
          typename PolylineInputIterator>
void
polylines_to_protect(const CGAL::Image_3& cgal_image,
                     std::vector<std::vector<P> >& polylines,
                     PolylineInputIterator existing_polylines_begin,
                     PolylineInputIterator existing_polylines_end)
{
  polylines_to_protect<P>
    (cgal_image,
     polylines,
     (Image_word_type*)0,
     CGAL::Null_subdomain_index(),
     existing_polylines_begin,
     existing_polylines_end);
}

} // namespace CGAL

#endif // CGAL_MESH_3_POLYLINES_TO_PROTECT_H
