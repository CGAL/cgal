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

#ifndef CGAL_MESH_3_POLYLINES_TO_PROTECT_IN_IMAGE_H
#define CGAL_MESH_3_POLYLINES_TO_PROTECT_IN_IMAGE_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/Mesh_3/polylines_to_protect.h>
#include <CGAL/Image_3.h>

namespace CGAL {

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
 std::optional<Image_word_type> scalar_interpolation_value = std::nullopt,
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
            {{ {{ { pix00, Point_3(), Domain_type(), 0, false, false },
                  { pix01, Point_3(), Domain_type(), 0, false, false } }},
               {{ { pix10, Point_3(), Domain_type(), 0, false, false },
                  { pix11, Point_3(), Domain_type(), 0, false, false } }} }};

          std::map<Domain_type, int> pixel_values_set;
          for(int ii = 0; ii < 2; ++ii) {
            for(int jj = 0; jj < 2; ++jj) {
              const Pixel& pixel = square[ii][jj].pixel;
              double x = pixel[0] * vx + tx;
              double y = pixel[1] * vy + ty;
              double z = pixel[2] * vz + tz;
              short sum_faces = ((0 == pixel[0] || (xdim - 1) == pixel[0]) ? 1 : 0)
                              + ((0 == pixel[1] || (ydim - 1) == pixel[1]) ? 1 : 0)
                              + ((0 == pixel[2] || (zdim - 1) == pixel[2]) ? 1 : 0);
              square[ii][jj].on_edge_of_the_cube = (sum_faces > 1);
              square[ii][jj].on_corner_of_the_cube = (sum_faces > 2);

#ifdef CGAL_MESH_3_DEBUG_POLYLINES_TO_PROTECT
              if(square[ii][jj].on_edge_of_the_cube) {
                std::cerr << "  Pixel(" << pixel[0] << ", " << pixel[1] << ", "
                          << pixel[2] << ") is on edge\n";
              }
              if (square[ii][jj].on_corner_of_the_cube) {
                std::cerr << "  Pixel(" << pixel[0] << ", " << pixel[1] << ", "
                  << pixel[2] << ") is on corner\n";
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
              if(scalar_interpolation_value != std::nullopt) {
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

            bool is_corner00 = square[0][0].on_corner_of_the_cube;
            bool is_corner10 = square[1][0].on_corner_of_the_cube;
            bool is_corner01 = square[0][1].on_corner_of_the_cube;
            bool is_corner11 = square[1][1].on_corner_of_the_cube;

            //
            // Protect the edges of the cube
            //
            if(pix00[axis_xx] == 0 &&
               ! ( out00 && out01 ) )
            {
              g_manip.try_add_edge(g_manip.get_vertex(p00, is_corner00),
                                   g_manip.get_vertex(p01, is_corner01));
            }
            if(pix11[axis_xx] == image_dims[axis_xx]-1 &&
               ! ( out10 && out11 ) )
            {
              g_manip.try_add_edge(g_manip.get_vertex(p10, is_corner10),
                                   g_manip.get_vertex(p11, is_corner11));
            }
            if(pix00[axis_yy] == 0 &&
               ! ( out00 && out10 ) )
            {
              g_manip.try_add_edge(g_manip.get_vertex(p00, is_corner00),
                                   g_manip.get_vertex(p10, is_corner10));
            }
            if(pix11[axis_yy] == image_dims[axis_yy]-1 &&
               ! ( out01 && out11 ) )
            {
              g_manip.try_add_edge(g_manip.get_vertex(p01, is_corner01),
                                   g_manip.get_vertex(p11, is_corner11));
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
                (scalar_interpolation_value == std::nullopt) ?
                Isoline_equation(1, -1, -1, 0) :
                Isoline_equation(v00, v10, v01, v11);
              insert_curve_in_graph.insert_curve(equation,
                                                 left, bottom,
                                                 Coords(0., interpolate(v00, v01) ),
                                                 Coords(interpolate(v00, v10), 0. ),
                                                 p00,
                                                 p10 - p00,
                                                 p01 - p00);
              if(scalar_interpolation_value == std::nullopt) {
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
                if(scalar_interpolation_value == std::nullopt) {
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
                // Else diagonal case 2-2
                // Same as the case with 4 colors
                CGAL_assertion(square[0][0].domain==square[1][1].domain);
                CGAL_assertion(square[1][0].domain==square[0][1].domain);
                CGAL_assertion(square[0][0].domain!=square[0][1].domain);

                if(scalar_interpolation_value != std::nullopt) {
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
                (scalar_interpolation_value == std::nullopt) ?
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

template <typename P,
          typename Image_word_type,
          typename PolylineInputIterator>
void
polylines_to_protect_on_bbox(const CGAL::Image_3& cgal_image,
                     std::vector<std::vector<P> >& polylines,
                     PolylineInputIterator existing_polylines_begin,
                     PolylineInputIterator existing_polylines_end)
{
  polylines_to_protect<P, Image_word_type>(cgal_image,
                                           polylines,
                                           existing_polylines_begin,
                                           existing_polylines_end);
}



template <typename PolylineRange1, typename PolylineRange2>
void
merge_and_snap_polylines(const CGAL::Image_3& image,
                         PolylineRange1& polylines_to_snap,
                         const PolylineRange2& existing_polylines)
{
  static_assert(std::is_same<typename PolylineRange1::value_type::value_type,
                             typename PolylineRange2::value_type::value_type>::value,
                "Polyline ranges should have same point type");
  using P = typename PolylineRange1::value_type::value_type;
  using K = typename Kernel_traits<P>::Kernel;

  using CGAL::internal::polylines_to_protect_namespace::Vertex_info;
  using Graph = boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS,
                                      Vertex_info<P> >;
  using vertex_descriptor = typename boost::graph_traits<Graph>::vertex_descriptor;

  // build graph of polylines_to_snap
  Graph graph;
  typedef Mesh_3::internal::Returns_midpoint<K, int> Midpoint_fct;
  Mesh_3::internal::Graph_manipulations<Graph,
    P,
    int,
    Midpoint_fct> g_manip(graph);

  for (const auto& polyline : polylines_to_snap)
  {
    if (polyline.size() < 2)
      continue;

    auto pit = polyline.begin();
    while (std::next(pit) != polyline.end())
    {
      vertex_descriptor v = g_manip.get_vertex(*pit, false);
      vertex_descriptor w = g_manip.get_vertex(*std::next(pit), false);
      g_manip.try_add_edge(v, w);
      ++pit;
    }
  }

  // snap graph to existing_polylines
  snap_graph_vertices(graph,
    image.vx(), image.vy(), image.vz(),
    std::begin(existing_polylines), std::end(existing_polylines),
    K());

  // rebuild polylines_to_snap
  polylines_to_snap.clear();
  Mesh_3::Polyline_visitor<P, Graph> visitor(polylines_to_snap, graph);
  Less_for_Graph_vertex_descriptors<Graph> less(graph);
  const Graph& const_graph = graph;
  Mesh_3::Angle_tester<K> angle_tester(90.);
  split_graph_into_polylines(const_graph, visitor, angle_tester, less);
}

} // CGAL namespace

#endif
