// Copyright (c) 2015 GeometryFactory
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
// Author(s)     : Laurent Rineau

#ifndef CGAL_POLYLINES_TO_PROTECT_H
#define CGAL_POLYLINES_TO_PROTECT_H

#include <vector>
#include <map>
#include <CGAL/tuple.h>
#include <CGAL/Image_3.h>
#include <CGAL/internal/Mesh_3/split_in_polylines.h>
#include <CGAL/internal/Mesh_3/Graph_manipulations.h>
#include <boost/graph/adjacency_list.hpp>

namespace CGAL {


  template <typename P>
  void
  polylines_to_protect(const CGAL::Image_3& cgal_image,
                       const double, const double, const double,
                       std::vector<std::vector<P> >& polylines)
  {
    typedef typename Kernel_traits<P>::Kernel K;
    typedef P Point_3;
    typedef boost::adjacency_list<boost::setS, boost::setS, boost::undirectedS, Point_3> Graph;
    typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex_descriptor;
   // typedef typename boost::graph_traits<Graph>::edge_iterator edge_iterator;

    const int xdim = static_cast<int>(cgal_image.xdim());
    const int ydim = static_cast<int>(cgal_image.ydim());
    const int zdim = static_cast<int>(cgal_image.zdim());

  const int image_dims[3] = { xdim, ydim, zdim };

  Graph graph;
  internal::Mesh_3::Graph_manipulations<Graph, Point_3> g_manip(graph);

  std::size_t
    case4 = 0, // 4 colors
    case211 = 0, case121 = 0, // 3 colors
    case31 = 0, case22 = 0, // 2 colors
    case1 = 0; // 1 color

  for(int axis = 0; axis < 3; ++axis)
  {
    for(int i = 0; i < xdim; i+= (axis == 0 ? (xdim-1) : 1 ) )
      for(int j = 0; j < ydim; j+= (axis == 1 ? (ydim-1) : 1 ) )
        for(int k = 0; k < zdim; k+= (axis == 2 ? (zdim-1) : 1 ) )
        {

          using CGAL::cpp11::array;
          using CGAL::cpp11::tuple;

          typedef array<int, 3> Pixel;

          Pixel pix00 = {{i  , j  , k  }},
            pix10 = pix00, pix01 = pix00, pix11 = pix00;

          const int axis_xx = (axis + 1) % 3;
          const int axis_yy = (axis + 2) % 3;

          ++pix10[axis_xx];
          ++pix11[axis_xx];
          ++pix01[axis_yy];
          ++pix11[axis_yy];
          if(pix11[0] >= xdim || pix11[1] >= ydim || pix11[2] >= zdim) {
            // we have gone too far
            continue;
          }
          typedef unsigned char Image_word_type;
          typedef tuple<Pixel, Point_3, Image_word_type> Enriched_pixel;

          array<array<Enriched_pixel, 2>, 2> square =
            {{ {{ Enriched_pixel(pix00, Point_3(), Image_word_type()),
                  Enriched_pixel(pix01, Point_3(), Image_word_type()) }},
               {{ Enriched_pixel(pix10, Point_3(), Image_word_type()),
                  Enriched_pixel(pix11, Point_3(), Image_word_type()) }} }};

          std::map<Image_word_type, int> pixel_values_set;
          for(int ii = 0; ii < 2; ++ii)
            for(int jj = 0; jj < 2; ++jj)
            {
              const Pixel& pixel = get<0>(square[ii][jj]);
              double x = pixel[0] * cgal_image.vx();
              double y = pixel[1] * cgal_image.vy();
              double z = pixel[2] * cgal_image.vz();
              get<1>(square[ii][jj]) = Point_3(x, y, z);
              get<2>(square[ii][jj]) = static_cast<Image_word_type>(cgal_image.value(pixel[0],
                                                                                     pixel[1],
                                                                                     pixel[2]));
              ++pixel_values_set[get<2>(square[ii][jj])];
            }
          const Point_3& p00 = get<1>(square[0][0]);
          const Point_3& p10 = get<1>(square[1][0]);
          const Point_3& p01 = get<1>(square[0][1]);
          const Point_3& p11 = get<1>(square[1][1]);

          //
          // Protect the edges of the cube
          //
          if(pix00[axis_xx] == 0) {
            g_manip.try_add_edge(g_manip.get_vertex(p00),
                                 g_manip.get_vertex(p01));
          }
          if(pix11[axis_xx] == image_dims[axis_xx]-1) {
            g_manip.try_add_edge(g_manip.get_vertex(p10),
                                 g_manip.get_vertex(p11));
          }
          if(pix00[axis_yy] == 0) {
            g_manip.try_add_edge(g_manip.get_vertex(p00),
                                 g_manip.get_vertex(p10));
          }
          if(pix11[axis_yy] == image_dims[axis_yy]-1) {
            g_manip.try_add_edge(g_manip.get_vertex(p01),
                                 g_manip.get_vertex(p11));
          }

          //
          // Protect lines inside the square
          //
          switch(pixel_values_set.size()) {
          case 4: {
            assert(get<2>(square[0][0]) != get<2>(square[0][1]));
            assert(get<2>(square[0][0]) != get<2>(square[1][0]));
            assert(get<2>(square[0][0]) != get<2>(square[1][1]));
            assert(get<2>(square[1][0]) != get<2>(square[1][1]));
            assert(get<2>(square[0][1]) != get<2>(square[1][1]));
            assert(get<2>(square[0][1]) != get<2>(square[1][0]));
case_4:
            // case 4 or case 2-2
            ++case4;
            vertex_descriptor left   = g_manip.split(p00, p01);
            vertex_descriptor right  = g_manip.split(p10, p11);
            vertex_descriptor top    = g_manip.split(p01, p11);
            vertex_descriptor bottom = g_manip.split(p00, p10);
            vertex_descriptor vmid   = g_manip.get_vertex(midpoint(p00, p11));
            g_manip.try_add_edge(left   , vmid);
            g_manip.try_add_edge(right  , vmid);
            g_manip.try_add_edge(top    , vmid);
            g_manip.try_add_edge(bottom , vmid);
          }
            break;
          case 3: {
            if(get<2>(square[0][0]) == get<2>(square[1][1])) {
              // Diagonal, but the wrong one.
              // Vertical swap
              std::swap(square[0][1], square[0][0]);
              std::swap(square[1][1], square[1][0]);
            }
            if(get<2>(square[0][1]) == get<2>(square[1][0])) {
              // diagonal case 1-2-1
              assert(get<2>(square[0][1]) == get<2>(square[1][0]));
              assert(get<2>(square[1][1]) != get<2>(square[0][0]));
              assert(get<2>(square[0][1]) != get<2>(square[0][0]));
              assert(get<2>(square[0][1]) != get<2>(square[1][1]));
              ++case121;
              vertex_descriptor left   = g_manip.split(p00, p01);
              vertex_descriptor right  = g_manip.split(p10, p11);
              vertex_descriptor top    = g_manip.split(p01, p11);
              vertex_descriptor bottom = g_manip.split(p00, p10);
              Point_3 inter_left = p00 + (1./3.) * (p11 - p00);
              Point_3 inter_right = p00 + (2./3.) * (p11 - p00);
              vertex_descriptor v_int_left = g_manip.get_vertex(inter_left);
              vertex_descriptor v_int_right = g_manip.get_vertex(inter_right);
              g_manip.try_add_edge(left, v_int_left);
              g_manip.try_add_edge(v_int_left, bottom);
              g_manip.try_add_edge(top, v_int_right);
              g_manip.try_add_edge(v_int_right, right);
            } else {
              // case 2-1-1
              if(get<2>(square[0][0]) == get<2>(square[1][0])) {
                // Diagonal swap
                std::swap(square[0][1], square[1][0]);
              } else
              if(get<2>(square[0][1]) == get<2>(square[1][1])) {
                // The other diagonal swap
                std::swap(square[0][0], square[1][1]);
              } else
              if(get<2>(square[1][0]) == get<2>(square[1][1])) {
                // Vertical swap
                std::swap(square[0][0], square[1][0]);
                std::swap(square[0][1], square[1][1]);
              }
              assert(get<2>(square[0][0]) == get<2>(square[0][1]));
              assert(get<2>(square[0][0]) != get<2>(square[1][0]));
              assert(get<2>(square[0][0]) != get<2>(square[1][1]));
              ++case211;
              Point_3 midleft =  midpoint(p00, p01);
              Point_3 midright = midpoint(p10, p11);
              Point_3 inter = midleft + (2./3)*(midright-midleft);
              vertex_descriptor v_inter = g_manip.get_vertex(inter);
              vertex_descriptor right  = g_manip.split(p10, p11);
              vertex_descriptor top    = g_manip.split(p01, p11);
              vertex_descriptor bottom = g_manip.split(p00, p10);
              g_manip.try_add_edge(top,    v_inter);
              g_manip.try_add_edge(bottom, v_inter);
              g_manip.try_add_edge(right,  v_inter);
            }
          }
            break;
          case 2: {
            if(pixel_values_set.begin()->second ==
               pixel_values_set.rbegin()->second)
            {
              // Case of two colors with two pixels each.

              if(get<2>(square[0][0])==get<2>(square[1][0])) {
                // case 2-2, diagonal swap
                std::swap(square[0][1], square[1][0]);
                assert(get<2>(square[0][0])==get<2>(square[0][1]));
              }
              if(get<2>(square[1][0])==get<2>(square[1][1])) {
                // case 2-2, vertical swap
                std::swap(square[0][1], square[1][1]);
                std::swap(square[0][0], square[1][0]);
                assert(get<2>(square[0][0])==get<2>(square[0][1]));
              }
              if(get<2>(square[0][1])==get<2>(square[1][1])) {
                // case 2-2, diagonal swap
                std::swap(square[0][0], square[1][1]);
                assert(get<2>(square[0][0])==get<2>(square[0][1]));
              }

              if(get<2>(square[0][0])==get<2>(square[0][1])) {
                // vertical case 2-2
                ++case22;
                assert(get<2>(square[1][0])==get<2>(square[1][1]));
                assert(get<2>(square[1][0])!=get<2>(square[0][1]));
                vertex_descriptor top    = g_manip.split(p01, p11);
                vertex_descriptor bottom = g_manip.split(p00, p10);
                g_manip.try_add_edge(top, bottom);
              } else {
                // Else diagonal case case 2-2
                // Same as the case with 4 colors
                assert(get<2>(square[0][0])==get<2>(square[1][1]));
                assert(get<2>(square[1][0])==get<2>(square[0][1]));
                assert(get<2>(square[0][0])!=get<2>(square[0][1]));
                goto case_4;
              }
            }
            else {
              // case of two colors with one pixel green and three red
              Image_word_type value_alone;
              if(pixel_values_set.begin()->second == 1) {
                value_alone = pixel_values_set.begin()->first;
              } else {
                assert(pixel_values_set.begin()->second == 3);
                assert(pixel_values_set.rbegin()->second ==1);
                value_alone = pixel_values_set.rbegin()->first;
              }
              if(get<2>(square[0][1]) == value_alone) {
                // central symmetry
                std::swap(square[0][1], square[1][0]);
                std::swap(square[0][0], square[1][1]);
              assert(get<2>(square[1][0]) == value_alone);
              }
              if(get<2>(square[1][1]) == value_alone) {
                // vertical swap
                std::swap(square[0][0], square[0][1]);
                std::swap(square[1][0], square[1][1]);
              assert(get<2>(square[1][0]) == value_alone);
              }
              if(get<2>(square[0][0]) == value_alone) {
                // horizontal swap
                std::swap(square[0][1], square[1][1]);
                std::swap(square[0][0], square[1][0]);
              assert(get<2>(square[1][0]) == value_alone);
              }
              ++case31;
              assert(get<2>(square[1][0]) == value_alone);
              assert(get<2>(square[1][0]) != get<2>(square[0][0]));
              assert(get<2>(square[1][1]) == get<2>(square[0][0]));
              assert(get<2>(square[0][1]) == get<2>(square[0][0]));
              vertex_descriptor bottom = g_manip.split(p00, p10);
              vertex_descriptor old = bottom;

              Point_3 inter;

              vertex_descriptor v_int;
              for(double x = 0.55; x < 1.; x+= 0.05)
              {
                inter = p00
                  +      x         * (p10 - p00)  // x
                  + (1.-1./(2.*x)) * (p01 - p00); // y
                v_int = g_manip.get_vertex(inter);
                g_manip.try_add_edge(old, v_int);
                old = v_int;
              }

              vertex_descriptor right  = g_manip.split(p10, p11);
              g_manip.try_add_edge(v_int, right);
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
    2*((xdim-1)*(ydim-1) + (ydim-1)*(zdim-1) + (xdim-1)*(zdim-1));

  // std::cerr << "nb of facets:           " << nb_facets << std::endl
  //           << " expected nb of facets: " << expected_nb_facets << std::endl;

  CGAL_assertion(nb_facets == expected_nb_facets);
  CGAL_USE(nb_facets); CGAL_USE(expected_nb_facets);


  internal::Mesh_3::split_in_polylines(graph, polylines, K());
  }


  template <typename P>
  void
  polylines_to_protect(const CGAL::Image_3& cgal_image,
                       std::vector<std::vector<P> >& polylines)
  {
    polylines_to_protect<P>
      (cgal_image,
       cgal_image.vx(), cgal_image.vy(),cgal_image.vz(),
       polylines);
  }


} // namespace CGAL


#endif // CGAL_POLYLINES_TO_PROTECT_H

