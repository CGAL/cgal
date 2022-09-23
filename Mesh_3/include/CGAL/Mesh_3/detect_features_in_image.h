// Copyright (c) 2022 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sebastien Loriot, Jane Tournois
//
//******************************************************************************
//
//******************************************************************************

#ifndef CGAL_MESH_3_DETECT_FEATURES_IN_IMAGE_H
#define CGAL_MESH_3_DETECT_FEATURES_IN_IMAGE_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/ImageIO.h>

#include <CGAL/Mesh_3/features_detection/coordinates.h>
#include <CGAL/Mesh_3/features_detection/combinations.h>
#include <CGAL/Mesh_3/features_detection/cases_table.h>
#include <CGAL/Mesh_3/features_detection/features_detection.h>
#include <CGAL/Mesh_3/features_detection/cube_isometries.h>
#include <CGAL/Mesh_3/features_detection/features_detection_helpers.h>

#include <CGAL/Mesh_3/polylines_to_protect.h>

#include <vector>
#include <array>

#ifdef CGAL_DEBUG_TRIPLE_LINES
#include <boost/range/join.hpp>
#endif

namespace CGAL
{
namespace Mesh_3
{

// Protect the intersection of the object with the box of the image,
// by declaring 1D-features. Note that `CGAL::polylines_to_protect` is
// not documented.
template<typename Word_type, typename Mesh_domain>
bool detect_features_in_image_with_know_word_type(const CGAL::Image_3& image,
                                                  Mesh_domain& domain)
{
  using Gt = typename Mesh_domain::R;
  using Point_3 = typename Gt::Point_3;
  using Vector_3 = typename Gt::Vector_3;
  using Polyline_type = std::vector<Point_3>;
  using Polylines = std::vector<Polyline_type>;

  CGAL::Mesh_3::Triple_line_extractor<Point_3> lines;

  Polylines features_inside;

  const double vx = image.vx();
  const double vy = image.vy();
  const double vz = image.vz();
  const double dist_bound = (std::min)(vx, (std::min)(vy, vz)) / 256;
  const double sq_dist_bound = dist_bound * dist_bound;

  const std::size_t xdim = image.xdim();
  const std::size_t ydim = image.ydim();
  const std::size_t zdim = image.zdim();

  const float tx = image.tx();
  const float ty = image.ty();
  const float tz = image.tz();

  using CGAL::IMAGEIO::static_evaluate;

  using Del = CGAL::Delaunay_triangulation_3<Gt>;
  using Cell_handle = typename Del::Cell_handle;
  using Vertex_handle = typename Del::Vertex_handle;
  Del triangulation;
  Cell_handle start_cell;

  using Word //use unsigned integral Word type to use it as an index
    = typename CGAL::IMAGEIO::Word_type_generator<WK_FIXED, SGN_UNSIGNED, sizeof(Word_type)>::type;

  using Color_transform = internal::Color_transformation_helper<Word>;
  typename Color_transform::type color_transformation;
  std::array<Word, 8>            inv_color_transformation;

  using Permutation = internal::Permutation;
  using Coord = internal::Coordinates;

  for (std::size_t k = 0, end_k = zdim - 1; k < end_k; ++k)
    for (std::size_t j = 0, end_j = ydim - 1; j < end_j; ++j)
      for (std::size_t i = 0, end_i = xdim - 1; i < end_i; ++i)
      {
        Vector_3 translation{ i * vx + tx,
                              j * vy + ty,
                              k * vz + tz };

        const std::array<Word, 8> cube = {
          static_evaluate<Word>(image.image(), i  , j  , k),
          static_evaluate<Word>(image.image(), i + 1, j  , k),
          static_evaluate<Word>(image.image(), i  , j + 1, k),
          static_evaluate<Word>(image.image(), i + 1, j + 1, k),
          static_evaluate<Word>(image.image(), i  , j  , k + 1),
          static_evaluate<Word>(image.image(), i + 1, j  , k + 1),
          static_evaluate<Word>(image.image(), i  , j + 1, k + 1),
          static_evaluate<Word>(image.image(), i + 1, j + 1, k + 1),
        }; /// TODO: optimize the access to the image data
        bool monocolor = (cube[0] == cube[1]);
        for (int i = 2; i < 8; ++i) monocolor = monocolor && (cube[0] == cube[i]);
        if (monocolor) continue;

        Color_transform::reset(color_transformation);

        std::uint8_t nb_color = 0;
        for (int i = 0; i < 8; ++i) {
          if (!Color_transform::is_valid(color_transformation, cube[i]))
          {
            color_transformation[cube[i]] = nb_color;
            inv_color_transformation[nb_color] = cube[i];
            ++nb_color;
          }
        }
        std::array<std::uint8_t, 8> reference_cube = {
            color_transformation[cube[0]],
            color_transformation[cube[1]],
            color_transformation[cube[2]],
            color_transformation[cube[3]],
            color_transformation[cube[4]],
            color_transformation[cube[5]],
            color_transformation[cube[6]],
            color_transformation[cube[7]]
        };
        auto case_it = internal::find_case(internal::cases, reference_cube);
        const bool case_found = (case_it != std::end(internal::cases));
        if (case_found) reference_cube = internal::combinations[(*case_it)[8]];
        else {
          //std::cerr << "Warning: case not found: " << reference_cube << '\n';
          CGAL_error();
        };
#ifdef CGAL_DEBUG_TRIPLE_LINES
        CGAL::Mesh_3::internal::debug_cerr("Cube", cube);
        CGAL::Mesh_3::internal::debug_cerr("reference cube", reference_cube);
        CGAL::Mesh_3::internal::debug_cerr("with transformation", internal::cube_isometries[(*case_it)[9]]);
#endif // CGAL_DEBUG_TRIPLE_LINES

        auto fct_it = lines.create_polylines_fcts.find(reference_cube);
        if (fct_it != lines.create_polylines_fcts.end())
        {
#ifdef CGAL_DEBUG_TRIPLE_LINES
          CGAL::Mesh_3::internal::debug_cerr("Using the function of", Cube(fct_it->first));
#endif // CGAL_DEBUG_TRIPLE_LINES

          Polylines cube_features = (fct_it->second)(10);
          if (case_found)
          {
            const Permutation& transformation = internal::cube_isometries[(*case_it)[9]];

            Coord a1 = internal::coordinates[transformation[0]];
            Coord u = internal::minus(internal::coordinates[transformation[1]], a1);
            Coord v = internal::minus(internal::coordinates[transformation[2]], a1);
            Coord w = internal::minus(internal::coordinates[transformation[4]], a1);

            const Point_3  pa{ a1[0], a1[1], a1[2] };
            const Vector_3 vu{ u[0], u[1], u[2] };
            const Vector_3 vv{ v[0], v[1], v[2] };
            const Vector_3 vw{ w[0], w[1], w[2] };
#ifdef CGAL_DEBUG_TRIPLE_LINES
            std::cerr << "pa: " << pa << "\n";
            std::cerr << "vu: " << vu << "\n";
            std::cerr << "vv: " << vv << "\n";
            std::cerr << "vw: " << vw << "\n";
#endif // CGAL_DEBUG_TRIPLE_LINES
            for (auto& polyline : cube_features) {
              for (auto& point : polyline) {
                point = pa
                  + point.x() * vu
                  + point.y() * vv
                  + point.z() * vw;
                point = { vx * point.x(),
                          vy * point.y(),
                          vz * point.z(), };
                point = point + translation;
              }
              for (int i = 0; i < 2; ++i) {
                Point_3& extremity = (i == 0) ? polyline.front() : polyline.back();
                Vertex_handle vh = triangulation.nearest_vertex(extremity, start_cell);
                if (Vertex_handle() != vh) {
                  if (squared_distance(vh->point(), extremity) < sq_dist_bound) {
                    extremity = vh->point();
                  }
                }
                vh = triangulation.insert(extremity, start_cell);
                start_cell = vh->cell();
              }
              features_inside.push_back(std::move(polyline));
            } // end loop on polylines
          } // end case where the transformation is not the identity
        } // end if the reference_cube has polylines
      }

  // call the split_graph_into_polylines, to create long polylines from the
  // short polylines that were generated per voxel.
  Polylines new_polylines_inside;
  CGAL::polylines_to_protect(new_polylines_inside,
                             features_inside.begin(),
                             features_inside.end());

  std::vector<std::vector<Point_3> > polylines_on_bbox;
  CGAL::polylines_to_protect<Point_3, Word_type>(image, polylines_on_bbox,
    new_polylines_inside.begin(),
    new_polylines_inside.end());

  domain.add_features(polylines_on_bbox.begin(), polylines_on_bbox.end());

  // It is very important that the polylines from the file `lines_fname`
  // contain only polylines in the inside of the box of the image.
  domain.add_features(new_polylines_inside.begin(), new_polylines_inside.end());

#ifdef CGAL_DEBUG_TRIPLE_LINES
  std::ofstream output_polylines("out-generated.polylines.txt");
  output_polylines.precision(17);
  for (auto poly : boost::range::join(polylines_on_bbox, new_polylines_inside)) {
    output_polylines << poly.size();
    for (auto p : poly) output_polylines << " " << p;
    output_polylines << std::endl;
  }
#endif

  return true;
}

template<typename Mesh_domain>
bool detect_features_in_image(const CGAL::Image_3& image, Mesh_domain& domain)
{
  CGAL_IMAGE_IO_CASE(image.image(),
    return detect_features_in_image_with_know_word_type<Word>(image, domain)
  );
  CGAL_error_msg("This place should never be reached, because it would mean "
                 "the image word type is a type that is not handled by "
                 "CGAL_ImageIO.");
  return false;
}


}//end namespace Mesh_3
}//end namespace CGAL


#endif //CGAL_MESH_3_DETECT_FEATURES_IN_IMAGE_H
