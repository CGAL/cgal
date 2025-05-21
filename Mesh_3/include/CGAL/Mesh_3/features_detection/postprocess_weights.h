// Copyright (c) 2023 GeometryFactory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jane Tournois

#ifndef CGAL_MESH_3_POSTPROCESS_LABEL_WEIGHTS_H
#define CGAL_MESH_3_POSTPROCESS_LABEL_WEIGHTS_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/Image_3.h>

#include <vector>

namespace CGAL {
namespace Mesh_3 {
namespace internal {

  template<typename Word_type>
  void set_voxel(CGAL::Image_3& img,
                 const std::size_t& i,
                 const std::size_t& j,
                 const std::size_t& k,
                 const Word_type& w)
  {
    using CGAL::IMAGEIO::static_evaluate;

    if (i == std::size_t(-1) || j == std::size_t(-1) || k == std::size_t(-1))
      return;
    else if (i > img.xdim() - 1 || j > img.ydim() - 1 || k > img.zdim() - 1)
      return;
    else
      static_evaluate<Word_type>(img.image(), i, j, k) = w;
  }

  template<typename Word_type>
  void set_voxels(CGAL::Image_3& weights,
                  const std::vector<std::array<std::size_t, 3>>& black_voxels,
                  const Word_type& wblack)
  {
    for (auto v : black_voxels)
    {
      const std::size_t& i = v[0];
      const std::size_t& j = v[1];
      const std::size_t& k = v[2];

      // i - 1 : 9 voxels
      internal::set_voxel(weights, i - 1, j - 1, k - 1, wblack);
      internal::set_voxel(weights, i - 1, j - 1, k, wblack);
      internal::set_voxel(weights, i - 1, j - 1, k + 1, wblack);
      internal::set_voxel(weights, i - 1, j, k - 1, wblack);
      internal::set_voxel(weights, i - 1, j, k, wblack);
      internal::set_voxel(weights, i - 1, j, k + 1, wblack);
      internal::set_voxel(weights, i - 1, j + 1, k - 1, wblack);
      internal::set_voxel(weights, i - 1, j + 1, k, wblack);
      internal::set_voxel(weights, i - 1, j + 1, k + 1, wblack);

      // i : 9 voxels
      internal::set_voxel(weights, i, j - 1, k - 1, wblack);
      internal::set_voxel(weights, i, j - 1, k, wblack);
      internal::set_voxel(weights, i, j - 1, k + 1, wblack);
      internal::set_voxel(weights, i, j, k - 1, wblack);
      internal::set_voxel(weights, i, j, k, wblack);
      internal::set_voxel(weights, i, j, k + 1, wblack);
      internal::set_voxel(weights, i, j + 1, k - 1, wblack);
      internal::set_voxel(weights, i, j + 1, k, wblack);
      internal::set_voxel(weights, i, j + 1, k + 1, wblack);

      // i + 1 : 9 voxels
      internal::set_voxel(weights, i + 1, j - 1, k - 1, wblack);
      internal::set_voxel(weights, i + 1, j - 1, k, wblack);
      internal::set_voxel(weights, i + 1, j - 1, k + 1, wblack);
      internal::set_voxel(weights, i + 1, j, k - 1, wblack);
      internal::set_voxel(weights, i + 1, j, k, wblack);
      internal::set_voxel(weights, i + 1, j, k + 1, wblack);
      internal::set_voxel(weights, i + 1, j + 1, k - 1, wblack);
      internal::set_voxel(weights, i + 1, j + 1, k, wblack);
      internal::set_voxel(weights, i + 1, j + 1, k + 1, wblack);
    }


#ifdef CGAL_MESH_3_WEIGHTED_IMAGES_DEBUG
    _writeImage(weights.image(), "weights-image_postprocessed.inr.gz");
#endif
  }

  // POLYLINES ON CUBE BOUNDARY
  template<typename Image_word_type>
  void feature_voxels_on_image_bbox(const CGAL::Image_3& image,
                                    std::vector<std::array<std::size_t, 3>>& black_voxels)
  {
    const std::size_t& xdim = image.xdim();
    const std::size_t& ydim = image.ydim();
    const std::size_t& zdim = image.zdim();

    const std::size_t wx = (std::max)(xdim - 1, std::size_t(1));
    const std::size_t wy = (std::max)(ydim - 1, std::size_t(1));
    const std::size_t wz = (std::max)(zdim - 1, std::size_t(1));

    for (int axis = 0; axis < 3; ++axis)
    {
      for (std::size_t i = 0; i < xdim; i += (axis == 0 ? wx : 1))
        for (std::size_t j = 0; j < ydim; j += (axis == 1 ? wy : 1))
          for (std::size_t k = 0; k < zdim; k += (axis == 2 ? wz : 1))
          {
            typedef std::array<std::size_t, 3> Pixel;

            Pixel pix00 = { {i  , j  , k  } },
              pix10 = pix00, pix01 = pix00, pix11 = pix00;

            const int axis_xx = (axis + 1) % 3;
            const int axis_yy = (axis + 2) % 3;

            ++pix10[axis_xx];
            ++pix11[axis_xx]; ++pix11[axis_yy];
            ++pix01[axis_yy];
            if (pix11[0] >= xdim || pix11[1] >= ydim || pix11[2] >= zdim) {
              // we have gone too far
              continue;
            }

            struct Enriched_pixel {
              Pixel pixel;
              Image_word_type word;
            };

            std::array<std::array<Enriched_pixel, 2>, 2> square =
            { { {{ { pix00, Image_word_type() },
                   { pix01, Image_word_type() } }},
                {{ { pix10, Image_word_type() },
                   { pix11, Image_word_type() } }} } };

            std::map<Image_word_type, int> pixel_values_set;
            for (int ii = 0; ii < 2; ++ii) {
              for (int jj = 0; jj < 2; ++jj) {
                const Pixel& pixel = square[ii][jj].pixel;
                short sum_faces =
                  ((0 == pixel[0] || (xdim - 1) == pixel[0]) ? 1 : 0)
                  + ((0 == pixel[1] || (ydim - 1) == pixel[1]) ? 1 : 0)
                  + ((0 == pixel[2] || (zdim - 1) == pixel[2]) ? 1 : 0);

                square[ii][jj].word = CGAL::IMAGEIO::static_evaluate<Image_word_type>
                  (image.image(), pixel[0], pixel[1], pixel[2]);
                ++pixel_values_set[square[ii][jj].word];

                if (pixel_values_set.size() > 1 || sum_faces > 1/*on edge of bbox*/)
                  black_voxels.push_back({ i, j, k });
              }
            }//end for loops on ii, jj
          }//end for loops on i,j,k
    }//end for loop on axis
  }

}//end namespace internal
}//namespace Mesh_3
}//namespace CGAL

#endif // CGAL_MESH_3_POSTPROCESS_LABEL_WEIGHTS_H
