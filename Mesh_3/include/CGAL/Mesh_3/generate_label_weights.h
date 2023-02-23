// Copyright (c) 2021 GeometryFactory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent Rineau, Jane Tournois

#ifndef CGAL_MESH_3_GENERATE_LABEL_WEIGHTS_H
#define CGAL_MESH_3_GENERATE_LABEL_WEIGHTS_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/Image_3.h>
#include <CGAL/ImageIO.h>

#include <itkImage.h>
#include <itkImageDuplicator.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkDiscreteGaussianImageFilter.h>
#include <itkMaximumImageFilter.h>

#include <CGAL/Mesh_3/features_detection/features_detection.h>
#include <CGAL/Mesh_3/features_detection/coordinates.h>
#include <CGAL/Mesh_3/features_detection/combinations.h>
#include <CGAL/Mesh_3/features_detection/cases_table.h>
#include <CGAL/Mesh_3/features_detection/cube_isometries.h>
#include <CGAL/Mesh_3/features_detection/features_detection_helpers.h>

#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <iostream>
#include <vector>
#include <set>
#include <type_traits>
#include <algorithm>

namespace CGAL {
namespace Mesh_3 {
namespace internal {

template<typename Image_word_type, typename LabelsSet>
void convert_image_3_to_itk(const CGAL::Image_3& image,
                            itk::Image<Image_word_type, 3>* const itk_img,
                            LabelsSet& labels)
{
  const double spacing[3] = {image.vx(), image.vy(), image.vz()};
  itk_img->SetSpacing(spacing);

  const double origin[3] =  {image.tx(), image.ty(), image.tz()};
  itk_img->SetOrigin(origin);

  using ImageType = itk::Image<Image_word_type, 3/*Dimension*/>;
  typename ImageType::IndexType  corner = {{0, 0, 0 }};
  typename ImageType::SizeType   size = {{image.xdim(), image.ydim(), image.zdim()}};
  typename ImageType::RegionType region(corner, size);
  itk_img->SetRegions(region);

  itk_img->Allocate();

  const Image_word_type* img_begin = static_cast<const Image_word_type*>(image.data());
  std::copy(img_begin, img_begin + image.size(), itk_img->GetBufferPointer());

  labels.insert(img_begin, img_begin + image.size());
}

#ifdef CGAL_MESH_3_WEIGHTED_IMAGES_DEBUG
template<typename Image_word_type>
int count_non_white_pixels(const CGAL::Image_3& image)
{
  auto diff255 = [&](const Image_word_type p)
  {
    return p != 255;
  };
  const Image_word_type* img_begin = static_cast<const Image_word_type*>(image.data());
  return std::count_if(img_begin,
                       img_begin + image.size(),
                       diff255);
}

template<typename Image_word_type>
int count_non_white_pixels(itk::Image<Image_word_type, 3>* itk_img)
{
  auto diff255 = [&](const Image_word_type p)
  {
    return p != 255;
  };
  auto size = itk_img->GetLargestPossibleRegion().GetSize();
  return std::count_if(itk_img->GetBufferPointer(),
                       itk_img->GetBufferPointer() + size[0]*size[1]*size[2],
                       diff255);
}
#endif //CGAL_MESH_3_WEIGHTED_IMAGES_DEBUG

template<typename Image_word_type>
WORD_KIND get_wordkind()
{
  if (std::is_floating_point<Image_word_type>::value)
    return WK_FLOAT;
  else
    return WK_FIXED;
/** unknown (uninitialized) */
//    WK_UNKNOWN
}

template<typename Image_word_type>
SIGN get_sign()
{
  if (std::is_signed<Image_word_type>::value)
    return SGN_SIGNED;
  else
    return SGN_UNSIGNED;
/** unknown (uninitialized or floating point words) */
//    SGN_UNKNOWN
}

#ifdef CGAL_MESH_3_WEIGHTED_IMAGES_DEBUG
template<typename Image_word_type>
void convert_itk_to_image_3(itk::Image<Image_word_type, 3>* const itk_img,
                            const char* filename = "")
{
  auto t = itk_img->GetOrigin();
  auto v = itk_img->GetSpacing();
  auto region = itk_img->GetRequestedRegion();

  _image* img
    = _createImage(region.GetSize(0), region.GetSize(1), region.GetSize(2),
      1,                                        //vectorial dimension
      v[0], v[1], v[2],
      sizeof(Image_word_type),                     //image word size in bytes
      internal::get_wordkind<Image_word_type>(),   //image word kind WK_FIXED, WK_FLOAT, WK_UNKNOWN
      internal::get_sign<Image_word_type>());      //image word sign
  Image_word_type* img_ptr = (Image_word_type*)(img->data);

  const int size = region.GetSize(0) * region.GetSize(1) * region.GetSize(2);
  std::fill(img_ptr,
            img_ptr + size,
            Image_word_type(0));
  img->tx = t[0];
  img->ty = t[1];
  img->tz = t[2];

  std::copy(itk_img->GetBufferPointer(),
            itk_img->GetBufferPointer() + size,
            img_ptr);

  if(filename != "")
    _writeImage(img, filename);
}
#endif

}//namespace internal

/// @cond INTERNAL
template<typename Image_word_type>
CGAL::Image_3 generate_label_weights_with_known_word_type(const CGAL::Image_3& image,
                                                          const float& sigma,
                                                          const bool with_features)
{
  typedef unsigned char Weights_type; //from 0 t 255
  const std::size_t img_size = image.size();

  //create weights image
  _image* weights
    = _createImage(image.xdim(), image.ydim(), image.zdim(),
                   1,                                        //vectorial dimension
                   image.vx(), image.vy(), image.vz(),
                   sizeof(Weights_type),                     //image word size in bytes
                   internal::get_wordkind<Weights_type>(),   //image word kind WK_FIXED, WK_FLOAT, WK_UNKNOWN
                   internal::get_sign<Weights_type>());      //image word sign
  Weights_type* weights_ptr = (Weights_type*)(weights->data);
  std::fill(weights_ptr,
            weights_ptr + img_size,
            Weights_type(0));
  weights->tx = image.tx();
  weights->ty = image.ty();
  weights->tz = image.tz();

  //convert image to itkImage
  using ImageType = itk::Image<Image_word_type, 3/*Dimension*/>;
  using WeightsType = itk::Image<Weights_type, 3>;
  typename ImageType::Pointer itk_img = ImageType::New();
  std::set<Image_word_type> labels;
  internal::convert_image_3_to_itk(image, itk_img.GetPointer(), labels);

#ifdef CGAL_MESH_3_WEIGHTED_IMAGES_DEBUG
  CGAL_assertion(internal::count_non_white_pixels<Image_word_type>(image)
              == internal::count_non_white_pixels<Image_word_type>(itk_img.GetPointer()));
#endif

  using DuplicatorType = itk::ImageDuplicator<ImageType>;
  using IndicatorFilter = itk::BinaryThresholdImageFilter<ImageType, WeightsType>;
  using GaussianFilterType = itk::DiscreteGaussianImageFilter<WeightsType, WeightsType>;
  using MaximumImageFilterType = itk::MaximumImageFilter<WeightsType>;

  std::vector<typename ImageType::Pointer> indicators(labels.size());
  typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage(itk_img);
  duplicator->Update();

  for (std::size_t id = 0; id < labels.size(); ++id)
  {
    if (id > 0)
    {
      duplicator->SetInputImage(indicators[id - 1]);
      duplicator->Update();
    }
    indicators[id] = duplicator->GetOutput();
  }

  int id = 0;
  typename WeightsType::Pointer blured_max = WeightsType::New();
  for (Image_word_type label : labels)
  {
#ifdef CGAL_MESH_3_WEIGHTED_IMAGES_DEBUG
    std::cout << "\nLABEL = " << label << std::endl;
#endif

    //compute "indicator image" for "label"
    typename IndicatorFilter::Pointer indicator = IndicatorFilter::New();
    indicator->SetInput(indicators[id]);
    indicator->SetOutsideValue(0);
    indicator->SetInsideValue(255);
    indicator->SetLowerThreshold(label);
    indicator->SetUpperThreshold(label);
    indicator->Update();

#ifdef CGAL_MESH_3_WEIGHTED_IMAGES_DEBUG
    std::ostringstream oss;
    oss << "indicator_" << id << ".inr.gz";
    std::cout << "filename = " << oss.str().c_str() << std::endl;
    internal::convert_itk_to_image_3(indicator->GetOutput(), oss.str().c_str());
#endif

    //perform gaussian smoothing
    typename GaussianFilterType::Pointer smoother = GaussianFilterType::New();
    smoother->SetUseImageSpacing(true);//variance/std deviation is counted real world distances
    smoother->SetInput(indicator->GetOutput());
    smoother->SetVariance(sigma*sigma);
    smoother->Update();

#ifdef CGAL_MESH_3_WEIGHTED_IMAGES_DEBUG
    std::ostringstream oss1;
    oss1 << "smooth_" << id << ".inr.gz";
    std::cout << "filename = " << oss1.str().c_str() << std::endl;
    internal::convert_itk_to_image_3(smoother->GetOutput(), oss1.str().c_str());
#endif

    //take the max of smoothed indicator functions
    if (id == 0)
      blured_max = smoother->GetOutput();
    else
    {
      typename MaximumImageFilterType::Pointer maximumImageFilter = MaximumImageFilterType::New();
      maximumImageFilter->SetInput(0, blured_max);
      maximumImageFilter->SetInput(1, smoother->GetOutput());
      maximumImageFilter->Update();
      blured_max = maximumImageFilter->GetOutput();
    }

    id++;
  }


#ifdef CGAL_MESH_3_WEIGHTED_IMAGES_DEBUG
    std::ostringstream oss2;
    oss2 << "max_" << "all" << ".inr.gz";
    std::cout << "filename = " << oss2.str().c_str() << std::endl;
    internal::convert_itk_to_image_3(blured_max.GetPointer(), oss2.str().c_str());
#endif

#ifdef CGAL_MESH_3_WEIGHTED_IMAGES_DEBUG
//    std::cout << "AFTER MAX (label = " << label << ") : " <<  std::endl;
    std::cout << "\tnon zero in max ("
      << id << ")\t= " << internal::count_non_white_pixels(blured_max.GetPointer()) << std::endl;
#endif

  //copy pixels to weights
  std::copy(blured_max->GetBufferPointer(),
            blured_max->GetBufferPointer() + img_size,
            weights_ptr);

  CGAL::Image_3 weights_img(weights);

#ifdef CGAL_MESH_3_WEIGHTED_IMAGES_DEBUG
  std::cout << "non white in image \t= "
    << internal::count_non_white_pixels<Image_word_type>(image) << std::endl;
  std::cout << "non white in weights \t= "
    << internal::count_non_white_pixels<Weights_type>(weights_img) << std::endl;
  std::cout << "non white in itkWeights \t= "
    << internal::count_non_white_pixels<Weights_type>(blured_max.GetPointer()) << std::endl;
  _writeImage(weights, "weights-image.inr.gz");
#endif

  if (with_features)
  {
    postprocess_weights_for_feature_protection(image, weights_img);
  }

  return weights_img;
}
/// @endcond

/*!
* \ingroup PkgMesh3Functions
* Free function that generates a `CGAL::Image_3` of weights associated to each
* voxel of `image`, to make the output mesh surfaces smoother.
* The weights image is generated using the algorithm described by Stalling et al
* in \cgalCite{stalling1998weighted}.
* The [Insight toolkit](https://itk.org/) is needed to compile this function.
*
* @param image the input labeled image from which the weights image is computed.
*   Both will then be used to construct a `Labeled_mesh_domain_3`.
* @param sigma the standard deviation parameter of the internal Gaussian filter,
*   measured in real-world distances. The size of a voxel (e.g. shortest length
*   or longest length) usually is a good value for this parameter.
*   Note that if `sigma` is too small, the "stair-effect" of meshing from
*   a voxel image can appear. On the other side, if `sigma` is too large,
*   thin volumes (basically one voxel thick) may be lost in the meshing process
*   because the computed weights are too blurry.
*
* @returns a `CGAL::Image_3` of weights used to build a quality `Labeled_mesh_domain_3`,
* with the same dimensions as `image`
*/
template<typename CGAL_NP_TEMPLATE_PARAMETERS>
CGAL::Image_3 generate_label_weights(const CGAL::Image_3& image,
    const float& sigma,
    const CGAL_NP_CLASS& np = parameters::default_values())
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  const bool with_features = choose_parameter(get_parameter(np, internal_np::with_features_param), false);
  CGAL_IMAGE_IO_CASE(image.image(),
    return generate_label_weights_with_known_word_type<Word>(image, sigma, with_features);
  );
  CGAL_error_msg("This place should never be reached, because it would mean "
    "the image word type is a type that is not handled by "
    "CGAL_ImageIO.");
  return CGAL::Image_3();
}

namespace internal
{
  template<typename Word_type>
  void set_voxel(CGAL::Image_3& img,
                 const std::size_t& i,
                 const std::size_t& j,
                 const std::size_t& k,
                 const Word_type& w)
  {
    using CGAL::IMAGEIO::static_evaluate;

    if (i < 0 || j < 0 || k < 0)
      return;
    else if (i > img.xdim() - 1 || j > img.ydim() - 1 || k > img.zdim() - 1)
      return;
    else
      static_evaluate<Word_type>(img.image(), i, j, k) = w;
  }
}
inline
void postprocess_weights_for_feature_protection(const CGAL::Image_3& image,
                                                CGAL::Image_3& weights)
{
  typedef unsigned char Image_word_type; //todo
  typedef unsigned char Weights_type; //from 0 t 255

  using CGAL::IMAGEIO::static_evaluate;

  using Word //use unsigned integral Word type to use it as an index
    = typename CGAL::IMAGEIO::Word_type_generator<WK_FIXED, SGN_UNSIGNED, sizeof(Image_word_type)>::type;

  using Color_transform = internal::Color_transformation_helper<Word>;
  typename Color_transform::type color_transformation;

  struct Dummy_point {
    double x, y, z;
    Dummy_point(double, double, double) {}
  };
  CGAL::Mesh_3::Triple_line_extractor<Dummy_point> lines;

  const std::size_t xdim = image.xdim();
  const std::size_t ydim = image.ydim();
  const std::size_t zdim = image.zdim();

  std::vector<std::array<std::size_t, 3>> black_voxels;

  // POLYLINES INSIDE THE CUBE
  for (std::size_t k = 0, end_k = zdim - 1; k < end_k; ++k)
  {
    for (std::size_t j = 0, end_j = ydim - 1; j < end_j; ++j)
    {
      for (std::size_t i = 0, end_i = xdim - 1; i < end_i; ++i)
      {
        using Cube = CGAL::Mesh_3::Cube;
        const Cube cube = {
          static_evaluate<Word>(image.image(), i  , j  , k),
          static_evaluate<Word>(image.image(), i + 1, j  , k),
          static_evaluate<Word>(image.image(), i  , j + 1, k),
          static_evaluate<Word>(image.image(), i + 1, j + 1, k),
          static_evaluate<Word>(image.image(), i  , j  , k + 1),
          static_evaluate<Word>(image.image(), i + 1, j  , k + 1),
          static_evaluate<Word>(image.image(), i  , j + 1, k + 1),
          static_evaluate<Word>(image.image(), i + 1, j + 1, k + 1),
        };

        bool monocolor = (cube[0] == cube[1]);
        for (int i = 2; i < 8; ++i) monocolor = monocolor && (cube[0] == cube[i]);
        if (monocolor) continue;

        Color_transform::reset(color_transformation);

        std::uint8_t nb_color = 0;
        for (int i = 0; i < 8; ++i) {
          if (!Color_transform::is_valid(color_transformation, cube[i]))
          {
            color_transformation[cube[i]] = nb_color;
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
        if (!case_found)
          CGAL_error();
        else
          reference_cube = internal::combinations[(*case_it)[8]];

        auto fct_it = lines.create_polylines_fcts.find(reference_cube);
        if (fct_it != lines.create_polylines_fcts.end())
          black_voxels.push_back({ i, j, k });
      }
    }
  }

  // POLYLINES ON CUBE BOUNDARY
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

  Weights_type wblack(0);
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
  std::cout << "non white in post-processed image \t= "
    << internal::count_non_white_pixels<Weights_type>(weights) << std::endl;
  _writeImage(weights.image(), "weights-image_postprocessed.inr.gz");
#endif
}

}//namespace Mesh_3
}//namespace CGAL

#endif // CGAL_MESH_3_GENERATE_LABEL_WEIGHTS_H
