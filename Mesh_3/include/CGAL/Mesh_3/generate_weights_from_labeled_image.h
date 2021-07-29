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

#ifndef CGAL_MESH_3_GENERATE_WEIGHTS_FROM_LABELED_IMAGE_H
#define CGAL_MESH_3_GENERATE_WEIGHTS_FROM_LABELED_IMAGE_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/Image_3.h>
#include <CGAL/ImageIO.h>

#include <itkImage.h>
#include <itkImageDuplicator.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkRecursiveGaussianImageFilter.h>
#include <itkMaximumImageFilter.h>

#include <iostream>
#include <vector>
#include <set>
#include <type_traits>

namespace CGAL {
namespace Mesh_3 {
namespace internal {

template<typename Image_word_type, typename LabelsSet>
void convert_image_3_to_itk(const CGAL::Image_3& image,
                            itk::Image<Image_word_type, 3>* const itk_img,
                            LabelsSet& labels)
{
  using ImageType = itk::Image<Image_word_type, 3/*Dimension*/>;

  typename ImageType::SpacingType spacing;
  spacing[0] = image.vx();
  spacing[1] = image.vy();
  spacing[2] = image.vz();
  itk_img->SetSpacing(spacing);

  typename ImageType::PointType origin;
  origin[0] = image.tx();
  origin[1] = image.ty();
  origin[2] = image.tz();
  itk_img->SetOrigin(origin);

  typename ImageType::IndexType  corner = {{0, 0, 0 }};
  typename ImageType::SizeType   size = {{image.xdim(), image.ydim(), image.zdim()}};
  typename ImageType::RegionType region(corner, size);
  itk_img->SetRegions(region);

//  itk_img->Allocate(image.xdim() * image.ydim() * image.zdim());
  itk_img->Allocate();

  using Index = itk::Index<3>::IndexValueType;
  for (std::size_t i = 0; i < image.xdim(); ++i)
  {
    for (std::size_t j = 0; j < image.ydim(); ++j)
    {
      for (std::size_t k = 0; k < image.zdim(); ++k)
      {
        typename ImageType::IndexType index = {(Index)i, (Index)j, (Index)k};
        const Image_word_type label = image.value(i, j, k);
        labels.insert(label);
        itk_img->SetPixel(index, label);
      }
    }
  }
}

#ifdef CGAL_MESH_3_WEIGHTED_IMAGES_DEBUG
int count_non_white_pixels(const CGAL::Image_3& image)
{
  int nb_nonzero = 0;
  for (std::size_t i = 0; i < image.xdim(); ++i)
  {
    for (std::size_t j = 0; j < image.ydim(); ++j)
    {
      for (std::size_t k = 0; k < image.zdim(); ++k)
      {
        if (image.value(i, j, k) != 255)
          nb_nonzero++;
      }
    }
  }
  return nb_nonzero;
}

template<typename Image_word_type>
int count_non_white_pixels(const itk::Image<Image_word_type, 3>* image)
{
  int nb_nonzero = 0;
  const auto sizeOfImage = image->GetLargestPossibleRegion().GetSize();
  using Index = itk::Index<3>::IndexValueType;

  for (std::size_t i = 0; i < sizeOfImage[0]; ++i)
  {
    for (std::size_t j = 0; j < sizeOfImage[1]; ++j)
    {
      for (std::size_t k = 0; k < sizeOfImage[2]; ++k)
      {
        const itk::Index<3>  index = { (Index)i, (Index)j, (Index)k };
        if (image->GetPixel(index) != 255)
          nb_nonzero++;
      }
    }
  }
  return nb_nonzero;
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

}//namespace internal

/// @cond INTERNAL
template<typename Image_word_type>
CGAL::Image_3 generate_weights_with_known_word_type(const CGAL::Image_3& image,
                                                    const float& sigma)
{
  typedef unsigned char Weights_type; //from 0 t 255

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
            weights_ptr + image.xdim() * image.ydim() * image.zdim(),
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

  using DuplicatorType = itk::ImageDuplicator<ImageType>;
  using IndicatorFilter = itk::BinaryThresholdImageFilter<ImageType, WeightsType>;
  using GaussianFilterType = itk::RecursiveGaussianImageFilter<WeightsType, WeightsType>;
  using MaximumImageFilterType = itk::MaximumImageFilter<WeightsType>;

  std::vector<typename ImageType::Pointer> indicators(labels.size());
  typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage(itk_img);
  duplicator->Update();

  int id = 0;
  for (Image_word_type label : labels)
  {
    if (id > 0)
    {
      duplicator->SetInputImage(indicators[id - 1]);
      duplicator->Update();
    }
    indicators[id++] = duplicator->GetOutput();
  }

  id = 0;
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

    //perform gaussian smoothing
    typename GaussianFilterType::Pointer smoother = GaussianFilterType::New();
    smoother->SetInput(indicator->GetOutput());
    smoother->SetSigma(sigma);
    smoother->Update();

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

#ifdef CGAL_MESH_3_WEIGHTED_IMAGES_DEBUG
    std::cout << "AFTER MAX (label = " << label << ") : " <<  std::endl;
    std::cout << "\tnon zero in max ("
      << label << ")\t= " << internal::count_non_white_pixels(blured_max.GetPointer()) << std::endl;
#endif
  }

  //copy pixels to weights
  using Index = itk::Index<3>::IndexValueType;
  for (std::size_t i = 0; i < image.xdim(); ++i)
  {
    for (std::size_t j = 0; j < image.ydim(); ++j)
    {
      for (std::size_t k = 0; k < image.zdim(); ++k)
      {
        typename ImageType::IndexType index = { (Index)i, (Index)j, (Index)k };
        using CGAL::IMAGEIO::static_evaluate;
        static_evaluate<Weights_type>(weights, i, j, k) = blured_max->GetPixel(index);
      }
    }
  }

#ifdef CGAL_MESH_3_WEIGHTED_IMAGES_DEBUG
  std::cout << "non zero in image \t= " << internal::count_non_white_pixels(image) << std::endl;
  std::cout << "non zero in weights \t= " << internal::count_non_white_pixels(blured_max.GetPointer()) << std::endl;
#endif

  _writeImage(weights, "weights-image.inr.gz");
  return CGAL::Image_3(weights);
}
/// @endcond

/*!
* Free function that generates a `CGAL::Image_3` of weights associated to each
* voxel of `image`, to make the output mesh surfaces smoother.
* The weights image is generated using the algorithm described by Stalling et al
* in \cgalCite{stalling1998weighted}.
*
* @param image the input labeled image from which the weights image is computed.
*   Both will then be used to construct a `Labeled_mesh_domain_3`.
* @param sigma the standard deviation parameter of the internal Gaussian filter
*
* @returns a `CGAL::Image_3` of weights used to build a quality `Labeled_mesh_domain_3`
*/

CGAL::Image_3 generate_weights(const CGAL::Image_3& image,
                               const float& sigma)
{
  CGAL_IMAGE_IO_CASE(image.image(),
    return generate_weights_with_known_word_type<Word>(image, sigma);
  );
  CGAL_error_msg("This place should never be reached, because it would mean "
    "the image word type is a type that is not handled by "
    "CGAL_ImageIO.");
  return CGAL::Image_3();
}

}//namespace Mesh_3
}//namespace CGAL

#endif // CGAL_MESH_3_GENERATE_WEIGHTS_FROM_LABELED_IMAGE_H
