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

#include <itkImage.h>
#include <itkThresholdImageFilter.h>
#include <itkRescaleIntensityImageFilter.h>

//#include <itkMaximumImageFilter.h>
//#include <itkSmoothingRecursiveGaussianImageFilter.h>
//#include <itkConnectedComponentImageFilter.h>

#include <boost/container/flat_set.hpp>

#include <iostream>
#include <vector>

namespace CGAL {
namespace Mesh_3 {
namespace internal {

template<typename Image_word_type>
void convert_image_3_to_itk(const CGAL::Image_3& image,
                            itk::Image<Image_word_type, 3>& itk_img,
                            boost::container::flat_set<Image_word_type>& labels)
{
  using PixelType = Image_word_type;
  using ImageType = itk::Image<PixelType, 3/*Dimension*/>;

  typename ImageType::SpacingType spacing;
  spacing[0] = image.vx();
  spacing[1] = image.vy();
  spacing[2] = image.vz();
  itk_img.SetSpacing(spacing);

  typename ImageType::PointType origin;
  origin[0] = image.tx();
  origin[1] = image.ty();
  origin[2] = image.tz();
  itk_img.SetOrigin(origin);

  using Index = itk::Index<3>::IndexValueType;
  for (std::size_t i = 0; i < image.xdim(); ++i)
  {
    for (std::size_t j = 0; j < image.ydim(); ++j)
    {
      for (std::size_t k = 0; j < image.zdim(); ++k)
      {
        typename ImageType::IndexType index = {(Index)i, (Index)j, (Index)k};
        const Image_word_type& label = image.value(i, j, k);
        itk_img.SetPixel(index, label);
        labels.insert(label);
      }
    }
  }
}

}//namespace internal

template<typename Image_word_type>
void generate_weights(const CGAL::Image_3& image,
                      CGAL::Image_3& weights,
                      const float& sigma,
                      Image_word_type)
{
  using PixelType = Image_word_type;
  using ImageType = itk::Image<PixelType, 3/*Dimension*/>;

  typename ImageType::Pointer itk_img = ImageType::New();
  boost::container::flat_set<Image_word_type> labels;
  internal::convert_image_3_to_itk(image, *itk_img, labels);

  //create indicator images
  using IndicatorFilter = itk::ThresholdImageFilter<ImageType>;
  using RescaleFilterType = itk::RescaleIntensityImageFilter<ImageType, ImageType>;

  std::vector<typename ImageType::Pointer> indicators(labels.size());

  for (Image_word_type label : labels)
  {
    //compute "indicator image" for "label"
    IndicatorFilter::Pointer indicator = IndicatorFilter::New();
    indicator->SetInput(itk_img);
    indicator->SetOutsideValue(0);
    indicator->ThresholdOutside(label, label);
    indicator->Update();

    //rescale it "* 255"
    RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
    rescaler->SetOutputMinimum(0);
    rescaler->SetOutputMaximum(255);
    rescaler->SetInput(indicator->GetOutput());
    rescaler->Update();

    //save for later use
    indicators.push_back(rescaler->GetOutput());
  }

  //  using ConnectedComponentImageFilterType = itk::ConnectedComponentImageFilter<ImageType, ImageType>;
  //
  //  ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New();
  //  //connected->SetInput(reader->GetOutput());
  //  connected->Update();
  //
  //  std::size_t nb_labels = connected->GetObjectCount();
  //  std::cout << "Number of objects: " << nb_labels << std::endl;
  //
  //  std::vector<ImageType> indicators(nb_labels);
  //
  //  //smooth indicator functions
  //  using FilterType = itk::SmoothingRecursiveGaussianImageFilter<ImageType, ImageType>;
  //  FilterType::Pointer smoothFilter = FilterType::New();
  //
  //  smoothFilter->SetSigma(sigma);
  ////  smoothFilter->SetInput(reader->GetOutput());
  ////
  //////  writer->SetInput(smoothFilter->GetOutput());
  ////
  //
//    //take the max of indicator functions
//  using MaximumImageFilterType = itk::MaximumImageFilter<ImageType>;
//  MaximumImageFilterType::Pointer maximumImageFilter = MaximumImageFilterType::New();
//  //  for(itk::ProcessObject::DataObjectIdentifierType i = 0; i < connected->GetObjectCount(); ++i)
//  //    maximumImageFilter->SetInput(i, indicators[i]);
//  //  maximumImageFilter->Execute();
//  //  maximumImageFilter->Update();
//
//    //output
//  weights.set_data(maximumImageFilter->GetOutput());
}

}//namespace Mesh_3
}//namespace CGAL

#endif // CGAL_MESH_3_GENERATE_WEIGHTS_FROM_LABELED_IMAGE_H
