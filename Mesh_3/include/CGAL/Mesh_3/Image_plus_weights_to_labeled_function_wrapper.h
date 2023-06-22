// Copyright (c) 2016,2021  GeometryFactory Sarl (France).
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
//

#ifndef CGAL_MESH_3_IMAGE_PLUS_WEIGHTS_TO_LABELED_FUNCTION_WRAPPER_H
#define CGAL_MESH_3_IMAGE_PLUS_WEIGHTS_TO_LABELED_FUNCTION_WRAPPER_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/function_objects.h>
#include <CGAL/Image_3.h>

namespace CGAL {

namespace ImageIO {

template <typename Image_word_type,
          typename Weights_type>
class Weighted_indicator_factory
{
  const Image_3& image;
  const Image_3& weights;

public:

  Weighted_indicator_factory(const Image_3& image,
                             const Image_3& weights)
    : image(image)
    , weights(weights)
  {}

  class Indicator
  {
    const Weighted_indicator_factory& f;
    const Image_word_type label;
  public:
    Indicator(const Weighted_indicator_factory& f,
              const Image_word_type label) : f(f), label(label) {}

    double operator()(const Image_word_type& x) const
    {
      const std::ptrdiff_t offset = &x - (Image_word_type *)f.image.data();
      const Weights_type w = (std::max)(
          Weights_type(128), ((Weights_type *)f.weights.data())[offset]);
      return (x == label) ? w : (255 - w);
    }
  }; // end nested class Indicator

  Indicator indicator(const Image_word_type label) const
  {
    return Indicator(*this, label);
  }
}; // end template class ImageIO::Weighted_indicator_factory

} // end namespace CGAL::ImageIO

namespace Mesh_3 {


/**
 * @class Image_plus_weights_to_labeled_function_wrapper
 *
 * Wraps a pair of images into a labeled function which takes his values into
 * N. Uses weighted trilinear interpolation.
 */
template<typename Image_word_type = unsigned char,
         typename Interpolation_type = Image_word_type,
         typename Weights_type = unsigned char,
         typename Return_type = int>
class Image_plus_weights_to_labeled_function_wrapper
{
public:
  typedef std::function<Return_type(Interpolation_type)>
                                                    Image_values_to_labels;

  // Types
  typedef Return_type return_type;
  typedef Image_word_type word_type;
  typedef CGAL::Image_3 Image_;

  /// Constructor
  Image_plus_weights_to_labeled_function_wrapper
  (const Image_& image,
   const Image_& weights_image,
   Image_values_to_labels transform = Identity<Return_type>(),
   const Interpolation_type value_outside = Interpolation_type())
    : r_im_(image)
    , r_weights_im_(weights_image)
    , transform(transform)
    , value_outside(value_outside)
    , indicator_factory(image, weights_image)
  {
  }

  // Default copy constructor and assignment operator are ok

  /// Destructor
  ~Image_plus_weights_to_labeled_function_wrapper() {}

  /**
   * Returns an int corresponding to the label at point `p`.
   * @param p the input point
   * @return the label at point `p`
   */
  template <typename Point_3>
  return_type operator()(const Point_3& p) const
  {
    return static_cast<return_type>(transform(
      r_im_.template labellized_trilinear_interpolation<Image_word_type>(
          CGAL::to_double(p.x() - r_im_.image()->tx),
          CGAL::to_double(p.y() - r_im_.image()->ty),
          CGAL::to_double(p.z() - r_im_.image()->tz),
          value_outside,
          indicator_factory)));
  }

private:
  /// Labeled image to wrap
  const Image_& r_im_;
  const Image_& r_weights_im_;
  const Image_values_to_labels transform;
  const Interpolation_type value_outside;
  CGAL::ImageIO::Weighted_indicator_factory<Image_word_type,
                                            Weights_type> indicator_factory;
};  // end class Image_plus_weights_to_labeled_function_wrapper


}  // end namespace Mesh_3

}  // end namespace CGAL

#endif // CGAL_MESH_3_IMAGE_PLUS_WEIGHTS_TO_LABELED_FUNCTION_WRAPPER_H
