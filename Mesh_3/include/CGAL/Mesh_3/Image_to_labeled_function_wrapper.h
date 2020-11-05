// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Stéphane Tayeb
//
//******************************************************************************
// File Description :
// Image_to_labeled_function_wrapper declaration and implementation. See
// class description.
//******************************************************************************

#ifndef CGAL_MESH_3_IMAGE_TO_LABELED_FUNCTION_WRAPPER_H
#define CGAL_MESH_3_IMAGE_TO_LABELED_FUNCTION_WRAPPER_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Image_3.h>
#include <CGAL/function_objects.h>
#include <functional>
#include <boost/mpl/if.hpp>

namespace CGAL {

namespace Mesh_3 {

/**
 * @class Image_to_labeled_function_wrapper
 *
 * Wraps a labeled image into a labeled function which takes his values into
 * N. Uses trilinear interpolation.
 * Note: Image has to be labeled with unsigned char
 */
template<typename Image_word_type = unsigned char,
         typename Interpolation_type = Image_word_type,
         typename Return_type = int,
         bool labeled_image = true,
         bool use_trilinear_interpolation=true>
class Image_to_labeled_function_wrapper
{
public:
  typedef std::function<Return_type(Interpolation_type)>
                                                    Image_values_to_labels;

  // Types
  typedef Return_type return_type;
  typedef Image_word_type word_type;
  typedef CGAL::Image_3 Image_;

  /// Constructor
  Image_to_labeled_function_wrapper
  (const Image_& image,
   Image_values_to_labels transform = Identity<Return_type>(),
   const Interpolation_type value_outside = Interpolation_type())
    : r_im_(image)
    , transform(transform)
    , value_outside(value_outside)
  {
  }

  // Default copy constructor and assignment operator are ok

  /// Destructor
  ~Image_to_labeled_function_wrapper() {}

  /**
   * Returns an int corresponding to the label at point \c p
   * @param p the input point
   * @return the label at point \c p
   */
  template <typename Point_3>
  return_type operator()(const Point_3& p) const
  {
    return eval(p,
                CGAL::Boolean_tag<use_trilinear_interpolation>(),
                CGAL::Boolean_tag<labeled_image>());
  }

private:
  template <typename Point_3>
  return_type eval(const Point_3& p,
                   CGAL::Tag_true /*trilinear*/,
                   CGAL::Tag_true /*labeled*/) const
  {
    return static_cast<return_type>(transform(
       r_im_.template labellized_trilinear_interpolation<Image_word_type>(
          CGAL::to_double(p.x()-r_im_.image()->tx),
          CGAL::to_double(p.y()-r_im_.image()->ty),
          CGAL::to_double(p.z()-r_im_.image()->tz),
          value_outside)));
  }

  template <typename Point_3>
  return_type eval(const Point_3& p,
                   CGAL::Tag_true /*trilinear*/,
                   CGAL::Tag_false /*labeled*/) const
  {
    return transform(
        r_im_.template trilinear_interpolation<Image_word_type,
                                               Interpolation_type>(
          CGAL::to_double(p.x()-r_im_.image()->tx),
          CGAL::to_double(p.y()-r_im_.image()->ty),
          CGAL::to_double(p.z()-r_im_.image()->tz),
          value_outside));
  }

  template <typename Labeled_tag, typename Point_3>
  return_type eval(const Point_3& p,
                   CGAL::Tag_false /*trilinear*/,
                   Labeled_tag /*labeled*/) const
  {
    const std::ptrdiff_t px = static_cast<std::ptrdiff_t>((p.x()-r_im_.image()->tx)/r_im_.vx());
    const std::ptrdiff_t py = static_cast<std::ptrdiff_t>((p.y()-r_im_.image()->ty)/r_im_.vy());
    const std::ptrdiff_t pz = static_cast<std::ptrdiff_t>((p.z()-r_im_.image()->tz)/r_im_.vz());

    const std::ptrdiff_t dimx = static_cast<std::ptrdiff_t>(r_im_.xdim());
    const std::ptrdiff_t dimy = static_cast<std::ptrdiff_t>(r_im_.ydim());
    const std::ptrdiff_t dimz = static_cast<std::ptrdiff_t>(r_im_.zdim());

    if(px < 0 ||
       py < 0 ||
       pz < 0 ||
       px+1 >= dimx ||
       py+1 >= dimy ||
       pz+1 >= dimz)
    {
      return 0;
    }

    const word_type* data = static_cast<const word_type*>(r_im_.data());
    return static_cast<return_type>(transform(
              data[pz*dimy*dimx + py*dimx + px]));
  }

  /// Labeled image to wrap
  const Image_& r_im_;
  const Image_values_to_labels transform;
  const Interpolation_type value_outside;

};  // end class Image_to_labeled_function_wrapper


}  // end namespace Mesh_3

}  // end namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_MESH_3_IMAGE_TO_LABELED_FUNCTION_WRAPPER_H
