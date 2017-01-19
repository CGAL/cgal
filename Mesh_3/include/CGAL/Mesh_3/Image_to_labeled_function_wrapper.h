// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
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


#include <CGAL/Image_3.h>
#include <CGAL/function_objects.h>


namespace CGAL {

namespace Mesh_3 {

/**
 * @class Image_to_labeled_function_wrapper
 *
 * Wraps a labeled image into a labeled function which takes his values into
 * N. Uses trilinear interpolation.
 * Note: Image has to be labeled with unsigned char
 */
template<class Image_,
         class BGT,
         typename Image_word_type = unsigned char,
         typename Return_type = int,
         typename Transform = Identity<Return_type>,
         bool labeled_image = true,
         bool use_trilinear_interpolation=true>
class Image_to_labeled_function_wrapper
{
public:
  // Types
  typedef Return_type return_type;
  typedef Image_word_type word_type;
  typedef typename BGT::Point_3   Point_3;

  /// Constructor
  Image_to_labeled_function_wrapper(const Image_& image, 
                                    const Transform& transform = Transform(),
                                    const Return_type value_outside = 0)
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
  return_type operator()(const Point_3& p, const bool = true) const
  {
    return eval(p,
                CGAL::Boolean_tag<use_trilinear_interpolation>(),
                CGAL::Boolean_tag<labeled_image>());
  }

private:
  return_type eval(const Point_3& p,
                   CGAL::Tag_true /*trilinear*/,
                   CGAL::Tag_true /*labeled*/) const
  {
    return static_cast<return_type>(transform(
      r_im_.template labellized_trilinear_interpolation<Image_word_type>(
          CGAL::to_double(p.x()),
          CGAL::to_double(p.y()),
          CGAL::to_double(p.z()),
          value_outside)));
  }

  return_type eval(const Point_3& p,
                   CGAL::Tag_true /*trilinear*/,
                   CGAL::Tag_false /*labeled*/) const
  {
    return transform(
        r_im_.template trilinear_interpolation<Image_word_type, double>(
          CGAL::to_double(p.x()),
          CGAL::to_double(p.y()),
          CGAL::to_double(p.z()),
          value_outside));
  }

  template <typename Labeled_tag>
  return_type eval(const Point_3& p,
                   CGAL::Tag_false /*trilinear*/,
                   Labeled_tag /*labeled*/) const
  {
    const std::ptrdiff_t px = static_cast<std::ptrdiff_t>(p.x()/r_im_.vx());
    const std::ptrdiff_t py = static_cast<std::ptrdiff_t>(p.y()/r_im_.vy());
    const std::ptrdiff_t pz = static_cast<std::ptrdiff_t>(p.z()/r_im_.vz());

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
  const Transform transform;
  const Return_type value_outside;

};  // end class Image_to_labeled_function_wrapper


}  // end namespace Mesh_3

}  // end namespace CGAL

#endif // CGAL_MESH_3_IMAGE_TO_LABELED_FUNCTION_WRAPPER_H
