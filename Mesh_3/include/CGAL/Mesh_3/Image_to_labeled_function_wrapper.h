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



#include <CGAL/Image_3.h>

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
         bool use_trilinear_interpolation=true>
class Image_to_labeled_function_wrapper
{
public:
  // Types
  typedef Return_type return_type;
  typedef Image_word_type word_type;
  typedef typename BGT::Point_3   Point_3;

  /// Constructor
  Image_to_labeled_function_wrapper(const Image_& image)
    : r_im_(image) {}

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
    if ( use_trilinear_interpolation )
    {
      return static_cast<return_type>(
          r_im_.labellized_trilinear_interpolation(
              CGAL::to_double(p.x()),
              CGAL::to_double(p.y()),
              CGAL::to_double(p.z()),
              word_type(0)));
    }
    else
    {
      const int px = static_cast<int>(p.x()/r_im_.vx());
      const int py = static_cast<int>(p.y()/r_im_.vy());
      const int pz = static_cast<int>(p.z()/r_im_.vz());

      const int dimx = r_im_.xdim();
      const int dimy = r_im_.ydim();
      const int dimz = r_im_.zdim();

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
      return data[ pz*dimy*dimx + py*dimx + px ];
    }
  }

private:
  /// Labeled image to wrap
  const Image_& r_im_;

};  // end class Image_to_labeled_function_wrapper


}  // end namespace Mesh_3

}  // end namespace CGAL

#endif // CGAL_MESH_3_IMAGE_TO_LABELED_FUNCTION_WRAPPER_H
