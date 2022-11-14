// Copyright (c) 2005-2007  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Laurent RINEAU

#ifndef CGAL_MESH_3_GRAY_LEVEL_IMAGE_3_H
#define CGAL_MESH_3_GRAY_LEVEL_IMAGE_3_H

#include <CGAL/license/Surface_mesher.h>


#include <CGAL/basic.h>

#include <CGAL/auto_link/ImageIO.h>

#include <memory>
#include <string>

#ifdef CGAL_SURFACE_MESHER_DEBUG_GRAY_LEVEL_IMAGE_3_CONSTRUCTOR
#include <boost/format.hpp>
#endif

#include <CGAL/Image_3.h>

namespace CGAL {

template <typename FT, typename Point>
class Gray_level_image_3 : public Image_3
{
  float isovalue;
  bool positive_inside;
  float value_outside;
public:
  Gray_level_image_3(const Image_3& image, float isoval)
    : Image_3(image),
      isovalue(isoval),
      positive_inside(true),
      value_outside(0.f)
  {
  }

  Gray_level_image_3(const std::string& file, float isoval, bool positive_inside_=true, float value_outside = 0.f)
    : Image_3(),
      isovalue(isoval),
      positive_inside(positive_inside_),
      value_outside(value_outside)
  {
#ifdef CGAL_SURFACE_MESHER_DEBUG_GRAY_LEVEL_IMAGE_3_CONSTRUCTOR
    std::cerr <<
      ::boost::format("Constructing a Gray_level_image_3(\"%1%\")... ") % file;
#endif
    Image_3::read(file);
#ifdef CGAL_SURFACE_MESHER_DEBUG_GRAY_LEVEL_IMAGE_3_CONSTRUCTOR
    if( image_ptr.get() != 0 )
    {
      std::cerr << ::boost::format(" = %1%\n") % image_ptr.get();
    }
#endif
  }

  ~Gray_level_image_3()
  {
#ifdef CGAL_SURFACE_MESHER_DEBUG_GRAY_LEVEL_IMAGE_3_CONSTRUCTOR
      std::cerr << ::boost::format("~Gray_level_image_3() image=%1%\n") % image;
#endif
  }

  static void print_supported_file_format()
  {
    ::printSupportedFileFormat();
  }

  FT operator()(Point p) const
  {
    const float X=static_cast<float>(to_double(p.x()));
    const float Y=static_cast<float>(to_double(p.y()));
    const float Z=static_cast<float>(to_double(p.z()));

    float value = ::triLinInterp(this->image_ptr.get(), X, Y, Z, value_outside);
    if (positive_inside)
    {
      if (value > isovalue) // inside
        return FT(-1);
      else if (value < isovalue) // outside
        return FT(1);
      else
        return FT(0);
    }
    else
    {
      if (value < isovalue) // inside
        return FT(-1);
      else if (value > isovalue) // outside
        return FT(1);
      else
        return FT(0);
    }
  }
}; // end Gray_level_image_3

} // end namespace CGAL

#endif // CGAL_MESH_3_GRAY_LEVEL_IMAGE_3_H
