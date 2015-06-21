// Copyright (c) 2005-2008  INRIA Sophia-Antipolis (France).
// Copyright (c) 2008  GeometryFactory (France)
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
// Author(s)     : Laurent RINEAU, Pierre ALLIEZ

#ifndef BINARY_IMAGE_3
#define BINARY_IMAGE_3

#include <CGAL/basic.h>

#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>
#include <CGAL/algorithm.h>

#include <boost/format.hpp>
#include <CGAL/ImageIO.h>

#include <limits>

#include <CGAL/Image_3.h>

template <typename FT_, typename Point>
class CBinary_image_3 : public CGAL::Image_3
{
  bool interpolate_;
  bool labellized_;

public:
  double min_value;
  double max_value;

  typedef FT_ FT;

public:
  CBinary_image_3() : Image_3(), interpolate_(true)
  {
  }

  CBinary_image_3(const CBinary_image_3& bi)
    : Image_3(bi), interpolate_(bi.interpolate_),labellized_(bi.labellized_)
  {
    std::cerr << "CBinary_image_3::copy_constructor\n";
    min_value = bi.min_value;
    max_value = bi.max_value;
  }

  ~CBinary_image_3()
  {
  }

  void finish_open() {
    CGAL_IMAGE_IO_CASE(image_ptr.get(),
                       Word *min; Word *max;
                       (boost::tie(min, max)) = 
                         (CGAL::min_max_element((Word*)(data()),
                                                (Word*)(data()) + 
                                                xdim() * ydim() * zdim()));
                        min_value = *min;
                        max_value = *max;)
  }

  float xmax() const
  {
    return (float)(((image_ptr->xdim) - 1.0)*(image_ptr->vx));
  }

  float ymax() const
  {
    return (float)(((image_ptr->ydim) - 1.0)*(image_ptr->vy));
  }

  float zmax() const
  {
    return (float)(((image_ptr->zdim) - 1.0)*(image_ptr->vz));
  }

  Point center() 
  {
    FT cx = 0.5 * xmax();
    FT cy = 0.5 * ymax();
    FT cz = 0.5 * zmax();
    return Point(cx,cy,cz);
  }

  FT radius()
  {
    return (std::max)((std::max)(xmax(),ymax()),zmax());
  }

  Point point(const std::size_t i,
              const std::size_t j,
              const std::size_t k) const
  {
    return Point(i * (image_ptr->vx),
                 j * (image_ptr->vy),
                 k * (image_ptr->vz));
  }

public:
  bool inside(const float x,
              const float y, 
              const float z) const
  {
    return ( x >= 0.0f && 
             y >= 0.0f && 
             z >= 0.0f && 
             x <= xmax() &&
             y <= ymax() &&
             z <= zmax() );
  }

  float rand_x() { return (float)rand() / (float)RAND_MAX * xmax(); }
  float rand_y() { return (float)rand() / (float)RAND_MAX * ymax(); }
  float rand_z() { return (float)rand() / (float)RAND_MAX * zmax(); }

  void set_interpolation(const bool b)
  {
    interpolate_ = b;
  }

  bool interpolation() const {
    return interpolate_;
  }

  void set_labellized(const bool b)
  {
    labellized_ = b;
  }

  bool labellized() const {
    return labellized_;
  }

  FT operator()(Point p) const
  {
    const float x = static_cast<float>(CGAL::to_double(p.x()));
    const float y = static_cast<float>(CGAL::to_double(p.y()));
    const float z = static_cast<float>(CGAL::to_double(p.z()));
      
    if(interpolation()) {
      if(labellized()) {
	CGAL_IMAGE_IO_CASE(image_ptr.get(),
			   return (this->labellized_trilinear_interpolation<Word, double>(x, y, z, min_value));)
      }
      else {
	CGAL_IMAGE_IO_CASE(image_ptr.get(),
 			   return (this->trilinear_interpolation<Word, double>(x, y, z, min_value));)
      }
    }
    else {
      const std::ptrdiff_t i = static_cast<std::ptrdiff_t>(x/image()->vx + 0.5f);
      const std::ptrdiff_t j = static_cast<std::ptrdiff_t>(y/image()->vy + 0.5f);
      const std::ptrdiff_t k = static_cast<std::ptrdiff_t>(z/image()->vz + 0.5f);
      if( i < 0 ||
	  j < 0 ||
	  k < 0 )
      {
	return 0;
      }
      else
      {    
	const std::size_t ui = static_cast<std::size_t>(i);
	const std::size_t uj = static_cast<std::size_t>(j);
	const std::size_t uk = static_cast<std::size_t>(k);
	if( ui >= image()->xdim ||
	    uj >= image()->ydim ||
	    uk >= image()->zdim )
	{
	  return 0;
	}
	else
	{
	  return this->value(ui, uj, uk);
	}
      }
    }
    return FT();
  }
}; // end CBinary_image_3
 
#endif // BINARY_IMAGE_3
