// Copyright (c) 2005-2008  INRIA Sophia-Antipolis (France).
//               2008 GeometryFactory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// Author(s)     : Laurent RINEAU, Pierre ALLIEZ

#ifndef CGAL_IMAGE_3
#define CGAL_IMAGE_3

#include <CGAL/basic.h>

#include <boost/shared_ptr.hpp>

#include <boost/format.hpp>
#include <CGAL/ImageIO.h>

#include <limits>

class vtkImageData;

namespace CGAL {

class Image_3
{
  struct Image_deleter {
    void operator()(_image* image)
    {
      ::_freeImage(image);
    }
  };
public:
  typedef boost::shared_ptr<_image> Image_shared_ptr;
  typedef Image_shared_ptr Pointer;

protected:
  Image_shared_ptr image_ptr;

  bool private_read(_image* im); // implementation in src/CGALimageIO/Image_3.cpp

public:
  Image_3()
    : image_ptr()
  {
  }

  Image_3(const Image_3& bi)
    : image_ptr(bi.image_ptr)
  {
    std::cerr << "Image_3::copy_constructor\n";
  }

  ~Image_3()
  {
  }

  const _image* image() const
  {
    return image_ptr.get();
  }

  _image* image()
  {
    return image_ptr.get();
  }

  const void* data() const
  {
    return image()->data;
  }

  void* data()
  {
    return image()->data;
  }

  unsigned int xdim() const { return image_ptr->xdim; }
  unsigned int ydim() const { return image_ptr->ydim; }
  unsigned int zdim() const { return image_ptr->zdim; }

  unsigned int size() const { return xdim() * ydim() * zdim(); }

  float vx() const { return image_ptr->vx; }
  float vy() const { return image_ptr->vy; }
  float vz() const { return image_ptr->vz; }

  float value(const unsigned int i,
              const unsigned int j,
              const unsigned int k) const
  {
    return ::evaluate(image(),i,j,k);
  }

public:

  bool read(const char* file)
  {
    return private_read(::_readImage(file));
  }

  bool read_raw(const char* file,
                const unsigned int rx,
                const unsigned int ry,
                const unsigned int rz,
                const float vx = 1,
                const float vy = 1,
                const float vz = 1)
  {
    return private_read(::_readImage_raw(file,
                                         rx,ry,rz,
                                         vx,vy,vz));
  }

#ifdef CGAL_USE_VTK
  bool read_vtk_image_data(vtkImageData*);
#endif // CGAL_USE_VTK

  void gl_draw(const float point_size,
               const unsigned char r,
               const unsigned char g,
               const unsigned char b); // implementation in src/CGALimageIO/Image_3.cpp

  void gl_draw_bbox(const float line_width,
                    const unsigned char red,
                    const unsigned char green,
                    const unsigned char blue); // implementation in src/CGALimageIO/Image_3.cpp

}; // end Image_3

} // end namespace CGAL

 
#endif // CGAL_IMAGE_3
