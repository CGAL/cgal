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
#include <GL/gl.h>

#include <limits>

namespace CGAL {

template <typename FT_, typename Point>
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
  Image_shared_ptr image_ptr;
  float max_x, max_y, max_z;
  float min_value;
  float max_value;

  typedef FT_ FT;

  bool private_read(_image* im)
  {
    if(im != 0)
    {
      if(image() != 0)
      {
        ::_freeImage(image());
      }
      image_ptr = Image_shared_ptr(im, Image_deleter());

      std::cerr << 
        boost::format("image=%1% (xdim=%2%, ydim=%3%, zdim=%4%)\n")
        % image_ptr.get() % image_ptr->xdim % image_ptr->ydim % image_ptr->zdim;

      max_x = (float)(((image_ptr->xdim) - 1.0)*(image_ptr->vx));
      max_y = (float)(((image_ptr->ydim) - 1.0)*(image_ptr->vy));
      max_z = (float)(((image_ptr->zdim) - 1.0)*(image_ptr->vz));
    }
    compute_min_max();
    return im != 0;
  }

  void compute_min_max()
  {
    if(image() == 0) {
      min_value = 0;
      max_value = 0;
      return;
    }
    min_value = max_value = ::evaluate(image(),0,0,0);
    for(unsigned int i = 0; i < xdim(); ++i) {
      for(unsigned int j = 0; j < ydim(); ++j) {
        for(unsigned int k = 0; k < zdim(); ++k) {
          const float v = ::evaluate(image(), i, j, k);
          if(v > max_value) max_value = v;
          if(v < min_value) min_value = v;
        }
      }
    }
  }

public:
  Image_3()
  {
    image_ptr = Image_shared_ptr();
    max_x = max_y = max_z = 0.0f;
  }

  Image_3(const Image_3& bi)
  {
    std::cerr << "Image_3::copy_constructor\n";
    image_ptr = bi.image_ptr;
    max_x = bi.max_x;
    max_y = bi.max_y;
    max_z = bi.max_z;
    min_value = bi.min_value;
    max_value = bi.max_value;
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

  float xmax() const { return max_x; }
  float ymax() const { return max_y; }
  float zmax() const { return max_z; }

  unsigned int xdim() const { return image_ptr->xdim; }
  unsigned int ydim() const { return image_ptr->ydim; }
  unsigned int zdim() const { return image_ptr->zdim; }

  float value(const unsigned int i,
              const unsigned int j,
              const unsigned int k) const
  {
    return ::evaluate(image(),i,j,k);
  }

  Point point(const unsigned int i,
              const unsigned int j,
              const unsigned int k) const
  {
    return Point(i * (image_ptr->vx),
                 j * (image_ptr->vy),
                 k * (image_ptr->vz));
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

  void gl_draw(const float point_size,
               const unsigned char r,
               const unsigned char g,
               const unsigned char b)
  {
    if(image_ptr.get() == NULL)
      return;

    ::glPointSize(point_size);
    ::glColor3ub(r,g,b);
    ::glBegin(GL_POINTS);
    unsigned char *pData = (unsigned char*)image_ptr->data;
    unsigned int xy = image_ptr->xdim * image_ptr->ydim;
    for(unsigned int i=0;i<image_ptr->xdim;i+=5)
      for(unsigned int j=0;j<image_ptr->ydim;j+=5)
        for(unsigned int k=0;k<image_ptr->zdim;k+=5)
        {
          unsigned char value = pData[xy*k + j*image_ptr->xdim + i];
          if(value > 0)
          {
            double x = image_ptr->vx * i;
            double y = image_ptr->vy * j;
            double z = image_ptr->vz * k;
            ::glVertex3d(x,y,z);
          }
        }
    ::glEnd();
  }

  void gl_draw_bbox(const float line_width,
                    const unsigned char red,
                    const unsigned char green,
                    const unsigned char blue)
  {
    ::glLineWidth(line_width);
    ::glColor3ub(red,green,blue);
    ::glBegin(GL_LINES);

    Point a(0.0, 0.0,    0.0);
    Point b(0.0, ymax(), 0.0);
    Point c(0.0, ymax(), zmax());
    Point d(0.0, 0.0,    zmax());
    Point e(xmax(), 0.0,    0.0);
    Point f(xmax(), ymax(), 0.0);
    Point g(xmax(), ymax(), zmax());
    Point h(xmax(), 0.0,    zmax());

    ::glVertex3d(a.x(),a.y(),a.z());
    ::glVertex3d(b.x(),b.y(),b.z());

    ::glVertex3d(b.x(),b.y(),b.z());
    ::glVertex3d(c.x(),c.y(),c.z());

    ::glVertex3d(c.x(),c.y(),c.z());
    ::glVertex3d(d.x(),d.y(),d.z());

    ::glVertex3d(d.x(),d.y(),d.z());
    ::glVertex3d(a.x(),a.y(),a.z());

    ::glVertex3d(e.x(),e.y(),e.z());
    ::glVertex3d(f.x(),f.y(),f.z());

    ::glVertex3d(f.x(),f.y(),f.z());
    ::glVertex3d(g.x(),g.y(),g.z());

    ::glVertex3d(g.x(),g.y(),g.z());
    ::glVertex3d(h.x(),h.y(),h.z());

    ::glVertex3d(h.x(),h.y(),h.z());
    ::glVertex3d(e.x(),e.y(),e.z());

    ::glVertex3d(a.x(),a.y(),a.z());
    ::glVertex3d(e.x(),e.y(),e.z());

    ::glVertex3d(d.x(),d.y(),d.z());
    ::glVertex3d(h.x(),h.y(),h.z());

    ::glVertex3d(c.x(),c.y(),c.z());
    ::glVertex3d(g.x(),g.y(),g.z());

    ::glVertex3d(b.x(),b.y(),b.z());
    ::glVertex3d(f.x(),f.y(),f.z());

    ::glEnd();
  }

//   FT operator()(Point p) const
//   {
//     const float x = static_cast<float>(CGAL::to_double(p.x()));
//     const float y = static_cast<float>(CGAL::to_double(p.y()));
//     const float z = static_cast<float>(CGAL::to_double(p.z()));

//     if(inside(x,y,z))
//       return FT(::trilinear_interpolation(image_ptr.get(),x,y,z));
//     else
//       return 0;
//   }
}; // end Image_3

} // end namespace CGAL
 
#endif // CGAL_IMAGE_3
