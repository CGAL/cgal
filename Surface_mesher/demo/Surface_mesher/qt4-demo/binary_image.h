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

#ifndef BINARY_IMAGE_3
#define BINARY_IMAGE_3

#include <CGAL/basic.h>

#include <boost/shared_ptr.hpp>

#include <boost/format.hpp>
#include <CGAL/ImageIO.h>

#include <limits>

template <typename FT, typename Point>
class CBinary_image_3
{
  struct Image_deleter {
#ifdef CGAL_SURFACE_MESHER_DEBUG_GRAY_LEVEL_IMAGE_3_CONSTRUCTOR
    void operator()(_image* image)
    {
      std::cerr << ::boost::format("Deletion of image %1%.\n") % image;
      ::_freeImage(image);
    }
#else
    void operator()(_image* ) {}
#endif
  };
public:
  typedef boost::shared_ptr<_image> Image_shared_ptr;
  Image_shared_ptr image_ptr;
  float m_isovalue;
  float min_x, min_y, min_z;
  float max_x, max_y, max_z;
  Point m_sink;
  float min_value;
  float max_value;

public:
  CBinary_image_3()
  {
    image_ptr = Image_shared_ptr();
    m_isovalue = 0.0f;
    max_x = max_y = max_z = 0.0f;
  }

  CBinary_image_3(const CBinary_image_3& bi)
  {
    std::cerr << "CBinary_image_3::copy_constructor\n";
    image_ptr = bi.image_ptr;
    m_isovalue = bi.m_isovalue;
    max_x = bi.max_x;
    max_y = bi.max_y;
    max_z = bi.max_z;
    m_sink = bi.m_sink;
    min_value = bi.min_value;
    max_value = bi.max_value;
  }

  float xmax() const { return max_x; }
  float ymax() const { return max_y; }
  float zmax() const { return max_z; }

  void compute_min_max()
  {
    if(image_ptr.get() == 0) {
      min_value = 0;
      max_value = 0;
      return;
    }
    min_value = max_value = evaluate(image_ptr.get(),0,0,0);
    for(int i = 0; i < xmax(); ++i) {
      for(int j = 0; j < ymax(); ++j)
        for(int k = 0; k < zmax(); ++k) {
          const float v = evaluate(image_ptr.get(), i, j, k);
          if(v > max_value) max_value = v;
          if(v < min_value) min_value = v;
        }
    }
  }

  Point center() 
  {
    FT cx = 0.5 * max_x;
    FT cy = 0.5 * max_y;
    FT cz = 0.5 * max_z;
    return Point(cx,cy,cz);
  }

  Point& sink() { return m_sink; }
  const Point& sink() const { return m_sink; }

  float& isovalue() { return m_isovalue; }
  float isovalue() const { return m_isovalue; }

  FT radius()
  {
    return std::max(std::max(max_x,max_y),max_z);
  }

  unsigned int xdim() { return image_ptr->xdim; }
  unsigned int ydim() { return image_ptr->ydim; }
  unsigned int zdim() { return image_ptr->zdim; }

  float value(const unsigned int i,
              const unsigned int j,
              const unsigned int k)
  {
    return evaluate(image_ptr.get(),i,j,k);
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

  bool read(const char* file,
            const float isoval)
  {
    if(image_ptr.get() != NULL)
      delete image_ptr.get();

    image_ptr = Image_shared_ptr(::_readImage(file), Image_deleter());
    std::cerr << boost::format("image=%1% (xdim=%2%, ydim=%3%, zdim=%4%)\n")
      % image_ptr.get() % image_ptr->xdim % image_ptr->ydim % image_ptr->zdim;
    if(image_ptr.get() != NULL)
    {
      m_isovalue = isoval;
      max_x = (float)(((image_ptr->xdim) - 1.0)*(image_ptr->vx));
      max_y = (float)(((image_ptr->ydim) - 1.0)*(image_ptr->vy));
      max_z = (float)(((image_ptr->zdim) - 1.0)*(image_ptr->vz));
      compute_min_max();
      return true;
    }
    compute_min_max();
    return false;
  }

  bool read_raw(const char* file,
                const float isoval,
                const unsigned int rx,
                const unsigned int ry,
                const unsigned int rz)
  {
    if(image_ptr.get() != NULL)
      delete image_ptr.get();

    image_ptr = Image_shared_ptr(::_readImage_raw(file,rx,ry,rz), Image_deleter());
    if(image_ptr.get() != NULL)
    {
      m_isovalue = isoval;
      max_x = (float)(((image_ptr->xdim) - 1.0)*(image_ptr->vx));
      max_y = (float)(((image_ptr->ydim) - 1.0)*(image_ptr->vy));
      max_z = (float)(((image_ptr->zdim) - 1.0)*(image_ptr->vz));
      compute_min_max();
      return true;
    }
    compute_min_max();
    return false;
  }

  ~CBinary_image_3()
  {
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

  unsigned int threshold(const unsigned char value,
                         const unsigned char equal,
                         const unsigned char diff)
  {
    if(image_ptr.get() == NULL)
      return 0;

    unsigned int nb = 0;
    unsigned char *pData = (unsigned char*)image_ptr->data;
    unsigned int xy = image_ptr->xdim * image_ptr->ydim;
    for(unsigned int i=0;i<image_ptr->xdim;i++)
      for(unsigned int j=0;j<image_ptr->ydim;j++)
        for(unsigned int k=0;k<image_ptr->zdim;k++)
        {
          unsigned char voxel = pData[xy*k + j*image_ptr->xdim + i];
          if(voxel == value)
          {
            pData[xy*k + j*image_ptr->xdim + i] = equal;
            nb++;
          }
          else
            pData[xy*k + j*image_ptr->xdim + i] = diff;
        }
    return nb;
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

  bool inside(const float x,
              const float y, 
              const float z) const
  {
    return ( x >= 0.0 && 
             y >= 0.0 && 
             z >= 0.0 && 
             x <= max_x &&
             y <= max_y &&
             z <= max_z );
  }

  bool probe_sink(const unsigned int nb_max)
  {
    unsigned int nb = 0;
    while(!probe(m_sink)) 
    {
      nb++;
      if(nb > nb_max)
        return false;
    }
    return true;
  }

  bool probe(Point& center)
  {
    float rx = rand_x();
    float ry = rand_y();
    float rz = rand_z();
    float value = ::trilinear_interpolation(image_ptr.get(),rx,ry,rz); 
    if(value > m_isovalue)
    {
      center = Point(rx,ry,rz);
      return true;
    }
    return false;
  }

  float rand_x() { return (float)rand() / (float)RAND_MAX * max_x; }
  float rand_y() { return (float)rand() / (float)RAND_MAX * max_y; }
  float rand_z() { return (float)rand() / (float)RAND_MAX * max_z; }

  FT operator()(Point p) const
  {
    float x = static_cast<float>(CGAL::to_double(p.x()));
    float y = static_cast<float>(CGAL::to_double(p.y()));
    float z = static_cast<float>(CGAL::to_double(p.z()));

    if(!inside(x,y,z))
      return FT(1);
    else
    {
      float value = ::trilinear_interpolation(image_ptr.get(),x,y,z); 
      if (value > m_isovalue) // inside
        return FT(-1);
      else
        if (value < m_isovalue) // outside
          return FT(1);
        else
          return FT(0);
    }
  }
}; // end CBinary_image_3
 
#endif // BINARY_IMAGE_3
