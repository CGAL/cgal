// Copyright (c) 2005-2008  INRIA Sophia-Antipolis (France).
//               2008 GeometryFactory, Sophia Antipolis (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s)     : Laurent Rineau, Pierre Alliez

#include <CGAL/Image_3.h>
#include <CGAL/gl.h>

namespace CGAL {

bool Image_3::private_read(_image* im)
{
  if(im != 0)
  {
    if(image() != 0)
    {
      ::_freeImage(image());
    }
    image_ptr = Image_shared_ptr(im, Image_deleter());

//     std::cerr << 
//       boost::format("image=%1% (xdim=%2%, ydim=%3%, zdim=%4%)\n")
//       % image_ptr.get() % image_ptr->xdim % image_ptr->ydim % image_ptr->zdim;

  }
  return im != 0;
}

void Image_3::gl_draw(const float point_size,
                      const unsigned char r,
                      const unsigned char g,
                      const unsigned char b)
{
  if(image_ptr.get() == NULL)
    return;

  glPointSize(point_size);
  glColor3ub(r,g,b);
  glBegin(GL_POINTS);
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
          glVertex3d(x,y,z);
        }
      }
  glEnd();
} // end Image_3::gl_draw


void Image_3::gl_draw_bbox(const float line_width,
                           const unsigned char red,
                           const unsigned char green,
                           const unsigned char blue)
{
  if(!image_ptr)
    return;

  glLineWidth(line_width);
  glColor3ub(red,green,blue);
  glBegin(GL_LINES);

  struct Point {
    double x_;
    double y_;
    double z_;
    Point(double x, double y, double z) : x_(x), y_(y), z_(z) {};

    double x() const { return x_; }
    double y() const { return y_; }
    double z() const { return z_; }
  };

  const double xmax = (image_ptr->xdim - 1.0)*(image_ptr->vx);
  const double ymax = (image_ptr->ydim - 1.0)*(image_ptr->vy);
  const double zmax = (image_ptr->zdim - 1.0)*(image_ptr->vz);

  Point a(0.0, 0.0,    0.0);
  Point b(0.0, ymax, 0.0);
  Point c(0.0, ymax, zmax);
  Point d(0.0, 0.0,    zmax);
  Point e(xmax, 0.0,    0.0);
  Point f(xmax, ymax, 0.0);
  Point g(xmax, ymax, zmax);
  Point h(xmax, 0.0,    zmax);

  glVertex3d(a.x(),a.y(),a.z());
  glVertex3d(b.x(),b.y(),b.z());

  glVertex3d(b.x(),b.y(),b.z());
  glVertex3d(c.x(),c.y(),c.z());

  glVertex3d(c.x(),c.y(),c.z());
  glVertex3d(d.x(),d.y(),d.z());

  glVertex3d(d.x(),d.y(),d.z());
  glVertex3d(a.x(),a.y(),a.z());

  glVertex3d(e.x(),e.y(),e.z());
  glVertex3d(f.x(),f.y(),f.z());

  glVertex3d(f.x(),f.y(),f.z());
  glVertex3d(g.x(),g.y(),g.z());

  glVertex3d(g.x(),g.y(),g.z());
  glVertex3d(h.x(),h.y(),h.z());

  glVertex3d(h.x(),h.y(),h.z());
  glVertex3d(e.x(),e.y(),e.z());

  glVertex3d(a.x(),a.y(),a.z());
  glVertex3d(e.x(),e.y(),e.z());

  glVertex3d(d.x(),d.y(),d.z());
  glVertex3d(h.x(),h.y(),h.z());

  glVertex3d(c.x(),c.y(),c.z());
  glVertex3d(g.x(),g.y(),g.z());

  glVertex3d(b.x(),b.y(),b.z());
  glVertex3d(f.x(),f.y(),f.z());

  glEnd();
} // end Image_3::gl_draw_bbox

} // end namespace CGAL

#ifdef CGAL_USE_VTK

#include <vtkImageData.h>
#include <CGAL/Image_3_vtk_interface.h>

namespace CGAL {

namespace {

struct VTK_to_ImageIO_type_mapper {
  WORD_KIND wordKind;
  SIGN sign;
  unsigned int wdim;
};

static const VTK_to_ImageIO_type_mapper VTK_to_ImageIO_type[VTK_ID_TYPE] = 
  { { WK_UNKNOWN, SGN_UNKNOWN,  0}, //  0=VTK_VOID
    { WK_UNKNOWN, SGN_UNKNOWN,  0}, //  1=VTK_BIT
    { WK_FIXED,   SGN_SIGNED,   1}, //  2=VTK_CHAR
    { WK_FIXED,   SGN_UNSIGNED, 1}, //  3=VTK_UNSIGNED_CHAR
    { WK_FIXED,   SGN_SIGNED,   2}, //  4=VTK_SHORT 
    { WK_FIXED,   SGN_UNSIGNED, 2}, //  5=VTK_UNSIGNED_SHORT
    { WK_FIXED,   SGN_SIGNED,   4}, //  6=VTK_INT
    { WK_FIXED,   SGN_UNSIGNED, 4}, //  7=VTK_UNSIGNED_INT
    { WK_FIXED,   SGN_SIGNED,   8}, //  8=VTK_LONG
    { WK_FIXED,   SGN_UNSIGNED, 8}, //  9=VTK_UNSIGNED_LONG
    { WK_FLOAT,   SGN_SIGNED,   4}, // 10=VTK_FLOAT
    { WK_FIXED,   SGN_SIGNED,   8}  // 11=VTK_DOUBLE
 }; 

} //end anonymous namespace

bool
Image_3::read_vtk_image_data(vtkImageData* vtk_image)
{
  if(!vtk_image)
    return false;

  _image* image = ::_initImage();
  const int* dims = vtk_image->GetDimensions();
  const double* spacing = vtk_image->GetSpacing();
  image->vectMode = VM_SCALAR;
  image->xdim = dims[0];
  image->ydim = dims[1];
  image->zdim = dims[2];
  image->vdim = 1;
  image->vx = spacing[0];
  image->vy = spacing[1];
  image->vz = spacing[2];
  vtk_image->Update();
  image->endianness = ::_getEndianness();
  int vtk_type = vtk_image->GetScalarType();
  if(vtk_type == VTK_SIGNED_CHAR) vtk_type = VTK_CHAR;
  if(vtk_type < 0 || vtk_type > VTK_DOUBLE)
    vtk_type = VTK_DOUBLE;
  const VTK_to_ImageIO_type_mapper& imageio_type = 
    VTK_to_ImageIO_type[vtk_type];
  image->wdim = imageio_type.wdim;
  image->wordKind = imageio_type.wordKind;
  image->sign = imageio_type.sign;
  image->data = ::ImageIO_alloc(dims[0]*dims[1]*dims[2]*image->wdim);
  std::cerr << "GetNumberOfTuples()=" << vtk_image->GetPointData()->GetScalars()->GetNumberOfTuples()
            << "\nimage->size()=" << dims[0]*dims[1]*dims[2]
            << "\nwdim=" << image->wdim << '\n';
  assert(vtk_image->GetPointData()->GetScalars()->GetNumberOfTuples() == dims[0]*dims[1]*dims[2]);
  vtk_image->GetPointData()->GetScalars()->ExportToVoidPointer(image->data);

  return this->private_read(image);
}

} // end namespace CGAL

#endif // CGAL_USE_VTK
