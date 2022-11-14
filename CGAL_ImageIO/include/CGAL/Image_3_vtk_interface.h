// Copyright (c) 2008 GeometryFactory, Sophia-Antipolis (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Laurent Rineau


#ifndef CGAL_IMAGE_3_VTK_INTERFACE_H
#define CGAL_IMAGE_3_VTK_INTERFACE_H

#include <CGAL/config.h>

#include <CGAL/Image_3.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCharArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkShortArray.h>
#include <vtkUnsignedShortArray.h>
#include <vtkIntArray.h>
#include <vtkUnsignedIntArray.h>

#include <boost/cstdint.hpp> // for uint32_t, etc.

namespace CGAL {

template <typename Word>
struct VTK_type_generator {
};

template <>
struct VTK_type_generator<double> {
  static const int type = VTK_DOUBLE;
  typedef vtkDoubleArray ArrayType;
};

template <>
struct VTK_type_generator<float> {
  static const int type = VTK_FLOAT;
  typedef vtkFloatArray ArrayType;
};

template <>
struct VTK_type_generator<char> {
  static const int type = VTK_CHAR;
  typedef vtkCharArray ArrayType;
};

template <>
struct VTK_type_generator<boost::uint8_t> {
  static const int type = VTK_UNSIGNED_CHAR;
  typedef vtkUnsignedCharArray ArrayType;
};

template <>
struct VTK_type_generator<boost::int16_t> {
  static const int type = VTK_SHORT;
  typedef vtkShortArray ArrayType;
};

template <>
struct VTK_type_generator<boost::uint16_t> {
  static const int type = VTK_UNSIGNED_SHORT;
  typedef vtkUnsignedShortArray ArrayType;
};

template <>
struct VTK_type_generator<boost::int32_t> {
  static const int type = VTK_INT;
  typedef vtkIntArray ArrayType;
};

template <>
struct VTK_type_generator<boost::uint32_t> {
  static const int type = VTK_UNSIGNED_INT;
  typedef vtkUnsignedIntArray ArrayType;
};

inline
::vtkImageData* vtk_image_sharing_same_data_pointer(Image_3& image)
{
  vtkImageData* vtk_image = vtkImageData::New();
  vtkDataArray* data_array = 0;
  int type = 0;

  _image* image_ptr = image.image();

  CGAL_IMAGE_IO_CASE
    (image_ptr,
     type = VTK_type_generator<Word>::type;
     typedef VTK_type_generator<Word>::ArrayType VTKArray;
     VTKArray* array = VTKArray::New();
     array->SetArray(static_cast<Word*>(image.data()), image.size(), 1);
     data_array = array;
     );

  vtk_image->SetDimensions((int)image.xdim(),
                           (int)image.ydim(),
                           (int)image.zdim());
  vtk_image->SetExtent(0, (int)(image.xdim() - 1),
                       0, (int)(image.ydim() - 1),
                       0, (int)(image.zdim() - 1));
  vtk_image->SetSpacing(image.vx(),
                        image.vy(),
                        image.vz());
  vtk_image->SetOrigin(image.tx(),
                       image.ty(),
                       image.tz());
  vtk_image->AllocateScalars(type, 1);
  vtk_image->GetPointData()->SetScalars(data_array);
  return vtk_image;
} // end vtk_image_sharing_same_data_pointer

} //end namespace CGAL


#endif // CGAL_IMAGE_3_VTK_INTERFACE_H
