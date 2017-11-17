// Copyright (c) 2005-2008  INRIA Sophia-Antipolis (France).
//               2008 GeometryFactory
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
// SPDX-License-Identifier: LGPL-3.0+
//
//
// Author(s)     : Laurent Rineau, Pierre Alliez

#ifndef CGAL_READ_VTK_IMAGE_DATA_H
#define CGAL_READ_VTK_IMAGE_DATA_H

#include <CGAL/Image_3.h>
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

inline 
Image_3
read_vtk_image_data(vtkImageData* vtk_image)
{
  if(!vtk_image)
    return Image_3();

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
  CGAL_assertion(vtk_image->GetPointData()->GetScalars()->GetNumberOfTuples() == dims[0]*dims[1]*dims[2]);
  vtk_image->GetPointData()->GetScalars()->ExportToVoidPointer(image->data);

  return Image_3(image);
}

} // namespace CGAL

#endif // CGAL_READ_VTK_IMAGE_DATA_H
