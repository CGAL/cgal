// Copyright (c) 2005-2008  INRIA Sophia-Antipolis (France).
//               2008 GeometryFactory, Sophia Antipolis (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Laurent Rineau, Pierre Alliez

#ifdef CGAL_HEADER_ONLY
#define CGAL_INLINE_FUNCTION inline
#else
#define CGAL_INLINE_FUNCTION
#endif

#include <CGAL/assertions.h>

namespace CGAL {

CGAL_INLINE_FUNCTION
bool Image_3::private_read(_image* im, Own own)
{
  if(im != 0)
  {
    if(image() != 0)
    {
      ::_freeImage(image());
    }
    image_ptr = Image_shared_ptr(im, Image_deleter(own == OWN_THE_DATA));

//     std::cerr <<
//       boost::format("image=%1% (xdim=%2%, ydim=%3%, zdim=%4%)\n")
//       % image_ptr.get() % image_ptr->xdim % image_ptr->ydim % image_ptr->zdim;

  }
  return im != 0;
}


} // end namespace CGAL

#if 0

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

CGAL_INLINE_FUNCTION
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

  return this->private_read(image);
}

} // end namespace CGAL

#endif // CGAL_USE_VTK
