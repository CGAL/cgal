// Copyright (c) 2005-2006  INRIA Sophia-Antipolis (France).
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
// $Source: /proj/geometrica/home/CGAL/Local/cvsroot/InrImage/include/CGAL/Isosurface.h,v $
// $Revision: 1.2 $ $Date: 2005/08/01 17:49:26 $
// $Name:  $
//
// Author(s)     : Laurent RINEAU

#ifndef CGAL_MESH_3_GRAY_LEVEL_IMAGE_3_H
#define CGAL_MESH_3_GRAY_LEVEL_IMAGE_3_H

#include "imageio/ImageIO.h"

namespace CGAL {

template <typename FT>
class Gray_level_image_3
{
  _image *image;
  float isovalue;
  float max_x, max_y, max_z;

public:
  Gray_level_image_3(const char* file, float isoval=0)
  {
    image = ::_readImage(file);
    ::convertImageTypeToFloat(image);
    isovalue=isoval;
    max_x = ((image->xdim) - 1.0)*(image->vx);
    max_y = ((image->ydim) - 1.0)*(image->vy);
    max_z = ((image->zdim) - 1.0)*(image->vz);
  }

  bool inside(float X, float Y, float Z) const
  {
    return ( X>=0. && Y>=0. && Z>=0. && 
             X<=max_x && Y<=max_y && Z<=max_z );
  }

  int operator()(double x, double y, double z) const
  {
    float X=static_cast<float>(x);
    float Y=static_cast<float>(y);
    float Z=static_cast<float>(z);

    if (!inside(X,Y,Z))
      return 1;
    else{
      float value = ::triLinInterp(image, X, Y, Z); 

      if (value > isovalue)
	return -1;
      else if (value < isovalue)
	return 1;
      else
	return 0;
    }
  }
}; // end Gray_level_image_3
 
} // end namespace CGAL

#endif // CGAL_MESH_3_GRAY_LEVEL_IMAGE_3_H
