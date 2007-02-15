// Copyright (c) 2005-2007  INRIA Sophia-Antipolis (France).
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
// $Revision$ $Date$
// $Name:  $
//
// Author(s)     : Laurent RINEAU

#ifndef CGAL_MESH_3_GRAY_LEVEL_IMAGE_3_H
#define CGAL_MESH_3_GRAY_LEVEL_IMAGE_3_H

#include <CGAL/basic.h>

struct _image;

/* Copy-paste from <imageio/ImageIO.h> */
_image* _readImage(const char *name);
void convertImageTypeToFloat(_image* image);
float triLinInterp(_image* image,float posx, float posy, float posz);
/* End of copy-paste from <imageio/ImageIO.h> */

namespace CGAL {

  template <typename FT, typename Point>
class Gray_level_image_3
{
  _image *image;
  float isovalue;
  float max_x, max_y, max_z;

public:
  Gray_level_image_3(const char* file, float isoval)
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

  FT operator()(Point p) const
  {
    float X=static_cast<float>(to_double(p.x()));
    float Y=static_cast<float>(to_double(p.y()));
    float Z=static_cast<float>(to_double(p.z()));

    if (!inside(X,Y,Z))
      return FT(1);
    else{
      float value = ::triLinInterp(image, X, Y, Z); 

      if (value > isovalue) // inside
	return FT(-1);
      else if (value < isovalue) // outside
	return FT(1);
      else
	return FT(0);
    }
  }
}; // end Gray_level_image_3
 
} // end namespace CGAL

#endif // CGAL_MESH_3_GRAY_LEVEL_IMAGE_3_H
