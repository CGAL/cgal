#ifndef CGAL_MESH_3_ISOSURFACE_H
#define CGAL_MESH_3_ISOSURFACE_H

#include "imageio/ImageIO.h"

namespace CGAL {

class Inrimage_isosurface
{
  _image *image;
  float isovalue;
  float max_x, max_y, max_z;

public:
  Inrimage_isosurface(const char* file, float isoval=0)
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
      float value = triLinInterp(image, X, Y, Z); 

      if (value > isovalue)
	return -1;
      else if (value < isovalue)
	return 1;
      else
	return 0;
    }
  }
}; // end Inrimage_isosurface
 
} // end namespace CGAL

#endif // CGAL_MESH_3_ISOSURFACE_H
