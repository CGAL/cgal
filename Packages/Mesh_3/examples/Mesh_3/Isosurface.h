#ifndef CGAL_MESH_3_ISOSURFACE_H
#define CGAL_MESH_3_ISOSURFACE_H

#include <inrimage/Inrimage.h>

namespace CGAL {

template <class FT> 
class Inrimage_isosurface
{
  yav::Inrimage *entree;
  yav::Voxel *voxl;
  double isovalue;

public:
  enum SURFACE_LOCATION {IN = -1, ON = 0, OUT = 1};

  Inrimage_isosurface(char* file, double isoval=0)
  {
    entree = new yav::Inrimage(file);
    yav::Inrimage::WORD_TYPE wt = entree->getType();
    voxl = yav::Voxel::voxelByType(wt);
    isovalue=isoval;
  }

  SURFACE_LOCATION operator()(FT x, FT y, FT z) const
  {
    unsigned X=static_cast<int>(floor(CGAL::to_double (x)));
    unsigned Y=static_cast<int>(floor(CGAL::to_double (y)));
    unsigned Z=static_cast<int>(floor(CGAL::to_double (z)));

    bool inside = entree->inside(X,Y,Z);

    if (!inside)
      return IN;
    else{
      entree->getRealVoxel(x,y,z,*voxl);
      double value = voxl->getValueAsDouble();    

      if (value < isovalue)
	return IN;
      else if (value > isovalue)
	return OUT;
      else
	return ON;
    }
  }
}; // end Inrimage_isosurface
 
} // end namespace CGAL

#endif // CGAL_MESH_3_ISOSURFACE_H
