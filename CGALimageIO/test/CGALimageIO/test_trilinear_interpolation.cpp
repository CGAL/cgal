#include <CGAL/config.h>
#include <CGAL/Image_3.h>

#include <cassert>
#include <iostream>

typedef unsigned char Word;

int main() {
  CGAL::Image_3 image(_createImage(2, 2, 2, 1,
				   1., 1., 1.,
				   1, WK_FIXED, SGN_UNSIGNED));
  
  Word data[8] = { 0, 0, 0, 0, 0, 0, 0, 0};

  ImageIO_free(image.data());
  image.set_data(&data[0]);
  _image* im = image.image();

//   im->xdim = 2;
//   im->ydim = 2;
//   im->zdim = 2;
//   im->vdim = 1;
//   im->vx = 1;
//   im->vy = 1;
//   im->vz = 1;

//   im->cx = im->cy = im->cz = 0;
//   im->tx = im->ty = im->tz = 0.0;
//   im->rx = im->ry = im->rz = 0.0;


//   im->fd = NULL;
//   im->openMode = OM_CLOSE;
//   im->endianness = END_UNKNOWN;

//   im->dataMode = DM_BINARY;

//   // word type (unsigned char)
//   im->wdim = 1;
//   im->wordKind = WK_FIXED;
//   im->vectMode = VM_SCALAR;
//   im->sign = SGN_UNSIGNED;
//   im->imageFormat = NULL;
  std::cerr << std::setprecision(2) << std::fixed;

  for(int x = 0; x <= 1; ++x)
    for(int y = 0; y <= 1; ++y)
      for(int z = 0; z <= 1; ++z)
      {
	data[x * 4 + y * 2 + z] = 1;

	std::cerr << "#### data"
		  << "[" << x << "]"
		  << "[" << y << "]"
		  << "[" << z << "]"
		  << " = 1\n";
	for(double d_x = 0.; d_x <= 1.; d_x += 0.499999)
	  for(double d_y = 0.; d_y <= 1.; d_y += 0.499999)
	    for(double d_z = 0.; d_z <= 1.; d_z += 0.499999)
	    {
	      assert((int)d_x == 0);
	      assert((int)d_y == 0);
	      assert((int)d_z == 0);
	      const double value = 
		image.trilinear_interpolation<Word, double, float>(d_x,
								   d_y,
								   d_z,
								   255);
	      std::cerr << "val(" << d_x << ", " << d_y << " , " << d_z << ") = "
			<< value << std::endl;
	      const double sq_dist = 
		(d_x - x) * (d_x - x) +
		(d_y - y) * (d_y - y) +
		(d_z - z) * (d_z - z);
	      if(sq_dist <= 0.001)
		assert(value >= 0.9);
	      else if(sq_dist <= 0.51*0.51)
		assert(value >= 0.49);
	      else if(sq_dist <= 2*0.51*0.51)
		assert(value >= 0.24);
	      else if(sq_dist <= 3*0.51*0.51)
		assert(value >= 0.12);
	      else
		assert(value <= 0.001);
	    }

	data[x * 4 + y * 2 + z] = 0;
      }
  image.set_data(0);
}
