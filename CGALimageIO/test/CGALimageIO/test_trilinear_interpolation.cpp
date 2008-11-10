#include <CGAL/config.h>
#include <CGAL/Image_3.h>
#include <CGAL/Random.h>
#include <CGAL/Timer.h>

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

  std::cerr << std::setprecision(2) << std::fixed;

  for(int x = 0; x <= 1; ++x)
    for(int y = 0; y <= 1; ++y)
      for(int z = 0; z <= 1; ++z)
      {
	data[z * 4 + y * 2 + x] = 1;

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
		image.trilinear_interpolation<Word, double, double>(d_x,
								    d_y,
								    d_z,
								    255);

	      const Word label = 
		image.labellized_trilinear_interpolation<Word,
	          double>(d_x, d_y, d_z, 255);

	      std::cerr << "val(" << d_x << ", " << d_y << " , " << d_z << ") = "
			<< value << std::endl;
	      const double sq_dist = 
		(d_x - x) * (d_x - x) +
		(d_y - y) * (d_y - y) +
		(d_z - z) * (d_z - z);

	      if(sq_dist <= 0.001) {
		assert(value >= 0.9);
		assert(label == 1);
	      }
	      else if(sq_dist <= 0.51*0.51/2.) {
		assert(value >= 0.79);
		assert(label == 1);
	      }
	      else if(sq_dist <= 0.51*0.51) {
		assert(value >= 0.49);
	      }
	      else if(sq_dist <= 2*0.51*0.51) {
		assert(value >= 0.24);
		assert(label == 0);
	      }
	      else if(sq_dist <= 3*0.51*0.51) {
		assert(value >= 0.12);
		assert(label == 0);
	      }
	      else {
		assert(value <= 0.001);
		assert(label == 0);
	      }

	      const float value2 = 
		image.trilinear_interpolation<Word, float, float>(d_x,
								  d_y,
								  d_z,
								  0);

	      const float value3 = triLinInterp(image.image(),
						d_x, 
						d_y,
						d_z,
						0.f);
	      std::cerr << "tri(" << d_x << ", " << d_y << " , " << d_z << ") = "
			<< value3 << std::endl;
	      if(value2 != value3)
		std::cerr << std::setprecision(30)
			  << "   " << value2 << "\n!= " << value3 << std::endl;
	      assert(value2 == value3);
	    }

	data[z * 4 + y * 2 + x] = 0;
      }

  // BENCH
  std::cerr.precision(10);

  float max_diff = 0.f;
  CGAL::Timer timer_new_implementation;
  CGAL::Timer timer_old_implementation;

  for(int image_nb = 0; image_nb < 1000; ++image_nb)
  {
    // fill the 2x2x2 image
    for(int i = 0; i < 8; ++i) {
      data[i] = CGAL::default_random.get_int(0,255);
    }

    // test the difference between the two implementations
    for(double d_x = 0.; d_x < 0.9; d_x += 0.1)
      for(double d_y = 0.; d_y < 0.9; d_y += 0.1)
	for(double d_z = 0.; d_z < 0.9; d_z += 0.1)
	{

	  const float value1 = 
	    image.trilinear_interpolation<Word, float, float>(d_x,
							      d_y,
							      d_z,
							      0);
	  const float value2 = triLinInterp(image.image(),
					    d_x, 
					    d_y,
					    d_z,
					    0.f);
	  float diff = value1 - value2;
	  if(diff < 0) {
	    diff =-diff;
	  }
	  assert(diff < 0.1f);
	  if(diff > max_diff)
	    max_diff = diff;
	}
    // bench new implementation
    float sum = 0.f;
    timer_new_implementation.start();
    for(double d_x = 0.; d_x < 0.9; d_x += 0.05)
      for(double d_y = 0.; d_y < 0.9; d_y += 0.05)
	for(double d_z = 0.; d_z < 0.9; d_z += 0.05)
	{

	  sum += 
	    image.trilinear_interpolation<Word, float, float>(d_x,
							      d_y,
							      d_z,
							      0);
	}
    timer_new_implementation.stop();

    sum = 0.f;
    timer_old_implementation.start();
    // bench old implementation    
    for(double d_x = 0.; d_x < 0.9; d_x += 0.05)
      for(double d_y = 0.; d_y < 0.9; d_y += 0.05)
	for(double d_z = 0.; d_z < 0.9; d_z += 0.05)
	{
	  sum += triLinInterp(image.image(),
			      d_x, 
			      d_y,
			      d_z,
			      0.f);
	}
    timer_old_implementation.stop();
  }
  std::cerr << "max difference = " << max_diff << "\n"
	    << "timer new implementation: " << timer_new_implementation.time()
	    << "\ntimer old implementation: " << timer_old_implementation.time()
	    << "\n";
  image.set_data(0); // trick to avoid ~Image_3 segfault.
}
