#include <CGAL/Image_3.h>
#include <CGAL/Random.h>
#include <algorithm>

CGAL::Image_3 random_labeled_image()
{
  const int dim = 400;
  const unsigned char number_of_spheres = 50;
  const int max_radius_of_spheres = 10;
  const int radius_of_big_sphere = 80;
  _image* image = _createImage(dim, dim, dim, 1,
                               1.f, 1.f, 1.f, 1,
                               WK_FIXED, SGN_UNSIGNED);
  unsigned char* ptr = (unsigned char*)(image->data);
  std::fill(ptr, ptr+dim*dim*dim, '\0');

  std::ptrdiff_t center = dim / 2;
  CGAL::Random rand(0);
  for(unsigned char n = number_of_spheres; n > 0 ; --n) {
    std::size_t i, j, k;
    do {
      i = rand.uniform_smallint(1 + max_radius_of_spheres,
                                dim-2 - max_radius_of_spheres);
      j = rand.uniform_smallint(1 + max_radius_of_spheres,
                                dim-2 - max_radius_of_spheres);
      k = rand.uniform_smallint(1 + max_radius_of_spheres,
                                dim-2 - max_radius_of_spheres);
    } while ( ( CGAL::square(double(center) - double(i)) +
                CGAL::square(double(center) - double(j)) +
                CGAL::square(double(center) - double(k)) )
              <
              CGAL::square(double(radius_of_big_sphere) + 4 * max_radius_of_spheres) );
    std::ptrdiff_t radius = max_radius_of_spheres;
    if(n==1) {
      i = j = k = center;
      radius = radius_of_big_sphere;
    }
    for(std::ptrdiff_t ii = - radius; ii <= radius; ++ii)
    {
      for(std::ptrdiff_t jj = - radius; jj <= radius; ++jj)
      {
        for(std::ptrdiff_t kk = - radius; kk <= radius; ++kk)
        {
          if(ii*ii + jj*jj + kk*kk > radius * radius) continue;
          using CGAL::IMAGEIO::static_evaluate;
          static_evaluate<unsigned char>(image, i+ii, j+jj, k+kk) = n;
        }
      }
    }
  }
  _writeImage(image, "random-image.inr");
  return CGAL::Image_3(image);
}
