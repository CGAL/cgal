#include <CGAL/Image_3.h>
#include <CGAL/Random.h>
#include <algorithm>

CGAL::Image_3 random_labeled_image()
{
  const int dim = 400;
  const unsigned char number_of_spheres = 50;
  const int max_radius_of_spheres = 10;
  _image* image = _createImage(dim, dim, dim, 1,
                               1.f, 1.f, 1.f, 1,
                               WK_FIXED, SGN_UNSIGNED);
  unsigned char* ptr = (unsigned char*)(image->data);
  std::fill(ptr, ptr+dim*dim*dim, '\0');

  CGAL::Random rand(0);
  for(unsigned char n = number_of_spheres; n > 0 ; --n) {
    std::size_t i = rand.uniform_smallint(    1 + max_radius_of_spheres,
                                          dim-2 - max_radius_of_spheres);
    std::size_t j = rand.uniform_smallint(    1 + max_radius_of_spheres,
                                          dim-2 - max_radius_of_spheres);
    std::size_t k = rand.uniform_smallint(    1 + max_radius_of_spheres,
                                          dim-2 - max_radius_of_spheres);
    for(std::ptrdiff_t ii = - max_radius_of_spheres;
        ii <= max_radius_of_spheres; ++ii)
    {
      for(std::ptrdiff_t jj = - max_radius_of_spheres;
          jj <= max_radius_of_spheres; ++jj)
      {
        for(std::ptrdiff_t kk = - max_radius_of_spheres;
            kk <= max_radius_of_spheres; ++kk)
        {
          if(ii*ii + jj*jj + kk*kk >
             max_radius_of_spheres * max_radius_of_spheres) continue;
          using CGAL::IMAGEIO::static_evaluate;
          static_evaluate<unsigned char>(image, i+ii, j+jj, k+kk) = n;
        }
      }
    }
  }
  return CGAL::Image_3(image);
}
