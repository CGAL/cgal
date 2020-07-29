#include <CGAL/ImageIO.h>
#include <iostream>
#include <string>

int main(int argc, char **argv) {
  if (argc != 9) {
    std::cerr << "Usage: extract_a_sub_image <input> <output> <xmin> <xmax> "
                 "<ymin> <ymax> <zmin> <zmax>\n";
    return argc != 1;
  }
  _image *image = ::_readImage(argv[1]);
  if (!image)
    return 2;
  const auto xmin = std::stoul(argv[3]);
  const auto xmax = std::stoul(argv[4]);
  const auto ymin = std::stoul(argv[5]);
  const auto ymax = std::stoul(argv[6]);
  const auto zmin = std::stoul(argv[7]);
  const auto zmax = std::stoul(argv[8]);
  assert(xmax < image->xdim);
  assert(ymax < image->ydim);
  assert(zmax < image->zdim);

  const auto new_xdim = xmax + 1 - xmin;
  const auto new_ydim = ymax + 1 - ymin;
  const auto new_zdim = zmax + 1 - zmin;

  auto *new_image = ::_createImage(new_xdim, new_ydim, new_zdim, 1,
                                   image->vx, image->vy, image->vz, image->wdim,
                                   image->wordKind, image->sign);
  const auto* const data = static_cast<char*>(image->data);
  auto* new_data = static_cast<char*>(new_image->data);
  for (auto k = zmin; k < zmax; ++k)
    for (auto j = ymin; j <= ymax; ++j)
      for (auto i = xmin; i <= xmax; ++i) {
        auto pos = data + image->wdim * (i + image->xdim * (j + image->zdim * k));
        std::copy(pos, pos + image->wdim, new_data);
        new_data += image->wdim;
      }
  auto r = ::_writeImage(new_image, argv[2]);
  if(r != ImageIO_NO_ERROR) return 3;
  ::_freeImage(image);
  ::_freeImage(new_image);
}
