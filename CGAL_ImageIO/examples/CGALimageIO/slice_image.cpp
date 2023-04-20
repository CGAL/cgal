#include <CGAL/ImageIO.h>
#include <iostream>
#include <string>

int main(int argc, char **argv) {
  if (argc != 3) {
    std::cerr << "Usage: slice_size <input> <output>\n";
    return argc != 1;
  }
  _image *image = ::_readImage(argv[1]);
  if (!image)
    return 2;
  auto *new_image = ::_createImage(image->xdim, image->ydim, image->zdim / 2 + 1, 1,
                                   image->vx, image->vy, image->vz*2, image->wdim,
                                   image->wordKind, image->sign);
  const auto* const data = static_cast<char*>(image->data);
  auto* new_data = static_cast<char*>(new_image->data);
  const auto slice_size = image->wdim * image->xdim * image->ydim;
  for (auto k = 0ul; k < image->zdim; k+=2) {
    auto pos = data + slice_size * k;
    new_data = std::copy(pos, pos + slice_size, new_data);
  }
  auto r = ::_writeImage(new_image, argv[2]);
  if(r != ImageIO_NO_ERROR) return 3;
  ::_freeImage(image);
  ::_freeImage(new_image);
}
