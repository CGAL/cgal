// Use something defined not in headers but in the CGAL library to test that is was indeed properly built and linked to,

#include <CGAL/ImageIO.h>

int main()
{

  volatile _image* i = _initImage();

  bool ok = (i != 0);
  _freeImage(const_cast<_image*>(i));
  return ok ? 0 : 1;
}
