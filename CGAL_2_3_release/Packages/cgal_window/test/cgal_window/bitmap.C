#include <CGAL/basic.h>

#if defined(CGAL_USE_CGAL_WINDOW)

#include <CGAL/LEDA/window.h>
#include <CGAL/LEDA/bitmaps/button32/voro.xbm>
#include <list>

using namespace std;

/*
#if defined(__KCC)
list<CGAL::window_point> L1;
list<std::string>        L2;
#endif
*/

int main(int argc, char *argv[])
{ 
  if (argc >= 2) { return 0; }

  CGAL::window W(400,400);
  W.set_bg_color(CGAL::yellow);
  W.display();

  // construct bitmap from the bitmap data in 
  // <CGAL/LEDA/bitmaps/button32/voro.xbm>

  char* bm = W.create_bitmap(32,
                             32, 
                             voro_bits);

  // copy copies of bm into the window

  CGAL::window_point p;
  while (W.read_mouse(p) != MOUSE_BUTTON(3)) {
    W.put_bitmap(p.xcoord(),p.ycoord(),bm,CGAL::blue);
  }

  W.del_bitmap(bm);

  //W.screenshot("bitmap");
  return 0;
}

#else
#include <iostream>

int main()
{
 std::cout << "CGAL::window is not available !\n";

 return 0;
}

#endif


