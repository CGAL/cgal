#include <CGAL/basic.h>

#if defined(CGAL_USE_CGAL_WINDOW)

#include <CGAL/LEDA/window.h>

int main(int argc, char *argv[])
{
  if (argc >= 2) { return 0; }

  int count = 0;

  CGAL::window W(400,100);
  W.set_item_width(300);
  W.int_item("progress",count,0,1000);
  W.display(CGAL::window::center,CGAL::window::center);

  for(;;)
  { count = 0;
    while (count < 1000)
    { W.redraw_panel();
      W.flush();
      count++;
     }
    if (W.read_mouse() == MOUSE_BUTTON(3)) break;
   }

  W.screenshot("progress");
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

