#include <CGAL/basic.h>

#if defined(CGAL_USE_CGAL_WINDOW)

#include <CGAL/LEDA/window.h>
#include <cmath>

int main(int argc, char *argv[])
{
  if (argc >= 2) { return 0; }

  CGAL::window W(450,500,"Event Demo");
  W.display();

  W.start_buffering();

  for(;;)
  { 
    // read the first corner p of the rectangle
    // terminate if the right button was clicked 

    CGAL::window_point p;
    if (W.read_mouse(p) == MOUSE_BUTTON(3)) break;  
  
    // draw rectangle from p to current position while button down 

    int  val;
    double x,y;

    char* win_buf = W.get_window_pixrect();

    while (W.read_event(val,x,y) != CGAL::button_release_event) 
    { CGAL::window_point q(x,y);
      W.put_pixrect(win_buf);
      W.draw_box(p,q,CGAL::yellow);
      W.draw_rectangle(p,q,CGAL::black);
      W.flush_buffer();
     }

    W.del_pixrect(win_buf);
  }

  W.stop_buffering();

  W.screenshot("event");
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

