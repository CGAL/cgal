#include <CGAL/basic.h>

#if defined(CGAL_USE_CGAL_WINDOW)

#include <CGAL/LEDA/window.h>
#include <CGAL/LEDA/bitmaps/button32.h>

using namespace std;

int main(int argc, char *argv[])
{ 
  if (argc >= 2) { return 0; }

  CGAL::panel P("Bitmap Buttons");
  P.buttons_per_line(8);
  P.set_button_space(3);

  for(int i=0; i < num_xbm_button32; i++) 
    P.button(32,32,xbm_button32[i],string(name_xbm_button32[i]));
  
  P.open();
 
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

