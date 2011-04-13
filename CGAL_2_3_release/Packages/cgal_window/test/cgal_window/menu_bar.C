#include <CGAL/basic.h>

#if defined(CGAL_USE_CGAL_WINDOW)

#include <CGAL/LEDA/window.h>
#include <CGAL/LEDA/file.h>


int main(int argc, char *argv[]) 
{
  if (argc >= 2) { return 0; }

  CGAL::menu M;
  M.button("button 1");
  M.button("button 2");
  M.button("button 3");
  M.button("button 4");
  M.button("button 5");

  CGAL::window W(400,300,"Menu Demo");

  W.button("File",M);
  W.button("Edit",M);
  W.button("Help",M);
  W.button("exit");

  W.make_menu_bar();

  W.display();

  W.read_mouse();

  W.screenshot("menu_bar");
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
