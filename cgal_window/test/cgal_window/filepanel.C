#include <CGAL/basic.h>

#if defined(CGAL_USE_CGAL_WINDOW)

#include <CGAL/LEDA/file_panel.h>

using namespace std;

void save_handler(string fn)
{
 std::cout << "save_handler got parameter " << fn << "\n"; 
}

void dummy(string fn)
{
 std::cout << "dummy got parameter " << fn << "\n"; 
}

int main(int argc, char *argv[])
{
 if (argc >= 2) { return 0; }

 CGAL::window Win(600,600);
 Win.init(0,600,0);
 
 Win.display();
 
 // open a file panel ...
 std::string fname("test"), dname(".");
 
 CGAL::file_panel P(Win, fname, dname);
 P.set_save_handler(save_handler);
 P.set_cancel_handler(dummy); 
 P.set_pattern("C Files (*.C)","*.C");

 P.open();
 
 std::cout << fname << "  " << dname << "\n";

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
