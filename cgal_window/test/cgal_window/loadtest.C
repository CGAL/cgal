#include <CGAL/basic.h>

#if defined(CGAL_USE_CGAL_WINDOW)

#include <CGAL/LEDA/file_panel.h>
#include <CGAL/LEDA/file.h>
#include <fstream>

using namespace std;

void load_handler(string fn)
{
 std::cout << "load_handler got parameter " << fn << "\n"; 
 
 if (CGAL::is_file(fn)) {
   std::cout << fn << " is a file!\n";
   
   std::ifstream I(fn.c_str());
   
   char ch;
   
   while (I.get(ch)) cout << ch;
   
   I.close();
 }
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
 P.set_load_handler(load_handler);
 P.set_cancel_handler(dummy); 
 P.set_pattern("C Files (*.C)","*.C");
 

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
