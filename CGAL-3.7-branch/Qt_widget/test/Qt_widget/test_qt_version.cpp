#include <CGAL/basic.h>
  #include <iostream>
  #include <qglobal.h>
  int main(int, char**)
  {
    std::cout << "QT_VERSION= " << QT_VERSION << std::endl
	      << "QT_VERSION_STR= " << QT_VERSION_STR << std::endl;
    return 0;
  }
