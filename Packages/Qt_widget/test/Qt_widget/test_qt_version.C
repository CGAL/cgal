#ifndef CGAL_USE_QT
  #include <iostream>
  int main(int, char*)
  {
    std::cout << "This platform does not have QT installed.";
    std::cout << std::endl;
    return 0;
  }
#else
  #include <iostream>
  #include <qglobal.h>
  int main(char*, char**)
  {
    std::cout << QT_VERSION << std::endl;
    return 0;
  }
#endif
