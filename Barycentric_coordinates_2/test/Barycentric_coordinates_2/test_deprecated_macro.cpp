#include <iostream>
#include <CGAL/config.h>


template <class T>
struct
CGAL_DEPRECATED_MSG("Toto is deprecated")
 Toto {
  T t;
};


int main()
{
  std::cout << "done" << std::endl;
  return 0;
}
