// Use something defined not in headers but in the CGAL library to test that is was indeed properly built and linked to,

#include <CGAL/IO/Color.h>

int main()
{
  volatile const CGAL::Color* c = &CGAL::BLACK;
  
  return 0;
}
