// Use something defined not in headers but in the CGAL library to test that is was indeed properly built and linked to,

#include <CGAL/CORE/BigFloat.h>

int main()
{
  CORE::BigFloat n(1.2);
  
  volatile double d = n.doubleValue();
  
  return (d > 1) ? 0 : 1;
}
