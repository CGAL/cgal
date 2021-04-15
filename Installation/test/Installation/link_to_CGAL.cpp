// Use something defined not in headers but in the CGAL library to test that is was indeed properly built and linked to,

#include <CGAL/Random.h>

int main()
{
  volatile const CGAL::Random* r = &CGAL::get_default_random();
  return (r != 0) ? 0 : 1;
}
