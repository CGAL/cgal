#define CGAL_PROFILE

#include <CGAL/Profile_timer.h>

int main()
{
  CGAL_TIME_PROFILER("seconds spent in this for loop");
  for (int i=0; i<10; ++i)
  {
    // do something
    double d = 1+1;
    (void) d;
  }
  return 0;
}
