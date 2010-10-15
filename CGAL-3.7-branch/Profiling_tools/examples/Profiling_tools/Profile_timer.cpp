#define CGAL_PROFILE

#include <CGAL/Profile_timer.h>

void f()
{
  CGAL_TIME_PROFILER("seconds spent in this function (total time)");
  // do something
  double d = 1+1;
  (void) d;
}

int main()
{
  for (int i=0; i<100000; ++i)
    f();
  return 0;
}
