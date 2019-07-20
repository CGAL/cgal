#define CGAL_PROFILE

#include <CGAL/Profile_counter.h>

int main()
{
  for (int i=0; i<10000; ++i)
  {
    CGAL_PROFILER("iterations of the for-loop");
  }
  return 0;
}
