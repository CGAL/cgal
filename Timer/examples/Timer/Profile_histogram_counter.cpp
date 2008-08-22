#define CGAL_PROFILE

#include <CGAL/Profile_counter.h>

int main()
{
  for (int i=0; i<10; ++i)
  {
    for (int j=0; j<i; ++j)
      CGAL_HISTOGRAM_PROFILER("[iterations of the for-loop]", j);
  }
  return 0;
}
