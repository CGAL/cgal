#ifndef TIMER_HXX
#define TIMER_HXX

#include <unistd.h>
#include <sys/time.h>

struct Timer
{
 public:
  struct timeval t1,t2;

  double t;

  void start()
    {
      gettimeofday (&t1, NULL);
    }

  void stop()
    {
      gettimeofday (&t2, NULL);
      t = (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec) / 1000000.0;
    }

  void reset()
    {
      t     = 0.0;
    }
};

#endif

