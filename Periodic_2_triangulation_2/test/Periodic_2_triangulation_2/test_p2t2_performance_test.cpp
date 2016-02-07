// Author(s)     : Nico Kruithof  <Nico@nghk.nl>

#include <CGAL/config.h>

#ifndef CGAL_NDEBUG

#include <iostream>
int main()
{
  std::cerr << "No performance test for a Debug build" << std::endl;
  return 0;
}

#else

#include "./types.h"
#include <CGAL/Triangulation_2.h>
#include <ctime>

typedef Triangulation_2<Gt>              EuclideanTriangulation;

const int N_RUNS = 2;
const int N_PTS  = 50000;
//#define PERIODIC


/* Timings for N_PTS=50000, iMac, 2nd run, real time:
 * Normal triangulation:                      6.29806s
 * Normal triangulation with periodic copies: 31.6s
 * Periodic triangulation:
 * r58602 : 153.909s
 * r58618 : 154.491s
 */

double total_time = 0.0;
double locate_time = 0.0;
double insert_time = 0.0;
double periodic_locate_time = 0.0;
double periodic_insert_time = 0.0;

int main()
{
  Random random(1284141159);
  Random_points_in_square g(0.495, random);
  Vector midpoint(0.5, 0.5);

#if 0
  // Should take 5 seconds on the iMac
  std::vector<Point> pts;
  pts.resize(500000);
  for (size_t i = 0; i < pts.size(); ++i) pts[i] = *(++g) + midpoint;

  Gt gt;
  for (size_t i = 0; i < pts.size() - 2; ++i) gt.orientation_2_object()(pts[i], pts[i + 1], pts[i + 2]);

  return 0;
#endif

  std::cout << "i" << ", \t"
            << "total time"  << ", \t"
            << "total_time cumm" << ", \t"
            << "locate_time" << ", \t" << "insert_time" << ", \t"
            << "periodic_insert_time" <<  ", \t" << "periodic_insert_time" << std::endl;


  // First run is for heating up the CPU, don't output the stats
  for (int run = 0; run < N_RUNS; ++run)
    {
      // Reset timings
      total_time = 0.0;
      locate_time = 0.0;
      insert_time = 0.0;
      periodic_locate_time = 0.0;
      periodic_insert_time = 0.0;

#ifdef PERIODIC
      Triangulation t;
      const bool insert_periodic_copies = false;
#else
      EuclideanTriangulation t;
      const bool insert_periodic_copies = false;
#endif

      std::clock_t start_time = std::clock();
      // Do one additional point to get the statistics right
      for (int i = 0; i <= N_PTS; ++i)
        {
          Point p = *(++g) + midpoint;
          t.insert(p);

          if (insert_periodic_copies)
            {
              for (int x = 0; x < 3; ++x)
                {
                  for (int y = 0; y < 3; ++y)
                    {
                      if (x + y > 0)
                        {
                          t.insert(p + Vector(x, y));
                        }
                    }
                }
            }

          if (i % 500 == 0)
            {
              std::cout << i << ", \t"
                        << (std::clock() - start_time) / (double)CLOCKS_PER_SEC << ", \t"
                        << total_time << ", \t"
                        << locate_time << ", \t" << insert_time << ", \t"
                        << periodic_insert_time <<  ", \t" << periodic_insert_time << std::endl;
            }
        }
      CGAL_assertion(t.is_valid());

      std::cout << std::endl;
    }

  return 0;
}

#endif
