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
#include <CGAL/Delaunay_triangulation_2.h>

#include <iostream>
#include <fstream>

#include <ctime>
#include <algorithm>

const bool pre_run = false;
const bool do_remove = true;
const int n_runs = 2;

bool load_data(const char *filename, Iso_rectangle &domain, std::vector<Point> &pts)
{
  std::ifstream file (filename, std::ios::in | std::ios::binary);
  if (!file.is_open()) {
    std::cout << "Couldn't load " << filename << std::endl;
    return false;
  }

  float dom[2];
  file.read((char *)&dom[0], 2 * sizeof(float));
  domain = Iso_rectangle(0, 0, dom[0], dom[1]);

  float coords[2];
  while (!file.eof()) {
    file.read((char *)&coords[0], 2 * sizeof(float));
    while (coords[0] < 0) coords[0] += dom[0];
    while (coords[1] < 0) coords[1] += dom[1];
    while (coords[0] >= dom[0]) coords[0] -= dom[0];
    while (coords[1] >= dom[1]) coords[1] -= dom[1];
    
    pts.push_back(Point(coords[0], coords[1]));
  }
  return true;
}

template <class T>
void test(const std::vector<Point> &input, T &t)
{
  t.insert(input.begin(), input.end(), true);

  if (do_remove)
    {
      std::vector<typename T::Vertex_handle> vhs;
      for (typename T::Vertex_iterator it = t.vertices_begin(); it != t.vertices_end(); ++it)
        {
          vhs.push_back(it);
        }

      std::random_shuffle(vhs.begin(), vhs.end());
      vhs.resize(vhs.size() / 2);
      for (size_t i = 0; i < vhs.size(); ++i)
        t.remove(vhs[i]);
    }
}

template <class T>
bool test(const char *filename) {
  Iso_rectangle domain;
  std::vector<Point> pts;
  if (!load_data(filename, domain, pts))
    return false;

  for (int run = 0; run < 3; ++run)
    {
      // if (true) {
      //   if (pre_run) {
      //     Delaunay_triangulation_2<Gt> t;
      //     test(pts, t);
      //   }

      //   std::clock_t total_start = std::clock();
      //   for (int i=0; i<n_runs; ++i) {
      //     Delaunay_triangulation_2<Gt> t;
      //     test(pts, t);
      //   }
      //   double total_time = (std::clock()-total_start)/(double)CLOCKS_PER_SEC;

      //   std::cout << "Euclidean space, " << filename << ", " << total_time << std::endl;
      // }

      if (true)
        {
          if (pre_run)
            {
              T t(domain);
              test(pts, t);
            }

          std::clock_t total_start = std::clock();
          for (int i = 0; i < n_runs; ++i)
            {
              T t(domain);
              test(pts, t);
            }
          double total_time = (std::clock() - total_start) / (double)CLOCKS_PER_SEC;

          std::cout << "Periodic  space, " << filename << ", " << total_time << std::endl;
        }
    }
  return true;
}

int main(int argc, char * argv[])
{
  typedef Periodic_2_Delaunay_triangulation_2<Gt> T;
  srand(42);

  int result = 0;
  if (argc == 2) {
    if (!test<T>(argv[1])) result++;
  }
  else {
    if (!test<T>("data/gredner000.dat")) result++;
    if (!test<T>("data/gredner050.dat")) result++;
    if (!test<T>("data/gredner100.dat")) result++;
    if (!test<T>("data/gredner150.dat")) result++;
    if (!test<T>("data/gredner200.dat")) result++;
  }

  return result;
}

#endif
