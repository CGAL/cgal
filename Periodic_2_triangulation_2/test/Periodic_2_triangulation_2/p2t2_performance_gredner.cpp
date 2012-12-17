#define CGAL_PROFILE

#include "./types.h"
#include <CGAL/Delaunay_triangulation_2.h>

#include <iostream>
#include <fstream>

#include <ctime>
#include <algorithm>

template <class T>
void test(const std::vector<Point> &input, 
          T &t,
          std::vector< std::pair<int, double> > &timings)
{
  const size_t interval = 100;
  const size_t n_pts = input.size();
  timings.reserve(n_pts / interval);

  if (true) {
    std::clock_t total_start = std::clock();
    t.insert(input.begin(), input.end());
    timings.push_back(std::make_pair(t.number_of_vertices(),
                                     (std::clock()-total_start)/(double)CLOCKS_PER_SEC));
  } else {
    std::clock_t total_start = std::clock();
    for (size_t i=0; i<input.size(); ++i) {
      // std::cout << i << std::endl;
      t.insert(input[i]);
      CGAL_assertion(t.is_valid());
    
      if (i%1000 == 0) {
        std::cout << t.number_of_vertices() << " \t"
                  << ((std::clock()-total_start)/(double)CLOCKS_PER_SEC) << std::endl;
        timings.push_back(std::make_pair(t.number_of_vertices(),
                                         (std::clock()-total_start)/(double)CLOCKS_PER_SEC));
      }
    }
  }
}

int main(int argc, char * argv[]) {
  srand(42);
  bool periodic_triangulation = true;
  const char *filename = "/home/nico/Code/periodic_data_sets/512000_000.dat";
  if (argc == 3) {
    periodic_triangulation = (argv[1][0] == 'p');
    filename = argv[2];
  }
  std::cout << "testing file: " << filename << std::endl;
  std::ifstream file (filename, std::ios::in|std::ios::binary);
  if (!file.is_open()) return 0;

  float domain[2];
  file.read((char *)&domain[0], 2 * sizeof(float));

  std::vector<Point> pts;
  float coords[2];
  while (!file.eof()) {
    file.read((char *)&coords[0], 2 * sizeof(float));
    pts.push_back(Point(coords[0], coords[1]));
  }
  std::random_shuffle(pts.begin(), pts.end());

  std::vector< std::pair<int, double> > timings;
  if (periodic_triangulation) {
    Periodic_2_Delaunay_triangulation_2<Gt> t(Iso_rectangle(0,0,domain[0],domain[1]));
    test(pts, t, timings);

    std::cout << "Periodic  space, " << filename << ", ";
  } else {
    Delaunay_triangulation_2<Gt> t;
    test(pts, t, timings);

    std::cout << "Euclidean space, " << filename << ", ";
  }

  std::cout << timings.back().first << ", " << timings.back().second << std::endl;
  return 0;
}
