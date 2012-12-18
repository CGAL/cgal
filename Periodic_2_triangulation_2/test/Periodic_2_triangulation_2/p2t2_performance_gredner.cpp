// #define CGAL_PROFILE

#include "./types.h"
#include <CGAL/Delaunay_triangulation_2.h>

#include <iostream>
#include <fstream>

#include <ctime>
#include <algorithm>

template <class T>
double test(const std::vector<Point> &input, T &t)
{
  std::clock_t total_start = std::clock();
  t.insert(input.begin(), input.end());
  return (std::clock()-total_start)/(double)CLOCKS_PER_SEC;
}

int main(int argc, char * argv[]) {
  srand(42);
  const char *filename = "/home/nico/Code/periodic_data_sets/512000_000.dat";
  if (argc == 2) {
    filename = argv[1];
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
    while (coords[0] < 0) coords[0] += domain[0];
    while (coords[1] < 0) coords[1] += domain[1];
    while (coords[0] >= domain[0]) coords[0] -= domain[0];
    while (coords[1] >= domain[1]) coords[1] -= domain[1];

    pts.push_back(Point(coords[0], coords[1]));
  }

  if (false) {
    // Warming up ...
    std::random_shuffle(pts.begin(), pts.end());
    Delaunay_triangulation_2<Gt> t;
    test(pts, t);

    std::random_shuffle(pts.begin(), pts.end());
    Periodic_2_Delaunay_triangulation_2<Gt> t2(Iso_rectangle(0,0,domain[0],domain[1]));
    test(pts, t2);
  }
  if (true) {
    std::random_shuffle(pts.begin(), pts.end());
    Delaunay_triangulation_2<Gt> t;

    std::cout << "Euclidean space, " << filename << ", ";
    std::cout << test(pts, t) << std::endl;
  }
  if (true) {
    std::random_shuffle(pts.begin(), pts.end());
    Periodic_2_Delaunay_triangulation_2<Gt> t(Iso_rectangle(0,0,domain[0],domain[1]));

    std::cout << "Periodic  space, " << filename << ", ";
    std::cout << test(pts, t) << std::endl;
  }

  return 0;
}
