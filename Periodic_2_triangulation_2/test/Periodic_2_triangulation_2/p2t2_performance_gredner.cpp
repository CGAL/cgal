// #define CGAL_PROFILE

#include "./types.h"
#include <CGAL/Delaunay_triangulation_2.h>

#include <iostream>
#include <fstream>

#include <ctime>
#include <algorithm>

const bool pre_run = true;
const bool do_remove = true;
const int n_runs = 3;

template <class T>
void test(const std::vector<Point> &input, T &t)
{
  t.insert(input.begin(), input.end());

  if (do_remove) {
    std::vector<typename T::Vertex_handle> vhs;
    for (typename T::Vertex_iterator it = t.vertices_begin(); it != t.vertices_end(); ++it) {
      vhs.push_back(it);
    }
    
    std::random_shuffle(vhs.begin(), vhs.end());
    vhs.resize(vhs.size()/2);
    for (size_t i=0; i<vhs.size(); ++i)
      t.remove(vhs[i]);
  }
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

  if (true) {
    if (pre_run) {
      Delaunay_triangulation_2<Gt> t;
      test(pts, t);
    }

    std::clock_t total_start = std::clock();
    for (int i=0; i<n_runs; ++i) {
      Delaunay_triangulation_2<Gt> t;
      test(pts, t);
    }
    double total_time = (std::clock()-total_start)/(double)CLOCKS_PER_SEC;

    std::cout << "Euclidean space, " << filename << ", " << total_time << std::endl;
  }

  if (true) {
    if (pre_run) {
      Periodic_2_Delaunay_triangulation_2<Gt> t(Iso_rectangle(0,0,domain[0],domain[1]));
      test(pts, t);
    }

    std::clock_t total_start = std::clock();
    for (int i=0; i<n_runs; ++i) {
      Periodic_2_Delaunay_triangulation_2<Gt> t(Iso_rectangle(0,0,domain[0],domain[1]));
      test(pts, t);
    }
    double total_time = (std::clock()-total_start)/(double)CLOCKS_PER_SEC;

    std::cout << "Periodic  space, " << filename << ", " << total_time << std::endl;

  }

  return 0;
}
