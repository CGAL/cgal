#include <iostream>
#include <cassert>
#include <CGAL/algorithm.h>

template <int pivot_, int expected_value>
int test() {
  int nbs[10] = { 2, 6, 1, 0, 3, 4, 6, 9, 8, 1};

  int pivot = pivot_;

  std::less<int> cmp;
  CGAL::nth_element(nbs, nbs+pivot, nbs+sizeof(nbs) / sizeof(int), cmp);
  std::cerr << "After nth_element, nbs["<< pivot << "] = " << nbs[pivot] << std::endl;
  assert(nbs[pivot]==expected_value);
  for(std::size_t i = 0; i < sizeof(nbs) / sizeof(int); ++i) {
    std::cerr << " " << nbs[i];
  }
  std::cerr << std::endl;
  for(pivot=0; pivot < 10; ++pivot)
    CGAL::nth_element(nbs, nbs+pivot, nbs+sizeof(nbs) / sizeof(int), cmp);
  std::cerr << "After sort:\n";
  for(std::size_t i = 0; i < sizeof(nbs) / sizeof(int); ++i) {
    std::cerr << " " << nbs[i];
  }
  std::cerr << std::endl;
  for(std::size_t i = 1; i < sizeof(nbs) / sizeof(int); ++i) {
    assert(nbs[i]>=nbs[i-1]);
  }
  return 0;
}

int main() {
  test<0, 0>();
  test<1, 1>();
  test<2, 1>();
  test<3, 2>();
  test<4, 3>();
  test<5, 4>();
  test<6, 6>();
  test<7, 6>();
  test<8, 8>();
  test<9, 9>();
}
