#include <CGAL/Random.h>
#include <array>
#include <iostream>

int main() {
  CGAL::Random rng;
  std::cout.precision(17);
  std ::cout << 6 << std::endl;
  std::cout << 100000 << std::endl;
  for(int i = 0; i < 100000; ++i) {
    std::array<double, 6> arr;
    double slength = 0;
    for(int i = 0; i < arr.size(); ++i) {
      arr[i] = rng.get_double(-1, 1);
      slength += arr[i] * arr[i];
    }
    slength = std::sqrt(slength);
    for(int i = 0; i < arr.size(); ++i) {
      arr[i] /= slength;
    }
    std::cout  << "6 "  << arr[0] << " " << arr[1] << " " << arr[2] << " " << arr[3] << " " << arr[4]<< " " << arr[5] << std::endl;
  }
  return 0;
}
