#include <ctime>
#include <cstdlib>
#include <iostream>

#define N 500

int main() {
  double v[N];
  std::clock_t start = std::clock();
  for(int i = 0; i < N; ++i) {
    v[i] = (0. + std::rand()) / RAND_MAX;
  }
  std::clock_t end = std::clock();
  std::cout << "init time = " << (0. + end - start) / CLOCKS_PER_SEC << " seconds"
            << " (" << (end - start) << " clock ticks)" << std::endl;
  int counter = 0;
  start = std::clock();
  for(int k = 0; k < 10000; ++k) {
    for(int i = 0; i < N; ++i) {
      for(int j = i + 1; j < N; ++j) {
        counter += (v[i] < v[j]);
      }
    }
  }
  end = std::clock();
  std::cout << "counter = " << counter << std::endl;
  std::cout << "% (should be around 0.5) = " << 0.0002 * (0. + counter) / N / (N-1) << std::endl;
  std::cout << "time = " << (0. + end - start) / CLOCKS_PER_SEC << " seconds"
            << " (" << (end - start) << " clock ticks)" << std::endl;
  
}
