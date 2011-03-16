#include <ctime>
#include <cstdlib>
#include <vector>
#include <iostream>

int main() {
  const int N = 50000;
  std::vector<double> v(N);
  for(int i = 0; i < N; ++i) {
    v[i] = (0. + std::rand()) / RAND_MAX;
  }
  int counter = 0;
  std::clock_t start = std::clock();
  for(int i = 0; i < N; ++i) {
    for(int j = i + 1; j < N; ++j) {
      counter += (v[i] < v[j]);
    }
  }
  std::clock_t end = std::clock();
  std::cout << "counter = " << counter << std::endl;
  std::cout << "% = " << 2 * (0. + counter) / N / (N-1) << std::endl;
  std::cout << "time = " << (0. + end - start) / CLOCKS_PER_SEC << " seconds\n";
}
