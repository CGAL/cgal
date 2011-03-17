#include <ctime>
#include <cstdlib>
#include <iostream>

#define N 500

template <bool adds>
void test_1 () {
  double v[N];
  std::clock_t start = std::clock();
  for(int i = 0; i < N; ++i) {
    v[i] = (0. + std::rand()) / RAND_MAX;
  }
  std::clock_t end = std::clock();
  // std::cout << "init time = " << (0. + end - start) / CLOCKS_PER_SEC << " seconds"
  //           << " (" << (end - start) << " clock ticks)" << std::endl;
  int counter = 0;
  start = std::clock();
  for(int k = 0; k < 10000; ++k) {
    for(int i = 0; i < N; ++i) {
      for(int j = i + 1; j < N; ++j) {
        if(adds)
          counter += (v[i] < v[j]);
        else
          if(v[i] < v[j]) counter = 1;
      }
    }
  }
  end = std::clock();
  std::cout << "counter = " << counter << std::endl;
  // std::cout << "% (should be around 0.5) = " << 0.0002 * (0. + counter) / N / (N-1) << std::endl;
  std::cout << "time = " << (0. + end - start) / CLOCKS_PER_SEC << " seconds"
            << " (" << (end - start) << " clock ticks)" << std::endl;
}

template <bool adds>
void test_2() {
  volatile double d1 = 1.;
  volatile double d2 = 2.;
  const unsigned long int iter = -1;
  std::clock_t start = std::clock();
  int result = 0;
  for(unsigned long int i = 0; i < iter; ++i) {
    if(adds) {
      result += (d1 < d2);
    } else {
      if(d1 < d2) result = 1;
    }
  }
  std::clock_t end = std::clock();
  std::cout << "result = " << result << std::endl;
  //           << "nb iter = " << iter << std::endl;
  std::cout << "time = " << (0. + end - start) / CLOCKS_PER_SEC << " seconds"
            << " (" << (end - start) << " clock ticks)" << std::endl;
}

int main() {
  test_1<false>();
  test_1<true>();
  test_2<false>();
  test_2<true>();
}
