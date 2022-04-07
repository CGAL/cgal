#include <iostream>

// This file displays the value of `__cplusplus` on `stderr`, and then
// return 0 if and only if C++ `thread_local` feature can be used.
// The tests are the same as for the definition of
// `CGAL_CAN_USE_CXX11_THREAD_LOCAL` in `<CGAL/config.h>`.

#ifndef __has_feature
  #define __has_feature(x) 0  // Compatibility with non-clang compilers.
#endif

int main() {
#if ! defined(__cplusplus)
  std::cout << "Undefined";
  return 1;
#else
  std::cout << __cplusplus;
#endif

#if __has_feature(cxx_thread_local) || \
    ( __GNUC__ > 0 && __cplusplus >= 201103L )
  return 0;
#else
  return 1;
#endif
}
