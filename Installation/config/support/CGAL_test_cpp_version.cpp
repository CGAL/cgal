#include <iostream>

#ifndef __has_feature
  #define __has_feature(x) 0  // Compatibility with non-clang compilers.
#endif

int main() {
#if ! defined(__cplusplus)
  std::cout << "Undefined";
  return 1;
#else
  std::cout << __cplusplus;
#  if __has_feature(cxx_thread_local) || __cplusplus >= 201103L
  return 0;
#  else
  return 1;
#  endif
#endif
}
