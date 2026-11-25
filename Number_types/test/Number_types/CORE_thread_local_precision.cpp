// Tests that the precisions of CORE (the default relative and absolute precisions, and those of BigFloat division and
// square root) are thread-local: a precision set in one thread does not affect other threads, each thread starts with
// the default precisions, and threads that compute concurrently, each with precisions of its own, obtain the results
// that they obtain sequentially. (The bigfloat intervals of the algebraic kernels set the precisions temporarily;
// shared precisions made concurrent computations use the precisions of each other.)

#include <iostream>

#include <CGAL/config.h>

#if ! defined(CGAL_USE_CORE) || ! defined(CGAL_HAS_THREADS)

int main() {
  std::cout << "This test requires CORE and threads; it is skipped." << std::endl;
  return 0;
}

#else

#include <thread>
#include <vector>

#include <CGAL/CORE_BigFloat.h>

using BigFloat = CORE::BigFloat;
using BFI = CGAL::Bigfloat_interval_traits<BigFloat>::Type;

/*! checks that the precisions of the calling thread are the given ones.
 */
bool check(long rel, long div, const char* where) {
  const long r = CORE::get_static_defRelPrec().toLong();
  const long d = CORE::get_static_defBFdivRelPrec().toLong();
  if (r == rel && d == div) return true;
  std::cerr << "Error: " << where << ": relative precision " << r << " (expected " << rel << "), division precision "
            << d << " (expected " << div << ")" << std::endl;
  return false;
}

/*! computes quotients and square roots at the given precision (set as the bigfloat intervals set it), repeatedly.
 */
std::vector<BigFloat> compute(long prec, int rounds) {
  std::vector<BigFloat> res;
  for (int i = 0; i < rounds; ++i) {
    const long old = CGAL::set_precision(BFI(), prec);
    BigFloat x(1 + i);
    BigFloat y(3 + 2 * i);
    res.push_back(x / y);
    res.push_back(CORE::sqrt(x / y));
    CGAL::set_precision(BFI(), old);
  }
  return res;
}

bool same(const std::vector<BigFloat>& a, const std::vector<BigFloat>& b) {
  if (a.size() != b.size()) return false;
  for (std::size_t i = 0; i < a.size(); ++i) {
    // the same mantissa, error and exponent
    if (a[i].m() != b[i].m() || a[i].err() != b[i].err() || a[i].exp() != b[i].exp()) return false;
  }
  return true;
}

int main() {
  bool ok = true;
  const long def_rel = CORE::get_static_defRelPrec().toLong();
  const long def_div = CORE::get_static_defBFdivRelPrec().toLong();

  // a precision set in one thread does not affect other threads, which start with the default precisions
  CGAL::set_precision(BFI(), 300);
  ok &= check(300, 300, "the main thread, after setting the precision");
  std::thread t([&ok, def_rel, def_div]() {
    ok &= check(def_rel, def_div, "a new thread");
    CGAL::set_precision(BFI(), 500);
    ok &= check(500, 500, "a new thread, after setting the precision");
  });
  t.join();
  ok &= check(300, 300, "the main thread, after another thread set the precision");
  CGAL::set_precision(BFI(), def_rel);

  // concurrent computations, each with precisions of its own, obtain their sequential results
  const int num_threads = 8;
  const int rounds = 2000;
  std::vector<std::vector<BigFloat>> expected(num_threads);
  for (int i = 0; i < num_threads; ++i) expected[i] = compute(40 + 60 * i, rounds);
  std::vector<std::vector<BigFloat>> results(num_threads);
  std::vector<std::thread> threads;
  for (int i = 0; i < num_threads; ++i)
    threads.emplace_back([&results, i, rounds]() { results[i] = compute(40 + 60 * i, rounds); });
  for (auto& th : threads) th.join();
  for (int i = 0; i < num_threads; ++i) {
    if (same(results[i], expected[i])) continue;
    std::cerr << "Error: thread " << i << " obtains results that differ from the sequential ones" << std::endl;
    ok = false;
  }

  std::cout << (ok ? "passed" : "failed") << std::endl;
  return ok ? 0 : 1;
}

#endif
