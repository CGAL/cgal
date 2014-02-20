#include <cstdio>
#include <iterator>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <list>

#include <boost/timer.hpp>
#include <boost/lexical_cast.hpp>

#include <CGAL/algorithm.h>
#include <CGAL/compiler_config.h>

//class with non-trivial copy ctor
class non_trivial_cctor
{
public:
  non_trivial_cctor() { data = new char[n]; std::fill_n(data, n, 0); }
  non_trivial_cctor(const non_trivial_cctor& x) { data = new char[n]; std::copy(x.data, x.data + n, this->data);  }
  ~non_trivial_cctor() { delete[] data; }
private:
  static const int n;
  char* data;
};

const int non_trivial_cctor::n = 64;

int format_output(const char* lib, const char* container, const char* to, int n, float time) {
  return std::printf("| %s || %s || %s || %d || %.4fM items/sec\n", lib, container, to, n, time);
}

struct std_tag {};
struct cgal_tag {};

#ifndef CGAL_CFG_NO_CPP0X_COPY_N
template <typename ForwardIterator, typename Size, typename OutputIterator>
inline double test(ForwardIterator it, Size n, OutputIterator result, int repeats, std_tag) {
  boost::timer timer;
  timer.restart();
  for (int i = 0; i < repeats; ++i) { std::copy_n(it, n, result); }
  return (double)n*repeats/timer.elapsed()/1.0E6;
}
#endif

template <typename ForwardIterator, typename Size, typename OutputIterator>
inline double test(ForwardIterator it, Size n, OutputIterator result, int repeats, cgal_tag) {
  boost::timer timer;
  timer.restart();
  for (int i = 0; i < repeats; ++i) { CGAL::copy_n(it, n, result); }
  return (double)n*repeats/timer.elapsed()/1.0E6;
}

int main(int argc, char* argv[]) {
  int n = 20000000;
  int repeats = 20;

  if(argc > 1)
    n = boost::lexical_cast<int>(argv[1]);

  if(argc > 2)
    repeats = boost::lexical_cast<int>(argv[2]);

  typedef std::vector<int> vector;
  typedef std::list<int> list;
  typedef std::vector<non_trivial_cctor> vector2;

  vector v(n);

  typedef int* copy_mem;
  typedef non_trivial_cctor* copy_mem2;

  copy_mem copy_m = new int[n];

  //wiki markup header
  std::cout << 
    "{| \n"
    "! Library !! From Container !! To !! #Elements !! items/sec \n"
    "|- \n";
  float item_sec;

#ifndef CGAL_CFG_NO_CPP0X_COPY_N
  item_sec = test(v.begin(), n, copy_m, repeats, std_tag());
  format_output("stdlib", "vector<int>", "int*", n, item_sec);
  std::cout << "|- \n";
#endif

  item_sec = test(v.begin(), n, copy_m, repeats, cgal_tag());
  format_output("CGAL", "vector<int>", "int*", n, item_sec);
  std::cout << "|- \n";

  //purge
  v.clear();
  std::fill_n(copy_m, n, 0);

  list l(n);

#ifndef CGAL_CFG_NO_CPP0X_COPY_N
  item_sec = test(l.begin(), n, copy_m, repeats, std_tag());
  format_output("stdlib", "list<int>", "int*", n, item_sec);
  std::cout << "|- \n";
#endif

  item_sec = test(l.begin(), n, copy_m, repeats, cgal_tag());
  format_output("CGAL", "list<int>", "int*", n, item_sec);
  std::cout << "|- \n";

  delete[] copy_m;

  vector2 v2(n);
  copy_mem2 copy_m2 = new non_trivial_cctor[n];
  
#ifndef CGAL_CFG_NO_CPP0X_COPY_N
  item_sec = test(v2.begin(), n, copy_m2, repeats, std_tag());
  format_output("stdlib", "vector<non_trivial_cctor>", "non_trivial_cctor*", n, item_sec);
  std::cout << "|- \n";
#endif

  item_sec = test(v2.begin(), n, copy_m2, repeats, cgal_tag());
  format_output("CGAL", "vector<non_trivial_cctor>", "non_trivial_cctor*", n, item_sec);

  //wiki markup footer
  std::cout << "|}" << std::endl;
  


  return EXIT_SUCCESS;
}
