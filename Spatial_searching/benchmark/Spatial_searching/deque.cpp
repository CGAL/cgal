#include <vector>
#include <deque>
#include <CGAL/Memory_sizer.h>
#include <CGAL/array.h>
#include <CGAL/Block_list.h>
#include <CGAL/Timer.h>
#include <boost/container/deque.hpp>

#if 1 // 
// leaf
struct X {
  bool b;
  int i;
  double d1;
};
#else
// internal
struct X {
  bool b;
  int i;
  double d1, d2, d3, d4,d5;
};
#endif
int main()
{
  CGAL::Memory_sizer ms;
  CGAL::Timer t;
  t.start();
#if 0
  std::cout <<"blocklist"<< std::endl;
  CGAL::Block_list<X, 256> de;
#elif 0
  std::cout <<"boost::container::deque"<< std::endl;
  boost::container::deque<X> de;
#elif 0
  std::cout <<"std::deque"<< std::endl;
  std::deque<X> de;
#else
  std::cout <<"vector"<< std::endl;
  std::vector<X> de;
  // de.reserve(100000000);
#endif

  for(int i=0; i < 10000000; i++){
    de.push_back(X());
  }
  t.stop();
  std::cout << t.time() << "sec"<< std::endl;
  std::cout << ms.virtual_size() << " " << ms.resident_size() << std::endl;

  return 0;
}
