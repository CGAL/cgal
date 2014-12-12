#include <vector>
#include <deque>
#include <CGAL/Memory_sizer.h>
#include <CGAL/array.h>
#include <CGAL/Block_list.h>
#include <CGAL/Timer.h>
#include <boost/container/deque.hpp>

const int N = 10000000;

struct Leaf {
  bool b;
  int i;
  double d1;
};


struct Internal {
  bool b;
  int i;
  double d1, d2, d3, d4,d5;
};

template <typename Container, typename T>
void bench_container()
{
  std::cout << typeid(Container).name() << std::endl;

  CGAL::Memory_sizer ms;
  CGAL::Timer t;
  t.start();
  
  T element;
  Container c;
  for(int i=0; i <N ; i++){
    c.push_back(element);
  }
  t.stop();
  std::cout << t.time() << "sec"<< std::endl;
  std::cout << "virtual : " << ms.virtual_size() << std::endl
            << "resident: " << ms.resident_size() << std::endl << std::endl;

}

template <typename T>
void
bench()
{
  std::cout << "sizeof(" << typeid(T).name() << ") = " << sizeof(T) <<std::endl;
  //bench_container<CGAL::Block_list<T, 256>,T>();
  bench_container<std::deque<T>,T>();
  bench_container<boost::container::deque<T>,T>();
  //bench_container<std::vector<T>,T>();
}

int main()
{
  bench<double>();
  bench<Leaf>();
  bench<Internal>();
  return 0;
}
