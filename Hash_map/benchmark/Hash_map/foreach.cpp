#include <iostream>
#include <vector>

#include<boost/range/iterator_range.hpp> 
#include <boost/foreach.hpp>

#include <CGAL/Iterator_range.h>
#include <CGAL/Timer.h>
#include <CGAL/Simple_cartesian.h>

typedef CGAL::Timer Timer;

int main()
{
  int N = 100000;
  std::vector<int> V(N), V2(N);
  
  Timer t;

  t.start();
  for(int k=0; k < N; k++){
    boost::iterator_range<std::vector<int>::iterator> bir(V.begin(), V.end());
    int j = 0;
    BOOST_FOREACH(int i, bir){
      V2[j++] = i;
    }
  }
  t.stop();
  std::cerr << "boost::iterator_range: " << t.time() << "sec.\n";

  t.reset();
  t.start();
  for(int k=0; k < N; k++){
       int j = 0;
    BOOST_FOREACH(int i, CGAL::make_range(V.begin(), V.end())){
      V2[j++] = i;
    }
  }
  t.stop();
  std::cerr << "CGAL::iterator_range: " << t.time() << "sec.\n";

  return 0;
}
