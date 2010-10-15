#include <CGAL/basic.h>
#include <CGAL/Double_map.h>
#include <utility>
#include <CGAL/Memory_sizer.h>
#include <CGAL/Timer.h>

#include <list>
#include <iostream>

#include <cstdlib>

typedef CGAL::Double_map<int, int> Map;
typedef Map::size_type size_type;

int main(int argc, char** argv)
{
  size_type number_of_elements = 100000;
  int number_of_loops = 10;

  if(argc > 1)
    number_of_elements = std::atoi(argv[1]);
  if(argc > 2)
    number_of_loops = std::atoi(argv[2]);

  Map f;
  CGAL::Timer time_insert;
  CGAL::Timer time_pop;
  CGAL::Timer time;
  CGAL::Memory_sizer memory;

  time.start();
  for(int loop = 0; loop < number_of_loops; ++loop)
  {
    time_pop.start();
    while(!f.empty())
      f.pop_front();
    time_pop.stop();
    time_insert.start();
    for(size_type i = 0; i < number_of_elements; ++i)
      f.insert(i, i);
    time_insert.stop();
  }
  time.stop();
  std::cerr << "Total time:           " << time.time()
	    << "\nTime for 'insert':    " << time_insert.time() 
	    << "\nTime for 'pop_front': " << time_pop.time() 
	    << "\nResident memory:      " << memory.resident_size() 
	    << "\nVirtual  memory:      " << memory.virtual_size()
	    << "\n";
  return 0;
}
