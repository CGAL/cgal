// Copyright (c) 2012 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:  $
// $Id:  $
//
// Author(s)     : Clement Jamin

#include <CGAL/Mesh_3/Profiling_tools.h>
#include <CGAL/Compact_container.h>
//#include <CGAL/Compact_container - spion.h>
#include <CGAL/Concurrent_compact_container.h>

#include <tbb/tbb.h>

#include <iostream>


const size_t NUM_ITERATIONS = 20000000;

using namespace CGAL;

struct T
{
  T(int ii) : i(ii), p(0) {}

  int i;
  int *p;
  
  void *   for_compact_container() const { return (void*)p; }
  void * & for_compact_container()       { return (void*&)p; }
};

//typedef Compact_container<T, CGAL_ALLOCATOR(T), CGAL::Parallel_tag> Conc_container;
typedef Concurrent_compact_container<T> Conc_container;

int main()
{
  std::cerr << "== Compact_container == Iterator size: " << sizeof(Compact_container<T>::iterator) << std::endl;
  
  Compact_container<T> cc;
  WallClockTimer t1;
  for( size_t i = 0 ; i != NUM_ITERATIONS ; ++i)
  {
    std::vector<Compact_container<T>::iterator> listOfIts;
    if (listOfIts.empty())
    {
      Compact_container<T>::iterator it = cc.emplace(rand());
      if (rand() % 10 < 7)
        listOfIts.push_back(it);
    }
    else
    {
      cc.erase(listOfIts.back());
      listOfIts.pop_back();
    }
  }
  std::cerr << t1.elapsed() << " seconds." << std::endl;

  std::cerr << "== Concurrent_compact_container == Iterator size: " << sizeof(Conc_container::iterator) << std::endl;
  for (int num_threads = 1 ; num_threads <= 8 ; ++num_threads)
  {
    tbb::task_scheduler_init init(num_threads);
    Conc_container ccc;
  
    WallClockTimer t2;
    tbb::parallel_for(tbb::blocked_range<size_t>(0, NUM_ITERATIONS),
      [&]( const tbb::blocked_range<size_t>& r ) { // CJTODO: lambdas ok?
        for( size_t i = r.begin() ; i != r.end() ; ++i)
        {
          std::vector<Conc_container::iterator> listOfIts;
          if (listOfIts.empty())
          {
            Conc_container::iterator it = ccc.emplace(rand());
            if (rand() % 10 < 7)
              listOfIts.push_back(it);
          }
          else
          {
            ccc.erase(listOfIts.back());
            listOfIts.pop_back();
          }
        }
    });
    std::cerr << "Num_threads=" << num_threads << ": " << t2.elapsed() << " seconds." << std::endl;
  }

  /*ccc.emplace(12);
  Concurrent_compact_container<T>::iterator it2 = ccc.emplace(1);
  ccc.erase(it2);
  ccc.emplace(13);
  ccc.emplace(13);

  Concurrent_compact_container<T>::const_iterator it = ccc.begin();
  Concurrent_compact_container<T>::const_iterator it_end = ccc.end();
  for( ; it != it_end ; ++it)
    std::cerr << it->i << " " << std::endl;
*/
  return 0;
}
