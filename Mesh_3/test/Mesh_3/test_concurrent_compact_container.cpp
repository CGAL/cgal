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
#include "XML_exporter.h"
//#include <CGAL/Compact_container.h>
//#include <CGAL/Compact_container - spion.h>
#include <CGAL/Concurrent_compact_container.h>
#include <CGAL/Concurrent_compact_container - TBB based.h>

#include <tbb/tbb.h>
/*
// AN ATTEMPT TO REPRODUCE THE C1060 ERROR => COULDN'T DO IT
#include <CGAL/Triangulation_ds_vertex_base_3.h>
#include <CGAL/Triangulation_data_structure_3.h>

#include <boost/bimap/bimap.hpp>
#include <boost/bimap/multiset_of.hpp>

#include <libQGLViewer/libQ

#include <tbb/tbb.h>

#include <iostream>

namespace CGAL
{
template <typename Tr,
          typename CornerIndex = int,
          typename CurveSegmentIndex = int,
          class Vb = Triangulation_ds_vertex_base_3<>,
          class Cb = Triangulation_ds_cell_base_3<> >
class Mesh_complex_3_in_triangulation_3
{
public:
  typedef Triangulation_data_structure_3<Vb,Cb > Tds;
  typedef typename Vb::template Rebind_TDS<Tds>::Other  Vertex;
  typedef typename Cb::template Rebind_TDS<Tds>::Other  Cell;

  typedef Concurrent_compact_container<Vertex>     Vertex_range;
  typedef typename Vertex_range::iterator          Vertex_iterator;
  typedef Vertex_iterator                          Vertex_handle;


  typedef int               CurveSegmentIndex;
  typedef CurveSegmentIndex Curve_segment_index;

  typedef boost::bimaps::bimap< 
    boost::bimaps::multiset_of<Vertex_handle>,
    boost::bimaps::multiset_of<Vertex_handle>,
    boost::bimaps::set_of_relation<>,
    boost::bimaps::with_info<Curve_segment_index> >   Edge_map;

  typedef typename Edge_map::value_type               Internal_edge; 
  
 
  Mesh_complex_3_in_triangulation_3() 
    : edges_()
  {
    Vertex_handle v1, v2;
    Internal_edge e = make_internal_edge(v1, v2);
  }

private:
  Internal_edge make_internal_edge(const Vertex_handle& v1,
                                   const Vertex_handle& v2) const
  {
    if ( v1 < v2 ) { return Internal_edge(v1,v2); }
    else { return Internal_edge(v2,v1); }
  }

  Edge_map edges_;
};


} // namespace CGAL

int main()
{
  CGAL::Mesh_complex_3_in_triangulation_3<int> m;
  return 0;
}*/

#ifdef _DEBUG
  const size_t NUM_ITERATIONS = 2000000;
#else
  const size_t NUM_ITERATIONS = 200000000;
#endif

using namespace CGAL;

struct T
{
  T(int ii) : i(ii), p(0) {}

  int i;
  int *p;
  
  void *   for_compact_container() const { return (void*)p; }
  void * & for_compact_container()       { return (void*&)p; }
};


template <typename Container>
void testContainer(const char *containerTypeName)
{
  std::cerr << "== " << containerTypeName << " ==" << std::endl;
  std::cerr << "* Iterator size = " << sizeof(Container::iterator) << std::endl;

  std::vector<std::string> subelements;
  subelements.push_back("NumThreads");
  subelements.push_back("EraseRatio");
  subelements.push_back("Time");
  Simple_XML_exporter<int> xml("ContainerPerformance", "Perf", subelements);

  for (int num_threads = 1 ; num_threads <= 8 ; ++num_threads)
  {
    for (int erase_ratio = 0 ; erase_ratio <= 10 ; ++erase_ratio)
    {
      tbb::task_scheduler_init init(num_threads);
      Container ccc;
  
      WallClockTimer t2;
      tbb::parallel_for(tbb::blocked_range<size_t>(0, NUM_ITERATIONS),
        [&]( const tbb::blocked_range<size_t>& r ) { // CJTODO: lambdas ok?
          for( size_t i = r.begin() ; i != r.end() ; ++i)
          {
            Container::iterator it = ccc.emplace(rand());
            if (rand() % 10 < erase_ratio)
            {
              ccc.erase(it);
            }
          }
      });
      std::vector<int> value;
      value.push_back( num_threads );
      value.push_back( erase_ratio );
      value.push_back( static_cast<int>(t2.elapsed()*1000) );
      xml.add_element(value);
      std::cerr << "Num_threads=" << num_threads << " / Erase_ratio=" << (float)erase_ratio/10.f << " => " << t2.elapsed() << " seconds." << std::endl;
    }
  }

  xml.export_to_xml(std::string("C:\\INRIA\\Reports\\CGAL\\Concurrent_compact_container\\") + containerTypeName + ".xml");
}

typedef Concurrent_compact_container<T> Conc_container_CCbased;
typedef Concurrent_compact_container_TBB_based<T> Conc_container_TBBbased;

int main()
{
  // PERFORMANCE TEST
  testContainer<Conc_container_CCbased>("Concurrent_compact_container - CC based");
  testContainer<Conc_container_TBBbased>("Concurrent_compact_container - TBB based");
  /*
  std::cerr << "== Compact_container == Iterator size: " << sizeof(Compact_container<T>::iterator) << std::endl;
  
  Compact_container<T> cc;
  WallClockTimer t1;
  for( size_t i = 0 ; i != NUM_ITERATIONS ; ++i)
  {
    Compact_container<T>::iterator it = cc.emplace(rand());
    if (rand() % 10 < 7)
    {
      cc.erase(it);
    }
  }
  std::cerr << t1.elapsed() << " seconds." << std::endl;*/

  /*Conc_container ccc;
  ccc.emplace(12);
  Conc_container::iterator it2 = ccc.emplace(1);
  ccc.erase(it2);
  ccc.emplace(13);
  ccc.emplace(13);

  Conc_container::const_iterator it = ccc.begin();
  Conc_container::const_iterator it_end = ccc.end();
  for( ; it != it_end ; ++it)
    std::cerr << it->i << " " << std::endl;*/

  return 0;
}
