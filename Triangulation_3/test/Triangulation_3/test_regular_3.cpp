#include <CGAL/Regular_triangulation_3.h>

bool del = true;

#include <CGAL/_test_types.h>
#include <CGAL/_test_cls_regular_3.h>
#include <CGAL/_test_cls_parallel_triangulation_3.h>

#include <cassert>
#include <iostream>
#include <list>

typedef CGAL::Exact_predicates_inexact_constructions_kernel EPIC;

// Explicit instantiation of the whole class :
template class CGAL::Regular_triangulation_3<EPIC>;

template<typename K>
void test_kernel()
{
  using Cls = CGAL::Regular_triangulation_3<K>;

  _test_cls_regular_3(Cls());

#ifdef CGAL_LINKED_WITH_TBB
  typedef CGAL::Spatial_lock_grid_3<
    CGAL::Tag_priority_blocking>                      Lock_ds;
  typedef CGAL::Triangulation_data_structure_3<
    CGAL::Regular_triangulation_vertex_base_3<K>,
    CGAL::Regular_triangulation_cell_base_3<K>,
    CGAL::Parallel_tag >                              Tds_parallel;
  typedef CGAL::Regular_triangulation_3<
    K, Tds_parallel, Lock_ds>                         RT_parallel;

  // The following test won't do things in parallel since it doesn't provide a lock data structure
  _test_cls_regular_3(RT_parallel());

  // This test performs parallel operations
  _test_cls_parallel_triangulation_3(RT_parallel());
#endif
}

int main()
{
  test_kernel<EPIC>();
  test_kernel<EPEC>();

  std::cout << "Done!" << std::endl;
  return EXIT_SUCCESS;
}

