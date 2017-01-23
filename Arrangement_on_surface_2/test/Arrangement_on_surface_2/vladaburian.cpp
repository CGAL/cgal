#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel   Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>                  Traits_2;
typedef CGAL::Arrangement_2<Traits_2>                       Arrangement_2;
typedef Traits_2::Segment_2                                 Segment_2;

void try_arr_insert(std::vector<Segment_2> &segments)
{
  Arrangement_2 arr;
  insert(arr, segments.begin(), segments.end());
  std::cout << arr << std::endl;
  if (! CGAL::is_valid(arr)) {
    std::cerr << "ERROR :  failed!" << std::endl;
    return;
  }
}

int main()
{
  std::vector<Segment_2> segments;

  segments = {
    {{1,3},{4,0}},
    {{4,0},{5,0}},
    {{5,0},{7,2}},
    {{0,0},{9,0}},
    {{3,0},{6,0}},
    //{{7,1},{8,0}},
    {{2,0},{8,0}},
  };

  try_arr_insert(segments);

  segments = {
    {{1,3},{4,0}},
    {{4,0},{5,0}},
    //{{5,0},{7,2}},
    {{0,0},{9,0}},
    {{3,0},{6,0}},
    //{{7,1},{8,0}},
    {{2,0},{8,0}},
  };

  try_arr_insert(segments);

  segments = {
    {{1,3},{4,0}},
    {{4,0},{5,0}},
    //{{5,0},{7,2}},
    {{0,0},{9,0}},
    {{3,0},{6,0}},
    {{7,1},{8,0}},
    {{2,0},{8,0}},
  };

  try_arr_insert(segments);

  return 0;
}
