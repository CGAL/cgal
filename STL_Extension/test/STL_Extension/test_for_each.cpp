#include <CGAL/for_each.h>

#include <vector>
#include <list>

int main()
{
  std::vector<int> vec { 0, 1, 2, 3, 4, 5 };
  std::list<int> list { 0, 1, 2, 3, 4, 5 };

  std::cerr << "Testing sequential random access" << std::endl;
  CGAL::for_each<CGAL::Sequential_tag>
    (vec, [](int& i) -> bool { i *= 2; return true; });

  for (const int& i : vec)
    std::cerr << i << " ";
  std::cerr << std::endl;

#ifdef CGAL_LINKED_WITH_TBB
  std::cerr << "Testing parallel random access" << std::endl;
  CGAL::for_each<CGAL::Parallel_tag>
    (vec, [](int& i) -> bool { i *= 2; return true; });

  for (const int& i : vec)
    std::cerr << i << " ";
  std::cerr << std::endl;
#endif

  std::cerr << "Testing sequential non-random access" << std::endl;
  CGAL::for_each<CGAL::Sequential_tag>
    (list, [](int& i) -> bool { i *= 2; return true; });

  for (const int& i : list)
    std::cerr << i << " ";
  std::cerr << std::endl;

#ifdef CGAL_LINKED_WITH_TBB
  std::cerr << "Testing parallel non-random access" << std::endl;
  CGAL::for_each<CGAL::Parallel_tag>
    (list, [](int& i) -> bool { i *= 2; return true; });

  for (const int& i : list)
    std::cerr << i << " ";
  std::cerr << std::endl;
#endif

  return 0;
}
