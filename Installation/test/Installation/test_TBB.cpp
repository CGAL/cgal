//#undef CGAL_LINKED_WITH_TBB

#include <CGAL/tags.h>

#ifdef CGAL_LINKED_WITH_TBB
# include <tbb/blocked_range.h>
# include <tbb/parallel_for.h>
#endif

#include <iostream>
#include <vector>

#ifdef CGAL_LINKED_WITH_TBB
class Change_array_functor
{
public:
  Change_array_functor(std::vector<int> &v)
    : m_v(v) {}

  void operator() (const tbb::blocked_range<size_t>& r) const
  {
    for(size_t i = r.begin(); i != r.end(); i++ )
      m_v[i] *= 10;
  }

private:
  std::vector<int> &m_v;
};

void change_array(std::vector<int> &v, CGAL::Parallel_tag)
{
  std::cout << "[Parallel] Using a tbb::parallel_for loop...";
  tbb::parallel_for(
    tbb::blocked_range<std::size_t>(0, v.size()),
    Change_array_functor(v)
  );
  std::cout << " done." << std::endl;
}
#endif

void change_array(std::vector<int> &v, CGAL::Sequential_tag)
{
  std::cout << "[Sequential] Using a sequential for loop...";
  std::vector<int>::iterator it = v.begin(), it_end = v.end();
  for ( ; it != it_end ; ++it)
    *it *= 10;
  std::cout << " done." << std::endl;
}

int main()
{
  std::vector<int> v;
  for (int i = 0 ; i < 10000 ; ++i)
    v.push_back(i);

  std::cout << "Trying to call the sequential algorithm => ";
  change_array(v, CGAL::Sequential_tag());
  std::cout << "Trying to call the parallel algorithm => ";
  change_array(v, CGAL::Parallel_tag());

  /*for (int i = 0 ; i < 100 ; ++i)
    std::cerr << v[i] << " ";*/

  return 0;
}
