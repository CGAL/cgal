//#undef CGAL_LINKED_WITH_TBB

#include <CGAL/tags.h>
#include <CGAL/Compact_container.h>
#include <CGAL/Real_timer.h>

#ifdef CGAL_LINKED_WITH_TBB
# include <tbb/blocked_range.h>
# include <tbb/parallel_for.h>
# include <tbb/parallel_do.h>
# include <tbb/cache_aligned_allocator.h>
#endif

//#define TEST_CCSS

#include <iostream>
#include <vector>

struct Sequential_with_forward_access_tag : public CGAL::Sequential_tag {};
struct Sequential_with_random_access_tag  : public CGAL::Sequential_tag {};

struct Parallel_for_tag : public CGAL::Parallel_tag {};
struct Parallel_do_tag  : public CGAL::Parallel_tag {};

struct Truc
{
  Truc(int v = 0) : value(v), /*value2(v), */p(NULL) {}
  void *   for_compact_container() const { return p; }
  void * & for_compact_container()       { return p; }

  int value;
  int value2;
  void * p;
};


inline void compute_the_thing(int &thing)
{
  thing = static_cast<int>(10.*sqrt(0.3 + thing*thing/1.5)/2.6);
  //thing = 10*thing;
}

inline void compute_the_thing(Truc &thing)
{
  thing.value2=thing.value;
  compute_the_thing(thing.value2);
}

#ifdef CGAL_LINKED_WITH_TBB
// For parallel_for
template <typename Array_t>
class Change_array_functor
{
public:
  Change_array_functor(Array_t &v)
    : m_v(v) {}

  Change_array_functor(const Change_array_functor &caf)
    : m_v(caf.m_v) {}

  void operator() (const tbb::blocked_range<size_t>& r) const
  {
    for(size_t i = r.begin(); i != r.end(); i++ )
      if (m_v.is_used(i))
        compute_the_thing(m_v[i]);
  }

private:
  Array_t &m_v;
};

// For parallel_for
// Specialization for std::vector
template <>
class Change_array_functor<std::vector<Truc> >
{
public:
  Change_array_functor(std::vector<Truc> &v)
    : m_v(v) {}
  
  Change_array_functor(const Change_array_functor &caf)
    : m_v(caf.m_v) {}

  void operator() (const tbb::blocked_range<size_t>& r) const
  {
    for(size_t i = r.begin(); i != r.end(); i++ )
      compute_the_thing(m_v[i]);
  }

private:
  std::vector<Truc> &m_v;
};

// For parallel_do
template <typename Array_t>
class Change_array_functor_2
{
public:
  void operator() (typename Array_t::value_type& truc) const
  {
    compute_the_thing(truc);
  }
};

// Parallel_for
// For vector only
double change_array(std::vector<Truc> &v, Parallel_for_tag)
{
  std::cout << "[Parallel] Using a tbb::parallel_for loop...";
  CGAL::Real_timer t;
  t.start();
  tbb::parallel_for(
    tbb::blocked_range<std::size_t>(0, v.size()),
    Change_array_functor<std::vector<Truc> >(v)
  );
  t.stop();
  std::cout << " done in " << t.time() << " seconds." << std::endl;
  return t.time();
}

// Parallel_for
// For Compact_container only
template <typename Array_t>
double change_array(Array_t &v, Parallel_for_tag)
{
  std::cout << "[Parallel] Using a tbb::parallel_for loop...";
  CGAL::Real_timer t;
  t.start();
  tbb::parallel_for(
    tbb::blocked_range<std::size_t>(0, v.capacity()),
    Change_array_functor<Array_t>(v)
  );
  t.stop();
  std::cout << " done in " << t.time() << " seconds." << std::endl;
  return t.time();
}

// Parallel_do
// For vector and Compact_container
template <typename Array_t>
double change_array(Array_t &v, Parallel_do_tag)
{
  std::cout << "[Parallel] Using a tbb::parallel_do loop...";
  CGAL::Real_timer t;
  t.start();
  tbb::parallel_do(
    v.begin(), v.end(),
    Change_array_functor_2<Array_t>());
  t.stop();
  std::cout << " done in " << t.time() << " seconds." << std::endl;
  return t.time();
}
#endif

// Forward access
// For vector & Compact_container
template <typename Array_t>
double change_array(Array_t &v, Sequential_with_forward_access_tag)
{
  std::cout << "[Sequential] Using a sequential for loop (forward access)...";
  CGAL::Real_timer t;
  t.start();
  typename Array_t::iterator it = v.begin(), it_end = v.end();
  
  for ( ; it != it_end ; ++it)
    compute_the_thing(*it);

  t.stop();
  std::cout << " done in " << t.time() << " seconds." << std::endl;
  return t.time();
}

// Random access
// For vector only
double change_array(std::vector<Truc> &v, Sequential_with_random_access_tag)
{
  std::cout << "[Sequential] Using a sequential for loop (random access)...";
  CGAL::Real_timer t;
  t.start();
  std::vector<Truc>::iterator it = v.begin(), it_end = v.end();

  for (int i = 0 ; i < v.size() ; ++i)
    compute_the_thing(v[i]);

  t.stop();
  std::cout << " done in " << t.time() << " seconds." << std::endl;
  return t.time();
}

// Random access
// For Compact_container (use of "is_used")
template <typename Array_t>
double change_array(Array_t &v, Sequential_with_random_access_tag)
{
  std::cout << "[Sequential] Using a sequential for loop (random access)...";
  CGAL::Real_timer t;
  t.start();
  for (int i = 0 ; i < v.capacity() ; ++i)
  {
    if (v.is_used(i))
      compute_the_thing(v[i]);
  }
  t.stop();
  std::cout << " done in " << t.time() << " seconds." << std::endl;
  return t.time();
}


template <typename Array_t>
void benchmark(Array_t &v)
{
  std::cout << "* Sequential algorithm (random access) => ";
  double seq_time_random = change_array(v, Sequential_with_random_access_tag());
  std::cout << "* Parallel_for algorithm => ";
  double parallel_for_time = change_array(v, Parallel_for_tag());
  std::cout << "Speed-up parallel_for (operator[]) = " 
    << seq_time_random/parallel_for_time << std::endl << std::endl;
  
  std::cout << "* Sequential algorithm (forward access) => ";
  double seq_time_fw = change_array(v, Sequential_with_forward_access_tag());
  std::cout << "* Parallel_do algorithm => ";
  double parallel_do_time = change_array(v, Parallel_do_tag());
  std::cout << "Speed-up parallel_do (iterators) = " 
    << seq_time_fw/parallel_do_time << std::endl << std::endl;
}

int main()
{
#ifdef _DEBUG
  const int NUM_THINGS = 1000;
#else
  const int NUM_THINGS = 10000000;
#endif

  std::cout << std::endl << "*** VECTOR ***" << std::endl;
  CGAL::Real_timer t;
  t.start();
  std::vector<Truc> v;
  for (int i = 0 ; i < NUM_THINGS ; ++i)
    v.push_back(i);
  t.stop();
  std::cout << "vector created in " << t.time() << " seconds." << std::endl;
  benchmark(v);

  std::cout << std::endl << "*** CC<Addition> ***" << std::endl;
  t.reset();
  t.start();
  CGAL::Compact_container<Truc> v2;
  for (int i = 0 ; i < NUM_THINGS ; ++i)
    v2.emplace(i);
  t.stop();
  std::cout << "Compact_container created in " << t.time() << " seconds. Capacity = " << v2.capacity() << std::endl;
  benchmark(v2);

  std::cout << std::endl << "*** CC<Constant> ***" << std::endl;
  t.reset();
  t.start();
  CGAL::Compact_container<Truc, CGAL::Default, CGAL::Constant_size_policy<1024> > v3;
  for (int i = 0 ; i < NUM_THINGS ; ++i)
    v3.emplace(i);
  t.stop();
  std::cout << "Compact_container<Constant_size_policy> created in " << t.time() << " seconds. Capacity = " << v3.capacity() << std::endl;
  benchmark(v3);

  return 0;
}
