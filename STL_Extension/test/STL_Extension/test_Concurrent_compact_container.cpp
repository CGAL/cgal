// test program for Concurrent_compact_container.

#include <iostream>

#ifndef CGAL_LINKED_WITH_TBB

int main()
{
  std::cout << 
    "NOTICE: this test needs CGAL_LINKED_WITH_TBB, and will not be tested."
    << std::endl;
  return 0;
}

#else

#include <CGAL/basic.h>
#include <cassert>
#include <cstddef>
#include <list>
#include <vector>
#include <CGAL/Compact_container.h>
#include <CGAL/Concurrent_compact_container.h>
#include <CGAL/Random.h>
#include <CGAL/Testsuite/use.h>

# include <tbb/parallel_for.h>
# include <tbb/atomic.h>

struct Node_1
: public CGAL::Compact_container_base
{
  bool operator==(const Node_1 &) const { return true; }
  bool operator!=(const Node_1 &) const { return false; }
  bool operator< (const Node_1 &) const { return false; }
};

class Node_2
{
  union {
    Node_2 * p;
    void * p_cc;
  };

public:
  
  int      rnd;

  Node_2()
  : p(NULL), rnd(CGAL::get_default_random().get_int(0, 100)) {}

  bool operator==(const Node_2 &n) const { return rnd == n.rnd; }
  bool operator!=(const Node_2 &n) const { return rnd != n.rnd; }
  bool operator< (const Node_2 &n) const { return rnd <  n.rnd; }

  void *   for_compact_container() const { return p_cc; }
  void * & for_compact_container()       { return p_cc; }
};

template < class Cont >
inline bool check_empty(const Cont &c)
{
  return c.empty() && c.size() == 0 && c.begin() == c.end();
}

// For parallel_for
template <typename Values_vec, typename Cont>
class Insert_in_CCC_functor
{
  typedef std::vector<typename Cont::iterator> Iterators_vec;

public:
  Insert_in_CCC_functor(
    const Values_vec &values, Cont &cont, Iterators_vec &iterators)
    : m_values(values), m_cont(cont), m_iterators(iterators) 
  {}
  
  Insert_in_CCC_functor(const Insert_in_CCC_functor &other)
    : m_values(other.m_values), m_cont(other.m_cont), 
      m_iterators(other.m_iterators)
  {}

  void operator() (const tbb::blocked_range<size_t>& r) const
  {
    for( size_t i = r.begin() ; i != r.end() ; ++i)
      m_iterators[i] = m_cont.insert(m_values[i]);
  }

private:
  const Values_vec    & m_values;
  Cont                & m_cont;
  Iterators_vec       & m_iterators;
};

// For parallel_for
template <typename Cont>
class Erase_in_CCC_functor
{
  typedef std::vector<typename Cont::iterator> Iterators_vec;

public:
  Erase_in_CCC_functor(
    Cont &cont, Iterators_vec &iterators)
    : m_cont(cont), m_iterators(iterators) 
  {}
  
  Erase_in_CCC_functor(const Erase_in_CCC_functor &other)
    : m_cont(other.m_cont), 
      m_iterators(other.m_iterators)
  {}

  void operator() (const tbb::blocked_range<size_t>& r) const
  {
    for( size_t i = r.begin() ; i != r.end() ; ++i)
      m_cont.erase(m_iterators[i]);
  }

private:
  Cont                & m_cont;
  Iterators_vec       & m_iterators;
};

// For parallel_for
template <typename Values_vec, typename Cont>
class Insert_and_erase_in_CCC_functor
{
  typedef std::vector<typename Cont::iterator>  Iterators_vec;
  typedef std::vector<tbb::atomic<bool> >       Free_elts_vec;

public:
  Insert_and_erase_in_CCC_functor(
    const Values_vec &values, Cont &cont, Iterators_vec &iterators,
    Free_elts_vec &free_elements, tbb::atomic<unsigned int> &num_erasures)
  : m_values(values), m_cont(cont), m_iterators(iterators),
    m_free_elements(free_elements), m_num_erasures(num_erasures)
  {}
  
  Insert_and_erase_in_CCC_functor(const Insert_and_erase_in_CCC_functor &other)
    : m_values(other.m_values), m_cont(other.m_cont), 
      m_iterators(other.m_iterators), m_free_elements(other.m_free_elements),
      m_num_erasures(other.m_num_erasures)
  {}

  void operator() (const tbb::blocked_range<size_t>& r) const
  {
    for( size_t i = r.begin() ; i != r.end() ; ++i)
    {
      m_iterators[i] = m_cont.insert(m_values[i]);
      // Random-pick an element to erase
      int index_to_erase = rand() % m_values.size();
      // If it exists
      if (m_free_elements[index_to_erase].compare_and_swap(true, false) == false)
      {
        m_cont.erase(m_iterators[index_to_erase]);
        ++m_num_erasures;
      }
    }
  }

private:
  const Values_vec          & m_values;
  Cont                      & m_cont;
  Iterators_vec             & m_iterators;
  Free_elts_vec             & m_free_elements;
  tbb::atomic<unsigned int> & m_num_erasures;
};

template < class Cont >
void test(const Cont &)
{
  // Testing if all types are provided.

  typename Cont::value_type              t0;
  typename Cont::reference               t1 = t0;      CGAL_USE(t1);
  typename Cont::const_reference         t2 = t0;      CGAL_USE(t2);
  typename Cont::pointer                 t3 = &t0;
  typename Cont::const_pointer           t4 = &t0;     CGAL_USE(t4);
  typename Cont::size_type               t5 = 0;       CGAL_USE(t5);
  typename Cont::difference_type         t6 = t3-t3;   CGAL_USE(t6);
  typename Cont::iterator                t7;           CGAL_USE(t7);
  typename Cont::const_iterator          t8;           CGAL_USE(t8);
  typename Cont::reverse_iterator        t9;           CGAL_USE(t9);
  typename Cont::const_reverse_iterator  t10;          CGAL_USE(t10);
  typename Cont::allocator_type          t15;

  std::cout << "Testing empty containers." << std::endl;

  Cont c0, c1;
  Cont c2(t15);
  Cont c3(c2);
  Cont c4;
  c4 = c2;

  typedef std::vector<typename Cont::value_type> Vect;
  Vect v0;
  const Cont c5(v0.begin(), v0.end());
  Cont c6(c5.begin(), c5.end());
  typename Cont::allocator_type Al;
  Cont c7(c0.begin(), c0.end(), Al);
  Cont c8;
  c8.insert(c0.rbegin(), c0.rend());

  // test conversion iterator-> const_iterator.
  typename Cont::const_iterator t16 = c5.begin();  CGAL_USE(t16);
  assert(t16 == c5.begin());

  assert(c0 == c1);
  assert(! (c0 < c1));

  assert(check_empty(c0));
  assert(check_empty(c1));
  assert(check_empty(c2));
  assert(check_empty(c3));
  assert(check_empty(c4));
  assert(check_empty(c5));
  assert(check_empty(c6));
  assert(check_empty(c7));
  assert(check_empty(c8));

  c1.swap(c0);

  assert(check_empty(c0));
  assert(check_empty(c1));

  c1.merge(c0);

  assert(check_empty(c0));
  assert(check_empty(c1));

  typename Cont::allocator_type  t20 = c0.get_allocator();

  std::cout << "Now filling some containers" << std::endl;

  Vect v1(10000);
  Cont c9(v1.begin(), v1.end());

  assert(c9.size() == v1.size());
  assert(c9.max_size() >= v1.size());
  assert(c9.capacity() >= c9.size());

  Cont c10 = c9;

  assert(c10 == c9);
  assert(c10.size() == v1.size());
  assert(c10.max_size() >= v1.size());
  assert(c10.capacity() >= c10.size());

  c9.clear();

  assert(check_empty(c9));
  assert(c9.capacity() >= c9.size());
  assert(c0 == c9);

  c9.merge(c10);
  c10.swap(c9);

  assert(check_empty(c9));
  assert(c9.capacity() >= c9.size());

  assert(c10.size() == v1.size());
  assert(c10.max_size() >= v1.size());
  assert(c10.capacity() >= c10.size());

  std::cout << "Testing insertion methods" << std::endl;

  c9.assign(c10.begin(), c10.end());

  assert(c9 == c10);

  c10.assign(c9.begin(), c9.end());

  assert(c9 == c10);

  c9.insert(c10.begin(), c10.end());

  assert(c9.size() == 2*v1.size());

  c9.clear();

  assert(c9 != c10);

  c9.insert(c10.begin(), c10.end());

  assert(c9.size() == v1.size());
  assert(c9 == c10);


  typename Cont::iterator it = c9.iterator_to(*c9.begin());
  assert(it == c9.begin());
  typename Cont::const_iterator cit = c9.iterator_to(const_cast<typename Cont::const_reference>(*c9.begin()));
  assert(cit == c9.begin());

  typename Cont::iterator s_it = Cont::s_iterator_to(*c9.begin());
  assert(s_it == c9.begin());
  typename Cont::const_iterator s_cit = Cont::s_iterator_to(const_cast<typename Cont::const_reference>(*c9.begin()));
  assert(s_cit == c9.begin());


  c10 = Cont();

  assert(check_empty(c10));

  for(typename Vect::const_iterator it = v1.begin(); it != v1.end(); ++it)
    c10.insert(*it);

  assert(c10.size() == v1.size());
  assert(c9 == c10);

  c9.erase(c9.begin());
  c9.erase(c9.begin());

  assert(c9.size() == v1.size() - 2);

  // test reserve
  /*Cont c11;
  c11.reserve(v1.size());
  for(typename Vect::const_iterator it = v1.begin(); it != v1.end(); ++it)
    c11.insert(*it);
  
  assert(c11.size() == v1.size());
  assert(c10 == c11);*/

  // owns() and owns_dereferencable().
  for(typename Cont::const_iterator it = c9.begin(), end = c9.end(); it != end; ++it) {
    assert(c9.owns(it));
    assert(c9.owns_dereferencable(it));
    assert(! c10.owns(it));
    assert(! c10.owns_dereferencable(it));
  }
  assert(c9.owns(c9.end()));
  assert(! c9.owns_dereferencable(c9.end()));


  c9.erase(c9.begin(), c9.end());

  assert(check_empty(c9));
  
  std::cout << "Testing parallel insertion" << std::endl;
  {
  Cont c11;
  Vect v11(1000000);
  std::vector<typename Cont::iterator> iterators(v11.size());
  tbb::parallel_for(
    tbb::blocked_range<size_t>( 0, v11.size() ),
    Insert_in_CCC_functor<Vect, Cont>(v11, c11, iterators)
  );
  assert(c11.size() == v11.size());
  
  std::cout << "Testing parallel erasure" << std::endl;
  tbb::parallel_for(
    tbb::blocked_range<size_t>( 0, v11.size() ),
    Erase_in_CCC_functor<Cont>(c11, iterators)
  );
  assert(c11.empty());
  }

  std::cout << "Testing parallel insertion AND erasure" << std::endl;
  {
  Cont c12;
  Vect v12(1000000);
  std::vector<tbb::atomic<bool> > free_elements(v12.size());
  for(typename std::vector<tbb::atomic<bool> >::iterator 
    it = free_elements.begin(), end = free_elements.end(); it != end; ++it) 
  {
    *it = true;
  }
    
  tbb::atomic<unsigned int> num_erasures; 
  num_erasures = 0;
  std::vector<typename Cont::iterator> iterators(v12.size());
  tbb::parallel_for(
    tbb::blocked_range<size_t>( 0, v12.size() ),
    Insert_and_erase_in_CCC_functor<Vect, Cont>(
      v12, c12, iterators, free_elements, num_erasures)
  );
  assert(c12.size() == v12.size() - num_erasures);
  }
}


int main()
{
  CGAL::Concurrent_compact_container<Node_1> C1;
  CGAL::Concurrent_compact_container<Node_2> C2;
  test(C1);
  test(C2);

  /*
  // Verbose merging test
  typedef CGAL::Concurrent_compact_container<Node_2> CCC;
  CCC cc1, cc2;
  tbb::parallel_for(
    tbb::blocked_range<size_t>( 0, 100, 1 ),
    [&] (const tbb::blocked_range<size_t>& r) // TODO: lambdas ok?
    {
      for( size_t i = r.begin() ; i != r.end() ; ++i)
      {
        Node_2 n;
        n.rnd = i;
        cc1.insert(n);
      }
    });

  for (int i = 10 ; i < 14 ; ++i)
  {
    Node_2 n;
    n.rnd = i;
    cc2.insert(n);
  }
  
  std::cout << "cc1 capacity: " << cc1.capacity() << std::endl;
  std::cout << "cc1 size: " << cc1.size() << std::endl;
  for(CCC::const_iterator it = cc1.begin(), end = cc1.end(); it != end; ++it) {
    std::cout << "cc1: " << it->rnd << " / " << std::endl;
  }
  std::cout << "cc2 capacity: " << cc2.capacity() << std::endl;
  std::cout << "cc2 size: " << cc2.size() << std::endl;
  for(CCC::const_iterator it = cc2.begin(), end = cc2.end(); it != end; ++it) {
    std::cout << "cc2: " << it->rnd << " / " << std::endl;
  }
  std::cout << "Merging...";
  //cc1.merge(cc2);
  cc2.merge(cc1);
  std::cout << " done." << std::endl;
  std::cout << "cc1 capacity: " << cc1.capacity() << std::endl;
  std::cout << "cc1 size: " << cc1.size() << std::endl;
  for(CCC::const_iterator it = cc1.begin(), end = cc1.end(); it != end; ++it) {
    std::cout << "cc1: " << it->rnd << " / " << std::endl;
  }
  std::cout << "cc2 capacity: " << cc2.capacity() << std::endl;
  std::cout << "cc2 size: " << cc2.size() << std::endl;
  for(CCC::const_iterator it = cc2.begin(), end = cc2.end(); it != end; ++it) {
    std::cout << "cc2: " << it->rnd << " / " << std::endl;
  }*/

  return 0;
}

#endif // CGAL_LINKED_WITH_TBB

// EOF //
