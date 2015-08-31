// test program for Compact_container.

#include <CGAL/basic.h>
#include <cassert>
#include <cstddef>
#include <list>
#include <vector>
#include <CGAL/Compact_container.h>
#include <CGAL/Random.h>
#include <CGAL/Testsuite/use.h>

#include <CGAL/tags.h>
#include <CGAL/use.h>
#include <CGAL/assertions.h>

class Node_1
{
  union {
    Node_1 * p;
    void * p_cc;
  };

public:

  Node_1() : p(NULL)
  {}

  void *   for_compact_container() const { return p_cc; }
  void * & for_compact_container()       { return p_cc; }
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
  CGAL_USE(t20);

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
  Cont c11;
  c11.reserve(v1.size());
  for(typename Vect::const_iterator it = v1.begin(); it != v1.end(); ++it)
    c11.insert(*it);

  assert(c11.size() == v1.size());
  assert(c10 == c11);

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
}

int main()
{
  typedef CGAL::Compact_container<Node_1> C1;

  C1 c1;
  if (!c1.empty())
  {
    std::cout<<"PB new container is not empty."<<std::endl;
    return EXIT_FAILURE;
  }

  for (int i = 0 ; i < 1000 ; ++i)
  {
    C1::iterator it=c1.emplace();
    if ( !c1.is_used(it) )
    {
      std::cout<<"PB new emplace element is not used."<<std::endl;
      return EXIT_FAILURE;
    }
  }

  std::size_t nb=0;
  for (C1::iterator it = c1.begin(), itend=c1.end(); it!=itend; ++it, ++nb)
  {
    c1.erase(it);
    if ( c1.is_used(it) )
    {
      std::cout<<"PB erase element is used."<<std::endl;
      return EXIT_FAILURE;
    }
    if ( c1.is_used(nb) )
    {
      std::cout<<"PB erase element is used (2)."<<std::endl;
      return EXIT_FAILURE;
    }
  }

  if (!c1.empty())
  {
    std::cout<<"PB container at the end is not empty."<<std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
