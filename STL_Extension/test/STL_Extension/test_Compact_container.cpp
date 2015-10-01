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

#include <boost/type_traits/is_base_of.hpp>

template <typename Has_timestamp_ = CGAL::Tag_true>
struct Node_1
: public CGAL::Compact_container_base
{
  bool operator==(const Node_1 &) const { return true; }
  bool operator!=(const Node_1 &) const { return false; }
  bool operator< (const Node_1 &) const { return false; }

  // Erase counter (cf. Compact_container)
  unsigned int erase_counter() const
  {
    return this->m_erase_counter;
  }
  void set_erase_counter(unsigned int c)
  {
    this->m_erase_counter = c;
  }
  void increment_erase_counter()
  {
    ++this->m_erase_counter;
  }

  /// For the determinism of Compact_container iterators
  ///@{
  typedef Has_timestamp_ Has_timestamp;

  std::size_t time_stamp() const {
    return time_stamp_;
  }
  void set_time_stamp(const std::size_t& ts) {
    time_stamp_ = ts;
  }
  ///@}
  int m_erase_counter;
  std::size_t time_stamp_;
};

class Node_2
{
  union {
    Node_2 * p;
    void * p_cc;
  };
  int      rnd;

public:

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

template < class Cont >
void test_index(const Cont &C)
{
  test(C);

  Cont c1;
  for (int i = 0 ; i < 1000000 ; ++i)
    c1.emplace();

  typename Cont::iterator it = c1.begin();
  for (int i=0; i < 1000000 ; ++i, ++it)
  {
    assert( c1[i]==*it ); // test the contents
    if ( i%1000==0 )
    {
      assert( Cont::s_iterator_to(c1[i])==it );
      assert( c1.iterator_to(c1[i])==it );
      assert( Cont::s_iterator_to(c1[i])==it );
    }
  }
}

struct Incomplete_struct;

int main()
{

  typedef Node_1<CGAL::Tag_true > T1;
  typedef CGAL::Compact_container<T1> C1; //    with timestamps

  typedef Node_1<CGAL::Tag_false> T2;
  typedef CGAL::Compact_container<T2> C2; // without timestamps

  typedef CGAL::Compact_container<T2,
                                  CGAL::Default,
                                  CGAL::Default,
                                  CGAL::Time_stamper<T2> > C4;
                                          //    with timestamps

  typedef Node_2 T3;
  typedef CGAL::Compact_container<T3> C3; // without timestamps

  C1 c1;  test(c1);
  C2 c2;  test(c2);
  C3 c3;  test(c3);
  C4 c4;  test(c4);

  // Check the time stamper policies
  if(! boost::is_base_of<CGAL::Time_stamper<T1>,
     C1::Time_stamper_impl>::value)
  {
    std::cerr << "Error timestamper of C1\n"; return 1;
  }
  if(! boost::is_base_of<CGAL::No_time_stamp<T2>,
     C2::Time_stamper_impl>::value)
  {
    std::cerr << "Error timestamper of C2\n"; return 1;
  }
  if(! boost::is_base_of<CGAL::No_time_stamp<T3>,
     C3::Time_stamper_impl>::value)
  {
    std::cerr << "Error timestamper of C3\n"; return 1;
  }
  if(! boost::is_base_of<CGAL::Time_stamper<T2>,
     C4::Time_stamper_impl>::value)
  {
    std::cerr << "Error timestamper of C4\n"; return 1;
  }

  // Check that Compact_container does not require a complete type.
  CGAL_static_assertion(sizeof(CGAL::Compact_container<Incomplete_struct>) > 0);

  // Test increment policy
  CGAL::Compact_container<Node_2, CGAL::Default, CGAL::Constant_size_policy<1024> > C5;
  CGAL::Compact_container<Node_2, CGAL::Default, CGAL::Addition_size_policy<14,16> > C6;

  test_index(C5);
  test_index(C6);
  return 0;
}
// EOF //
