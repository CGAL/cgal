// ============================================================================
//
// Copyright (c) 2003 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: $
// release_date  : $CGAL_Date: $
//
// file          : test_Compact_container.C
// chapter       : $CGAL_Chapter: STL Extensions for CGAL $
// package       : $CGAL_Package: STL_Extension $
// source        : stl_extension.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//                 Lutz Kettner <kettner@mpi-sb.mpg.de>
//                 Sylvain Pion <Sylvain.Pion@mpi-sb.mpg.de>
//
// maintainer    : Michael Hoffmann <hoffmann@inf.ethz.ch>
// coordinator   : ETH
//
// Stl_Extensions: Compact container.
// ============================================================================


#include <CGAL/basic.h>
#include <cassert>
#include <cstddef>
#include <list>
#include <vector>
#include <CGAL/Compact_container.h>
#include <CGAL/Random.h>

template < typename T >
void use(const T&) {}

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
  int      rnd;

public:

  Node_2()
  : p(NULL), rnd(CGAL::default_random.get_int(0, 100)) {}

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
  typename Cont::reference               t1 = t0;      use(t1);
  typename Cont::const_reference         t2 = t0;      use(t2);
  typename Cont::pointer                 t3 = &t0;
  typename Cont::const_pointer           t4 = &t0;     use(t4);
  typename Cont::size_type               t5 = 0;       use(t5);
  typename Cont::difference_type         t6 = t3-t3;   use(t6);
  typename Cont::iterator                t7;           use(t7);
  typename Cont::const_iterator          t8;           use(t8);
  typename Cont::reverse_iterator        t9;           use(t9);
  typename Cont::const_reverse_iterator  t10;          use(t10);
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
  Cont c7(c0.begin(), c0.end(), typename Cont::allocator_type());
  Cont c8;
#ifndef __SUNPRO_CC
  c8.insert(c0.rbegin(), c0.rend());
#endif

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

  c10 = Cont();

  assert(check_empty(c10));

  for(typename Vect::const_iterator it = v1.begin(); it != v1.end(); ++it)
    c10.insert(*it);

  assert(c10.size() == v1.size());
  assert(c9 == c10);

  c9.erase(c9.begin());
  c9.erase(c9.begin());

  assert(c9.size() == v1.size() - 2);

  c9.erase(c9.begin(), c9.end());

  assert(check_empty(c9));
}


int main()
{
  CGAL::Compact_container<Node_1> C1;
  CGAL::Compact_container<Node_2> C2;
  test(C1);
  test(C2);
  return 0;
}
// EOF //
