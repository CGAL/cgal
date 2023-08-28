#define BOOST_TEST_MODULE Skiplist
#include <boost/test/included/unit_test.hpp>

#include <iostream>
#include <string>
#include <boost/range/algorithm.hpp>
#include <boost/assign.hpp>

#include <CGAL/Skiplist.h>

typedef CGAL::Skiplist<int> skip;
typedef skip::skip_iterator skip_iterator;
typedef skip::all_iterator all_iterator;

using namespace boost::assign;

struct Fixture {
  Fixture() {
    skips += 1, 2, 3, 4, 5, 6, 7;
    all = skips;
    l.insert(l.all_begin(), skips.begin(), skips.end());
  }
  std::vector<int> all;
  std::vector<int> skips;
  skip l;
};

BOOST_AUTO_TEST_CASE( test_construction )
{
  std::vector<int> v;
  v += 1, 2, 3, 4, 5, 6, 7;
  skip l(v.begin(), v.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(l.all_begin(), l.all_end(),
                                v.begin(), v.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(l.skip_begin(), l.skip_end(),
                                v.begin(), v.end());
}

BOOST_FIXTURE_TEST_CASE( test_insert, Fixture )
{
  l.insert(l.all_end(), 8);
  l.insert(l.all_end(), 9);
  all += 8, 9;
  // 8 and 9 should be included in the all range
  BOOST_CHECK_EQUAL_COLLECTIONS(l.all_begin(), l.all_end(),
                                all.begin(), all.end());
  // but should not be part of the skipped set
  BOOST_CHECK_EQUAL_COLLECTIONS(l.skip_begin(), l.skip_end(),
                                skips.begin(), skips.end());

  // clear and try again
  l.clear();
  BOOST_CHECK_EQUAL(l.all_size(), std::size_t(0));
  BOOST_CHECK_EQUAL(l.skip_size(), std::size_t(0));
  l.insert(l.all_begin(), all.begin(), all.end());
  skips += 8, 9;
  BOOST_CHECK_EQUAL_COLLECTIONS(l.all_begin(), l.all_end(),
                                all.begin(), all.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(l.skip_begin(), l.skip_end(),
                                skips.begin(), skips.end());

  // the same goes for inserting at an arbitrary position
  l.insert(std::next(l.all_begin(), 3)
           , 20);
  all.insert(std::next(all.begin(), 3)
             , 20);
  BOOST_CHECK_EQUAL_COLLECTIONS(l.skip_begin(), l.skip_end(),
                                skips.begin(), skips.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(l.all_begin(), l.all_end(),
                                all.begin(), all.end());
}

BOOST_FIXTURE_TEST_CASE( test_single_skip, Fixture )
{
  for(skip::all_iterator it = l.all_begin(); it != l.all_end(); ++it) {
    BOOST_CHECK(!l.is_skipped(it));
  }

  // skip somewhere in between and at the end
  // skip 2 and 7
  l.skip(std::next(l.all_begin()));
  BOOST_CHECK(l.is_skipped(std::next(l.all_begin())));
  skips.erase(std::remove(skips.begin(), skips.end(), 2), skips.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(l.skip_begin(), l.skip_end(),
                                skips.begin(), skips.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(l.all_begin(), l.all_end(),
                                all.begin(), all.end());
  l.skip(std::prev(l.all_end()));
  BOOST_CHECK(l.is_skipped(std::prev(l.all_end())));
  skips.erase(std::remove(skips.begin(), skips.end(), 7), skips.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(l.skip_begin(), l.skip_end(),
                                skips.begin(), skips.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(l.all_begin(), l.all_end(),
                                all.begin(), all.end());
  // skip the beginning
  l.skip(l.all_begin());
  BOOST_CHECK(l.is_skipped(l.all_begin()));
  skips.erase(std::remove(skips.begin(), skips.end(), 1), skips.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(l.skip_begin(), l.skip_end(),
                                skips.begin(), skips.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(l.all_begin(), l.all_end(),
                                all.begin(), all.end());
}

BOOST_FIXTURE_TEST_CASE( test_range_skip, Fixture )
{
  // drop all from 2 up to 4
  l.skip(std::next(l.all_begin()), std::next(l.all_begin(), 3));
  skips.erase(std::remove(skips.begin(), skips.end(), 2), skips.end());
  skips.erase(std::remove(skips.begin(), skips.end(), 3), skips.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(l.skip_begin(), l.skip_end(),
                                skips.begin(), skips.end());
  // drop 6 and 7
  l.skip(std::prev(l.all_end(), 2), l.all_end());
  skips.erase(std::remove(skips.begin(), skips.end(), 6), skips.end());
  skips.erase(std::remove(skips.begin(), skips.end(), 7), skips.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(l.skip_begin(), l.skip_end(),
                                skips.begin(), skips.end());
}

BOOST_FIXTURE_TEST_CASE( skip_all_case, Fixture )
{
  l.skip(l.all_begin(), l.all_end());
  skips.clear();
  BOOST_CHECK_EQUAL(l.all_size(), all.size());
  BOOST_CHECK_EQUAL(l.skip_size(), std::size_t(0));
  BOOST_CHECK_EQUAL_COLLECTIONS(l.skip_begin(), l.skip_end(),
                                skips.begin(), skips.end());
}

BOOST_AUTO_TEST_CASE( test_continous_insert )
{

}

BOOST_FIXTURE_TEST_CASE( test_unskip, Fixture )
{
  // skip 2 and 3
  l.skip(std::next(l.all_begin()), std::next(l.all_begin(), 3));
  skips.erase(std::remove(skips.begin(), skips.end(), 2), skips.end());
  skips.erase(std::remove(skips.begin(), skips.end(), 3), skips.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(l.skip_begin(), l.skip_end(),
                                skips.begin(), skips.end());

  // unskip 2
  l.unskip(std::next(l.skip_begin()), std::next(l.all_begin()));
  skips.insert(std::next(skips.begin()), 2);
  BOOST_CHECK_EQUAL_COLLECTIONS(l.skip_begin(), l.skip_end(),
                                skips.begin(), skips.end());
  // unskip 3
  l.unskip(std::next(l.skip_begin(), 2), std::next(l.all_begin(), 2));
  skips.insert(std::next(skips.begin(), 2), 3);
  BOOST_CHECK_EQUAL_COLLECTIONS(l.skip_begin(), l.skip_end(),
                                skips.begin(), skips.end());

}

BOOST_FIXTURE_TEST_CASE( test_push, Fixture )
{
  using namespace boost::assign;

  l.push_back(8);
  l.push_back(9);
  all.push_back(8); all.push_back(9);
  skips.push_back(8); skips.push_back(9);
  BOOST_CHECK_EQUAL_COLLECTIONS(l.skip_begin(), l.skip_end(),
                                skips.begin(), skips.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(l.all_begin(), l.all_end(),
                                all.begin(), all.end());

  l.push_front(0); l.push_front(-1);
  all.insert(all.begin(), 0);
  all.insert(all.begin(), -1);
  skips.insert(skips.begin(), 0);
  skips.insert(skips.begin(), -1);

  BOOST_CHECK_EQUAL_COLLECTIONS(l.skip_begin(), l.skip_end(),
                                skips.begin(), skips.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(l.all_begin(), l.all_end(),
                                all.begin(), all.end());

  // check how it works with skipped elements
  l.skip(l.all_begin(), std::next(l.all_begin(), 4));
  skips.erase(skips.begin(), std::next(skips.begin(), 4));
  l.push_front(20);
  l.push_front(21);
  skips.insert(skips.begin(), 20);
  skips.insert(skips.begin(), 21);
  all.insert(all.begin(), 20);
  all.insert(all.begin(), 21);
  BOOST_CHECK_EQUAL_COLLECTIONS(l.skip_begin(), l.skip_end(),
                                skips.begin(), skips.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(l.all_begin(), l.all_end(),
                                all.begin(), all.end());
}

BOOST_FIXTURE_TEST_CASE( test_implicit_conversion, Fixture )
{
  all_iterator all;
  skip_iterator skip = l.skip_begin();
  all = skip;
  BOOST_CHECK(all == l.all_begin());
}

BOOST_FIXTURE_TEST_CASE( test_erase, Fixture )
{
  // erase 3
  l.erase(std::next(l.all_begin(), 2));
  skips.erase(std::remove(skips.begin(), skips.end(), 3), skips.end());
  all.erase(std::remove(all.begin(), all.end(), 3), all.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(l.skip_begin(), l.skip_end(),
                                skips.begin(), skips.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(l.all_begin(), l.all_end(),
                                all.begin(), all.end());

  // skip 2 first and then erase it
  l.skip(std::next(l.all_begin()));
  skips.erase(std::remove(skips.begin(), skips.end(), 2), skips.end());
  all.erase(std::remove(all.begin(), all.end(), 2), all.end());
  l.erase(std::next(l.all_begin()));
  BOOST_CHECK_EQUAL_COLLECTIONS(l.skip_begin(), l.skip_end(),
                                skips.begin(), skips.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(l.all_begin(), l.all_end(),
                                all.begin(), all.end());
}

BOOST_AUTO_TEST_CASE( test_swap )
{
  using std::swap;
  Fixture a, b;
  all_iterator it = std::prev(b.l.all_end());
  b.l.push_back(8); b.l.push_back(9); b.l.push_back(10);
  swap(a.l, b.l);

  // b should be equal to the default
  BOOST_CHECK_EQUAL_COLLECTIONS(b.l.skip_begin(), b.l.skip_end(),
                                b.all.begin(), b.all.end());

  a.all += 8, 9, 10;
  a.skips += 8, 9, 10;
  // a should have this shape
  BOOST_CHECK_EQUAL_COLLECTIONS(a.l.all_begin(), a.l.all_end(),
                                a.all.begin(), a.all.end());
  // this iterator should still be valid and now point into a
  BOOST_CHECK_EQUAL_COLLECTIONS(it, a.l.all_end(),
                                std::prev(a.all.end(), 4), a.all.end());
}

BOOST_AUTO_TEST_CASE( test_splice )
{
  Fixture a, b;
  a.all.insert(std::next(a.all.begin()), b.all.begin(), b.all.end());
  a.skips.insert(std::next(a.skips.begin()), b.skips.begin(), b.skips.end());

  a.l.splice(std::next(a.l.skip_begin()), b.l, b.l.skip_begin(), b.l.skip_end());
  BOOST_CHECK_EQUAL_COLLECTIONS(a.l.all_begin(), a.l.all_end(),
                                a.all.begin(), a.all.end());

  BOOST_CHECK_EQUAL_COLLECTIONS(a.l.skip_begin(), a.l.skip_end(),
                                a.skips.begin(), a.skips.end());
}

// trick cgal_create_cmake_script
// int main()
// {
// }
