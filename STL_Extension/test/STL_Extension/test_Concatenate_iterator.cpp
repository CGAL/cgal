#include <iostream>
#include <list>
#include <vector>
#include <cassert>

#include <CGAL/Concatenate_iterator.h>
#include <CGAL/use.h>
#include <CGAL/Simple_cartesian.h>

template<class C1, class C2>
class Concatenate_container
{
public:
  typedef
  CGAL::Concatenate_iterator<typename C1::iterator,
                             typename C2::iterator>  Iterator;

  typedef typename C1::size_type  size_type;

  Concatenate_container(const C1& c1, const C2& c2)
    : c1(c1), c2(c2) {}

  Iterator begin() const {
    return Iterator(c1.end(), c2.begin(), c1.begin());
  }

  Iterator end() const {
    return Iterator(c1.end(), c2.begin(), c2.end(), 0);
  }

  size_type size() const { return c1.size() + c2.size(); }
  size_type lazy_size()
  {
    size_type k(0);
    for (Iterator it = begin(); it != end(); it++) { k++; }
    return k;
  }


private:
  mutable C1  c1;
  mutable C2  c2;
};



template<class C, class I>
void print_container(C c, bool is_concatenated = true)
{
  typedef    I       Iterator;

  std::cout << "Contents of ";
  std::cout << (is_concatenated ? "concatenated " : "input ");
  std::cout << "container in forward order:" << std::endl;
  for (Iterator it = c.begin(); it != c.end(); it++) {
    std::cout << " " << (*it);
  }
  std::cout << std::endl << std::endl;

  std::cout << "Contents of ";
  std::cout << (is_concatenated ? "concatenated " : "input ");
  std::cout << "container in reverse order:" << std::endl;
  if ( c.begin() != c.end() ) {
    Iterator it = c.end();
    --it;
    for (; it != c.begin(); it--) {
      std::cout << " " << (*it);
    }
    std::cout << " " << *it;
  }
  std::cout << std::endl << std::endl;
}


template<class C, class A, class B>
void copy(C& c, A a, B b)
{
  c.clear();

  for (typename A::iterator it = a.begin(); it != a.end(); it++) {
    c.push_back(*it);
  }

  for (typename B::iterator it = b.begin(); it != b.end(); it++) {
    c.push_back(*it);
  }
}


template<class C, class AB>
bool test_creation(C c, AB ab)
{
  typedef typename C::iterator   It_c;
  typedef typename AB::Iterator  It_ab;

  It_c  it_c;
  It_ab it_ab;
  for (it_c = c.begin(), it_ab = ab.begin(); it_ab != ab.end();
       it_c++, it_ab++) {
    assert( *it_c == *it_ab );
  }

  return true;
}

template<class A, class B>
int test(A a, B b)
{
  typedef Concatenate_container<A,B>        AB_container;
  typedef Concatenate_container<A,A>        AA_container;
  typedef Concatenate_container<B,A>        BA_container;

  typedef typename AB_container::Iterator   AB_iterator;
  typedef typename AA_container::Iterator   AA_iterator;
  typedef typename BA_container::Iterator   BA_iterator;

  CGAL_USE_TYPE(typename AB_iterator::reference);
  CGAL_USE_TYPE(typename AB_iterator::pointer);
  typedef typename AB_iterator::value_type         value_type;
  CGAL_USE_TYPE(typename AB_iterator::difference_type);
  CGAL_USE_TYPE(typename AB_iterator::iterator_category);

  std::vector<value_type>  c;


  AB_container ab(a, b);

  ::copy(c, a, b);
  assert( test_creation(c, ab) );
  assert( ab.size() == ab.lazy_size() && ab.size() == c.size() );

  AA_container aa(a, a);

  ::copy(c, a, a);
  assert( test_creation(c, aa) );
  assert( aa.size() == aa.lazy_size() && aa.size() == c.size() );

  BA_container ba(b, a);

  ::copy(c, b, a);
  assert( test_creation(c, ba) );
  assert( ba.size() == ba.lazy_size() && ba.size() == c.size() );


  std::cout << "==========================================="
            << std::endl << std::endl;
  print_container<A,typename A::iterator>(a, false);
  print_container<B,typename B::iterator>(b, false);
  std::cout << "AB container:" << std::endl;
  print_container<AB_container,AB_iterator>(ab);
  std::cout << "AA container:" << std::endl;
  print_container<AA_container,AA_iterator>(aa);
  std::cout << "BA container:" << std::endl;
  print_container<BA_container,BA_iterator>(ba);
  std::cout << "==========================================="
            << std::endl << std::endl << std::endl;

  return 1;
}


struct Integer
{
  Integer(int i) : i(i) {}
  int i;
};

bool operator==(const Integer& i, const Integer& j)
{
  return i.i == j.i;
}

std::ostream& operator<<(std::ostream& os, const Integer& i)
{
  return os << i.i;
}


int main()
{
  typedef int Data;

  std::vector<Data> v1;
  std::vector<Data> v2;
  std::vector<Data> v3;
  std::vector<Data> v4;
  std::list<Data> l1;
  std::list<Data> l2;


  // generic case
  {
    int n = 10;
    for (int i = 0; i < n; i++) {
      v1.push_back(i);
      v2.push_back(i+n);
      l1.push_back(i+2*n);
    }

    std::cout << "testing generic cases..." << std::endl;

    assert( test(v1, l1) );
    assert( test(l1, v1) );
    assert( test(v1, v2) );
    assert( test(v1, v1) );
  }

  // one is empty
  {
    l1.clear();
    v2.clear();

    std::cout << "testing cases where at least one container "
              << "is empty..." << std::endl;

    assert( test(v1, l1) );
    assert( test(l1, v1) );
    assert( test(v1, v2) );
    assert( test(v2, v2) );
  }

  // Cartesian_const_iterator
  {
    using Kernel = CGAL::Simple_cartesian<double>;
    using Point = Kernel::Point_3;
    using Cartesian_const_iterator = typename Kernel::Cartesian_const_iterator_3;
    using Construct_cartesian_const_iterator = typename Kernel::Construct_cartesian_const_iterator_3;
    using Concatenate_iterator =CGAL::Concatenate_iterator<Cartesian_const_iterator,Cartesian_const_iterator>;

    Point p0(0,1,2), p1(10,11,12);
    Construct_cartesian_const_iterator ccc;
    // return Iterator(c1.end(), c2.begin(), c1.begin());
    Concatenate_iterator ppb(ccc(p0,0), ccc(p1), ccc(p0));

    // return Iterator(c1.end(), c2.begin(), c2.end(), 0);
    Concatenate_iterator ppe(ccc(p0,0), ccc(p1), ccc(p1,0),0);
    assert(ppb != ppe);
    ++ppb;
    assert(ppb != ppe);
    ++ppb;
    assert(ppb != ppe);
     ++ppb;
    assert(ppb != ppe);
    for(; ppb != ppe; ++ppb){
      std::cout << *ppb << std::endl;
    }
  }
  {
    using Kernel = CGAL::Simple_cartesian<double>;
    using Point = Kernel::Point_3;
    using Cartesian_const_iterator = typename Kernel::Cartesian_const_iterator_3;
    using Construct_cartesian_const_iterator = typename Kernel::Construct_cartesian_const_iterator_3;
    using Concatenate_iterator =CGAL::Concatenate_iterator<Cartesian_const_iterator,Cartesian_const_iterator>;
    Point p[2] = {Point(0,1,2), Point(10,11,12)};

    Construct_cartesian_const_iterator ccc;
    Concatenate_iterator ppb(ccc(p[0],0), ccc(p[1]), ccc(p[0]));
    Concatenate_iterator ppe(ccc(p[0],0), ccc(p[1]), ccc(p[1],0),0);
    assert(ppb != ppe);
    ++ppb;
    assert(ppb != ppe);
     ++ppb;
     assert(ppb != ppe);
     ++ppb;
     assert(ppb != ppe);
  }

  // Test O(1) operator+ and operator- for random access iterators (vectors)
  {
    std::cout << "testing O(1) operator+ and operator-..." << std::endl;

    std::vector<int> v_a = {0, 1, 2, 3, 4};
    std::vector<int> v_b = {10, 11, 12, 13, 14};

    using VV_Iterator = CGAL::Concatenate_iterator<
      std::vector<int>::iterator, std::vector<int>::iterator>;

    VV_Iterator begin(v_a.end(), v_b.begin(), v_a.begin());
    VV_Iterator end(v_a.end(), v_b.begin(), v_b.end(), 0);

    // Test operator+ within first range
    VV_Iterator it = begin + 2;
    assert(*it == 2);

    // Test operator+ crossing from first to second range
    it = begin + 6;
    assert(*it == 11);

    // Test operator+ within second range
    it = begin + 8;
    assert(*it == 13);

    // Test operator- within second range
    it = end - 2;
    assert(*it == 13);

    // Test operator- crossing from second to first range
    it = end - 7;
    assert(*it == 3);

    // Test operator- within first range
    VV_Iterator mid = begin + 3;
    it = mid - 2;
    assert(*it == 1);

    // Test difference operator
    assert(end - begin == 10);
    assert((begin + 7) - begin == 7);
    assert(end - (end - 3) == 3);

    // Test operator+= and operator-=
    it = begin;
    it += 5;
    assert(*it == 10);
    it -= 3;
    assert(*it == 2);

    // Test negative offsets
    it = begin + 5;
    it = it + (-2);
    assert(*it == 3);
    it = it - (-3);
    assert(*it == 11);

    // Test negative differences
    assert(begin - end == -10);
    assert((begin + 3) - (begin + 7) == -4);
    assert((begin + 5) - (begin + 2) == 3);

    std::cout << "O(1) operator tests passed!" << std::endl;
  }


  //------------------------------------------------------------------
#if 0
  {
    vv_iterator vv_begin(v1.end(), v2.begin(), v1.begin());
    vv_iterator vv_end(v1.end(), v2.begin(), v2.end(),0);

    std::cout << "Vector-vector combined iterator:" << std::endl;
    for (vv_iterator vv_it = vv_begin; vv_it != vv_end; ++vv_it) {
      std::cout << " " << (*vv_it);
    }
    std::cout << std::endl << std::endl;
  }

  //------------------------------------------------------------------

  {
    vl_iterator vl_begin(v1.end(), l2.begin(), v1.begin());
    vl_iterator vl_end(v1.end(), l2.begin(), l2.end(),0);


    std::cout << "Vector-list combined iterator:" << std::endl;
    for (vl_iterator vl_it = vl_begin; vl_it != vl_end; ++vl_it) {
      std::cout << " " << (*vl_it);
    }
    std::cout << std::endl << std::endl;
  }
  //------------------------------------------------------------------
  {
    vv_iterator vv_begin(v3.end(), v4.begin(), v3.begin());
    vv_iterator vv_end(v3.end(), v4.begin(), v4.end(),0);

    std::cout << "Vector-vector combined iterator:" << std::endl;
    for (vv_iterator vv_it = vv_begin; vv_it != vv_end; ++vv_it) {
      std::cout << " " << (*vv_it);
    }
    std::cout << std::endl << std::endl;
  }

  //------------------------------------------------------------------
  {
    vv_iterator vv_begin(v3.end(), v2.begin(), v3.begin());
    vv_iterator vv_end(v3.end(), v2.begin(), v2.end(),0);

    std::cout << "Vector-vector combined iterator:" << std::endl;
    for (vv_iterator vv_it = vv_begin; vv_it != vv_end; ++vv_it) {
      std::cout << " " << (*vv_it);
    }
    std::cout << std::endl << std::endl;
  }
#endif

  // Test mixed iterator categories (Random Access + Bidirectional) regression
  {
      std::cout << "testing mixed iterator categories (Random Access + Bidirectional) regression..." << std::endl;
      std::vector<int> v = {1, 2, 3};
      std::list<int> l = {4, 5, 6};
      using Iterator = CGAL::Concatenate_iterator<std::vector<int>::iterator, std::list<int>::iterator>;
      Iterator begin(v.end(), l.begin(), v.begin());

      Iterator it = begin + 4;
      assert(*it == 5);
      std::cout << "Mixed iterator category test passed!" << std::endl;
  }

  return 0;
}

