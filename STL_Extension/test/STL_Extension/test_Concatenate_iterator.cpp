#include <iostream>
#include <list>
#include <vector>
#include <cassert>

#include <CGAL/Concatenate_iterator.h>
#include <CGAL/use.h>

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

  return 0;
}

