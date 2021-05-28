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
// file          : test_stl_extension.C
// chapter       : $CGAL_Chapter: STL Extensions for CGAL $
// package       : $CGAL_Package: STL_Extension $
// source        : stl_extension.fw
// revision      : $Id$
// revision_date : $Date$
// author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//                 Lutz Kettner <kettner@mpi-sb.mpg.de>
//                 Sylvain Pion
//
// maintainer    : Michael Hoffmann <hoffmann@inf.ethz.ch>
// coordinator   : ETH
//
// Stl_Extensions: Iterator Adaptor.
// ============================================================================


#include <cstddef>
#include <list>
#include <vector>
#include <cassert>
#include <CGAL/use.h>
#include <CGAL/utility.h>
#include <CGAL/iterator.h>
#include <CGAL/algorithm.h>
#include <CGAL/In_place_list.h>
#include <CGAL/Circulator_identity.h>
#include <CGAL/Iterator_project.h>
#include <CGAL/Iterator_transform.h>
#include <CGAL/function_objects.h>
#include <CGAL/Circulator_project.h>
#include <CGAL/Circulator_on_node.h>
#include <CGAL/Circulator/Circulator_adapters.h>
#include <CGAL/function_objects.h>
#include <CGAL/tuple.h>
#include <CGAL/assertions.h>
#include <boost/type_traits/is_same.hpp>
#include <boost/typeof/typeof.hpp>

#include <CGAL/internal/disable_deprecation_warnings_and_errors.h>
#include <CGAL/result_of.h>
using namespace CGAL;

struct item : public In_place_list_base<item> {
  int key;
  item() {}
  item( int i) : key(i) {}
  item( const item& i)
  : In_place_list_base<item>(i), key(i.key) {}
  bool operator== (const item& i) const { return key == i.key;}
  bool operator== (int i) const         { return key == i;}
  bool operator!= (int i) const         { return key != i;}
  bool operator<  (const item& i) const { return key < i.key;}
};

struct item_adaptor : public item
{
  item_adaptor() {}
  item_adaptor( int i) : item(i) {}
  item_adaptor( const item_adaptor& i) : item(i) {}

  int const& get_key() const { return this->key ; }
  int&       get_key()       { return this->key ; }
} ;

int test_value_type( item*)           { return 1;}
int test_value_type( item**)          { return 1;}
int test_value_type( item_adaptor*)   { return 1;}

// Global data structures.
std::list<int>   L;
std::vector<int> V;
// int*          C_array;  // Gnu gcov produces a core here!

// Workaround
int* new_C_array() {
    int* p = new int[5];
    for ( int i = 1; i <= 5; i++) {
        p[i-1] = i;
    }
    return p;
}

// Build a simple 'n'-element circular structure using struct's.
struct Node {
    int   key;
    Node* next;
    Node* prev;
    Node() : key(0), next(nullptr), prev(nullptr) { next = prev = this; }
    Node(int n) : key(n), next(nullptr), prev(nullptr) { next = prev = this; }
    Node(Node* nx_, Node* pv_, int n) : key(n), next(nx_), prev(pv_) {}
};
Node* new_node( Node* nx_, Node* pv_, int n) {
    return new Node( nx_, pv_, n);
}
void append_node( Node* p, int n) {
    Node* q = new_node( p, p->prev, n);
    p->prev->next = q;
    p->prev = q;
}
Node* generate_nodes( int n) {
    assert( n > 0);
    Node* p = new Node(1);
    for ( int i = 2; i <= n; i++)
        append_node( p, i);
    return p;
}
void delete_nodes( Node* p) {
    Node* end = p;
    Node* q   = p;
    p = p->next;
    while ( p != end) {
        delete q;
        q = p;
        p = p->next;
    }
    delete q;
}

typedef CGAL::Forward_circulator_over_struct< Node>
    Struct_circulator;
typedef CGAL::Forward_const_circulator_over_struct< Node>
    Struct_const_circulator;
typedef CGAL::Bidirectional_circulator_over_struct< Node>
    Struct_bi_circulator;
typedef CGAL::Bidirectional_const_circulator_over_struct< Node>
    Struct_bi_const_circulator;

// Build a simple 'n'-element circular structure using struct's.
class CNode {
    CNode* next_;
    CNode* prev_;
  public:
    int   key;
    CNode*       next()       { return next_;}
    const CNode* next() const { return next_;}
    CNode*       prev()       { return prev_;}
    const CNode* prev() const { return prev_;}
    CNode() : next_(nullptr), prev_(nullptr), key(0) { next_ = prev_ = this; }
    CNode( int n) : next_(nullptr), prev_(nullptr), key(n) { next_ = prev_ = this; }
    CNode( CNode* nx_, CNode* pv_, int n)
        : next_(nx_), prev_( pv_), key(n) {}
    friend CNode* new_cnode( CNode* nx_, CNode* pv_, int n);
    friend void append_cnode( CNode* p, int n);
    friend void delete_cnodes( CNode* p);
};
CNode* new_cnode( CNode* nx_, CNode* pv_, int n) {
    return new CNode( nx_, pv_, n);
}
void append_cnode( CNode* p, int n) {
    CNode* q = new_cnode( p, p->prev_, n);
    p->prev_->next_ = q;
    p->prev_ = q;
}
CNode* generate_cnodes( int n) {
    assert( n > 0);
    CNode* p = new CNode(1);
    for ( int i = 2; i <= n; i++)
        append_cnode( p, i);
    return p;
}
void delete_cnodes( CNode* p) {
    CNode* end = p;
    CNode* q   = p;
    p = p->next_;
    while ( p != end) {
        delete q;
        q = p;
        p = p->next_;
    }
    delete q;
}

typedef CGAL::Forward_circulator_over_class< CNode>
    Class_circulator;
typedef CGAL::Forward_const_circulator_over_class< CNode>
    Class_const_circulator;
typedef CGAL::Bidirectional_circulator_over_class< CNode>
    Class_bi_circulator;
typedef CGAL::Bidirectional_const_circulator_over_class< CNode>
    Class_bi_const_circulator;


struct Pair
{
  Pair ( int f_, int s_ ) : f(f_), s(s_) {}
  int f,s;
} ;

struct Point
{
  Point( int x_, int y_ ) : x(x_), y(y_) {}
  int x,y ;

  friend bool operator==( Point const& a, Point const& b ) { return a.x == b.x && a.y == b.y ; }
  friend bool operator!=( Point const& a, Point const& b ) { return !(a==b); }
} ;

int test_value_type( Point*)  { return 1;}

struct PairToPoint
{
  typedef Pair  argument_type ;
  typedef Point result_type ;

  Point operator() ( Pair const& p ) const { return Point(p.f,p.s); }
} ;

void init_global_data() {
    // C_array = new int[5];
    for ( int i = 1; i <= 5; i++) {
        L.push_back(i);
        V.push_back(i);
        // C_array[i-1] = i;
    }
}

void clean_global_data() {
    //delete[] C_array;
}

// Test value type and distance type.
int test_value_type( int*)                 { return 1;}
int test_value_type( Node*)                { return 1;}
int test_value_type( CNode*)               { return 1;}
int test_value_type( char*)                { return 2;}
int test_value_type( double*)              { return 3;}

int test_distance_type( std::ptrdiff_t*)   { return 1;}
int test_distance_type( char*)             { return 2;}
int test_distance_type( double*)           { return 3;}
void test_Circulator_identity() {
  {
    typedef In_place_list<item,false> List;
    List l;
    l.push_back( *new item(1));
    l.push_back( *new item(2));
    l.push_back( *new item(3));
    l.push_back( *new item(4));
    l.push_back( *new item(5));
    assert( l.size() == 5);
    typedef List::iterator IterBase;
    typedef Bidirectional_circulator_from_iterator<IterBase,item,
      std::size_t,std::ptrdiff_t> CircBase;
    typedef Circulator_identity<CircBase,item&,item*> Circulator;
    Circulator begin(CircBase( l.begin(),l.end()));
    Assert_bidirectional_category(begin);
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_is_at_least_forward_category(begin);
        CGAL::Assert_is_at_least_forward_category(begin);
        typedef std::iterator_traits< Circulator >::value_type      VT;
        typedef std::iterator_traits< Circulator >::difference_type DT;
        assert(1==test_value_type(static_cast< VT* >(nullptr)));
        assert(1==test_distance_type(static_cast< DT* >(nullptr)));

        // Default constructor.
        Circulator z = Circulator();
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Circulator i = begin;

        // Check general support for circulators and iterators.
        assert( CGAL::is_empty_range( z, z));
        assert( ! CGAL::is_empty_range( i, begin));

        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                assert( k == (*i).key);
                su += (*i).key;
                ++k;
                Circulator j = ++i;
                assert( i ==  j);
                if ( i != begin) {
                    assert( (*i).key == (*j).key);
                }
            } while (i != begin);  // Inequality and equality checked.
        }
        assert( i == begin);  // Equality checked.
        assert( su == 15);

        // Assignment.
        i = begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                assert( k == (*i).key);
                su += (*i).key;
                ++k;
                Circulator j = i++;
                assert(  i !=  j);
                if ( i != begin) {
                    assert( (*i).key == (*j).key + 1);
                }
            } while (i != begin);
        }
        assert( i == begin);
        assert( su == 15);
    }
    { // Open own scope to hide local variables.
        // Change three elements and check post-/pre-increment.
        Circulator i = begin;
        (*i++).key = 4;
        assert( 4 == (*begin).key);
        assert( 2 == (*i).key);
        (*i++).key = 3;
        assert( 3 == (*i).key);
        (*++i).key = 7;
        assert( 7 == (*i).key);

        // Check the setting and reset these elements
        // to their original values.
        i = begin;
        assert( 4 == (*i).key);
        (*i).key = 1;
        i++;
        assert( 3 == (*i).key);
        (*i++).key = 2;
        assert( 3 == (*i).key);
        i++;
        assert( 7 == (*i).key);
        (*i).key = 4;

        // Check the resetting.
        i = begin;
        int k = 1;
        do {
            assert( k == (*i).key);
            ++i;
            ++k;
        } while (i != begin);
    }
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_is_at_least_forward_category(begin);
        CGAL::Assert_is_at_least_forward_category(begin);
        typedef std::iterator_traits< Circulator >::value_type      VT;
        typedef std::iterator_traits< Circulator >::difference_type DT;
        assert(1==test_value_type(static_cast< VT* >(nullptr)));
        assert(1==test_distance_type(static_cast< DT* >(nullptr)));

        // Default constructor.
        Circulator z = Circulator();
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Circulator i = begin;

        // Check general support for circulators and iterators.
        assert( CGAL::is_empty_range( z, z));
        assert( ! CGAL::is_empty_range( i, begin));

        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                assert( k == (*i).key);
                su += (*i).key;
                ++k;
                Circulator j = ++i;
                assert( i ==  j);
                if ( i != begin) {
                    assert( (*i).key == (*j).key);
                }
            } while (i != begin);  // Inequality and equality checked.
        }
        assert( i == begin);  // Equality checked.
        assert( su == 15);

        // Assignment.
        i = begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                assert( k == (*i).key);
                su += (*i).key;
                ++k;
                Circulator j = i++;
                assert(  i !=  j);
                if ( i != begin) {
                    assert( (*i).key == (*j).key + 1);
                }
            } while (i != begin);
        }
        assert( i == begin);
        assert( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_bidirectional_category(begin);
        CGAL::Assert_is_at_least_bidirectional_category(begin);
        // Loop backwards and pre-decrement.
        Circulator i = begin;
        int su = 0;
        int k  = 5;
        do {
            Circulator j = --i;
            assert(  i ==  j);
            assert( (*i).key == (*j).key);
            assert( k == (*i).key);
            su += (*i).key;
            --k;
        } while (i != begin);
        assert( i == begin);
        assert( su == 15);

        // Assignment.
        i = begin;
        su = 0;
        k  = 5;
        // Loop with post-decrement.
        do {
            Circulator j = i--;
            assert(  i !=  j);
            if ( j != begin) {
                assert( (*i).key == (*j).key - 1);
            }
            assert( k == (*i).key);
            su += (*i).key;
            --k;
        } while (i != begin);
        assert( i == begin);
        assert( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_circulator( begin);
        CGAL::Assert_circulator( begin);

        // Check the local type parameters.
        Circulator::value_type      k1;
        k1.key = 1;
        Circulator::reference       k2 = k1;
        (void)k2;
        assert( k2.key == 1);
        Circulator::pointer         p1 = &k1;
        (void)p1;
        assert( (*p1).key == 1);
        k1.key = 3;
        assert( k1.key == 3);
        assert( k2.key == 3);
        assert( (*p1).key == 3);
        k1.key = 6;
        assert( k1.key == 6);
        assert( k2.key == 6);
        assert( (*p1).key == 6);
        Circulator::size_type s = 5;
        assert( s == 5);
        Circulator::difference_type d = -5;
        assert( d == -5);

        // Check tests for empty data structures.
        Circulator z = Circulator();
        assert(   z == nullptr);
        assert( ! (z != nullptr));
        Circulator i = begin;
        assert( ! (i == nullptr));
        assert(   i != nullptr);
        assert( i == begin);
        assert( i == begin);
        // Do I reach myself.
        ++i;
        Circulator j = i;
        int k = 0;
        do {
            assert( k < 5);
            ++k;
            ++i;
        } while( i != j);
        assert( k == 5);
    }
    { // Open own scope to hide local variables.
        // Do I reach myself backwards.
        Circulator i = begin;
        ++i;
        Circulator j = i;
        int k = 0;
        do {
            assert( k < 5);
            ++k;
            --i;
        } while( i != j);
        assert( k == 5);
    }

    List l2 = l;
    const List& l1 = l2;
    typedef List::const_iterator ConstIterBase;
    typedef Bidirectional_const_circulator_from_iterator<
      ConstIterBase, item, std::size_t, std::ptrdiff_t>
    ConstCircBase;
typedef Circulator_identity<ConstCircBase,const item&,
  const item*> C_Circulator;
C_Circulator c_begin(ConstCircBase(l1.begin(),l1.end()));
Assert_bidirectional_category(c_begin);
{ // Open own scope to hide local variables.
    // Check generally correct parameter properties.
    CGAL::Assert_circulator_or_iterator(c_begin);
    CGAL::Assert_circulator_or_iterator(c_begin);
    CGAL::Assert_is_at_least_forward_category(c_begin);
    CGAL::Assert_is_at_least_forward_category(c_begin);
    typedef std::iterator_traits< C_Circulator >::value_type      VT;
    typedef std::iterator_traits< C_Circulator >::difference_type DT;
    assert(1==test_value_type(static_cast< VT* >(nullptr)));
    assert(1==test_distance_type(static_cast< DT* >(nullptr)));

    // Default constructor.
    C_Circulator z = C_Circulator();
    CGAL::Assert_circulator_or_iterator(z);
    // Copy constructor.
    C_Circulator i = c_begin;

    // Check general support for circulators and iterators.
    assert( CGAL::is_empty_range( z, z));
    assert( ! CGAL::is_empty_range( i, c_begin));

    int su = 0;
    int k  = 1;
    // Check general loop, pre-increment, dereference.
    if (! CGAL::is_empty_range( i, c_begin)) {   // superfluous
        do {
            assert( k == (*i).key);
            su += (*i).key;
            ++k;
            C_Circulator j = ++i;
            assert( i ==  j);
            if ( i != c_begin) {
                assert( (*i).key == (*j).key);
            }
        } while (i != c_begin);  // Inequality and equality checked.
    }
    assert( i == c_begin);  // Equality checked.
    assert( su == 15);

    // Assignment.
    i = c_begin;
    su = 0;
    k  = 1;
    // Loop with post increment.
    if (! CGAL::is_empty_range( i, c_begin)) {   // superfluous
        do {
            assert( k == (*i).key);
            su += (*i).key;
            ++k;
            C_Circulator j = i++;
            assert(  i !=  j);
            if ( i != c_begin) {
                assert( (*i).key == (*j).key + 1);
            }
        } while (i != c_begin);
    }
    assert( i == c_begin);
    assert( su == 15);
}
{ // Open own scope to hide local variables.
    CGAL::Assert_is_at_least_bidirectional_category(c_begin);
    CGAL::Assert_is_at_least_bidirectional_category(c_begin);
    // Loop backwards and pre-decrement.
    C_Circulator i = c_begin;
    int su = 0;
    int k  = 5;
    do {
        C_Circulator j = --i;
        assert(  i ==  j);
        assert( (*i).key == (*j).key);
        assert( k == (*i).key);
        su += (*i).key;
        --k;
    } while (i != c_begin);
    assert( i == c_begin);
    assert( su == 15);

    // Assignment.
    i = c_begin;
    su = 0;
    k  = 5;
    // Loop with post-decrement.
    do {
        C_Circulator j = i--;
        assert(  i !=  j);
        if ( j != c_begin) {
            assert( (*i).key == (*j).key - 1);
        }
        assert( k == (*i).key);
        su += (*i).key;
        --k;
    } while (i != c_begin);
    assert( i == c_begin);
    assert( su == 15);
}
{ // Open own scope to hide local variables.
    CGAL::Assert_circulator( c_begin);
    CGAL::Assert_circulator( c_begin);

    // Check the local type parameters.
    C_Circulator::value_type      k1;
    k1.key = 1;
    C_Circulator::reference       k2 = k1;
    (void)k2;
    assert( k2.key == 1);
    C_Circulator::pointer         p1 = &k1;
    (void)p1;
    assert( (*p1).key == 1);
    k1.key = 3;
    assert( k1.key == 3);
    assert( k2.key == 3);
    assert( (*p1).key == 3);
    k1.key = 6;
    assert( k1.key == 6);
    assert( k2.key == 6);
    assert( (*p1).key == 6);
    C_Circulator::size_type s = 5;
    assert( s == 5);
    C_Circulator::difference_type d = -5;
    assert( d == -5);

    // Check tests for empty data structures.
    C_Circulator z = C_Circulator();
    assert(   z == nullptr);
    assert( ! (z != nullptr));
    C_Circulator i = c_begin;
    assert( ! (i == nullptr));
    assert(   i != nullptr);
    assert( i == c_begin);
    assert( i == c_begin);
    // Do I reach myself.
    ++i;
    C_Circulator j = i;
    int k = 0;
    do {
        assert( k < 5);
        ++k;
        ++i;
    } while( i != j);
    assert( k == 5);
}
{ // Open own scope to hide local variables.
    // Do I reach myself backwards.
    C_Circulator i = c_begin;
    ++i;
    C_Circulator j = i;
    int k = 0;
    do {
        assert( k < 5);
        ++k;
        --i;
    } while( i != j);
    assert( k == 5);
}
l.destroy();
l2.destroy();
  }
  {
    typedef std::vector<int> Vector;
    Vector v;
    v.push_back(1);
    v.push_back(2);
    v.push_back(3);
    v.push_back(4);
    v.push_back(5);
    assert( v.size() == 5);
    typedef Vector::iterator IterBase;
    typedef Random_access_circulator_from_iterator<IterBase,int,
      std::size_t,std::ptrdiff_t> CircBase;
    typedef Circulator_identity<CircBase,int&,int*> Circulator;
    Circulator begin(CircBase( v.begin(),v.end()));
    Assert_random_access_category(begin);
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_is_at_least_forward_category(begin);
        CGAL::Assert_is_at_least_forward_category(begin);
        typedef std::iterator_traits< Circulator >::value_type      VT;
        typedef std::iterator_traits< Circulator >::difference_type DT;
        assert(1==test_value_type(static_cast< VT* >(nullptr)));
        assert(1==test_distance_type(static_cast< DT* >(nullptr)));

        // Default constructor.
        Circulator z;
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Circulator i = begin;

        // Check general support for circulators and iterators.
        assert( CGAL::is_empty_range( z, z));
        assert( ! CGAL::is_empty_range( i, begin));

        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                assert( k == (*i));
                su += (*i);
                ++k;
                Circulator j = ++i;
                assert( i ==  j);
                if ( i != begin) {
                    assert( (*i) == (*j));
                }
            } while (i != begin);  // Inequality and equality checked.
        }
        assert( i == begin);  // Equality checked.
        assert( su == 15);

        // Assignment.
        i = begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                assert( k == (*i));
                su += (*i);
                ++k;
                Circulator j = i++;
                assert(  i !=  j);
                if ( i != begin) {
                    assert( (*i) == (*j) + 1);
                }
            } while (i != begin);
        }
        assert( i == begin);
        assert( su == 15);
    }
    { // Open own scope to hide local variables.
        // Change three elements and check post-/pre-increment.
        Circulator i = begin;
        (*i++) = 4;
        assert( 4 == (*begin));
        assert( 2 == (*i));
        (*i++) = 3;
        assert( 3 == (*i));
        (*++i) = 7;
        assert( 7 == (*i));

        // Check the setting and reset these elements
        // to their original values.
        i = begin;
        assert( 4 == (*i));
        (*i) = 1;
        i++;
        assert( 3 == (*i));
        (*i++) = 2;
        assert( 3 == (*i));
        i++;
        assert( 7 == (*i));
        (*i) = 4;

        // Check the resetting.
        i = begin;
        int k = 1;
        do {
            assert( k == (*i));
            ++i;
            ++k;
        } while (i != begin);
    }
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_is_at_least_forward_category(begin);
        CGAL::Assert_is_at_least_forward_category(begin);
        typedef std::iterator_traits< Circulator >::value_type      VT;
        typedef std::iterator_traits< Circulator >::difference_type DT;
        assert(1==test_value_type(static_cast< VT* >(nullptr)));
        assert(1==test_distance_type(static_cast< DT* >(nullptr)));

        // Default constructor.
        Circulator z = Circulator();
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Circulator i = begin;

        // Check general support for circulators and iterators.
        assert( CGAL::is_empty_range( z, z));
        assert( ! CGAL::is_empty_range( i, begin));

        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                assert( k == (*i));
                su += (*i);
                ++k;
                Circulator j = ++i;
                assert( i ==  j);
                if ( i != begin) {
                    assert( (*i) == (*j));
                }
            } while (i != begin);  // Inequality and equality checked.
        }
        assert( i == begin);  // Equality checked.
        assert( su == 15);

        // Assignment.
        i = begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                assert( k == (*i));
                su += (*i);
                ++k;
                Circulator j = i++;
                assert(  i !=  j);
                if ( i != begin) {
                    assert( (*i) == (*j) + 1);
                }
            } while (i != begin);
        }
        assert( i == begin);
        assert( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_bidirectional_category(begin);
        CGAL::Assert_is_at_least_bidirectional_category(begin);
        // Loop backwards and pre-decrement.
        Circulator i = begin;
        int su = 0;
        int k  = 5;
        do {
            Circulator j = --i;
            assert(  i ==  j);
            assert( (*i) == (*j));
            assert( k == (*i));
            su += (*i);
            --k;
        } while (i != begin);
        assert( i == begin);
        assert( su == 15);

        // Assignment.
        i = begin;
        su = 0;
        k  = 5;
        // Loop with post-decrement.
        do {
            Circulator j = i--;
            assert(  i !=  j);
            if ( j != begin) {
                assert( (*i) == (*j) - 1);
            }
            assert( k == (*i));
            su += (*i);
            --k;
        } while (i != begin);
        assert( i == begin);
        assert( su == 15);
    }
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_is_at_least_forward_category(begin);
        CGAL::Assert_is_at_least_forward_category(begin);
        typedef std::iterator_traits< Circulator >::value_type      VT;
        typedef std::iterator_traits< Circulator >::difference_type DT;
        assert(1==test_value_type(static_cast< VT* >(nullptr)));
        assert(1==test_distance_type(static_cast< DT* >(nullptr)));

        // Default constructor.
        Circulator z = Circulator();
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Circulator i = begin;

        // Check general support for circulators and iterators.
        assert( CGAL::is_empty_range( z, z));
        assert( ! CGAL::is_empty_range( i, begin));

        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                assert( k == (*i));
                su += (*i);
                ++k;
                Circulator j = ++i;
                assert( i ==  j);
                if ( i != begin) {
                    assert( (*i) == (*j));
                }
            } while (i != begin);  // Inequality and equality checked.
        }
        assert( i == begin);  // Equality checked.
        assert( su == 15);

        // Assignment.
        i = begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                assert( k == (*i));
                su += (*i);
                ++k;
                Circulator j = i++;
                assert(  i !=  j);
                if ( i != begin) {
                    assert( (*i) == (*j) + 1);
                }
            } while (i != begin);
        }
        assert( i == begin);
        assert( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_bidirectional_category(begin);
        CGAL::Assert_is_at_least_bidirectional_category(begin);
        // Loop backwards and pre-decrement.
        Circulator i = begin;
        int su = 0;
        int k  = 5;
        do {
            Circulator j = --i;
            assert(  i ==  j);
            assert( (*i) == (*j));
            assert( k == (*i));
            su += (*i);
            --k;
        } while (i != begin);
        assert( i == begin);
        assert( su == 15);

        // Assignment.
        i = begin;
        su = 0;
        k  = 5;
        // Loop with post-decrement.
        do {
            Circulator j = i--;
            assert(  i !=  j);
            if ( j != begin) {
                assert( (*i) == (*j) - 1);
            }
            assert( k == (*i));
            su += (*i);
            --k;
        } while (i != begin);
        assert( i == begin);
        assert( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_random_access_category(begin);
        CGAL::Assert_is_at_least_random_access_category(begin);
        // Random access.
        int k;
        for( k = 0; k < 5; k++) {
            assert( 1+k == begin[k]);
        }

        int su = begin[0]
               + begin[1]
               + begin[2]
               + begin[3]
               + begin[4];
        assert( su == 15);

        // Jump around.
        Circulator i = begin;
        i += 3;
        assert( 4 == (*i));
        i -= 2;
        assert( 2 == (*i));
        i += 3;
        assert( 5 == (*i));
        i -= 4;
        assert( 1 == (*i));
        assert( i == begin);
        Circulator j = i + 3;
        assert( 4 == (*j));
        Circulator jj = j - 2;
        assert( 2 == (*jj));
        typedef std::ptrdiff_t PT;
        jj = PT(4) + jj;
        assert( jj == begin);
        Circulator ij = jj - 5;
        ij = jj - 5; // avoids warning with NDEBUG
        assert( ij == begin);

        // Difference test.
        assert( jj - i == 5  ||  jj - i == 0);
        assert( i + (j-i) == j);
        assert( (j-i) + i == j);
    }
    { // Open own scope to hide local variables.
        Circulator i = begin;
        i[2] = 18;
        i[4] = 9;
        i[3] = 12;
        assert( i[2] == 18);
        assert( i[4] == 9);
        assert( i[3] == 12);
        i[2] = 3;
        i[3] = 4;
        i[4] = 5;
        // Check the resetting.
        i = begin;
        int k = 1;
        do {
            assert( k == (*i));
            ++i;
            ++k;
        } while (i != begin);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_circulator( begin);
        CGAL::Assert_circulator( begin);

        // Check the local type parameters.
        Circulator::value_type      k1;
        k1 = 1;
        Circulator::reference       k2 = k1;
        (void)k2;
        assert( k2 == 1);
        Circulator::pointer         p1 = &k1;
        (void)p1;
        assert( (*p1) == 1);
        k1 = 3;
        assert( k1 == 3);
        assert( k2 == 3);
        assert( (*p1) == 3);
        k1 = 6;
        assert( k1 == 6);
        assert( k2 == 6);
        assert( (*p1) == 6);
        Circulator::size_type s = 5;
        assert( s == 5);
        Circulator::difference_type d = -5;
        assert( d == -5);

        // Check tests for empty data structures.
        Circulator z = Circulator();
        assert(   z == nullptr);
        assert( ! (z != nullptr));
        Circulator i = begin;
        assert( ! (i == nullptr));
        assert(   i != nullptr);
        assert( i == begin);
        assert( i == begin);
        // Do I reach myself.
        ++i;
        Circulator j = i;
        int k = 0;
        do {
            assert( k < 5);
            ++k;
            ++i;
        } while( i != j);
        assert( k == 5);
    }
    { // Open own scope to hide local variables.
        // Do I reach myself backwards.
        Circulator i = begin;
        ++i;
        Circulator j = i;
        int k = 0;
        do {
            assert( k < 5);
            ++k;
            --i;
        } while( i != j);
        assert( k == 5);
    }
    { // Open own scope to hide local variables.
        Circulator::difference_type d = begin - begin;
        assert( d == 0);
        d = begin - begin;
        assert( d == 0);
        Circulator i = begin + 1;
        assert( begin - i == 1 ||  begin - i == -1);
        assert( i - begin == 1 ||  i - begin == -1);
        // Check minimal circulator properties.
        i = i.min_circulator();
        Circulator j = i;
        assert( j - i == 0);
        j++;
        assert( j - i == 1);
        j++;
        assert( j - i == 2);
        j++;
        assert( j - i == 3);
        j++;
        assert( j - i == 4);
        j++;
        assert( j - i == 0);
    }

    Vector v2 = v;
    const Vector& v1 = v2;
    typedef Vector::const_iterator ConstIterBase;
    typedef Random_access_const_circulator_from_iterator<
      ConstIterBase, int, std::size_t, std::ptrdiff_t>
    ConstCircBase;
typedef Circulator_identity<ConstCircBase,const int&,
  const int*> C_Circulator;
C_Circulator c_begin(ConstCircBase(v1.begin(),v1.end()));
Assert_random_access_category(c_begin);
{ // Open own scope to hide local variables.
    // Check generally correct parameter properties.
    CGAL::Assert_circulator_or_iterator(c_begin);
    CGAL::Assert_circulator_or_iterator(c_begin);
    CGAL::Assert_is_at_least_forward_category(c_begin);
    CGAL::Assert_is_at_least_forward_category(c_begin);
    typedef std::iterator_traits< C_Circulator >::value_type      VT;
    typedef std::iterator_traits< C_Circulator >::difference_type DT;
    assert(1==test_value_type(static_cast< VT* >(nullptr)));
    assert(1==test_distance_type(static_cast< DT* >(nullptr)));

    // Default constructor.
    C_Circulator z = C_Circulator();
    CGAL::Assert_circulator_or_iterator(z);
    // Copy constructor.
    C_Circulator i = c_begin;

    // Check general support for circulators and iterators.
    assert( CGAL::is_empty_range( z, z));
    assert( ! CGAL::is_empty_range( i, c_begin));

    int su = 0;
    int k  = 1;
    // Check general loop, pre-increment, dereference.
    if (! CGAL::is_empty_range( i, c_begin)) {   // superfluous
        do {
            assert( k == (*i));
            su += (*i);
            ++k;
            C_Circulator j = ++i;
            assert( i ==  j);
            if ( i != c_begin) {
                assert( (*i) == (*j));
            }
        } while (i != c_begin);  // Inequality and equality checked.
    }
    assert( i == c_begin);  // Equality checked.
    assert( su == 15);

    // Assignment.
    i = c_begin;
    su = 0;
    k  = 1;
    // Loop with post increment.
    if (! CGAL::is_empty_range( i, c_begin)) {   // superfluous
        do {
            assert( k == (*i));
            su += (*i);
            ++k;
            C_Circulator j = i++;
            assert(  i !=  j);
            if ( i != c_begin) {
                assert( (*i) == (*j) + 1);
            }
        } while (i != c_begin);
    }
    assert( i == c_begin);
    assert( su == 15);
}
{ // Open own scope to hide local variables.
    CGAL::Assert_is_at_least_bidirectional_category(c_begin);
    CGAL::Assert_is_at_least_bidirectional_category(c_begin);
    // Loop backwards and pre-decrement.
    C_Circulator i = c_begin;
    int su = 0;
    int k  = 5;
    do {
        C_Circulator j = --i;
        assert(  i ==  j);
        assert( (*i) == (*j));
        assert( k == (*i));
        su += (*i);
        --k;
    } while (i != c_begin);
    assert( i == c_begin);
    assert( su == 15);

    // Assignment.
    i = c_begin;
    su = 0;
    k  = 5;
    // Loop with post-decrement.
    do {
        C_Circulator j = i--;
        assert(  i !=  j);
        if ( j != c_begin) {
            assert( (*i) == (*j) - 1);
        }
        assert( k == (*i));
        su += (*i);
        --k;
    } while (i != c_begin);
    assert( i == c_begin);
    assert( su == 15);
}
{ // Open own scope to hide local variables.
    CGAL::Assert_is_at_least_random_access_category(c_begin);
    CGAL::Assert_is_at_least_random_access_category(c_begin);
    // Random access.
    int k;
    for( k = 0; k < 5; k++) {
        assert( 1+k == c_begin[k]);
    }

      int su = c_begin[0]
             + c_begin[1]
             + c_begin[2]
             + c_begin[3]
             + c_begin[4];
    assert( su == 15);

    // Jump around.
    C_Circulator i = c_begin;
    i += 3;
    assert( 4 == (*i));
    i -= 2;
    assert( 2 == (*i));
    i += 3;
    assert( 5 == (*i));
    i -= 4;
    assert( 1 == (*i));
    assert( i == c_begin);
    C_Circulator j = i + 3;
    assert( 4 == (*j));
    C_Circulator jj = j - 2;
    assert( 2 == (*jj));
    typedef std::ptrdiff_t PT;
    jj = PT(4) + jj;
    assert( jj == c_begin);
    C_Circulator ij = jj - 5;
    ij = jj - 5; // avoids warning with NDEBUG
    assert( ij == c_begin);

    // Difference test.
    assert( jj - i == 5  ||  jj - i == 0);
    assert( i + (j-i) == j);
    assert( (j-i) + i == j);
}
{ // Open own scope to hide local variables.
    CGAL::Assert_circulator( c_begin);
    CGAL::Assert_circulator( c_begin);

    // Check the local type parameters.
    C_Circulator::value_type      k1;
    k1 = 1;
    C_Circulator::reference       k2 = k1;
    (void)k2;
    assert( k2 == 1);
    C_Circulator::pointer         p1 = &k1;
    (void)p1;
    assert( (*p1) == 1);
    k1 = 3;
    assert( k1 == 3);
    assert( k2 == 3);
    assert( (*p1) == 3);
    k1 = 6;
    assert( k1 == 6);
    assert( k2 == 6);
    assert( (*p1) == 6);
    C_Circulator::size_type s = 5;
    assert( s == 5);
    C_Circulator::difference_type d = -5;
    assert( d == -5);

    // Check tests for empty data structures.
    C_Circulator z = C_Circulator();
    assert(   z == nullptr);
    assert( ! (z != nullptr));
    C_Circulator i = c_begin;
    assert( ! (i == nullptr));
    assert(   i != nullptr);
    assert( i == c_begin);
    assert( i == c_begin);
    // Do I reach myself.
    ++i;
    C_Circulator j = i;
    int k = 0;
    do {
        assert( k < 5);
        ++k;
        ++i;
    } while( i != j);
    assert( k == 5);
}
{ // Open own scope to hide local variables.
    // Do I reach myself backwards.
    C_Circulator i = c_begin;
    ++i;
    C_Circulator j = i;
    int k = 0;
    do {
        assert( k < 5);
        ++k;
        --i;
    } while( i != j);
    assert( k == 5);
}
{ // Open own scope to hide local variables.
    C_Circulator::difference_type d = c_begin - c_begin;
    assert( d == 0);
    d = c_begin - c_begin;
    assert( d == 0);
    C_Circulator i = c_begin + 1;
    assert( c_begin - i == 1 ||  c_begin - i == -1);
    assert( i - c_begin == 1 ||  i - c_begin == -1);
    // Check minimal circulator properties.
    i = i.min_circulator();
    C_Circulator j = i;
    assert( j - i == 0);
    j++;
    assert( j - i == 1);
    j++;
    assert( j - i == 2);
    j++;
    assert( j - i == 3);
    j++;
    assert( j - i == 4);
    j++;
    assert( j - i == 0);
}
  }
}
void test_Iterator_project()
{
  {
    typedef In_place_list<item,false> List;
    List l;
    l.push_back( *new item(1));
    l.push_back( *new item(2));
    l.push_back( *new item(3));
    l.push_back( *new item(4));
    l.push_back( *new item(5));
    assert( l.size() == 5);
    typedef Identity<item> Ident;
    typedef List::iterator IterBase;
    typedef Iterator_project<IterBase,Ident,item&,item*,std::ptrdiff_t,
      std::bidirectional_iterator_tag> Iterator;
    Iterator begin(l.begin());
    Iterator end(l.end());
    Assert_bidirectional_category(begin);
    Assert_bidirectional_category(end);
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_circulator_or_iterator(end);
        CGAL::Assert_is_at_least_forward_category(begin);
        CGAL::Assert_is_at_least_forward_category(end);
        typedef std::iterator_traits< Iterator >::value_type      VT;
        typedef std::iterator_traits< Iterator >::difference_type DT;
        assert(1==test_value_type(static_cast< VT* >(nullptr)));
        assert(1==test_distance_type(static_cast< DT* >(nullptr)));

        // Default constructor.
        Iterator z ;
        CGAL::Assert_circulator_or_iterator(z);
        z = end ;
        // Copy constructor.
        Iterator i = begin;

        // Check general support for circulators and iterators.
        assert( CGAL::is_empty_range( z, z));
        assert( ! CGAL::is_empty_range( i, end));

        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, end)) {   // superfluous
            do {
                assert( k == (*i).key);
                su += (*i).key;
                ++k;
                Iterator j = ++i;
                assert( i ==  j);
                if ( i != end) {
                    assert( (*i).key == (*j).key);
                }
            } while (i != end);  // Inequality and equality checked.
        }
        assert( i == end);  // Equality checked.
        assert( su == 15);

        // Assignment.
        i = begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, end)) {   // superfluous
            do {
                assert( k == (*i).key);
                su += (*i).key;
                ++k;
                Iterator j = i++;
                assert(  i !=  j);
                if ( i != end) {
                    assert( (*i).key == (*j).key + 1);
                }
            } while (i != end);
        }
        assert( i == end);
        assert( su == 15);
    }
    { // Open own scope to hide local variables.
        // Change three elements and check post-/pre-increment.
        Iterator i = begin;
        (*i++).key = 4;
        assert( 4 == (*begin).key);
        assert( 2 == (*i).key);
        (*i++).key = 3;
        assert( 3 == (*i).key);
        (*++i).key = 7;
        assert( 7 == (*i).key);

        // Check the setting and reset these elements
        // to their original values.
        i = begin;
        assert( 4 == (*i).key);
        (*i).key = 1;
        i++;
        assert( 3 == (*i).key);
        (*i++).key = 2;
        assert( 3 == (*i).key);
        i++;
        assert( 7 == (*i).key);
        (*i).key = 4;

        // Check the resetting.
        i = begin;
        int k = 1;
        do {
            assert( k == (*i).key);
            ++i;
            ++k;
        } while (i != end);
    }
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_circulator_or_iterator(end);
        CGAL::Assert_is_at_least_forward_category(begin);
        CGAL::Assert_is_at_least_forward_category(end);
        typedef std::iterator_traits< Iterator >::value_type      VT;
        typedef std::iterator_traits< Iterator >::difference_type DT;
        assert(1==test_value_type(static_cast< VT* >(nullptr)));
        assert(1==test_distance_type(static_cast< DT* >(nullptr)));

        // Default constructor.
        Iterator z ;
        CGAL::Assert_circulator_or_iterator(z);
        z = end ;
        // Copy constructor.
        Iterator i = begin;

        // Check general support for circulators and iterators.
        assert( CGAL::is_empty_range( z, z));
        assert( ! CGAL::is_empty_range( i, end));

        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, end)) {   // superfluous
            do {
                assert( k == (*i).key);
                su += (*i).key;
                ++k;
                Iterator j = ++i;
                assert( i ==  j);
                if ( i != end) {
                    assert( (*i).key == (*j).key);
                }
            } while (i != end);  // Inequality and equality checked.
        }
        assert( i == end);  // Equality checked.
        assert( su == 15);

        // Assignment.
        i = begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, end)) {   // superfluous
            do {
                assert( k == (*i).key);
                su += (*i).key;
                ++k;
                Iterator j = i++;
                assert(  i !=  j);
                if ( i != end) {
                    assert( (*i).key == (*j).key + 1);
                }
            } while (i != end);
        }
        assert( i == end);
        assert( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_bidirectional_category(begin);
        CGAL::Assert_is_at_least_bidirectional_category(end);
        // Loop backwards and pre-decrement.
        Iterator i = end;
        int su = 0;
        int k  = 5;
        do {
            Iterator j = --i;
            assert(  i ==  j);
            assert( (*i).key == (*j).key);
            assert( k == (*i).key);
            su += (*i).key;
            --k;
        } while (i != begin);
        assert( i == begin);
        assert( su == 15);

        // Assignment.
        i = end;
        su = 0;
        k  = 5;
        // Loop with post-decrement.
        do {
            Iterator j = i--;
            assert(  i !=  j);
            if ( j != end) {
                assert( (*i).key == (*j).key - 1);
            }
            assert( k == (*i).key);
            su += (*i).key;
            --k;
        } while (i != begin);
        assert( i == begin);
        assert( su == 15);
    }
    CGAL::Assert_iterator( begin);
    CGAL::Assert_iterator( end);

    List l2 = l;
    const List& l1 = l2;
    typedef List::const_iterator ConstIterBase;
    typedef Iterator_project< ConstIterBase,Ident>
    C_Iterator;
    C_Iterator c_begin(l1.begin());
    C_Iterator c_end(l1.end());
    Assert_bidirectional_category(c_begin);
    Assert_bidirectional_category(c_end);
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(c_begin);
        CGAL::Assert_circulator_or_iterator(c_end);
        CGAL::Assert_is_at_least_forward_category(c_begin);
        CGAL::Assert_is_at_least_forward_category(c_end);
        typedef std::iterator_traits< C_Iterator >::value_type      VT;
        typedef std::iterator_traits< C_Iterator >::difference_type DT;
        assert(1==test_value_type(static_cast< VT* >(nullptr)));
        assert(1==test_distance_type(static_cast< DT* >(nullptr)));

        // Default constructor.
        C_Iterator z ;
        CGAL::Assert_circulator_or_iterator(z);
        z = c_end ;
        // Copy constructor.
        C_Iterator i = c_begin;

        // Check general support for circulators and iterators.
        assert( CGAL::is_empty_range( z, z));
        assert( ! CGAL::is_empty_range( i, c_end));

        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, c_end)) {   // superfluous
            do {
                assert( k == (*i).key);
                su += (*i).key;
                ++k;
                C_Iterator j = ++i;
                assert( i ==  j);
                if ( i != c_end) {
                    assert( (*i).key == (*j).key);
                }
            } while (i != c_end);  // Inequality and equality checked.
        }
        assert( i == c_end);  // Equality checked.
        assert( su == 15);

        // Assignment.
        i = c_begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, c_end)) {   // superfluous
            do {
                assert( k == (*i).key);
                su += (*i).key;
                ++k;
                C_Iterator j = i++;
                assert(  i !=  j);
                if ( i != c_end) {
                    assert( (*i).key == (*j).key + 1);
                }
            } while (i != c_end);
        }
        assert( i == c_end);
        assert( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_bidirectional_category(c_begin);
        CGAL::Assert_is_at_least_bidirectional_category(c_end);
        // Loop backwards and pre-decrement.
        C_Iterator i = c_end;
        int su = 0;
        int k  = 5;
        do {
            C_Iterator j = --i;
            assert(  i ==  j);
            assert( (*i).key == (*j).key);
            assert( k == (*i).key);
            su += (*i).key;
            --k;
        } while (i != c_begin);
        assert( i == c_begin);
        assert( su == 15);

        // Assignment.
        i = c_end;
        su = 0;
        k  = 5;
        // Loop with post-decrement.
        do {
            C_Iterator j = i--;
            assert(  i !=  j);
            if ( j != c_end) {
                assert( (*i).key == (*j).key - 1);
            }
            assert( k == (*i).key);
            su += (*i).key;
            --k;
        } while (i != c_begin);
        assert( i == c_begin);
        assert( su == 15);
    }
    CGAL::Assert_iterator( c_begin);
    CGAL::Assert_iterator( c_end);
    l.destroy();
    l2.destroy();
    }

    //
    //
    //============================================
    //
    //

    {
    typedef In_place_list<item,false> List;
    List l;
    l.push_back( *new item(1));
    l.push_back( *new item(2));
    l.push_back( *new item(3));
    l.push_back( *new item(4));
    l.push_back( *new item(5));
    assert( l.size() == 5);
    typedef Cast_function_object<item,item> Ident;
    typedef List::iterator IterBase;
    typedef Iterator_project<IterBase,Ident> Iterator;
    Iterator begin(l.begin());
    Iterator end(l.end());
    Assert_bidirectional_category(begin);
    Assert_bidirectional_category(end);
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_circulator_or_iterator(end);
        CGAL::Assert_is_at_least_forward_category(begin);
        CGAL::Assert_is_at_least_forward_category(end);
        typedef std::iterator_traits< Iterator >::value_type      VT;
        typedef std::iterator_traits< Iterator >::difference_type DT;
        assert(1==test_value_type(static_cast< VT* >(nullptr)));
        assert(1==test_distance_type(static_cast< DT* >(nullptr)));

        // Default constructor.
        Iterator z ;
        CGAL::Assert_circulator_or_iterator(z);
        z = end ;
        // Copy constructor.
        Iterator i = begin;

        // Check general support for circulators and iterators.
        assert( CGAL::is_empty_range( z, z));
        assert( ! CGAL::is_empty_range( i, end));

        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, end)) {   // superfluous
            do {
                assert( k == (*i).key);
                su += (*i).key;
                ++k;
                Iterator j = ++i;
                assert( i ==  j);
                if ( i != end) {
                    assert( (*i).key == (*j).key);
                }
            } while (i != end);  // Inequality and equality checked.
        }
        assert( i == end);  // Equality checked.
        assert( su == 15);

        // Assignment.
        i = begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, end)) {   // superfluous
            do {
                assert( k == (*i).key);
                su += (*i).key;
                ++k;
                Iterator j = i++;
                assert(  i !=  j);
                if ( i != end) {
                    assert( (*i).key == (*j).key + 1);
                }
            } while (i != end);
        }
        assert( i == end);
        assert( su == 15);
    }
    { // Open own scope to hide local variables.
        // Change three elements and check post-/pre-increment.
        Iterator i = begin;
        (*i++).key = 4;
        assert( 4 == (*begin).key);
        assert( 2 == (*i).key);
        (*i++).key = 3;
        assert( 3 == (*i).key);
        (*++i).key = 7;
        assert( 7 == (*i).key);

        // Check the setting and reset these elements
        // to their original values.
        i = begin;
        assert( 4 == (*i).key);
        (*i).key = 1;
        i++;
        assert( 3 == (*i).key);
        (*i++).key = 2;
        assert( 3 == (*i).key);
        i++;
        assert( 7 == (*i).key);
        (*i).key = 4;

        // Check the resetting.
        i = begin;
        int k = 1;
        do {
            assert( k == (*i).key);
            ++i;
            ++k;
        } while (i != end);
    }
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_circulator_or_iterator(end);
        CGAL::Assert_is_at_least_forward_category(begin);
        CGAL::Assert_is_at_least_forward_category(end);
        typedef std::iterator_traits< Iterator >::value_type      VT;
        typedef std::iterator_traits< Iterator >::difference_type DT;
        assert(1==test_value_type(static_cast< VT* >(nullptr)));
        assert(1==test_distance_type(static_cast< DT* >(nullptr)));

        // Default constructor.
        Iterator z ;
        CGAL::Assert_circulator_or_iterator(z);
        z = end ;
        // Copy constructor.
        Iterator i = begin;

        // Check general support for circulators and iterators.
        assert( CGAL::is_empty_range( z, z));
        assert( ! CGAL::is_empty_range( i, end));

        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, end)) {   // superfluous
            do {
                assert( k == (*i).key);
                su += (*i).key;
                ++k;
                Iterator j = ++i;
                assert( i ==  j);
                if ( i != end) {
                    assert( (*i).key == (*j).key);
                }
            } while (i != end);  // Inequality and equality checked.
        }
        assert( i == end);  // Equality checked.
        assert( su == 15);

        // Assignment.
        i = begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, end)) {   // superfluous
            do {
                assert( k == (*i).key);
                su += (*i).key;
                ++k;
                Iterator j = i++;
                assert(  i !=  j);
                if ( i != end) {
                    assert( (*i).key == (*j).key + 1);
                }
            } while (i != end);
        }
        assert( i == end);
        assert( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_bidirectional_category(begin);
        CGAL::Assert_is_at_least_bidirectional_category(end);
        // Loop backwards and pre-decrement.
        Iterator i = end;
        int su = 0;
        int k  = 5;
        do {
            Iterator j = --i;
            assert(  i ==  j);
            assert( (*i).key == (*j).key);
            assert( k == (*i).key);
            su += (*i).key;
            --k;
        } while (i != begin);
        assert( i == begin);
        assert( su == 15);

        // Assignment.
        i = end;
        su = 0;
        k  = 5;
        // Loop with post-decrement.
        do {
            Iterator j = i--;
            assert(  i !=  j);
            if ( j != end) {
                assert( (*i).key == (*j).key - 1);
            }
            assert( k == (*i).key);
            su += (*i).key;
            --k;
        } while (i != begin);
        assert( i == begin);
        assert( su == 15);
    }
    CGAL::Assert_iterator( begin);
    CGAL::Assert_iterator( end);

    List l2 = l;
    const List& l1 = l2;
    typedef List::const_iterator ConstIterBase;
    typedef Iterator_project< ConstIterBase,Ident,const item&,
      const item*, std::ptrdiff_t,
                              std::bidirectional_iterator_tag>
    C_Iterator;
    C_Iterator c_begin(l1.begin());
    C_Iterator c_end(l1.end());
    Assert_bidirectional_category(c_begin);
    Assert_bidirectional_category(c_end);
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(c_begin);
        CGAL::Assert_circulator_or_iterator(c_end);
        CGAL::Assert_is_at_least_forward_category(c_begin);
        CGAL::Assert_is_at_least_forward_category(c_end);
        typedef std::iterator_traits< C_Iterator >::value_type      VT;
        typedef std::iterator_traits< C_Iterator >::difference_type DT;
        assert(1==test_value_type(static_cast< VT* >(nullptr)));
        assert(1==test_distance_type(static_cast< DT* >(nullptr)));

        // Default constructor.
        C_Iterator z ;
        CGAL::Assert_circulator_or_iterator(z);
        z = c_end ;
        // Copy constructor.
        C_Iterator i = c_begin;

        // Check general support for circulators and iterators.
        assert( CGAL::is_empty_range( z, z));
        assert( ! CGAL::is_empty_range( i, c_end));

        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, c_end)) {   // superfluous
            do {
                assert( k == (*i).key);
                su += (*i).key;
                ++k;
                C_Iterator j = ++i;
                assert( i ==  j);
                if ( i != c_end) {
                    assert( (*i).key == (*j).key);
                }
            } while (i != c_end);  // Inequality and equality checked.
        }
        assert( i == c_end);  // Equality checked.
        assert( su == 15);

        // Assignment.
        i = c_begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, c_end)) {   // superfluous
            do {
                assert( k == (*i).key);
                su += (*i).key;
                ++k;
                C_Iterator j = i++;
                assert(  i !=  j);
                if ( i != c_end) {
                    assert( (*i).key == (*j).key + 1);
                }
            } while (i != c_end);
        }
        assert( i == c_end);
        assert( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_bidirectional_category(c_begin);
        CGAL::Assert_is_at_least_bidirectional_category(c_end);
        // Loop backwards and pre-decrement.
        C_Iterator i = c_end;
        int su = 0;
        int k  = 5;
        do {
            C_Iterator j = --i;
            assert(  i ==  j);
            assert( (*i).key == (*j).key);
            assert( k == (*i).key);
            su += (*i).key;
            --k;
        } while (i != c_begin);
        assert( i == c_begin);
        assert( su == 15);

        // Assignment.
        i = c_end;
        su = 0;
        k  = 5;
        // Loop with post-decrement.
        do {
            C_Iterator j = i--;
            assert(  i !=  j);
            if ( j != c_end) {
                assert( (*i).key == (*j).key - 1);
            }
            assert( k == (*i).key);
            su += (*i).key;
            --k;
        } while (i != c_begin);
        assert( i == c_begin);
        assert( su == 15);
    }
    CGAL::Assert_iterator( c_begin);
    CGAL::Assert_iterator( c_end);
    l.destroy();
    l2.destroy();
    }

    //
    //
    //============================================
    //
    //

    {
    typedef std::vector<int> Vector;
    Vector v;
    v.push_back(1);
    v.push_back(2);
    v.push_back(3);
    v.push_back(4);
    v.push_back(5);
    assert( v.size() == 5);
    typedef Vector::iterator IterBase;
    typedef Identity<int> Ident;
    typedef Iterator_project<IterBase,Ident,int&,int*,std::ptrdiff_t,
      std::random_access_iterator_tag> Iterator;
    Iterator begin(v.begin());
    Iterator end(v.end());
    Assert_random_access_category(begin);
    Assert_random_access_category(end);
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_circulator_or_iterator(end);
        CGAL::Assert_is_at_least_forward_category(begin);
        CGAL::Assert_is_at_least_forward_category(end);
        typedef std::iterator_traits< Iterator >::value_type      VT;
        typedef std::iterator_traits< Iterator >::difference_type DT;
        assert(1==test_value_type(static_cast< VT* >(nullptr)));
        assert(1==test_distance_type(static_cast< DT* >(nullptr)));

        // Default constructor.
        Iterator z ;
        CGAL::Assert_circulator_or_iterator(z);
        z = end ;
        // Copy constructor.
        Iterator i = begin;

        // Check general support for circulators and iterators.
        assert( CGAL::is_empty_range( z, z));
        assert( ! CGAL::is_empty_range( i, end));

        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, end)) {   // superfluous
            do {
                assert( k == (*i));
                su += (*i);
                ++k;
                Iterator j = ++i;
                assert( i ==  j);
                if ( i != end) {
                    assert( (*i) == (*j));
                }
            } while (i != end);  // Inequality and equality checked.
        }
        assert( i == end);  // Equality checked.
        assert( su == 15);

        // Assignment.
        i = begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, end)) {   // superfluous
            do {
                assert( k == (*i));
                su += (*i);
                ++k;
                Iterator j = i++;
                assert(  i !=  j);
                if ( i != end) {
                    assert( (*i) == (*j) + 1);
                }
            } while (i != end);
        }
        assert( i == end);
        assert( su == 15);
    }
    { // Open own scope to hide local variables.
        // Change three elements and check post-/pre-increment.
        Iterator i = begin;
        (*i++) = 4;
        assert( 4 == (*begin));
        assert( 2 == (*i));
        (*i++) = 3;
        assert( 3 == (*i));
        (*++i) = 7;
        assert( 7 == (*i));

        // Check the setting and reset these elements
        // to their original values.
        i = begin;
        assert( 4 == (*i));
        (*i) = 1;
        i++;
        assert( 3 == (*i));
        (*i++) = 2;
        assert( 3 == (*i));
        i++;
        assert( 7 == (*i));
        (*i) = 4;

        // Check the resetting.
        i = begin;
        int k = 1;
        do {
            assert( k == (*i));
            ++i;
            ++k;
        } while (i != end);
    }
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_circulator_or_iterator(end);
        CGAL::Assert_is_at_least_forward_category(begin);
        CGAL::Assert_is_at_least_forward_category(end);
        typedef std::iterator_traits< Iterator >::value_type      VT;
        typedef std::iterator_traits< Iterator >::difference_type DT;
        assert(1==test_value_type(static_cast< VT* >(nullptr)));
        assert(1==test_distance_type(static_cast< DT* >(nullptr)));

        // Default constructor.
        Iterator z ;
        CGAL::Assert_circulator_or_iterator(z);
        z = end ;
        // Copy constructor.
        Iterator i = begin;


        // Check general support for circulators and iterators.
        assert( CGAL::is_empty_range( z, z));
        assert( ! CGAL::is_empty_range( i, end));

        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, end)) {   // superfluous
            do {
                assert( k == (*i));
                su += (*i);
                ++k;
                Iterator j = ++i;
                assert( i ==  j);
                if ( i != end) {
                    assert( (*i) == (*j));
                }
            } while (i != end);  // Inequality and equality checked.
        }
        assert( i == end);  // Equality checked.
        assert( su == 15);

        // Assignment.
        i = begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, end)) {   // superfluous
            do {
                assert( k == (*i));
                su += (*i);
                ++k;
                Iterator j = i++;
                assert(  i !=  j);
                if ( i != end) {
                    assert( (*i) == (*j) + 1);
                }
            } while (i != end);
        }
        assert( i == end);
        assert( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_bidirectional_category(begin);
        CGAL::Assert_is_at_least_bidirectional_category(end);
        // Loop backwards and pre-decrement.
        Iterator i = end;
        int su = 0;
        int k  = 5;
        do {
            Iterator j = --i;
            assert(  i ==  j);
            assert( (*i) == (*j));
            assert( k == (*i));
            su += (*i);
            --k;
        } while (i != begin);
        assert( i == begin);
        assert( su == 15);

        // Assignment.
        i = end;
        su = 0;
        k  = 5;
        // Loop with post-decrement.
        do {
            Iterator j = i--;
            assert(  i !=  j);
            if ( j != end) {
                assert( (*i) == (*j) - 1);
            }
            assert( k == (*i));
            su += (*i);
            --k;
        } while (i != begin);
        assert( i == begin);
        assert( su == 15);
    }
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_circulator_or_iterator(end);
        CGAL::Assert_is_at_least_forward_category(begin);
        CGAL::Assert_is_at_least_forward_category(end);
        typedef std::iterator_traits< Iterator >::value_type      VT;
        typedef std::iterator_traits< Iterator >::difference_type DT;
        assert(1==test_value_type(static_cast< VT* >(nullptr)));
        assert(1==test_distance_type(static_cast< DT* >(nullptr)));

        // Default constructor.
        Iterator z ;
        CGAL::Assert_circulator_or_iterator(z);
        z = end ;
        // Copy constructor.
        Iterator i = begin;

        // Check general support for circulators and iterators.
        assert( CGAL::is_empty_range( z, z));
        assert( ! CGAL::is_empty_range( i, end));

        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, end)) {   // superfluous
            do {
                assert( k == (*i));
                su += (*i);
                ++k;
                Iterator j = ++i;
                assert( i ==  j);
                if ( i != end) {
                    assert( (*i) == (*j));
                }
            } while (i != end);  // Inequality and equality checked.
        }
        assert( i == end);  // Equality checked.
        assert( su == 15);

        // Assignment.
        i = begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, end)) {   // superfluous
            do {
                assert( k == (*i));
                su += (*i);
                ++k;
                Iterator j = i++;
                assert(  i !=  j);
                if ( i != end) {
                    assert( (*i) == (*j) + 1);
                }
            } while (i != end);
        }
        assert( i == end);
        assert( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_bidirectional_category(begin);
        CGAL::Assert_is_at_least_bidirectional_category(end);
        // Loop backwards and pre-decrement.
        Iterator i = end;
        int su = 0;
        int k  = 5;
        do {
            Iterator j = --i;
            assert(  i ==  j);
            assert( (*i) == (*j));
            assert( k == (*i));
            su += (*i);
            --k;
        } while (i != begin);
        assert( i == begin);
        assert( su == 15);

        // Assignment.
        i = end;
        su = 0;
        k  = 5;
        // Loop with post-decrement.
        do {
            Iterator j = i--;
            assert(  i !=  j);
            if ( j != end) {
                assert( (*i) == (*j) - 1);
            }
            assert( k == (*i));
            su += (*i);
            --k;
        } while (i != begin);
        assert( i == begin);
        assert( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_random_access_category(begin);
        CGAL::Assert_is_at_least_random_access_category(end);
        // Random access.
        int k;
        for( k = 0; k < 5; k++) {
            assert( 1+k == begin[k]);
        }

          int su = begin[0]
                 + begin[1]
                 + begin[2]
                 + begin[3]
                 + begin[4];
        assert( su == 15);

        // Jump around.
        Iterator i = begin;
        i += 3;
        assert( 4 == (*i));
        i -= 2;
        assert( 2 == (*i));
        i += 3;
        assert( 5 == (*i));
        i -= 4;
        assert( 1 == (*i));
        assert( i == begin);
        Iterator j = i + 3;
        assert( 4 == (*j));
        Iterator jj = j - 2;
        assert( 2 == (*jj));
        typedef std::ptrdiff_t PT;
        jj = PT(4) + jj;
        assert( jj == end);
        Iterator ij = jj - 5;
        ij = jj - 5; // avoids warning with NDEBUG
        assert( ij == begin);

        // Difference test.
        assert( jj - i == 5  ||  jj - i == 0);
        assert( i + (j-i) == j);
        assert( (j-i) + i == j);
    }
    { // Open own scope to hide local variables.
        Iterator i = begin;
        i[2] = 18;
        i[4] = 9;
        i[3] = 12;
        assert( i[2] == 18);
        assert( i[4] == 9);
        assert( i[3] == 12);
        i[2] = 3;
        i[3] = 4;
        i[4] = 5;
        // Check the resetting.
        i = begin;
        int k = 1;
        do {
            assert( k == (*i));
            ++i;
            ++k;
        } while (i != end);
    }
    CGAL::Assert_iterator( begin);
    CGAL::Assert_iterator( end);
    { // Open own scope to hide local variables.
        assert( end - begin ==  5);
        assert( begin - end == -5);
        // Relational operator.
        Iterator i = begin;
        ++i;
        Iterator j = i;
        ++j;
        assert( begin < i);
        assert( i < j);
        assert( j < end);
        assert( j > i);
        assert( i <= j);
        assert( j >= i);
        assert( i <= i);
        assert( i >= i);

        assert( !( i >= j));
        assert( !( j <= i));
        assert( !( i > j));
        assert( !( j < i));
        assert( !( i > i));
        assert( !( i < i));
    }

    Vector v2 = v;
    const Vector& v1 = v2;
    typedef Vector::const_iterator ConstIterBase;
    typedef Iterator_project< ConstIterBase,Ident,const int&,
      const int*, std::ptrdiff_t,
                              std::random_access_iterator_tag>
    C_Iterator;
    C_Iterator c_begin(v1.begin());
    C_Iterator c_end(v1.end());
    Assert_random_access_category(c_begin);
    Assert_random_access_category(c_end);
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(c_begin);
        CGAL::Assert_circulator_or_iterator(c_end);
        CGAL::Assert_is_at_least_forward_category(c_begin);
        CGAL::Assert_is_at_least_forward_category(c_end);
        typedef std::iterator_traits< C_Iterator >::value_type      VT;
        typedef std::iterator_traits< C_Iterator >::difference_type DT;
        assert(1==test_value_type(static_cast< VT* >(nullptr)));
        assert(1==test_distance_type(static_cast< DT* >(nullptr)));

        // Default constructor.
        C_Iterator z ;
        CGAL::Assert_circulator_or_iterator(z);
        z = c_end ;
        // Copy constructor.
        C_Iterator i = c_begin;

        // Check general support for circulators and iterators.
        assert( CGAL::is_empty_range( z, z));
        assert( ! CGAL::is_empty_range( i, c_end));

        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, c_end)) {   // superfluous
            do {
                assert( k == (*i));
                su += (*i);
                ++k;
                C_Iterator j = ++i;
                assert( i ==  j);
                if ( i != c_end) {
                    assert( (*i) == (*j));
                }
            } while (i != c_end);  // Inequality and equality checked.
        }
        assert( i == c_end);  // Equality checked.
        assert( su == 15);

        // Assignment.
        i = c_begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, c_end)) {   // superfluous
            do {
                assert( k == (*i));
                su += (*i);
                ++k;
                C_Iterator j = i++;
                assert(  i !=  j);
                if ( i != c_end) {
                    assert( (*i) == (*j) + 1);
                }
            } while (i != c_end);
        }
        assert( i == c_end);
        assert( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_bidirectional_category(c_begin);
        CGAL::Assert_is_at_least_bidirectional_category(c_end);
        // Loop backwards and pre-decrement.
        C_Iterator i = c_end;
        int su = 0;
        int k  = 5;
        do {
            C_Iterator j = --i;
            assert(  i ==  j);
            assert( (*i) == (*j));
            assert( k == (*i));
            su += (*i);
            --k;
        } while (i != c_begin);
        assert( i == c_begin);
        assert( su == 15);

        // Assignment.
        i = c_end;
        su = 0;
        k  = 5;
        // Loop with post-decrement.
        do {
            C_Iterator j = i--;
            assert(  i !=  j);
            if ( j != c_end) {
                assert( (*i) == (*j) - 1);
            }
            assert( k == (*i));
            su += (*i);
            --k;
        } while (i != c_begin);
        assert( i == c_begin);
        assert( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_random_access_category(c_begin);
        CGAL::Assert_is_at_least_random_access_category(c_end);
        // Random access.
        int k;
        for( k = 0; k < 5; k++) {
            assert( 1+k == c_begin[k]);
        }

          int su = c_begin[0]
                 + c_begin[1]
                 + c_begin[2]
                 + c_begin[3]
                 + c_begin[4];
        assert( su == 15);

        // Jump around.
        C_Iterator i = c_begin;
        i += 3;
        assert( 4 == (*i));
        i -= 2;
        assert( 2 == (*i));
        i += 3;
        assert( 5 == (*i));
        i -= 4;
        assert( 1 == (*i));
        assert( i == c_begin);
        C_Iterator j = i + 3;
        assert( 4 == (*j));
        C_Iterator jj = j - 2;
        assert( 2 == (*jj));
        typedef std::ptrdiff_t PT;
        jj = PT(4) + jj;
        assert( jj == c_end);
        C_Iterator ij = jj - 5;
        ij = jj - 5; // avoids warning with NDEBUG
        assert( ij == c_begin);

        // Difference test.
        assert( jj - i == 5  ||  jj - i == 0);
        assert( i + (j-i) == j);
        assert( (j-i) + i == j);
    }
    CGAL::Assert_iterator( c_begin);
    CGAL::Assert_iterator( c_end);
    { // Open own scope to hide local variables.
        assert( c_end - c_begin ==  5);
        assert( c_begin - c_end == -5);
        // Relational operator.
        C_Iterator i = c_begin;
        ++i;
        C_Iterator j = i;
        ++j;
        assert( c_begin < i);
        assert( i < j);
        assert( j < c_end);
        assert( j > i);
        assert( i <= j);
        assert( j >= i);
        assert( i <= i);
        assert( i >= i);

        assert( !( i >= j));
        assert( !( j <= i));
        assert( !( i > j));
        assert( !( j < i));
        assert( !( i > i));
        assert( !( i < i));
    }
  }

    //
    //
    //============================================
    //
    //

    {
    typedef std::list<void*> PtrList;
    PtrList l;
    l.push_back( new item(1));
    l.push_back( new item(2));
    l.push_back( new item(3));
    l.push_back( new item(4));
    l.push_back( new item(5));
    assert( l.size() == 5);
    typedef Cast_function_object<void*,item*> Ident;
    typedef PtrList::iterator IterBase;
    typedef Iterator_project<IterBase,Ident> Iterator;
    Iterator begin(l.begin());
    Iterator end(l.end());
    Assert_bidirectional_category(begin);
    Assert_bidirectional_category(end);
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_circulator_or_iterator(end);
        CGAL::Assert_is_at_least_forward_category(begin);
        CGAL::Assert_is_at_least_forward_category(end);
        typedef std::iterator_traits< Iterator >::value_type      VT;
        typedef std::iterator_traits< Iterator >::difference_type DT;
        assert(1==test_value_type(static_cast< VT* >(nullptr)));
        assert(1==test_distance_type(static_cast< DT* >(nullptr)));

        // Default constructor.
        Iterator z ;
        CGAL::Assert_circulator_or_iterator(z);
        z = end ;
        // Copy constructor.
        Iterator i = begin;

        // Check general support for circulators and iterators.
        assert( CGAL::is_empty_range( z, z));
        assert( ! CGAL::is_empty_range( i, end));

        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, end)) {   // superfluous
            do {
                assert( k == (*i)->key);
                su += (*i)->key;
                ++k;
                Iterator j = ++i;
                assert( i ==  j);
                if ( i != end) {
                    assert( (*i)->key == (*j)->key);
                }
            } while (i != end);  // Inequality and equality checked.
        }
        assert( i == end);  // Equality checked.
        assert( su == 15);

        // Assignment.
        i = begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, end)) {   // superfluous
            do {
                assert( k == (*i)->key);
                su += (*i)->key;
                ++k;
                Iterator j = i++;
                assert(  i !=  j);
                if ( i != end) {
                    assert( (*i)->key == (*j)->key + 1);
                }
            } while (i != end);
        }
        assert( i == end);
        assert( su == 15);
    }
    { // Open own scope to hide local variables.
        // Change three elements and check post-/pre-increment.
        Iterator i = begin;
        (*i++)->key = 4;
        assert( 4 == (*begin)->key);
        assert( 2 == (*i)->key);
        (*i++)->key = 3;
        assert( 3 == (*i)->key);
        (*++i)->key = 7;
        assert( 7 == (*i)->key);

        // Check the setting and reset these elements
        // to their original values.
        i = begin;
        assert( 4 == (*i)->key);
        (*i)->key = 1;
        i++;
        assert( 3 == (*i)->key);
        (*i++)->key = 2;
        assert( 3 == (*i)->key);
        i++;
        assert( 7 == (*i)->key);
        (*i)->key = 4;

        // Check the resetting.
        i = begin;
        int k = 1;
        do {
            assert( k == (*i)->key);
            ++i;
            ++k;
        } while (i != end);
    }
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_circulator_or_iterator(end);
        CGAL::Assert_is_at_least_forward_category(begin);
        CGAL::Assert_is_at_least_forward_category(end);
        typedef std::iterator_traits< Iterator >::value_type      VT;
        typedef std::iterator_traits< Iterator >::difference_type DT;
        assert(1==test_value_type(static_cast< VT* >(nullptr)));
        assert(1==test_distance_type(static_cast< DT* >(nullptr)));

        // Default constructor.
        Iterator z ;
        CGAL::Assert_circulator_or_iterator(z);
        z = end ;
        // Copy constructor.
        Iterator i = begin;

        // Check general support for circulators and iterators.
        assert( CGAL::is_empty_range( z, z));
        assert( ! CGAL::is_empty_range( i, end));

        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, end)) {   // superfluous
            do {
                assert( k == (*i)->key);
                su += (*i)->key;
                ++k;
                Iterator j = ++i;
                assert( i ==  j);
                if ( i != end) {
                    assert( (*i)->key == (*j)->key);
                }
            } while (i != end);  // Inequality and equality checked.
        }
        assert( i == end);  // Equality checked.
        assert( su == 15);

        // Assignment.
        i = begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, end)) {   // superfluous
            do {
                assert( k == (*i)->key);
                su += (*i)->key;
                ++k;
                Iterator j = i++;
                assert(  i !=  j);
                if ( i != end) {
                    assert( (*i)->key == (*j)->key + 1);
                }
            } while (i != end);
        }
        assert( i == end);
        assert( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_bidirectional_category(begin);
        CGAL::Assert_is_at_least_bidirectional_category(end);
        // Loop backwards and pre-decrement.
        Iterator i = end;
        int su = 0;
        int k  = 5;
        do {
            Iterator j = --i;
            assert(  i ==  j);
            assert( (*i)->key == (*j)->key);
            assert( k == (*i)->key);
            su += (*i)->key;
            --k;
        } while (i != begin);
        assert( i == begin);
        assert( su == 15);

        // Assignment.
        i = end;
        su = 0;
        k  = 5;
        // Loop with post-decrement.
        do {
            Iterator j = i--;
            assert(  i !=  j);
            if ( j != end) {
                assert( (*i)->key == (*j)->key - 1);
            }
            assert( k == (*i)->key);
            su += (*i)->key;
            --k;
        } while (i != begin);
        assert( i == begin);
        assert( su == 15);
    }
    CGAL::Assert_iterator( begin);
    CGAL::Assert_iterator( end);

    PtrList l2 = l;
    const PtrList& l1 = l2;
    typedef PtrList::const_iterator ConstIterBase;
    typedef Iterator_project< ConstIterBase,Ident>  C_Iterator;
    C_Iterator c_begin(l1.begin());
    C_Iterator c_end(l1.end());
    Assert_bidirectional_category(c_begin);
    Assert_bidirectional_category(c_end);
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(c_begin);
        CGAL::Assert_circulator_or_iterator(c_end);
        CGAL::Assert_is_at_least_forward_category(c_begin);
        CGAL::Assert_is_at_least_forward_category(c_end);
        typedef std::iterator_traits< C_Iterator >::value_type      VT;
        typedef std::iterator_traits< C_Iterator >::difference_type DT;
        assert(1==test_value_type(static_cast< VT* >(nullptr)));
        assert(1==test_distance_type(static_cast< DT* >(nullptr)));

        // Default constructor.
        C_Iterator z ;
        CGAL::Assert_circulator_or_iterator(z);
        z = c_end ;
        // Copy constructor.
        C_Iterator i = c_begin;

        // Check general support for circulators and iterators.
        assert( CGAL::is_empty_range( z, z));
        assert( ! CGAL::is_empty_range( i, c_end));

        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, c_end)) {   // superfluous
            do {
                assert( k == (*i)->key);
                su += (*i)->key;
                ++k;
                C_Iterator j = ++i;
                assert( i ==  j);
                if ( i != c_end) {
                    assert( (*i)->key == (*j)->key);
                }
            } while (i != c_end);  // Inequality and equality checked.
        }
        assert( i == c_end);  // Equality checked.
        assert( su == 15);

        // Assignment.
        i = c_begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, c_end)) {   // superfluous
            do {
                assert( k == (*i)->key);
                su += (*i)->key;
                ++k;
                C_Iterator j = i++;
                assert(  i !=  j);
                if ( i != c_end) {
                    assert( (*i)->key == (*j)->key + 1);
                }
            } while (i != c_end);
        }
        assert( i == c_end);
        assert( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_bidirectional_category(c_begin);
        CGAL::Assert_is_at_least_bidirectional_category(c_end);
        // Loop backwards and pre-decrement.
        C_Iterator i = c_end;
        int su = 0;
        int k  = 5;
        do {
            C_Iterator j = --i;
            assert(  i ==  j);
            assert( (*i)->key == (*j)->key);
            assert( k == (*i)->key);
            su += (*i)->key;
            --k;
        } while (i != c_begin);
        assert( i == c_begin);
        assert( su == 15);

        // Assignment.
        i = c_end;
        su = 0;
        k  = 5;
        // Loop with post-decrement.
        do {
            C_Iterator j = i--;
            assert(  i !=  j);
            if ( j != c_end) {
                assert( (*i)->key == (*j)->key - 1);
            }
            assert( k == (*i)->key);
            su += (*i)->key;
            --k;
        } while (i != c_begin);
        assert( i == c_begin);
        assert( su == 15);
    }
    CGAL::Assert_iterator( c_begin);
    CGAL::Assert_iterator( c_end);

    while(! l.empty()){
      delete static_cast<item*>(l.back());
      l.pop_back();
    }

    }

    //
    //
    //============================================
    //
    //

    {
    typedef In_place_list<item,false> List;
    List l;
    l.push_back( *new item(1));
    l.push_back( *new item(2));
    l.push_back( *new item(3));
    l.push_back( *new item(4));
    l.push_back( *new item(5));
    assert( l.size() == 5);
    typedef Cast_function_object<item,item_adaptor> Ident;
    typedef List::iterator IterBase;
    typedef Iterator_project<IterBase,Ident> Iterator;
    Iterator begin(l.begin());
    Iterator end(l.end());
    Assert_bidirectional_category(begin);
    Assert_bidirectional_category(end);
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_circulator_or_iterator(end);
        CGAL::Assert_is_at_least_forward_category(begin);
        CGAL::Assert_is_at_least_forward_category(end);
        typedef std::iterator_traits< Iterator >::value_type      VT;
        typedef std::iterator_traits< Iterator >::difference_type DT;
        assert(1==test_value_type(static_cast< VT* >(nullptr)));
        assert(1==test_distance_type(static_cast< DT* >(nullptr)));

        // Default constructor.
        Iterator z ;
        CGAL::Assert_circulator_or_iterator(z);
        z = end ;
        // Copy constructor.
        Iterator i = begin;

        // Check general support for circulators and iterators.
        assert( CGAL::is_empty_range( z, z));
        assert( ! CGAL::is_empty_range( i, end));

        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, end)) {   // superfluous
            do {
                assert( k == (*i).get_key());
                su += (*i).get_key();
                ++k;
                Iterator j = ++i;
                assert( i ==  j);
                if ( i != end) {
                    assert( (*i).get_key() == (*j).get_key());
                }
            } while (i != end);  // Inequality and equality checked.
        }
        assert( i == end);  // Equality checked.
        assert( su == 15);

        // Assignment.
        i = begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, end)) {   // superfluous
            do {
                assert( k == (*i).get_key());
                su += (*i).get_key();
                ++k;
                Iterator j = i++;
                assert(  i !=  j);
                if ( i != end) {
                    assert( (*i).get_key() == (*j).get_key() + 1);
                }
            } while (i != end);
        }
        assert( i == end);
        assert( su == 15);
    }
    { // Open own scope to hide local variables.
        // Change three elements and check post-/pre-increment.
        Iterator i = begin;
        (*i++).get_key() = 4;
        assert( 4 == (*begin).get_key());
        assert( 2 == (*i).get_key());
        (*i++).get_key() = 3;
        assert( 3 == (*i).get_key());
        (*++i).get_key() = 7;
        assert( 7 == (*i).get_key());

        // Check the setting and reset these elements
        // to their original values.
        i = begin;
        assert( 4 == (*i).get_key());
        (*i).get_key() = 1;
        i++;
        assert( 3 == (*i).get_key());
        (*i++).get_key() = 2;
        assert( 3 == (*i).get_key());
        i++;
        assert( 7 == (*i).get_key());
        (*i).get_key() = 4;

        // Check the resetting.
        i = begin;
        int k = 1;
        do {
            assert( k == (*i).get_key());
            ++i;
            ++k;
        } while (i != end);
    }
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_circulator_or_iterator(end);
        CGAL::Assert_is_at_least_forward_category(begin);
        CGAL::Assert_is_at_least_forward_category(end);
        typedef std::iterator_traits< Iterator >::value_type      VT;
        typedef std::iterator_traits< Iterator >::difference_type DT;
        assert(1==test_value_type(static_cast< VT* >(nullptr)));
        assert(1==test_distance_type(static_cast< DT* >(nullptr)));

        // Default constructor.
        Iterator z ;
        CGAL::Assert_circulator_or_iterator(z);
        z = end ;
        // Copy constructor.
        Iterator i = begin;

        // Check general support for circulators and iterators.
        assert( CGAL::is_empty_range( z, z));
        assert( ! CGAL::is_empty_range( i, end));

        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, end)) {   // superfluous
            do {
                assert( k == (*i).get_key());
                su += (*i).get_key();
                ++k;
                Iterator j = ++i;
                assert( i ==  j);
                if ( i != end) {
                    assert( (*i).get_key() == (*j).get_key());
                }
            } while (i != end);  // Inequality and equality checked.
        }
        assert( i == end);  // Equality checked.
        assert( su == 15);

        // Assignment.
        i = begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, end)) {   // superfluous
            do {
                assert( k == (*i).get_key());
                su += (*i).get_key();
                ++k;
                Iterator j = i++;
                assert(  i !=  j);
                if ( i != end) {
                    assert( (*i).get_key() == (*j).get_key() + 1);
                }
            } while (i != end);
        }
        assert( i == end);
        assert( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_bidirectional_category(begin);
        CGAL::Assert_is_at_least_bidirectional_category(end);
        // Loop backwards and pre-decrement.
        Iterator i = end;
        int su = 0;
        int k  = 5;
        do {
            Iterator j = --i;
            assert(  i ==  j);
            assert( (*i).get_key() == (*j).get_key());
            assert( k == (*i).get_key());
            su += (*i).get_key();
            --k;
        } while (i != begin);
        assert( i == begin);
        assert( su == 15);

        // Assignment.
        i = end;
        su = 0;
        k  = 5;
        // Loop with post-decrement.
        do {
            Iterator j = i--;
            assert(  i !=  j);
            if ( j != end) {
                assert( (*i).get_key() == (*j).get_key() - 1);
            }
            assert( k == (*i).get_key());
            su += (*i).get_key();
            --k;
        } while (i != begin);
        assert( i == begin);
        assert( su == 15);
    }
    CGAL::Assert_iterator( begin);
    CGAL::Assert_iterator( end);

    List l2 = l;
    const List& l1 = l2;
    typedef List::const_iterator ConstIterBase;
    typedef Iterator_project< ConstIterBase,Ident,const item&,
      const item*, std::ptrdiff_t,
                              std::bidirectional_iterator_tag>
    C_Iterator;
    C_Iterator c_begin(l1.begin());
    C_Iterator c_end(l1.end());
    Assert_bidirectional_category(c_begin);
    Assert_bidirectional_category(c_end);
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(c_begin);
        CGAL::Assert_circulator_or_iterator(c_end);
        CGAL::Assert_is_at_least_forward_category(c_begin);
        CGAL::Assert_is_at_least_forward_category(c_end);
        typedef std::iterator_traits< C_Iterator >::value_type      VT;
        typedef std::iterator_traits< C_Iterator >::difference_type DT;
        assert(1==test_value_type(static_cast< VT* >(nullptr)));
        assert(1==test_distance_type(static_cast< DT* >(nullptr)));

        // Default constructor.
        C_Iterator z ;
        CGAL::Assert_circulator_or_iterator(z);
        z = c_end ;
        // Copy constructor.
        C_Iterator i = c_begin;

        // Check general support for circulators and iterators.
        assert( CGAL::is_empty_range( z, z));
        assert( ! CGAL::is_empty_range( i, c_end));

        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, c_end)) {   // superfluous
            do {
                assert( k == (*i).get_key());
                su += (*i).get_key();
                ++k;
                C_Iterator j = ++i;
                assert( i ==  j);
                if ( i != c_end) {
                    assert( (*i).get_key() == (*j).get_key());
                }
            } while (i != c_end);  // Inequality and equality checked.
        }
        assert( i == c_end);  // Equality checked.
        assert( su == 15);

        // Assignment.
        i = c_begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, c_end)) {   // superfluous
            do {
                assert( k == (*i).get_key());
                su += (*i).get_key();
                ++k;
                C_Iterator j = i++;
                assert(  i !=  j);
                if ( i != c_end) {
                    assert( (*i).get_key() == (*j).get_key() + 1);
                }
            } while (i != c_end);
        }
        assert( i == c_end);
        assert( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_bidirectional_category(c_begin);
        CGAL::Assert_is_at_least_bidirectional_category(c_end);
        // Loop backwards and pre-decrement.
        C_Iterator i = c_end;
        int su = 0;
        int k  = 5;
        do {
            C_Iterator j = --i;
            assert(  i ==  j);
            assert( (*i).get_key() == (*j).get_key());
            assert( k == (*i).get_key());
            su += (*i).get_key();
            --k;
        } while (i != c_begin);
        assert( i == c_begin);
        assert( su == 15);

        // Assignment.
        i = c_end;
        su = 0;
        k  = 5;
        // Loop with post-decrement.
        do {
            C_Iterator j = i--;
            assert(  i !=  j);
            if ( j != c_end) {
                assert( (*i).get_key() == (*j).get_key() - 1);
            }
            assert( k == (*i).get_key());
            su += (*i).get_key();
            --k;
        } while (i != c_begin);
        assert( i == c_begin);
        assert( su == 15);
    }
    CGAL::Assert_iterator( c_begin);
    CGAL::Assert_iterator( c_end);
    l.destroy();
    l2.destroy();
    }
}

void test_Iterator_transform()
{
  typedef std::vector<Pair> Vector;
  Vector v;
  v.push_back( Pair(0,1) );
  v.push_back( Pair(1,2) );
  v.push_back( Pair(2,3) );
  v.push_back( Pair(3,4) );
  v.push_back( Pair(4,5) );
  assert( v.size() == 5);

  typedef Vector::iterator IterBase;
  typedef Iterator_transform<IterBase,PairToPoint> Iterator;
  Iterator begin(v.begin());
  Iterator end  (v.end  ());
  Iterator last (v.end  () - 1 );

  assert( begin < end   ) ;
  assert( last  < end   ) ;
  assert( end   > begin ) ;
  assert( last  > begin ) ;

  assert( (last + 1 == end  ) ) ;
  assert( (end - 1  == last ) ) ;

  Assert_random_access_category(begin);
  CGAL::Assert_circulator_or_iterator(begin);
  CGAL::Assert_is_at_least_forward_category(begin);

  typedef std::iterator_traits< Iterator >::value_type      VT;
  typedef std::iterator_traits< Iterator >::difference_type DT;
  assert(1==test_value_type(static_cast< VT* >(nullptr)));
  assert(1==test_distance_type(static_cast< DT* >(nullptr)));

  assert(   CGAL::is_empty_range( end  , end) );
  assert( ! CGAL::is_empty_range( begin, end) );

  // Default constructor.
  Iterator z ;

  // Assignment
  z = end ;

  assert( z == end ) ;

  // Copy constructor.
  Iterator f1 = begin;
  Iterator f2 = begin;
  Iterator b1 = last ;
  Iterator b2 = last ;

  assert( f1 == begin ) ;
  assert( b1 == last  ) ;

  Point fp(0,1);
  Point bp(4,5);

  // Check pre/post-increment and decrement
  // Check dereference.
  // Check equality/inequality
  bool done = false ;

  do
  {
    assert( fp == (*f1) );
    assert( fp == (*f2) );
    assert( bp == (*b1) );
    assert( bp == (*b2) );

    if ( f1 < last )
    {
      Iterator of1 = f1 ;
      Iterator of2 = f2 ;
      Iterator ob1 = b1 ;
      Iterator ob2 = b2 ;

      Iterator fn1 = ++f1;
      assert( f1 != of1 );
      assert( fn1 == f1 );

      Iterator fn2 = f2++ ;
      assert( f2 != of2 );
      assert( fn2 == of2 );

      Iterator bn1 = --b1;
      assert( b1 != ob1 );
      assert( bn1 == b1 );

      Iterator bn2 = b2-- ;
      assert( b2 != ob2 );
      assert( bn2 == ob2 );

      Iterator next = of1 + 1 ;
      assert( f1 == next );

      assert (  f1 >   of1 ) ;
      assert (  f1 >=  of1 ) ;
      assert (  of1 <  f1  ) ;
      assert (  of1 <= f1  ) ;
      assert ( !(f1 < of1) ) ;
      assert ( !(of1 > f1) ) ;

      assert( of1[0] == *of1 ) ;
      assert( of1[1] == *f1  ) ;

      of1 += 1 ;
      assert( f1 == of1 );

      Iterator prev = ob1 - 1 ;
      assert( b1 == prev );

      ob1 -= 1 ;
      assert( b1 == ob1 );

      fp.x ++ ; fp.y ++ ;
      bp.x -- ; bp.y -- ;
    }

    bool done1 = ( f1 == last ) ;
    bool done2 = ( b1 == begin ) ;

    assert ( done1 == done2 ) ;

    done = done1 || done2 ;
  }
  while ( !done );

}

void test_Circulator_project() {
  {
    typedef In_place_list<item,false> List;
    List l;
    l.push_back( *new item(1));
    l.push_back( *new item(2));
    l.push_back( *new item(3));
    l.push_back( *new item(4));
    l.push_back( *new item(5));
    assert( l.size() == 5);
    typedef List::iterator IterBase;
    typedef Identity<item> Ident;
    typedef Bidirectional_circulator_from_iterator<IterBase,item,
      std::size_t,std::ptrdiff_t> CircBase;
    typedef Circulator_project<CircBase,Ident,item&,item*>
      Circulator;
    Circulator begin(CircBase( l.begin(),l.end()));
    Assert_bidirectional_category(begin);
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_is_at_least_forward_category(begin);
        CGAL::Assert_is_at_least_forward_category(begin);
        typedef std::iterator_traits< Circulator >::value_type      VT;
        typedef std::iterator_traits< Circulator >::difference_type DT;
        assert(1==test_value_type(static_cast< VT* >(nullptr)));
        assert(1==test_distance_type(static_cast< DT* >(nullptr)));

        // Default constructor.
        Circulator z = Circulator();
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Circulator i = begin;

        // Check general support for circulators and iterators.
        assert( CGAL::is_empty_range( z, z));
        assert( ! CGAL::is_empty_range( i, begin));

        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                assert( k == (*i).key);
                su += (*i).key;
                ++k;
                Circulator j = ++i;
                assert( i ==  j);
                if ( i != begin) {
                    assert( (*i).key == (*j).key);
                }
            } while (i != begin);  // Inequality and equality checked.
        }
        assert( i == begin);  // Equality checked.
        assert( su == 15);

        // Assignment.
        i = begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                assert( k == (*i).key);
                su += (*i).key;
                ++k;
                Circulator j = i++;
                assert(  i !=  j);
                if ( i != begin) {
                    assert( (*i).key == (*j).key + 1);
                }
            } while (i != begin);
        }
        assert( i == begin);
        assert( su == 15);
    }
    { // Open own scope to hide local variables.
        // Change three elements and check post-/pre-increment.
        Circulator i = begin;
        (*i++).key = 4;
        assert( 4 == (*begin).key);
        assert( 2 == (*i).key);
        (*i++).key = 3;
        assert( 3 == (*i).key);
        (*++i).key = 7;
        assert( 7 == (*i).key);

        // Check the setting and reset these elements
        // to their original values.
        i = begin;
        assert( 4 == (*i).key);
        (*i).key = 1;
        i++;
        assert( 3 == (*i).key);
        (*i++).key = 2;
        assert( 3 == (*i).key);
        i++;
        assert( 7 == (*i).key);
        (*i).key = 4;

        // Check the resetting.
        i = begin;
        int k = 1;
        do {
            assert( k == (*i).key);
            ++i;
            ++k;
        } while (i != begin);
    }
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_is_at_least_forward_category(begin);
        CGAL::Assert_is_at_least_forward_category(begin);
        typedef std::iterator_traits< Circulator >::value_type      VT;
        typedef std::iterator_traits< Circulator >::difference_type DT;
        assert(1==test_value_type(static_cast< VT* >(nullptr)));
        assert(1==test_distance_type(static_cast< DT* >(nullptr)));

        // Default constructor.
        Circulator z = Circulator();
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Circulator i = begin;

        // Check general support for circulators and iterators.
        assert( CGAL::is_empty_range( z, z));
        assert( ! CGAL::is_empty_range( i, begin));

        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                assert( k == (*i).key);
                su += (*i).key;
                ++k;
                Circulator j = ++i;
                assert( i ==  j);
                if ( i != begin) {
                    assert( (*i).key == (*j).key);
                }
            } while (i != begin);  // Inequality and equality checked.
        }
        assert( i == begin);  // Equality checked.
        assert( su == 15);

        // Assignment.
        i = begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                assert( k == (*i).key);
                su += (*i).key;
                ++k;
                Circulator j = i++;
                assert(  i !=  j);
                if ( i != begin) {
                    assert( (*i).key == (*j).key + 1);
                }
            } while (i != begin);
        }
        assert( i == begin);
        assert( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_bidirectional_category(begin);
        CGAL::Assert_is_at_least_bidirectional_category(begin);
        // Loop backwards and pre-decrement.
        Circulator i = begin;
        int su = 0;
        int k  = 5;
        do {
            Circulator j = --i;
            assert(  i ==  j);
            assert( (*i).key == (*j).key);
            assert( k == (*i).key);
            su += (*i).key;
            --k;
        } while (i != begin);
        assert( i == begin);
        assert( su == 15);

        // Assignment.
        i = begin;
        su = 0;
        k  = 5;
        // Loop with post-decrement.
        do {
            Circulator j = i--;
            assert(  i !=  j);
            if ( j != begin) {
                assert( (*i).key == (*j).key - 1);
            }
            assert( k == (*i).key);
            su += (*i).key;
            --k;
        } while (i != begin);
        assert( i == begin);
        assert( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_circulator( begin);
        CGAL::Assert_circulator( begin);

        // Check the local type parameters.
        Circulator::value_type      k1;
        k1.key = 1;
        Circulator::reference       k2 = k1;
        (void)k2;
        assert( k2.key == 1);
        Circulator::pointer         p1 = &k1;
        (void)p1;
        assert( (*p1).key == 1);
        k1.key = 3;
        assert( k1.key == 3);
        assert( k2.key == 3);
        assert( (*p1).key == 3);
        k1.key = 6;
        assert( k1.key == 6);
        assert( k2.key == 6);
        assert( (*p1).key == 6);
        Circulator::size_type s = 5;
        assert( s == 5);
        Circulator::difference_type d = -5;
        assert( d == -5);

        // Check tests for empty data structures.
        Circulator z = Circulator();
        assert(   z == nullptr);
        assert( ! (z != nullptr));
        Circulator i = begin;
        assert( ! (i == nullptr));
        assert(   i != nullptr);
        assert( i == begin);
        assert( i == begin);
        // Do I reach myself.
        ++i;
        Circulator j = i;
        int k = 0;
        do {
            assert( k < 5);
            ++k;
            ++i;
        } while( i != j);
        assert( k == 5);
    }
    { // Open own scope to hide local variables.
        // Do I reach myself backwards.
        Circulator i = begin;
        ++i;
        Circulator j = i;
        int k = 0;
        do {
            assert( k < 5);
            ++k;
            --i;
        } while( i != j);
        assert( k == 5);
    }

    List l2 = l;
    const List& l1 = l2;
    typedef List::const_iterator ConstIterBase;
    typedef Bidirectional_const_circulator_from_iterator<
      ConstIterBase, item, std::size_t, std::ptrdiff_t>
    ConstCircBase;
typedef Circulator_project<ConstCircBase,Ident,const item&,
  const item*> C_Circulator;
C_Circulator c_begin(ConstCircBase(l1.begin(),l1.end()));
Assert_bidirectional_category(c_begin);
{ // Open own scope to hide local variables.
    // Check generally correct parameter properties.
    CGAL::Assert_circulator_or_iterator(c_begin);
    CGAL::Assert_circulator_or_iterator(c_begin);
    CGAL::Assert_is_at_least_forward_category(c_begin);
    CGAL::Assert_is_at_least_forward_category(c_begin);
    typedef std::iterator_traits< C_Circulator >::value_type      VT;
    typedef std::iterator_traits< C_Circulator >::difference_type DT;
    assert(1==test_value_type(static_cast< VT* >(nullptr)));
    assert(1==test_distance_type(static_cast< DT* >(nullptr)));

    // Default constructor.
    C_Circulator z = C_Circulator();
    CGAL::Assert_circulator_or_iterator(z);
    // Copy constructor.
    C_Circulator i = c_begin;

    // Check general support for circulators and iterators.
    assert( CGAL::is_empty_range( z, z));
    assert( ! CGAL::is_empty_range( i, c_begin));

    int su = 0;
    int k  = 1;
    // Check general loop, pre-increment, dereference.
    if (! CGAL::is_empty_range( i, c_begin)) {   // superfluous
        do {
            assert( k == (*i).key);
            su += (*i).key;
            ++k;
            C_Circulator j = ++i;
            assert( i ==  j);
            if ( i != c_begin) {
                assert( (*i).key == (*j).key);
            }
        } while (i != c_begin);  // Inequality and equality checked.
    }
    assert( i == c_begin);  // Equality checked.
    assert( su == 15);

    // Assignment.
    i = c_begin;
    su = 0;
    k  = 1;
    // Loop with post increment.
    if (! CGAL::is_empty_range( i, c_begin)) {   // superfluous
        do {
            assert( k == (*i).key);
            su += (*i).key;
            ++k;
            C_Circulator j = i++;
            assert(  i !=  j);
            if ( i != c_begin) {
                assert( (*i).key == (*j).key + 1);
            }
        } while (i != c_begin);
    }
    assert( i == c_begin);
    assert( su == 15);
}
{ // Open own scope to hide local variables.
    CGAL::Assert_is_at_least_bidirectional_category(c_begin);
    CGAL::Assert_is_at_least_bidirectional_category(c_begin);
    // Loop backwards and pre-decrement.
    C_Circulator i = c_begin;
    int su = 0;
    int k  = 5;
    do {
        C_Circulator j = --i;
        assert(  i ==  j);
        assert( (*i).key == (*j).key);
        assert( k == (*i).key);
        su += (*i).key;
        --k;
    } while (i != c_begin);
    assert( i == c_begin);
    assert( su == 15);

    // Assignment.
    i = c_begin;
    su = 0;
    k  = 5;
    // Loop with post-decrement.
    do {
        C_Circulator j = i--;
        assert(  i !=  j);
        if ( j != c_begin) {
            assert( (*i).key == (*j).key - 1);
        }
        assert( k == (*i).key);
        su += (*i).key;
        --k;
    } while (i != c_begin);
    assert( i == c_begin);
    assert( su == 15);
}
{ // Open own scope to hide local variables.
    CGAL::Assert_circulator( c_begin);
    CGAL::Assert_circulator( c_begin);

    // Check the local type parameters.
    C_Circulator::value_type      k1;
    k1.key = 1;
    C_Circulator::reference       k2 = k1;
    (void)k2;
    assert( k2.key == 1);
    C_Circulator::pointer         p1 = &k1;
    (void)p1;
    assert( (*p1).key == 1);
    k1.key = 3;
    assert( k1.key == 3);
    assert( k2.key == 3);
    assert( (*p1).key == 3);
    k1.key = 6;
    assert( k1.key == 6);
    assert( k2.key == 6);
    assert( (*p1).key == 6);
    C_Circulator::size_type s = 5;
    assert( s == 5);
    C_Circulator::difference_type d = -5;
    assert( d == -5);

    // Check tests for empty data structures.
    C_Circulator z = C_Circulator();
    assert(   z == nullptr);
    assert( ! (z != nullptr));
    C_Circulator i = c_begin;
    assert( ! (i == nullptr));
    assert(   i != nullptr);
    assert( i == c_begin);
    assert( i == c_begin);
    // Do I reach myself.
    ++i;
    C_Circulator j = i;
    int k = 0;
    do {
        assert( k < 5);
        ++k;
        ++i;
    } while( i != j);
    assert( k == 5);
}
{ // Open own scope to hide local variables.
    // Do I reach myself backwards.
    C_Circulator i = c_begin;
    ++i;
    C_Circulator j = i;
    int k = 0;
    do {
        assert( k < 5);
        ++k;
        --i;
    } while( i != j);
    assert( k == 5);
}
l.destroy();
l2.destroy();
  }
  {
    typedef In_place_list<item,false> List;
    List l;
    l.push_back( *new item(1));
    l.push_back( *new item(2));
    l.push_back( *new item(3));
    l.push_back( *new item(4));
    l.push_back( *new item(5));
    assert( l.size() == 5);
    typedef List::iterator IterBase;
    typedef Cast_function_object<item,item> Ident;
    typedef Bidirectional_circulator_from_iterator<IterBase,item,
      std::size_t,std::ptrdiff_t> CircBase;
    typedef Circulator_project<CircBase,Ident,item&,item*>
      Circulator;
    Circulator begin(CircBase( l.begin(),l.end()));
    Assert_bidirectional_category(begin);
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_is_at_least_forward_category(begin);
        CGAL::Assert_is_at_least_forward_category(begin);
        typedef std::iterator_traits< Circulator >::value_type      VT;
        typedef std::iterator_traits< Circulator >::difference_type DT;
        assert(1==test_value_type(static_cast< VT* >(nullptr)));
        assert(1==test_distance_type(static_cast< DT* >(nullptr)));

        // Default constructor.
        Circulator z = Circulator();
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Circulator i = begin;

        // Check general support for circulators and iterators.
        assert( CGAL::is_empty_range( z, z));
        assert( ! CGAL::is_empty_range( i, begin));

        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                assert( k == (*i).key);
                su += (*i).key;
                ++k;
                Circulator j = ++i;
                assert( i ==  j);
                if ( i != begin) {
                    assert( (*i).key == (*j).key);
                }
            } while (i != begin);  // Inequality and equality checked.
        }
        assert( i == begin);  // Equality checked.
        assert( su == 15);

        // Assignment.
        i = begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                assert( k == (*i).key);
                su += (*i).key;
                ++k;
                Circulator j = i++;
                assert(  i !=  j);
                if ( i != begin) {
                    assert( (*i).key == (*j).key + 1);
                }
            } while (i != begin);
        }
        assert( i == begin);
        assert( su == 15);
    }
    { // Open own scope to hide local variables.
        // Change three elements and check post-/pre-increment.
        Circulator i = begin;
        (*i++).key = 4;
        assert( 4 == (*begin).key);
        assert( 2 == (*i).key);
        (*i++).key = 3;
        assert( 3 == (*i).key);
        (*++i).key = 7;
        assert( 7 == (*i).key);

        // Check the setting and reset these elements
        // to their original values.
        i = begin;
        assert( 4 == (*i).key);
        (*i).key = 1;
        i++;
        assert( 3 == (*i).key);
        (*i++).key = 2;
        assert( 3 == (*i).key);
        i++;
        assert( 7 == (*i).key);
        (*i).key = 4;

        // Check the resetting.
        i = begin;
        int k = 1;
        do {
            assert( k == (*i).key);
            ++i;
            ++k;
        } while (i != begin);
    }
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_is_at_least_forward_category(begin);
        CGAL::Assert_is_at_least_forward_category(begin);
        typedef std::iterator_traits< Circulator >::value_type      VT;
        typedef std::iterator_traits< Circulator >::difference_type DT;
        assert(1==test_value_type(static_cast< VT* >(nullptr)));
        assert(1==test_distance_type(static_cast< DT* >(nullptr)));

        // Default constructor.
        Circulator z = Circulator();
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Circulator i = begin;

        // Check general support for circulators and iterators.
        assert( CGAL::is_empty_range( z, z));
        assert( ! CGAL::is_empty_range( i, begin));

        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                assert( k == (*i).key);
                su += (*i).key;
                ++k;
                Circulator j = ++i;
                assert( i ==  j);
                if ( i != begin) {
                    assert( (*i).key == (*j).key);
                }
            } while (i != begin);  // Inequality and equality checked.
        }
        assert( i == begin);  // Equality checked.
        assert( su == 15);

        // Assignment.
        i = begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                assert( k == (*i).key);
                su += (*i).key;
                ++k;
                Circulator j = i++;
                assert(  i !=  j);
                if ( i != begin) {
                    assert( (*i).key == (*j).key + 1);
                }
            } while (i != begin);
        }
        assert( i == begin);
        assert( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_bidirectional_category(begin);
        CGAL::Assert_is_at_least_bidirectional_category(begin);
        // Loop backwards and pre-decrement.
        Circulator i = begin;
        int su = 0;
        int k  = 5;
        do {
            Circulator j = --i;
            assert(  i ==  j);
            assert( (*i).key == (*j).key);
            assert( k == (*i).key);
            su += (*i).key;
            --k;
        } while (i != begin);
        assert( i == begin);
        assert( su == 15);

        // Assignment.
        i = begin;
        su = 0;
        k  = 5;
        // Loop with post-decrement.
        do {
            Circulator j = i--;
            assert(  i !=  j);
            if ( j != begin) {
                assert( (*i).key == (*j).key - 1);
            }
            assert( k == (*i).key);
            su += (*i).key;
            --k;
        } while (i != begin);
        assert( i == begin);
        assert( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_circulator( begin);
        CGAL::Assert_circulator( begin);

        // Check the local type parameters.
        Circulator::value_type      k1;
        k1.key = 1;
        Circulator::reference       k2 = k1;
        (void)k2;
        assert( k2.key == 1);
        Circulator::pointer         p1 = &k1;
        (void)p1;
        assert( (*p1).key == 1);
        k1.key = 3;
        assert( k1.key == 3);
        assert( k2.key == 3);
        assert( (*p1).key == 3);
        k1.key = 6;
        assert( k1.key == 6);
        assert( k2.key == 6);
        assert( (*p1).key == 6);
        Circulator::size_type s = 5;
        assert( s == 5);
        Circulator::difference_type d = -5;
        assert( d == -5);

        // Check tests for empty data structures.
        Circulator z = Circulator();
        assert(   z == nullptr);
        assert( ! (z != nullptr));
        Circulator i = begin;
        assert( ! (i == nullptr));
        assert(   i != nullptr);
        assert( i == begin);
        assert( i == begin);
        // Do I reach myself.
        ++i;
        Circulator j = i;
        int k = 0;
        do {
            assert( k < 5);
            ++k;
            ++i;
        } while( i != j);
        assert( k == 5);
    }
    { // Open own scope to hide local variables.
        // Do I reach myself backwards.
        Circulator i = begin;
        ++i;
        Circulator j = i;
        int k = 0;
        do {
            assert( k < 5);
            ++k;
            --i;
        } while( i != j);
        assert( k == 5);
    }

    List l2 = l;
    const List& l1 = l2;
    typedef List::const_iterator ConstIterBase;
    typedef Bidirectional_const_circulator_from_iterator<
      ConstIterBase, item, std::size_t, std::ptrdiff_t>
    ConstCircBase;
typedef Circulator_project<ConstCircBase,Ident,const item&,
  const item*> C_Circulator;
C_Circulator c_begin(ConstCircBase(l1.begin(),l1.end()));
Assert_bidirectional_category(c_begin);
{ // Open own scope to hide local variables.
    // Check generally correct parameter properties.
    CGAL::Assert_circulator_or_iterator(c_begin);
    CGAL::Assert_circulator_or_iterator(c_begin);
    CGAL::Assert_is_at_least_forward_category(c_begin);
    CGAL::Assert_is_at_least_forward_category(c_begin);
    typedef std::iterator_traits< C_Circulator >::value_type      VT;
    typedef std::iterator_traits< C_Circulator >::difference_type DT;
    assert(1==test_value_type(static_cast< VT* >(nullptr)));
    assert(1==test_distance_type(static_cast< DT* >(nullptr)));

    // Default constructor.
    C_Circulator z = C_Circulator();
    CGAL::Assert_circulator_or_iterator(z);
    // Copy constructor.
    C_Circulator i = c_begin;

    // Check general support for circulators and iterators.
    assert( CGAL::is_empty_range( z, z));
    assert( ! CGAL::is_empty_range( i, c_begin));

    int su = 0;
    int k  = 1;
    // Check general loop, pre-increment, dereference.
    if (! CGAL::is_empty_range( i, c_begin)) {   // superfluous
        do {
            assert( k == (*i).key);
            su += (*i).key;
            ++k;
            C_Circulator j = ++i;
            assert( i ==  j);
            if ( i != c_begin) {
                assert( (*i).key == (*j).key);
            }
        } while (i != c_begin);  // Inequality and equality checked.
    }
    assert( i == c_begin);  // Equality checked.
    assert( su == 15);

    // Assignment.
    i = c_begin;
    su = 0;
    k  = 1;
    // Loop with post increment.
    if (! CGAL::is_empty_range( i, c_begin)) {   // superfluous
        do {
            assert( k == (*i).key);
            su += (*i).key;
            ++k;
            C_Circulator j = i++;
            assert(  i !=  j);
            if ( i != c_begin) {
                assert( (*i).key == (*j).key + 1);
            }
        } while (i != c_begin);
    }
    assert( i == c_begin);
    assert( su == 15);
}
{ // Open own scope to hide local variables.
    CGAL::Assert_is_at_least_bidirectional_category(c_begin);
    CGAL::Assert_is_at_least_bidirectional_category(c_begin);
    // Loop backwards and pre-decrement.
    C_Circulator i = c_begin;
    int su = 0;
    int k  = 5;
    do {
        C_Circulator j = --i;
        assert(  i ==  j);
        assert( (*i).key == (*j).key);
        assert( k == (*i).key);
        su += (*i).key;
        --k;
    } while (i != c_begin);
    assert( i == c_begin);
    assert( su == 15);

    // Assignment.
    i = c_begin;
    su = 0;
    k  = 5;
    // Loop with post-decrement.
    do {
        C_Circulator j = i--;
        assert(  i !=  j);
        if ( j != c_begin) {
            assert( (*i).key == (*j).key - 1);
        }
        assert( k == (*i).key);
        su += (*i).key;
        --k;
    } while (i != c_begin);
    assert( i == c_begin);
    assert( su == 15);
}
{ // Open own scope to hide local variables.
    CGAL::Assert_circulator( c_begin);
    CGAL::Assert_circulator( c_begin);

    // Check the local type parameters.
    C_Circulator::value_type      k1;
    k1.key = 1;
    C_Circulator::reference       k2 = k1;
    (void)k2;
    assert( k2.key == 1);
    C_Circulator::pointer         p1 = &k1;
    (void)p1;
    assert( (*p1).key == 1);
    k1.key = 3;
    assert( k1.key == 3);
    assert( k2.key == 3);
    assert( (*p1).key == 3);
    k1.key = 6;
    assert( k1.key == 6);
    assert( k2.key == 6);
    assert( (*p1).key == 6);
    C_Circulator::size_type s = 5;
    assert( s == 5);
    C_Circulator::difference_type d = -5;
    assert( d == -5);

    // Check tests for empty data structures.
    C_Circulator z = C_Circulator();
    assert(   z == nullptr);
    assert( ! (z != nullptr));
    C_Circulator i = c_begin;
    assert( ! (i == nullptr));
    assert(   i != nullptr);
    assert( i == c_begin);
    assert( i == c_begin);
    // Do I reach myself.
    ++i;
    C_Circulator j = i;
    int k = 0;
    do {
        assert( k < 5);
        ++k;
        ++i;
    } while( i != j);
    assert( k == 5);
}
{ // Open own scope to hide local variables.
    // Do I reach myself backwards.
    C_Circulator i = c_begin;
    ++i;
    C_Circulator j = i;
    int k = 0;
    do {
        assert( k < 5);
        ++k;
        --i;
    } while( i != j);
    assert( k == 5);
}
l.destroy();
l2.destroy();
  }
  {
    typedef std::vector<int> Vector;
    Vector v;
    v.push_back(1);
    v.push_back(2);
    v.push_back(3);
    v.push_back(4);
    v.push_back(5);
    assert( v.size() == 5);
    typedef Vector::iterator IterBase;
    typedef Identity<int> Ident;
    typedef Random_access_circulator_from_iterator<IterBase,int,
      std::size_t,std::ptrdiff_t> CircBase;
    typedef Circulator_project<CircBase,Ident,int&,int*>
      Circulator;
    Circulator begin(CircBase( v.begin(),v.end()));
    Assert_random_access_category(begin);
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_is_at_least_forward_category(begin);
        CGAL::Assert_is_at_least_forward_category(begin);
        typedef std::iterator_traits< Circulator >::value_type      VT;
        typedef std::iterator_traits< Circulator >::difference_type DT;
        assert(1==test_value_type(static_cast< VT* >(nullptr)));
        assert(1==test_distance_type(static_cast< DT* >(nullptr)));

        // Default constructor.
        Circulator z = Circulator();
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Circulator i = begin;

        // Check general support for circulators and iterators.
        assert( CGAL::is_empty_range( z, z));
        assert( ! CGAL::is_empty_range( i, begin));

        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                assert( k == (*i));
                su += (*i);
                ++k;
                Circulator j = ++i;
                assert( i ==  j);
                if ( i != begin) {
                    assert( (*i) == (*j));
                }
            } while (i != begin);  // Inequality and equality checked.
        }
        assert( i == begin);  // Equality checked.
        assert( su == 15);

        // Assignment.
        i = begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                assert( k == (*i));
                su += (*i);
                ++k;
                Circulator j = i++;
                assert(  i !=  j);
                if ( i != begin) {
                    assert( (*i) == (*j) + 1);
                }
            } while (i != begin);
        }
        assert( i == begin);
        assert( su == 15);
    }
    { // Open own scope to hide local variables.
        // Change three elements and check post-/pre-increment.
        Circulator i = begin;
        (*i++) = 4;
        assert( 4 == (*begin));
        assert( 2 == (*i));
        (*i++) = 3;
        assert( 3 == (*i));
        (*++i) = 7;
        assert( 7 == (*i));

        // Check the setting and reset these elements
        // to their original values.
        i = begin;
        assert( 4 == (*i));
        (*i) = 1;
        i++;
        assert( 3 == (*i));
        (*i++) = 2;
        assert( 3 == (*i));
        i++;
        assert( 7 == (*i));
        (*i) = 4;

        // Check the resetting.
        i = begin;
        int k = 1;
        do {
            assert( k == (*i));
            ++i;
            ++k;
        } while (i != begin);
    }
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_is_at_least_forward_category(begin);
        CGAL::Assert_is_at_least_forward_category(begin);
        typedef std::iterator_traits< Circulator >::value_type      VT;
        typedef std::iterator_traits< Circulator >::difference_type DT;
        assert(1==test_value_type(static_cast< VT* >(nullptr)));
        assert(1==test_distance_type(static_cast< DT* >(nullptr)));

        // Default constructor.
        Circulator z = Circulator();
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Circulator i = begin;

        // Check general support for circulators and iterators.
        assert( CGAL::is_empty_range( z, z));
        assert( ! CGAL::is_empty_range( i, begin));

        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                assert( k == (*i));
                su += (*i);
                ++k;
                Circulator j = ++i;
                assert( i ==  j);
                if ( i != begin) {
                    assert( (*i) == (*j));
                }
            } while (i != begin);  // Inequality and equality checked.
        }
        assert( i == begin);  // Equality checked.
        assert( su == 15);

        // Assignment.
        i = begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                assert( k == (*i));
                su += (*i);
                ++k;
                Circulator j = i++;
                assert(  i !=  j);
                if ( i != begin) {
                    assert( (*i) == (*j) + 1);
                }
            } while (i != begin);
        }
        assert( i == begin);
        assert( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_bidirectional_category(begin);
        CGAL::Assert_is_at_least_bidirectional_category(begin);
        // Loop backwards and pre-decrement.
        Circulator i = begin;
        int su = 0;
        int k  = 5;
        do {
            Circulator j = --i;
            assert(  i ==  j);
            assert( (*i) == (*j));
            assert( k == (*i));
            su += (*i);
            --k;
        } while (i != begin);
        assert( i == begin);
        assert( su == 15);

        // Assignment.
        i = begin;
        su = 0;
        k  = 5;
        // Loop with post-decrement.
        do {
            Circulator j = i--;
            assert(  i !=  j);
            if ( j != begin) {
                assert( (*i) == (*j) - 1);
            }
            assert( k == (*i));
            su += (*i);
            --k;
        } while (i != begin);
        assert( i == begin);
        assert( su == 15);
    }
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_is_at_least_forward_category(begin);
        CGAL::Assert_is_at_least_forward_category(begin);
        typedef std::iterator_traits< Circulator >::value_type      VT;
        typedef std::iterator_traits< Circulator >::difference_type DT;
        assert(1==test_value_type(static_cast< VT* >(nullptr)));
        assert(1==test_distance_type(static_cast< DT* >(nullptr)));

        // Default constructor.
        Circulator z = Circulator();
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Circulator i = begin;

        // Check general support for circulators and iterators.
        assert( CGAL::is_empty_range( z, z));
        assert( ! CGAL::is_empty_range( i, begin));

        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                assert( k == (*i));
                su += (*i);
                ++k;
                Circulator j = ++i;
                assert( i ==  j);
                if ( i != begin) {
                    assert( (*i) == (*j));
                }
            } while (i != begin);  // Inequality and equality checked.
        }
        assert( i == begin);  // Equality checked.
        assert( su == 15);

        // Assignment.
        i = begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                assert( k == (*i));
                su += (*i);
                ++k;
                Circulator j = i++;
                assert(  i !=  j);
                if ( i != begin) {
                    assert( (*i) == (*j) + 1);
                }
            } while (i != begin);
        }
        assert( i == begin);
        assert( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_bidirectional_category(begin);
        CGAL::Assert_is_at_least_bidirectional_category(begin);
        // Loop backwards and pre-decrement.
        Circulator i = begin;
        int su = 0;
        int k  = 5;
        do {
            Circulator j = --i;
            assert(  i ==  j);
            assert( (*i) == (*j));
            assert( k == (*i));
            su += (*i);
            --k;
        } while (i != begin);
        assert( i == begin);
        assert( su == 15);

        // Assignment.
        i = begin;
        su = 0;
        k  = 5;
        // Loop with post-decrement.
        do {
            Circulator j = i--;
            assert(  i !=  j);
            if ( j != begin) {
                assert( (*i) == (*j) - 1);
            }
            assert( k == (*i));
            su += (*i);
            --k;
        } while (i != begin);
        assert( i == begin);
        assert( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_random_access_category(begin);
        CGAL::Assert_is_at_least_random_access_category(begin);
        // Random access.
        int k;
        for( k = 0; k < 5; k++) {
            assert( 1+k == begin[k]);
        }

          int su = begin[0]
                 + begin[1]
                 + begin[2]
                 + begin[3]
                 + begin[4];
        assert( su == 15);

        // Jump around.
        Circulator i = begin;
        i += 3;
        assert( 4 == (*i));
        i -= 2;
        assert( 2 == (*i));
        i += 3;
        assert( 5 == (*i));
        i -= 4;
        assert( 1 == (*i));
        assert( i == begin);
        Circulator j = i + 3;
        assert( 4 == (*j));
        Circulator jj = j - 2;
        assert( 2 == (*jj));
        typedef std::ptrdiff_t PT;
        jj = PT(4) + jj;
        assert( jj == begin);
        Circulator ij = jj - 5;
        ij = jj - 5; // avoids warning with NDEBUG
        assert( ij == begin);

        // Difference test.
        assert( jj - i == 5  ||  jj - i == 0);
        assert( i + (j-i) == j);
        assert( (j-i) + i == j);
    }
    { // Open own scope to hide local variables.
        Circulator i = begin;
        i[2] = 18;
        i[4] = 9;
        i[3] = 12;
        assert( i[2] == 18);
        assert( i[4] == 9);
        assert( i[3] == 12);
        i[2] = 3;
        i[3] = 4;
        i[4] = 5;
        // Check the resetting.
        i = begin;
        int k = 1;
        do {
            assert( k == (*i));
            ++i;
            ++k;
        } while (i != begin);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_circulator( begin);
        CGAL::Assert_circulator( begin);

        // Check the local type parameters.
        Circulator::value_type      k1;
        k1 = 1;
        Circulator::reference       k2 = k1;
        (void)k2;
        assert( k2 == 1);
        Circulator::pointer         p1 = &k1;
        (void)p1;
        assert( (*p1) == 1);
        k1 = 3;
        assert( k1 == 3);
        assert( k2 == 3);
        assert( (*p1) == 3);
        k1 = 6;
        assert( k1 == 6);
        assert( k2 == 6);
        assert( (*p1) == 6);
        Circulator::size_type s = 5;
        assert( s == 5);
        Circulator::difference_type d = -5;
        assert( d == -5);

        // Check tests for empty data structures.
        Circulator z = Circulator();
        assert(   z == nullptr);
        assert( ! (z != nullptr));
        Circulator i = begin;
        assert( ! (i == nullptr));
        assert(   i != nullptr);
        assert( i == begin);
        assert( i == begin);
        // Do I reach myself.
        ++i;
        Circulator j = i;
        int k = 0;
        do {
            assert( k < 5);
            ++k;
            ++i;
        } while( i != j);
        assert( k == 5);
    }
    { // Open own scope to hide local variables.
        // Do I reach myself backwards.
        Circulator i = begin;
        ++i;
        Circulator j = i;
        int k = 0;
        do {
            assert( k < 5);
            ++k;
            --i;
        } while( i != j);
        assert( k == 5);
    }
    { // Open own scope to hide local variables.
        Circulator::difference_type d = begin - begin;
        assert( d == 0);
        d = begin - begin;
        assert( d == 0);
        Circulator i = begin + 1;
        assert( begin - i == 1 ||  begin - i == -1);
        assert( i - begin == 1 ||  i - begin == -1);
        // Check minimal circulator properties.
        i = i.min_circulator();
        Circulator j = i;
        assert( j - i == 0);
        j++;
        assert( j - i == 1);
        j++;
        assert( j - i == 2);
        j++;
        assert( j - i == 3);
        j++;
        assert( j - i == 4);
        j++;
        assert( j - i == 0);
    }

    Vector v2 = v;
    const Vector& v1 = v2;
    typedef Vector::const_iterator ConstIterBase;
    typedef Random_access_const_circulator_from_iterator<
      ConstIterBase, int, std::size_t, std::ptrdiff_t>
    ConstCircBase;
typedef Circulator_project<ConstCircBase,Ident,const int&,
  const int*> C_Circulator;
C_Circulator c_begin(ConstCircBase(v1.begin(),v1.end()));
Assert_random_access_category(c_begin);
{ // Open own scope to hide local variables.
    // Check generally correct parameter properties.
    CGAL::Assert_circulator_or_iterator(c_begin);
    CGAL::Assert_circulator_or_iterator(c_begin);
    CGAL::Assert_is_at_least_forward_category(c_begin);
    CGAL::Assert_is_at_least_forward_category(c_begin);
    typedef std::iterator_traits< C_Circulator >::value_type      VT;
    typedef std::iterator_traits< C_Circulator >::difference_type DT;
    assert(1==test_value_type(static_cast< VT* >(nullptr)));
    assert(1==test_distance_type(static_cast< DT* >(nullptr)));

    // Default constructor.
    C_Circulator z = C_Circulator();
    CGAL::Assert_circulator_or_iterator(z);
    // Copy constructor.
    C_Circulator i = c_begin;

    // Check general support for circulators and iterators.
    assert( CGAL::is_empty_range( z, z));
    assert( ! CGAL::is_empty_range( i, c_begin));

    int su = 0;
    int k  = 1;
    // Check general loop, pre-increment, dereference.
    if (! CGAL::is_empty_range( i, c_begin)) {   // superfluous
        do {
            assert( k == (*i));
            su += (*i);
            ++k;
            C_Circulator j = ++i;
            assert( i ==  j);
            if ( i != c_begin) {
                assert( (*i) == (*j));
            }
        } while (i != c_begin);  // Inequality and equality checked.
    }
    assert( i == c_begin);  // Equality checked.
    assert( su == 15);

    // Assignment.
    i = c_begin;
    su = 0;
    k  = 1;
    // Loop with post increment.
    if (! CGAL::is_empty_range( i, c_begin)) {   // superfluous
        do {
            assert( k == (*i));
            su += (*i);
            ++k;
            C_Circulator j = i++;
            assert(  i !=  j);
            if ( i != c_begin) {
                assert( (*i) == (*j) + 1);
            }
        } while (i != c_begin);
    }
    assert( i == c_begin);
    assert( su == 15);
}
{ // Open own scope to hide local variables.
    CGAL::Assert_is_at_least_bidirectional_category(c_begin);
    CGAL::Assert_is_at_least_bidirectional_category(c_begin);
    // Loop backwards and pre-decrement.
    C_Circulator i = c_begin;
    int su = 0;
    int k  = 5;
    do {
        C_Circulator j = --i;
        assert(  i ==  j);
        assert( (*i) == (*j));
        assert( k == (*i));
        su += (*i);
        --k;
    } while (i != c_begin);
    assert( i == c_begin);
    assert( su == 15);

    // Assignment.
    i = c_begin;
    su = 0;
    k  = 5;
    // Loop with post-decrement.
    do {
        C_Circulator j = i--;
        assert(  i !=  j);
        if ( j != c_begin) {
            assert( (*i) == (*j) - 1);
        }
        assert( k == (*i));
        su += (*i);
        --k;
    } while (i != c_begin);
    assert( i == c_begin);
    assert( su == 15);
}
{ // Open own scope to hide local variables.
    CGAL::Assert_is_at_least_random_access_category(c_begin);
    CGAL::Assert_is_at_least_random_access_category(c_begin);
    // Random access.
    int k;
    for( k = 0; k < 5; k++) {
        assert( 1+k == c_begin[k]);
    }

      int su = c_begin[0]
             + c_begin[1]
             + c_begin[2]
             + c_begin[3]
             + c_begin[4];
    assert( su == 15);

    // Jump around.
    C_Circulator i = c_begin;
    i += 3;
    assert( 4 == (*i));
    i -= 2;
    assert( 2 == (*i));
    i += 3;
    assert( 5 == (*i));
    i -= 4;
    assert( 1 == (*i));
    assert( i == c_begin);
    C_Circulator j = i + 3;
    assert( 4 == (*j));
    C_Circulator jj = j - 2;
    assert( 2 == (*jj));
    typedef std::ptrdiff_t PT;
    jj = PT(4) + jj;
    assert( jj == c_begin);
    C_Circulator ij = jj - 5;
    ij = jj - 5; // avoids warning with NDEBUG
    assert( ij == c_begin);

    // Difference test.
    assert( jj - i == 5  ||  jj - i == 0);
    assert( i + (j-i) == j);
    assert( (j-i) + i == j);
}
{ // Open own scope to hide local variables.
    CGAL::Assert_circulator( c_begin);
    CGAL::Assert_circulator( c_begin);

    // Check the local type parameters.
    C_Circulator::value_type      k1;
    k1 = 1;
    C_Circulator::reference       k2 = k1;
    (void)k2;
    assert( k2 == 1);
    C_Circulator::pointer         p1 = &k1;
    (void)p1;
    assert( (*p1) == 1);
    k1 = 3;
    assert( k1 == 3);
    assert( k2 == 3);
    assert( (*p1) == 3);
    k1 = 6;
    assert( k1 == 6);
    assert( k2 == 6);
    assert( (*p1) == 6);
    C_Circulator::size_type s = 5;
    assert( s == 5);
    C_Circulator::difference_type d = -5;
    assert( d == -5);

    // Check tests for empty data structures.
    C_Circulator z = C_Circulator();
    assert(   z == nullptr);
    assert( ! (z != nullptr));
    C_Circulator i = c_begin;
    assert( ! (i == nullptr));
    assert(   i != nullptr);
    assert( i == c_begin);
    assert( i == c_begin);
    // Do I reach myself.
    ++i;
    C_Circulator j = i;
    int k = 0;
    do {
        assert( k < 5);
        ++k;
        ++i;
    } while( i != j);
    assert( k == 5);
}
{ // Open own scope to hide local variables.
    // Do I reach myself backwards.
    C_Circulator i = c_begin;
    ++i;
    C_Circulator j = i;
    int k = 0;
    do {
        assert( k < 5);
        ++k;
        --i;
    } while( i != j);
    assert( k == 5);
}
{ // Open own scope to hide local variables.
    C_Circulator::difference_type d = c_begin - c_begin;
    assert( d == 0);
    d = c_begin - c_begin;
    assert( d == 0);
    C_Circulator i = c_begin + 1;
    assert( c_begin - i == 1 ||  c_begin - i == -1);
    assert( i - c_begin == 1 ||  i - c_begin == -1);
    // Check minimal circulator properties.
    i = i.min_circulator();
    C_Circulator j = i;
    assert( j - i == 0);
    j++;
    assert( j - i == 1);
    j++;
    assert( j - i == 2);
    j++;
    assert( j - i == 3);
    j++;
    assert( j - i == 4);
    j++;
    assert( j - i == 0);
}
  }
}
struct NN {
  NN* nn;
  int key;
  NN() : nn(nullptr), key(-1) {}
  NN( int k, NN* p) : nn(p), key(k) {}
  NN*       next()       { return nn; }
  const NN* next() const { return nn; }
};
int test_value_type( NN*) { return 1;}

void test_Circulator_on_node() {
  {
    NN* end   = new NN( 5, nullptr);
    NN* p     = new NN( 4, end);
    NN* start = new NN( 3, p);
    p         = new NN( 2, start);
    start     = new NN( 1, p);
    end->nn   = start;
    typedef Circulator_on_node< NN,
      Project_next<NN>,
                                Project_next<NN>,
                                NN&,
                                NN*,
                                Forward_circulator_tag>
    Circulator;
    Circulator begin( start);
    Assert_forward_category(begin);
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_is_at_least_forward_category(begin);
        CGAL::Assert_is_at_least_forward_category(begin);
        typedef std::iterator_traits< Circulator >::value_type      VT;
        typedef std::iterator_traits< Circulator >::difference_type DT;
        assert(1==test_value_type(static_cast< VT* >(nullptr)));
        assert(1==test_distance_type(static_cast< DT* >(nullptr)));

        // Default constructor.
        Circulator z = Circulator();
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Circulator i = begin;

        // Check general support for circulators and iterators.
        assert( CGAL::is_empty_range( z, z));
        assert( ! CGAL::is_empty_range( i, begin));

        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                assert( k == (*i).key);
                su += (*i).key;
                ++k;
                Circulator j = ++i;
                assert( i ==  j);
                if ( i != begin) {
                    assert( (*i).key == (*j).key);
                }
            } while (i != begin);  // Inequality and equality checked.
        }
        assert( i == begin);  // Equality checked.
        assert( su == 15);

        // Assignment.
        i = begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                assert( k == (*i).key);
                su += (*i).key;
                ++k;
                Circulator j = i++;
                assert(  i !=  j);
                if ( i != begin) {
                    assert( (*i).key == (*j).key + 1);
                }
            } while (i != begin);
        }
        assert( i == begin);
        assert( su == 15);
    }
    { // Open own scope to hide local variables.
        // Change three elements and check post-/pre-increment.
        Circulator i = begin;
        (*i++).key = 4;
        assert( 4 == (*begin).key);
        assert( 2 == (*i).key);
        (*i++).key = 3;
        assert( 3 == (*i).key);
        (*++i).key = 7;
        assert( 7 == (*i).key);

        // Check the setting and reset these elements
        // to their original values.
        i = begin;
        assert( 4 == (*i).key);
        (*i).key = 1;
        i++;
        assert( 3 == (*i).key);
        (*i++).key = 2;
        assert( 3 == (*i).key);
        i++;
        assert( 7 == (*i).key);
        (*i).key = 4;

        // Check the resetting.
        i = begin;
        int k = 1;
        do {
            assert( k == (*i).key);
            ++i;
            ++k;
        } while (i != begin);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_circulator( begin);
        CGAL::Assert_circulator( begin);

        // Check the local type parameters.
        Circulator::value_type      k1;
        k1.key = 1;
        Circulator::reference       k2 = k1;
        (void)k2;
        assert( k2.key == 1);
        Circulator::pointer         p1 = &k1;
        (void)p1;
        assert( (*p1).key == 1);
        k1.key = 3;
        assert( k1.key == 3);
        assert( k2.key == 3);
        assert( (*p1).key == 3);
        k1.key = 6;
        assert( k1.key == 6);
        assert( k2.key == 6);
        assert( (*p1).key == 6);
        Circulator::size_type s = 5;
        assert( s == 5);
        Circulator::difference_type d = -5;
        assert( d == -5);

        // Check tests for empty data structures.
        Circulator z = Circulator();
        assert(   z == nullptr);
        assert( ! (z != nullptr));
        Circulator i = begin;
        assert( ! (i == nullptr));
        assert(   i != nullptr);
        assert( i == begin);
        assert( i == begin);
        // Do I reach myself.
        ++i;
        Circulator j = i;
        int k = 0;
        do {
            assert( k < 5);
            ++k;
            ++i;
        } while( i != j);
        assert( k == 5);
    }

    typedef Circulator_on_node< NN,
      Project_next<NN>,
                                Project_next<NN>,
                                const NN&,
                                const NN*,
                                Forward_circulator_tag>
    C_Circulator;
    C_Circulator c_begin( start);
    Assert_forward_category(c_begin);
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(c_begin);
        CGAL::Assert_circulator_or_iterator(c_begin);
        CGAL::Assert_is_at_least_forward_category(c_begin);
        CGAL::Assert_is_at_least_forward_category(c_begin);
        typedef std::iterator_traits< C_Circulator >::value_type      VT;
        typedef std::iterator_traits< C_Circulator >::difference_type DT;
        assert(1==test_value_type(static_cast< VT* >(nullptr)));
        assert(1==test_distance_type(static_cast< DT* >(nullptr)));

        // Default constructor.
        C_Circulator z = C_Circulator();
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        C_Circulator i = c_begin;

        // Check general support for circulators and iterators.
        assert( CGAL::is_empty_range( z, z));
        assert( ! CGAL::is_empty_range( i, c_begin));

        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, c_begin)) {   // superfluous
            do {
                assert( k == (*i).key);
                su += (*i).key;
                ++k;
                C_Circulator j = ++i;
                assert( i ==  j);
                if ( i != c_begin) {
                    assert( (*i).key == (*j).key);
                }
            } while (i != c_begin);  // Inequality and equality checked.
        }
        assert( i == c_begin);  // Equality checked.
        assert( su == 15);

        // Assignment.
        i = c_begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, c_begin)) {   // superfluous
            do {
                assert( k == (*i).key);
                su += (*i).key;
                ++k;
                C_Circulator j = i++;
                assert(  i !=  j);
                if ( i != c_begin) {
                    assert( (*i).key == (*j).key + 1);
                }
            } while (i != c_begin);
        }
        assert( i == c_begin);
        assert( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_circulator( c_begin);
        CGAL::Assert_circulator( c_begin);

        // Check the local type parameters.
        C_Circulator::value_type      k1;
        k1.key = 1;
        C_Circulator::reference       k2 = k1;
        (void)k2;
        assert( k2.key == 1);
        C_Circulator::pointer         p1 = &k1;
        (void)p1;
        assert( (*p1).key == 1);
        k1.key = 3;
        assert( k1.key == 3);
        assert( k2.key == 3);
        assert( (*p1).key == 3);
        k1.key = 6;
        assert( k1.key == 6);
        assert( k2.key == 6);
        assert( (*p1).key == 6);
        C_Circulator::size_type s = 5;
        assert( s == 5);
        C_Circulator::difference_type d = -5;
        assert( d == -5);

        // Check tests for empty data structures.
        C_Circulator z = C_Circulator();
        assert(   z == nullptr);
        assert( ! (z != nullptr));
        C_Circulator i = c_begin;
        assert( ! (i == nullptr));
        assert(   i != nullptr);
        assert( i == c_begin);
        assert( i == c_begin);
        // Do I reach myself.
        ++i;
        C_Circulator j = i;
        int k = 0;
        do {
            assert( k < 5);
            ++k;
            ++i;
        } while( i != j);
        assert( k == 5);
    }
    while ( start != end) {
      p = start->nn;
      delete start;
      start = p;
    }
    delete end;
  }
}
void test_N_step_adaptor() {
  {
    typedef In_place_list<item,false> List;
    List l;
    l.push_back( *new item(1));
    l.push_back( *new item(1));
    l.push_back( *new item(2));
    l.push_back( *new item(2));
    l.push_back( *new item(3));
    l.push_back( *new item(3));
    l.push_back( *new item(4));
    l.push_back( *new item(4));
    l.push_back( *new item(5));
    l.push_back( *new item(5));
    assert( l.size() == 10);
    typedef List::iterator IterBase;
    typedef N_step_adaptor<IterBase,2> Iterator;
    Iterator begin(l.begin());
    Iterator end(l.end());
    Assert_bidirectional_category(begin);
    Assert_bidirectional_category(end);
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_circulator_or_iterator(end);
        CGAL::Assert_is_at_least_forward_category(begin);
        CGAL::Assert_is_at_least_forward_category(end);
        typedef std::iterator_traits< Iterator >::value_type      VT;
        typedef std::iterator_traits< Iterator >::difference_type DT;
        assert(1==test_value_type(static_cast< VT* >(nullptr)));
        assert(1==test_distance_type(static_cast< DT* >(nullptr)));

        // Default constructor.
        Iterator z ;
        CGAL::Assert_circulator_or_iterator(z);
        z = end ;
        // Copy constructor.
        Iterator i = begin;

        // Check general support for circulators and iterators.
        assert( CGAL::is_empty_range( z, z));
        assert( ! CGAL::is_empty_range( i, end));

        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, end)) {   // superfluous
            do {
                assert( k == (*i).key);
                su += (*i).key;
                ++k;
                Iterator j = ++i;
                assert( i ==  j);
                if ( i != end) {
                    assert( (*i).key == (*j).key);
                }
            } while (i != end);  // Inequality and equality checked.
        }
        assert( i == end);  // Equality checked.
        assert( su == 15);

        // Assignment.
        i = begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, end)) {   // superfluous
            do {
                assert( k == (*i).key);
                su += (*i).key;
                ++k;
                Iterator j = i++;
                assert(  i !=  j);
                if ( i != end) {
                    assert( (*i).key == (*j).key + 1);
                }
            } while (i != end);
        }
        assert( i == end);
        assert( su == 15);
    }
    { // Open own scope to hide local variables.
        // Change three elements and check post-/pre-increment.
        Iterator i = begin;
        (*i++).key = 4;
        assert( 4 == (*begin).key);
        assert( 2 == (*i).key);
        (*i++).key = 3;
        assert( 3 == (*i).key);
        (*++i).key = 7;
        assert( 7 == (*i).key);

        // Check the setting and reset these elements
        // to their original values.
        i = begin;
        assert( 4 == (*i).key);
        (*i).key = 1;
        i++;
        assert( 3 == (*i).key);
        (*i++).key = 2;
        assert( 3 == (*i).key);
        i++;
        assert( 7 == (*i).key);
        (*i).key = 4;

        // Check the resetting.
        i = begin;
        int k = 1;
        do {
            assert( k == (*i).key);
            ++i;
            ++k;
        } while (i != end);
    }
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_circulator_or_iterator(end);
        CGAL::Assert_is_at_least_forward_category(begin);
        CGAL::Assert_is_at_least_forward_category(end);
        typedef std::iterator_traits< Iterator >::value_type      VT;
        typedef std::iterator_traits< Iterator >::difference_type DT;
        assert(1==test_value_type(static_cast< VT* >(nullptr)));
        assert(1==test_distance_type(static_cast< DT* >(nullptr)));

        // Default constructor.
        Iterator z ;
        CGAL::Assert_circulator_or_iterator(z);
        z = end ;
        // Copy constructor.
        Iterator i = begin;

        // Check general support for circulators and iterators.
        assert( CGAL::is_empty_range( z, z));
        assert( ! CGAL::is_empty_range( i, end));

        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, end)) {   // superfluous
            do {
                assert( k == (*i).key);
                su += (*i).key;
                ++k;
                Iterator j = ++i;
                assert( i ==  j);
                if ( i != end) {
                    assert( (*i).key == (*j).key);
                }
            } while (i != end);  // Inequality and equality checked.
        }
        assert( i == end);  // Equality checked.
        assert( su == 15);

        // Assignment.
        i = begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, end)) {   // superfluous
            do {
                assert( k == (*i).key);
                su += (*i).key;
                ++k;
                Iterator j = i++;
                assert(  i !=  j);
                if ( i != end) {
                    assert( (*i).key == (*j).key + 1);
                }
            } while (i != end);
        }
        assert( i == end);
        assert( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_bidirectional_category(begin);
        CGAL::Assert_is_at_least_bidirectional_category(end);
        // Loop backwards and pre-decrement.
        Iterator i = end;
        int su = 0;
        int k  = 5;
        do {
            Iterator j = --i;
            assert(  i ==  j);
            assert( (*i).key == (*j).key);
            assert( k == (*i).key);
            su += (*i).key;
            --k;
        } while (i != begin);
        assert( i == begin);
        assert( su == 15);

        // Assignment.
        i = end;
        su = 0;
        k  = 5;
        // Loop with post-decrement.
        do {
            Iterator j = i--;
            assert(  i !=  j);
            if ( j != end) {
                assert( (*i).key == (*j).key - 1);
            }
            assert( k == (*i).key);
            su += (*i).key;
            --k;
        } while (i != begin);
        assert( i == begin);
        assert( su == 15);
    }
    CGAL::Assert_iterator( begin);
    CGAL::Assert_iterator( end);

    List l2 = l;
    const List& l1 = l2;
    typedef List::const_iterator ConstIterBase;
    typedef N_step_adaptor<ConstIterBase,2> C_Iterator;
C_Iterator c_begin(l1.begin());
C_Iterator c_end(l1.end());
Assert_bidirectional_category(c_begin);
Assert_bidirectional_category(c_end);
{ // Open own scope to hide local variables.
    // Check generally correct parameter properties.
    CGAL::Assert_circulator_or_iterator(c_begin);
    CGAL::Assert_circulator_or_iterator(c_end);
    CGAL::Assert_is_at_least_forward_category(c_begin);
    CGAL::Assert_is_at_least_forward_category(c_end);
    typedef std::iterator_traits< C_Iterator >::value_type      VT;
    typedef std::iterator_traits< C_Iterator >::difference_type DT;
    assert(1==test_value_type(static_cast< VT* >(nullptr)));
    assert(1==test_distance_type(static_cast< DT* >(nullptr)));

    // Default constructor.
    C_Iterator z ;
    CGAL::Assert_circulator_or_iterator(z);
    z = c_end ;
    // Copy constructor.
    C_Iterator i = c_begin;

    // Check general support for circulators and iterators.
    assert( CGAL::is_empty_range( z, z));
    assert( ! CGAL::is_empty_range( i, c_end));

    int su = 0;
    int k  = 1;
    // Check general loop, pre-increment, dereference.
    if (! CGAL::is_empty_range( i, c_end)) {   // superfluous
        do {
            assert( k == (*i).key);
            su += (*i).key;
            ++k;
            C_Iterator j = ++i;
            assert( i ==  j);
            if ( i != c_end) {
                assert( (*i).key == (*j).key);
            }
        } while (i != c_end);  // Inequality and equality checked.
    }
    assert( i == c_end);  // Equality checked.
    assert( su == 15);

    // Assignment.
    i = c_begin;
    su = 0;
    k  = 1;
    // Loop with post increment.
    if (! CGAL::is_empty_range( i, c_end)) {   // superfluous
        do {
            assert( k == (*i).key);
            su += (*i).key;
            ++k;
            C_Iterator j = i++;
            assert(  i !=  j);
            if ( i != c_end) {
                assert( (*i).key == (*j).key + 1);
            }
        } while (i != c_end);
    }
    assert( i == c_end);
    assert( su == 15);
}
{ // Open own scope to hide local variables.
    CGAL::Assert_is_at_least_bidirectional_category(c_begin);
    CGAL::Assert_is_at_least_bidirectional_category(c_end);
    // Loop backwards and pre-decrement.
    C_Iterator i = c_end;
    int su = 0;
    int k  = 5;
    do {
        C_Iterator j = --i;
        assert(  i ==  j);
        assert( (*i).key == (*j).key);
        assert( k == (*i).key);
        su += (*i).key;
        --k;
    } while (i != c_begin);
    assert( i == c_begin);
    assert( su == 15);

    // Assignment.
    i = c_end;
    su = 0;
    k  = 5;
    // Loop with post-decrement.
    do {
        C_Iterator j = i--;
        assert(  i !=  j);
        if ( j != c_end) {
            assert( (*i).key == (*j).key - 1);
        }
        assert( k == (*i).key);
        su += (*i).key;
        --k;
    } while (i != c_begin);
    assert( i == c_begin);
    assert( su == 15);
}
CGAL::Assert_iterator( c_begin);
CGAL::Assert_iterator( c_end);
l.destroy();
l2.destroy();
  }
  {
    typedef std::vector<int> Vector;
    Vector v;
    v.push_back(1);
    v.push_back(1);
    v.push_back(2);
    v.push_back(2);
    v.push_back(3);
    v.push_back(3);
    v.push_back(4);
    v.push_back(4);
    v.push_back(5);
    v.push_back(5);
    assert( v.size() == 10);
    typedef Vector::iterator IterBase;
    typedef N_step_adaptor<IterBase,2> Iterator;
    Iterator begin(v.begin());
    Iterator end(v.end());
    Assert_random_access_category(begin);
    Assert_random_access_category(end);
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_circulator_or_iterator(end);
        CGAL::Assert_is_at_least_forward_category(begin);
        CGAL::Assert_is_at_least_forward_category(end);
        typedef std::iterator_traits< Iterator >::value_type      VT;
        typedef std::iterator_traits< Iterator >::difference_type DT;
        assert(1==test_value_type(static_cast< VT* >(nullptr)));
        assert(1==test_distance_type(static_cast< DT* >(nullptr)));

        // Default constructor.
        Iterator z ;
        CGAL::Assert_circulator_or_iterator(z);
        z = end ;
        // Copy constructor.
        Iterator i = begin;

        // Check general support for circulators and iterators.
        assert( CGAL::is_empty_range( z, z));
        assert( ! CGAL::is_empty_range( i, end));

        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, end)) {   // superfluous
            do {
                assert( k == (*i));
                su += (*i);
                ++k;
                Iterator j = ++i;
                assert( i ==  j);
                if ( i != end) {
                    assert( (*i) == (*j));
                }
            } while (i != end);  // Inequality and equality checked.
        }
        assert( i == end);  // Equality checked.
        assert( su == 15);

        // Assignment.
        i = begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, end)) {   // superfluous
            do {
                assert( k == (*i));
                su += (*i);
                ++k;
                Iterator j = i++;
                assert(  i !=  j);
                if ( i != end) {
                    assert( (*i) == (*j) + 1);
                }
            } while (i != end);
        }
        assert( i == end);
        assert( su == 15);
    }
    { // Open own scope to hide local variables.
        // Change three elements and check post-/pre-increment.
        Iterator i = begin;
        (*i++) = 4;
        assert( 4 == (*begin));
        assert( 2 == (*i));
        (*i++) = 3;
        assert( 3 == (*i));
        (*++i) = 7;
        assert( 7 == (*i));

        // Check the setting and reset these elements
        // to their original values.
        i = begin;
        assert( 4 == (*i));
        (*i) = 1;
        i++;
        assert( 3 == (*i));
        (*i++) = 2;
        assert( 3 == (*i));
        i++;
        assert( 7 == (*i));
        (*i) = 4;

        // Check the resetting.
        i = begin;
        int k = 1;
        do {
            assert( k == (*i));
            ++i;
            ++k;
        } while (i != end);
    }
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_circulator_or_iterator(end);
        CGAL::Assert_is_at_least_forward_category(begin);
        CGAL::Assert_is_at_least_forward_category(end);
        typedef std::iterator_traits< Iterator >::value_type      VT;
        typedef std::iterator_traits< Iterator >::difference_type DT;
        assert(1==test_value_type(static_cast< VT* >(nullptr)));
        assert(1==test_distance_type(static_cast< DT* >(nullptr)));

        // Default constructor.
        Iterator z ;
        CGAL::Assert_circulator_or_iterator(z);
        z = end ;
        // Copy constructor.
        Iterator i = begin;

        // Check general support for circulators and iterators.
        assert( CGAL::is_empty_range( z, z));
        assert( ! CGAL::is_empty_range( i, end));

        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, end)) {   // superfluous
            do {
                assert( k == (*i));
                su += (*i);
                ++k;
                Iterator j = ++i;
                assert( i ==  j);
                if ( i != end) {
                    assert( (*i) == (*j));
                }
            } while (i != end);  // Inequality and equality checked.
        }
        assert( i == end);  // Equality checked.
        assert( su == 15);

        // Assignment.
        i = begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, end)) {   // superfluous
            do {
                assert( k == (*i));
                su += (*i);
                ++k;
                Iterator j = i++;
                assert(  i !=  j);
                if ( i != end) {
                    assert( (*i) == (*j) + 1);
                }
            } while (i != end);
        }
        assert( i == end);
        assert( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_bidirectional_category(begin);
        CGAL::Assert_is_at_least_bidirectional_category(end);
        // Loop backwards and pre-decrement.
        Iterator i = end;
        int su = 0;
        int k  = 5;
        do {
            Iterator j = --i;
            assert(  i ==  j);
            assert( (*i) == (*j));
            assert( k == (*i));
            su += (*i);
            --k;
        } while (i != begin);
        assert( i == begin);
        assert( su == 15);

        // Assignment.
        i = end;
        su = 0;
        k  = 5;
        // Loop with post-decrement.
        do {
            Iterator j = i--;
            assert(  i !=  j);
            if ( j != end) {
                assert( (*i) == (*j) - 1);
            }
            assert( k == (*i));
            su += (*i);
            --k;
        } while (i != begin);
        assert( i == begin);
        assert( su == 15);
    }
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_circulator_or_iterator(end);
        CGAL::Assert_is_at_least_forward_category(begin);
        CGAL::Assert_is_at_least_forward_category(end);
        typedef std::iterator_traits< Iterator >::value_type      VT;
        typedef std::iterator_traits< Iterator >::difference_type DT;
        assert(1==test_value_type(static_cast< VT* >(nullptr)));
        assert(1==test_distance_type(static_cast< DT* >(nullptr)));

        // Default constructor.
        Iterator z ;
        CGAL::Assert_circulator_or_iterator(z);
        z = end ;
        // Copy constructor.
        Iterator i = begin;

        // Check general support for circulators and iterators.
        assert( CGAL::is_empty_range( z, z));
        assert( ! CGAL::is_empty_range( i, end));

        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, end)) {   // superfluous
            do {
                assert( k == (*i));
                su += (*i);
                ++k;
                Iterator j = ++i;
                assert( i ==  j);
                if ( i != end) {
                    assert( (*i) == (*j));
                }
            } while (i != end);  // Inequality and equality checked.
        }
        assert( i == end);  // Equality checked.
        assert( su == 15);

        // Assignment.
        i = begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, end)) {   // superfluous
            do {
                assert( k == (*i));
                su += (*i);
                ++k;
                Iterator j = i++;
                assert(  i !=  j);
                if ( i != end) {
                    assert( (*i) == (*j) + 1);
                }
            } while (i != end);
        }
        assert( i == end);
        assert( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_bidirectional_category(begin);
        CGAL::Assert_is_at_least_bidirectional_category(end);
        // Loop backwards and pre-decrement.
        Iterator i = end;
        int su = 0;
        int k  = 5;
        do {
            Iterator j = --i;
            assert(  i ==  j);
            assert( (*i) == (*j));
            assert( k == (*i));
            su += (*i);
            --k;
        } while (i != begin);
        assert( i == begin);
        assert( su == 15);

        // Assignment.
        i = end;
        su = 0;
        k  = 5;
        // Loop with post-decrement.
        do {
            Iterator j = i--;
            assert(  i !=  j);
            if ( j != end) {
                assert( (*i) == (*j) - 1);
            }
            assert( k == (*i));
            su += (*i);
            --k;
        } while (i != begin);
        assert( i == begin);
        assert( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_random_access_category(begin);
        CGAL::Assert_is_at_least_random_access_category(end);
        // Random access.
        int k;
        for( k = 0; k < 5; k++) {
            assert( 1+k == begin[k]);
        }

          int su = begin[0]
                 + begin[1]
                 + begin[2]
                 + begin[3]
                 + begin[4];
        assert( su == 15);

        // Jump around.
        Iterator i = begin;
        i += 3;
        assert( 4 == (*i));
        i -= 2;
        assert( 2 == (*i));
        i += 3;
        assert( 5 == (*i));
        i -= 4;
        assert( 1 == (*i));
        assert( i == begin);
        Iterator j = i + 3;
        assert( 4 == (*j));
        Iterator jj = j - 2;
        assert( 2 == (*jj));
        typedef std::ptrdiff_t PT;
        jj = PT(4) + jj;
        assert( jj == end);
        Iterator ij = jj - 5;
        ij = jj - 5; // avoids warning with NDEBUG
        assert( ij == begin);

        // Difference test.
        assert( jj - i == 5  ||  jj - i == 0);
        assert( i + (j-i) == j);
        assert( (j-i) + i == j);
    }
    { // Open own scope to hide local variables.
        Iterator i = begin;
        i[2] = 18;
        i[4] = 9;
        i[3] = 12;
        assert( i[2] == 18);
        assert( i[4] == 9);
        assert( i[3] == 12);
        i[2] = 3;
        i[3] = 4;
        i[4] = 5;
        // Check the resetting.
        i = begin;
        int k = 1;
        do {
            assert( k == (*i));
            ++i;
            ++k;
        } while (i != end);
    }
    CGAL::Assert_iterator( begin);
    CGAL::Assert_iterator( end);
    { // Open own scope to hide local variables.
        assert( end - begin ==  5);
        assert( begin - end == -5);
        // Relational operator.
        Iterator i = begin;
        ++i;
        Iterator j = i;
        ++j;
        assert( begin < i);
        assert( i < j);
        assert( j < end);
        assert( j > i);
        assert( i <= j);
        assert( j >= i);
        assert( i <= i);
        assert( i >= i);

        assert( !( i >= j));
        assert( !( j <= i));
        assert( !( i > j));
        assert( !( j < i));
        assert( !( i > i));
        assert( !( i < i));
    }

    Vector v2 = v;
    const Vector& v1 = v2;
    typedef Vector::const_iterator ConstIterBase;
    typedef N_step_adaptor< ConstIterBase,2> C_Iterator;
    C_Iterator c_begin(v1.begin());
    C_Iterator c_end(v1.end());
    Assert_random_access_category(c_begin);
    Assert_random_access_category(c_end);
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(c_begin);
        CGAL::Assert_circulator_or_iterator(c_end);
        CGAL::Assert_is_at_least_forward_category(c_begin);
        CGAL::Assert_is_at_least_forward_category(c_end);
        typedef std::iterator_traits< C_Iterator >::value_type      VT;
        typedef std::iterator_traits< C_Iterator >::difference_type DT;
        assert(1==test_value_type(static_cast< VT* >(nullptr)));
        assert(1==test_distance_type(static_cast< DT* >(nullptr)));

        // Default constructor.
        C_Iterator z ;
        CGAL::Assert_circulator_or_iterator(z);
        z = c_end ;
        // Copy constructor.
        C_Iterator i = c_begin;

        // Check general support for circulators and iterators.
        assert( CGAL::is_empty_range( z, z));
        assert( ! CGAL::is_empty_range( i, c_end));

        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, c_end)) {   // superfluous
            do {
                assert( k == (*i));
                su += (*i);
                ++k;
                C_Iterator j = ++i;
                assert( i ==  j);
                if ( i != c_end) {
                    assert( (*i) == (*j));
                }
            } while (i != c_end);  // Inequality and equality checked.
        }
        assert( i == c_end);  // Equality checked.
        assert( su == 15);

        // Assignment.
        i = c_begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, c_end)) {   // superfluous
            do {
                assert( k == (*i));
                su += (*i);
                ++k;
                C_Iterator j = i++;
                assert(  i !=  j);
                if ( i != c_end) {
                    assert( (*i) == (*j) + 1);
                }
            } while (i != c_end);
        }
        assert( i == c_end);
        assert( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_bidirectional_category(c_begin);
        CGAL::Assert_is_at_least_bidirectional_category(c_end);
        // Loop backwards and pre-decrement.
        C_Iterator i = c_end;
        int su = 0;
        int k  = 5;
        do {
            C_Iterator j = --i;
            assert(  i ==  j);
            assert( (*i) == (*j));
            assert( k == (*i));
            su += (*i);
            --k;
        } while (i != c_begin);
        assert( i == c_begin);
        assert( su == 15);

        // Assignment.
        i = c_end;
        su = 0;
        k  = 5;
        // Loop with post-decrement.
        do {
            C_Iterator j = i--;
            assert(  i !=  j);
            if ( j != c_end) {
                assert( (*i) == (*j) - 1);
            }
            assert( k == (*i));
            su += (*i);
            --k;
        } while (i != c_begin);
        assert( i == c_begin);
        assert( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_random_access_category(c_begin);
        CGAL::Assert_is_at_least_random_access_category(c_end);
        // Random access.
        int k;
        for( k = 0; k < 5; k++) {
            assert( 1+k == c_begin[k]);
        }

          int su = c_begin[0]
                 + c_begin[1]
                 + c_begin[2]
                 + c_begin[3]
                 + c_begin[4];
        assert( su == 15);

        // Jump around.
        C_Iterator i = c_begin;
        i += 3;
        assert( 4 == (*i));
        i -= 2;
        assert( 2 == (*i));
        i += 3;
        assert( 5 == (*i));
        i -= 4;
        assert( 1 == (*i));
        assert( i == c_begin);
        C_Iterator j = i + 3;
        assert( 4 == (*j));
        C_Iterator jj = j - 2;
        assert( 2 == (*jj));
        typedef std::ptrdiff_t PT;
        jj = PT(4) + jj;
        assert( jj == c_end);
        C_Iterator ij = jj - 5;
        ij = jj - 5; // avoids warning with NDEBUG
        assert( ij == c_begin);

        // Difference test.
        assert( jj - i == 5  ||  jj - i == 0);
        assert( i + (j-i) == j);
        assert( (j-i) + i == j);
    }
    CGAL::Assert_iterator( c_begin);
    CGAL::Assert_iterator( c_end);
    { // Open own scope to hide local variables.
        assert( c_end - c_begin ==  5);
        assert( c_begin - c_end == -5);
        // Relational operator.
        C_Iterator i = c_begin;
        ++i;
        C_Iterator j = i;
        ++j;
        assert( c_begin < i);
        assert( i < j);
        assert( j < c_end);
        assert( j > i);
        assert( i <= j);
        assert( j >= i);
        assert( i <= i);
        assert( i >= i);

        assert( !( i >= j));
        assert( !( j <= i));
        assert( !( i > j));
        assert( !( j < i));
        assert( !( i > i));
        assert( !( i < i));
    }
  }
  {
    typedef std::vector<int> Vector;
    Vector v;
    v.push_back(1);
    v.push_back(1);
    v.push_back(2);
    v.push_back(2);
    v.push_back(3);
    v.push_back(3);
    v.push_back(4);
    v.push_back(4);
    v.push_back(5);
    v.push_back(5);
    assert( v.size() == 10);
    typedef Vector::iterator IterBase;
    typedef Random_access_circulator_from_iterator<IterBase,int,
      std::size_t,std::ptrdiff_t> CircBase;
    typedef N_step_adaptor<CircBase,2> Circulator;
    Circulator begin(CircBase( v.begin(),v.end()));
    Assert_random_access_category(begin);
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_is_at_least_forward_category(begin);
        CGAL::Assert_is_at_least_forward_category(begin);
        typedef std::iterator_traits< Circulator >::value_type      VT;
        typedef std::iterator_traits< Circulator >::difference_type DT;
        assert(1==test_value_type(static_cast< VT* >(nullptr)));
        assert(1==test_distance_type(static_cast< DT* >(nullptr)));

        // Default constructor.
        Circulator z = Circulator();
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Circulator i = begin;

        // Check general support for circulators and iterators.
        assert( CGAL::is_empty_range( z, z));
        assert( ! CGAL::is_empty_range( i, begin));

        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                assert( k == (*i));
                su += (*i);
                ++k;
                Circulator j = ++i;
                assert( i ==  j);
                if ( i != begin) {
                    assert( (*i) == (*j));
                }
            } while (i != begin);  // Inequality and equality checked.
        }
        assert( i == begin);  // Equality checked.
        assert( su == 15);

        // Assignment.
        i = begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                assert( k == (*i));
                su += (*i);
                ++k;
                Circulator j = i++;
                assert(  i !=  j);
                if ( i != begin) {
                    assert( (*i) == (*j) + 1);
                }
            } while (i != begin);
        }
        assert( i == begin);
        assert( su == 15);
    }
    { // Open own scope to hide local variables.
        // Change three elements and check post-/pre-increment.
        Circulator i = begin;
        (*i++) = 4;
        assert( 4 == (*begin));
        assert( 2 == (*i));
        (*i++) = 3;
        assert( 3 == (*i));
        (*++i) = 7;
        assert( 7 == (*i));

        // Check the setting and reset these elements
        // to their original values.
        i = begin;
        assert( 4 == (*i));
        (*i) = 1;
        i++;
        assert( 3 == (*i));
        (*i++) = 2;
        assert( 3 == (*i));
        i++;
        assert( 7 == (*i));
        (*i) = 4;

        // Check the resetting.
        i = begin;
        int k = 1;
        do {
            assert( k == (*i));
            ++i;
            ++k;
        } while (i != begin);
    }
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_is_at_least_forward_category(begin);
        CGAL::Assert_is_at_least_forward_category(begin);
        typedef std::iterator_traits< Circulator >::value_type      VT;
        typedef std::iterator_traits< Circulator >::difference_type DT;
        assert(1==test_value_type(static_cast< VT* >(nullptr)));
        assert(1==test_distance_type(static_cast< DT* >(nullptr)));

        // Default constructor.
        Circulator z = Circulator();
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Circulator i = begin;

        // Check general support for circulators and iterators.
        assert( CGAL::is_empty_range( z, z));
        assert( ! CGAL::is_empty_range( i, begin));

        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                assert( k == (*i));
                su += (*i);
                ++k;
                Circulator j = ++i;
                assert( i ==  j);
                if ( i != begin) {
                    assert( (*i) == (*j));
                }
            } while (i != begin);  // Inequality and equality checked.
        }
        assert( i == begin);  // Equality checked.
        assert( su == 15);

        // Assignment.
        i = begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                assert( k == (*i));
                su += (*i);
                ++k;
                Circulator j = i++;
                assert(  i !=  j);
                if ( i != begin) {
                    assert( (*i) == (*j) + 1);
                }
            } while (i != begin);
        }
        assert( i == begin);
        assert( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_bidirectional_category(begin);
        CGAL::Assert_is_at_least_bidirectional_category(begin);
        // Loop backwards and pre-decrement.
        Circulator i = begin;
        int su = 0;
        int k  = 5;
        do {
            Circulator j = --i;
            assert(  i ==  j);
            assert( (*i) == (*j));
            assert( k == (*i));
            su += (*i);
            --k;
        } while (i != begin);
        assert( i == begin);
        assert( su == 15);

        // Assignment.
        i = begin;
        su = 0;
        k  = 5;
        // Loop with post-decrement.
        do {
            Circulator j = i--;
            assert(  i !=  j);
            if ( j != begin) {
                assert( (*i) == (*j) - 1);
            }
            assert( k == (*i));
            su += (*i);
            --k;
        } while (i != begin);
        assert( i == begin);
        assert( su == 15);
    }
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_is_at_least_forward_category(begin);
        CGAL::Assert_is_at_least_forward_category(begin);
        typedef std::iterator_traits< Circulator >::value_type      VT;
        typedef std::iterator_traits< Circulator >::difference_type DT;
        assert(1==test_value_type(static_cast< VT* >(nullptr)));
        assert(1==test_distance_type(static_cast< DT* >(nullptr)));

        // Default constructor.
        Circulator z = Circulator();
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Circulator i = begin;

        // Check general support for circulators and iterators.
        assert( CGAL::is_empty_range( z, z));
        assert( ! CGAL::is_empty_range( i, begin));

        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                assert( k == (*i));
                su += (*i);
                ++k;
                Circulator j = ++i;
                assert( i ==  j);
                if ( i != begin) {
                    assert( (*i) == (*j));
                }
            } while (i != begin);  // Inequality and equality checked.
        }
        assert( i == begin);  // Equality checked.
        assert( su == 15);

        // Assignment.
        i = begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                assert( k == (*i));
                su += (*i);
                ++k;
                Circulator j = i++;
                assert(  i !=  j);
                if ( i != begin) {
                    assert( (*i) == (*j) + 1);
                }
            } while (i != begin);
        }
        assert( i == begin);
        assert( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_bidirectional_category(begin);
        CGAL::Assert_is_at_least_bidirectional_category(begin);
        // Loop backwards and pre-decrement.
        Circulator i = begin;
        int su = 0;
        int k  = 5;
        do {
            Circulator j = --i;
            assert(  i ==  j);
            assert( (*i) == (*j));
            assert( k == (*i));
            su += (*i);
            --k;
        } while (i != begin);
        assert( i == begin);
        assert( su == 15);

        // Assignment.
        i = begin;
        su = 0;
        k  = 5;
        // Loop with post-decrement.
        do {
            Circulator j = i--;
            assert(  i !=  j);
            if ( j != begin) {
                assert( (*i) == (*j) - 1);
            }
            assert( k == (*i));
            su += (*i);
            --k;
        } while (i != begin);
        assert( i == begin);
        assert( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_random_access_category(begin);
        CGAL::Assert_is_at_least_random_access_category(begin);
        // Random access.
        int k;
        for( k = 0; k < 5; k++) {
            assert( 1+k == begin[k]);
        }

          int su = begin[0]
                 + begin[1]
                 + begin[2]
                 + begin[3]
                 + begin[4];
        assert( su == 15);

        // Jump around.
        Circulator i = begin;
        i += 3;
        assert( 4 == (*i));
        i -= 2;
        assert( 2 == (*i));
        i += 3;
        assert( 5 == (*i));
        i -= 4;
        assert( 1 == (*i));
        assert( i == begin);
        Circulator j = i + 3;
        assert( 4 == (*j));
        Circulator jj = j - 2;
        assert( 2 == (*jj));
        typedef std::ptrdiff_t PT;
        jj = PT(4) + jj;
        assert( jj == begin);
        Circulator ij = jj - 5;
        ij = jj - 5; // avoids warning with NDEBUG
        assert( ij == begin);

        // Difference test.
        assert( jj - i == 5  ||  jj - i == 0);
        assert( i + (j-i) == j);
        assert( (j-i) + i == j);
    }
    { // Open own scope to hide local variables.
        Circulator i = begin;
        i[2] = 18;
        i[4] = 9;
        i[3] = 12;
        assert( i[2] == 18);
        assert( i[4] == 9);
        assert( i[3] == 12);
        i[2] = 3;
        i[3] = 4;
        i[4] = 5;
        // Check the resetting.
        i = begin;
        int k = 1;
        do {
            assert( k == (*i));
            ++i;
            ++k;
        } while (i != begin);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_circulator( begin);
        CGAL::Assert_circulator( begin);

        // Check the local type parameters.
        Circulator::value_type      k1;
        k1 = 1;
        Circulator::reference       k2 = k1;
        (void)k2;
        assert( k2 == 1);
        Circulator::pointer         p1 = &k1;
        (void)p1;
        assert( (*p1) == 1);
        k1 = 3;
        assert( k1 == 3);
        assert( k2 == 3);
        assert( (*p1) == 3);
        k1 = 6;
        assert( k1 == 6);
        assert( k2 == 6);
        assert( (*p1) == 6);
        Circulator::size_type s = 5;
        assert( s == 5);
        Circulator::difference_type d = -5;
        assert( d == -5);

        // Check tests for empty data structures.
        Circulator z = Circulator();
        assert(   z == nullptr);
        assert( ! (z != nullptr));
        Circulator i = begin;
        assert( ! (i == nullptr));
        assert(   i != nullptr);
        assert( i == begin);
        assert( i == begin);
        // Do I reach myself.
        ++i;
        Circulator j = i;
        int k = 0;
        do {
            assert( k < 5);
            ++k;
            ++i;
        } while( i != j);
        assert( k == 5);
    }
    { // Open own scope to hide local variables.
        // Do I reach myself backwards.
        Circulator i = begin;
        ++i;
        Circulator j = i;
        int k = 0;
        do {
            assert( k < 5);
            ++k;
            --i;
        } while( i != j);
        assert( k == 5);
    }
    { // Open own scope to hide local variables.
        Circulator::difference_type d = begin - begin;
        assert( d == 0);
        d = begin - begin;
        assert( d == 0);
        Circulator i = begin + 1;
        assert( begin - i == 1 ||  begin - i == -1);
        assert( i - begin == 1 ||  i - begin == -1);
        // Check minimal circulator properties.
        i = i.min_circulator();
        Circulator j = i;
        assert( j - i == 0);
        j++;
        assert( j - i == 1);
        j++;
        assert( j - i == 2);
        j++;
        assert( j - i == 3);
        j++;
        assert( j - i == 4);
        j++;
        assert( j - i == 0);
    }

    Vector v2 = v;
    const Vector& v1 = v2;
    typedef Vector::const_iterator ConstIterBase;
    typedef Random_access_const_circulator_from_iterator<
      ConstIterBase, int, std::size_t, std::ptrdiff_t>
    ConstCircBase;
typedef N_step_adaptor<ConstCircBase,2> C_Circulator;
C_Circulator c_begin(ConstCircBase(v1.begin(),v1.end()));
Assert_random_access_category(c_begin);
{ // Open own scope to hide local variables.
    // Check generally correct parameter properties.
    CGAL::Assert_circulator_or_iterator(c_begin);
    CGAL::Assert_circulator_or_iterator(c_begin);
    CGAL::Assert_is_at_least_forward_category(c_begin);
    CGAL::Assert_is_at_least_forward_category(c_begin);
    typedef std::iterator_traits< C_Circulator >::value_type      VT;
    typedef std::iterator_traits< C_Circulator >::difference_type DT;
    assert(1==test_value_type(static_cast< VT* >(nullptr)));
    assert(1==test_distance_type(static_cast< DT* >(nullptr)));

    // Default constructor.
    C_Circulator z = C_Circulator();
    CGAL::Assert_circulator_or_iterator(z);
    // Copy constructor.
    C_Circulator i = c_begin;

    // Check general support for circulators and iterators.
    assert( CGAL::is_empty_range( z, z));
    assert( ! CGAL::is_empty_range( i, c_begin));

    int su = 0;
    int k  = 1;
    // Check general loop, pre-increment, dereference.
    if (! CGAL::is_empty_range( i, c_begin)) {   // superfluous
        do {
            assert( k == (*i));
            su += (*i);
            ++k;
            C_Circulator j = ++i;
            assert( i ==  j);
            if ( i != c_begin) {
                assert( (*i) == (*j));
            }
        } while (i != c_begin);  // Inequality and equality checked.
    }
    assert( i == c_begin);  // Equality checked.
    assert( su == 15);

    // Assignment.
    i = c_begin;
    su = 0;
    k  = 1;
    // Loop with post increment.
    if (! CGAL::is_empty_range( i, c_begin)) {   // superfluous
        do {
            assert( k == (*i));
            su += (*i);
            ++k;
            C_Circulator j = i++;
            assert(  i !=  j);
            if ( i != c_begin) {
                assert( (*i) == (*j) + 1);
            }
        } while (i != c_begin);
    }
    assert( i == c_begin);
    assert( su == 15);
}
{ // Open own scope to hide local variables.
    CGAL::Assert_is_at_least_bidirectional_category(c_begin);
    CGAL::Assert_is_at_least_bidirectional_category(c_begin);
    // Loop backwards and pre-decrement.
    C_Circulator i = c_begin;
    int su = 0;
    int k  = 5;
    do {
        C_Circulator j = --i;
        assert(  i ==  j);
        assert( (*i) == (*j));
        assert( k == (*i));
        su += (*i);
        --k;
    } while (i != c_begin);
    assert( i == c_begin);
    assert( su == 15);

    // Assignment.
    i = c_begin;
    su = 0;
    k  = 5;
    // Loop with post-decrement.
    do {
        C_Circulator j = i--;
        assert(  i !=  j);
        if ( j != c_begin) {
            assert( (*i) == (*j) - 1);
        }
        assert( k == (*i));
        su += (*i);
        --k;
    } while (i != c_begin);
    assert( i == c_begin);
    assert( su == 15);
}
{ // Open own scope to hide local variables.
    CGAL::Assert_is_at_least_random_access_category(c_begin);
    CGAL::Assert_is_at_least_random_access_category(c_begin);
    // Random access.
    int k;
    for( k = 0; k < 5; k++) {
        assert( 1+k == c_begin[k]);
    }

      int su = c_begin[0]
             + c_begin[1]
             + c_begin[2]
             + c_begin[3]
             + c_begin[4];
    assert( su == 15);

    // Jump around.
    C_Circulator i = c_begin;
    i += 3;
    assert( 4 == (*i));
    i -= 2;
    assert( 2 == (*i));
    i += 3;
    assert( 5 == (*i));
    i -= 4;
    assert( 1 == (*i));
    assert( i == c_begin);
    C_Circulator j = i + 3;
    assert( 4 == (*j));
    C_Circulator jj = j - 2;
    assert( 2 == (*jj));
    typedef std::ptrdiff_t PT;
    jj = PT(4) + jj;
    assert( jj == c_begin);
    C_Circulator ij = jj - 5;
    ij = jj - 5; // avoids warning with NDEBUG
    assert( ij == c_begin);

    // Difference test.
    assert( jj - i == 5  ||  jj - i == 0);
    assert( i + (j-i) == j);
    assert( (j-i) + i == j);
}
{ // Open own scope to hide local variables.
    CGAL::Assert_circulator( c_begin);
    CGAL::Assert_circulator( c_begin);

    // Check the local type parameters.
    C_Circulator::value_type      k1;
    k1 = 1;
    C_Circulator::reference       k2 = k1;
    (void)k2;
    assert( k2 == 1);
    C_Circulator::pointer         p1 = &k1;
    (void)p1;
    assert( (*p1) == 1);
    k1 = 3;
    assert( k1 == 3);
    assert( k2 == 3);
    assert( (*p1) == 3);
    k1 = 6;
    assert( k1 == 6);
    assert( k2 == 6);
    assert( (*p1) == 6);
    C_Circulator::size_type s = 5;
    assert( s == 5);
    C_Circulator::difference_type d = -5;
    assert( d == -5);

    // Check tests for empty data structures.
    C_Circulator z = C_Circulator();
    assert(   z == nullptr);
    assert( ! (z != nullptr));
    C_Circulator i = c_begin;
    assert( ! (i == nullptr));
    assert(   i != nullptr);
    assert( i == c_begin);
    assert( i == c_begin);
    // Do I reach myself.
    ++i;
    C_Circulator j = i;
    int k = 0;
    do {
        assert( k < 5);
        ++k;
        ++i;
    } while( i != j);
    assert( k == 5);
}
{ // Open own scope to hide local variables.
    // Do I reach myself backwards.
    C_Circulator i = c_begin;
    ++i;
    C_Circulator j = i;
    int k = 0;
    do {
        assert( k < 5);
        ++k;
        --i;
    } while( i != j);
    assert( k == 5);
}
{ // Open own scope to hide local variables.
    C_Circulator::difference_type d = c_begin - c_begin;
    assert( d == 0);
    d = c_begin - c_begin;
    assert( d == 0);
    C_Circulator i = c_begin + 1;
    assert( c_begin - i == 1 ||  c_begin - i == -1);
    assert( i - c_begin == 1 ||  i - c_begin == -1);
    // Check minimal circulator properties.
    i = i.min_circulator();
    C_Circulator j = i;
    assert( j - i == 0);
    j++;
    assert( j - i == 1);
    j++;
    assert( j - i == 2);
    j++;
    assert( j - i == 3);
    j++;
    assert( j - i == 4);
    j++;
    assert( j - i == 0);
}
  }
}
void test_N_step_adaptor_derived() {
    {
        typedef In_place_list<item,false> List;
        List l;
        l.push_back( *new item(1));
        l.push_back( *new item(1));
        l.push_back( *new item(2));
        l.push_back( *new item(2));
        l.push_back( *new item(3));
        l.push_back( *new item(3));
        l.push_back( *new item(4));
        l.push_back( *new item(4));
        l.push_back( *new item(5));
        l.push_back( *new item(5));
        assert( l.size() == 10);
        typedef List::iterator IterBase;
        typedef N_step_adaptor_derived<IterBase,2> Iterator;
        Iterator begin(l.begin());
        Iterator end(l.end());
        Assert_bidirectional_category(begin);
        Assert_bidirectional_category(end);
        { // Open own scope to hide local variables.
            // Check generally correct parameter properties.
            CGAL::Assert_circulator_or_iterator(begin);
            CGAL::Assert_circulator_or_iterator(end);
            CGAL::Assert_is_at_least_forward_category(begin);
            CGAL::Assert_is_at_least_forward_category(end);
            typedef std::iterator_traits< Iterator >::value_type      VT;
            typedef std::iterator_traits< Iterator >::difference_type DT;
            assert(1==test_value_type(static_cast< VT* >(nullptr)));
            assert(1==test_distance_type(static_cast< DT* >(nullptr)));

            // Default constructor.
            Iterator z ;
            CGAL::Assert_circulator_or_iterator(z);
            z = end ;
            // Copy constructor.
            Iterator i = begin;

            // Check general support for circulators and iterators.
            assert( CGAL::is_empty_range( z, z));
            assert( ! CGAL::is_empty_range( i, end));

            int su = 0;
            int k  = 1;
            // Check general loop, pre-increment, dereference.
            if (! CGAL::is_empty_range( i, end)) {   // superfluous
                do {
                    assert( k == (*i).key);
                    su += (*i).key;
                    ++k;
                    Iterator j = ++i;
                    assert( i ==  j);
                    if ( i != end) {
                        assert( (*i).key == (*j).key);
                    }
                } while (i != end);  // Inequality and equality checked.
            }
            assert( i == end);  // Equality checked.
            assert( su == 15);

            // Assignment.
            i = begin;
            su = 0;
            k  = 1;
            // Loop with post increment.
            if (! CGAL::is_empty_range( i, end)) {   // superfluous
                do {
                    assert( k == (*i).key);
                    su += (*i).key;
                    ++k;
                    Iterator j = i++;
                    assert(  i !=  j);
                    if ( i != end) {
                        assert( (*i).key == (*j).key + 1);
                    }
                } while (i != end);
            }
            assert( i == end);
            assert( su == 15);
        }
        { // Open own scope to hide local variables.
            // Change three elements and check post-/pre-increment.
            Iterator i = begin;
            (*i++).key = 4;
            assert( 4 == (*begin).key);
            assert( 2 == (*i).key);
            (*i++).key = 3;
            assert( 3 == (*i).key);
            (*++i).key = 7;
            assert( 7 == (*i).key);

            // Check the setting and reset these elements
            // to their original values.
            i = begin;
            assert( 4 == (*i).key);
            (*i).key = 1;
            i++;
            assert( 3 == (*i).key);
            (*i++).key = 2;
            assert( 3 == (*i).key);
            i++;
            assert( 7 == (*i).key);
            (*i).key = 4;

            // Check the resetting.
            i = begin;
            int k = 1;
            do {
                assert( k == (*i).key);
                ++i;
                ++k;
            } while (i != end);
        }
        { // Open own scope to hide local variables.
            // Check generally correct parameter properties.
            CGAL::Assert_circulator_or_iterator(begin);
            CGAL::Assert_circulator_or_iterator(end);
            CGAL::Assert_is_at_least_forward_category(begin);
            CGAL::Assert_is_at_least_forward_category(end);
            typedef std::iterator_traits< Iterator >::value_type      VT;
            typedef std::iterator_traits< Iterator >::difference_type DT;
            assert(1==test_value_type(static_cast< VT* >(nullptr)));
            assert(1==test_distance_type(static_cast< DT* >(nullptr)));

            // Default constructor.
            Iterator z ;
            CGAL::Assert_circulator_or_iterator(z);
            z = end ;
            // Copy constructor.
            Iterator i = begin;

            // Check general support for circulators and iterators.
            assert( CGAL::is_empty_range( z, z));
            assert( ! CGAL::is_empty_range( i, end));

            int su = 0;
            int k  = 1;
            // Check general loop, pre-increment, dereference.
            if (! CGAL::is_empty_range( i, end)) {   // superfluous
                do {
                    assert( k == (*i).key);
                    su += (*i).key;
                    ++k;
                    Iterator j = ++i;
                    assert( i ==  j);
                    if ( i != end) {
                        assert( (*i).key == (*j).key);
                    }
                } while (i != end);  // Inequality and equality checked.
            }
            assert( i == end);  // Equality checked.
            assert( su == 15);

            // Assignment.
            i = begin;
            su = 0;
            k  = 1;
            // Loop with post increment.
            if (! CGAL::is_empty_range( i, end)) {   // superfluous
                do {
                    assert( k == (*i).key);
                    su += (*i).key;
                    ++k;
                    Iterator j = i++;
                    assert(  i !=  j);
                    if ( i != end) {
                        assert( (*i).key == (*j).key + 1);
                    }
                } while (i != end);
            }
            assert( i == end);
            assert( su == 15);
        }
        { // Open own scope to hide local variables.
            CGAL::Assert_is_at_least_bidirectional_category(begin);
            CGAL::Assert_is_at_least_bidirectional_category(end);
            // Loop backwards and pre-decrement.
            Iterator i = end;
            int su = 0;
            int k  = 5;
            do {
                Iterator j = --i;
                assert(  i ==  j);
                assert( (*i).key == (*j).key);
                assert( k == (*i).key);
                su += (*i).key;
                --k;
            } while (i != begin);
            assert( i == begin);
            assert( su == 15);

            // Assignment.
            i = end;
            su = 0;
            k  = 5;
            // Loop with post-decrement.
            do {
                Iterator j = i--;
                assert(  i !=  j);
                if ( j != end) {
                    assert( (*i).key == (*j).key - 1);
                }
                assert( k == (*i).key);
                su += (*i).key;
                --k;
            } while (i != begin);
            assert( i == begin);
            assert( su == 15);
        }
        CGAL::Assert_iterator( begin);
        CGAL::Assert_iterator( end);

        List l2 = l;
        const List& l1 = l2;
        typedef List::const_iterator ConstIterBase;
        typedef N_step_adaptor_derived<ConstIterBase,2> C_Iterator;
        C_Iterator c_begin(l1.begin());
        C_Iterator c_end(l1.end());
        Assert_bidirectional_category(c_begin);
        Assert_bidirectional_category(c_end);
        { // Open own scope to hide local variables.
            // Check generally correct parameter properties.
            CGAL::Assert_circulator_or_iterator(c_begin);
            CGAL::Assert_circulator_or_iterator(c_end);
            CGAL::Assert_is_at_least_forward_category(c_begin);
            CGAL::Assert_is_at_least_forward_category(c_end);
            typedef std::iterator_traits< C_Iterator >::value_type      VT;
            typedef std::iterator_traits< C_Iterator >::difference_type DT;
            assert(1==test_value_type(static_cast< VT* >(nullptr)));
            assert(1==test_distance_type(static_cast< DT* >(nullptr)));

            // Default constructor.
            C_Iterator z ;
            CGAL::Assert_circulator_or_iterator(z);
            z = c_end ;
            // Copy constructor.
            C_Iterator i = c_begin;

            // Check general support for circulators and iterators.
            assert( CGAL::is_empty_range( z, z));
            assert( ! CGAL::is_empty_range( i, c_end));

            int su = 0;
            int k  = 1;
            // Check general loop, pre-increment, dereference.
            if (! CGAL::is_empty_range( i, c_end)) {   // superfluous
                do {
                    assert( k == (*i).key);
                    su += (*i).key;
                    ++k;
                    C_Iterator j = ++i;
                    assert( i ==  j);
                    if ( i != c_end) {
                        assert( (*i).key == (*j).key);
                    }
                } while (i != c_end);  // Inequality and equality checked.
            }
            assert( i == c_end);  // Equality checked.
            assert( su == 15);

            // Assignment.
            i = c_begin;
            su = 0;
            k  = 1;
            // Loop with post increment.
            if (! CGAL::is_empty_range( i, c_end)) {   // superfluous
                do {
                    assert( k == (*i).key);
                    su += (*i).key;
                    ++k;
                    C_Iterator j = i++;
                    assert(  i !=  j);
                    if ( i != c_end) {
                        assert( (*i).key == (*j).key + 1);
                    }
                } while (i != c_end);
            }
            assert( i == c_end);
            assert( su == 15);
        }
        { // Open own scope to hide local variables.
            CGAL::Assert_is_at_least_bidirectional_category(c_begin);
            CGAL::Assert_is_at_least_bidirectional_category(c_end);
            // Loop backwards and pre-decrement.
            C_Iterator i = c_end;
            int su = 0;
            int k  = 5;
            do {
                C_Iterator j = --i;
                assert(  i ==  j);
                assert( (*i).key == (*j).key);
                assert( k == (*i).key);
                su += (*i).key;
                --k;
            } while (i != c_begin);
            assert( i == c_begin);
            assert( su == 15);

            // Assignment.
            i = c_end;
            su = 0;
            k  = 5;
            // Loop with post-decrement.
            do {
                C_Iterator j = i--;
                assert(  i !=  j);
                if ( j != c_end) {
                    assert( (*i).key == (*j).key - 1);
                }
                assert( k == (*i).key);
                su += (*i).key;
                --k;
            } while (i != c_begin);
            assert( i == c_begin);
            assert( su == 15);
        }
        CGAL::Assert_iterator( c_begin);
        CGAL::Assert_iterator( c_end);
        l.destroy();
        l2.destroy();
    }
    {
        typedef std::vector<int> Vector;
        Vector v;
        v.push_back(1);
        v.push_back(1);
        v.push_back(2);
        v.push_back(2);
        v.push_back(3);
        v.push_back(3);
        v.push_back(4);
        v.push_back(4);
        v.push_back(5);
        v.push_back(5);
        assert( v.size() == 10);
        typedef Vector::iterator IterBase;
        typedef Random_access_circulator_from_iterator<IterBase,int,
                std::size_t,std::ptrdiff_t> CircBase;
        typedef N_step_adaptor_derived<CircBase,2> Circulator;
        Circulator begin(CircBase( v.begin(),v.end()));
        Assert_random_access_category(begin);
        { // Open own scope to hide local variables.
            // Check generally correct parameter properties.
            CGAL::Assert_circulator_or_iterator(begin);
            CGAL::Assert_circulator_or_iterator(begin);
            CGAL::Assert_is_at_least_forward_category(begin);
            CGAL::Assert_is_at_least_forward_category(begin);
            typedef std::iterator_traits< Circulator >::value_type      VT;
            typedef std::iterator_traits< Circulator >::difference_type DT;
            assert(1==test_value_type(static_cast< VT* >(nullptr)));
            assert(1==test_distance_type(static_cast< DT* >(nullptr)));

            // Default constructor.
            Circulator z = Circulator();
            CGAL::Assert_circulator_or_iterator(z);
            // Copy constructor.
            Circulator i = begin;

            // Check general support for circulators and iterators.
            assert( CGAL::is_empty_range( z, z));
            assert( ! CGAL::is_empty_range( i, begin));

            int su = 0;
            int k  = 1;
            // Check general loop, pre-increment, dereference.
            if (! CGAL::is_empty_range( i, begin)) {   // superfluous
                do {
                    assert( k == (*i));
                    su += (*i);
                    ++k;
                    Circulator j = ++i;
                    assert( i ==  j);
                    if ( i != begin) {
                        assert( (*i) == (*j));
                    }
                } while (i != begin);  // Inequality and equality checked.
            }
            assert( i == begin);  // Equality checked.
            assert( su == 15);

            // Assignment.
            i = begin;
            su = 0;
            k  = 1;
            // Loop with post increment.
            if (! CGAL::is_empty_range( i, begin)) {   // superfluous
                do {
                    assert( k == (*i));
                    su += (*i);
                    ++k;
                    Circulator j = i++;
                    assert(  i !=  j);
                    if ( i != begin) {
                        assert( (*i) == (*j) + 1);
                    }
                } while (i != begin);
            }
            assert( i == begin);
            assert( su == 15);
        }
        { // Open own scope to hide local variables.
            // Change three elements and check post-/pre-increment.
            Circulator i = begin;
            (*i++) = 4;
            assert( 4 == (*begin));
            assert( 2 == (*i));
            (*i++) = 3;
            assert( 3 == (*i));
            (*++i) = 7;
            assert( 7 == (*i));

            // Check the setting and reset these elements
            // to their original values.
            i = begin;
            assert( 4 == (*i));
            (*i) = 1;
            i++;
            assert( 3 == (*i));
            (*i++) = 2;
            assert( 3 == (*i));
            i++;
            assert( 7 == (*i));
            (*i) = 4;

            // Check the resetting.
            i = begin;
            int k = 1;
            do {
                assert( k == (*i));
                ++i;
                ++k;
            } while (i != begin);
        }
        { // Open own scope to hide local variables.
            // Check generally correct parameter properties.
            CGAL::Assert_circulator_or_iterator(begin);
            CGAL::Assert_circulator_or_iterator(begin);
            CGAL::Assert_is_at_least_forward_category(begin);
            CGAL::Assert_is_at_least_forward_category(begin);
            typedef std::iterator_traits< Circulator >::value_type      VT;
            typedef std::iterator_traits< Circulator >::difference_type DT;
            assert(1==test_value_type(static_cast< VT* >(nullptr)));
            assert(1==test_distance_type(static_cast< DT* >(nullptr)));

            // Default constructor.
            Circulator z = Circulator();
            CGAL::Assert_circulator_or_iterator(z);
            // Copy constructor.
            Circulator i = begin;

            // Check general support for circulators and iterators.
            assert( CGAL::is_empty_range( z, z));
            assert( ! CGAL::is_empty_range( i, begin));

            int su = 0;
            int k  = 1;
            // Check general loop, pre-increment, dereference.
            if (! CGAL::is_empty_range( i, begin)) {   // superfluous
                do {
                    assert( k == (*i));
                    su += (*i);
                    ++k;
                    Circulator j = ++i;
                    assert( i ==  j);
                    if ( i != begin) {
                        assert( (*i) == (*j));
                    }
                } while (i != begin);  // Inequality and equality checked.
            }
            assert( i == begin);  // Equality checked.
            assert( su == 15);

            // Assignment.
            i = begin;
            su = 0;
            k  = 1;
            // Loop with post increment.
            if (! CGAL::is_empty_range( i, begin)) {   // superfluous
                do {
                    assert( k == (*i));
                    su += (*i);
                    ++k;
                    Circulator j = i++;
                    assert(  i !=  j);
                    if ( i != begin) {
                        assert( (*i) == (*j) + 1);
                    }
                } while (i != begin);
            }
            assert( i == begin);
            assert( su == 15);
        }
        { // Open own scope to hide local variables.
            CGAL::Assert_is_at_least_bidirectional_category(begin);
            CGAL::Assert_is_at_least_bidirectional_category(begin);
            // Loop backwards and pre-decrement.
            Circulator i = begin;
            int su = 0;
            int k  = 5;
            do {
                Circulator j = --i;
                assert(  i ==  j);
                assert( (*i) == (*j));
                assert( k == (*i));
                su += (*i);
                --k;
            } while (i != begin);
            assert( i == begin);
            assert( su == 15);

            // Assignment.
            i = begin;
            su = 0;
            k  = 5;
            // Loop with post-decrement.
            do {
                Circulator j = i--;
                assert(  i !=  j);
                if ( j != begin) {
                    assert( (*i) == (*j) - 1);
                }
                assert( k == (*i));
                su += (*i);
                --k;
            } while (i != begin);
            assert( i == begin);
            assert( su == 15);
        }
        { // Open own scope to hide local variables.
            // Check generally correct parameter properties.
            CGAL::Assert_circulator_or_iterator(begin);
            CGAL::Assert_circulator_or_iterator(begin);
            CGAL::Assert_is_at_least_forward_category(begin);
            CGAL::Assert_is_at_least_forward_category(begin);
            typedef std::iterator_traits< Circulator >::value_type      VT;
            typedef std::iterator_traits< Circulator >::difference_type DT;
            assert(1==test_value_type(static_cast< VT* >(nullptr)));
            assert(1==test_distance_type(static_cast< DT* >(nullptr)));

            // Default constructor.
            Circulator z = Circulator();
            CGAL::Assert_circulator_or_iterator(z);
            // Copy constructor.
            Circulator i = begin;

            // Check general support for circulators and iterators.
            assert( CGAL::is_empty_range( z, z));
            assert( ! CGAL::is_empty_range( i, begin));

            int su = 0;
            int k  = 1;
            // Check general loop, pre-increment, dereference.
            if (! CGAL::is_empty_range( i, begin)) {   // superfluous
                do {
                    assert( k == (*i));
                    su += (*i);
                    ++k;
                    Circulator j = ++i;
                    assert( i ==  j);
                    if ( i != begin) {
                        assert( (*i) == (*j));
                    }
                } while (i != begin);  // Inequality and equality checked.
            }
            assert( i == begin);  // Equality checked.
            assert( su == 15);

            // Assignment.
            i = begin;
            su = 0;
            k  = 1;
            // Loop with post increment.
            if (! CGAL::is_empty_range( i, begin)) {   // superfluous
                do {
                    assert( k == (*i));
                    su += (*i);
                    ++k;
                    Circulator j = i++;
                    assert(  i !=  j);
                    if ( i != begin) {
                        assert( (*i) == (*j) + 1);
                    }
                } while (i != begin);
            }
            assert( i == begin);
            assert( su == 15);
        }
        { // Open own scope to hide local variables.
            CGAL::Assert_is_at_least_bidirectional_category(begin);
            CGAL::Assert_is_at_least_bidirectional_category(begin);
            // Loop backwards and pre-decrement.
            Circulator i = begin;
            int su = 0;
            int k  = 5;
            do {
                Circulator j = --i;
                assert(  i ==  j);
                assert( (*i) == (*j));
                assert( k == (*i));
                su += (*i);
                --k;
            } while (i != begin);
            assert( i == begin);
            assert( su == 15);

            // Assignment.
            i = begin;
            su = 0;
            k  = 5;
            // Loop with post-decrement.
            do {
                Circulator j = i--;
                assert(  i !=  j);
                if ( j != begin) {
                    assert( (*i) == (*j) - 1);
                }
                assert( k == (*i));
                su += (*i);
                --k;
            } while (i != begin);
            assert( i == begin);
            assert( su == 15);
        }
        { // Open own scope to hide local variables.
            CGAL::Assert_is_at_least_random_access_category(begin);
            CGAL::Assert_is_at_least_random_access_category(begin);
            // Random access.
            int k;
            for( k = 0; k < 5; k++) {
                assert( 1+k == begin[k]);
            }

              int su = begin[0]
                     + begin[1]
                     + begin[2]
                     + begin[3]
                     + begin[4];
            assert( su == 15);

            // Jump around.
            Circulator i = begin;
            i += 3;
            assert( 4 == (*i));
            i -= 2;
            assert( 2 == (*i));
            i += 3;
            assert( 5 == (*i));
            i -= 4;
            assert( 1 == (*i));
            assert( i == begin);
            Circulator j = i + 3;
            assert( 4 == (*j));
            Circulator jj = j - 2;
            assert( 2 == (*jj));
            typedef std::ptrdiff_t PT;
            jj = PT(4) + jj;
            assert( jj == begin);
            Circulator ij = jj - 5;
            ij = jj - 5; // avoids warning with NDEBUG
            assert( ij == begin);

            // Difference test.
            assert( jj - i == 5  ||  jj - i == 0);
            assert( i + (j-i) == j);
            assert( (j-i) + i == j);
        }
        { // Open own scope to hide local variables.
            Circulator i = begin;
            i[2] = 18;
            i[4] = 9;
            i[3] = 12;
            assert( i[2] == 18);
            assert( i[4] == 9);
            assert( i[3] == 12);
            i[2] = 3;
            i[3] = 4;
            i[4] = 5;
            // Check the resetting.
            i = begin;
            int k = 1;
            do {
                assert( k == (*i));
                ++i;
                ++k;
            } while (i != begin);
        }
        { // Open own scope to hide local variables.
            CGAL::Assert_circulator( begin);
            CGAL::Assert_circulator( begin);

            // Check the local type parameters.
            Circulator::value_type      k1;
            k1 = 1;
            Circulator::reference       k2 = k1;
            (void)k2;
            assert( k2 == 1);
            Circulator::pointer         p1 = &k1;
            (void)p1;
            assert( (*p1) == 1);
            k1 = 3;
            assert( k1 == 3);
            assert( k2 == 3);
            assert( (*p1) == 3);
            k1 = 6;
            assert( k1 == 6);
            assert( k2 == 6);
            assert( (*p1) == 6);
            Circulator::size_type s = 5;
            assert( s == 5);
            Circulator::difference_type d = -5;
            assert( d == -5);

            // Check tests for empty data structures.
            Circulator z = Circulator();
            assert(   z == nullptr);
            assert( ! (z != nullptr));
            Circulator i = begin;
            assert( ! (i == nullptr));
            assert(   i != nullptr);
            assert( i == begin);
            assert( i == begin);
            // Do I reach myself.
            ++i;
            Circulator j = i;
            int k = 0;
            do {
                assert( k < 5);
                ++k;
                ++i;
            } while( i != j);
            assert( k == 5);
        }
        { // Open own scope to hide local variables.
            // Do I reach myself backwards.
            Circulator i = begin;
            ++i;
            Circulator j = i;
            int k = 0;
            do {
                assert( k < 5);
                ++k;
                --i;
            } while( i != j);
            assert( k == 5);
        }
        { // Open own scope to hide local variables.
            Circulator::difference_type d = begin - begin;
            assert( d == 0);
            d = begin - begin;
            assert( d == 0);
            Circulator i = begin + 1;
            assert( begin - i == 1 ||  begin - i == -1);
            assert( i - begin == 1 ||  i - begin == -1);
            // Check minimal circulator properties.
            i = i.min_circulator();
            Circulator j = i;
            assert( j - i == 0);
            j++;
            assert( j - i == 1);
            j++;
            assert( j - i == 2);
            j++;
            assert( j - i == 3);
            j++;
            assert( j - i == 4);
            j++;
            assert( j - i == 0);
        }

        Vector v2 = v;
        const Vector& v1 = v2;
        typedef Vector::const_iterator ConstIterBase;
        typedef Random_access_const_circulator_from_iterator<
                    ConstIterBase, int, std::size_t, std::ptrdiff_t>
                ConstCircBase;
        typedef N_step_adaptor_derived<ConstCircBase,2> C_Circulator;
        C_Circulator c_begin(ConstCircBase(v1.begin(),v1.end()));
        Assert_random_access_category(c_begin);
        { // Open own scope to hide local variables.
            // Check generally correct parameter properties.
            CGAL::Assert_circulator_or_iterator(c_begin);
            CGAL::Assert_circulator_or_iterator(c_begin);
            CGAL::Assert_is_at_least_forward_category(c_begin);
            CGAL::Assert_is_at_least_forward_category(c_begin);
            typedef std::iterator_traits< C_Circulator >::value_type      VT;
            typedef std::iterator_traits< C_Circulator >::difference_type DT;
            assert(1==test_value_type(static_cast< VT* >(nullptr)));
            assert(1==test_distance_type(static_cast< DT* >(nullptr)));

            // Default constructor.
            C_Circulator z = C_Circulator();
            CGAL::Assert_circulator_or_iterator(z);
            // Copy constructor.
            C_Circulator i = c_begin;

            // Check general support for circulators and iterators.
            assert( CGAL::is_empty_range( z, z));
            assert( ! CGAL::is_empty_range( i, c_begin));

            int su = 0;
            int k  = 1;
            // Check general loop, pre-increment, dereference.
            if (! CGAL::is_empty_range( i, c_begin)) {   // superfluous
                do {
                    assert( k == (*i));
                    su += (*i);
                    ++k;
                    C_Circulator j = ++i;
                    assert( i ==  j);
                    if ( i != c_begin) {
                        assert( (*i) == (*j));
                    }
                } while (i != c_begin);  // Inequality and equality checked.
            }
            assert( i == c_begin);  // Equality checked.
            assert( su == 15);

            // Assignment.
            i = c_begin;
            su = 0;
            k  = 1;
            // Loop with post increment.
            if (! CGAL::is_empty_range( i, c_begin)) {   // superfluous
                do {
                    assert( k == (*i));
                    su += (*i);
                    ++k;
                    C_Circulator j = i++;
                    assert(  i !=  j);
                    if ( i != c_begin) {
                        assert( (*i) == (*j) + 1);
                    }
                } while (i != c_begin);
            }
            assert( i == c_begin);
            assert( su == 15);
        }
        { // Open own scope to hide local variables.
            CGAL::Assert_is_at_least_bidirectional_category(c_begin);
            CGAL::Assert_is_at_least_bidirectional_category(c_begin);
            // Loop backwards and pre-decrement.
            C_Circulator i = c_begin;
            int su = 0;
            int k  = 5;
            do {
                C_Circulator j = --i;
                assert(  i ==  j);
                assert( (*i) == (*j));
                assert( k == (*i));
                su += (*i);
                --k;
            } while (i != c_begin);
            assert( i == c_begin);
            assert( su == 15);

            // Assignment.
            i = c_begin;
            su = 0;
            k  = 5;
            // Loop with post-decrement.
            do {
                C_Circulator j = i--;
                assert(  i !=  j);
                if ( j != c_begin) {
                    assert( (*i) == (*j) - 1);
                }
                assert( k == (*i));
                su += (*i);
                --k;
            } while (i != c_begin);
            assert( i == c_begin);
            assert( su == 15);
        }
        { // Open own scope to hide local variables.
            CGAL::Assert_is_at_least_random_access_category(c_begin);
            CGAL::Assert_is_at_least_random_access_category(c_begin);
            // Random access.
            int k;
            for( k = 0; k < 5; k++) {
                assert( 1+k == c_begin[k]);
            }

              int su = c_begin[0]
                     + c_begin[1]
                     + c_begin[2]
                     + c_begin[3]
                     + c_begin[4];
            assert( su == 15);

            // Jump around.
            C_Circulator i = c_begin;
            i += 3;
            assert( 4 == (*i));
            i -= 2;
            assert( 2 == (*i));
            i += 3;
            assert( 5 == (*i));
            i -= 4;
            assert( 1 == (*i));
            assert( i == c_begin);
            C_Circulator j = i + 3;
            assert( 4 == (*j));
            C_Circulator jj = j - 2;
            assert( 2 == (*jj));
            typedef std::ptrdiff_t PT;
            jj = PT(4) + jj;
            assert( jj == c_begin);
            C_Circulator ij = jj - 5;
            ij = jj - 5; // avoids warning with NDEBUG
            assert( ij == c_begin);

            // Difference test.
            assert( jj - i == 5  ||  jj - i == 0);
            assert( i + (j-i) == j);
            assert( (j-i) + i == j);
        }
        { // Open own scope to hide local variables.
            CGAL::Assert_circulator( c_begin);
            CGAL::Assert_circulator( c_begin);

            // Check the local type parameters.
            C_Circulator::value_type      k1;
            k1 = 1;
            C_Circulator::reference       k2 = k1;
            (void)k2;
            assert( k2 == 1);
            C_Circulator::pointer         p1 = &k1;
            (void)p1;
            assert( (*p1) == 1);
            k1 = 3;
            assert( k1 == 3);
            assert( k2 == 3);
            assert( (*p1) == 3);
            k1 = 6;
            assert( k1 == 6);
            assert( k2 == 6);
            assert( (*p1) == 6);
            C_Circulator::size_type s = 5;
            assert( s == 5);
            C_Circulator::difference_type d = -5;
            assert( d == -5);

            // Check tests for empty data structures.
            C_Circulator z = C_Circulator();
            assert(   z == nullptr);
            assert( ! (z != nullptr));
            C_Circulator i = c_begin;
            assert( ! (i == nullptr));
            assert(   i != nullptr);
            assert( i == c_begin);
            assert( i == c_begin);
            // Do I reach myself.
            ++i;
            C_Circulator j = i;
            int k = 0;
            do {
                assert( k < 5);
                ++k;
                ++i;
            } while( i != j);
            assert( k == 5);
        }
        { // Open own scope to hide local variables.
            // Do I reach myself backwards.
            C_Circulator i = c_begin;
            ++i;
            C_Circulator j = i;
            int k = 0;
            do {
                assert( k < 5);
                ++k;
                --i;
            } while( i != j);
            assert( k == 5);
        }
        { // Open own scope to hide local variables.
            C_Circulator::difference_type d = c_begin - c_begin;
            assert( d == 0);
            d = c_begin - c_begin;
            assert( d == 0);
            C_Circulator i = c_begin + 1;
            assert( c_begin - i == 1 ||  c_begin - i == -1);
            assert( i - c_begin == 1 ||  i - c_begin == -1);
            // Check minimal circulator properties.
            i = i.min_circulator();
            C_Circulator j = i;
            assert( j - i == 0);
            j++;
            assert( j - i == 1);
            j++;
            assert( j - i == 2);
            j++;
            assert( j - i == 3);
            j++;
            assert( j - i == 4);
            j++;
            assert( j - i == 0);
        }
    }
}
struct Filt1 {
  template < class I >
  bool operator()(const I& i) const { return *i <= 5; }
};

struct Filt2 {
  template < class I >
  bool operator()(const I& i) const { return i->b <= 5; }
};

struct Filtest {
  Filtest(int a) : b(a) {}
  int b;
};

void test_Filter_iterator()
{
  {
    typedef std::list< int > Cont;
    typedef Cont::iterator iterator;
    typedef Filter_iterator< iterator, Filt1 > FI1;
    FI1 f;

    Cont l;
    l.push_back(3); l.push_back(8); l.push_back(1); l.push_back(6);
    l.push_back(9); l.push_back(5); l.push_back(3); l.push_back(2);
    l.push_back(6); l.push_back(4); l.push_back(7); l.push_back(9);
    Filt1 fi;

    FI1 f1(l.end(), fi, l.begin());
    ++f1;
    f1++;
    assert(*f1 == 9);
    --f1;
    f1--;
    assert(
      6 ==
      std::distance(filter_iterator(l.end(), fi, l.begin()),
                    filter_iterator(l.end(), fi, l.end())));
  }
  {
    typedef std::list< Filtest > Cont;
    typedef Cont::iterator iterator;
    typedef Filter_iterator< iterator, Filt2 > FI1;

    Cont l;
    Filtest f1(3), f2(4), f3(5), f4(6);
    l.push_back(f1); l.push_back(f2); l.push_back(f3); l.push_back(f4);
    Filt2 fi;
    FI1 f(l.end(), fi, l.begin());
    // next line just to get rid of "unused variable f" warning
    if (f == filter_iterator(l.end(), fi, l.begin())) f1 = 3;
    assert( 1 ==
      std::distance(f, filter_iterator(l.end(), fi, l.end())));
  }
}
void test_Inverse_index() {
  {
    typedef std::list<std::size_t> List;
    List l;
    l.push_back( 1);
    l.push_back( 2);
    l.push_back( 3);
    l.push_back( 4);
    l.push_back( 5);
    assert( l.size() == 5);
    typedef List::iterator Iterator;
    Inverse_index<Iterator> index(l.begin(),l.end());
    l.push_back( 6);
    index.push_back( -- (l.end()));
    l.push_back( 7);
    index.push_back( -- (l.end()));
    assert( l.size() == 7);
    for ( Iterator i = l.begin(); i != l.end(); ++i) {
      assert( *i == index[i] + 1);
    }
  }
  {
    typedef std::vector<std::size_t> Vector;
    Vector v;
    v.reserve(7);
    v.push_back( 1);
    v.push_back( 2);
    v.push_back( 3);
    v.push_back( 4);
    v.push_back( 5);
    assert( v.size() == 5);
    typedef Vector::iterator Iterator;
    Inverse_index<Iterator> index(v.begin(),v.end());
    v.push_back( 6);
    index.push_back( v.end() - 1);
    v.push_back( 7);
    index.push_back( v.end() - 1);
    assert( v.size() == 7);
    for ( Iterator i = v.begin(); i != v.end(); ++i) {
      assert( *i == index[i] + 1);
    }
  }
  {
    typedef std::list<std::size_t> List;
    List l;
    l.push_back( 1);
    l.push_back( 2);
    l.push_back( 3);
    l.push_back( 4);
    l.push_back( 5);
    assert( l.size() == 5);
    typedef List::iterator Iterator;
    typedef Bidirectional_circulator_from_iterator<Iterator,
      std::size_t, std::size_t, std::ptrdiff_t> Circulator;
    Inverse_index<Circulator> index(
      Circulator( l.begin(),l.end()),
      Circulator( l.begin(),l.end()));
    l.push_back( 6);
    Circulator cc( l.begin(),l.end(), -- (l.end()));
    index.push_back( cc);
    l.push_back( 7);
    cc = Circulator( l.begin(),l.end(), -- (l.end()));
    index.push_back( cc);
    assert( l.size() == 7);
    Circulator c( l.begin(),l.end());
    Circulator d = c;
    do {
      assert( *c == index[c] + 1);
    } while ( ++c != d);
  }
}
void test_Random_access_adaptor() {
  {
    typedef std::list<std::size_t> List;
    List l;
    l.push_back( 1);
    l.push_back( 2);
    l.push_back( 3);
    l.push_back( 4);
    l.push_back( 5);
    assert( l.size() == 5);
    typedef List::iterator Iterator;
    Random_access_adaptor<Iterator> index(l.begin(),l.end());
    l.push_back( 6);
    index.push_back( -- (l.end()));
    l.push_back( 7);
    index.push_back( -- (l.end()));
    assert( l.size() == 7);
    for ( Iterator i = l.begin(); i != l.end(); ++i) {
      assert( *i == *(index[*i - 1]));
    }
  }
  {
    typedef std::vector<std::size_t> Vector;
    Vector v;
    v.reserve(10);
    v.push_back( 1);
    v.push_back( 2);
    v.push_back( 3);
    v.push_back( 4);
    v.push_back( 5);
    assert( v.size() == 5);
    typedef Vector::iterator Iterator;
    Random_access_adaptor<Iterator> index(v.begin(),v.end());
    v.push_back( 6);
    index.push_back( v.end() - 1);
    v.push_back( 7);
    index.push_back( v.end() - 1);
    assert( v.size() == 7);
    for ( Iterator i = v.begin(); i != v.end(); ++i) {
      assert( *i == *(index[*i - 1]));
    }
  }
  {
    typedef std::list<std::size_t> List;
    List l;
    l.push_back( 1);
    l.push_back( 2);
    l.push_back( 3);
    l.push_back( 4);
    l.push_back( 5);
    assert( l.size() == 5);
    typedef List::iterator Iterator;
    typedef Bidirectional_circulator_from_iterator<Iterator,
      std::size_t, std::size_t, std::ptrdiff_t> Circulator;
    Random_access_adaptor<Circulator> index(
      Circulator( l.begin(),l.end()),
      Circulator( l.begin(),l.end()));
    l.push_back( 6);
    Circulator cc( l.begin(),l.end(), -- (l.end()));
    index.push_back( cc);
    l.push_back( 7);
    cc = Circulator( l.begin(),l.end(), -- (l.end()));
    index.push_back( cc);
    assert( l.size() == 7);
    Circulator c( l.begin(),l.end());
    Circulator d = c;
    do {
      assert( *c == *(index[*c - 1]));
    } while ( ++c != d);
  }
  {
    typedef std::list<std::size_t> List;
    List l;
    l.push_back( 1);
    l.push_back( 2);
    l.push_back( 3);
    l.push_back( 4);
    l.push_back( 5);
    assert( l.size() == 5);
    typedef List::iterator Iterator;
    Random_access_value_adaptor<Iterator,std::size_t>
    index(l.begin(),l.end());
  l.push_back( 6);
  index.push_back( -- (l.end()));
  l.push_back( 7);
  index.push_back( -- (l.end()));
  assert( l.size() == 7);
  for ( Iterator i = l.begin(); i != l.end(); ++i) {
    assert( *i == index[*i - 1]);
  }
  }
}
void test_Emptyset_iterator()
{
  Emptyset_iterator e;
  Assert_output_category(e);
  Emptyset_iterator f(e);
  Emptyset_iterator g = f;
  *g++ = 2;
  *g = 3.0;
  ++g;
}
struct A {};
struct B {
  B() {}
  B(int) {}
  operator A() { return A(); }
};

void test_Oneset_iterator()
{
  A a;
  B b;
  Oneset_iterator<A> e(a);
  Assert_bidirectional_category(e);
  Oneset_iterator<A> f(e);
  Oneset_iterator<A> g = f;
  *g++ = a;
  *g = b;
  ++g;
  --g;
  *g-- = a;
}
void test_Const_oneset_iterator()
{
  Const_oneset_iterator<B> e(B(1));
  Assert_random_access_category(e);
  Const_oneset_iterator<B> f(e);
  Const_oneset_iterator<B> g = f;
  A a = *g++;
  ++g;
  a = *g--;
  --g;
}
struct My_to_int    { operator int() const { return 2; } };
struct My_to_float  { operator float() const { return 3.25; } };
struct My_to_double { operator double() const { return 4.25; } };
struct My_to_bool   { operator bool() const { return true; } };

void test_Triple()
{
  typedef CGAL::Triple<int,double,bool>                    T1;
  typedef CGAL::Triple<T1,double,bool>                     T2;
  typedef CGAL::Triple<T1,T2,bool>                         T3;
  typedef CGAL::Triple<My_to_int,My_to_double,My_to_bool>  T4;

  T1 x1;
  x1 = CGAL::make_triple(1, 1.5, true);

  T1 xx1 = CGAL::make_tuple(1, 1.5, true);
  assert(x1 == xx1);
  assert(xx1.get<0>() == 1);
  assert(xx1.get<1>() == 1.5);
  assert(xx1.get<2>() == true);

  My_to_int mti;
  My_to_double mtd;
  My_to_bool mtb;
  T1 xx(mti, mtd, mtb);
  assert(xx.first == 2);
  T4 x4(mti, mtd, mtb);
  x1 = x4;
  assert(x1.first == 2);
  T2 x2;
  T3 x3;
  T1 y1(2, 2.2, true);
  T2 y2(x1, 2.2, true);
  T3 y3(y1, y2, false);

  T1 z1 = CGAL::make_triple(2, 2.0, false);
  T3 z3 = CGAL::make_triple(z1, CGAL::make_triple(z1, 1.0, true), false);
  x2 = CGAL::make_triple(x1, 2.2, true);

  assert(z3 < y3);
  if (z3 < y3) z3 = y3;
  assert(x2 == y2);
  if (x2 == y2) x3 = z3;
}

void test_Quadruple()
{
  typedef CGAL::Quadruple<int,float,double,bool>                         T1;
  typedef CGAL::Quadruple<T1,float,double,bool>                          T2;
  typedef CGAL::Quadruple<T1,T2,double,bool>                             T3;
  typedef CGAL::Quadruple<T1,T2,T3,bool>                                 T4;
  typedef CGAL::Quadruple<My_to_int,My_to_float,My_to_double,My_to_bool> T5;

  T1 x1;
  x1 = CGAL::make_quadruple(1, 2.5f, 1.5, true);

  T1 xx1 = CGAL::make_tuple(1, 2.5f, 1.5, true);
  assert(x1 == xx1);
  assert(xx1.get<0>() == 1);
  assert(xx1.get<1>() == 2.5f);
  assert(xx1.get<2>() == 1.5);
  assert(xx1.get<3>() == true);

  My_to_int mti;
  My_to_float mtf;
  My_to_double mtd;
  My_to_bool mtb;
  T1 xx(mti, mtf, mtd, mtb);
  assert(xx.first == 2);
  T5 x5(mti, mtf, mtd, mtb);
  x1 = x5;
  assert(x1.first == 2);
  x1 = CGAL::make_quadruple(1, 2.5f, 1.5, true);
  T2 x2;
  T3 x3;
  T4 x4;
  T1 y1(2, 1.5f, 2.2, true);
  T2 y2(x1, 2.5f, 2.2, true);
  T3 y3(y1, y2, 2.2, false);
  T4 y4(y1, y2, y3, true);

  T1 z1 = CGAL::make_quadruple(2, 1.0f, 2.0, false);
  T4 z4 = CGAL::make_quadruple(z1, y2,
                               CGAL::make_quadruple(z1, y2, 2.0, true),
                               false);

  assert(z4 < y4);
  if (z4 < y4) x2 = make_quadruple(x1, 2.5f , 2.2, true);
  assert(x2 == y2);
  if (x2 == y2) x4.third = x3;
}

void test_tuple(){
  typedef std::tuple<> T0;
  typedef std::tuple<int,int> T1;
  typedef std::tuple<int,My_to_int,int,int> T2;

  CGAL_USE_TYPE(T0);
  CGAL_USE_TYPE(T2);
  CGAL_static_assertion( std::tuple_size<T0>::value == 0 );
  CGAL_static_assertion( std::tuple_size<T1>::value == 2 );
  CGAL_static_assertion( std::tuple_size<T2>::value == 4 );
  CGAL_static_assertion( (boost::is_same<std::tuple_element<1,T2>::type,My_to_int>::value) );

  T1 t1=std::make_tuple(1,2);
  T1 t1_2=std::make_tuple(1,2);

  assert(t1==t1_2); // test the equality operator

  // T2 t2 = T2();
  // assert( t2 == T2() );
  //
  // Do not test equality between default initialized tuples, because
  // GNU/g++ version 4.1.2 does not default-initialize correctly
  // std::tr1::tuple.

  // Test std::tie
  int i1=-1,i2=-1;
  std::tie(i1,i2)=t1;
  assert( std::get<0>(t1)==i1 );
  assert( std::get<1>(t1)==i2 );

  // Test std::get for a pair
  double d = 1;
  std::pair<int, double *> pair(-3, &d);
  const std::pair<int, double *> const_pair(2, &d);

  assert(std::get<0>(pair) == -3);
  assert(std::get<1>(pair) == &d);
  assert(std::get<0>(const_pair) == 2);
  assert(std::get<1>(const_pair) == &d);
}

void test_prev_next()
{
  std::vector<int> V;
  V.push_back(1);
  V.push_back(2);
  V.push_back(3);

  assert(std::next(std::next(V.begin())) == std::prev(V.end()));
}

void test_copy_n() {
  std::vector<int> V;
  for(int i = 0; i < 10; ++i)
    V.push_back(i);

  std::vector<int> V2(5);
  std::copy_n(V.begin(), 5, V2.begin());

  assert(std::equal(V2.begin(), V2.end(), V.begin()));
}

struct SP_struct{
  SP_struct(int k):i(k){}
  int i;
  bool operator==(SP_struct other) const{
    return other.i==i;
  }
};

struct Cmp_SP_struct{
  bool operator()(SP_struct s1, SP_struct s2) const
  {
    return s1.i<s2.i;
  }
};

void test_make_sorted_pair() {
  std::pair<int,int> p1(1,2);
  std::pair<int,int> p2(2,1);
  assert( CGAL::make_sorted_pair(1,2)==p1 );
  assert( CGAL::make_sorted_pair(2,1)==p1 );
  assert( CGAL::make_sorted_pair(1,2,std::greater<int>())==p2 );
  assert( CGAL::make_sorted_pair(2,1,std::greater<int>())==p2 );

  SP_struct s1(1);
  SP_struct s2(2);
  assert( std::make_pair(s1,s2) ==
          CGAL::make_sorted_pair(s2,s1,Cmp_SP_struct()) );
  assert( std::make_pair(s2,s1) !=
          CGAL::make_sorted_pair(s2,s1,Cmp_SP_struct()) );
  std::pair<SP_struct,SP_struct> p3(s1,s2),
  p4=CGAL::make_sorted_pair(s2,s1,Cmp_SP_struct());
  assert(p3==p4);
  int i=2;
  assert( CGAL::make_sorted_pair(1,i) == std::make_pair(1,i) );
  CGAL_static_assertion( (boost::is_same<
                          BOOST_TYPEOF(CGAL::make_sorted_pair<long>(1L,i)),
                          std::pair<long,long> >::value) );
  assert( (CGAL::make_sorted_pair<long>(i,1L) == std::pair<long,long>(1L,2L)) );

  CGAL_static_assertion( (boost::is_same<
                          BOOST_TYPEOF(CGAL::make_sorted_pair<double>(1,2L)),
                          std::pair<double,double> >::value) );
  CGAL_static_assertion( (boost::is_same<
                          BOOST_TYPEOF(CGAL::make_sorted_pair<int>(1,2L)),
                          std::pair<int,int> >::value) );
}

void test_result_of() {
  struct Result_functor
  {
    int operator()()
    {
      return 0;
    }

    float operator()(const int&)
    {
      return 0.0f;
    }
  };

  typedef CGAL::cpp11::result_of<Result_functor(void)>::type result_type;
  typedef CGAL::cpp11::result_of<Result_functor(int)>::type result_type_float;
  CGAL_USE_TYPE(result_type);
  CGAL_USE_TYPE(result_type_float);
  CGAL_static_assertion((boost::is_same<result_type, int>::value));
  CGAL_static_assertion((boost::is_same<result_type_float, float>::value));

}

int main() {
  init_global_data();
  test_Circulator_identity();
  test_Iterator_project();
  test_Iterator_transform();
  test_Circulator_project();
  test_Circulator_on_node();
  test_N_step_adaptor();
  test_N_step_adaptor_derived();
  test_Filter_iterator();
  test_Inverse_index();
  test_Random_access_adaptor();
  test_Emptyset_iterator();
  test_Oneset_iterator();
  test_Const_oneset_iterator();
  test_Triple();
  test_Quadruple();
  clean_global_data();
  test_tuple();
  test_prev_next();
  test_copy_n();
  test_make_sorted_pair();
  test_result_of();
  return 0;
}
// EOF //
