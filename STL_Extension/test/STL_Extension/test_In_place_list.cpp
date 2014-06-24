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
// file          : test_In_place_list.C
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
// Stl_Extensions: In place list.
// ============================================================================


#include <CGAL/basic.h>
#include <cstddef>
#include <cassert>
#include <iterator>
#include <list>
#include <vector>
#include <CGAL/circulator.h>  // Needed for iterator category assertions.
#include <CGAL/Circulator/Circulator_adapters.h>  // Needed for test data structures.
#include <CGAL/In_place_list.h>

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
    Node() : key(0), next(0), prev(0) { next = prev = this; }
    Node(int n) : key(n), next(0), prev(0) { next = prev = this; }
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
    CNode() : next_(0), prev_(0), key(0) { next_ = prev_ = this; }
    CNode( int n) : next_(0), prev_(0), key(n) { next_ = prev_ = this; }
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

using namespace CGAL;

struct item : public In_place_list_base<item> {
  int key;
  item() {}
  item( int i) : key(i) {}
  item( const item& i)
  : In_place_list_base<item>(i), key(i.key) {}
  bool operator== (const item& i) const { return key == i.key;}
  bool operator!= (const item& i) const { return key != i.key;}
  bool operator== (int i) const         { return key == i;}
  bool operator!= (int i) const         { return key != i;}
  bool operator<  (const item& i) const { return key < i.key;}
};
int test_value_type( item*)           { return 1;}

struct Cont : public CGAL::In_place_list_base<Cont>
{
  int i_;
  Cont() : i_(4711) {}
  Cont(int i) : i_(i) {};
};

struct ContLt {
  bool operator()(const Cont& c1, const Cont& c2) const
  { return (c1.i_)<(c2.i_); }
};

void test_In_place_list() {
  {
    typedef In_place_list<item,false> List;
    typedef List::iterator            Iterator;
    typedef List::const_iterator      Const_iterator;
    typedef List::size_type           Size;
    List l1;
    assert( l1.empty());
    assert( l1.size() == 0);
    assert( l1.max_size() == Size(-1));
    assert( l1 == l1);
    assert( l1.begin() == l1.end());
    List l2(l1);
    assert( l1 == l2);
    assert( l2.begin() == l2.end());
    List l;
    l.push_back( *new item(1));
    l.push_back( *new item(2));
    l.push_back( *new item(3));
    l.push_back( *new item(4));
    l.push_back( *new item(5));
    l1 = l;
    assert( ! l1.empty());
    assert( l1.size() == 5);
    assert( l1 == l);
    assert( l1 != l2);
    Iterator begin = l.begin();
    Iterator end   = l.end();
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
        assert(1==test_value_type(static_cast< VT* >(0)));
        assert(1==test_distance_type(static_cast< DT* >(0)));
    
        // Default constructor.
        Iterator z = Iterator();
        z = Iterator(); // avoids warning with NDEBUG
        CGAL::Assert_circulator_or_iterator(z);
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
        assert(1==test_value_type(static_cast< VT* >(0)));
        assert(1==test_distance_type(static_cast< DT* >(0)));
    
        // Default constructor.
        Iterator z = Iterator();
        z = Iterator(); // avoids warning with NDEBUG
        CGAL::Assert_circulator_or_iterator(z);
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

    List l4 = l;
    const List& l3 = l4;
    Const_iterator c_begin = l3.begin();
    Const_iterator c_end   = l3.end();
    Assert_bidirectional_category(c_begin);
    Assert_bidirectional_category(c_end);
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(c_begin);
        CGAL::Assert_circulator_or_iterator(c_end);
        CGAL::Assert_is_at_least_forward_category(c_begin);
        CGAL::Assert_is_at_least_forward_category(c_end);
        typedef std::iterator_traits< Const_iterator >::value_type      VT;
        typedef std::iterator_traits< Const_iterator >::difference_type DT;
        assert(1==test_value_type(static_cast< VT* >(0)));
        assert(1==test_distance_type(static_cast< DT* >(0)));
    
        // Default constructor.
        Const_iterator z = Const_iterator();
        z = Const_iterator(); // avoids warning with NDEBUG
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Const_iterator i = c_begin;
    
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
                Const_iterator j = ++i;
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
                Const_iterator j = i++;
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
        Const_iterator i = c_end;
        int su = 0;
        int k  = 5;
        do {
            Const_iterator j = --i;
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
            Const_iterator j = i--;
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
    l1.destroy();
    l2.destroy();
    l4.destroy();
  }{
    typedef In_place_list<item,true> List;
    typedef List::iterator       Iterator;
    typedef List::const_iterator Const_iterator;
    typedef List::size_type      Size;
    List l1;
    assert( l1.empty());
    assert( l1.size() == 0);
    assert( l1.max_size() == Size(-1));
    assert( l1 == l1);
    assert( l1.begin() == l1.end());
    List l2(l1);
    assert( l1 == l2);
    assert( l2.begin() == l2.end());
    List l;
    l.push_back( *new item(1));
    l.push_back( *new item(2));
    l.push_back( *new item(3));
    l.push_back( *new item(4));
    l.push_back( *new item(5));
    l1 = l;
    assert( ! l1.empty());
    assert( l1.size() == 5);
    assert( l1 == l);
    assert( l1 != l2);
    Iterator begin = l.begin();
    Iterator end   = l.end();
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
        assert(1==test_value_type(static_cast< VT* >(0)));
        assert(1==test_distance_type(static_cast< DT* >(0)));
    
        // Default constructor.
        Iterator z = Iterator();
        z = Iterator(); // avoids warning with NDEBUG
        CGAL::Assert_circulator_or_iterator(z);
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
        assert(1==test_value_type(static_cast< VT* >(0)));
        assert(1==test_distance_type(static_cast< DT* >(0)));
    
        // Default constructor.
        Iterator z = Iterator();
        z = Iterator(); // avoids warning with NDEBUG
        CGAL::Assert_circulator_or_iterator(z);
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

    const List l3( l);
    Const_iterator c_begin = l3.begin();
    Const_iterator c_end   = l3.end();
    Assert_bidirectional_category(c_begin);
    Assert_bidirectional_category(c_end);
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(c_begin);
        CGAL::Assert_circulator_or_iterator(c_end);
        CGAL::Assert_is_at_least_forward_category(c_begin);
        CGAL::Assert_is_at_least_forward_category(c_end);
        typedef std::iterator_traits< Const_iterator >::value_type      VT;
        typedef std::iterator_traits< Const_iterator >::difference_type DT;
        assert(1==test_value_type(static_cast< VT* >(0)));
        assert(1==test_distance_type(static_cast< DT* >(0)));
    
        // Default constructor.
        Const_iterator z = Const_iterator();
        z = Const_iterator(); // avoids warning with NDEBUG
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Const_iterator i = c_begin;
    
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
                Const_iterator j = ++i;
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
                Const_iterator j = i++;
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
        Const_iterator i = c_end;
        int su = 0;
        int k  = 5;
        do {
            Const_iterator j = --i;
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
            Const_iterator j = i--;
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
  }{
    // in_list_prog.C
    typedef In_place_list<item,true> List;
    typedef List::iterator       Iterator;
    List  l;
    item* p = new item(1);
    l.push_back( *p);
    l.push_back( *new item(2));
    l.push_front( *new item(3));
    l.push_front( *new item(4));
    l.push_front( *new item(2));
    Iterator i = l.begin();
    ++i;
    l.insert(i,*new item(5));
    l.insert(p,*new item(5));
    int a[7] = {2,5,4,3,5,1,2};
    a[0] = 2; // avoids warning with NDEBUG
    assert( std::equal( l.begin(), l.end(), a));
    l.sort();
    l.unique();
    int b[5] = {1,2,3,4,5};
    b[0] = 1; // avoids warning with NDEBUG
    assert( l.size() == 5);
    assert( std::equal( l.begin(), l.end(), b));
  }

  {
    typedef CGAL::In_place_list<Cont,false> ContList;
    typedef CGAL::In_place_list<Cont,false>::iterator Iterator;
    typedef CGAL::In_place_list<Cont,false>::const_iterator Const_iterator;
  
    ContList L;
    L.push_back(* new Cont(3));
    L.push_back(* new Cont(5));
    L.push_back(* new Cont(-2));
    L.push_back(* new Cont(-2));
    L.push_back(* new Cont(-7));
    L.push_back(* new Cont(11));
    L.sort(ContLt());
    Iterator it1 = L.begin(), it2 = it1;
    for (++it2; it2 != L.end(); it1=it2, ++it2) {
      assert( (*it1).i_ <= (*it2).i_ );
    }

    // test remove_const
    Const_iterator cit=L.begin();
    Iterator it=cit.remove_const();

    CGAL_USE(cit);
    CGAL_USE(it);

    L.destroy();
  }

}

int main() {
  init_global_data();
  test_In_place_list();
  clean_global_data();
  return 0;
}
// EOF //
