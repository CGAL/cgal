// ============================================================================
//
// Copyright (c) 1997, 1998, 1999, 2000 The CGAL Consortium
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
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//                 Lutz Kettner <kettner@cs.unc.edu>
//
// maintainer    : Michael Hoffmann <hoffmann@inf.ethz.ch>
//
// Stl_Extensions: In place list.
// ============================================================================


#include <CGAL/basic.h>
#include <cstddef>
#include <iterator>
#include <list>
#include <vector>
#include <CGAL/circulator.h>  // Needed for iterator category assertions.
#include <CGAL/circulator_impl.h>  // Needed for test data structures.
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
    Node() : key(0), next(this), prev( this) {}
    Node( int n) : key(n), next(this), prev( this) {}
    Node( Node* nx_, Node* pv_, int n)
        : key(n), next(nx_), prev( pv_) {}
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
    CGAL_assertion( n > 0);
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
    CNode() : next_(this), prev_( this), key(0) {}
    CNode( int n) : next_(this), prev_( this), key(n) {}
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
    CGAL_assertion( n > 0);
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

void test_In_place_list() {
  {
    typedef In_place_list<item,false> List;
    typedef List::iterator            Iterator;
    typedef List::const_iterator      Const_iterator;
    typedef List::size_type           Size;
    List l1;
    CGAL_assertion( l1.empty());
    CGAL_assertion( l1.size() == 0);
    CGAL_assertion( l1.max_size() == Size(-1));
    CGAL_assertion( l1 == l1);
    CGAL_assertion( l1.begin() == l1.end());
    List l2(l1);
    CGAL_assertion( l1 == l2);
    CGAL_assertion( l2.begin() == l2.end());
    List l;
    l.push_back( *new item(1));
    l.push_back( *new item(2));
    l.push_back( *new item(3));
    l.push_back( *new item(4));
    l.push_back( *new item(5));
    l1 = l;
    CGAL_assertion( ! l1.empty());
    CGAL_assertion( l1.size() == 5);
    CGAL_assertion( l1 == l);
    CGAL_assertion( l1 != l2);
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
        CGAL_assertion(1==test_value_type(static_cast< VT* >(0)));
        CGAL_assertion(1==test_distance_type(static_cast< DT* >(0)));
    
        // Default constructor.
        Iterator z = Iterator();
        z = Iterator(); // avoids warning with NDEBUG
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Iterator i = begin;
    
        // Check general support for circulators and iterators.
        CGAL_assertion( CGAL::is_empty_range( z, z));
        CGAL_assertion( ! CGAL::is_empty_range( i, end));
    
        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, end)) {   // superfluous
            do {
                CGAL_assertion( k == (*i).key);
                su += (*i).key;
                ++k;
                CGAL_assertion_code( Iterator j =) ++i;
                CGAL_assertion( i ==  j);
                if ( i != end) {
                    CGAL_assertion( (*i).key == (*j).key);
                }
            } while (i != end);  // Inequality and equality checked.
        }
        CGAL_assertion( i == end);  // Equality checked.
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, end)) {   // superfluous
            do {
                CGAL_assertion( k == (*i).key);
                su += (*i).key;
                ++k;
                CGAL_assertion_code( Iterator j =) i++;
                CGAL_assertion(  i !=  j);
                if ( i != end) {
                    CGAL_assertion( (*i).key == (*j).key + 1);
                }
            } while (i != end);
        }
        CGAL_assertion( i == end);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        // Change three elements and check post-/pre-increment.
        Iterator i = begin;
        (*i++).key = 4;
        CGAL_assertion( 4 == (*begin).key);
        CGAL_assertion( 2 == (*i).key);
        (*i++).key = 3;
        CGAL_assertion( 3 == (*i).key);
        (*++i).key = 7;
        CGAL_assertion( 7 == (*i).key);
    
        // Check the setting and reset these elements
        // to their original values.
        i = begin;
        CGAL_assertion( 4 == (*i).key);
        (*i).key = 1;
        i++;
        CGAL_assertion( 3 == (*i).key);
        (*i++).key = 2;
        CGAL_assertion( 3 == (*i).key);
        i++;
        CGAL_assertion( 7 == (*i).key);
        (*i).key = 4;
    
        // Check the resetting.
        i = begin;
        int k = 1;
        do {
            CGAL_assertion( k == (*i).key);
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
        CGAL_assertion(1==test_value_type(static_cast< VT* >(0)));
        CGAL_assertion(1==test_distance_type(static_cast< DT* >(0)));
    
        // Default constructor.
        Iterator z = Iterator();
        z = Iterator(); // avoids warning with NDEBUG
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Iterator i = begin;
    
        // Check general support for circulators and iterators.
        CGAL_assertion( CGAL::is_empty_range( z, z));
        CGAL_assertion( ! CGAL::is_empty_range( i, end));
    
        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, end)) {   // superfluous
            do {
                CGAL_assertion( k == (*i).key);
                su += (*i).key;
                ++k;
                CGAL_assertion_code( Iterator j =) ++i;
                CGAL_assertion( i ==  j);
                if ( i != end) {
                    CGAL_assertion( (*i).key == (*j).key);
                }
            } while (i != end);  // Inequality and equality checked.
        }
        CGAL_assertion( i == end);  // Equality checked.
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, end)) {   // superfluous
            do {
                CGAL_assertion( k == (*i).key);
                su += (*i).key;
                ++k;
                CGAL_assertion_code( Iterator j =) i++;
                CGAL_assertion(  i !=  j);
                if ( i != end) {
                    CGAL_assertion( (*i).key == (*j).key + 1);
                }
            } while (i != end);
        }
        CGAL_assertion( i == end);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_bidirectional_category(begin);
        CGAL::Assert_is_at_least_bidirectional_category(end);
        // Loop backwards and pre-decrement.
        Iterator i = end;
        int su = 0;
        int k  = 5;
        do {
            CGAL_assertion_code(Iterator j =) --i;
            CGAL_assertion(  i ==  j);
            CGAL_assertion( (*i).key == (*j).key);
            CGAL_assertion( k == (*i).key);
            su += (*i).key;
            --k;
        } while (i != begin);
        CGAL_assertion( i == begin);
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = end;
        su = 0;
        k  = 5;
        // Loop with post-decrement.
        do {
            Iterator j = i--;
            CGAL_assertion(  i !=  j);
            if ( j != end) {
                CGAL_assertion( (*i).key == (*j).key - 1);
            }
            CGAL_assertion( k == (*i).key);
            su += (*i).key;
            --k;
        } while (i != begin);
        CGAL_assertion( i == begin);
        CGAL_assertion( su == 15);
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
        CGAL_assertion(1==test_value_type(static_cast< VT* >(0)));
        CGAL_assertion(1==test_distance_type(static_cast< DT* >(0)));
    
        // Default constructor.
        Const_iterator z = Const_iterator();
        z = Const_iterator(); // avoids warning with NDEBUG
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Const_iterator i = c_begin;
    
        // Check general support for circulators and iterators.
        CGAL_assertion( CGAL::is_empty_range( z, z));
        CGAL_assertion( ! CGAL::is_empty_range( i, c_end));
    
        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, c_end)) {   // superfluous
            do {
                CGAL_assertion( k == (*i).key);
                su += (*i).key;
                ++k;
                CGAL_assertion_code( Const_iterator j =) ++i;
                CGAL_assertion( i ==  j);
                if ( i != c_end) {
                    CGAL_assertion( (*i).key == (*j).key);
                }
            } while (i != c_end);  // Inequality and equality checked.
        }
        CGAL_assertion( i == c_end);  // Equality checked.
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = c_begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, c_end)) {   // superfluous
            do {
                CGAL_assertion( k == (*i).key);
                su += (*i).key;
                ++k;
                CGAL_assertion_code( Const_iterator j =) i++;
                CGAL_assertion(  i !=  j);
                if ( i != c_end) {
                    CGAL_assertion( (*i).key == (*j).key + 1);
                }
            } while (i != c_end);
        }
        CGAL_assertion( i == c_end);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_bidirectional_category(c_begin);
        CGAL::Assert_is_at_least_bidirectional_category(c_end);
        // Loop backwards and pre-decrement.
        Const_iterator i = c_end;
        int su = 0;
        int k  = 5;
        do {
            CGAL_assertion_code(Const_iterator j =) --i;
            CGAL_assertion(  i ==  j);
            CGAL_assertion( (*i).key == (*j).key);
            CGAL_assertion( k == (*i).key);
            su += (*i).key;
            --k;
        } while (i != c_begin);
        CGAL_assertion( i == c_begin);
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = c_end;
        su = 0;
        k  = 5;
        // Loop with post-decrement.
        do {
            Const_iterator j = i--;
            CGAL_assertion(  i !=  j);
            if ( j != c_end) {
                CGAL_assertion( (*i).key == (*j).key - 1);
            }
            CGAL_assertion( k == (*i).key);
            su += (*i).key;
            --k;
        } while (i != c_begin);
        CGAL_assertion( i == c_begin);
        CGAL_assertion( su == 15);
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
    CGAL_assertion( l1.empty());
    CGAL_assertion( l1.size() == 0);
    CGAL_assertion( l1.max_size() == Size(-1));
    CGAL_assertion( l1 == l1);
    CGAL_assertion( l1.begin() == l1.end());
    List l2(l1);
    CGAL_assertion( l1 == l2);
    CGAL_assertion( l2.begin() == l2.end());
    List l;
    l.push_back( *new item(1));
    l.push_back( *new item(2));
    l.push_back( *new item(3));
    l.push_back( *new item(4));
    l.push_back( *new item(5));
    l1 = l;
    CGAL_assertion( ! l1.empty());
    CGAL_assertion( l1.size() == 5);
    CGAL_assertion( l1 == l);
    CGAL_assertion( l1 != l2);
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
        CGAL_assertion(1==test_value_type(static_cast< VT* >(0)));
        CGAL_assertion(1==test_distance_type(static_cast< DT* >(0)));
    
        // Default constructor.
        Iterator z = Iterator();
        z = Iterator(); // avoids warning with NDEBUG
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Iterator i = begin;
    
        // Check general support for circulators and iterators.
        CGAL_assertion( CGAL::is_empty_range( z, z));
        CGAL_assertion( ! CGAL::is_empty_range( i, end));
    
        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, end)) {   // superfluous
            do {
                CGAL_assertion( k == (*i).key);
                su += (*i).key;
                ++k;
                CGAL_assertion_code( Iterator j =) ++i;
                CGAL_assertion( i ==  j);
                if ( i != end) {
                    CGAL_assertion( (*i).key == (*j).key);
                }
            } while (i != end);  // Inequality and equality checked.
        }
        CGAL_assertion( i == end);  // Equality checked.
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, end)) {   // superfluous
            do {
                CGAL_assertion( k == (*i).key);
                su += (*i).key;
                ++k;
                CGAL_assertion_code( Iterator j =) i++;
                CGAL_assertion(  i !=  j);
                if ( i != end) {
                    CGAL_assertion( (*i).key == (*j).key + 1);
                }
            } while (i != end);
        }
        CGAL_assertion( i == end);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        // Change three elements and check post-/pre-increment.
        Iterator i = begin;
        (*i++).key = 4;
        CGAL_assertion( 4 == (*begin).key);
        CGAL_assertion( 2 == (*i).key);
        (*i++).key = 3;
        CGAL_assertion( 3 == (*i).key);
        (*++i).key = 7;
        CGAL_assertion( 7 == (*i).key);
    
        // Check the setting and reset these elements
        // to their original values.
        i = begin;
        CGAL_assertion( 4 == (*i).key);
        (*i).key = 1;
        i++;
        CGAL_assertion( 3 == (*i).key);
        (*i++).key = 2;
        CGAL_assertion( 3 == (*i).key);
        i++;
        CGAL_assertion( 7 == (*i).key);
        (*i).key = 4;
    
        // Check the resetting.
        i = begin;
        int k = 1;
        do {
            CGAL_assertion( k == (*i).key);
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
        CGAL_assertion(1==test_value_type(static_cast< VT* >(0)));
        CGAL_assertion(1==test_distance_type(static_cast< DT* >(0)));
    
        // Default constructor.
        Iterator z = Iterator();
        z = Iterator(); // avoids warning with NDEBUG
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Iterator i = begin;
    
        // Check general support for circulators and iterators.
        CGAL_assertion( CGAL::is_empty_range( z, z));
        CGAL_assertion( ! CGAL::is_empty_range( i, end));
    
        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, end)) {   // superfluous
            do {
                CGAL_assertion( k == (*i).key);
                su += (*i).key;
                ++k;
                CGAL_assertion_code( Iterator j =) ++i;
                CGAL_assertion( i ==  j);
                if ( i != end) {
                    CGAL_assertion( (*i).key == (*j).key);
                }
            } while (i != end);  // Inequality and equality checked.
        }
        CGAL_assertion( i == end);  // Equality checked.
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, end)) {   // superfluous
            do {
                CGAL_assertion( k == (*i).key);
                su += (*i).key;
                ++k;
                CGAL_assertion_code( Iterator j =) i++;
                CGAL_assertion(  i !=  j);
                if ( i != end) {
                    CGAL_assertion( (*i).key == (*j).key + 1);
                }
            } while (i != end);
        }
        CGAL_assertion( i == end);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_bidirectional_category(begin);
        CGAL::Assert_is_at_least_bidirectional_category(end);
        // Loop backwards and pre-decrement.
        Iterator i = end;
        int su = 0;
        int k  = 5;
        do {
            CGAL_assertion_code(Iterator j =) --i;
            CGAL_assertion(  i ==  j);
            CGAL_assertion( (*i).key == (*j).key);
            CGAL_assertion( k == (*i).key);
            su += (*i).key;
            --k;
        } while (i != begin);
        CGAL_assertion( i == begin);
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = end;
        su = 0;
        k  = 5;
        // Loop with post-decrement.
        do {
            Iterator j = i--;
            CGAL_assertion(  i !=  j);
            if ( j != end) {
                CGAL_assertion( (*i).key == (*j).key - 1);
            }
            CGAL_assertion( k == (*i).key);
            su += (*i).key;
            --k;
        } while (i != begin);
        CGAL_assertion( i == begin);
        CGAL_assertion( su == 15);
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
        CGAL_assertion(1==test_value_type(static_cast< VT* >(0)));
        CGAL_assertion(1==test_distance_type(static_cast< DT* >(0)));
    
        // Default constructor.
        Const_iterator z = Const_iterator();
        z = Const_iterator(); // avoids warning with NDEBUG
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Const_iterator i = c_begin;
    
        // Check general support for circulators and iterators.
        CGAL_assertion( CGAL::is_empty_range( z, z));
        CGAL_assertion( ! CGAL::is_empty_range( i, c_end));
    
        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, c_end)) {   // superfluous
            do {
                CGAL_assertion( k == (*i).key);
                su += (*i).key;
                ++k;
                CGAL_assertion_code( Const_iterator j =) ++i;
                CGAL_assertion( i ==  j);
                if ( i != c_end) {
                    CGAL_assertion( (*i).key == (*j).key);
                }
            } while (i != c_end);  // Inequality and equality checked.
        }
        CGAL_assertion( i == c_end);  // Equality checked.
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = c_begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, c_end)) {   // superfluous
            do {
                CGAL_assertion( k == (*i).key);
                su += (*i).key;
                ++k;
                CGAL_assertion_code( Const_iterator j =) i++;
                CGAL_assertion(  i !=  j);
                if ( i != c_end) {
                    CGAL_assertion( (*i).key == (*j).key + 1);
                }
            } while (i != c_end);
        }
        CGAL_assertion( i == c_end);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_bidirectional_category(c_begin);
        CGAL::Assert_is_at_least_bidirectional_category(c_end);
        // Loop backwards and pre-decrement.
        Const_iterator i = c_end;
        int su = 0;
        int k  = 5;
        do {
            CGAL_assertion_code(Const_iterator j =) --i;
            CGAL_assertion(  i ==  j);
            CGAL_assertion( (*i).key == (*j).key);
            CGAL_assertion( k == (*i).key);
            su += (*i).key;
            --k;
        } while (i != c_begin);
        CGAL_assertion( i == c_begin);
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = c_end;
        su = 0;
        k  = 5;
        // Loop with post-decrement.
        do {
            Const_iterator j = i--;
            CGAL_assertion(  i !=  j);
            if ( j != c_end) {
                CGAL_assertion( (*i).key == (*j).key - 1);
            }
            CGAL_assertion( k == (*i).key);
            su += (*i).key;
            --k;
        } while (i != c_begin);
        CGAL_assertion( i == c_begin);
        CGAL_assertion( su == 15);
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
    CGAL_assertion( std::equal( l.begin(), l.end(), a));
    l.sort();
    l.unique();
    int b[5] = {1,2,3,4,5};
    b[0] = 1; // avoids warning with NDEBUG
    CGAL_assertion( l.size() == 5);
    CGAL_assertion( std::equal( l.begin(), l.end(), b));
  }
}

int main() {
  init_global_data();
  test_In_place_list();
  clean_global_data();
  return 0;
}
// EOF //
