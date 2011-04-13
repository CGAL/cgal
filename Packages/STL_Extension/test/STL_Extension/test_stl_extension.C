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
// file          : test_stl_extension.C
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
// Stl_Extensions: Iterator Adaptor.
// ============================================================================


#include <CGAL/basic.h>
#include <cstddef>
#include <list>
#include <vector>
#include <CGAL/In_place_list.h>
#include <CGAL/Iterator_identity.h>
#include <CGAL/Circulator_identity.h>
#include <CGAL/Iterator_project.h>
#include <CGAL/function_objects.h>
#include <CGAL/Circulator_project.h>
#include <CGAL/Circulator_on_node.h>
#include <CGAL/N_step_adaptor.h>
#include <CGAL/N_step_adaptor_derived.h>
#include <CGAL/function_objects.h>
#include <CGAL/Inverse_index.h>
#include <CGAL/Random_access_adaptor.h>
#include <CGAL/Random_access_value_adaptor.h>

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
int test_value_type( item*)           { return 1;}

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
void test_Iterator_identity() {
  {
    typedef In_place_list<item,false> List;
    List l;
    l.push_back( *new item(1));
    l.push_back( *new item(2));
    l.push_back( *new item(3));
    l.push_back( *new item(4));
    l.push_back( *new item(5));
    CGAL_assertion( l.size() == 5);
    typedef List::iterator IterBase;
    typedef Iterator_identity<IterBase,item&,item*,item,std::ptrdiff_t,
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

    List l2 = l;
    const List& l1 = l2;
    typedef List::const_iterator ConstIterBase;
    typedef Iterator_identity< ConstIterBase,const item&,
      const item*, item, std::ptrdiff_t,
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
        CGAL_assertion(1==test_value_type(static_cast< VT* >(0)));
        CGAL_assertion(1==test_distance_type(static_cast< DT* >(0)));
    
        // Default constructor.
        C_Iterator z = C_Iterator();
        z = C_Iterator(); // avoids warning with NDEBUG
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        C_Iterator i = c_begin;
    
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
                CGAL_assertion_code( C_Iterator j =) ++i;
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
                CGAL_assertion_code( C_Iterator j =) i++;
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
        C_Iterator i = c_end;
        int su = 0;
        int k  = 5;
        do {
            CGAL_assertion_code(C_Iterator j =) --i;
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
            C_Iterator j = i--;
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
    CGAL_assertion( v.size() == 5);
    typedef Vector::iterator IterBase;
    typedef Iterator_identity<IterBase,int&,int*,int,std::ptrdiff_t,
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
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                CGAL_assertion_code( Iterator j =) ++i;
                CGAL_assertion( i ==  j);
                if ( i != end) {
                    CGAL_assertion( (*i) == (*j));
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
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                CGAL_assertion_code( Iterator j =) i++;
                CGAL_assertion(  i !=  j);
                if ( i != end) {
                    CGAL_assertion( (*i) == (*j) + 1);
                }
            } while (i != end);
        }
        CGAL_assertion( i == end);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        // Change three elements and check post-/pre-increment.
        Iterator i = begin;
        (*i++) = 4;
        CGAL_assertion( 4 == (*begin));
        CGAL_assertion( 2 == (*i));
        (*i++) = 3;
        CGAL_assertion( 3 == (*i));
        (*++i) = 7;
        CGAL_assertion( 7 == (*i));
    
        // Check the setting and reset these elements
        // to their original values.
        i = begin;
        CGAL_assertion( 4 == (*i));
        (*i) = 1;
        i++;
        CGAL_assertion( 3 == (*i));
        (*i++) = 2;
        CGAL_assertion( 3 == (*i));
        i++;
        CGAL_assertion( 7 == (*i));
        (*i) = 4;
    
        // Check the resetting.
        i = begin;
        int k = 1;
        do {
            CGAL_assertion( k == (*i));
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
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                CGAL_assertion_code( Iterator j =) ++i;
                CGAL_assertion( i ==  j);
                if ( i != end) {
                    CGAL_assertion( (*i) == (*j));
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
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                CGAL_assertion_code( Iterator j =) i++;
                CGAL_assertion(  i !=  j);
                if ( i != end) {
                    CGAL_assertion( (*i) == (*j) + 1);
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
            CGAL_assertion( (*i) == (*j));
            CGAL_assertion( k == (*i));
            su += (*i);
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
                CGAL_assertion( (*i) == (*j) - 1);
            }
            CGAL_assertion( k == (*i));
            su += (*i);
            --k;
        } while (i != begin);
        CGAL_assertion( i == begin);
        CGAL_assertion( su == 15);
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
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                CGAL_assertion_code( Iterator j =) ++i;
                CGAL_assertion( i ==  j);
                if ( i != end) {
                    CGAL_assertion( (*i) == (*j));
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
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                CGAL_assertion_code( Iterator j =) i++;
                CGAL_assertion(  i !=  j);
                if ( i != end) {
                    CGAL_assertion( (*i) == (*j) + 1);
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
            CGAL_assertion( (*i) == (*j));
            CGAL_assertion( k == (*i));
            su += (*i);
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
                CGAL_assertion( (*i) == (*j) - 1);
            }
            CGAL_assertion( k == (*i));
            su += (*i);
            --k;
        } while (i != begin);
        CGAL_assertion( i == begin);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_random_access_category(begin);
        CGAL::Assert_is_at_least_random_access_category(end);
        // Random access.
        int k;
        for( k = 0; k < 5; k++) {
            CGAL_assertion( 1+k == begin[k]);
        }
        CGAL_assertion_code(
          int su = begin[0]
                 + begin[1]
                 + begin[2]
                 + begin[3]
                 + begin[4];)
        CGAL_assertion( su == 15);
    
        // Jump around.
        Iterator i = begin;
        i += 3;
        CGAL_assertion( 4 == (*i));
        i -= 2;
        CGAL_assertion( 2 == (*i));
        i += 3;
        CGAL_assertion( 5 == (*i));
        i -= 4;
        CGAL_assertion( 1 == (*i));
        CGAL_assertion( i == begin);
        Iterator j = i + 3;
        CGAL_assertion( 4 == (*j));
        Iterator jj = j - 2;
        CGAL_assertion( 2 == (*jj));
        jj = 4 + jj;
        CGAL_assertion( jj == end);
        Iterator ij = jj - 5;
        ij = jj - 5; // avoids warning with NDEBUG
        CGAL_assertion( ij == begin);
    
        // Difference test.
        CGAL_assertion( jj - i == 5  ||  jj - i == 0);
        CGAL_assertion( i + (j-i) == j);
        CGAL_assertion( (j-i) + i == j);
    }
    { // Open own scope to hide local variables.
        Iterator i = begin;
        i[2] = 18;
        i[4] = 9;
        i[3] = 12;
        CGAL_assertion( i[2] == 18);
        CGAL_assertion( i[4] == 9);
        CGAL_assertion( i[3] == 12);
        i[2] = 3;
        i[3] = 4;
        i[4] = 5;
        // Check the resetting.
        i = begin;
        int k = 1;
        do {
            CGAL_assertion( k == (*i));
            ++i;
            ++k;
        } while (i != end);
    }
    CGAL::Assert_iterator( begin);
    CGAL::Assert_iterator( end);
    { // Open own scope to hide local variables.
        CGAL_assertion( end - begin ==  5);
        CGAL_assertion( begin - end == -5);
        // Relational operator.
        Iterator i = begin;
        ++i;
        Iterator j = i;
        ++j;
        CGAL_assertion( begin < i);
        CGAL_assertion( i < j);
        CGAL_assertion( j < end);
        CGAL_assertion( j > i);
        CGAL_assertion( i <= j);
        CGAL_assertion( j >= i);
        CGAL_assertion( i <= i);
        CGAL_assertion( i >= i);
    
        CGAL_assertion( !( i >= j));
        CGAL_assertion( !( j <= i));
        CGAL_assertion( !( i > j));
        CGAL_assertion( !( j < i));
        CGAL_assertion( !( i > i));
        CGAL_assertion( !( i < i));
    }

    Vector v2 = v;
    const Vector& v1 = v2;
    typedef Vector::const_iterator ConstIterBase;
    typedef Iterator_identity< ConstIterBase,const int&,
      const int*, int, std::ptrdiff_t,
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
        CGAL_assertion(1==test_value_type(static_cast< VT* >(0)));
        CGAL_assertion(1==test_distance_type(static_cast< DT* >(0)));
    
        // Default constructor.
        C_Iterator z = C_Iterator();
        z = C_Iterator(); // avoids warning with NDEBUG
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        C_Iterator i = c_begin;
    
        // Check general support for circulators and iterators.
        CGAL_assertion( CGAL::is_empty_range( z, z));
        CGAL_assertion( ! CGAL::is_empty_range( i, c_end));
    
        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, c_end)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                CGAL_assertion_code( C_Iterator j =) ++i;
                CGAL_assertion( i ==  j);
                if ( i != c_end) {
                    CGAL_assertion( (*i) == (*j));
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
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                CGAL_assertion_code( C_Iterator j =) i++;
                CGAL_assertion(  i !=  j);
                if ( i != c_end) {
                    CGAL_assertion( (*i) == (*j) + 1);
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
        C_Iterator i = c_end;
        int su = 0;
        int k  = 5;
        do {
            CGAL_assertion_code(C_Iterator j =) --i;
            CGAL_assertion(  i ==  j);
            CGAL_assertion( (*i) == (*j));
            CGAL_assertion( k == (*i));
            su += (*i);
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
            C_Iterator j = i--;
            CGAL_assertion(  i !=  j);
            if ( j != c_end) {
                CGAL_assertion( (*i) == (*j) - 1);
            }
            CGAL_assertion( k == (*i));
            su += (*i);
            --k;
        } while (i != c_begin);
        CGAL_assertion( i == c_begin);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_random_access_category(c_begin);
        CGAL::Assert_is_at_least_random_access_category(c_end);
        // Random access.
        int k;
        for( k = 0; k < 5; k++) {
            CGAL_assertion( 1+k == c_begin[k]);
        }
        CGAL_assertion_code(
          int su = c_begin[0]
                 + c_begin[1]
                 + c_begin[2]
                 + c_begin[3]
                 + c_begin[4];)
        CGAL_assertion( su == 15);
    
        // Jump around.
        C_Iterator i = c_begin;
        i += 3;
        CGAL_assertion( 4 == (*i));
        i -= 2;
        CGAL_assertion( 2 == (*i));
        i += 3;
        CGAL_assertion( 5 == (*i));
        i -= 4;
        CGAL_assertion( 1 == (*i));
        CGAL_assertion( i == c_begin);
        C_Iterator j = i + 3;
        CGAL_assertion( 4 == (*j));
        C_Iterator jj = j - 2;
        CGAL_assertion( 2 == (*jj));
        jj = 4 + jj;
        CGAL_assertion( jj == c_end);
        C_Iterator ij = jj - 5;
        ij = jj - 5; // avoids warning with NDEBUG
        CGAL_assertion( ij == c_begin);
    
        // Difference test.
        CGAL_assertion( jj - i == 5  ||  jj - i == 0);
        CGAL_assertion( i + (j-i) == j);
        CGAL_assertion( (j-i) + i == j);
    }
    CGAL::Assert_iterator( c_begin);
    CGAL::Assert_iterator( c_end);
    { // Open own scope to hide local variables.
        CGAL_assertion( c_end - c_begin ==  5);
        CGAL_assertion( c_begin - c_end == -5);
        // Relational operator.
        C_Iterator i = c_begin;
        ++i;
        C_Iterator j = i;
        ++j;
        CGAL_assertion( c_begin < i);
        CGAL_assertion( i < j);
        CGAL_assertion( j < c_end);
        CGAL_assertion( j > i);
        CGAL_assertion( i <= j);
        CGAL_assertion( j >= i);
        CGAL_assertion( i <= i);
        CGAL_assertion( i >= i);
    
        CGAL_assertion( !( i >= j));
        CGAL_assertion( !( j <= i));
        CGAL_assertion( !( i > j));
        CGAL_assertion( !( j < i));
        CGAL_assertion( !( i > i));
        CGAL_assertion( !( i < i));
    }
  }
}
void test_Circulator_identity() {
  {
    typedef In_place_list<item,false> List;
    List l;
    l.push_back( *new item(1));
    l.push_back( *new item(2));
    l.push_back( *new item(3));
    l.push_back( *new item(4));
    l.push_back( *new item(5));
    CGAL_assertion( l.size() == 5);
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
        CGAL_assertion(1==test_value_type(static_cast< VT* >(0)));
        CGAL_assertion(1==test_distance_type(static_cast< DT* >(0)));
    
        // Default constructor.
        Circulator z = Circulator();
        z = Circulator(); // avoids warning with NDEBUG
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Circulator i = begin;
    
        // Check general support for circulators and iterators.
        CGAL_assertion( CGAL::is_empty_range( z, z));
        CGAL_assertion( ! CGAL::is_empty_range( i, begin));
    
        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                CGAL_assertion( k == (*i).key);
                su += (*i).key;
                ++k;
                CGAL_assertion_code( Circulator j =) ++i;
                CGAL_assertion( i ==  j);
                if ( i != begin) {
                    CGAL_assertion( (*i).key == (*j).key);
                }
            } while (i != begin);  // Inequality and equality checked.
        }
        CGAL_assertion( i == begin);  // Equality checked.
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                CGAL_assertion( k == (*i).key);
                su += (*i).key;
                ++k;
                CGAL_assertion_code( Circulator j =) i++;
                CGAL_assertion(  i !=  j);
                if ( i != begin) {
                    CGAL_assertion( (*i).key == (*j).key + 1);
                }
            } while (i != begin);
        }
        CGAL_assertion( i == begin);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        // Change three elements and check post-/pre-increment.
        Circulator i = begin;
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
        CGAL_assertion(1==test_value_type(static_cast< VT* >(0)));
        CGAL_assertion(1==test_distance_type(static_cast< DT* >(0)));
    
        // Default constructor.
        Circulator z = Circulator();
        z = Circulator(); // avoids warning with NDEBUG
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Circulator i = begin;
    
        // Check general support for circulators and iterators.
        CGAL_assertion( CGAL::is_empty_range( z, z));
        CGAL_assertion( ! CGAL::is_empty_range( i, begin));
    
        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                CGAL_assertion( k == (*i).key);
                su += (*i).key;
                ++k;
                CGAL_assertion_code( Circulator j =) ++i;
                CGAL_assertion( i ==  j);
                if ( i != begin) {
                    CGAL_assertion( (*i).key == (*j).key);
                }
            } while (i != begin);  // Inequality and equality checked.
        }
        CGAL_assertion( i == begin);  // Equality checked.
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                CGAL_assertion( k == (*i).key);
                su += (*i).key;
                ++k;
                CGAL_assertion_code( Circulator j =) i++;
                CGAL_assertion(  i !=  j);
                if ( i != begin) {
                    CGAL_assertion( (*i).key == (*j).key + 1);
                }
            } while (i != begin);
        }
        CGAL_assertion( i == begin);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_bidirectional_category(begin);
        CGAL::Assert_is_at_least_bidirectional_category(begin);
        // Loop backwards and pre-decrement.
        Circulator i = begin;
        int su = 0;
        int k  = 5;
        do {
            CGAL_assertion_code(Circulator j =) --i;
            CGAL_assertion(  i ==  j);
            CGAL_assertion( (*i).key == (*j).key);
            CGAL_assertion( k == (*i).key);
            su += (*i).key;
            --k;
        } while (i != begin);
        CGAL_assertion( i == begin);
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = begin;
        su = 0;
        k  = 5;
        // Loop with post-decrement.
        do {
            Circulator j = i--;
            CGAL_assertion(  i !=  j);
            if ( j != begin) {
                CGAL_assertion( (*i).key == (*j).key - 1);
            }
            CGAL_assertion( k == (*i).key);
            su += (*i).key;
            --k;
        } while (i != begin);
        CGAL_assertion( i == begin);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_circulator( begin);
        CGAL::Assert_circulator( begin);
    
        // Check the local type parameters.
        Circulator::value_type      k1;
        k1.key = 1;
        Circulator::reference       k2 = k1;
        (void)k2;
        CGAL_assertion( k2.key == 1);
        Circulator::pointer         p1 = &k1;
        (void)p1;
        CGAL_assertion( (*p1).key == 1);
        k1.key = 3;
        CGAL_assertion( k1.key == 3);
        CGAL_assertion( k2.key == 3);
        CGAL_assertion( (*p1).key == 3);
        k1.key = 6;
        CGAL_assertion( k1.key == 6);
        CGAL_assertion( k2.key == 6);
        CGAL_assertion( (*p1).key == 6);
        CGAL_assertion_code( Circulator::size_type s = 5;)
        CGAL_assertion( s == 5);
        CGAL_assertion_code(Circulator::difference_type d = -5;)
        CGAL_assertion( d == -5);
    
        // Check tests for empty data structures.
        Circulator z = Circulator();
        z = Circulator(); // avoids warning with NDEBUG
        CGAL_assertion(   z == CGAL_CIRC_NULL);
        CGAL_assertion( ! (z != CGAL_CIRC_NULL));
        Circulator i = begin;
        CGAL_assertion( ! (i == CGAL_CIRC_NULL));
        CGAL_assertion(   i != CGAL_CIRC_NULL);
        CGAL_assertion( i == begin);
        CGAL_assertion( i == begin);
        // Do I reach myself.
        ++i;
        Circulator j = i;
        int k = 0;
        do {
            CGAL_assertion( k < 5);
            ++k;
            ++i;
        } while( i != j);
        CGAL_assertion( k == 5);
    }
    { // Open own scope to hide local variables.
        // Do I reach myself backwards.
        Circulator i = begin;
        ++i;
        Circulator j = i;
        int k = 0;
        do {
            CGAL_assertion( k < 5);
            ++k;
            --i;
        } while( i != j);
        CGAL_assertion( k == 5);
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
    CGAL_assertion(1==test_value_type(static_cast< VT* >(0)));
    CGAL_assertion(1==test_distance_type(static_cast< DT* >(0)));

    // Default constructor.
    C_Circulator z = C_Circulator();
    z = C_Circulator(); // avoids warning with NDEBUG
    CGAL::Assert_circulator_or_iterator(z);
    // Copy constructor.
    C_Circulator i = c_begin;

    // Check general support for circulators and iterators.
    CGAL_assertion( CGAL::is_empty_range( z, z));
    CGAL_assertion( ! CGAL::is_empty_range( i, c_begin));

    int su = 0;
    int k  = 1;
    // Check general loop, pre-increment, dereference.
    if (! CGAL::is_empty_range( i, c_begin)) {   // superfluous
        do {
            CGAL_assertion( k == (*i).key);
            su += (*i).key;
            ++k;
            CGAL_assertion_code( C_Circulator j =) ++i;
            CGAL_assertion( i ==  j);
            if ( i != c_begin) {
                CGAL_assertion( (*i).key == (*j).key);
            }
        } while (i != c_begin);  // Inequality and equality checked.
    }
    CGAL_assertion( i == c_begin);  // Equality checked.
    CGAL_assertion( su == 15);

    // Assignment.
    i = c_begin;
    su = 0;
    k  = 1;
    // Loop with post increment.
    if (! CGAL::is_empty_range( i, c_begin)) {   // superfluous
        do {
            CGAL_assertion( k == (*i).key);
            su += (*i).key;
            ++k;
            CGAL_assertion_code( C_Circulator j =) i++;
            CGAL_assertion(  i !=  j);
            if ( i != c_begin) {
                CGAL_assertion( (*i).key == (*j).key + 1);
            }
        } while (i != c_begin);
    }
    CGAL_assertion( i == c_begin);
    CGAL_assertion( su == 15);
}
{ // Open own scope to hide local variables.
    CGAL::Assert_is_at_least_bidirectional_category(c_begin);
    CGAL::Assert_is_at_least_bidirectional_category(c_begin);
    // Loop backwards and pre-decrement.
    C_Circulator i = c_begin;
    int su = 0;
    int k  = 5;
    do {
        CGAL_assertion_code(C_Circulator j =) --i;
        CGAL_assertion(  i ==  j);
        CGAL_assertion( (*i).key == (*j).key);
        CGAL_assertion( k == (*i).key);
        su += (*i).key;
        --k;
    } while (i != c_begin);
    CGAL_assertion( i == c_begin);
    CGAL_assertion( su == 15);

    // Assignment.
    i = c_begin;
    su = 0;
    k  = 5;
    // Loop with post-decrement.
    do {
        C_Circulator j = i--;
        CGAL_assertion(  i !=  j);
        if ( j != c_begin) {
            CGAL_assertion( (*i).key == (*j).key - 1);
        }
        CGAL_assertion( k == (*i).key);
        su += (*i).key;
        --k;
    } while (i != c_begin);
    CGAL_assertion( i == c_begin);
    CGAL_assertion( su == 15);
}
{ // Open own scope to hide local variables.
    CGAL::Assert_circulator( c_begin);
    CGAL::Assert_circulator( c_begin);

    // Check the local type parameters.
    C_Circulator::value_type      k1;
    k1.key = 1;
    C_Circulator::reference       k2 = k1;
    (void)k2;
    CGAL_assertion( k2.key == 1);
    C_Circulator::pointer         p1 = &k1;
    (void)p1;
    CGAL_assertion( (*p1).key == 1);
    k1.key = 3;
    CGAL_assertion( k1.key == 3);
    CGAL_assertion( k2.key == 3);
    CGAL_assertion( (*p1).key == 3);
    k1.key = 6;
    CGAL_assertion( k1.key == 6);
    CGAL_assertion( k2.key == 6);
    CGAL_assertion( (*p1).key == 6);
    CGAL_assertion_code( C_Circulator::size_type s = 5;)
    CGAL_assertion( s == 5);
    CGAL_assertion_code(C_Circulator::difference_type d = -5;)
    CGAL_assertion( d == -5);

    // Check tests for empty data structures.
    C_Circulator z = C_Circulator();
    z = C_Circulator(); // avoids warning with NDEBUG
    CGAL_assertion(   z == CGAL_CIRC_NULL);
    CGAL_assertion( ! (z != CGAL_CIRC_NULL));
    C_Circulator i = c_begin;
    CGAL_assertion( ! (i == CGAL_CIRC_NULL));
    CGAL_assertion(   i != CGAL_CIRC_NULL);
    CGAL_assertion( i == c_begin);
    CGAL_assertion( i == c_begin);
    // Do I reach myself.
    ++i;
    C_Circulator j = i;
    int k = 0;
    do {
        CGAL_assertion( k < 5);
        ++k;
        ++i;
    } while( i != j);
    CGAL_assertion( k == 5);
}
{ // Open own scope to hide local variables.
    // Do I reach myself backwards.
    C_Circulator i = c_begin;
    ++i;
    C_Circulator j = i;
    int k = 0;
    do {
        CGAL_assertion( k < 5);
        ++k;
        --i;
    } while( i != j);
    CGAL_assertion( k == 5);
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
    CGAL_assertion( v.size() == 5);
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
        CGAL_assertion(1==test_value_type(static_cast< VT* >(0)));
        CGAL_assertion(1==test_distance_type(static_cast< DT* >(0)));
    
        // Default constructor.
        Circulator z = Circulator();
        z = Circulator(); // avoids warning with NDEBUG
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Circulator i = begin;
    
        // Check general support for circulators and iterators.
        CGAL_assertion( CGAL::is_empty_range( z, z));
        CGAL_assertion( ! CGAL::is_empty_range( i, begin));
    
        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                CGAL_assertion_code( Circulator j =) ++i;
                CGAL_assertion( i ==  j);
                if ( i != begin) {
                    CGAL_assertion( (*i) == (*j));
                }
            } while (i != begin);  // Inequality and equality checked.
        }
        CGAL_assertion( i == begin);  // Equality checked.
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                CGAL_assertion_code( Circulator j =) i++;
                CGAL_assertion(  i !=  j);
                if ( i != begin) {
                    CGAL_assertion( (*i) == (*j) + 1);
                }
            } while (i != begin);
        }
        CGAL_assertion( i == begin);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        // Change three elements and check post-/pre-increment.
        Circulator i = begin;
        (*i++) = 4;
        CGAL_assertion( 4 == (*begin));
        CGAL_assertion( 2 == (*i));
        (*i++) = 3;
        CGAL_assertion( 3 == (*i));
        (*++i) = 7;
        CGAL_assertion( 7 == (*i));
    
        // Check the setting and reset these elements
        // to their original values.
        i = begin;
        CGAL_assertion( 4 == (*i));
        (*i) = 1;
        i++;
        CGAL_assertion( 3 == (*i));
        (*i++) = 2;
        CGAL_assertion( 3 == (*i));
        i++;
        CGAL_assertion( 7 == (*i));
        (*i) = 4;
    
        // Check the resetting.
        i = begin;
        int k = 1;
        do {
            CGAL_assertion( k == (*i));
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
        CGAL_assertion(1==test_value_type(static_cast< VT* >(0)));
        CGAL_assertion(1==test_distance_type(static_cast< DT* >(0)));
    
        // Default constructor.
        Circulator z = Circulator();
        z = Circulator(); // avoids warning with NDEBUG
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Circulator i = begin;
    
        // Check general support for circulators and iterators.
        CGAL_assertion( CGAL::is_empty_range( z, z));
        CGAL_assertion( ! CGAL::is_empty_range( i, begin));
    
        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                CGAL_assertion_code( Circulator j =) ++i;
                CGAL_assertion( i ==  j);
                if ( i != begin) {
                    CGAL_assertion( (*i) == (*j));
                }
            } while (i != begin);  // Inequality and equality checked.
        }
        CGAL_assertion( i == begin);  // Equality checked.
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                CGAL_assertion_code( Circulator j =) i++;
                CGAL_assertion(  i !=  j);
                if ( i != begin) {
                    CGAL_assertion( (*i) == (*j) + 1);
                }
            } while (i != begin);
        }
        CGAL_assertion( i == begin);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_bidirectional_category(begin);
        CGAL::Assert_is_at_least_bidirectional_category(begin);
        // Loop backwards and pre-decrement.
        Circulator i = begin;
        int su = 0;
        int k  = 5;
        do {
            CGAL_assertion_code(Circulator j =) --i;
            CGAL_assertion(  i ==  j);
            CGAL_assertion( (*i) == (*j));
            CGAL_assertion( k == (*i));
            su += (*i);
            --k;
        } while (i != begin);
        CGAL_assertion( i == begin);
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = begin;
        su = 0;
        k  = 5;
        // Loop with post-decrement.
        do {
            Circulator j = i--;
            CGAL_assertion(  i !=  j);
            if ( j != begin) {
                CGAL_assertion( (*i) == (*j) - 1);
            }
            CGAL_assertion( k == (*i));
            su += (*i);
            --k;
        } while (i != begin);
        CGAL_assertion( i == begin);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_is_at_least_forward_category(begin);
        CGAL::Assert_is_at_least_forward_category(begin);
        typedef std::iterator_traits< Circulator >::value_type      VT;
        typedef std::iterator_traits< Circulator >::difference_type DT;
        CGAL_assertion(1==test_value_type(static_cast< VT* >(0)));
        CGAL_assertion(1==test_distance_type(static_cast< DT* >(0)));
    
        // Default constructor.
        Circulator z = Circulator();
        z = Circulator(); // avoids warning with NDEBUG
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Circulator i = begin;
    
        // Check general support for circulators and iterators.
        CGAL_assertion( CGAL::is_empty_range( z, z));
        CGAL_assertion( ! CGAL::is_empty_range( i, begin));
    
        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                CGAL_assertion_code( Circulator j =) ++i;
                CGAL_assertion( i ==  j);
                if ( i != begin) {
                    CGAL_assertion( (*i) == (*j));
                }
            } while (i != begin);  // Inequality and equality checked.
        }
        CGAL_assertion( i == begin);  // Equality checked.
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                CGAL_assertion_code( Circulator j =) i++;
                CGAL_assertion(  i !=  j);
                if ( i != begin) {
                    CGAL_assertion( (*i) == (*j) + 1);
                }
            } while (i != begin);
        }
        CGAL_assertion( i == begin);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_bidirectional_category(begin);
        CGAL::Assert_is_at_least_bidirectional_category(begin);
        // Loop backwards and pre-decrement.
        Circulator i = begin;
        int su = 0;
        int k  = 5;
        do {
            CGAL_assertion_code(Circulator j =) --i;
            CGAL_assertion(  i ==  j);
            CGAL_assertion( (*i) == (*j));
            CGAL_assertion( k == (*i));
            su += (*i);
            --k;
        } while (i != begin);
        CGAL_assertion( i == begin);
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = begin;
        su = 0;
        k  = 5;
        // Loop with post-decrement.
        do {
            Circulator j = i--;
            CGAL_assertion(  i !=  j);
            if ( j != begin) {
                CGAL_assertion( (*i) == (*j) - 1);
            }
            CGAL_assertion( k == (*i));
            su += (*i);
            --k;
        } while (i != begin);
        CGAL_assertion( i == begin);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_random_access_category(begin);
        CGAL::Assert_is_at_least_random_access_category(begin);
        // Random access.
        int k;
        for( k = 0; k < 5; k++) {
            CGAL_assertion( 1+k == begin[k]);
        }
        CGAL_assertion_code(
          int su = begin[0]
                 + begin[1]
                 + begin[2]
                 + begin[3]
                 + begin[4];)
        CGAL_assertion( su == 15);
    
        // Jump around.
        Circulator i = begin;
        i += 3;
        CGAL_assertion( 4 == (*i));
        i -= 2;
        CGAL_assertion( 2 == (*i));
        i += 3;
        CGAL_assertion( 5 == (*i));
        i -= 4;
        CGAL_assertion( 1 == (*i));
        CGAL_assertion( i == begin);
        Circulator j = i + 3;
        CGAL_assertion( 4 == (*j));
        Circulator jj = j - 2;
        CGAL_assertion( 2 == (*jj));
        jj = 4 + jj;
        CGAL_assertion( jj == begin);
        Circulator ij = jj - 5;
        ij = jj - 5; // avoids warning with NDEBUG
        CGAL_assertion( ij == begin);
    
        // Difference test.
        CGAL_assertion( jj - i == 5  ||  jj - i == 0);
        CGAL_assertion( i + (j-i) == j);
        CGAL_assertion( (j-i) + i == j);
    }
    { // Open own scope to hide local variables.
        Circulator i = begin;
        i[2] = 18;
        i[4] = 9;
        i[3] = 12;
        CGAL_assertion( i[2] == 18);
        CGAL_assertion( i[4] == 9);
        CGAL_assertion( i[3] == 12);
        i[2] = 3;
        i[3] = 4;
        i[4] = 5;
        // Check the resetting.
        i = begin;
        int k = 1;
        do {
            CGAL_assertion( k == (*i));
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
        CGAL_assertion( k2 == 1);
        Circulator::pointer         p1 = &k1;
        (void)p1;
        CGAL_assertion( (*p1) == 1);
        k1 = 3;
        CGAL_assertion( k1 == 3);
        CGAL_assertion( k2 == 3);
        CGAL_assertion( (*p1) == 3);
        k1 = 6;
        CGAL_assertion( k1 == 6);
        CGAL_assertion( k2 == 6);
        CGAL_assertion( (*p1) == 6);
        CGAL_assertion_code( Circulator::size_type s = 5;)
        CGAL_assertion( s == 5);
        CGAL_assertion_code(Circulator::difference_type d = -5;)
        CGAL_assertion( d == -5);
    
        // Check tests for empty data structures.
        Circulator z = Circulator();
        z = Circulator(); // avoids warning with NDEBUG
        CGAL_assertion(   z == CGAL_CIRC_NULL);
        CGAL_assertion( ! (z != CGAL_CIRC_NULL));
        Circulator i = begin;
        CGAL_assertion( ! (i == CGAL_CIRC_NULL));
        CGAL_assertion(   i != CGAL_CIRC_NULL);
        CGAL_assertion( i == begin);
        CGAL_assertion( i == begin);
        // Do I reach myself.
        ++i;
        Circulator j = i;
        int k = 0;
        do {
            CGAL_assertion( k < 5);
            ++k;
            ++i;
        } while( i != j);
        CGAL_assertion( k == 5);
    }
    { // Open own scope to hide local variables.
        // Do I reach myself backwards.
        Circulator i = begin;
        ++i;
        Circulator j = i;
        int k = 0;
        do {
            CGAL_assertion( k < 5);
            ++k;
            --i;
        } while( i != j);
        CGAL_assertion( k == 5);
    }
    { // Open own scope to hide local variables.
        Circulator::difference_type d = begin - begin;
        CGAL_assertion( d == 0);
        d = begin - begin;
        CGAL_assertion( d == 0);
        Circulator i = begin + 1;
        CGAL_assertion( begin - i == 1 ||  begin - i == -1);
        CGAL_assertion( i - begin == 1 ||  i - begin == -1);
        // Check minimal circulator properties.
        i = i.min_circulator();
        Circulator j = i;
        CGAL_assertion( j - i == 0);
        j++;
        CGAL_assertion( j - i == 1);
        j++;
        CGAL_assertion( j - i == 2);
        j++;
        CGAL_assertion( j - i == 3);
        j++;
        CGAL_assertion( j - i == 4);
        j++;
        CGAL_assertion( j - i == 0);
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
    CGAL_assertion(1==test_value_type(static_cast< VT* >(0)));
    CGAL_assertion(1==test_distance_type(static_cast< DT* >(0)));

    // Default constructor.
    C_Circulator z = C_Circulator();
    z = C_Circulator(); // avoids warning with NDEBUG
    CGAL::Assert_circulator_or_iterator(z);
    // Copy constructor.
    C_Circulator i = c_begin;

    // Check general support for circulators and iterators.
    CGAL_assertion( CGAL::is_empty_range( z, z));
    CGAL_assertion( ! CGAL::is_empty_range( i, c_begin));

    int su = 0;
    int k  = 1;
    // Check general loop, pre-increment, dereference.
    if (! CGAL::is_empty_range( i, c_begin)) {   // superfluous
        do {
            CGAL_assertion( k == (*i));
            su += (*i);
            ++k;
            CGAL_assertion_code( C_Circulator j =) ++i;
            CGAL_assertion( i ==  j);
            if ( i != c_begin) {
                CGAL_assertion( (*i) == (*j));
            }
        } while (i != c_begin);  // Inequality and equality checked.
    }
    CGAL_assertion( i == c_begin);  // Equality checked.
    CGAL_assertion( su == 15);

    // Assignment.
    i = c_begin;
    su = 0;
    k  = 1;
    // Loop with post increment.
    if (! CGAL::is_empty_range( i, c_begin)) {   // superfluous
        do {
            CGAL_assertion( k == (*i));
            su += (*i);
            ++k;
            CGAL_assertion_code( C_Circulator j =) i++;
            CGAL_assertion(  i !=  j);
            if ( i != c_begin) {
                CGAL_assertion( (*i) == (*j) + 1);
            }
        } while (i != c_begin);
    }
    CGAL_assertion( i == c_begin);
    CGAL_assertion( su == 15);
}
{ // Open own scope to hide local variables.
    CGAL::Assert_is_at_least_bidirectional_category(c_begin);
    CGAL::Assert_is_at_least_bidirectional_category(c_begin);
    // Loop backwards and pre-decrement.
    C_Circulator i = c_begin;
    int su = 0;
    int k  = 5;
    do {
        CGAL_assertion_code(C_Circulator j =) --i;
        CGAL_assertion(  i ==  j);
        CGAL_assertion( (*i) == (*j));
        CGAL_assertion( k == (*i));
        su += (*i);
        --k;
    } while (i != c_begin);
    CGAL_assertion( i == c_begin);
    CGAL_assertion( su == 15);

    // Assignment.
    i = c_begin;
    su = 0;
    k  = 5;
    // Loop with post-decrement.
    do {
        C_Circulator j = i--;
        CGAL_assertion(  i !=  j);
        if ( j != c_begin) {
            CGAL_assertion( (*i) == (*j) - 1);
        }
        CGAL_assertion( k == (*i));
        su += (*i);
        --k;
    } while (i != c_begin);
    CGAL_assertion( i == c_begin);
    CGAL_assertion( su == 15);
}
{ // Open own scope to hide local variables.
    CGAL::Assert_is_at_least_random_access_category(c_begin);
    CGAL::Assert_is_at_least_random_access_category(c_begin);
    // Random access.
    int k;
    for( k = 0; k < 5; k++) {
        CGAL_assertion( 1+k == c_begin[k]);
    }
    CGAL_assertion_code(
      int su = c_begin[0]
             + c_begin[1]
             + c_begin[2]
             + c_begin[3]
             + c_begin[4];)
    CGAL_assertion( su == 15);

    // Jump around.
    C_Circulator i = c_begin;
    i += 3;
    CGAL_assertion( 4 == (*i));
    i -= 2;
    CGAL_assertion( 2 == (*i));
    i += 3;
    CGAL_assertion( 5 == (*i));
    i -= 4;
    CGAL_assertion( 1 == (*i));
    CGAL_assertion( i == c_begin);
    C_Circulator j = i + 3;
    CGAL_assertion( 4 == (*j));
    C_Circulator jj = j - 2;
    CGAL_assertion( 2 == (*jj));
    jj = 4 + jj;
    CGAL_assertion( jj == c_begin);
    C_Circulator ij = jj - 5;
    ij = jj - 5; // avoids warning with NDEBUG
    CGAL_assertion( ij == c_begin);

    // Difference test.
    CGAL_assertion( jj - i == 5  ||  jj - i == 0);
    CGAL_assertion( i + (j-i) == j);
    CGAL_assertion( (j-i) + i == j);
}
{ // Open own scope to hide local variables.
    CGAL::Assert_circulator( c_begin);
    CGAL::Assert_circulator( c_begin);

    // Check the local type parameters.
    C_Circulator::value_type      k1;
    k1 = 1;
    C_Circulator::reference       k2 = k1;
    (void)k2;
    CGAL_assertion( k2 == 1);
    C_Circulator::pointer         p1 = &k1;
    (void)p1;
    CGAL_assertion( (*p1) == 1);
    k1 = 3;
    CGAL_assertion( k1 == 3);
    CGAL_assertion( k2 == 3);
    CGAL_assertion( (*p1) == 3);
    k1 = 6;
    CGAL_assertion( k1 == 6);
    CGAL_assertion( k2 == 6);
    CGAL_assertion( (*p1) == 6);
    CGAL_assertion_code( C_Circulator::size_type s = 5;)
    CGAL_assertion( s == 5);
    CGAL_assertion_code(C_Circulator::difference_type d = -5;)
    CGAL_assertion( d == -5);

    // Check tests for empty data structures.
    C_Circulator z = C_Circulator();
    z = C_Circulator(); // avoids warning with NDEBUG
    CGAL_assertion(   z == CGAL_CIRC_NULL);
    CGAL_assertion( ! (z != CGAL_CIRC_NULL));
    C_Circulator i = c_begin;
    CGAL_assertion( ! (i == CGAL_CIRC_NULL));
    CGAL_assertion(   i != CGAL_CIRC_NULL);
    CGAL_assertion( i == c_begin);
    CGAL_assertion( i == c_begin);
    // Do I reach myself.
    ++i;
    C_Circulator j = i;
    int k = 0;
    do {
        CGAL_assertion( k < 5);
        ++k;
        ++i;
    } while( i != j);
    CGAL_assertion( k == 5);
}
{ // Open own scope to hide local variables.
    // Do I reach myself backwards.
    C_Circulator i = c_begin;
    ++i;
    C_Circulator j = i;
    int k = 0;
    do {
        CGAL_assertion( k < 5);
        ++k;
        --i;
    } while( i != j);
    CGAL_assertion( k == 5);
}
{ // Open own scope to hide local variables.
    C_Circulator::difference_type d = c_begin - c_begin;
    CGAL_assertion( d == 0);
    d = c_begin - c_begin;
    CGAL_assertion( d == 0);
    C_Circulator i = c_begin + 1;
    CGAL_assertion( c_begin - i == 1 ||  c_begin - i == -1);
    CGAL_assertion( i - c_begin == 1 ||  i - c_begin == -1);
    // Check minimal circulator properties.
    i = i.min_circulator();
    C_Circulator j = i;
    CGAL_assertion( j - i == 0);
    j++;
    CGAL_assertion( j - i == 1);
    j++;
    CGAL_assertion( j - i == 2);
    j++;
    CGAL_assertion( j - i == 3);
    j++;
    CGAL_assertion( j - i == 4);
    j++;
    CGAL_assertion( j - i == 0);
}
  }
}
void test_Iterator_project() {
  {
    typedef In_place_list<item,false> List;
    List l;
    l.push_back( *new item(1));
    l.push_back( *new item(2));
    l.push_back( *new item(3));
    l.push_back( *new item(4));
    l.push_back( *new item(5));
    CGAL_assertion( l.size() == 5);
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
        CGAL_assertion(1==test_value_type(static_cast< VT* >(0)));
        CGAL_assertion(1==test_distance_type(static_cast< DT* >(0)));
    
        // Default constructor.
        C_Iterator z = C_Iterator();
        z = C_Iterator(); // avoids warning with NDEBUG
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        C_Iterator i = c_begin;
    
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
                CGAL_assertion_code( C_Iterator j =) ++i;
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
                CGAL_assertion_code( C_Iterator j =) i++;
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
        C_Iterator i = c_end;
        int su = 0;
        int k  = 5;
        do {
            CGAL_assertion_code(C_Iterator j =) --i;
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
            C_Iterator j = i--;
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
    CGAL_assertion( l.size() == 5);
    typedef Cast_function_object<item,item> Ident;
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
        CGAL_assertion(1==test_value_type(static_cast< VT* >(0)));
        CGAL_assertion(1==test_distance_type(static_cast< DT* >(0)));
    
        // Default constructor.
        C_Iterator z = C_Iterator();
        z = C_Iterator(); // avoids warning with NDEBUG
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        C_Iterator i = c_begin;
    
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
                CGAL_assertion_code( C_Iterator j =) ++i;
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
                CGAL_assertion_code( C_Iterator j =) i++;
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
        C_Iterator i = c_end;
        int su = 0;
        int k  = 5;
        do {
            CGAL_assertion_code(C_Iterator j =) --i;
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
            C_Iterator j = i--;
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
    CGAL_assertion( v.size() == 5);
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
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                CGAL_assertion_code( Iterator j =) ++i;
                CGAL_assertion( i ==  j);
                if ( i != end) {
                    CGAL_assertion( (*i) == (*j));
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
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                CGAL_assertion_code( Iterator j =) i++;
                CGAL_assertion(  i !=  j);
                if ( i != end) {
                    CGAL_assertion( (*i) == (*j) + 1);
                }
            } while (i != end);
        }
        CGAL_assertion( i == end);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        // Change three elements and check post-/pre-increment.
        Iterator i = begin;
        (*i++) = 4;
        CGAL_assertion( 4 == (*begin));
        CGAL_assertion( 2 == (*i));
        (*i++) = 3;
        CGAL_assertion( 3 == (*i));
        (*++i) = 7;
        CGAL_assertion( 7 == (*i));
    
        // Check the setting and reset these elements
        // to their original values.
        i = begin;
        CGAL_assertion( 4 == (*i));
        (*i) = 1;
        i++;
        CGAL_assertion( 3 == (*i));
        (*i++) = 2;
        CGAL_assertion( 3 == (*i));
        i++;
        CGAL_assertion( 7 == (*i));
        (*i) = 4;
    
        // Check the resetting.
        i = begin;
        int k = 1;
        do {
            CGAL_assertion( k == (*i));
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
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                CGAL_assertion_code( Iterator j =) ++i;
                CGAL_assertion( i ==  j);
                if ( i != end) {
                    CGAL_assertion( (*i) == (*j));
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
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                CGAL_assertion_code( Iterator j =) i++;
                CGAL_assertion(  i !=  j);
                if ( i != end) {
                    CGAL_assertion( (*i) == (*j) + 1);
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
            CGAL_assertion( (*i) == (*j));
            CGAL_assertion( k == (*i));
            su += (*i);
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
                CGAL_assertion( (*i) == (*j) - 1);
            }
            CGAL_assertion( k == (*i));
            su += (*i);
            --k;
        } while (i != begin);
        CGAL_assertion( i == begin);
        CGAL_assertion( su == 15);
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
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                CGAL_assertion_code( Iterator j =) ++i;
                CGAL_assertion( i ==  j);
                if ( i != end) {
                    CGAL_assertion( (*i) == (*j));
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
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                CGAL_assertion_code( Iterator j =) i++;
                CGAL_assertion(  i !=  j);
                if ( i != end) {
                    CGAL_assertion( (*i) == (*j) + 1);
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
            CGAL_assertion( (*i) == (*j));
            CGAL_assertion( k == (*i));
            su += (*i);
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
                CGAL_assertion( (*i) == (*j) - 1);
            }
            CGAL_assertion( k == (*i));
            su += (*i);
            --k;
        } while (i != begin);
        CGAL_assertion( i == begin);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_random_access_category(begin);
        CGAL::Assert_is_at_least_random_access_category(end);
        // Random access.
        int k;
        for( k = 0; k < 5; k++) {
            CGAL_assertion( 1+k == begin[k]);
        }
        CGAL_assertion_code(
          int su = begin[0]
                 + begin[1]
                 + begin[2]
                 + begin[3]
                 + begin[4];)
        CGAL_assertion( su == 15);
    
        // Jump around.
        Iterator i = begin;
        i += 3;
        CGAL_assertion( 4 == (*i));
        i -= 2;
        CGAL_assertion( 2 == (*i));
        i += 3;
        CGAL_assertion( 5 == (*i));
        i -= 4;
        CGAL_assertion( 1 == (*i));
        CGAL_assertion( i == begin);
        Iterator j = i + 3;
        CGAL_assertion( 4 == (*j));
        Iterator jj = j - 2;
        CGAL_assertion( 2 == (*jj));
        jj = 4 + jj;
        CGAL_assertion( jj == end);
        Iterator ij = jj - 5;
        ij = jj - 5; // avoids warning with NDEBUG
        CGAL_assertion( ij == begin);
    
        // Difference test.
        CGAL_assertion( jj - i == 5  ||  jj - i == 0);
        CGAL_assertion( i + (j-i) == j);
        CGAL_assertion( (j-i) + i == j);
    }
    { // Open own scope to hide local variables.
        Iterator i = begin;
        i[2] = 18;
        i[4] = 9;
        i[3] = 12;
        CGAL_assertion( i[2] == 18);
        CGAL_assertion( i[4] == 9);
        CGAL_assertion( i[3] == 12);
        i[2] = 3;
        i[3] = 4;
        i[4] = 5;
        // Check the resetting.
        i = begin;
        int k = 1;
        do {
            CGAL_assertion( k == (*i));
            ++i;
            ++k;
        } while (i != end);
    }
    CGAL::Assert_iterator( begin);
    CGAL::Assert_iterator( end);
    { // Open own scope to hide local variables.
        CGAL_assertion( end - begin ==  5);
        CGAL_assertion( begin - end == -5);
        // Relational operator.
        Iterator i = begin;
        ++i;
        Iterator j = i;
        ++j;
        CGAL_assertion( begin < i);
        CGAL_assertion( i < j);
        CGAL_assertion( j < end);
        CGAL_assertion( j > i);
        CGAL_assertion( i <= j);
        CGAL_assertion( j >= i);
        CGAL_assertion( i <= i);
        CGAL_assertion( i >= i);
    
        CGAL_assertion( !( i >= j));
        CGAL_assertion( !( j <= i));
        CGAL_assertion( !( i > j));
        CGAL_assertion( !( j < i));
        CGAL_assertion( !( i > i));
        CGAL_assertion( !( i < i));
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
        CGAL_assertion(1==test_value_type(static_cast< VT* >(0)));
        CGAL_assertion(1==test_distance_type(static_cast< DT* >(0)));
    
        // Default constructor.
        C_Iterator z = C_Iterator();
        z = C_Iterator(); // avoids warning with NDEBUG
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        C_Iterator i = c_begin;
    
        // Check general support for circulators and iterators.
        CGAL_assertion( CGAL::is_empty_range( z, z));
        CGAL_assertion( ! CGAL::is_empty_range( i, c_end));
    
        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, c_end)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                CGAL_assertion_code( C_Iterator j =) ++i;
                CGAL_assertion( i ==  j);
                if ( i != c_end) {
                    CGAL_assertion( (*i) == (*j));
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
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                CGAL_assertion_code( C_Iterator j =) i++;
                CGAL_assertion(  i !=  j);
                if ( i != c_end) {
                    CGAL_assertion( (*i) == (*j) + 1);
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
        C_Iterator i = c_end;
        int su = 0;
        int k  = 5;
        do {
            CGAL_assertion_code(C_Iterator j =) --i;
            CGAL_assertion(  i ==  j);
            CGAL_assertion( (*i) == (*j));
            CGAL_assertion( k == (*i));
            su += (*i);
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
            C_Iterator j = i--;
            CGAL_assertion(  i !=  j);
            if ( j != c_end) {
                CGAL_assertion( (*i) == (*j) - 1);
            }
            CGAL_assertion( k == (*i));
            su += (*i);
            --k;
        } while (i != c_begin);
        CGAL_assertion( i == c_begin);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_random_access_category(c_begin);
        CGAL::Assert_is_at_least_random_access_category(c_end);
        // Random access.
        int k;
        for( k = 0; k < 5; k++) {
            CGAL_assertion( 1+k == c_begin[k]);
        }
        CGAL_assertion_code(
          int su = c_begin[0]
                 + c_begin[1]
                 + c_begin[2]
                 + c_begin[3]
                 + c_begin[4];)
        CGAL_assertion( su == 15);
    
        // Jump around.
        C_Iterator i = c_begin;
        i += 3;
        CGAL_assertion( 4 == (*i));
        i -= 2;
        CGAL_assertion( 2 == (*i));
        i += 3;
        CGAL_assertion( 5 == (*i));
        i -= 4;
        CGAL_assertion( 1 == (*i));
        CGAL_assertion( i == c_begin);
        C_Iterator j = i + 3;
        CGAL_assertion( 4 == (*j));
        C_Iterator jj = j - 2;
        CGAL_assertion( 2 == (*jj));
        jj = 4 + jj;
        CGAL_assertion( jj == c_end);
        C_Iterator ij = jj - 5;
        ij = jj - 5; // avoids warning with NDEBUG
        CGAL_assertion( ij == c_begin);
    
        // Difference test.
        CGAL_assertion( jj - i == 5  ||  jj - i == 0);
        CGAL_assertion( i + (j-i) == j);
        CGAL_assertion( (j-i) + i == j);
    }
    CGAL::Assert_iterator( c_begin);
    CGAL::Assert_iterator( c_end);
    { // Open own scope to hide local variables.
        CGAL_assertion( c_end - c_begin ==  5);
        CGAL_assertion( c_begin - c_end == -5);
        // Relational operator.
        C_Iterator i = c_begin;
        ++i;
        C_Iterator j = i;
        ++j;
        CGAL_assertion( c_begin < i);
        CGAL_assertion( i < j);
        CGAL_assertion( j < c_end);
        CGAL_assertion( j > i);
        CGAL_assertion( i <= j);
        CGAL_assertion( j >= i);
        CGAL_assertion( i <= i);
        CGAL_assertion( i >= i);
    
        CGAL_assertion( !( i >= j));
        CGAL_assertion( !( j <= i));
        CGAL_assertion( !( i > j));
        CGAL_assertion( !( j < i));
        CGAL_assertion( !( i > i));
        CGAL_assertion( !( i < i));
    }
  }
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
    CGAL_assertion( l.size() == 5);
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
        CGAL_assertion(1==test_value_type(static_cast< VT* >(0)));
        CGAL_assertion(1==test_distance_type(static_cast< DT* >(0)));
    
        // Default constructor.
        Circulator z = Circulator();
        z = Circulator(); // avoids warning with NDEBUG
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Circulator i = begin;
    
        // Check general support for circulators and iterators.
        CGAL_assertion( CGAL::is_empty_range( z, z));
        CGAL_assertion( ! CGAL::is_empty_range( i, begin));
    
        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                CGAL_assertion( k == (*i).key);
                su += (*i).key;
                ++k;
                CGAL_assertion_code( Circulator j =) ++i;
                CGAL_assertion( i ==  j);
                if ( i != begin) {
                    CGAL_assertion( (*i).key == (*j).key);
                }
            } while (i != begin);  // Inequality and equality checked.
        }
        CGAL_assertion( i == begin);  // Equality checked.
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                CGAL_assertion( k == (*i).key);
                su += (*i).key;
                ++k;
                CGAL_assertion_code( Circulator j =) i++;
                CGAL_assertion(  i !=  j);
                if ( i != begin) {
                    CGAL_assertion( (*i).key == (*j).key + 1);
                }
            } while (i != begin);
        }
        CGAL_assertion( i == begin);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        // Change three elements and check post-/pre-increment.
        Circulator i = begin;
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
        CGAL_assertion(1==test_value_type(static_cast< VT* >(0)));
        CGAL_assertion(1==test_distance_type(static_cast< DT* >(0)));
    
        // Default constructor.
        Circulator z = Circulator();
        z = Circulator(); // avoids warning with NDEBUG
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Circulator i = begin;
    
        // Check general support for circulators and iterators.
        CGAL_assertion( CGAL::is_empty_range( z, z));
        CGAL_assertion( ! CGAL::is_empty_range( i, begin));
    
        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                CGAL_assertion( k == (*i).key);
                su += (*i).key;
                ++k;
                CGAL_assertion_code( Circulator j =) ++i;
                CGAL_assertion( i ==  j);
                if ( i != begin) {
                    CGAL_assertion( (*i).key == (*j).key);
                }
            } while (i != begin);  // Inequality and equality checked.
        }
        CGAL_assertion( i == begin);  // Equality checked.
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                CGAL_assertion( k == (*i).key);
                su += (*i).key;
                ++k;
                CGAL_assertion_code( Circulator j =) i++;
                CGAL_assertion(  i !=  j);
                if ( i != begin) {
                    CGAL_assertion( (*i).key == (*j).key + 1);
                }
            } while (i != begin);
        }
        CGAL_assertion( i == begin);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_bidirectional_category(begin);
        CGAL::Assert_is_at_least_bidirectional_category(begin);
        // Loop backwards and pre-decrement.
        Circulator i = begin;
        int su = 0;
        int k  = 5;
        do {
            CGAL_assertion_code(Circulator j =) --i;
            CGAL_assertion(  i ==  j);
            CGAL_assertion( (*i).key == (*j).key);
            CGAL_assertion( k == (*i).key);
            su += (*i).key;
            --k;
        } while (i != begin);
        CGAL_assertion( i == begin);
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = begin;
        su = 0;
        k  = 5;
        // Loop with post-decrement.
        do {
            Circulator j = i--;
            CGAL_assertion(  i !=  j);
            if ( j != begin) {
                CGAL_assertion( (*i).key == (*j).key - 1);
            }
            CGAL_assertion( k == (*i).key);
            su += (*i).key;
            --k;
        } while (i != begin);
        CGAL_assertion( i == begin);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_circulator( begin);
        CGAL::Assert_circulator( begin);
    
        // Check the local type parameters.
        Circulator::value_type      k1;
        k1.key = 1;
        Circulator::reference       k2 = k1;
        (void)k2;
        CGAL_assertion( k2.key == 1);
        Circulator::pointer         p1 = &k1;
        (void)p1;
        CGAL_assertion( (*p1).key == 1);
        k1.key = 3;
        CGAL_assertion( k1.key == 3);
        CGAL_assertion( k2.key == 3);
        CGAL_assertion( (*p1).key == 3);
        k1.key = 6;
        CGAL_assertion( k1.key == 6);
        CGAL_assertion( k2.key == 6);
        CGAL_assertion( (*p1).key == 6);
        CGAL_assertion_code( Circulator::size_type s = 5;)
        CGAL_assertion( s == 5);
        CGAL_assertion_code(Circulator::difference_type d = -5;)
        CGAL_assertion( d == -5);
    
        // Check tests for empty data structures.
        Circulator z = Circulator();
        z = Circulator(); // avoids warning with NDEBUG
        CGAL_assertion(   z == CGAL_CIRC_NULL);
        CGAL_assertion( ! (z != CGAL_CIRC_NULL));
        Circulator i = begin;
        CGAL_assertion( ! (i == CGAL_CIRC_NULL));
        CGAL_assertion(   i != CGAL_CIRC_NULL);
        CGAL_assertion( i == begin);
        CGAL_assertion( i == begin);
        // Do I reach myself.
        ++i;
        Circulator j = i;
        int k = 0;
        do {
            CGAL_assertion( k < 5);
            ++k;
            ++i;
        } while( i != j);
        CGAL_assertion( k == 5);
    }
    { // Open own scope to hide local variables.
        // Do I reach myself backwards.
        Circulator i = begin;
        ++i;
        Circulator j = i;
        int k = 0;
        do {
            CGAL_assertion( k < 5);
            ++k;
            --i;
        } while( i != j);
        CGAL_assertion( k == 5);
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
    CGAL_assertion(1==test_value_type(static_cast< VT* >(0)));
    CGAL_assertion(1==test_distance_type(static_cast< DT* >(0)));

    // Default constructor.
    C_Circulator z = C_Circulator();
    z = C_Circulator(); // avoids warning with NDEBUG
    CGAL::Assert_circulator_or_iterator(z);
    // Copy constructor.
    C_Circulator i = c_begin;

    // Check general support for circulators and iterators.
    CGAL_assertion( CGAL::is_empty_range( z, z));
    CGAL_assertion( ! CGAL::is_empty_range( i, c_begin));

    int su = 0;
    int k  = 1;
    // Check general loop, pre-increment, dereference.
    if (! CGAL::is_empty_range( i, c_begin)) {   // superfluous
        do {
            CGAL_assertion( k == (*i).key);
            su += (*i).key;
            ++k;
            CGAL_assertion_code( C_Circulator j =) ++i;
            CGAL_assertion( i ==  j);
            if ( i != c_begin) {
                CGAL_assertion( (*i).key == (*j).key);
            }
        } while (i != c_begin);  // Inequality and equality checked.
    }
    CGAL_assertion( i == c_begin);  // Equality checked.
    CGAL_assertion( su == 15);

    // Assignment.
    i = c_begin;
    su = 0;
    k  = 1;
    // Loop with post increment.
    if (! CGAL::is_empty_range( i, c_begin)) {   // superfluous
        do {
            CGAL_assertion( k == (*i).key);
            su += (*i).key;
            ++k;
            CGAL_assertion_code( C_Circulator j =) i++;
            CGAL_assertion(  i !=  j);
            if ( i != c_begin) {
                CGAL_assertion( (*i).key == (*j).key + 1);
            }
        } while (i != c_begin);
    }
    CGAL_assertion( i == c_begin);
    CGAL_assertion( su == 15);
}
{ // Open own scope to hide local variables.
    CGAL::Assert_is_at_least_bidirectional_category(c_begin);
    CGAL::Assert_is_at_least_bidirectional_category(c_begin);
    // Loop backwards and pre-decrement.
    C_Circulator i = c_begin;
    int su = 0;
    int k  = 5;
    do {
        CGAL_assertion_code(C_Circulator j =) --i;
        CGAL_assertion(  i ==  j);
        CGAL_assertion( (*i).key == (*j).key);
        CGAL_assertion( k == (*i).key);
        su += (*i).key;
        --k;
    } while (i != c_begin);
    CGAL_assertion( i == c_begin);
    CGAL_assertion( su == 15);

    // Assignment.
    i = c_begin;
    su = 0;
    k  = 5;
    // Loop with post-decrement.
    do {
        C_Circulator j = i--;
        CGAL_assertion(  i !=  j);
        if ( j != c_begin) {
            CGAL_assertion( (*i).key == (*j).key - 1);
        }
        CGAL_assertion( k == (*i).key);
        su += (*i).key;
        --k;
    } while (i != c_begin);
    CGAL_assertion( i == c_begin);
    CGAL_assertion( su == 15);
}
{ // Open own scope to hide local variables.
    CGAL::Assert_circulator( c_begin);
    CGAL::Assert_circulator( c_begin);

    // Check the local type parameters.
    C_Circulator::value_type      k1;
    k1.key = 1;
    C_Circulator::reference       k2 = k1;
    (void)k2;
    CGAL_assertion( k2.key == 1);
    C_Circulator::pointer         p1 = &k1;
    (void)p1;
    CGAL_assertion( (*p1).key == 1);
    k1.key = 3;
    CGAL_assertion( k1.key == 3);
    CGAL_assertion( k2.key == 3);
    CGAL_assertion( (*p1).key == 3);
    k1.key = 6;
    CGAL_assertion( k1.key == 6);
    CGAL_assertion( k2.key == 6);
    CGAL_assertion( (*p1).key == 6);
    CGAL_assertion_code( C_Circulator::size_type s = 5;)
    CGAL_assertion( s == 5);
    CGAL_assertion_code(C_Circulator::difference_type d = -5;)
    CGAL_assertion( d == -5);

    // Check tests for empty data structures.
    C_Circulator z = C_Circulator();
    z = C_Circulator(); // avoids warning with NDEBUG
    CGAL_assertion(   z == CGAL_CIRC_NULL);
    CGAL_assertion( ! (z != CGAL_CIRC_NULL));
    C_Circulator i = c_begin;
    CGAL_assertion( ! (i == CGAL_CIRC_NULL));
    CGAL_assertion(   i != CGAL_CIRC_NULL);
    CGAL_assertion( i == c_begin);
    CGAL_assertion( i == c_begin);
    // Do I reach myself.
    ++i;
    C_Circulator j = i;
    int k = 0;
    do {
        CGAL_assertion( k < 5);
        ++k;
        ++i;
    } while( i != j);
    CGAL_assertion( k == 5);
}
{ // Open own scope to hide local variables.
    // Do I reach myself backwards.
    C_Circulator i = c_begin;
    ++i;
    C_Circulator j = i;
    int k = 0;
    do {
        CGAL_assertion( k < 5);
        ++k;
        --i;
    } while( i != j);
    CGAL_assertion( k == 5);
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
    CGAL_assertion( l.size() == 5);
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
        CGAL_assertion(1==test_value_type(static_cast< VT* >(0)));
        CGAL_assertion(1==test_distance_type(static_cast< DT* >(0)));
    
        // Default constructor.
        Circulator z = Circulator();
        z = Circulator(); // avoids warning with NDEBUG
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Circulator i = begin;
    
        // Check general support for circulators and iterators.
        CGAL_assertion( CGAL::is_empty_range( z, z));
        CGAL_assertion( ! CGAL::is_empty_range( i, begin));
    
        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                CGAL_assertion( k == (*i).key);
                su += (*i).key;
                ++k;
                CGAL_assertion_code( Circulator j =) ++i;
                CGAL_assertion( i ==  j);
                if ( i != begin) {
                    CGAL_assertion( (*i).key == (*j).key);
                }
            } while (i != begin);  // Inequality and equality checked.
        }
        CGAL_assertion( i == begin);  // Equality checked.
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                CGAL_assertion( k == (*i).key);
                su += (*i).key;
                ++k;
                CGAL_assertion_code( Circulator j =) i++;
                CGAL_assertion(  i !=  j);
                if ( i != begin) {
                    CGAL_assertion( (*i).key == (*j).key + 1);
                }
            } while (i != begin);
        }
        CGAL_assertion( i == begin);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        // Change three elements and check post-/pre-increment.
        Circulator i = begin;
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
        CGAL_assertion(1==test_value_type(static_cast< VT* >(0)));
        CGAL_assertion(1==test_distance_type(static_cast< DT* >(0)));
    
        // Default constructor.
        Circulator z = Circulator();
        z = Circulator(); // avoids warning with NDEBUG
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Circulator i = begin;
    
        // Check general support for circulators and iterators.
        CGAL_assertion( CGAL::is_empty_range( z, z));
        CGAL_assertion( ! CGAL::is_empty_range( i, begin));
    
        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                CGAL_assertion( k == (*i).key);
                su += (*i).key;
                ++k;
                CGAL_assertion_code( Circulator j =) ++i;
                CGAL_assertion( i ==  j);
                if ( i != begin) {
                    CGAL_assertion( (*i).key == (*j).key);
                }
            } while (i != begin);  // Inequality and equality checked.
        }
        CGAL_assertion( i == begin);  // Equality checked.
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                CGAL_assertion( k == (*i).key);
                su += (*i).key;
                ++k;
                CGAL_assertion_code( Circulator j =) i++;
                CGAL_assertion(  i !=  j);
                if ( i != begin) {
                    CGAL_assertion( (*i).key == (*j).key + 1);
                }
            } while (i != begin);
        }
        CGAL_assertion( i == begin);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_bidirectional_category(begin);
        CGAL::Assert_is_at_least_bidirectional_category(begin);
        // Loop backwards and pre-decrement.
        Circulator i = begin;
        int su = 0;
        int k  = 5;
        do {
            CGAL_assertion_code(Circulator j =) --i;
            CGAL_assertion(  i ==  j);
            CGAL_assertion( (*i).key == (*j).key);
            CGAL_assertion( k == (*i).key);
            su += (*i).key;
            --k;
        } while (i != begin);
        CGAL_assertion( i == begin);
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = begin;
        su = 0;
        k  = 5;
        // Loop with post-decrement.
        do {
            Circulator j = i--;
            CGAL_assertion(  i !=  j);
            if ( j != begin) {
                CGAL_assertion( (*i).key == (*j).key - 1);
            }
            CGAL_assertion( k == (*i).key);
            su += (*i).key;
            --k;
        } while (i != begin);
        CGAL_assertion( i == begin);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_circulator( begin);
        CGAL::Assert_circulator( begin);
    
        // Check the local type parameters.
        Circulator::value_type      k1;
        k1.key = 1;
        Circulator::reference       k2 = k1;
        (void)k2;
        CGAL_assertion( k2.key == 1);
        Circulator::pointer         p1 = &k1;
        (void)p1;
        CGAL_assertion( (*p1).key == 1);
        k1.key = 3;
        CGAL_assertion( k1.key == 3);
        CGAL_assertion( k2.key == 3);
        CGAL_assertion( (*p1).key == 3);
        k1.key = 6;
        CGAL_assertion( k1.key == 6);
        CGAL_assertion( k2.key == 6);
        CGAL_assertion( (*p1).key == 6);
        CGAL_assertion_code( Circulator::size_type s = 5;)
        CGAL_assertion( s == 5);
        CGAL_assertion_code(Circulator::difference_type d = -5;)
        CGAL_assertion( d == -5);
    
        // Check tests for empty data structures.
        Circulator z = Circulator();
        z = Circulator(); // avoids warning with NDEBUG
        CGAL_assertion(   z == CGAL_CIRC_NULL);
        CGAL_assertion( ! (z != CGAL_CIRC_NULL));
        Circulator i = begin;
        CGAL_assertion( ! (i == CGAL_CIRC_NULL));
        CGAL_assertion(   i != CGAL_CIRC_NULL);
        CGAL_assertion( i == begin);
        CGAL_assertion( i == begin);
        // Do I reach myself.
        ++i;
        Circulator j = i;
        int k = 0;
        do {
            CGAL_assertion( k < 5);
            ++k;
            ++i;
        } while( i != j);
        CGAL_assertion( k == 5);
    }
    { // Open own scope to hide local variables.
        // Do I reach myself backwards.
        Circulator i = begin;
        ++i;
        Circulator j = i;
        int k = 0;
        do {
            CGAL_assertion( k < 5);
            ++k;
            --i;
        } while( i != j);
        CGAL_assertion( k == 5);
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
    CGAL_assertion(1==test_value_type(static_cast< VT* >(0)));
    CGAL_assertion(1==test_distance_type(static_cast< DT* >(0)));

    // Default constructor.
    C_Circulator z = C_Circulator();
    z = C_Circulator(); // avoids warning with NDEBUG
    CGAL::Assert_circulator_or_iterator(z);
    // Copy constructor.
    C_Circulator i = c_begin;

    // Check general support for circulators and iterators.
    CGAL_assertion( CGAL::is_empty_range( z, z));
    CGAL_assertion( ! CGAL::is_empty_range( i, c_begin));

    int su = 0;
    int k  = 1;
    // Check general loop, pre-increment, dereference.
    if (! CGAL::is_empty_range( i, c_begin)) {   // superfluous
        do {
            CGAL_assertion( k == (*i).key);
            su += (*i).key;
            ++k;
            CGAL_assertion_code( C_Circulator j =) ++i;
            CGAL_assertion( i ==  j);
            if ( i != c_begin) {
                CGAL_assertion( (*i).key == (*j).key);
            }
        } while (i != c_begin);  // Inequality and equality checked.
    }
    CGAL_assertion( i == c_begin);  // Equality checked.
    CGAL_assertion( su == 15);

    // Assignment.
    i = c_begin;
    su = 0;
    k  = 1;
    // Loop with post increment.
    if (! CGAL::is_empty_range( i, c_begin)) {   // superfluous
        do {
            CGAL_assertion( k == (*i).key);
            su += (*i).key;
            ++k;
            CGAL_assertion_code( C_Circulator j =) i++;
            CGAL_assertion(  i !=  j);
            if ( i != c_begin) {
                CGAL_assertion( (*i).key == (*j).key + 1);
            }
        } while (i != c_begin);
    }
    CGAL_assertion( i == c_begin);
    CGAL_assertion( su == 15);
}
{ // Open own scope to hide local variables.
    CGAL::Assert_is_at_least_bidirectional_category(c_begin);
    CGAL::Assert_is_at_least_bidirectional_category(c_begin);
    // Loop backwards and pre-decrement.
    C_Circulator i = c_begin;
    int su = 0;
    int k  = 5;
    do {
        CGAL_assertion_code(C_Circulator j =) --i;
        CGAL_assertion(  i ==  j);
        CGAL_assertion( (*i).key == (*j).key);
        CGAL_assertion( k == (*i).key);
        su += (*i).key;
        --k;
    } while (i != c_begin);
    CGAL_assertion( i == c_begin);
    CGAL_assertion( su == 15);

    // Assignment.
    i = c_begin;
    su = 0;
    k  = 5;
    // Loop with post-decrement.
    do {
        C_Circulator j = i--;
        CGAL_assertion(  i !=  j);
        if ( j != c_begin) {
            CGAL_assertion( (*i).key == (*j).key - 1);
        }
        CGAL_assertion( k == (*i).key);
        su += (*i).key;
        --k;
    } while (i != c_begin);
    CGAL_assertion( i == c_begin);
    CGAL_assertion( su == 15);
}
{ // Open own scope to hide local variables.
    CGAL::Assert_circulator( c_begin);
    CGAL::Assert_circulator( c_begin);

    // Check the local type parameters.
    C_Circulator::value_type      k1;
    k1.key = 1;
    C_Circulator::reference       k2 = k1;
    (void)k2;
    CGAL_assertion( k2.key == 1);
    C_Circulator::pointer         p1 = &k1;
    (void)p1;
    CGAL_assertion( (*p1).key == 1);
    k1.key = 3;
    CGAL_assertion( k1.key == 3);
    CGAL_assertion( k2.key == 3);
    CGAL_assertion( (*p1).key == 3);
    k1.key = 6;
    CGAL_assertion( k1.key == 6);
    CGAL_assertion( k2.key == 6);
    CGAL_assertion( (*p1).key == 6);
    CGAL_assertion_code( C_Circulator::size_type s = 5;)
    CGAL_assertion( s == 5);
    CGAL_assertion_code(C_Circulator::difference_type d = -5;)
    CGAL_assertion( d == -5);

    // Check tests for empty data structures.
    C_Circulator z = C_Circulator();
    z = C_Circulator(); // avoids warning with NDEBUG
    CGAL_assertion(   z == CGAL_CIRC_NULL);
    CGAL_assertion( ! (z != CGAL_CIRC_NULL));
    C_Circulator i = c_begin;
    CGAL_assertion( ! (i == CGAL_CIRC_NULL));
    CGAL_assertion(   i != CGAL_CIRC_NULL);
    CGAL_assertion( i == c_begin);
    CGAL_assertion( i == c_begin);
    // Do I reach myself.
    ++i;
    C_Circulator j = i;
    int k = 0;
    do {
        CGAL_assertion( k < 5);
        ++k;
        ++i;
    } while( i != j);
    CGAL_assertion( k == 5);
}
{ // Open own scope to hide local variables.
    // Do I reach myself backwards.
    C_Circulator i = c_begin;
    ++i;
    C_Circulator j = i;
    int k = 0;
    do {
        CGAL_assertion( k < 5);
        ++k;
        --i;
    } while( i != j);
    CGAL_assertion( k == 5);
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
    CGAL_assertion( v.size() == 5);
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
        CGAL_assertion(1==test_value_type(static_cast< VT* >(0)));
        CGAL_assertion(1==test_distance_type(static_cast< DT* >(0)));
    
        // Default constructor.
        Circulator z = Circulator();
        z = Circulator(); // avoids warning with NDEBUG
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Circulator i = begin;
    
        // Check general support for circulators and iterators.
        CGAL_assertion( CGAL::is_empty_range( z, z));
        CGAL_assertion( ! CGAL::is_empty_range( i, begin));
    
        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                CGAL_assertion_code( Circulator j =) ++i;
                CGAL_assertion( i ==  j);
                if ( i != begin) {
                    CGAL_assertion( (*i) == (*j));
                }
            } while (i != begin);  // Inequality and equality checked.
        }
        CGAL_assertion( i == begin);  // Equality checked.
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                CGAL_assertion_code( Circulator j =) i++;
                CGAL_assertion(  i !=  j);
                if ( i != begin) {
                    CGAL_assertion( (*i) == (*j) + 1);
                }
            } while (i != begin);
        }
        CGAL_assertion( i == begin);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        // Change three elements and check post-/pre-increment.
        Circulator i = begin;
        (*i++) = 4;
        CGAL_assertion( 4 == (*begin));
        CGAL_assertion( 2 == (*i));
        (*i++) = 3;
        CGAL_assertion( 3 == (*i));
        (*++i) = 7;
        CGAL_assertion( 7 == (*i));
    
        // Check the setting and reset these elements
        // to their original values.
        i = begin;
        CGAL_assertion( 4 == (*i));
        (*i) = 1;
        i++;
        CGAL_assertion( 3 == (*i));
        (*i++) = 2;
        CGAL_assertion( 3 == (*i));
        i++;
        CGAL_assertion( 7 == (*i));
        (*i) = 4;
    
        // Check the resetting.
        i = begin;
        int k = 1;
        do {
            CGAL_assertion( k == (*i));
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
        CGAL_assertion(1==test_value_type(static_cast< VT* >(0)));
        CGAL_assertion(1==test_distance_type(static_cast< DT* >(0)));
    
        // Default constructor.
        Circulator z = Circulator();
        z = Circulator(); // avoids warning with NDEBUG
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Circulator i = begin;
    
        // Check general support for circulators and iterators.
        CGAL_assertion( CGAL::is_empty_range( z, z));
        CGAL_assertion( ! CGAL::is_empty_range( i, begin));
    
        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                CGAL_assertion_code( Circulator j =) ++i;
                CGAL_assertion( i ==  j);
                if ( i != begin) {
                    CGAL_assertion( (*i) == (*j));
                }
            } while (i != begin);  // Inequality and equality checked.
        }
        CGAL_assertion( i == begin);  // Equality checked.
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                CGAL_assertion_code( Circulator j =) i++;
                CGAL_assertion(  i !=  j);
                if ( i != begin) {
                    CGAL_assertion( (*i) == (*j) + 1);
                }
            } while (i != begin);
        }
        CGAL_assertion( i == begin);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_bidirectional_category(begin);
        CGAL::Assert_is_at_least_bidirectional_category(begin);
        // Loop backwards and pre-decrement.
        Circulator i = begin;
        int su = 0;
        int k  = 5;
        do {
            CGAL_assertion_code(Circulator j =) --i;
            CGAL_assertion(  i ==  j);
            CGAL_assertion( (*i) == (*j));
            CGAL_assertion( k == (*i));
            su += (*i);
            --k;
        } while (i != begin);
        CGAL_assertion( i == begin);
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = begin;
        su = 0;
        k  = 5;
        // Loop with post-decrement.
        do {
            Circulator j = i--;
            CGAL_assertion(  i !=  j);
            if ( j != begin) {
                CGAL_assertion( (*i) == (*j) - 1);
            }
            CGAL_assertion( k == (*i));
            su += (*i);
            --k;
        } while (i != begin);
        CGAL_assertion( i == begin);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_is_at_least_forward_category(begin);
        CGAL::Assert_is_at_least_forward_category(begin);
        typedef std::iterator_traits< Circulator >::value_type      VT;
        typedef std::iterator_traits< Circulator >::difference_type DT;
        CGAL_assertion(1==test_value_type(static_cast< VT* >(0)));
        CGAL_assertion(1==test_distance_type(static_cast< DT* >(0)));
    
        // Default constructor.
        Circulator z = Circulator();
        z = Circulator(); // avoids warning with NDEBUG
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Circulator i = begin;
    
        // Check general support for circulators and iterators.
        CGAL_assertion( CGAL::is_empty_range( z, z));
        CGAL_assertion( ! CGAL::is_empty_range( i, begin));
    
        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                CGAL_assertion_code( Circulator j =) ++i;
                CGAL_assertion( i ==  j);
                if ( i != begin) {
                    CGAL_assertion( (*i) == (*j));
                }
            } while (i != begin);  // Inequality and equality checked.
        }
        CGAL_assertion( i == begin);  // Equality checked.
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                CGAL_assertion_code( Circulator j =) i++;
                CGAL_assertion(  i !=  j);
                if ( i != begin) {
                    CGAL_assertion( (*i) == (*j) + 1);
                }
            } while (i != begin);
        }
        CGAL_assertion( i == begin);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_bidirectional_category(begin);
        CGAL::Assert_is_at_least_bidirectional_category(begin);
        // Loop backwards and pre-decrement.
        Circulator i = begin;
        int su = 0;
        int k  = 5;
        do {
            CGAL_assertion_code(Circulator j =) --i;
            CGAL_assertion(  i ==  j);
            CGAL_assertion( (*i) == (*j));
            CGAL_assertion( k == (*i));
            su += (*i);
            --k;
        } while (i != begin);
        CGAL_assertion( i == begin);
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = begin;
        su = 0;
        k  = 5;
        // Loop with post-decrement.
        do {
            Circulator j = i--;
            CGAL_assertion(  i !=  j);
            if ( j != begin) {
                CGAL_assertion( (*i) == (*j) - 1);
            }
            CGAL_assertion( k == (*i));
            su += (*i);
            --k;
        } while (i != begin);
        CGAL_assertion( i == begin);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_random_access_category(begin);
        CGAL::Assert_is_at_least_random_access_category(begin);
        // Random access.
        int k;
        for( k = 0; k < 5; k++) {
            CGAL_assertion( 1+k == begin[k]);
        }
        CGAL_assertion_code(
          int su = begin[0]
                 + begin[1]
                 + begin[2]
                 + begin[3]
                 + begin[4];)
        CGAL_assertion( su == 15);
    
        // Jump around.
        Circulator i = begin;
        i += 3;
        CGAL_assertion( 4 == (*i));
        i -= 2;
        CGAL_assertion( 2 == (*i));
        i += 3;
        CGAL_assertion( 5 == (*i));
        i -= 4;
        CGAL_assertion( 1 == (*i));
        CGAL_assertion( i == begin);
        Circulator j = i + 3;
        CGAL_assertion( 4 == (*j));
        Circulator jj = j - 2;
        CGAL_assertion( 2 == (*jj));
        jj = 4 + jj;
        CGAL_assertion( jj == begin);
        Circulator ij = jj - 5;
        ij = jj - 5; // avoids warning with NDEBUG
        CGAL_assertion( ij == begin);
    
        // Difference test.
        CGAL_assertion( jj - i == 5  ||  jj - i == 0);
        CGAL_assertion( i + (j-i) == j);
        CGAL_assertion( (j-i) + i == j);
    }
    { // Open own scope to hide local variables.
        Circulator i = begin;
        i[2] = 18;
        i[4] = 9;
        i[3] = 12;
        CGAL_assertion( i[2] == 18);
        CGAL_assertion( i[4] == 9);
        CGAL_assertion( i[3] == 12);
        i[2] = 3;
        i[3] = 4;
        i[4] = 5;
        // Check the resetting.
        i = begin;
        int k = 1;
        do {
            CGAL_assertion( k == (*i));
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
        CGAL_assertion( k2 == 1);
        Circulator::pointer         p1 = &k1;
        (void)p1;
        CGAL_assertion( (*p1) == 1);
        k1 = 3;
        CGAL_assertion( k1 == 3);
        CGAL_assertion( k2 == 3);
        CGAL_assertion( (*p1) == 3);
        k1 = 6;
        CGAL_assertion( k1 == 6);
        CGAL_assertion( k2 == 6);
        CGAL_assertion( (*p1) == 6);
        CGAL_assertion_code( Circulator::size_type s = 5;)
        CGAL_assertion( s == 5);
        CGAL_assertion_code(Circulator::difference_type d = -5;)
        CGAL_assertion( d == -5);
    
        // Check tests for empty data structures.
        Circulator z = Circulator();
        z = Circulator(); // avoids warning with NDEBUG
        CGAL_assertion(   z == CGAL_CIRC_NULL);
        CGAL_assertion( ! (z != CGAL_CIRC_NULL));
        Circulator i = begin;
        CGAL_assertion( ! (i == CGAL_CIRC_NULL));
        CGAL_assertion(   i != CGAL_CIRC_NULL);
        CGAL_assertion( i == begin);
        CGAL_assertion( i == begin);
        // Do I reach myself.
        ++i;
        Circulator j = i;
        int k = 0;
        do {
            CGAL_assertion( k < 5);
            ++k;
            ++i;
        } while( i != j);
        CGAL_assertion( k == 5);
    }
    { // Open own scope to hide local variables.
        // Do I reach myself backwards.
        Circulator i = begin;
        ++i;
        Circulator j = i;
        int k = 0;
        do {
            CGAL_assertion( k < 5);
            ++k;
            --i;
        } while( i != j);
        CGAL_assertion( k == 5);
    }
    { // Open own scope to hide local variables.
        Circulator::difference_type d = begin - begin;
        CGAL_assertion( d == 0);
        d = begin - begin;
        CGAL_assertion( d == 0);
        Circulator i = begin + 1;
        CGAL_assertion( begin - i == 1 ||  begin - i == -1);
        CGAL_assertion( i - begin == 1 ||  i - begin == -1);
        // Check minimal circulator properties.
        i = i.min_circulator();
        Circulator j = i;
        CGAL_assertion( j - i == 0);
        j++;
        CGAL_assertion( j - i == 1);
        j++;
        CGAL_assertion( j - i == 2);
        j++;
        CGAL_assertion( j - i == 3);
        j++;
        CGAL_assertion( j - i == 4);
        j++;
        CGAL_assertion( j - i == 0);
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
    CGAL_assertion(1==test_value_type(static_cast< VT* >(0)));
    CGAL_assertion(1==test_distance_type(static_cast< DT* >(0)));

    // Default constructor.
    C_Circulator z = C_Circulator();
    z = C_Circulator(); // avoids warning with NDEBUG
    CGAL::Assert_circulator_or_iterator(z);
    // Copy constructor.
    C_Circulator i = c_begin;

    // Check general support for circulators and iterators.
    CGAL_assertion( CGAL::is_empty_range( z, z));
    CGAL_assertion( ! CGAL::is_empty_range( i, c_begin));

    int su = 0;
    int k  = 1;
    // Check general loop, pre-increment, dereference.
    if (! CGAL::is_empty_range( i, c_begin)) {   // superfluous
        do {
            CGAL_assertion( k == (*i));
            su += (*i);
            ++k;
            CGAL_assertion_code( C_Circulator j =) ++i;
            CGAL_assertion( i ==  j);
            if ( i != c_begin) {
                CGAL_assertion( (*i) == (*j));
            }
        } while (i != c_begin);  // Inequality and equality checked.
    }
    CGAL_assertion( i == c_begin);  // Equality checked.
    CGAL_assertion( su == 15);

    // Assignment.
    i = c_begin;
    su = 0;
    k  = 1;
    // Loop with post increment.
    if (! CGAL::is_empty_range( i, c_begin)) {   // superfluous
        do {
            CGAL_assertion( k == (*i));
            su += (*i);
            ++k;
            CGAL_assertion_code( C_Circulator j =) i++;
            CGAL_assertion(  i !=  j);
            if ( i != c_begin) {
                CGAL_assertion( (*i) == (*j) + 1);
            }
        } while (i != c_begin);
    }
    CGAL_assertion( i == c_begin);
    CGAL_assertion( su == 15);
}
{ // Open own scope to hide local variables.
    CGAL::Assert_is_at_least_bidirectional_category(c_begin);
    CGAL::Assert_is_at_least_bidirectional_category(c_begin);
    // Loop backwards and pre-decrement.
    C_Circulator i = c_begin;
    int su = 0;
    int k  = 5;
    do {
        CGAL_assertion_code(C_Circulator j =) --i;
        CGAL_assertion(  i ==  j);
        CGAL_assertion( (*i) == (*j));
        CGAL_assertion( k == (*i));
        su += (*i);
        --k;
    } while (i != c_begin);
    CGAL_assertion( i == c_begin);
    CGAL_assertion( su == 15);

    // Assignment.
    i = c_begin;
    su = 0;
    k  = 5;
    // Loop with post-decrement.
    do {
        C_Circulator j = i--;
        CGAL_assertion(  i !=  j);
        if ( j != c_begin) {
            CGAL_assertion( (*i) == (*j) - 1);
        }
        CGAL_assertion( k == (*i));
        su += (*i);
        --k;
    } while (i != c_begin);
    CGAL_assertion( i == c_begin);
    CGAL_assertion( su == 15);
}
{ // Open own scope to hide local variables.
    CGAL::Assert_is_at_least_random_access_category(c_begin);
    CGAL::Assert_is_at_least_random_access_category(c_begin);
    // Random access.
    int k;
    for( k = 0; k < 5; k++) {
        CGAL_assertion( 1+k == c_begin[k]);
    }
    CGAL_assertion_code(
      int su = c_begin[0]
             + c_begin[1]
             + c_begin[2]
             + c_begin[3]
             + c_begin[4];)
    CGAL_assertion( su == 15);

    // Jump around.
    C_Circulator i = c_begin;
    i += 3;
    CGAL_assertion( 4 == (*i));
    i -= 2;
    CGAL_assertion( 2 == (*i));
    i += 3;
    CGAL_assertion( 5 == (*i));
    i -= 4;
    CGAL_assertion( 1 == (*i));
    CGAL_assertion( i == c_begin);
    C_Circulator j = i + 3;
    CGAL_assertion( 4 == (*j));
    C_Circulator jj = j - 2;
    CGAL_assertion( 2 == (*jj));
    jj = 4 + jj;
    CGAL_assertion( jj == c_begin);
    C_Circulator ij = jj - 5;
    ij = jj - 5; // avoids warning with NDEBUG
    CGAL_assertion( ij == c_begin);

    // Difference test.
    CGAL_assertion( jj - i == 5  ||  jj - i == 0);
    CGAL_assertion( i + (j-i) == j);
    CGAL_assertion( (j-i) + i == j);
}
{ // Open own scope to hide local variables.
    CGAL::Assert_circulator( c_begin);
    CGAL::Assert_circulator( c_begin);

    // Check the local type parameters.
    C_Circulator::value_type      k1;
    k1 = 1;
    C_Circulator::reference       k2 = k1;
    (void)k2;
    CGAL_assertion( k2 == 1);
    C_Circulator::pointer         p1 = &k1;
    (void)p1;
    CGAL_assertion( (*p1) == 1);
    k1 = 3;
    CGAL_assertion( k1 == 3);
    CGAL_assertion( k2 == 3);
    CGAL_assertion( (*p1) == 3);
    k1 = 6;
    CGAL_assertion( k1 == 6);
    CGAL_assertion( k2 == 6);
    CGAL_assertion( (*p1) == 6);
    CGAL_assertion_code( C_Circulator::size_type s = 5;)
    CGAL_assertion( s == 5);
    CGAL_assertion_code(C_Circulator::difference_type d = -5;)
    CGAL_assertion( d == -5);

    // Check tests for empty data structures.
    C_Circulator z = C_Circulator();
    z = C_Circulator(); // avoids warning with NDEBUG
    CGAL_assertion(   z == CGAL_CIRC_NULL);
    CGAL_assertion( ! (z != CGAL_CIRC_NULL));
    C_Circulator i = c_begin;
    CGAL_assertion( ! (i == CGAL_CIRC_NULL));
    CGAL_assertion(   i != CGAL_CIRC_NULL);
    CGAL_assertion( i == c_begin);
    CGAL_assertion( i == c_begin);
    // Do I reach myself.
    ++i;
    C_Circulator j = i;
    int k = 0;
    do {
        CGAL_assertion( k < 5);
        ++k;
        ++i;
    } while( i != j);
    CGAL_assertion( k == 5);
}
{ // Open own scope to hide local variables.
    // Do I reach myself backwards.
    C_Circulator i = c_begin;
    ++i;
    C_Circulator j = i;
    int k = 0;
    do {
        CGAL_assertion( k < 5);
        ++k;
        --i;
    } while( i != j);
    CGAL_assertion( k == 5);
}
{ // Open own scope to hide local variables.
    C_Circulator::difference_type d = c_begin - c_begin;
    CGAL_assertion( d == 0);
    d = c_begin - c_begin;
    CGAL_assertion( d == 0);
    C_Circulator i = c_begin + 1;
    CGAL_assertion( c_begin - i == 1 ||  c_begin - i == -1);
    CGAL_assertion( i - c_begin == 1 ||  i - c_begin == -1);
    // Check minimal circulator properties.
    i = i.min_circulator();
    C_Circulator j = i;
    CGAL_assertion( j - i == 0);
    j++;
    CGAL_assertion( j - i == 1);
    j++;
    CGAL_assertion( j - i == 2);
    j++;
    CGAL_assertion( j - i == 3);
    j++;
    CGAL_assertion( j - i == 4);
    j++;
    CGAL_assertion( j - i == 0);
}
  }
}
struct NN {
  NN* nn;
  int key;
  NN() : nn(0), key(-1) {}
  NN( int k, NN* p) : nn(p), key(k) {}
  NN*       next()       { return nn; }
  const NN* next() const { return nn; }
};
int test_value_type( NN*) { return 1;}

void test_Circulator_on_node() {
  {
    NN* end   = new NN( 5, 0);
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
        CGAL_assertion(1==test_value_type(static_cast< VT* >(0)));
        CGAL_assertion(1==test_distance_type(static_cast< DT* >(0)));
    
        // Default constructor.
        Circulator z = Circulator();
        z = Circulator(); // avoids warning with NDEBUG
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Circulator i = begin;
    
        // Check general support for circulators and iterators.
        CGAL_assertion( CGAL::is_empty_range( z, z));
        CGAL_assertion( ! CGAL::is_empty_range( i, begin));
    
        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                CGAL_assertion( k == (*i).key);
                su += (*i).key;
                ++k;
                CGAL_assertion_code( Circulator j =) ++i;
                CGAL_assertion( i ==  j);
                if ( i != begin) {
                    CGAL_assertion( (*i).key == (*j).key);
                }
            } while (i != begin);  // Inequality and equality checked.
        }
        CGAL_assertion( i == begin);  // Equality checked.
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                CGAL_assertion( k == (*i).key);
                su += (*i).key;
                ++k;
                CGAL_assertion_code( Circulator j =) i++;
                CGAL_assertion(  i !=  j);
                if ( i != begin) {
                    CGAL_assertion( (*i).key == (*j).key + 1);
                }
            } while (i != begin);
        }
        CGAL_assertion( i == begin);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        // Change three elements and check post-/pre-increment.
        Circulator i = begin;
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
        CGAL_assertion( k2.key == 1);
        Circulator::pointer         p1 = &k1;
        (void)p1;
        CGAL_assertion( (*p1).key == 1);
        k1.key = 3;
        CGAL_assertion( k1.key == 3);
        CGAL_assertion( k2.key == 3);
        CGAL_assertion( (*p1).key == 3);
        k1.key = 6;
        CGAL_assertion( k1.key == 6);
        CGAL_assertion( k2.key == 6);
        CGAL_assertion( (*p1).key == 6);
        CGAL_assertion_code( Circulator::size_type s = 5;)
        CGAL_assertion( s == 5);
        CGAL_assertion_code(Circulator::difference_type d = -5;)
        CGAL_assertion( d == -5);
    
        // Check tests for empty data structures.
        Circulator z = Circulator();
        z = Circulator(); // avoids warning with NDEBUG
        CGAL_assertion(   z == CGAL_CIRC_NULL);
        CGAL_assertion( ! (z != CGAL_CIRC_NULL));
        Circulator i = begin;
        CGAL_assertion( ! (i == CGAL_CIRC_NULL));
        CGAL_assertion(   i != CGAL_CIRC_NULL);
        CGAL_assertion( i == begin);
        CGAL_assertion( i == begin);
        // Do I reach myself.
        ++i;
        Circulator j = i;
        int k = 0;
        do {
            CGAL_assertion( k < 5);
            ++k;
            ++i;
        } while( i != j);
        CGAL_assertion( k == 5);
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
        CGAL_assertion(1==test_value_type(static_cast< VT* >(0)));
        CGAL_assertion(1==test_distance_type(static_cast< DT* >(0)));
    
        // Default constructor.
        C_Circulator z = C_Circulator();
        z = C_Circulator(); // avoids warning with NDEBUG
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        C_Circulator i = c_begin;
    
        // Check general support for circulators and iterators.
        CGAL_assertion( CGAL::is_empty_range( z, z));
        CGAL_assertion( ! CGAL::is_empty_range( i, c_begin));
    
        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, c_begin)) {   // superfluous
            do {
                CGAL_assertion( k == (*i).key);
                su += (*i).key;
                ++k;
                CGAL_assertion_code( C_Circulator j =) ++i;
                CGAL_assertion( i ==  j);
                if ( i != c_begin) {
                    CGAL_assertion( (*i).key == (*j).key);
                }
            } while (i != c_begin);  // Inequality and equality checked.
        }
        CGAL_assertion( i == c_begin);  // Equality checked.
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = c_begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, c_begin)) {   // superfluous
            do {
                CGAL_assertion( k == (*i).key);
                su += (*i).key;
                ++k;
                CGAL_assertion_code( C_Circulator j =) i++;
                CGAL_assertion(  i !=  j);
                if ( i != c_begin) {
                    CGAL_assertion( (*i).key == (*j).key + 1);
                }
            } while (i != c_begin);
        }
        CGAL_assertion( i == c_begin);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_circulator( c_begin);
        CGAL::Assert_circulator( c_begin);
    
        // Check the local type parameters.
        C_Circulator::value_type      k1;
        k1.key = 1;
        C_Circulator::reference       k2 = k1;
        (void)k2;
        CGAL_assertion( k2.key == 1);
        C_Circulator::pointer         p1 = &k1;
        (void)p1;
        CGAL_assertion( (*p1).key == 1);
        k1.key = 3;
        CGAL_assertion( k1.key == 3);
        CGAL_assertion( k2.key == 3);
        CGAL_assertion( (*p1).key == 3);
        k1.key = 6;
        CGAL_assertion( k1.key == 6);
        CGAL_assertion( k2.key == 6);
        CGAL_assertion( (*p1).key == 6);
        CGAL_assertion_code( C_Circulator::size_type s = 5;)
        CGAL_assertion( s == 5);
        CGAL_assertion_code(C_Circulator::difference_type d = -5;)
        CGAL_assertion( d == -5);
    
        // Check tests for empty data structures.
        C_Circulator z = C_Circulator();
        z = C_Circulator(); // avoids warning with NDEBUG
        CGAL_assertion(   z == CGAL_CIRC_NULL);
        CGAL_assertion( ! (z != CGAL_CIRC_NULL));
        C_Circulator i = c_begin;
        CGAL_assertion( ! (i == CGAL_CIRC_NULL));
        CGAL_assertion(   i != CGAL_CIRC_NULL);
        CGAL_assertion( i == c_begin);
        CGAL_assertion( i == c_begin);
        // Do I reach myself.
        ++i;
        C_Circulator j = i;
        int k = 0;
        do {
            CGAL_assertion( k < 5);
            ++k;
            ++i;
        } while( i != j);
        CGAL_assertion( k == 5);
    }
    while ( start != end) {
      p = start->nn;
      delete start;
      start = p;
    }
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
    CGAL_assertion( l.size() == 10);
    typedef List::iterator IterBase;
    typedef N_step_adaptor<IterBase,2,item&,item*,item,std::ptrdiff_t,
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

    List l2 = l;
    const List& l1 = l2;
    typedef List::const_iterator ConstIterBase;
    typedef N_step_adaptor<ConstIterBase,2,const item&,const item*,
      item,std::ptrdiff_t, std::bidirectional_iterator_tag>
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
    CGAL_assertion(1==test_value_type(static_cast< VT* >(0)));
    CGAL_assertion(1==test_distance_type(static_cast< DT* >(0)));

    // Default constructor.
    C_Iterator z = C_Iterator();
    z = C_Iterator(); // avoids warning with NDEBUG
    CGAL::Assert_circulator_or_iterator(z);
    // Copy constructor.
    C_Iterator i = c_begin;

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
            CGAL_assertion_code( C_Iterator j =) ++i;
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
            CGAL_assertion_code( C_Iterator j =) i++;
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
    C_Iterator i = c_end;
    int su = 0;
    int k  = 5;
    do {
        CGAL_assertion_code(C_Iterator j =) --i;
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
        C_Iterator j = i--;
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
    CGAL_assertion( v.size() == 10);
    typedef Vector::iterator IterBase;
    typedef N_step_adaptor<IterBase,2,int&,int*,int,std::ptrdiff_t,
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
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                CGAL_assertion_code( Iterator j =) ++i;
                CGAL_assertion( i ==  j);
                if ( i != end) {
                    CGAL_assertion( (*i) == (*j));
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
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                CGAL_assertion_code( Iterator j =) i++;
                CGAL_assertion(  i !=  j);
                if ( i != end) {
                    CGAL_assertion( (*i) == (*j) + 1);
                }
            } while (i != end);
        }
        CGAL_assertion( i == end);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        // Change three elements and check post-/pre-increment.
        Iterator i = begin;
        (*i++) = 4;
        CGAL_assertion( 4 == (*begin));
        CGAL_assertion( 2 == (*i));
        (*i++) = 3;
        CGAL_assertion( 3 == (*i));
        (*++i) = 7;
        CGAL_assertion( 7 == (*i));
    
        // Check the setting and reset these elements
        // to their original values.
        i = begin;
        CGAL_assertion( 4 == (*i));
        (*i) = 1;
        i++;
        CGAL_assertion( 3 == (*i));
        (*i++) = 2;
        CGAL_assertion( 3 == (*i));
        i++;
        CGAL_assertion( 7 == (*i));
        (*i) = 4;
    
        // Check the resetting.
        i = begin;
        int k = 1;
        do {
            CGAL_assertion( k == (*i));
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
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                CGAL_assertion_code( Iterator j =) ++i;
                CGAL_assertion( i ==  j);
                if ( i != end) {
                    CGAL_assertion( (*i) == (*j));
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
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                CGAL_assertion_code( Iterator j =) i++;
                CGAL_assertion(  i !=  j);
                if ( i != end) {
                    CGAL_assertion( (*i) == (*j) + 1);
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
            CGAL_assertion( (*i) == (*j));
            CGAL_assertion( k == (*i));
            su += (*i);
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
                CGAL_assertion( (*i) == (*j) - 1);
            }
            CGAL_assertion( k == (*i));
            su += (*i);
            --k;
        } while (i != begin);
        CGAL_assertion( i == begin);
        CGAL_assertion( su == 15);
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
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                CGAL_assertion_code( Iterator j =) ++i;
                CGAL_assertion( i ==  j);
                if ( i != end) {
                    CGAL_assertion( (*i) == (*j));
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
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                CGAL_assertion_code( Iterator j =) i++;
                CGAL_assertion(  i !=  j);
                if ( i != end) {
                    CGAL_assertion( (*i) == (*j) + 1);
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
            CGAL_assertion( (*i) == (*j));
            CGAL_assertion( k == (*i));
            su += (*i);
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
                CGAL_assertion( (*i) == (*j) - 1);
            }
            CGAL_assertion( k == (*i));
            su += (*i);
            --k;
        } while (i != begin);
        CGAL_assertion( i == begin);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_random_access_category(begin);
        CGAL::Assert_is_at_least_random_access_category(end);
        // Random access.
        int k;
        for( k = 0; k < 5; k++) {
            CGAL_assertion( 1+k == begin[k]);
        }
        CGAL_assertion_code(
          int su = begin[0]
                 + begin[1]
                 + begin[2]
                 + begin[3]
                 + begin[4];)
        CGAL_assertion( su == 15);
    
        // Jump around.
        Iterator i = begin;
        i += 3;
        CGAL_assertion( 4 == (*i));
        i -= 2;
        CGAL_assertion( 2 == (*i));
        i += 3;
        CGAL_assertion( 5 == (*i));
        i -= 4;
        CGAL_assertion( 1 == (*i));
        CGAL_assertion( i == begin);
        Iterator j = i + 3;
        CGAL_assertion( 4 == (*j));
        Iterator jj = j - 2;
        CGAL_assertion( 2 == (*jj));
        jj = 4 + jj;
        CGAL_assertion( jj == end);
        Iterator ij = jj - 5;
        ij = jj - 5; // avoids warning with NDEBUG
        CGAL_assertion( ij == begin);
    
        // Difference test.
        CGAL_assertion( jj - i == 5  ||  jj - i == 0);
        CGAL_assertion( i + (j-i) == j);
        CGAL_assertion( (j-i) + i == j);
    }
    { // Open own scope to hide local variables.
        Iterator i = begin;
        i[2] = 18;
        i[4] = 9;
        i[3] = 12;
        CGAL_assertion( i[2] == 18);
        CGAL_assertion( i[4] == 9);
        CGAL_assertion( i[3] == 12);
        i[2] = 3;
        i[3] = 4;
        i[4] = 5;
        // Check the resetting.
        i = begin;
        int k = 1;
        do {
            CGAL_assertion( k == (*i));
            ++i;
            ++k;
        } while (i != end);
    }
    CGAL::Assert_iterator( begin);
    CGAL::Assert_iterator( end);
    { // Open own scope to hide local variables.
        CGAL_assertion( end - begin ==  5);
        CGAL_assertion( begin - end == -5);
        // Relational operator.
        Iterator i = begin;
        ++i;
        Iterator j = i;
        ++j;
        CGAL_assertion( begin < i);
        CGAL_assertion( i < j);
        CGAL_assertion( j < end);
        CGAL_assertion( j > i);
        CGAL_assertion( i <= j);
        CGAL_assertion( j >= i);
        CGAL_assertion( i <= i);
        CGAL_assertion( i >= i);
    
        CGAL_assertion( !( i >= j));
        CGAL_assertion( !( j <= i));
        CGAL_assertion( !( i > j));
        CGAL_assertion( !( j < i));
        CGAL_assertion( !( i > i));
        CGAL_assertion( !( i < i));
    }

    Vector v2 = v;
    const Vector& v1 = v2;
    typedef Vector::const_iterator ConstIterBase;
    typedef N_step_adaptor< ConstIterBase,2,const int&,
      const int*, int, std::ptrdiff_t,
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
        CGAL_assertion(1==test_value_type(static_cast< VT* >(0)));
        CGAL_assertion(1==test_distance_type(static_cast< DT* >(0)));
    
        // Default constructor.
        C_Iterator z = C_Iterator();
        z = C_Iterator(); // avoids warning with NDEBUG
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        C_Iterator i = c_begin;
    
        // Check general support for circulators and iterators.
        CGAL_assertion( CGAL::is_empty_range( z, z));
        CGAL_assertion( ! CGAL::is_empty_range( i, c_end));
    
        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, c_end)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                CGAL_assertion_code( C_Iterator j =) ++i;
                CGAL_assertion( i ==  j);
                if ( i != c_end) {
                    CGAL_assertion( (*i) == (*j));
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
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                CGAL_assertion_code( C_Iterator j =) i++;
                CGAL_assertion(  i !=  j);
                if ( i != c_end) {
                    CGAL_assertion( (*i) == (*j) + 1);
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
        C_Iterator i = c_end;
        int su = 0;
        int k  = 5;
        do {
            CGAL_assertion_code(C_Iterator j =) --i;
            CGAL_assertion(  i ==  j);
            CGAL_assertion( (*i) == (*j));
            CGAL_assertion( k == (*i));
            su += (*i);
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
            C_Iterator j = i--;
            CGAL_assertion(  i !=  j);
            if ( j != c_end) {
                CGAL_assertion( (*i) == (*j) - 1);
            }
            CGAL_assertion( k == (*i));
            su += (*i);
            --k;
        } while (i != c_begin);
        CGAL_assertion( i == c_begin);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_random_access_category(c_begin);
        CGAL::Assert_is_at_least_random_access_category(c_end);
        // Random access.
        int k;
        for( k = 0; k < 5; k++) {
            CGAL_assertion( 1+k == c_begin[k]);
        }
        CGAL_assertion_code(
          int su = c_begin[0]
                 + c_begin[1]
                 + c_begin[2]
                 + c_begin[3]
                 + c_begin[4];)
        CGAL_assertion( su == 15);
    
        // Jump around.
        C_Iterator i = c_begin;
        i += 3;
        CGAL_assertion( 4 == (*i));
        i -= 2;
        CGAL_assertion( 2 == (*i));
        i += 3;
        CGAL_assertion( 5 == (*i));
        i -= 4;
        CGAL_assertion( 1 == (*i));
        CGAL_assertion( i == c_begin);
        C_Iterator j = i + 3;
        CGAL_assertion( 4 == (*j));
        C_Iterator jj = j - 2;
        CGAL_assertion( 2 == (*jj));
        jj = 4 + jj;
        CGAL_assertion( jj == c_end);
        C_Iterator ij = jj - 5;
        ij = jj - 5; // avoids warning with NDEBUG
        CGAL_assertion( ij == c_begin);
    
        // Difference test.
        CGAL_assertion( jj - i == 5  ||  jj - i == 0);
        CGAL_assertion( i + (j-i) == j);
        CGAL_assertion( (j-i) + i == j);
    }
    CGAL::Assert_iterator( c_begin);
    CGAL::Assert_iterator( c_end);
    { // Open own scope to hide local variables.
        CGAL_assertion( c_end - c_begin ==  5);
        CGAL_assertion( c_begin - c_end == -5);
        // Relational operator.
        C_Iterator i = c_begin;
        ++i;
        C_Iterator j = i;
        ++j;
        CGAL_assertion( c_begin < i);
        CGAL_assertion( i < j);
        CGAL_assertion( j < c_end);
        CGAL_assertion( j > i);
        CGAL_assertion( i <= j);
        CGAL_assertion( j >= i);
        CGAL_assertion( i <= i);
        CGAL_assertion( i >= i);
    
        CGAL_assertion( !( i >= j));
        CGAL_assertion( !( j <= i));
        CGAL_assertion( !( i > j));
        CGAL_assertion( !( j < i));
        CGAL_assertion( !( i > i));
        CGAL_assertion( !( i < i));
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
    CGAL_assertion( v.size() == 10);
    typedef Vector::iterator IterBase;
    typedef Random_access_circulator_from_iterator<IterBase,int,
      std::size_t,std::ptrdiff_t> CircBase;
    typedef N_step_adaptor<CircBase,2,int&,int*,int,std::ptrdiff_t,
      Random_access_circulator_tag> Circulator;
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
        CGAL_assertion(1==test_value_type(static_cast< VT* >(0)));
        CGAL_assertion(1==test_distance_type(static_cast< DT* >(0)));
    
        // Default constructor.
        Circulator z = Circulator();
        z = Circulator(); // avoids warning with NDEBUG
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Circulator i = begin;
    
        // Check general support for circulators and iterators.
        CGAL_assertion( CGAL::is_empty_range( z, z));
        CGAL_assertion( ! CGAL::is_empty_range( i, begin));
    
        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                CGAL_assertion_code( Circulator j =) ++i;
                CGAL_assertion( i ==  j);
                if ( i != begin) {
                    CGAL_assertion( (*i) == (*j));
                }
            } while (i != begin);  // Inequality and equality checked.
        }
        CGAL_assertion( i == begin);  // Equality checked.
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                CGAL_assertion_code( Circulator j =) i++;
                CGAL_assertion(  i !=  j);
                if ( i != begin) {
                    CGAL_assertion( (*i) == (*j) + 1);
                }
            } while (i != begin);
        }
        CGAL_assertion( i == begin);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        // Change three elements and check post-/pre-increment.
        Circulator i = begin;
        (*i++) = 4;
        CGAL_assertion( 4 == (*begin));
        CGAL_assertion( 2 == (*i));
        (*i++) = 3;
        CGAL_assertion( 3 == (*i));
        (*++i) = 7;
        CGAL_assertion( 7 == (*i));
    
        // Check the setting and reset these elements
        // to their original values.
        i = begin;
        CGAL_assertion( 4 == (*i));
        (*i) = 1;
        i++;
        CGAL_assertion( 3 == (*i));
        (*i++) = 2;
        CGAL_assertion( 3 == (*i));
        i++;
        CGAL_assertion( 7 == (*i));
        (*i) = 4;
    
        // Check the resetting.
        i = begin;
        int k = 1;
        do {
            CGAL_assertion( k == (*i));
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
        CGAL_assertion(1==test_value_type(static_cast< VT* >(0)));
        CGAL_assertion(1==test_distance_type(static_cast< DT* >(0)));
    
        // Default constructor.
        Circulator z = Circulator();
        z = Circulator(); // avoids warning with NDEBUG
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Circulator i = begin;
    
        // Check general support for circulators and iterators.
        CGAL_assertion( CGAL::is_empty_range( z, z));
        CGAL_assertion( ! CGAL::is_empty_range( i, begin));
    
        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                CGAL_assertion_code( Circulator j =) ++i;
                CGAL_assertion( i ==  j);
                if ( i != begin) {
                    CGAL_assertion( (*i) == (*j));
                }
            } while (i != begin);  // Inequality and equality checked.
        }
        CGAL_assertion( i == begin);  // Equality checked.
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                CGAL_assertion_code( Circulator j =) i++;
                CGAL_assertion(  i !=  j);
                if ( i != begin) {
                    CGAL_assertion( (*i) == (*j) + 1);
                }
            } while (i != begin);
        }
        CGAL_assertion( i == begin);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_bidirectional_category(begin);
        CGAL::Assert_is_at_least_bidirectional_category(begin);
        // Loop backwards and pre-decrement.
        Circulator i = begin;
        int su = 0;
        int k  = 5;
        do {
            CGAL_assertion_code(Circulator j =) --i;
            CGAL_assertion(  i ==  j);
            CGAL_assertion( (*i) == (*j));
            CGAL_assertion( k == (*i));
            su += (*i);
            --k;
        } while (i != begin);
        CGAL_assertion( i == begin);
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = begin;
        su = 0;
        k  = 5;
        // Loop with post-decrement.
        do {
            Circulator j = i--;
            CGAL_assertion(  i !=  j);
            if ( j != begin) {
                CGAL_assertion( (*i) == (*j) - 1);
            }
            CGAL_assertion( k == (*i));
            su += (*i);
            --k;
        } while (i != begin);
        CGAL_assertion( i == begin);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_is_at_least_forward_category(begin);
        CGAL::Assert_is_at_least_forward_category(begin);
        typedef std::iterator_traits< Circulator >::value_type      VT;
        typedef std::iterator_traits< Circulator >::difference_type DT;
        CGAL_assertion(1==test_value_type(static_cast< VT* >(0)));
        CGAL_assertion(1==test_distance_type(static_cast< DT* >(0)));
    
        // Default constructor.
        Circulator z = Circulator();
        z = Circulator(); // avoids warning with NDEBUG
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Circulator i = begin;
    
        // Check general support for circulators and iterators.
        CGAL_assertion( CGAL::is_empty_range( z, z));
        CGAL_assertion( ! CGAL::is_empty_range( i, begin));
    
        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                CGAL_assertion_code( Circulator j =) ++i;
                CGAL_assertion( i ==  j);
                if ( i != begin) {
                    CGAL_assertion( (*i) == (*j));
                }
            } while (i != begin);  // Inequality and equality checked.
        }
        CGAL_assertion( i == begin);  // Equality checked.
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = begin;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, begin)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                CGAL_assertion_code( Circulator j =) i++;
                CGAL_assertion(  i !=  j);
                if ( i != begin) {
                    CGAL_assertion( (*i) == (*j) + 1);
                }
            } while (i != begin);
        }
        CGAL_assertion( i == begin);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_bidirectional_category(begin);
        CGAL::Assert_is_at_least_bidirectional_category(begin);
        // Loop backwards and pre-decrement.
        Circulator i = begin;
        int su = 0;
        int k  = 5;
        do {
            CGAL_assertion_code(Circulator j =) --i;
            CGAL_assertion(  i ==  j);
            CGAL_assertion( (*i) == (*j));
            CGAL_assertion( k == (*i));
            su += (*i);
            --k;
        } while (i != begin);
        CGAL_assertion( i == begin);
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = begin;
        su = 0;
        k  = 5;
        // Loop with post-decrement.
        do {
            Circulator j = i--;
            CGAL_assertion(  i !=  j);
            if ( j != begin) {
                CGAL_assertion( (*i) == (*j) - 1);
            }
            CGAL_assertion( k == (*i));
            su += (*i);
            --k;
        } while (i != begin);
        CGAL_assertion( i == begin);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_random_access_category(begin);
        CGAL::Assert_is_at_least_random_access_category(begin);
        // Random access.
        int k;
        for( k = 0; k < 5; k++) {
            CGAL_assertion( 1+k == begin[k]);
        }
        CGAL_assertion_code(
          int su = begin[0]
                 + begin[1]
                 + begin[2]
                 + begin[3]
                 + begin[4];)
        CGAL_assertion( su == 15);
    
        // Jump around.
        Circulator i = begin;
        i += 3;
        CGAL_assertion( 4 == (*i));
        i -= 2;
        CGAL_assertion( 2 == (*i));
        i += 3;
        CGAL_assertion( 5 == (*i));
        i -= 4;
        CGAL_assertion( 1 == (*i));
        CGAL_assertion( i == begin);
        Circulator j = i + 3;
        CGAL_assertion( 4 == (*j));
        Circulator jj = j - 2;
        CGAL_assertion( 2 == (*jj));
        jj = 4 + jj;
        CGAL_assertion( jj == begin);
        Circulator ij = jj - 5;
        ij = jj - 5; // avoids warning with NDEBUG
        CGAL_assertion( ij == begin);
    
        // Difference test.
        CGAL_assertion( jj - i == 5  ||  jj - i == 0);
        CGAL_assertion( i + (j-i) == j);
        CGAL_assertion( (j-i) + i == j);
    }
    { // Open own scope to hide local variables.
        Circulator i = begin;
        i[2] = 18;
        i[4] = 9;
        i[3] = 12;
        CGAL_assertion( i[2] == 18);
        CGAL_assertion( i[4] == 9);
        CGAL_assertion( i[3] == 12);
        i[2] = 3;
        i[3] = 4;
        i[4] = 5;
        // Check the resetting.
        i = begin;
        int k = 1;
        do {
            CGAL_assertion( k == (*i));
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
        CGAL_assertion( k2 == 1);
        Circulator::pointer         p1 = &k1;
        (void)p1;
        CGAL_assertion( (*p1) == 1);
        k1 = 3;
        CGAL_assertion( k1 == 3);
        CGAL_assertion( k2 == 3);
        CGAL_assertion( (*p1) == 3);
        k1 = 6;
        CGAL_assertion( k1 == 6);
        CGAL_assertion( k2 == 6);
        CGAL_assertion( (*p1) == 6);
        CGAL_assertion_code( Circulator::size_type s = 5;)
        CGAL_assertion( s == 5);
        CGAL_assertion_code(Circulator::difference_type d = -5;)
        CGAL_assertion( d == -5);
    
        // Check tests for empty data structures.
        Circulator z = Circulator();
        z = Circulator(); // avoids warning with NDEBUG
        CGAL_assertion(   z == CGAL_CIRC_NULL);
        CGAL_assertion( ! (z != CGAL_CIRC_NULL));
        Circulator i = begin;
        CGAL_assertion( ! (i == CGAL_CIRC_NULL));
        CGAL_assertion(   i != CGAL_CIRC_NULL);
        CGAL_assertion( i == begin);
        CGAL_assertion( i == begin);
        // Do I reach myself.
        ++i;
        Circulator j = i;
        int k = 0;
        do {
            CGAL_assertion( k < 5);
            ++k;
            ++i;
        } while( i != j);
        CGAL_assertion( k == 5);
    }
    { // Open own scope to hide local variables.
        // Do I reach myself backwards.
        Circulator i = begin;
        ++i;
        Circulator j = i;
        int k = 0;
        do {
            CGAL_assertion( k < 5);
            ++k;
            --i;
        } while( i != j);
        CGAL_assertion( k == 5);
    }
    { // Open own scope to hide local variables.
        Circulator::difference_type d = begin - begin;
        CGAL_assertion( d == 0);
        d = begin - begin;
        CGAL_assertion( d == 0);
        Circulator i = begin + 1;
        CGAL_assertion( begin - i == 1 ||  begin - i == -1);
        CGAL_assertion( i - begin == 1 ||  i - begin == -1);
        // Check minimal circulator properties.
        i = i.min_circulator();
        Circulator j = i;
        CGAL_assertion( j - i == 0);
        j++;
        CGAL_assertion( j - i == 1);
        j++;
        CGAL_assertion( j - i == 2);
        j++;
        CGAL_assertion( j - i == 3);
        j++;
        CGAL_assertion( j - i == 4);
        j++;
        CGAL_assertion( j - i == 0);
    }

    Vector v2 = v;
    const Vector& v1 = v2;
    typedef Vector::const_iterator ConstIterBase;
    typedef Random_access_const_circulator_from_iterator<
      ConstIterBase, int, std::size_t, std::ptrdiff_t>
    ConstCircBase;
typedef N_step_adaptor<ConstCircBase,2,const int&,const int*,
  int,std::ptrdiff_t, Random_access_circulator_tag> C_Circulator;
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
    CGAL_assertion(1==test_value_type(static_cast< VT* >(0)));
    CGAL_assertion(1==test_distance_type(static_cast< DT* >(0)));

    // Default constructor.
    C_Circulator z = C_Circulator();
    z = C_Circulator(); // avoids warning with NDEBUG
    CGAL::Assert_circulator_or_iterator(z);
    // Copy constructor.
    C_Circulator i = c_begin;

    // Check general support for circulators and iterators.
    CGAL_assertion( CGAL::is_empty_range( z, z));
    CGAL_assertion( ! CGAL::is_empty_range( i, c_begin));

    int su = 0;
    int k  = 1;
    // Check general loop, pre-increment, dereference.
    if (! CGAL::is_empty_range( i, c_begin)) {   // superfluous
        do {
            CGAL_assertion( k == (*i));
            su += (*i);
            ++k;
            CGAL_assertion_code( C_Circulator j =) ++i;
            CGAL_assertion( i ==  j);
            if ( i != c_begin) {
                CGAL_assertion( (*i) == (*j));
            }
        } while (i != c_begin);  // Inequality and equality checked.
    }
    CGAL_assertion( i == c_begin);  // Equality checked.
    CGAL_assertion( su == 15);

    // Assignment.
    i = c_begin;
    su = 0;
    k  = 1;
    // Loop with post increment.
    if (! CGAL::is_empty_range( i, c_begin)) {   // superfluous
        do {
            CGAL_assertion( k == (*i));
            su += (*i);
            ++k;
            CGAL_assertion_code( C_Circulator j =) i++;
            CGAL_assertion(  i !=  j);
            if ( i != c_begin) {
                CGAL_assertion( (*i) == (*j) + 1);
            }
        } while (i != c_begin);
    }
    CGAL_assertion( i == c_begin);
    CGAL_assertion( su == 15);
}
{ // Open own scope to hide local variables.
    CGAL::Assert_is_at_least_bidirectional_category(c_begin);
    CGAL::Assert_is_at_least_bidirectional_category(c_begin);
    // Loop backwards and pre-decrement.
    C_Circulator i = c_begin;
    int su = 0;
    int k  = 5;
    do {
        CGAL_assertion_code(C_Circulator j =) --i;
        CGAL_assertion(  i ==  j);
        CGAL_assertion( (*i) == (*j));
        CGAL_assertion( k == (*i));
        su += (*i);
        --k;
    } while (i != c_begin);
    CGAL_assertion( i == c_begin);
    CGAL_assertion( su == 15);

    // Assignment.
    i = c_begin;
    su = 0;
    k  = 5;
    // Loop with post-decrement.
    do {
        C_Circulator j = i--;
        CGAL_assertion(  i !=  j);
        if ( j != c_begin) {
            CGAL_assertion( (*i) == (*j) - 1);
        }
        CGAL_assertion( k == (*i));
        su += (*i);
        --k;
    } while (i != c_begin);
    CGAL_assertion( i == c_begin);
    CGAL_assertion( su == 15);
}
{ // Open own scope to hide local variables.
    CGAL::Assert_is_at_least_random_access_category(c_begin);
    CGAL::Assert_is_at_least_random_access_category(c_begin);
    // Random access.
    int k;
    for( k = 0; k < 5; k++) {
        CGAL_assertion( 1+k == c_begin[k]);
    }
    CGAL_assertion_code(
      int su = c_begin[0]
             + c_begin[1]
             + c_begin[2]
             + c_begin[3]
             + c_begin[4];)
    CGAL_assertion( su == 15);

    // Jump around.
    C_Circulator i = c_begin;
    i += 3;
    CGAL_assertion( 4 == (*i));
    i -= 2;
    CGAL_assertion( 2 == (*i));
    i += 3;
    CGAL_assertion( 5 == (*i));
    i -= 4;
    CGAL_assertion( 1 == (*i));
    CGAL_assertion( i == c_begin);
    C_Circulator j = i + 3;
    CGAL_assertion( 4 == (*j));
    C_Circulator jj = j - 2;
    CGAL_assertion( 2 == (*jj));
    jj = 4 + jj;
    CGAL_assertion( jj == c_begin);
    C_Circulator ij = jj - 5;
    ij = jj - 5; // avoids warning with NDEBUG
    CGAL_assertion( ij == c_begin);

    // Difference test.
    CGAL_assertion( jj - i == 5  ||  jj - i == 0);
    CGAL_assertion( i + (j-i) == j);
    CGAL_assertion( (j-i) + i == j);
}
{ // Open own scope to hide local variables.
    CGAL::Assert_circulator( c_begin);
    CGAL::Assert_circulator( c_begin);

    // Check the local type parameters.
    C_Circulator::value_type      k1;
    k1 = 1;
    C_Circulator::reference       k2 = k1;
    (void)k2;
    CGAL_assertion( k2 == 1);
    C_Circulator::pointer         p1 = &k1;
    (void)p1;
    CGAL_assertion( (*p1) == 1);
    k1 = 3;
    CGAL_assertion( k1 == 3);
    CGAL_assertion( k2 == 3);
    CGAL_assertion( (*p1) == 3);
    k1 = 6;
    CGAL_assertion( k1 == 6);
    CGAL_assertion( k2 == 6);
    CGAL_assertion( (*p1) == 6);
    CGAL_assertion_code( C_Circulator::size_type s = 5;)
    CGAL_assertion( s == 5);
    CGAL_assertion_code(C_Circulator::difference_type d = -5;)
    CGAL_assertion( d == -5);

    // Check tests for empty data structures.
    C_Circulator z = C_Circulator();
    z = C_Circulator(); // avoids warning with NDEBUG
    CGAL_assertion(   z == CGAL_CIRC_NULL);
    CGAL_assertion( ! (z != CGAL_CIRC_NULL));
    C_Circulator i = c_begin;
    CGAL_assertion( ! (i == CGAL_CIRC_NULL));
    CGAL_assertion(   i != CGAL_CIRC_NULL);
    CGAL_assertion( i == c_begin);
    CGAL_assertion( i == c_begin);
    // Do I reach myself.
    ++i;
    C_Circulator j = i;
    int k = 0;
    do {
        CGAL_assertion( k < 5);
        ++k;
        ++i;
    } while( i != j);
    CGAL_assertion( k == 5);
}
{ // Open own scope to hide local variables.
    // Do I reach myself backwards.
    C_Circulator i = c_begin;
    ++i;
    C_Circulator j = i;
    int k = 0;
    do {
        CGAL_assertion( k < 5);
        ++k;
        --i;
    } while( i != j);
    CGAL_assertion( k == 5);
}
{ // Open own scope to hide local variables.
    C_Circulator::difference_type d = c_begin - c_begin;
    CGAL_assertion( d == 0);
    d = c_begin - c_begin;
    CGAL_assertion( d == 0);
    C_Circulator i = c_begin + 1;
    CGAL_assertion( c_begin - i == 1 ||  c_begin - i == -1);
    CGAL_assertion( i - c_begin == 1 ||  i - c_begin == -1);
    // Check minimal circulator properties.
    i = i.min_circulator();
    C_Circulator j = i;
    CGAL_assertion( j - i == 0);
    j++;
    CGAL_assertion( j - i == 1);
    j++;
    CGAL_assertion( j - i == 2);
    j++;
    CGAL_assertion( j - i == 3);
    j++;
    CGAL_assertion( j - i == 4);
    j++;
    CGAL_assertion( j - i == 0);
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
        CGAL_assertion( l.size() == 10);
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
            CGAL_assertion(1==test_value_type(static_cast< VT* >(0)));
            CGAL_assertion(1==test_distance_type(static_cast< DT* >(0)));
        
            // Default constructor.
            C_Iterator z = C_Iterator();
            z = C_Iterator(); // avoids warning with NDEBUG
            CGAL::Assert_circulator_or_iterator(z);
            // Copy constructor.
            C_Iterator i = c_begin;
        
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
                    CGAL_assertion_code( C_Iterator j =) ++i;
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
                    CGAL_assertion_code( C_Iterator j =) i++;
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
            C_Iterator i = c_end;
            int su = 0;
            int k  = 5;
            do {
                CGAL_assertion_code(C_Iterator j =) --i;
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
                C_Iterator j = i--;
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
        CGAL_assertion( v.size() == 10);
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
            CGAL_assertion(1==test_value_type(static_cast< VT* >(0)));
            CGAL_assertion(1==test_distance_type(static_cast< DT* >(0)));
        
            // Default constructor.
            Circulator z = Circulator();
            z = Circulator(); // avoids warning with NDEBUG
            CGAL::Assert_circulator_or_iterator(z);
            // Copy constructor.
            Circulator i = begin;
        
            // Check general support for circulators and iterators.
            CGAL_assertion( CGAL::is_empty_range( z, z));
            CGAL_assertion( ! CGAL::is_empty_range( i, begin));
        
            int su = 0;
            int k  = 1;
            // Check general loop, pre-increment, dereference.
            if (! CGAL::is_empty_range( i, begin)) {   // superfluous
                do {
                    CGAL_assertion( k == (*i));
                    su += (*i);
                    ++k;
                    CGAL_assertion_code( Circulator j =) ++i;
                    CGAL_assertion( i ==  j);
                    if ( i != begin) {
                        CGAL_assertion( (*i) == (*j));
                    }
                } while (i != begin);  // Inequality and equality checked.
            }
            CGAL_assertion( i == begin);  // Equality checked.
            CGAL_assertion( su == 15);
        
            // Assignment.
            i = begin;
            su = 0;
            k  = 1;
            // Loop with post increment.
            if (! CGAL::is_empty_range( i, begin)) {   // superfluous
                do {
                    CGAL_assertion( k == (*i));
                    su += (*i);
                    ++k;
                    CGAL_assertion_code( Circulator j =) i++;
                    CGAL_assertion(  i !=  j);
                    if ( i != begin) {
                        CGAL_assertion( (*i) == (*j) + 1);
                    }
                } while (i != begin);
            }
            CGAL_assertion( i == begin);
            CGAL_assertion( su == 15);
        }
        { // Open own scope to hide local variables.
            // Change three elements and check post-/pre-increment.
            Circulator i = begin;
            (*i++) = 4;
            CGAL_assertion( 4 == (*begin));
            CGAL_assertion( 2 == (*i));
            (*i++) = 3;
            CGAL_assertion( 3 == (*i));
            (*++i) = 7;
            CGAL_assertion( 7 == (*i));
        
            // Check the setting and reset these elements
            // to their original values.
            i = begin;
            CGAL_assertion( 4 == (*i));
            (*i) = 1;
            i++;
            CGAL_assertion( 3 == (*i));
            (*i++) = 2;
            CGAL_assertion( 3 == (*i));
            i++;
            CGAL_assertion( 7 == (*i));
            (*i) = 4;
        
            // Check the resetting.
            i = begin;
            int k = 1;
            do {
                CGAL_assertion( k == (*i));
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
            CGAL_assertion(1==test_value_type(static_cast< VT* >(0)));
            CGAL_assertion(1==test_distance_type(static_cast< DT* >(0)));
        
            // Default constructor.
            Circulator z = Circulator();
            z = Circulator(); // avoids warning with NDEBUG
            CGAL::Assert_circulator_or_iterator(z);
            // Copy constructor.
            Circulator i = begin;
        
            // Check general support for circulators and iterators.
            CGAL_assertion( CGAL::is_empty_range( z, z));
            CGAL_assertion( ! CGAL::is_empty_range( i, begin));
        
            int su = 0;
            int k  = 1;
            // Check general loop, pre-increment, dereference.
            if (! CGAL::is_empty_range( i, begin)) {   // superfluous
                do {
                    CGAL_assertion( k == (*i));
                    su += (*i);
                    ++k;
                    CGAL_assertion_code( Circulator j =) ++i;
                    CGAL_assertion( i ==  j);
                    if ( i != begin) {
                        CGAL_assertion( (*i) == (*j));
                    }
                } while (i != begin);  // Inequality and equality checked.
            }
            CGAL_assertion( i == begin);  // Equality checked.
            CGAL_assertion( su == 15);
        
            // Assignment.
            i = begin;
            su = 0;
            k  = 1;
            // Loop with post increment.
            if (! CGAL::is_empty_range( i, begin)) {   // superfluous
                do {
                    CGAL_assertion( k == (*i));
                    su += (*i);
                    ++k;
                    CGAL_assertion_code( Circulator j =) i++;
                    CGAL_assertion(  i !=  j);
                    if ( i != begin) {
                        CGAL_assertion( (*i) == (*j) + 1);
                    }
                } while (i != begin);
            }
            CGAL_assertion( i == begin);
            CGAL_assertion( su == 15);
        }
        { // Open own scope to hide local variables.
            CGAL::Assert_is_at_least_bidirectional_category(begin);
            CGAL::Assert_is_at_least_bidirectional_category(begin);
            // Loop backwards and pre-decrement.
            Circulator i = begin;
            int su = 0;
            int k  = 5;
            do {
                CGAL_assertion_code(Circulator j =) --i;
                CGAL_assertion(  i ==  j);
                CGAL_assertion( (*i) == (*j));
                CGAL_assertion( k == (*i));
                su += (*i);
                --k;
            } while (i != begin);
            CGAL_assertion( i == begin);
            CGAL_assertion( su == 15);
        
            // Assignment.
            i = begin;
            su = 0;
            k  = 5;
            // Loop with post-decrement.
            do {
                Circulator j = i--;
                CGAL_assertion(  i !=  j);
                if ( j != begin) {
                    CGAL_assertion( (*i) == (*j) - 1);
                }
                CGAL_assertion( k == (*i));
                su += (*i);
                --k;
            } while (i != begin);
            CGAL_assertion( i == begin);
            CGAL_assertion( su == 15);
        }
        { // Open own scope to hide local variables.
            // Check generally correct parameter properties.
            CGAL::Assert_circulator_or_iterator(begin);
            CGAL::Assert_circulator_or_iterator(begin);
            CGAL::Assert_is_at_least_forward_category(begin);
            CGAL::Assert_is_at_least_forward_category(begin);
            typedef std::iterator_traits< Circulator >::value_type      VT;
            typedef std::iterator_traits< Circulator >::difference_type DT;
            CGAL_assertion(1==test_value_type(static_cast< VT* >(0)));
            CGAL_assertion(1==test_distance_type(static_cast< DT* >(0)));
        
            // Default constructor.
            Circulator z = Circulator();
            z = Circulator(); // avoids warning with NDEBUG
            CGAL::Assert_circulator_or_iterator(z);
            // Copy constructor.
            Circulator i = begin;
        
            // Check general support for circulators and iterators.
            CGAL_assertion( CGAL::is_empty_range( z, z));
            CGAL_assertion( ! CGAL::is_empty_range( i, begin));
        
            int su = 0;
            int k  = 1;
            // Check general loop, pre-increment, dereference.
            if (! CGAL::is_empty_range( i, begin)) {   // superfluous
                do {
                    CGAL_assertion( k == (*i));
                    su += (*i);
                    ++k;
                    CGAL_assertion_code( Circulator j =) ++i;
                    CGAL_assertion( i ==  j);
                    if ( i != begin) {
                        CGAL_assertion( (*i) == (*j));
                    }
                } while (i != begin);  // Inequality and equality checked.
            }
            CGAL_assertion( i == begin);  // Equality checked.
            CGAL_assertion( su == 15);
        
            // Assignment.
            i = begin;
            su = 0;
            k  = 1;
            // Loop with post increment.
            if (! CGAL::is_empty_range( i, begin)) {   // superfluous
                do {
                    CGAL_assertion( k == (*i));
                    su += (*i);
                    ++k;
                    CGAL_assertion_code( Circulator j =) i++;
                    CGAL_assertion(  i !=  j);
                    if ( i != begin) {
                        CGAL_assertion( (*i) == (*j) + 1);
                    }
                } while (i != begin);
            }
            CGAL_assertion( i == begin);
            CGAL_assertion( su == 15);
        }
        { // Open own scope to hide local variables.
            CGAL::Assert_is_at_least_bidirectional_category(begin);
            CGAL::Assert_is_at_least_bidirectional_category(begin);
            // Loop backwards and pre-decrement.
            Circulator i = begin;
            int su = 0;
            int k  = 5;
            do {
                CGAL_assertion_code(Circulator j =) --i;
                CGAL_assertion(  i ==  j);
                CGAL_assertion( (*i) == (*j));
                CGAL_assertion( k == (*i));
                su += (*i);
                --k;
            } while (i != begin);
            CGAL_assertion( i == begin);
            CGAL_assertion( su == 15);
        
            // Assignment.
            i = begin;
            su = 0;
            k  = 5;
            // Loop with post-decrement.
            do {
                Circulator j = i--;
                CGAL_assertion(  i !=  j);
                if ( j != begin) {
                    CGAL_assertion( (*i) == (*j) - 1);
                }
                CGAL_assertion( k == (*i));
                su += (*i);
                --k;
            } while (i != begin);
            CGAL_assertion( i == begin);
            CGAL_assertion( su == 15);
        }
        { // Open own scope to hide local variables.
            CGAL::Assert_is_at_least_random_access_category(begin);
            CGAL::Assert_is_at_least_random_access_category(begin);
            // Random access.
            int k;
            for( k = 0; k < 5; k++) {
                CGAL_assertion( 1+k == begin[k]);
            }
            CGAL_assertion_code(
              int su = begin[0]
                     + begin[1]
                     + begin[2]
                     + begin[3]
                     + begin[4];)
            CGAL_assertion( su == 15);
        
            // Jump around.
            Circulator i = begin;
            i += 3;
            CGAL_assertion( 4 == (*i));
            i -= 2;
            CGAL_assertion( 2 == (*i));
            i += 3;
            CGAL_assertion( 5 == (*i));
            i -= 4;
            CGAL_assertion( 1 == (*i));
            CGAL_assertion( i == begin);
            Circulator j = i + 3;
            CGAL_assertion( 4 == (*j));
            Circulator jj = j - 2;
            CGAL_assertion( 2 == (*jj));
            jj = 4 + jj;
            CGAL_assertion( jj == begin);
            Circulator ij = jj - 5;
            ij = jj - 5; // avoids warning with NDEBUG
            CGAL_assertion( ij == begin);
        
            // Difference test.
            CGAL_assertion( jj - i == 5  ||  jj - i == 0);
            CGAL_assertion( i + (j-i) == j);
            CGAL_assertion( (j-i) + i == j);
        }
        { // Open own scope to hide local variables.
            Circulator i = begin;
            i[2] = 18;
            i[4] = 9;
            i[3] = 12;
            CGAL_assertion( i[2] == 18);
            CGAL_assertion( i[4] == 9);
            CGAL_assertion( i[3] == 12);
            i[2] = 3;
            i[3] = 4;
            i[4] = 5;
            // Check the resetting.
            i = begin;
            int k = 1;
            do {
                CGAL_assertion( k == (*i));
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
            CGAL_assertion( k2 == 1);
            Circulator::pointer         p1 = &k1;
            (void)p1;
            CGAL_assertion( (*p1) == 1);
            k1 = 3;
            CGAL_assertion( k1 == 3);
            CGAL_assertion( k2 == 3);
            CGAL_assertion( (*p1) == 3);
            k1 = 6;
            CGAL_assertion( k1 == 6);
            CGAL_assertion( k2 == 6);
            CGAL_assertion( (*p1) == 6);
            CGAL_assertion_code( Circulator::size_type s = 5;)
            CGAL_assertion( s == 5);
            CGAL_assertion_code(Circulator::difference_type d = -5;)
            CGAL_assertion( d == -5);
        
            // Check tests for empty data structures.
            Circulator z = Circulator();
            z = Circulator(); // avoids warning with NDEBUG
            CGAL_assertion(   z == CGAL_CIRC_NULL);
            CGAL_assertion( ! (z != CGAL_CIRC_NULL));
            Circulator i = begin;
            CGAL_assertion( ! (i == CGAL_CIRC_NULL));
            CGAL_assertion(   i != CGAL_CIRC_NULL);
            CGAL_assertion( i == begin);
            CGAL_assertion( i == begin);
            // Do I reach myself.
            ++i;
            Circulator j = i;
            int k = 0;
            do {
                CGAL_assertion( k < 5);
                ++k;
                ++i;
            } while( i != j);
            CGAL_assertion( k == 5);
        }
        { // Open own scope to hide local variables.
            // Do I reach myself backwards.
            Circulator i = begin;
            ++i;
            Circulator j = i;
            int k = 0;
            do {
                CGAL_assertion( k < 5);
                ++k;
                --i;
            } while( i != j);
            CGAL_assertion( k == 5);
        }
        { // Open own scope to hide local variables.
            Circulator::difference_type d = begin - begin;
            CGAL_assertion( d == 0);
            d = begin - begin;
            CGAL_assertion( d == 0);
            Circulator i = begin + 1;
            CGAL_assertion( begin - i == 1 ||  begin - i == -1);
            CGAL_assertion( i - begin == 1 ||  i - begin == -1);
            // Check minimal circulator properties.
            i = i.min_circulator();
            Circulator j = i;
            CGAL_assertion( j - i == 0);
            j++;
            CGAL_assertion( j - i == 1);
            j++;
            CGAL_assertion( j - i == 2);
            j++;
            CGAL_assertion( j - i == 3);
            j++;
            CGAL_assertion( j - i == 4);
            j++;
            CGAL_assertion( j - i == 0);
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
            CGAL_assertion(1==test_value_type(static_cast< VT* >(0)));
            CGAL_assertion(1==test_distance_type(static_cast< DT* >(0)));
        
            // Default constructor.
            C_Circulator z = C_Circulator();
            z = C_Circulator(); // avoids warning with NDEBUG
            CGAL::Assert_circulator_or_iterator(z);
            // Copy constructor.
            C_Circulator i = c_begin;
        
            // Check general support for circulators and iterators.
            CGAL_assertion( CGAL::is_empty_range( z, z));
            CGAL_assertion( ! CGAL::is_empty_range( i, c_begin));
        
            int su = 0;
            int k  = 1;
            // Check general loop, pre-increment, dereference.
            if (! CGAL::is_empty_range( i, c_begin)) {   // superfluous
                do {
                    CGAL_assertion( k == (*i));
                    su += (*i);
                    ++k;
                    CGAL_assertion_code( C_Circulator j =) ++i;
                    CGAL_assertion( i ==  j);
                    if ( i != c_begin) {
                        CGAL_assertion( (*i) == (*j));
                    }
                } while (i != c_begin);  // Inequality and equality checked.
            }
            CGAL_assertion( i == c_begin);  // Equality checked.
            CGAL_assertion( su == 15);
        
            // Assignment.
            i = c_begin;
            su = 0;
            k  = 1;
            // Loop with post increment.
            if (! CGAL::is_empty_range( i, c_begin)) {   // superfluous
                do {
                    CGAL_assertion( k == (*i));
                    su += (*i);
                    ++k;
                    CGAL_assertion_code( C_Circulator j =) i++;
                    CGAL_assertion(  i !=  j);
                    if ( i != c_begin) {
                        CGAL_assertion( (*i) == (*j) + 1);
                    }
                } while (i != c_begin);
            }
            CGAL_assertion( i == c_begin);
            CGAL_assertion( su == 15);
        }
        { // Open own scope to hide local variables.
            CGAL::Assert_is_at_least_bidirectional_category(c_begin);
            CGAL::Assert_is_at_least_bidirectional_category(c_begin);
            // Loop backwards and pre-decrement.
            C_Circulator i = c_begin;
            int su = 0;
            int k  = 5;
            do {
                CGAL_assertion_code(C_Circulator j =) --i;
                CGAL_assertion(  i ==  j);
                CGAL_assertion( (*i) == (*j));
                CGAL_assertion( k == (*i));
                su += (*i);
                --k;
            } while (i != c_begin);
            CGAL_assertion( i == c_begin);
            CGAL_assertion( su == 15);
        
            // Assignment.
            i = c_begin;
            su = 0;
            k  = 5;
            // Loop with post-decrement.
            do {
                C_Circulator j = i--;
                CGAL_assertion(  i !=  j);
                if ( j != c_begin) {
                    CGAL_assertion( (*i) == (*j) - 1);
                }
                CGAL_assertion( k == (*i));
                su += (*i);
                --k;
            } while (i != c_begin);
            CGAL_assertion( i == c_begin);
            CGAL_assertion( su == 15);
        }
        { // Open own scope to hide local variables.
            CGAL::Assert_is_at_least_random_access_category(c_begin);
            CGAL::Assert_is_at_least_random_access_category(c_begin);
            // Random access.
            int k;
            for( k = 0; k < 5; k++) {
                CGAL_assertion( 1+k == c_begin[k]);
            }
            CGAL_assertion_code(
              int su = c_begin[0]
                     + c_begin[1]
                     + c_begin[2]
                     + c_begin[3]
                     + c_begin[4];)
            CGAL_assertion( su == 15);
        
            // Jump around.
            C_Circulator i = c_begin;
            i += 3;
            CGAL_assertion( 4 == (*i));
            i -= 2;
            CGAL_assertion( 2 == (*i));
            i += 3;
            CGAL_assertion( 5 == (*i));
            i -= 4;
            CGAL_assertion( 1 == (*i));
            CGAL_assertion( i == c_begin);
            C_Circulator j = i + 3;
            CGAL_assertion( 4 == (*j));
            C_Circulator jj = j - 2;
            CGAL_assertion( 2 == (*jj));
            jj = 4 + jj;
            CGAL_assertion( jj == c_begin);
            C_Circulator ij = jj - 5;
            ij = jj - 5; // avoids warning with NDEBUG
            CGAL_assertion( ij == c_begin);
        
            // Difference test.
            CGAL_assertion( jj - i == 5  ||  jj - i == 0);
            CGAL_assertion( i + (j-i) == j);
            CGAL_assertion( (j-i) + i == j);
        }
        { // Open own scope to hide local variables.
            CGAL::Assert_circulator( c_begin);
            CGAL::Assert_circulator( c_begin);
        
            // Check the local type parameters.
            C_Circulator::value_type      k1;
            k1 = 1;
            C_Circulator::reference       k2 = k1;
            (void)k2;
            CGAL_assertion( k2 == 1);
            C_Circulator::pointer         p1 = &k1;
            (void)p1;
            CGAL_assertion( (*p1) == 1);
            k1 = 3;
            CGAL_assertion( k1 == 3);
            CGAL_assertion( k2 == 3);
            CGAL_assertion( (*p1) == 3);
            k1 = 6;
            CGAL_assertion( k1 == 6);
            CGAL_assertion( k2 == 6);
            CGAL_assertion( (*p1) == 6);
            CGAL_assertion_code( C_Circulator::size_type s = 5;)
            CGAL_assertion( s == 5);
            CGAL_assertion_code(C_Circulator::difference_type d = -5;)
            CGAL_assertion( d == -5);
        
            // Check tests for empty data structures.
            C_Circulator z = C_Circulator();
            z = C_Circulator(); // avoids warning with NDEBUG
            CGAL_assertion(   z == CGAL_CIRC_NULL);
            CGAL_assertion( ! (z != CGAL_CIRC_NULL));
            C_Circulator i = c_begin;
            CGAL_assertion( ! (i == CGAL_CIRC_NULL));
            CGAL_assertion(   i != CGAL_CIRC_NULL);
            CGAL_assertion( i == c_begin);
            CGAL_assertion( i == c_begin);
            // Do I reach myself.
            ++i;
            C_Circulator j = i;
            int k = 0;
            do {
                CGAL_assertion( k < 5);
                ++k;
                ++i;
            } while( i != j);
            CGAL_assertion( k == 5);
        }
        { // Open own scope to hide local variables.
            // Do I reach myself backwards.
            C_Circulator i = c_begin;
            ++i;
            C_Circulator j = i;
            int k = 0;
            do {
                CGAL_assertion( k < 5);
                ++k;
                --i;
            } while( i != j);
            CGAL_assertion( k == 5);
        }
        { // Open own scope to hide local variables.
            C_Circulator::difference_type d = c_begin - c_begin;
            CGAL_assertion( d == 0);
            d = c_begin - c_begin;
            CGAL_assertion( d == 0);
            C_Circulator i = c_begin + 1;
            CGAL_assertion( c_begin - i == 1 ||  c_begin - i == -1);
            CGAL_assertion( i - c_begin == 1 ||  i - c_begin == -1);
            // Check minimal circulator properties.
            i = i.min_circulator();
            C_Circulator j = i;
            CGAL_assertion( j - i == 0);
            j++;
            CGAL_assertion( j - i == 1);
            j++;
            CGAL_assertion( j - i == 2);
            j++;
            CGAL_assertion( j - i == 3);
            j++;
            CGAL_assertion( j - i == 4);
            j++;
            CGAL_assertion( j - i == 0);
        }
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
    CGAL_assertion( l.size() == 5);
    typedef List::iterator Iterator;
    Inverse_index<Iterator> index(l.begin(),l.end());
    l.push_back( 6);
    index.push_back( -- (l.end()));
    l.push_back( 7);
    index.push_back( -- (l.end()));
    CGAL_assertion( l.size() == 7);
    for ( Iterator i = l.begin(); i != l.end(); ++i) {
      CGAL_assertion( *i == index[i] + 1);
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
    CGAL_assertion( v.size() == 5);
    typedef Vector::iterator Iterator;
    Inverse_index<Iterator> index(v.begin(),v.end());
    v.push_back( 6);
    index.push_back( v.end() - 1);
    v.push_back( 7);
    index.push_back( v.end() - 1);
    CGAL_assertion( v.size() == 7);
    for ( Iterator i = v.begin(); i != v.end(); ++i) {
      CGAL_assertion( *i == index[i] + 1);
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
    CGAL_assertion( l.size() == 5);
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
    CGAL_assertion( l.size() == 7);
    Circulator c( l.begin(),l.end());
    Circulator d = c;
    do {
      CGAL_assertion( *c == index[c] + 1);
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
    CGAL_assertion( l.size() == 5);
    typedef List::iterator Iterator;
    Random_access_adaptor<Iterator> index(l.begin(),l.end());
    l.push_back( 6);
    index.push_back( -- (l.end()));
    l.push_back( 7);
    index.push_back( -- (l.end()));
    CGAL_assertion( l.size() == 7);
    for ( Iterator i = l.begin(); i != l.end(); ++i) {
      CGAL_assertion( *i == *(index[*i - 1]));
    }
  }
  {
    typedef std::vector<std::size_t> Vector;
    Vector v;
    v.push_back( 1);
    v.push_back( 2);
    v.push_back( 3);
    v.push_back( 4);
    v.push_back( 5);
    CGAL_assertion( v.size() == 5);
    typedef Vector::iterator Iterator;
    Random_access_adaptor<Iterator> index(v.begin(),v.end());
    v.push_back( 6);
    index.push_back( v.end() - 1);
    v.push_back( 7);
    index.push_back( v.end() - 1);
    CGAL_assertion( v.size() == 7);
    for ( Iterator i = v.begin(); i != v.end(); ++i) {
      CGAL_assertion( *i == *(index[*i - 1]));
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
    CGAL_assertion( l.size() == 5);
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
    CGAL_assertion( l.size() == 7);
    Circulator c( l.begin(),l.end());
    Circulator d = c;
    do {
      CGAL_assertion( *c == *(index[*c - 1]));
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
    CGAL_assertion( l.size() == 5);
    typedef List::iterator Iterator;
    Random_access_value_adaptor<Iterator,std::size_t>
    index(l.begin(),l.end());
  l.push_back( 6);
  index.push_back( -- (l.end()));
  l.push_back( 7);
  index.push_back( -- (l.end()));
  CGAL_assertion( l.size() == 7);
  for ( Iterator i = l.begin(); i != l.end(); ++i) {
    CGAL_assertion( *i == index[*i - 1]);
  }
  }
}

int main() {
  init_global_data();
  test_Iterator_identity();
  test_Circulator_identity();
  test_Iterator_project();
  test_Circulator_project();
  test_Circulator_on_node();
  test_N_step_adaptor();
  test_N_step_adaptor_derived();
  test_Inverse_index();
  test_Random_access_adaptor();
  clean_global_data();
  return 0;
}
// EOF //
