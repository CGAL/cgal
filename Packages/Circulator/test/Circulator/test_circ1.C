// ============================================================================
//
// Copyright (c) 1997 The CGAL Consortium
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
// file          : test_circ1.C
// chapter       : $CGAL_Chapter: Circulators $
// package       : $CGAL_Package: Circulator 3.4 (02 Sep 1999) $
// source        : circulator.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : INRIA, Sophia Antipolis
//
// Test circulator support and adaptor between circulator and iterator.
// ============================================================================


#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif // CGAL_BASIC_H
#ifndef CGAL_PROTECT_CSTDDEF
#include <cstddef>
#define CGAL_PROTECT_CSTDDEF
#endif
#ifndef CGAL_PROTECT_ITERATOR
#include <iterator>
#define CGAL_PROTECT_ITERATOR
#endif
#ifndef CGAL_PROTECT_LIST
#include <list>
#define CGAL_PROTECT_LIST
#endif
#ifndef CGAL_PROTECT_VECTOR
#include <vector>
#define CGAL_PROTECT_VECTOR
#endif
#ifndef CGAL_CIRCULATOR_H
#include <CGAL/circulator.h>
#endif // CGAL_CIRCULATOR_H
// needed for test data structures
#ifndef CGAL_CIRCULATOR_IMPL_H
#include <CGAL/circulator_impl.h>
#endif // CGAL_CIRCULATOR_IMPL_H

using namespace CGAL;

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
    Node( Node* _nx, Node* _pv, int n)
        : key(n), next(_nx), prev( _pv) {}
};
Node* new_node( Node* _nx, Node* _pv, int n) {
    return new Node( _nx, _pv, n);
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
    CNode* _next;
    CNode* _prev;
  public:
    int   key;
    CNode*       next()       { return _next;}
    const CNode* next() const { return _next;}
    CNode*       prev()       { return _prev;}
    const CNode* prev() const { return _prev;}
    CNode() : _next(this), _prev( this), key(0) {}
    CNode( int n) : _next(this), _prev( this), key(n) {}
    CNode( CNode* _nx, CNode* _pv, int n)
        : _next(_nx), _prev( _pv), key(n) {}
    friend CNode* new_cnode( CNode* _nx, CNode* _pv, int n);
    friend void append_cnode( CNode* p, int n);
    friend void delete_cnodes( CNode* p);
};
CNode* new_cnode( CNode* _nx, CNode* _pv, int n) {
    return new CNode( _nx, _pv, n);
}
void append_cnode( CNode* p, int n) {
    CNode* q = new_cnode( p, p->_prev, n);
    p->_prev->_next = q;
    p->_prev = q;
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
    p = p->_next;
    while ( p != end) {
        delete q;
        q = p;
        p = p->_next;
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

int test_difference_type( std::ptrdiff_t*)   { return 1;}
int test_difference_type( char*)             { return 2;}
int test_difference_type( double*)           { return 3;}
template< class T> inline
int foo3( T t, Forward_circulator_tag) {
    Assert_circulator( t);
    Assert_forward_category( t);
    return 1;
}

template< class T> inline
int foo3( T t,  Bidirectional_circulator_tag) {
    Assert_circulator( t);
    Assert_bidirectional_category( t);
    return 2;
}

template< class T> inline
int foo3( T t,  Random_access_circulator_tag) {
    Assert_circulator( t);
    Assert_random_access_category( t);
    return 3;
}

template< class T> inline
int foo3( T t) { return -1; }  // never used

template< class T> inline
int foo2( T t, Circulator_tag) {
    Assert_circulator( t);
    typedef std::iterator_traits<T> Traits;
    typedef typename Traits::iterator_category iterator_category;
    return foo3( t, iterator_category());
}

template< class T> inline
int foo2( T t, Iterator_tag) {
    Assert_iterator( t);
    return 4;
}

template< class T> inline
int foo2( T t) { return -1; }  // never used

template< class T> inline
int foo( T t) { return foo2( t, query_circulator_or_iterator( t));}

int bar( std::size_t)    { return 1;}
int bar( std::ptrdiff_t) { return 2;}
int bar( char)           { return 3;}
int bar( double)         { return 4;}

void test_tags() {
  {
    std::list<int>   l;
    CGAL_assertion( 4 == foo( l.begin()));
    CGAL_assertion( 4 == foo( l.end()));
    std::vector<int> v;
    CGAL_assertion( 4 == foo( v.begin()));
    CGAL_assertion( 4 == foo( v.end()));

    int* p = NULL;
    CGAL_assertion( 4 == foo( p));
    {
        typedef Forward_circulator_base<char,std::ptrdiff_t,std::size_t>
            FC;
        typedef Bidirectional_circulator_base<double,char,std::size_t>
            BC;
        typedef Random_access_circulator_base<std::size_t,char,std::size_t>
            RC;

        FC f_c = FC();
        BC b_c = BC();
        RC r_c = RC();
        CGAL_assertion( 1 == foo( f_c));
        CGAL_assertion( 2 == foo( b_c));
        CGAL_assertion( 3 == foo( r_c));
        CGAL_assertion( 3 == bar( std::iterator_traits<FC>::value_type()));
        CGAL_assertion( 4 == bar( std::iterator_traits<BC>::value_type()));
        CGAL_assertion( 1 == bar( std::iterator_traits<RC>::value_type()));
        CGAL_assertion( 2 == bar( std::iterator_traits<FC>::
                                      difference_type()));
        CGAL_assertion( 3 == bar( std::iterator_traits<BC>::
                                      difference_type()));
        CGAL_assertion( 3 == bar( std::iterator_traits<RC>::
                                      difference_type()));
    }
    {
        typedef Forward_circulator_ptrbase<char,std::ptrdiff_t,std::size_t>
            FC;
        typedef Bidirectional_circulator_ptrbase<double,char,std::size_t>
            BC;
        typedef Random_access_circulator_ptrbase<std::size_t,char,
            std::size_t>   RC;

        FC f_c = FC();
        BC b_c = BC();
        RC r_c = RC();
        CGAL_assertion( 1 == foo( f_c));
        CGAL_assertion( 2 == foo( b_c));
        CGAL_assertion( 3 == foo( r_c));
        CGAL_assertion( 3 == bar( std::iterator_traits<FC>::value_type()));
        CGAL_assertion( 4 == bar( std::iterator_traits<BC>::value_type()));
        CGAL_assertion( 1 == bar( std::iterator_traits<RC>::value_type()));
        CGAL_assertion( 2 == bar( std::iterator_traits<FC>::
                                      difference_type()));
        CGAL_assertion( 3 == bar( std::iterator_traits<BC>::
                                      difference_type()));
        CGAL_assertion( 3 == bar( std::iterator_traits<RC>::
                                      difference_type()));
    }
    // a bit more complicated cases.
    typedef Random_access_circulator_from_iterator<
        std::vector<int>::iterator,
        int,
        std::vector<int>::size_type,
        std::vector<int>::difference_type
    > Circulator;
    typedef Random_access_container_from_circulator<Circulator> Container;
    typedef Container::iterator Iterator;
    Circulator ci( v.begin(), v.end());
    Container  Co( ci);
    CGAL_assertion( 3 == foo( ci));
    CGAL_assertion( 4 == foo( Co.begin()));
    CGAL_assertion( 4 == foo( Co.end()));
  }
}
void test_functions_for_circulators() {
  Node* data_struct = generate_nodes( 5);
  Struct_circulator           start1(data_struct);
  Struct_const_circulator     start2(data_struct);
  Struct_bi_circulator        start3(data_struct);
  Struct_bi_const_circulator  start4(data_struct);
  CGAL_assertion( circulator_size(start1) == 5);
  CGAL_assertion( circulator_size(start2) == 5);
  CGAL_assertion( circulator_size(start3) == 5);
  CGAL_assertion( circulator_size(start4) == 5);
  Struct_circulator           start1a(data_struct);
  Struct_const_circulator     start2a(data_struct);
  Struct_bi_circulator        start3a(data_struct);
  Struct_bi_const_circulator  start4a(data_struct);
  CGAL_assertion( circulator_distance(start1, start1a) == 5);
  CGAL_assertion( circulator_distance(start2, start2a) == 5);
  CGAL_assertion( circulator_distance(start3, start3a) == 5);
  CGAL_assertion( circulator_distance(start4, start4a) == 5);
  CGAL_assertion( iterator_distance(start1, start1a) == 5);
  CGAL_assertion( iterator_distance(start2, start2a) == 5);
  CGAL_assertion( iterator_distance(start3, start3a) == 5);
  CGAL_assertion( iterator_distance(start4, start4a) == 5);
  ++ start1a;
  ++ start2a;
  ++ start3a;
  ++ start4a;
  CGAL_assertion( circulator_distance(start1, start1a) == 1);
  CGAL_assertion( circulator_distance(start2, start2a) == 1);
  CGAL_assertion( circulator_distance(start3, start3a) == 1);
  CGAL_assertion( circulator_distance(start4, start4a) == 1);
  CGAL_assertion( iterator_distance(start1, start1a) == 1);
  CGAL_assertion( iterator_distance(start2, start2a) == 1);
  CGAL_assertion( iterator_distance(start3, start3a) == 1);
  CGAL_assertion( iterator_distance(start4, start4a) == 1);
  ++ start1a;
  ++ start2a;
  ++ start1a;
  ++ start2a;
  ++ start1a;
  ++ start2a;
  -- start3a;
  -- start4a;
  -- start3a;
  -- start4a;
  CGAL_assertion( circulator_distance(start1, start1a) == 4);
  CGAL_assertion( circulator_distance(start2, start2a) == 4);
  CGAL_assertion( circulator_distance(start3, start3a) == 4);
  CGAL_assertion( circulator_distance(start4, start4a) == 4);
  CGAL_assertion( iterator_distance(start1, start1a) == 4);
  CGAL_assertion( iterator_distance(start2, start2a) == 4);
  CGAL_assertion( iterator_distance(start3, start3a) == 4);
  CGAL_assertion( iterator_distance(start4, start4a) == 4);
  delete_nodes(data_struct);
  std::list<int> l;
  CGAL_assertion( iterator_distance(l.begin(), l.end()) == 0);
  l.push_back(3);
  CGAL_assertion( iterator_distance(l.begin(), l.end()) == 1);
  int* my_C_array  = new_C_array();
  CGAL_assertion( iterator_distance(my_C_array, my_C_array+5) == 5);
  typedef Random_access_circulator_from_iterator<
      int*, int, std::size_t, std::ptrdiff_t> Circulator1;
  Circulator1 c1( my_C_array, my_C_array+5);
  CGAL_assertion( circulator_size(c1) == 5);
  Circulator1 c3 = c1;
  ++c3;
  CGAL_assertion( c1 == c3.min_circulator());
  --c3;
  --c3;
  CGAL_assertion( c1 == c3.min_circulator());
  --c3;
  CGAL_assertion( c1 == c3.min_circulator());
  typedef Random_access_const_circulator_from_iterator<
      int*, int, std::size_t, std::ptrdiff_t> Circulator2;
  Circulator2 c2( my_C_array, my_C_array+5);
  CGAL_assertion( circulator_size(c2) == 5);
  c3 = c1;
  Circulator2 c4(c2);
  CGAL_assertion( circulator_distance(c1,c3) == 5);
  CGAL_assertion( circulator_distance(c2,c4) == 5);
  c3 ++;
  c4 ++;
  CGAL_assertion( circulator_distance(c1,c3) == 1);
  CGAL_assertion( circulator_distance(c2,c4) == 1);
  c3 --;
  c4 --;
  c3 --;
  c4 --;
  CGAL_assertion( circulator_distance(c1,c3) == 4);
  CGAL_assertion( circulator_distance(c2,c4) == 4);
  delete[] my_C_array;
  CGAL_assertion( 2 == non_negative_mod( -4, 3));
  CGAL_assertion( 0 == non_negative_mod( -3, 3));
  CGAL_assertion( 1 == non_negative_mod( -2, 3));
  CGAL_assertion( 2 == non_negative_mod( -1, 3));
  CGAL_assertion( 0 == non_negative_mod(  0, 3));
  CGAL_assertion( 1 == non_negative_mod(  1, 3));
  CGAL_assertion( 2 == non_negative_mod(  2, 3));
  CGAL_assertion( 0 == non_negative_mod(  3, 3));
  CGAL_assertion( 1 == non_negative_mod(  4, 3));
}
void test_iterator_and_circulators() {
    std::vector<int> v;
    CGAL_assertion( is_empty_range( v.begin(), v.end()));
    typedef Random_access_circulator_from_iterator<
        std::vector<int>::iterator,
        int,
        std::vector<int>::size_type,
        std::vector<int>::difference_type
    > Circulator;
    Circulator c( v.begin(), v.end());
    CGAL_assertion( is_empty_range( c, c));
    v.push_back( 5);
    CGAL_assertion( is_empty_range( v.begin(), v.begin()));
    CGAL_assertion( ! is_empty_range( v.begin(), v.end()));
    Circulator d( v.begin(), v.end());
    CGAL_assertion( ! is_empty_range( d, d));
    std::vector<int>::iterator i = v.begin();
    int j = 0;
    {
        CGAL_assertion( v.size() == 1);
        CGAL_For_all( i, v.end()) {
            CGAL_assertion( *i == 5);
            j++;
        }
        CGAL_assertion( j == 1);
    }{
        CGAL_For_all( d, d) {
            CGAL_assertion( *d == 5);
            j++;
        }
        CGAL_assertion( j == 2);
    }{
        CGAL_For_all_backwards( d, d) {
            CGAL_assertion( *d == 5);
            j++;
        }
        CGAL_assertion( j == 3);
    }
}
void test_container_from_circulator() {
  Node* data_struct = generate_nodes( 5);
  {
    Struct_circulator start(data_struct);
    typedef Forward_container_from_circulator<Struct_circulator>
        Container;
    Container X;
    CGAL_assertion( X.begin() == X.end());
    Container Y(X);
    CGAL_assertion( Y.begin() == Y.end());
    Container C( start);
    typedef Container::iterator Iterator;
    Iterator begin = C.begin();
    Iterator end   = C.end();
    Assert_forward_category(begin);
    Assert_forward_category(end);
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_circulator_or_iterator(end);
        CGAL::Assert_is_at_least_forward_category(begin);
        CGAL::Assert_is_at_least_forward_category(end);
        typedef std::iterator_traits<Iterator> I_Traits;
        typedef I_Traits::value_type      I_value_type;
        typedef I_Traits::difference_type I_difference_type;
        CGAL_assertion(1==test_value_type( (I_value_type*)(0)));
        CGAL_assertion(1==test_difference_type( (I_difference_type*)(0)));
    
        // Default constructor.
        Iterator z = Iterator();
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
                Iterator j = ++i;
                CGAL_assertion(  i ==  j);
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
                Iterator j = i++;
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
    CGAL::Assert_iterator( begin);
    CGAL::Assert_iterator( end);
    begin++;
    CGAL_assertion( (*(begin.current_circulator())).key == 2);
  }
  {
    Struct_const_circulator start(data_struct);
    typedef Forward_container_from_circulator<Struct_const_circulator>
        Container;
    const Container X;
    CGAL_assertion( X.begin() == X.end());
    const Container Y(X);
    CGAL_assertion( Y.begin() == Y.end());
    const Container C( start);
    typedef Container::const_iterator Const_iterator;
    Const_iterator begin = C.begin();
    Const_iterator end   = C.end();
    Assert_forward_category(begin);
    Assert_forward_category(end);
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_circulator_or_iterator(end);
        CGAL::Assert_is_at_least_forward_category(begin);
        CGAL::Assert_is_at_least_forward_category(end);
        typedef std::iterator_traits<Const_iterator> I_Traits;
        typedef I_Traits::value_type      I_value_type;
        typedef I_Traits::difference_type I_difference_type;
        CGAL_assertion(1==test_value_type( (I_value_type*)(0)));
        CGAL_assertion(1==test_difference_type( (I_difference_type*)(0)));
    
        // Default constructor.
        Const_iterator z = Const_iterator();
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Const_iterator i = begin;
    
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
                Const_iterator j = ++i;
                CGAL_assertion(  i ==  j);
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
                Const_iterator j = i++;
                CGAL_assertion(  i !=  j);
                if ( i != end) {
                    CGAL_assertion( (*i).key == (*j).key + 1);
                }
            } while (i != end);
        }
        CGAL_assertion( i == end);
        CGAL_assertion( su == 15);
    }
    CGAL::Assert_iterator( begin);
    CGAL::Assert_iterator( end);
    begin++;
    CGAL_assertion( (*(begin.current_circulator())).key == 2);
  }
  {
    Struct_bi_circulator start(data_struct);
    typedef Bidirectional_container_from_circulator<
        Struct_bi_circulator> Container;
    Container X;
    CGAL_assertion( X.begin() == X.end());
    Container Y(X);
    CGAL_assertion( Y.begin() == Y.end());
    Container C( start);
    typedef Container::iterator Iterator;
    Iterator begin = C.begin();
    Iterator end   = C.end();
    Assert_bidirectional_category(begin);
    Assert_bidirectional_category(end);
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_circulator_or_iterator(end);
        CGAL::Assert_is_at_least_forward_category(begin);
        CGAL::Assert_is_at_least_forward_category(end);
        typedef std::iterator_traits<Iterator> I_Traits;
        typedef I_Traits::value_type      I_value_type;
        typedef I_Traits::difference_type I_difference_type;
        CGAL_assertion(1==test_value_type( (I_value_type*)(0)));
        CGAL_assertion(1==test_difference_type( (I_difference_type*)(0)));
    
        // Default constructor.
        Iterator z = Iterator();
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
                Iterator j = ++i;
                CGAL_assertion(  i ==  j);
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
                Iterator j = i++;
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
        typedef std::iterator_traits<Iterator> I_Traits;
        typedef I_Traits::value_type      I_value_type;
        typedef I_Traits::difference_type I_difference_type;
        CGAL_assertion(1==test_value_type( (I_value_type*)(0)));
        CGAL_assertion(1==test_difference_type( (I_difference_type*)(0)));
    
        // Default constructor.
        Iterator z = Iterator();
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
                Iterator j = ++i;
                CGAL_assertion(  i ==  j);
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
                Iterator j = i++;
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
            Iterator j = --i;
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
    begin++;
    CGAL_assertion( (*(begin.current_circulator())).key == 2);
  }
  {
    Struct_bi_const_circulator start(data_struct);
    typedef Bidirectional_container_from_circulator<
        Struct_bi_const_circulator> Container;
    const Container X;
    CGAL_assertion( X.begin() == X.end());
    const Container Y(X);
    CGAL_assertion( Y.begin() == Y.end());
    const Container C( start);
    typedef Container::const_iterator Const_iterator;
    Const_iterator begin = C.begin();
    Const_iterator end   = C.end();
    Assert_bidirectional_category(begin);
    Assert_bidirectional_category(end);
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_circulator_or_iterator(end);
        CGAL::Assert_is_at_least_forward_category(begin);
        CGAL::Assert_is_at_least_forward_category(end);
        typedef std::iterator_traits<Const_iterator> I_Traits;
        typedef I_Traits::value_type      I_value_type;
        typedef I_Traits::difference_type I_difference_type;
        CGAL_assertion(1==test_value_type( (I_value_type*)(0)));
        CGAL_assertion(1==test_difference_type( (I_difference_type*)(0)));
    
        // Default constructor.
        Const_iterator z = Const_iterator();
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Const_iterator i = begin;
    
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
                Const_iterator j = ++i;
                CGAL_assertion(  i ==  j);
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
                Const_iterator j = i++;
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
        Const_iterator i = end;
        int su = 0;
        int k  = 5;
        do {
            Const_iterator j = --i;
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
            Const_iterator j = i--;
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
    begin++;
    CGAL_assertion( (*(begin.current_circulator())).key == 2);
  }
  int* my_C_array  = new_C_array();
  {
    typedef Random_access_circulator_from_iterator<
        int*, int, std::size_t, std::ptrdiff_t> Circulator;
    Circulator c( my_C_array, my_C_array+5);
    typedef Random_access_container_from_circulator<Circulator>
        Container;
    Container X;
    CGAL_assertion( X.begin() == X.end());
    Container Y(X);
    CGAL_assertion( Y.begin() == Y.end());
    Container C( c);
    typedef Container::iterator Iterator;
    Iterator begin = C.begin();
    Iterator end   = C.end();
    Assert_random_access_category(begin);
    Assert_random_access_category(end);
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_circulator_or_iterator(end);
        CGAL::Assert_is_at_least_forward_category(begin);
        CGAL::Assert_is_at_least_forward_category(end);
        typedef std::iterator_traits<Iterator> I_Traits;
        typedef I_Traits::value_type      I_value_type;
        typedef I_Traits::difference_type I_difference_type;
        CGAL_assertion(1==test_value_type( (I_value_type*)(0)));
        CGAL_assertion(1==test_difference_type( (I_difference_type*)(0)));
    
        // Default constructor.
        Iterator z = Iterator();
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
                Iterator j = ++i;
                CGAL_assertion(  i ==  j);
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
                Iterator j = i++;
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
        typedef std::iterator_traits<Iterator> I_Traits;
        typedef I_Traits::value_type      I_value_type;
        typedef I_Traits::difference_type I_difference_type;
        CGAL_assertion(1==test_value_type( (I_value_type*)(0)));
        CGAL_assertion(1==test_difference_type( (I_difference_type*)(0)));
    
        // Default constructor.
        Iterator z = Iterator();
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
                Iterator j = ++i;
                CGAL_assertion(  i ==  j);
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
                Iterator j = i++;
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
            Iterator j = --i;
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
        typedef std::iterator_traits<Iterator> I_Traits;
        typedef I_Traits::value_type      I_value_type;
        typedef I_Traits::difference_type I_difference_type;
        CGAL_assertion(1==test_value_type( (I_value_type*)(0)));
        CGAL_assertion(1==test_difference_type( (I_difference_type*)(0)));
    
        // Default constructor.
        Iterator z = Iterator();
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
                Iterator j = ++i;
                CGAL_assertion(  i ==  j);
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
                Iterator j = i++;
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
            Iterator j = --i;
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
        int su = begin[0]
               + begin[1]
               + begin[2]
               + begin[3]
               + begin[4];
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
    begin++;
    CGAL_assertion( (*(begin.current_circulator())) == 2);
  }
  {
    typedef Random_access_const_circulator_from_iterator<
        int*, int, std::size_t, std::ptrdiff_t> Circulator;
    Circulator c( my_C_array, my_C_array+5);
    typedef Random_access_container_from_circulator<Circulator>
        Container;
    const Container X;
    CGAL_assertion( X.begin() == X.end());
    const Container Y(X);
    CGAL_assertion( Y.begin() == Y.end());
    const Container C( c);
    typedef Container::const_iterator Const_iterator;
    Const_iterator begin = C.begin();
    Const_iterator end   = C.end();
    Assert_random_access_category(begin);
    Assert_random_access_category(end);
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_circulator_or_iterator(end);
        CGAL::Assert_is_at_least_forward_category(begin);
        CGAL::Assert_is_at_least_forward_category(end);
        typedef std::iterator_traits<Const_iterator> I_Traits;
        typedef I_Traits::value_type      I_value_type;
        typedef I_Traits::difference_type I_difference_type;
        CGAL_assertion(1==test_value_type( (I_value_type*)(0)));
        CGAL_assertion(1==test_difference_type( (I_difference_type*)(0)));
    
        // Default constructor.
        Const_iterator z = Const_iterator();
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Const_iterator i = begin;
    
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
                Const_iterator j = ++i;
                CGAL_assertion(  i ==  j);
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
                Const_iterator j = i++;
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
        Const_iterator i = end;
        int su = 0;
        int k  = 5;
        do {
            Const_iterator j = --i;
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
            Const_iterator j = i--;
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
        int su = begin[0]
               + begin[1]
               + begin[2]
               + begin[3]
               + begin[4];
        CGAL_assertion( su == 15);
    
        // Jump around.
        Const_iterator i = begin;
        i += 3;
        CGAL_assertion( 4 == (*i));
        i -= 2;
        CGAL_assertion( 2 == (*i));
        i += 3;
        CGAL_assertion( 5 == (*i));
        i -= 4;
        CGAL_assertion( 1 == (*i));
        CGAL_assertion( i == begin);
        Const_iterator j = i + 3;
        CGAL_assertion( 4 == (*j));
        Const_iterator jj = j - 2;
        CGAL_assertion( 2 == (*jj));
        jj = 4 + jj;
        CGAL_assertion( jj == end);
        Const_iterator ij = jj - 5;
        CGAL_assertion( ij == begin);
    
        // Difference test.
        CGAL_assertion( jj - i == 5  ||  jj - i == 0);
        CGAL_assertion( i + (j-i) == j);
        CGAL_assertion( (j-i) + i == j);
    }
    CGAL::Assert_iterator( begin);
    CGAL::Assert_iterator( end);
    { // Open own scope to hide local variables.
        CGAL_assertion( end - begin ==  5);
        CGAL_assertion( begin - end == -5);
        // Relational operator.
        Const_iterator i = begin;
        ++i;
        Const_iterator j = i;
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
    begin++;
    CGAL_assertion( (*(begin.current_circulator())) == 2);
  }
  {
    Struct_bi_circulator start(data_struct);
    typedef Container_from_circulator< Struct_bi_circulator> Container;
    Container X;
    CGAL_assertion( X.begin() == X.end());
    Container Y(X);
    CGAL_assertion( Y.begin() == Y.end());
    Container C( start);
    typedef Container::iterator Iterator;
    Iterator begin = C.begin();
    Iterator end   = C.end();
    Assert_bidirectional_category(begin);
    Assert_bidirectional_category(end);
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_circulator_or_iterator(end);
        CGAL::Assert_is_at_least_forward_category(begin);
        CGAL::Assert_is_at_least_forward_category(end);
        typedef std::iterator_traits<Iterator> I_Traits;
        typedef I_Traits::value_type      I_value_type;
        typedef I_Traits::difference_type I_difference_type;
        CGAL_assertion(1==test_value_type( (I_value_type*)(0)));
        CGAL_assertion(1==test_difference_type( (I_difference_type*)(0)));
    
        // Default constructor.
        Iterator z = Iterator();
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
                Iterator j = ++i;
                CGAL_assertion(  i ==  j);
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
                Iterator j = i++;
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
        typedef std::iterator_traits<Iterator> I_Traits;
        typedef I_Traits::value_type      I_value_type;
        typedef I_Traits::difference_type I_difference_type;
        CGAL_assertion(1==test_value_type( (I_value_type*)(0)));
        CGAL_assertion(1==test_difference_type( (I_difference_type*)(0)));
    
        // Default constructor.
        Iterator z = Iterator();
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
                Iterator j = ++i;
                CGAL_assertion(  i ==  j);
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
                Iterator j = i++;
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
            Iterator j = --i;
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
    begin++;
    CGAL_assertion( (*(begin.current_circulator())).key == 2);
  }
  {
    Struct_bi_const_circulator start(data_struct);
    typedef Container_from_circulator<
        Struct_bi_const_circulator> Container;
    const Container X;
    CGAL_assertion( X.begin() == X.end());
    const Container Y(X);
    CGAL_assertion( Y.begin() == Y.end());
    const Container C( start);
    typedef Container::const_iterator Const_iterator;
    Const_iterator begin = C.begin();
    Const_iterator end   = C.end();
    Assert_bidirectional_category(begin);
    Assert_bidirectional_category(end);
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_circulator_or_iterator(end);
        CGAL::Assert_is_at_least_forward_category(begin);
        CGAL::Assert_is_at_least_forward_category(end);
        typedef std::iterator_traits<Const_iterator> I_Traits;
        typedef I_Traits::value_type      I_value_type;
        typedef I_Traits::difference_type I_difference_type;
        CGAL_assertion(1==test_value_type( (I_value_type*)(0)));
        CGAL_assertion(1==test_difference_type( (I_difference_type*)(0)));
    
        // Default constructor.
        Const_iterator z = Const_iterator();
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Const_iterator i = begin;
    
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
                Const_iterator j = ++i;
                CGAL_assertion(  i ==  j);
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
                Const_iterator j = i++;
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
        Const_iterator i = end;
        int su = 0;
        int k  = 5;
        do {
            Const_iterator j = --i;
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
            Const_iterator j = i--;
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
    begin++;
    CGAL_assertion( (*(begin.current_circulator())).key == 2);
  }
  {
    typedef Random_access_circulator_from_iterator<
        int*, int, std::size_t, std::ptrdiff_t> Circulator;
    Circulator c( my_C_array, my_C_array+5);
    typedef Container_from_circulator<Circulator>
        Container;
    Container X;
    CGAL_assertion( X.begin() == X.end());
    Container Y(X);
    CGAL_assertion( Y.begin() == Y.end());
    Container C( c);
    typedef Container::iterator Iterator;
    Iterator begin = C.begin();
    Iterator end   = C.end();
    Assert_random_access_category(begin);
    Assert_random_access_category(end);
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(begin);
        CGAL::Assert_circulator_or_iterator(end);
        CGAL::Assert_is_at_least_forward_category(begin);
        CGAL::Assert_is_at_least_forward_category(end);
        typedef std::iterator_traits<Iterator> I_Traits;
        typedef I_Traits::value_type      I_value_type;
        typedef I_Traits::difference_type I_difference_type;
        CGAL_assertion(1==test_value_type( (I_value_type*)(0)));
        CGAL_assertion(1==test_difference_type( (I_difference_type*)(0)));
    
        // Default constructor.
        Iterator z = Iterator();
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
                Iterator j = ++i;
                CGAL_assertion(  i ==  j);
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
                Iterator j = i++;
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
        typedef std::iterator_traits<Iterator> I_Traits;
        typedef I_Traits::value_type      I_value_type;
        typedef I_Traits::difference_type I_difference_type;
        CGAL_assertion(1==test_value_type( (I_value_type*)(0)));
        CGAL_assertion(1==test_difference_type( (I_difference_type*)(0)));
    
        // Default constructor.
        Iterator z = Iterator();
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
                Iterator j = ++i;
                CGAL_assertion(  i ==  j);
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
                Iterator j = i++;
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
            Iterator j = --i;
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
        typedef std::iterator_traits<Iterator> I_Traits;
        typedef I_Traits::value_type      I_value_type;
        typedef I_Traits::difference_type I_difference_type;
        CGAL_assertion(1==test_value_type( (I_value_type*)(0)));
        CGAL_assertion(1==test_difference_type( (I_difference_type*)(0)));
    
        // Default constructor.
        Iterator z = Iterator();
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
                Iterator j = ++i;
                CGAL_assertion(  i ==  j);
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
                Iterator j = i++;
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
            Iterator j = --i;
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
        int su = begin[0]
               + begin[1]
               + begin[2]
               + begin[3]
               + begin[4];
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
    begin++;
    CGAL_assertion( (*(begin.current_circulator())) == 2);
  }
  delete_nodes(data_struct);
  delete[] my_C_array;
}
#ifndef CGAL_PROTECT_ALGORITHM
#include <algorithm>
#define CGAL_PROTECT_ALGORITHM
#endif

void test_circulator_from_iterator() {
  int* my_C_array  = new_C_array();
  {
    typedef Forward_circulator_from_iterator<
        int*, int, std::size_t, std::ptrdiff_t> Circulator;
    Circulator c( my_C_array, my_C_array+5);
    Assert_forward_category(c);
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(c);
        CGAL::Assert_circulator_or_iterator(c);
        CGAL::Assert_is_at_least_forward_category(c);
        CGAL::Assert_is_at_least_forward_category(c);
        typedef std::iterator_traits<Circulator> I_Traits;
        typedef I_Traits::value_type      I_value_type;
        typedef I_Traits::difference_type I_difference_type;
        CGAL_assertion(1==test_value_type( (I_value_type*)(0)));
        CGAL_assertion(1==test_difference_type( (I_difference_type*)(0)));
    
        // Default constructor.
        Circulator z = Circulator();
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Circulator i = c;
    
        // Check general support for circulators and iterators.
        CGAL_assertion( CGAL::is_empty_range( z, z));
        CGAL_assertion( ! CGAL::is_empty_range( i, c));
    
        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, c)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                Circulator j = ++i;
                CGAL_assertion(  i ==  j);
                if ( i != c) {
                    CGAL_assertion( (*i) == (*j));
                }
            } while (i != c);  // Inequality and equality checked.
        }
        CGAL_assertion( i == c);  // Equality checked.
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = c;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, c)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                Circulator j = i++;
                CGAL_assertion(  i !=  j);
                if ( i != c) {
                    CGAL_assertion( (*i) == (*j) + 1);
                }
            } while (i != c);
        }
        CGAL_assertion( i == c);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        // Change three elements and check post-/pre-increment.
        Circulator i = c;
        (*i++) = 4;
        CGAL_assertion( 4 == (*c));
        CGAL_assertion( 2 == (*i));
        (*i++) = 3;
        CGAL_assertion( 3 == (*i));
        (*++i) = 7;
        CGAL_assertion( 7 == (*i));
    
        // Check the setting and reset these elements
        // to their original values.
        i = c;
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
        i = c;
        int k = 1;
        do {
            CGAL_assertion( k == (*i));
            ++i;
            ++k;
        } while (i != c);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_circulator( c);
        CGAL::Assert_circulator( c);
    
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
        Circulator::size_type s = 5;
        CGAL_assertion( s == 5);
        Circulator::difference_type d = -5;
        CGAL_assertion( d == -5);
    
        // Check tests for empty data structures.
        Circulator z = Circulator();
        CGAL_assertion(   z == CGAL_CIRC_NULL);
        CGAL_assertion( ! (z != CGAL_CIRC_NULL));
        Circulator i = c;
        CGAL_assertion( ! (i == CGAL_CIRC_NULL));
        CGAL_assertion(   i != CGAL_CIRC_NULL);
        CGAL_assertion( i == c);
        CGAL_assertion( i == c);
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
    c++;
    CGAL_assertion( (*(c.current_iterator())) == 2);
  }{
    typedef Forward_const_circulator_from_iterator<
        int*, int, std::size_t, std::ptrdiff_t> Circulator;
    Circulator c( my_C_array, my_C_array+5);
    Assert_forward_category(c);
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(c);
        CGAL::Assert_circulator_or_iterator(c);
        CGAL::Assert_is_at_least_forward_category(c);
        CGAL::Assert_is_at_least_forward_category(c);
        typedef std::iterator_traits<Circulator> I_Traits;
        typedef I_Traits::value_type      I_value_type;
        typedef I_Traits::difference_type I_difference_type;
        CGAL_assertion(1==test_value_type( (I_value_type*)(0)));
        CGAL_assertion(1==test_difference_type( (I_difference_type*)(0)));
    
        // Default constructor.
        Circulator z = Circulator();
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Circulator i = c;
    
        // Check general support for circulators and iterators.
        CGAL_assertion( CGAL::is_empty_range( z, z));
        CGAL_assertion( ! CGAL::is_empty_range( i, c));
    
        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, c)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                Circulator j = ++i;
                CGAL_assertion(  i ==  j);
                if ( i != c) {
                    CGAL_assertion( (*i) == (*j));
                }
            } while (i != c);  // Inequality and equality checked.
        }
        CGAL_assertion( i == c);  // Equality checked.
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = c;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, c)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                Circulator j = i++;
                CGAL_assertion(  i !=  j);
                if ( i != c) {
                    CGAL_assertion( (*i) == (*j) + 1);
                }
            } while (i != c);
        }
        CGAL_assertion( i == c);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_circulator( c);
        CGAL::Assert_circulator( c);
    
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
        Circulator::size_type s = 5;
        CGAL_assertion( s == 5);
        Circulator::difference_type d = -5;
        CGAL_assertion( d == -5);
    
        // Check tests for empty data structures.
        Circulator z = Circulator();
        CGAL_assertion(   z == CGAL_CIRC_NULL);
        CGAL_assertion( ! (z != CGAL_CIRC_NULL));
        Circulator i = c;
        CGAL_assertion( ! (i == CGAL_CIRC_NULL));
        CGAL_assertion(   i != CGAL_CIRC_NULL);
        CGAL_assertion( i == c);
        CGAL_assertion( i == c);
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
    c++;
    CGAL_assertion( (*(c.current_iterator())) == 2);
  }{
    typedef Bidirectional_circulator_from_iterator<
        int*, int, std::size_t, std::ptrdiff_t> Circulator;
    Circulator c( my_C_array, my_C_array+5);
    Assert_bidirectional_category(c);
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(c);
        CGAL::Assert_circulator_or_iterator(c);
        CGAL::Assert_is_at_least_forward_category(c);
        CGAL::Assert_is_at_least_forward_category(c);
        typedef std::iterator_traits<Circulator> I_Traits;
        typedef I_Traits::value_type      I_value_type;
        typedef I_Traits::difference_type I_difference_type;
        CGAL_assertion(1==test_value_type( (I_value_type*)(0)));
        CGAL_assertion(1==test_difference_type( (I_difference_type*)(0)));
    
        // Default constructor.
        Circulator z = Circulator();
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Circulator i = c;
    
        // Check general support for circulators and iterators.
        CGAL_assertion( CGAL::is_empty_range( z, z));
        CGAL_assertion( ! CGAL::is_empty_range( i, c));
    
        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, c)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                Circulator j = ++i;
                CGAL_assertion(  i ==  j);
                if ( i != c) {
                    CGAL_assertion( (*i) == (*j));
                }
            } while (i != c);  // Inequality and equality checked.
        }
        CGAL_assertion( i == c);  // Equality checked.
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = c;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, c)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                Circulator j = i++;
                CGAL_assertion(  i !=  j);
                if ( i != c) {
                    CGAL_assertion( (*i) == (*j) + 1);
                }
            } while (i != c);
        }
        CGAL_assertion( i == c);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        // Change three elements and check post-/pre-increment.
        Circulator i = c;
        (*i++) = 4;
        CGAL_assertion( 4 == (*c));
        CGAL_assertion( 2 == (*i));
        (*i++) = 3;
        CGAL_assertion( 3 == (*i));
        (*++i) = 7;
        CGAL_assertion( 7 == (*i));
    
        // Check the setting and reset these elements
        // to their original values.
        i = c;
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
        i = c;
        int k = 1;
        do {
            CGAL_assertion( k == (*i));
            ++i;
            ++k;
        } while (i != c);
    }
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(c);
        CGAL::Assert_circulator_or_iterator(c);
        CGAL::Assert_is_at_least_forward_category(c);
        CGAL::Assert_is_at_least_forward_category(c);
        typedef std::iterator_traits<Circulator> I_Traits;
        typedef I_Traits::value_type      I_value_type;
        typedef I_Traits::difference_type I_difference_type;
        CGAL_assertion(1==test_value_type( (I_value_type*)(0)));
        CGAL_assertion(1==test_difference_type( (I_difference_type*)(0)));
    
        // Default constructor.
        Circulator z = Circulator();
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Circulator i = c;
    
        // Check general support for circulators and iterators.
        CGAL_assertion( CGAL::is_empty_range( z, z));
        CGAL_assertion( ! CGAL::is_empty_range( i, c));
    
        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, c)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                Circulator j = ++i;
                CGAL_assertion(  i ==  j);
                if ( i != c) {
                    CGAL_assertion( (*i) == (*j));
                }
            } while (i != c);  // Inequality and equality checked.
        }
        CGAL_assertion( i == c);  // Equality checked.
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = c;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, c)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                Circulator j = i++;
                CGAL_assertion(  i !=  j);
                if ( i != c) {
                    CGAL_assertion( (*i) == (*j) + 1);
                }
            } while (i != c);
        }
        CGAL_assertion( i == c);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_bidirectional_category(c);
        CGAL::Assert_is_at_least_bidirectional_category(c);
        // Loop backwards and pre-decrement.
        Circulator i = c;
        int su = 0;
        int k  = 5;
        do {
            Circulator j = --i;
            CGAL_assertion(  i ==  j);
            CGAL_assertion( (*i) == (*j));
            CGAL_assertion( k == (*i));
            su += (*i);
            --k;
        } while (i != c);
        CGAL_assertion( i == c);
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = c;
        su = 0;
        k  = 5;
        // Loop with post-decrement.
        do {
            Circulator j = i--;
            CGAL_assertion(  i !=  j);
            if ( j != c) {
                CGAL_assertion( (*i) == (*j) - 1);
            }
            CGAL_assertion( k == (*i));
            su += (*i);
            --k;
        } while (i != c);
        CGAL_assertion( i == c);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_circulator( c);
        CGAL::Assert_circulator( c);
    
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
        Circulator::size_type s = 5;
        CGAL_assertion( s == 5);
        Circulator::difference_type d = -5;
        CGAL_assertion( d == -5);
    
        // Check tests for empty data structures.
        Circulator z = Circulator();
        CGAL_assertion(   z == CGAL_CIRC_NULL);
        CGAL_assertion( ! (z != CGAL_CIRC_NULL));
        Circulator i = c;
        CGAL_assertion( ! (i == CGAL_CIRC_NULL));
        CGAL_assertion(   i != CGAL_CIRC_NULL);
        CGAL_assertion( i == c);
        CGAL_assertion( i == c);
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
        Circulator i = c;
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
    c++;
    CGAL_assertion( (*(c.current_iterator())) == 2);
  }{
    typedef Bidirectional_const_circulator_from_iterator<
        int*, int, std::size_t, std::ptrdiff_t> Circulator;
    Circulator c( my_C_array, my_C_array+5);
    Assert_bidirectional_category(c);
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(c);
        CGAL::Assert_circulator_or_iterator(c);
        CGAL::Assert_is_at_least_forward_category(c);
        CGAL::Assert_is_at_least_forward_category(c);
        typedef std::iterator_traits<Circulator> I_Traits;
        typedef I_Traits::value_type      I_value_type;
        typedef I_Traits::difference_type I_difference_type;
        CGAL_assertion(1==test_value_type( (I_value_type*)(0)));
        CGAL_assertion(1==test_difference_type( (I_difference_type*)(0)));
    
        // Default constructor.
        Circulator z = Circulator();
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Circulator i = c;
    
        // Check general support for circulators and iterators.
        CGAL_assertion( CGAL::is_empty_range( z, z));
        CGAL_assertion( ! CGAL::is_empty_range( i, c));
    
        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, c)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                Circulator j = ++i;
                CGAL_assertion(  i ==  j);
                if ( i != c) {
                    CGAL_assertion( (*i) == (*j));
                }
            } while (i != c);  // Inequality and equality checked.
        }
        CGAL_assertion( i == c);  // Equality checked.
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = c;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, c)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                Circulator j = i++;
                CGAL_assertion(  i !=  j);
                if ( i != c) {
                    CGAL_assertion( (*i) == (*j) + 1);
                }
            } while (i != c);
        }
        CGAL_assertion( i == c);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_bidirectional_category(c);
        CGAL::Assert_is_at_least_bidirectional_category(c);
        // Loop backwards and pre-decrement.
        Circulator i = c;
        int su = 0;
        int k  = 5;
        do {
            Circulator j = --i;
            CGAL_assertion(  i ==  j);
            CGAL_assertion( (*i) == (*j));
            CGAL_assertion( k == (*i));
            su += (*i);
            --k;
        } while (i != c);
        CGAL_assertion( i == c);
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = c;
        su = 0;
        k  = 5;
        // Loop with post-decrement.
        do {
            Circulator j = i--;
            CGAL_assertion(  i !=  j);
            if ( j != c) {
                CGAL_assertion( (*i) == (*j) - 1);
            }
            CGAL_assertion( k == (*i));
            su += (*i);
            --k;
        } while (i != c);
        CGAL_assertion( i == c);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_circulator( c);
        CGAL::Assert_circulator( c);
    
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
        Circulator::size_type s = 5;
        CGAL_assertion( s == 5);
        Circulator::difference_type d = -5;
        CGAL_assertion( d == -5);
    
        // Check tests for empty data structures.
        Circulator z = Circulator();
        CGAL_assertion(   z == CGAL_CIRC_NULL);
        CGAL_assertion( ! (z != CGAL_CIRC_NULL));
        Circulator i = c;
        CGAL_assertion( ! (i == CGAL_CIRC_NULL));
        CGAL_assertion(   i != CGAL_CIRC_NULL);
        CGAL_assertion( i == c);
        CGAL_assertion( i == c);
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
        Circulator i = c;
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
    c++;
    CGAL_assertion( (*(c.current_iterator())) == 2);
  }{
    typedef Random_access_circulator_from_iterator<
        int*, int, std::size_t, std::ptrdiff_t> Circulator;
    Circulator c( my_C_array, my_C_array+5);
    Assert_random_access_category(c);
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(c);
        CGAL::Assert_circulator_or_iterator(c);
        CGAL::Assert_is_at_least_forward_category(c);
        CGAL::Assert_is_at_least_forward_category(c);
        typedef std::iterator_traits<Circulator> I_Traits;
        typedef I_Traits::value_type      I_value_type;
        typedef I_Traits::difference_type I_difference_type;
        CGAL_assertion(1==test_value_type( (I_value_type*)(0)));
        CGAL_assertion(1==test_difference_type( (I_difference_type*)(0)));
    
        // Default constructor.
        Circulator z = Circulator();
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Circulator i = c;
    
        // Check general support for circulators and iterators.
        CGAL_assertion( CGAL::is_empty_range( z, z));
        CGAL_assertion( ! CGAL::is_empty_range( i, c));
    
        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, c)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                Circulator j = ++i;
                CGAL_assertion(  i ==  j);
                if ( i != c) {
                    CGAL_assertion( (*i) == (*j));
                }
            } while (i != c);  // Inequality and equality checked.
        }
        CGAL_assertion( i == c);  // Equality checked.
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = c;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, c)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                Circulator j = i++;
                CGAL_assertion(  i !=  j);
                if ( i != c) {
                    CGAL_assertion( (*i) == (*j) + 1);
                }
            } while (i != c);
        }
        CGAL_assertion( i == c);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        // Change three elements and check post-/pre-increment.
        Circulator i = c;
        (*i++) = 4;
        CGAL_assertion( 4 == (*c));
        CGAL_assertion( 2 == (*i));
        (*i++) = 3;
        CGAL_assertion( 3 == (*i));
        (*++i) = 7;
        CGAL_assertion( 7 == (*i));
    
        // Check the setting and reset these elements
        // to their original values.
        i = c;
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
        i = c;
        int k = 1;
        do {
            CGAL_assertion( k == (*i));
            ++i;
            ++k;
        } while (i != c);
    }
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(c);
        CGAL::Assert_circulator_or_iterator(c);
        CGAL::Assert_is_at_least_forward_category(c);
        CGAL::Assert_is_at_least_forward_category(c);
        typedef std::iterator_traits<Circulator> I_Traits;
        typedef I_Traits::value_type      I_value_type;
        typedef I_Traits::difference_type I_difference_type;
        CGAL_assertion(1==test_value_type( (I_value_type*)(0)));
        CGAL_assertion(1==test_difference_type( (I_difference_type*)(0)));
    
        // Default constructor.
        Circulator z = Circulator();
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Circulator i = c;
    
        // Check general support for circulators and iterators.
        CGAL_assertion( CGAL::is_empty_range( z, z));
        CGAL_assertion( ! CGAL::is_empty_range( i, c));
    
        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, c)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                Circulator j = ++i;
                CGAL_assertion(  i ==  j);
                if ( i != c) {
                    CGAL_assertion( (*i) == (*j));
                }
            } while (i != c);  // Inequality and equality checked.
        }
        CGAL_assertion( i == c);  // Equality checked.
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = c;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, c)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                Circulator j = i++;
                CGAL_assertion(  i !=  j);
                if ( i != c) {
                    CGAL_assertion( (*i) == (*j) + 1);
                }
            } while (i != c);
        }
        CGAL_assertion( i == c);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_bidirectional_category(c);
        CGAL::Assert_is_at_least_bidirectional_category(c);
        // Loop backwards and pre-decrement.
        Circulator i = c;
        int su = 0;
        int k  = 5;
        do {
            Circulator j = --i;
            CGAL_assertion(  i ==  j);
            CGAL_assertion( (*i) == (*j));
            CGAL_assertion( k == (*i));
            su += (*i);
            --k;
        } while (i != c);
        CGAL_assertion( i == c);
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = c;
        su = 0;
        k  = 5;
        // Loop with post-decrement.
        do {
            Circulator j = i--;
            CGAL_assertion(  i !=  j);
            if ( j != c) {
                CGAL_assertion( (*i) == (*j) - 1);
            }
            CGAL_assertion( k == (*i));
            su += (*i);
            --k;
        } while (i != c);
        CGAL_assertion( i == c);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(c);
        CGAL::Assert_circulator_or_iterator(c);
        CGAL::Assert_is_at_least_forward_category(c);
        CGAL::Assert_is_at_least_forward_category(c);
        typedef std::iterator_traits<Circulator> I_Traits;
        typedef I_Traits::value_type      I_value_type;
        typedef I_Traits::difference_type I_difference_type;
        CGAL_assertion(1==test_value_type( (I_value_type*)(0)));
        CGAL_assertion(1==test_difference_type( (I_difference_type*)(0)));
    
        // Default constructor.
        Circulator z = Circulator();
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Circulator i = c;
    
        // Check general support for circulators and iterators.
        CGAL_assertion( CGAL::is_empty_range( z, z));
        CGAL_assertion( ! CGAL::is_empty_range( i, c));
    
        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, c)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                Circulator j = ++i;
                CGAL_assertion(  i ==  j);
                if ( i != c) {
                    CGAL_assertion( (*i) == (*j));
                }
            } while (i != c);  // Inequality and equality checked.
        }
        CGAL_assertion( i == c);  // Equality checked.
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = c;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, c)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                Circulator j = i++;
                CGAL_assertion(  i !=  j);
                if ( i != c) {
                    CGAL_assertion( (*i) == (*j) + 1);
                }
            } while (i != c);
        }
        CGAL_assertion( i == c);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_bidirectional_category(c);
        CGAL::Assert_is_at_least_bidirectional_category(c);
        // Loop backwards and pre-decrement.
        Circulator i = c;
        int su = 0;
        int k  = 5;
        do {
            Circulator j = --i;
            CGAL_assertion(  i ==  j);
            CGAL_assertion( (*i) == (*j));
            CGAL_assertion( k == (*i));
            su += (*i);
            --k;
        } while (i != c);
        CGAL_assertion( i == c);
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = c;
        su = 0;
        k  = 5;
        // Loop with post-decrement.
        do {
            Circulator j = i--;
            CGAL_assertion(  i !=  j);
            if ( j != c) {
                CGAL_assertion( (*i) == (*j) - 1);
            }
            CGAL_assertion( k == (*i));
            su += (*i);
            --k;
        } while (i != c);
        CGAL_assertion( i == c);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_random_access_category(c);
        CGAL::Assert_is_at_least_random_access_category(c);
        // Random access.
        int k;
        for( k = 0; k < 5; k++) {
            CGAL_assertion( 1+k == c[k]);
        }
        int su = c[0]
               + c[1]
               + c[2]
               + c[3]
               + c[4];
        CGAL_assertion( su == 15);
    
        // Jump around.
        Circulator i = c;
        i += 3;
        CGAL_assertion( 4 == (*i));
        i -= 2;
        CGAL_assertion( 2 == (*i));
        i += 3;
        CGAL_assertion( 5 == (*i));
        i -= 4;
        CGAL_assertion( 1 == (*i));
        CGAL_assertion( i == c);
        Circulator j = i + 3;
        CGAL_assertion( 4 == (*j));
        Circulator jj = j - 2;
        CGAL_assertion( 2 == (*jj));
        jj = 4 + jj;
        CGAL_assertion( jj == c);
        Circulator ij = jj - 5;
        CGAL_assertion( ij == c);
    
        // Difference test.
        CGAL_assertion( jj - i == 5  ||  jj - i == 0);
        CGAL_assertion( i + (j-i) == j);
        CGAL_assertion( (j-i) + i == j);
    }
    { // Open own scope to hide local variables.
        Circulator i = c;
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
        i = c;
        int k = 1;
        do {
            CGAL_assertion( k == (*i));
            ++i;
            ++k;
        } while (i != c);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_circulator( c);
        CGAL::Assert_circulator( c);
    
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
        Circulator::size_type s = 5;
        CGAL_assertion( s == 5);
        Circulator::difference_type d = -5;
        CGAL_assertion( d == -5);
    
        // Check tests for empty data structures.
        Circulator z = Circulator();
        CGAL_assertion(   z == CGAL_CIRC_NULL);
        CGAL_assertion( ! (z != CGAL_CIRC_NULL));
        Circulator i = c;
        CGAL_assertion( ! (i == CGAL_CIRC_NULL));
        CGAL_assertion(   i != CGAL_CIRC_NULL);
        CGAL_assertion( i == c);
        CGAL_assertion( i == c);
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
        Circulator i = c;
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
        Circulator::difference_type d = c - c;
        CGAL_assertion( d == 0);
        d = c - c;
        CGAL_assertion( d == 0);
        Circulator i = c + 1;
        CGAL_assertion( c - i == 1 ||  c - i == -1);
        CGAL_assertion( i - c == 1 ||  i - c == -1);
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
    c++;
    CGAL_assertion( (*(c.current_iterator())) == 2);
  }{
    typedef Random_access_const_circulator_from_iterator<
        int*, int, std::size_t, std::ptrdiff_t> Circulator;
    Circulator c( my_C_array, my_C_array+5);
    Assert_random_access_category(c);
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(c);
        CGAL::Assert_circulator_or_iterator(c);
        CGAL::Assert_is_at_least_forward_category(c);
        CGAL::Assert_is_at_least_forward_category(c);
        typedef std::iterator_traits<Circulator> I_Traits;
        typedef I_Traits::value_type      I_value_type;
        typedef I_Traits::difference_type I_difference_type;
        CGAL_assertion(1==test_value_type( (I_value_type*)(0)));
        CGAL_assertion(1==test_difference_type( (I_difference_type*)(0)));
    
        // Default constructor.
        Circulator z = Circulator();
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Circulator i = c;
    
        // Check general support for circulators and iterators.
        CGAL_assertion( CGAL::is_empty_range( z, z));
        CGAL_assertion( ! CGAL::is_empty_range( i, c));
    
        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, c)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                Circulator j = ++i;
                CGAL_assertion(  i ==  j);
                if ( i != c) {
                    CGAL_assertion( (*i) == (*j));
                }
            } while (i != c);  // Inequality and equality checked.
        }
        CGAL_assertion( i == c);  // Equality checked.
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = c;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, c)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                Circulator j = i++;
                CGAL_assertion(  i !=  j);
                if ( i != c) {
                    CGAL_assertion( (*i) == (*j) + 1);
                }
            } while (i != c);
        }
        CGAL_assertion( i == c);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_bidirectional_category(c);
        CGAL::Assert_is_at_least_bidirectional_category(c);
        // Loop backwards and pre-decrement.
        Circulator i = c;
        int su = 0;
        int k  = 5;
        do {
            Circulator j = --i;
            CGAL_assertion(  i ==  j);
            CGAL_assertion( (*i) == (*j));
            CGAL_assertion( k == (*i));
            su += (*i);
            --k;
        } while (i != c);
        CGAL_assertion( i == c);
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = c;
        su = 0;
        k  = 5;
        // Loop with post-decrement.
        do {
            Circulator j = i--;
            CGAL_assertion(  i !=  j);
            if ( j != c) {
                CGAL_assertion( (*i) == (*j) - 1);
            }
            CGAL_assertion( k == (*i));
            su += (*i);
            --k;
        } while (i != c);
        CGAL_assertion( i == c);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_random_access_category(c);
        CGAL::Assert_is_at_least_random_access_category(c);
        // Random access.
        int k;
        for( k = 0; k < 5; k++) {
            CGAL_assertion( 1+k == c[k]);
        }
        int su = c[0]
               + c[1]
               + c[2]
               + c[3]
               + c[4];
        CGAL_assertion( su == 15);
    
        // Jump around.
        Circulator i = c;
        i += 3;
        CGAL_assertion( 4 == (*i));
        i -= 2;
        CGAL_assertion( 2 == (*i));
        i += 3;
        CGAL_assertion( 5 == (*i));
        i -= 4;
        CGAL_assertion( 1 == (*i));
        CGAL_assertion( i == c);
        Circulator j = i + 3;
        CGAL_assertion( 4 == (*j));
        Circulator jj = j - 2;
        CGAL_assertion( 2 == (*jj));
        jj = 4 + jj;
        CGAL_assertion( jj == c);
        Circulator ij = jj - 5;
        CGAL_assertion( ij == c);
    
        // Difference test.
        CGAL_assertion( jj - i == 5  ||  jj - i == 0);
        CGAL_assertion( i + (j-i) == j);
        CGAL_assertion( (j-i) + i == j);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_circulator( c);
        CGAL::Assert_circulator( c);
    
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
        Circulator::size_type s = 5;
        CGAL_assertion( s == 5);
        Circulator::difference_type d = -5;
        CGAL_assertion( d == -5);
    
        // Check tests for empty data structures.
        Circulator z = Circulator();
        CGAL_assertion(   z == CGAL_CIRC_NULL);
        CGAL_assertion( ! (z != CGAL_CIRC_NULL));
        Circulator i = c;
        CGAL_assertion( ! (i == CGAL_CIRC_NULL));
        CGAL_assertion(   i != CGAL_CIRC_NULL);
        CGAL_assertion( i == c);
        CGAL_assertion( i == c);
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
        Circulator i = c;
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
        Circulator::difference_type d = c - c;
        CGAL_assertion( d == 0);
        d = c - c;
        CGAL_assertion( d == 0);
        Circulator i = c + 1;
        CGAL_assertion( c - i == 1 ||  c - i == -1);
        CGAL_assertion( i - c == 1 ||  i - c == -1);
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
    c++;
    CGAL_assertion( (*(c.current_iterator())) == 2);
  }
  delete[] my_C_array;
  {
    // the example program `circulator_prog1.C'
    typedef Random_access_circulator_from_iterator<
        std::vector<int>::iterator,
        std::vector<int>::value_type,
        std::vector<int>::size_type,
        std::vector<int>::difference_type
    > Circulator;
    typedef Random_access_container_from_circulator<Circulator>
          Container;
    std::vector<int> W;
    W.push_back(5);
    W.push_back(2);
    W.push_back(9);
    Circulator ci( W.begin(), W.end());
    Container  Co( ci);
    std::sort( Co.begin(), Co.end());
    CGAL_assertion( W.begin()[0] == 2);
    CGAL_assertion( W.begin()[1] == 5);
    CGAL_assertion( W.begin()[2] == 9);
  }
}
#ifndef CGAL_PROTECT_ALGORITHM
#include <algorithm>
#define CGAL_PROTECT_ALGORITHM
#endif

void test_circulator_from_container() {
  {
    typedef Forward_circulator_from_container< std::vector<int> >
        Circulator;
    Circulator c( &V);
    Assert_forward_category(c);
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(c);
        CGAL::Assert_circulator_or_iterator(c);
        CGAL::Assert_is_at_least_forward_category(c);
        CGAL::Assert_is_at_least_forward_category(c);
        typedef std::iterator_traits<Circulator> I_Traits;
        typedef I_Traits::value_type      I_value_type;
        typedef I_Traits::difference_type I_difference_type;
        CGAL_assertion(1==test_value_type( (I_value_type*)(0)));
        CGAL_assertion(1==test_difference_type( (I_difference_type*)(0)));
    
        // Default constructor.
        Circulator z = Circulator();
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Circulator i = c;
    
        // Check general support for circulators and iterators.
        CGAL_assertion( CGAL::is_empty_range( z, z));
        CGAL_assertion( ! CGAL::is_empty_range( i, c));
    
        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, c)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                Circulator j = ++i;
                CGAL_assertion(  i ==  j);
                if ( i != c) {
                    CGAL_assertion( (*i) == (*j));
                }
            } while (i != c);  // Inequality and equality checked.
        }
        CGAL_assertion( i == c);  // Equality checked.
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = c;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, c)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                Circulator j = i++;
                CGAL_assertion(  i !=  j);
                if ( i != c) {
                    CGAL_assertion( (*i) == (*j) + 1);
                }
            } while (i != c);
        }
        CGAL_assertion( i == c);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        // Change three elements and check post-/pre-increment.
        Circulator i = c;
        (*i++) = 4;
        CGAL_assertion( 4 == (*c));
        CGAL_assertion( 2 == (*i));
        (*i++) = 3;
        CGAL_assertion( 3 == (*i));
        (*++i) = 7;
        CGAL_assertion( 7 == (*i));
    
        // Check the setting and reset these elements
        // to their original values.
        i = c;
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
        i = c;
        int k = 1;
        do {
            CGAL_assertion( k == (*i));
            ++i;
            ++k;
        } while (i != c);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_circulator( c);
        CGAL::Assert_circulator( c);
    
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
        Circulator::size_type s = 5;
        CGAL_assertion( s == 5);
        Circulator::difference_type d = -5;
        CGAL_assertion( d == -5);
    
        // Check tests for empty data structures.
        Circulator z = Circulator();
        CGAL_assertion(   z == CGAL_CIRC_NULL);
        CGAL_assertion( ! (z != CGAL_CIRC_NULL));
        Circulator i = c;
        CGAL_assertion( ! (i == CGAL_CIRC_NULL));
        CGAL_assertion(   i != CGAL_CIRC_NULL);
        CGAL_assertion( i == c);
        CGAL_assertion( i == c);
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
    c++;
    CGAL_assertion( (*(c.current_iterator())) == 2);
  }{
    typedef Forward_const_circulator_from_container<std::vector<int> >
        Circulator;
    Circulator c( &V);
    Assert_forward_category(c);
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(c);
        CGAL::Assert_circulator_or_iterator(c);
        CGAL::Assert_is_at_least_forward_category(c);
        CGAL::Assert_is_at_least_forward_category(c);
        typedef std::iterator_traits<Circulator> I_Traits;
        typedef I_Traits::value_type      I_value_type;
        typedef I_Traits::difference_type I_difference_type;
        CGAL_assertion(1==test_value_type( (I_value_type*)(0)));
        CGAL_assertion(1==test_difference_type( (I_difference_type*)(0)));
    
        // Default constructor.
        Circulator z = Circulator();
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Circulator i = c;
    
        // Check general support for circulators and iterators.
        CGAL_assertion( CGAL::is_empty_range( z, z));
        CGAL_assertion( ! CGAL::is_empty_range( i, c));
    
        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, c)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                Circulator j = ++i;
                CGAL_assertion(  i ==  j);
                if ( i != c) {
                    CGAL_assertion( (*i) == (*j));
                }
            } while (i != c);  // Inequality and equality checked.
        }
        CGAL_assertion( i == c);  // Equality checked.
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = c;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, c)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                Circulator j = i++;
                CGAL_assertion(  i !=  j);
                if ( i != c) {
                    CGAL_assertion( (*i) == (*j) + 1);
                }
            } while (i != c);
        }
        CGAL_assertion( i == c);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_circulator( c);
        CGAL::Assert_circulator( c);
    
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
        Circulator::size_type s = 5;
        CGAL_assertion( s == 5);
        Circulator::difference_type d = -5;
        CGAL_assertion( d == -5);
    
        // Check tests for empty data structures.
        Circulator z = Circulator();
        CGAL_assertion(   z == CGAL_CIRC_NULL);
        CGAL_assertion( ! (z != CGAL_CIRC_NULL));
        Circulator i = c;
        CGAL_assertion( ! (i == CGAL_CIRC_NULL));
        CGAL_assertion(   i != CGAL_CIRC_NULL);
        CGAL_assertion( i == c);
        CGAL_assertion( i == c);
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
    c++;
    CGAL_assertion( (*(c.current_iterator())) == 2);
  }{
    typedef Bidirectional_circulator_from_container< std::vector<int> >
        Circulator;
    Circulator c( &V);
    Assert_bidirectional_category(c);
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(c);
        CGAL::Assert_circulator_or_iterator(c);
        CGAL::Assert_is_at_least_forward_category(c);
        CGAL::Assert_is_at_least_forward_category(c);
        typedef std::iterator_traits<Circulator> I_Traits;
        typedef I_Traits::value_type      I_value_type;
        typedef I_Traits::difference_type I_difference_type;
        CGAL_assertion(1==test_value_type( (I_value_type*)(0)));
        CGAL_assertion(1==test_difference_type( (I_difference_type*)(0)));
    
        // Default constructor.
        Circulator z = Circulator();
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Circulator i = c;
    
        // Check general support for circulators and iterators.
        CGAL_assertion( CGAL::is_empty_range( z, z));
        CGAL_assertion( ! CGAL::is_empty_range( i, c));
    
        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, c)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                Circulator j = ++i;
                CGAL_assertion(  i ==  j);
                if ( i != c) {
                    CGAL_assertion( (*i) == (*j));
                }
            } while (i != c);  // Inequality and equality checked.
        }
        CGAL_assertion( i == c);  // Equality checked.
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = c;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, c)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                Circulator j = i++;
                CGAL_assertion(  i !=  j);
                if ( i != c) {
                    CGAL_assertion( (*i) == (*j) + 1);
                }
            } while (i != c);
        }
        CGAL_assertion( i == c);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        // Change three elements and check post-/pre-increment.
        Circulator i = c;
        (*i++) = 4;
        CGAL_assertion( 4 == (*c));
        CGAL_assertion( 2 == (*i));
        (*i++) = 3;
        CGAL_assertion( 3 == (*i));
        (*++i) = 7;
        CGAL_assertion( 7 == (*i));
    
        // Check the setting and reset these elements
        // to their original values.
        i = c;
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
        i = c;
        int k = 1;
        do {
            CGAL_assertion( k == (*i));
            ++i;
            ++k;
        } while (i != c);
    }
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(c);
        CGAL::Assert_circulator_or_iterator(c);
        CGAL::Assert_is_at_least_forward_category(c);
        CGAL::Assert_is_at_least_forward_category(c);
        typedef std::iterator_traits<Circulator> I_Traits;
        typedef I_Traits::value_type      I_value_type;
        typedef I_Traits::difference_type I_difference_type;
        CGAL_assertion(1==test_value_type( (I_value_type*)(0)));
        CGAL_assertion(1==test_difference_type( (I_difference_type*)(0)));
    
        // Default constructor.
        Circulator z = Circulator();
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Circulator i = c;
    
        // Check general support for circulators and iterators.
        CGAL_assertion( CGAL::is_empty_range( z, z));
        CGAL_assertion( ! CGAL::is_empty_range( i, c));
    
        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, c)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                Circulator j = ++i;
                CGAL_assertion(  i ==  j);
                if ( i != c) {
                    CGAL_assertion( (*i) == (*j));
                }
            } while (i != c);  // Inequality and equality checked.
        }
        CGAL_assertion( i == c);  // Equality checked.
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = c;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, c)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                Circulator j = i++;
                CGAL_assertion(  i !=  j);
                if ( i != c) {
                    CGAL_assertion( (*i) == (*j) + 1);
                }
            } while (i != c);
        }
        CGAL_assertion( i == c);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_bidirectional_category(c);
        CGAL::Assert_is_at_least_bidirectional_category(c);
        // Loop backwards and pre-decrement.
        Circulator i = c;
        int su = 0;
        int k  = 5;
        do {
            Circulator j = --i;
            CGAL_assertion(  i ==  j);
            CGAL_assertion( (*i) == (*j));
            CGAL_assertion( k == (*i));
            su += (*i);
            --k;
        } while (i != c);
        CGAL_assertion( i == c);
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = c;
        su = 0;
        k  = 5;
        // Loop with post-decrement.
        do {
            Circulator j = i--;
            CGAL_assertion(  i !=  j);
            if ( j != c) {
                CGAL_assertion( (*i) == (*j) - 1);
            }
            CGAL_assertion( k == (*i));
            su += (*i);
            --k;
        } while (i != c);
        CGAL_assertion( i == c);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_circulator( c);
        CGAL::Assert_circulator( c);
    
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
        Circulator::size_type s = 5;
        CGAL_assertion( s == 5);
        Circulator::difference_type d = -5;
        CGAL_assertion( d == -5);
    
        // Check tests for empty data structures.
        Circulator z = Circulator();
        CGAL_assertion(   z == CGAL_CIRC_NULL);
        CGAL_assertion( ! (z != CGAL_CIRC_NULL));
        Circulator i = c;
        CGAL_assertion( ! (i == CGAL_CIRC_NULL));
        CGAL_assertion(   i != CGAL_CIRC_NULL);
        CGAL_assertion( i == c);
        CGAL_assertion( i == c);
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
        Circulator i = c;
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
    c++;
    CGAL_assertion( (*(c.current_iterator())) == 2);
  }{
    typedef
        Bidirectional_const_circulator_from_container<std::vector<int> >
        Circulator;
    Circulator c( &V);
    Assert_bidirectional_category(c);
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(c);
        CGAL::Assert_circulator_or_iterator(c);
        CGAL::Assert_is_at_least_forward_category(c);
        CGAL::Assert_is_at_least_forward_category(c);
        typedef std::iterator_traits<Circulator> I_Traits;
        typedef I_Traits::value_type      I_value_type;
        typedef I_Traits::difference_type I_difference_type;
        CGAL_assertion(1==test_value_type( (I_value_type*)(0)));
        CGAL_assertion(1==test_difference_type( (I_difference_type*)(0)));
    
        // Default constructor.
        Circulator z = Circulator();
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Circulator i = c;
    
        // Check general support for circulators and iterators.
        CGAL_assertion( CGAL::is_empty_range( z, z));
        CGAL_assertion( ! CGAL::is_empty_range( i, c));
    
        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, c)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                Circulator j = ++i;
                CGAL_assertion(  i ==  j);
                if ( i != c) {
                    CGAL_assertion( (*i) == (*j));
                }
            } while (i != c);  // Inequality and equality checked.
        }
        CGAL_assertion( i == c);  // Equality checked.
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = c;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, c)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                Circulator j = i++;
                CGAL_assertion(  i !=  j);
                if ( i != c) {
                    CGAL_assertion( (*i) == (*j) + 1);
                }
            } while (i != c);
        }
        CGAL_assertion( i == c);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_bidirectional_category(c);
        CGAL::Assert_is_at_least_bidirectional_category(c);
        // Loop backwards and pre-decrement.
        Circulator i = c;
        int su = 0;
        int k  = 5;
        do {
            Circulator j = --i;
            CGAL_assertion(  i ==  j);
            CGAL_assertion( (*i) == (*j));
            CGAL_assertion( k == (*i));
            su += (*i);
            --k;
        } while (i != c);
        CGAL_assertion( i == c);
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = c;
        su = 0;
        k  = 5;
        // Loop with post-decrement.
        do {
            Circulator j = i--;
            CGAL_assertion(  i !=  j);
            if ( j != c) {
                CGAL_assertion( (*i) == (*j) - 1);
            }
            CGAL_assertion( k == (*i));
            su += (*i);
            --k;
        } while (i != c);
        CGAL_assertion( i == c);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_circulator( c);
        CGAL::Assert_circulator( c);
    
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
        Circulator::size_type s = 5;
        CGAL_assertion( s == 5);
        Circulator::difference_type d = -5;
        CGAL_assertion( d == -5);
    
        // Check tests for empty data structures.
        Circulator z = Circulator();
        CGAL_assertion(   z == CGAL_CIRC_NULL);
        CGAL_assertion( ! (z != CGAL_CIRC_NULL));
        Circulator i = c;
        CGAL_assertion( ! (i == CGAL_CIRC_NULL));
        CGAL_assertion(   i != CGAL_CIRC_NULL);
        CGAL_assertion( i == c);
        CGAL_assertion( i == c);
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
        Circulator i = c;
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
    c++;
    CGAL_assertion( (*(c.current_iterator())) == 2);
  }{
    typedef Random_access_circulator_from_container< std::vector<int> >
        Circulator;
    Circulator c( &V);
    Assert_random_access_category(c);
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(c);
        CGAL::Assert_circulator_or_iterator(c);
        CGAL::Assert_is_at_least_forward_category(c);
        CGAL::Assert_is_at_least_forward_category(c);
        typedef std::iterator_traits<Circulator> I_Traits;
        typedef I_Traits::value_type      I_value_type;
        typedef I_Traits::difference_type I_difference_type;
        CGAL_assertion(1==test_value_type( (I_value_type*)(0)));
        CGAL_assertion(1==test_difference_type( (I_difference_type*)(0)));
    
        // Default constructor.
        Circulator z = Circulator();
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Circulator i = c;
    
        // Check general support for circulators and iterators.
        CGAL_assertion( CGAL::is_empty_range( z, z));
        CGAL_assertion( ! CGAL::is_empty_range( i, c));
    
        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, c)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                Circulator j = ++i;
                CGAL_assertion(  i ==  j);
                if ( i != c) {
                    CGAL_assertion( (*i) == (*j));
                }
            } while (i != c);  // Inequality and equality checked.
        }
        CGAL_assertion( i == c);  // Equality checked.
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = c;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, c)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                Circulator j = i++;
                CGAL_assertion(  i !=  j);
                if ( i != c) {
                    CGAL_assertion( (*i) == (*j) + 1);
                }
            } while (i != c);
        }
        CGAL_assertion( i == c);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        // Change three elements and check post-/pre-increment.
        Circulator i = c;
        (*i++) = 4;
        CGAL_assertion( 4 == (*c));
        CGAL_assertion( 2 == (*i));
        (*i++) = 3;
        CGAL_assertion( 3 == (*i));
        (*++i) = 7;
        CGAL_assertion( 7 == (*i));
    
        // Check the setting and reset these elements
        // to their original values.
        i = c;
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
        i = c;
        int k = 1;
        do {
            CGAL_assertion( k == (*i));
            ++i;
            ++k;
        } while (i != c);
    }
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(c);
        CGAL::Assert_circulator_or_iterator(c);
        CGAL::Assert_is_at_least_forward_category(c);
        CGAL::Assert_is_at_least_forward_category(c);
        typedef std::iterator_traits<Circulator> I_Traits;
        typedef I_Traits::value_type      I_value_type;
        typedef I_Traits::difference_type I_difference_type;
        CGAL_assertion(1==test_value_type( (I_value_type*)(0)));
        CGAL_assertion(1==test_difference_type( (I_difference_type*)(0)));
    
        // Default constructor.
        Circulator z = Circulator();
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Circulator i = c;
    
        // Check general support for circulators and iterators.
        CGAL_assertion( CGAL::is_empty_range( z, z));
        CGAL_assertion( ! CGAL::is_empty_range( i, c));
    
        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, c)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                Circulator j = ++i;
                CGAL_assertion(  i ==  j);
                if ( i != c) {
                    CGAL_assertion( (*i) == (*j));
                }
            } while (i != c);  // Inequality and equality checked.
        }
        CGAL_assertion( i == c);  // Equality checked.
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = c;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, c)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                Circulator j = i++;
                CGAL_assertion(  i !=  j);
                if ( i != c) {
                    CGAL_assertion( (*i) == (*j) + 1);
                }
            } while (i != c);
        }
        CGAL_assertion( i == c);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_bidirectional_category(c);
        CGAL::Assert_is_at_least_bidirectional_category(c);
        // Loop backwards and pre-decrement.
        Circulator i = c;
        int su = 0;
        int k  = 5;
        do {
            Circulator j = --i;
            CGAL_assertion(  i ==  j);
            CGAL_assertion( (*i) == (*j));
            CGAL_assertion( k == (*i));
            su += (*i);
            --k;
        } while (i != c);
        CGAL_assertion( i == c);
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = c;
        su = 0;
        k  = 5;
        // Loop with post-decrement.
        do {
            Circulator j = i--;
            CGAL_assertion(  i !=  j);
            if ( j != c) {
                CGAL_assertion( (*i) == (*j) - 1);
            }
            CGAL_assertion( k == (*i));
            su += (*i);
            --k;
        } while (i != c);
        CGAL_assertion( i == c);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(c);
        CGAL::Assert_circulator_or_iterator(c);
        CGAL::Assert_is_at_least_forward_category(c);
        CGAL::Assert_is_at_least_forward_category(c);
        typedef std::iterator_traits<Circulator> I_Traits;
        typedef I_Traits::value_type      I_value_type;
        typedef I_Traits::difference_type I_difference_type;
        CGAL_assertion(1==test_value_type( (I_value_type*)(0)));
        CGAL_assertion(1==test_difference_type( (I_difference_type*)(0)));
    
        // Default constructor.
        Circulator z = Circulator();
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Circulator i = c;
    
        // Check general support for circulators and iterators.
        CGAL_assertion( CGAL::is_empty_range( z, z));
        CGAL_assertion( ! CGAL::is_empty_range( i, c));
    
        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, c)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                Circulator j = ++i;
                CGAL_assertion(  i ==  j);
                if ( i != c) {
                    CGAL_assertion( (*i) == (*j));
                }
            } while (i != c);  // Inequality and equality checked.
        }
        CGAL_assertion( i == c);  // Equality checked.
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = c;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, c)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                Circulator j = i++;
                CGAL_assertion(  i !=  j);
                if ( i != c) {
                    CGAL_assertion( (*i) == (*j) + 1);
                }
            } while (i != c);
        }
        CGAL_assertion( i == c);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_bidirectional_category(c);
        CGAL::Assert_is_at_least_bidirectional_category(c);
        // Loop backwards and pre-decrement.
        Circulator i = c;
        int su = 0;
        int k  = 5;
        do {
            Circulator j = --i;
            CGAL_assertion(  i ==  j);
            CGAL_assertion( (*i) == (*j));
            CGAL_assertion( k == (*i));
            su += (*i);
            --k;
        } while (i != c);
        CGAL_assertion( i == c);
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = c;
        su = 0;
        k  = 5;
        // Loop with post-decrement.
        do {
            Circulator j = i--;
            CGAL_assertion(  i !=  j);
            if ( j != c) {
                CGAL_assertion( (*i) == (*j) - 1);
            }
            CGAL_assertion( k == (*i));
            su += (*i);
            --k;
        } while (i != c);
        CGAL_assertion( i == c);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_random_access_category(c);
        CGAL::Assert_is_at_least_random_access_category(c);
        // Random access.
        int k;
        for( k = 0; k < 5; k++) {
            CGAL_assertion( 1+k == c[k]);
        }
        int su = c[0]
               + c[1]
               + c[2]
               + c[3]
               + c[4];
        CGAL_assertion( su == 15);
    
        // Jump around.
        Circulator i = c;
        i += 3;
        CGAL_assertion( 4 == (*i));
        i -= 2;
        CGAL_assertion( 2 == (*i));
        i += 3;
        CGAL_assertion( 5 == (*i));
        i -= 4;
        CGAL_assertion( 1 == (*i));
        CGAL_assertion( i == c);
        Circulator j = i + 3;
        CGAL_assertion( 4 == (*j));
        Circulator jj = j - 2;
        CGAL_assertion( 2 == (*jj));
        jj = 4 + jj;
        CGAL_assertion( jj == c);
        Circulator ij = jj - 5;
        CGAL_assertion( ij == c);
    
        // Difference test.
        CGAL_assertion( jj - i == 5  ||  jj - i == 0);
        CGAL_assertion( i + (j-i) == j);
        CGAL_assertion( (j-i) + i == j);
    }
    { // Open own scope to hide local variables.
        Circulator i = c;
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
        i = c;
        int k = 1;
        do {
            CGAL_assertion( k == (*i));
            ++i;
            ++k;
        } while (i != c);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_circulator( c);
        CGAL::Assert_circulator( c);
    
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
        Circulator::size_type s = 5;
        CGAL_assertion( s == 5);
        Circulator::difference_type d = -5;
        CGAL_assertion( d == -5);
    
        // Check tests for empty data structures.
        Circulator z = Circulator();
        CGAL_assertion(   z == CGAL_CIRC_NULL);
        CGAL_assertion( ! (z != CGAL_CIRC_NULL));
        Circulator i = c;
        CGAL_assertion( ! (i == CGAL_CIRC_NULL));
        CGAL_assertion(   i != CGAL_CIRC_NULL);
        CGAL_assertion( i == c);
        CGAL_assertion( i == c);
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
        Circulator i = c;
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
        Circulator::difference_type d = c - c;
        CGAL_assertion( d == 0);
        d = c - c;
        CGAL_assertion( d == 0);
        Circulator i = c + 1;
        CGAL_assertion( c - i == 1 ||  c - i == -1);
        CGAL_assertion( i - c == 1 ||  i - c == -1);
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
    c++;
    CGAL_assertion( (*(c.current_iterator())) == 2);
  }{
    typedef
        Random_access_const_circulator_from_container<std::vector<int> >
        Circulator;
    Circulator c( &V);
    Assert_random_access_category(c);
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(c);
        CGAL::Assert_circulator_or_iterator(c);
        CGAL::Assert_is_at_least_forward_category(c);
        CGAL::Assert_is_at_least_forward_category(c);
        typedef std::iterator_traits<Circulator> I_Traits;
        typedef I_Traits::value_type      I_value_type;
        typedef I_Traits::difference_type I_difference_type;
        CGAL_assertion(1==test_value_type( (I_value_type*)(0)));
        CGAL_assertion(1==test_difference_type( (I_difference_type*)(0)));
    
        // Default constructor.
        Circulator z = Circulator();
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Circulator i = c;
    
        // Check general support for circulators and iterators.
        CGAL_assertion( CGAL::is_empty_range( z, z));
        CGAL_assertion( ! CGAL::is_empty_range( i, c));
    
        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, c)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                Circulator j = ++i;
                CGAL_assertion(  i ==  j);
                if ( i != c) {
                    CGAL_assertion( (*i) == (*j));
                }
            } while (i != c);  // Inequality and equality checked.
        }
        CGAL_assertion( i == c);  // Equality checked.
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = c;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, c)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                Circulator j = i++;
                CGAL_assertion(  i !=  j);
                if ( i != c) {
                    CGAL_assertion( (*i) == (*j) + 1);
                }
            } while (i != c);
        }
        CGAL_assertion( i == c);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_bidirectional_category(c);
        CGAL::Assert_is_at_least_bidirectional_category(c);
        // Loop backwards and pre-decrement.
        Circulator i = c;
        int su = 0;
        int k  = 5;
        do {
            Circulator j = --i;
            CGAL_assertion(  i ==  j);
            CGAL_assertion( (*i) == (*j));
            CGAL_assertion( k == (*i));
            su += (*i);
            --k;
        } while (i != c);
        CGAL_assertion( i == c);
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = c;
        su = 0;
        k  = 5;
        // Loop with post-decrement.
        do {
            Circulator j = i--;
            CGAL_assertion(  i !=  j);
            if ( j != c) {
                CGAL_assertion( (*i) == (*j) - 1);
            }
            CGAL_assertion( k == (*i));
            su += (*i);
            --k;
        } while (i != c);
        CGAL_assertion( i == c);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_random_access_category(c);
        CGAL::Assert_is_at_least_random_access_category(c);
        // Random access.
        int k;
        for( k = 0; k < 5; k++) {
            CGAL_assertion( 1+k == c[k]);
        }
        int su = c[0]
               + c[1]
               + c[2]
               + c[3]
               + c[4];
        CGAL_assertion( su == 15);
    
        // Jump around.
        Circulator i = c;
        i += 3;
        CGAL_assertion( 4 == (*i));
        i -= 2;
        CGAL_assertion( 2 == (*i));
        i += 3;
        CGAL_assertion( 5 == (*i));
        i -= 4;
        CGAL_assertion( 1 == (*i));
        CGAL_assertion( i == c);
        Circulator j = i + 3;
        CGAL_assertion( 4 == (*j));
        Circulator jj = j - 2;
        CGAL_assertion( 2 == (*jj));
        jj = 4 + jj;
        CGAL_assertion( jj == c);
        Circulator ij = jj - 5;
        CGAL_assertion( ij == c);
    
        // Difference test.
        CGAL_assertion( jj - i == 5  ||  jj - i == 0);
        CGAL_assertion( i + (j-i) == j);
        CGAL_assertion( (j-i) + i == j);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_circulator( c);
        CGAL::Assert_circulator( c);
    
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
        Circulator::size_type s = 5;
        CGAL_assertion( s == 5);
        Circulator::difference_type d = -5;
        CGAL_assertion( d == -5);
    
        // Check tests for empty data structures.
        Circulator z = Circulator();
        CGAL_assertion(   z == CGAL_CIRC_NULL);
        CGAL_assertion( ! (z != CGAL_CIRC_NULL));
        Circulator i = c;
        CGAL_assertion( ! (i == CGAL_CIRC_NULL));
        CGAL_assertion(   i != CGAL_CIRC_NULL);
        CGAL_assertion( i == c);
        CGAL_assertion( i == c);
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
        Circulator i = c;
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
        Circulator::difference_type d = c - c;
        CGAL_assertion( d == 0);
        d = c - c;
        CGAL_assertion( d == 0);
        Circulator i = c + 1;
        CGAL_assertion( c - i == 1 ||  c - i == -1);
        CGAL_assertion( i - c == 1 ||  i - c == -1);
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
    c++;
    CGAL_assertion( (*(c.current_iterator())) == 2);
  }{
    // the example program `circulator_prog2.C'.
    typedef Random_access_circulator_from_container< std::vector<int> >
        Circulator;
    typedef Random_access_container_from_circulator<Circulator>
        Container;
    typedef Container::iterator Iterator;

    std::vector<int> v;
    v.push_back(5);
    v.push_back(2);
    v.push_back(9);
    Circulator c( &v);
    Container  container( c);
    std::sort( container.begin(), container.end());
    Iterator i = container.begin();
    CGAL_assertion( *i == 2);
    i++;
    CGAL_assertion( *i == 5);
    i++;
    CGAL_assertion( *i == 9);
    i++;
    CGAL_assertion( i == container.end());
  }{
    // An example program applying sort() through two adaptors.
    typedef Random_access_circulator_from_container< std::vector<int> >
        Circulator;
    typedef Random_access_container_from_circulator<Circulator>
        Container;
    std::vector<int> v;
    v.push_back(5);
    v.push_back(2);
    v.push_back(9);
    Circulator c( &v);
    Container  container( c);
    std::sort( container.begin(), container.end());
    CGAL_assertion( v.begin()[0] == 2);
    CGAL_assertion( v.begin()[1] == 5);
    CGAL_assertion( v.begin()[2] == 9);
  }
}


int main(){
    init_global_data();
    test_tags();
    test_iterator_and_circulators();
    test_functions_for_circulators();

    test_container_from_circulator();
    test_circulator_from_iterator();
    test_circulator_from_container();
    clean_global_data();
    return 0;
}
// EOF //
