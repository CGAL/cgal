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
// file          : test_circ2.C
// chapter       : $CGAL_Chapter: Circulators $
// package       : $CGAL_Package: Circulator 3.4 (02 Sep 1999) $
// source        : circulator.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : INRIA, Sophia Antipolis
//
// Test support to build own circulators.
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
void test_struct(){
  Node* data_struct = generate_nodes( 5);
  {
    Struct_circulator start(data_struct);
    Assert_forward_category(start);
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(start);
        CGAL::Assert_circulator_or_iterator(start);
        CGAL::Assert_is_at_least_forward_category(start);
        CGAL::Assert_is_at_least_forward_category(start);
        typedef std::iterator_traits<Struct_circulator> I_Traits;
        typedef I_Traits::value_type      I_value_type;
        typedef I_Traits::difference_type I_difference_type;
        CGAL_assertion(1==test_value_type( (I_value_type*)(0)));
        CGAL_assertion(1==test_difference_type( (I_difference_type*)(0)));
    
        // Default constructor.
        Struct_circulator z = Struct_circulator();
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Struct_circulator i = start;
    
        // Check general support for circulators and iterators.
        CGAL_assertion( CGAL::is_empty_range( z, z));
        CGAL_assertion( ! CGAL::is_empty_range( i, start));
    
        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, start)) {   // superfluous
            do {
                CGAL_assertion( k == (*i).key);
                su += (*i).key;
                ++k;
                Struct_circulator j = ++i;
                CGAL_assertion(  i ==  j);
                if ( i != start) {
                    CGAL_assertion( (*i).key == (*j).key);
                }
            } while (i != start);  // Inequality and equality checked.
        }
        CGAL_assertion( i == start);  // Equality checked.
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = start;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, start)) {   // superfluous
            do {
                CGAL_assertion( k == (*i).key);
                su += (*i).key;
                ++k;
                Struct_circulator j = i++;
                CGAL_assertion(  i !=  j);
                if ( i != start) {
                    CGAL_assertion( (*i).key == (*j).key + 1);
                }
            } while (i != start);
        }
        CGAL_assertion( i == start);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        // Change three elements and check post-/pre-increment.
        Struct_circulator i = start;
        (*i++).key = 4;
        CGAL_assertion( 4 == (*start).key);
        CGAL_assertion( 2 == (*i).key);
        (*i++).key = 3;
        CGAL_assertion( 3 == (*i).key);
        (*++i).key = 7;
        CGAL_assertion( 7 == (*i).key);
    
        // Check the setting and reset these elements
        // to their original values.
        i = start;
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
        i = start;
        int k = 1;
        do {
            CGAL_assertion( k == (*i).key);
            ++i;
            ++k;
        } while (i != start);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_circulator( start);
        CGAL::Assert_circulator( start);
    
        // Check the local type parameters.
        Struct_circulator::value_type      k1;
        k1.key = 1;
        Struct_circulator::reference       k2 = k1;
        (void)k2;
        CGAL_assertion( k2.key == 1);
        Struct_circulator::pointer         p1 = &k1;
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
        Struct_circulator::size_type s = 5;
        CGAL_assertion( s == 5);
        Struct_circulator::difference_type d = -5;
        CGAL_assertion( d == -5);
    
        // Check tests for empty data structures.
        Struct_circulator z = Struct_circulator();
        CGAL_assertion(   z == CGAL_CIRC_NULL);
        CGAL_assertion( ! (z != CGAL_CIRC_NULL));
        Struct_circulator i = start;
        CGAL_assertion( ! (i == CGAL_CIRC_NULL));
        CGAL_assertion(   i != CGAL_CIRC_NULL);
        CGAL_assertion( i == start);
        CGAL_assertion( i == start);
        // Do I reach myself.
        ++i;
        Struct_circulator j = i;
        int k = 0;
        do {
            CGAL_assertion( k < 5);
            ++k;
            ++i;
        } while( i != j);
        CGAL_assertion( k == 5);
    }
  }{
    Struct_const_circulator start(data_struct);
    Assert_forward_category(start);
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(start);
        CGAL::Assert_circulator_or_iterator(start);
        CGAL::Assert_is_at_least_forward_category(start);
        CGAL::Assert_is_at_least_forward_category(start);
        typedef std::iterator_traits<Struct_const_circulator> I_Traits;
        typedef I_Traits::value_type      I_value_type;
        typedef I_Traits::difference_type I_difference_type;
        CGAL_assertion(1==test_value_type( (I_value_type*)(0)));
        CGAL_assertion(1==test_difference_type( (I_difference_type*)(0)));
    
        // Default constructor.
        Struct_const_circulator z = Struct_const_circulator();
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Struct_const_circulator i = start;
    
        // Check general support for circulators and iterators.
        CGAL_assertion( CGAL::is_empty_range( z, z));
        CGAL_assertion( ! CGAL::is_empty_range( i, start));
    
        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, start)) {   // superfluous
            do {
                CGAL_assertion( k == (*i).key);
                su += (*i).key;
                ++k;
                Struct_const_circulator j = ++i;
                CGAL_assertion(  i ==  j);
                if ( i != start) {
                    CGAL_assertion( (*i).key == (*j).key);
                }
            } while (i != start);  // Inequality and equality checked.
        }
        CGAL_assertion( i == start);  // Equality checked.
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = start;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, start)) {   // superfluous
            do {
                CGAL_assertion( k == (*i).key);
                su += (*i).key;
                ++k;
                Struct_const_circulator j = i++;
                CGAL_assertion(  i !=  j);
                if ( i != start) {
                    CGAL_assertion( (*i).key == (*j).key + 1);
                }
            } while (i != start);
        }
        CGAL_assertion( i == start);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_circulator( start);
        CGAL::Assert_circulator( start);
    
        // Check the local type parameters.
        Struct_const_circulator::value_type      k1;
        k1.key = 1;
        Struct_const_circulator::reference       k2 = k1;
        (void)k2;
        CGAL_assertion( k2.key == 1);
        Struct_const_circulator::pointer         p1 = &k1;
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
        Struct_const_circulator::size_type s = 5;
        CGAL_assertion( s == 5);
        Struct_const_circulator::difference_type d = -5;
        CGAL_assertion( d == -5);
    
        // Check tests for empty data structures.
        Struct_const_circulator z = Struct_const_circulator();
        CGAL_assertion(   z == CGAL_CIRC_NULL);
        CGAL_assertion( ! (z != CGAL_CIRC_NULL));
        Struct_const_circulator i = start;
        CGAL_assertion( ! (i == CGAL_CIRC_NULL));
        CGAL_assertion(   i != CGAL_CIRC_NULL);
        CGAL_assertion( i == start);
        CGAL_assertion( i == start);
        // Do I reach myself.
        ++i;
        Struct_const_circulator j = i;
        int k = 0;
        do {
            CGAL_assertion( k < 5);
            ++k;
            ++i;
        } while( i != j);
        CGAL_assertion( k == 5);
    }
  }{
    Struct_bi_circulator start(data_struct);
    Assert_bidirectional_category(start);
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(start);
        CGAL::Assert_circulator_or_iterator(start);
        CGAL::Assert_is_at_least_forward_category(start);
        CGAL::Assert_is_at_least_forward_category(start);
        typedef std::iterator_traits<Struct_bi_circulator> I_Traits;
        typedef I_Traits::value_type      I_value_type;
        typedef I_Traits::difference_type I_difference_type;
        CGAL_assertion(1==test_value_type( (I_value_type*)(0)));
        CGAL_assertion(1==test_difference_type( (I_difference_type*)(0)));
    
        // Default constructor.
        Struct_bi_circulator z = Struct_bi_circulator();
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Struct_bi_circulator i = start;
    
        // Check general support for circulators and iterators.
        CGAL_assertion( CGAL::is_empty_range( z, z));
        CGAL_assertion( ! CGAL::is_empty_range( i, start));
    
        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, start)) {   // superfluous
            do {
                CGAL_assertion( k == (*i).key);
                su += (*i).key;
                ++k;
                Struct_bi_circulator j = ++i;
                CGAL_assertion(  i ==  j);
                if ( i != start) {
                    CGAL_assertion( (*i).key == (*j).key);
                }
            } while (i != start);  // Inequality and equality checked.
        }
        CGAL_assertion( i == start);  // Equality checked.
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = start;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, start)) {   // superfluous
            do {
                CGAL_assertion( k == (*i).key);
                su += (*i).key;
                ++k;
                Struct_bi_circulator j = i++;
                CGAL_assertion(  i !=  j);
                if ( i != start) {
                    CGAL_assertion( (*i).key == (*j).key + 1);
                }
            } while (i != start);
        }
        CGAL_assertion( i == start);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        // Change three elements and check post-/pre-increment.
        Struct_bi_circulator i = start;
        (*i++).key = 4;
        CGAL_assertion( 4 == (*start).key);
        CGAL_assertion( 2 == (*i).key);
        (*i++).key = 3;
        CGAL_assertion( 3 == (*i).key);
        (*++i).key = 7;
        CGAL_assertion( 7 == (*i).key);
    
        // Check the setting and reset these elements
        // to their original values.
        i = start;
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
        i = start;
        int k = 1;
        do {
            CGAL_assertion( k == (*i).key);
            ++i;
            ++k;
        } while (i != start);
    }
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(start);
        CGAL::Assert_circulator_or_iterator(start);
        CGAL::Assert_is_at_least_forward_category(start);
        CGAL::Assert_is_at_least_forward_category(start);
        typedef std::iterator_traits<Struct_bi_circulator> I_Traits;
        typedef I_Traits::value_type      I_value_type;
        typedef I_Traits::difference_type I_difference_type;
        CGAL_assertion(1==test_value_type( (I_value_type*)(0)));
        CGAL_assertion(1==test_difference_type( (I_difference_type*)(0)));
    
        // Default constructor.
        Struct_bi_circulator z = Struct_bi_circulator();
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Struct_bi_circulator i = start;
    
        // Check general support for circulators and iterators.
        CGAL_assertion( CGAL::is_empty_range( z, z));
        CGAL_assertion( ! CGAL::is_empty_range( i, start));
    
        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, start)) {   // superfluous
            do {
                CGAL_assertion( k == (*i).key);
                su += (*i).key;
                ++k;
                Struct_bi_circulator j = ++i;
                CGAL_assertion(  i ==  j);
                if ( i != start) {
                    CGAL_assertion( (*i).key == (*j).key);
                }
            } while (i != start);  // Inequality and equality checked.
        }
        CGAL_assertion( i == start);  // Equality checked.
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = start;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, start)) {   // superfluous
            do {
                CGAL_assertion( k == (*i).key);
                su += (*i).key;
                ++k;
                Struct_bi_circulator j = i++;
                CGAL_assertion(  i !=  j);
                if ( i != start) {
                    CGAL_assertion( (*i).key == (*j).key + 1);
                }
            } while (i != start);
        }
        CGAL_assertion( i == start);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_bidirectional_category(start);
        CGAL::Assert_is_at_least_bidirectional_category(start);
        // Loop backwards and pre-decrement.
        Struct_bi_circulator i = start;
        int su = 0;
        int k  = 5;
        do {
            Struct_bi_circulator j = --i;
            CGAL_assertion(  i ==  j);
            CGAL_assertion( (*i).key == (*j).key);
            CGAL_assertion( k == (*i).key);
            su += (*i).key;
            --k;
        } while (i != start);
        CGAL_assertion( i == start);
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = start;
        su = 0;
        k  = 5;
        // Loop with post-decrement.
        do {
            Struct_bi_circulator j = i--;
            CGAL_assertion(  i !=  j);
            if ( j != start) {
                CGAL_assertion( (*i).key == (*j).key - 1);
            }
            CGAL_assertion( k == (*i).key);
            su += (*i).key;
            --k;
        } while (i != start);
        CGAL_assertion( i == start);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_circulator( start);
        CGAL::Assert_circulator( start);
    
        // Check the local type parameters.
        Struct_bi_circulator::value_type      k1;
        k1.key = 1;
        Struct_bi_circulator::reference       k2 = k1;
        (void)k2;
        CGAL_assertion( k2.key == 1);
        Struct_bi_circulator::pointer         p1 = &k1;
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
        Struct_bi_circulator::size_type s = 5;
        CGAL_assertion( s == 5);
        Struct_bi_circulator::difference_type d = -5;
        CGAL_assertion( d == -5);
    
        // Check tests for empty data structures.
        Struct_bi_circulator z = Struct_bi_circulator();
        CGAL_assertion(   z == CGAL_CIRC_NULL);
        CGAL_assertion( ! (z != CGAL_CIRC_NULL));
        Struct_bi_circulator i = start;
        CGAL_assertion( ! (i == CGAL_CIRC_NULL));
        CGAL_assertion(   i != CGAL_CIRC_NULL);
        CGAL_assertion( i == start);
        CGAL_assertion( i == start);
        // Do I reach myself.
        ++i;
        Struct_bi_circulator j = i;
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
        Struct_bi_circulator i = start;
        ++i;
        Struct_bi_circulator j = i;
        int k = 0;
        do {
            CGAL_assertion( k < 5);
            ++k;
            --i;
        } while( i != j);
        CGAL_assertion( k == 5);
    }
  }{
    Struct_bi_const_circulator start(data_struct);
    Assert_bidirectional_category(start);
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(start);
        CGAL::Assert_circulator_or_iterator(start);
        CGAL::Assert_is_at_least_forward_category(start);
        CGAL::Assert_is_at_least_forward_category(start);
        typedef std::iterator_traits<Struct_bi_const_circulator> I_Traits;
        typedef I_Traits::value_type      I_value_type;
        typedef I_Traits::difference_type I_difference_type;
        CGAL_assertion(1==test_value_type( (I_value_type*)(0)));
        CGAL_assertion(1==test_difference_type( (I_difference_type*)(0)));
    
        // Default constructor.
        Struct_bi_const_circulator z = Struct_bi_const_circulator();
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Struct_bi_const_circulator i = start;
    
        // Check general support for circulators and iterators.
        CGAL_assertion( CGAL::is_empty_range( z, z));
        CGAL_assertion( ! CGAL::is_empty_range( i, start));
    
        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, start)) {   // superfluous
            do {
                CGAL_assertion( k == (*i).key);
                su += (*i).key;
                ++k;
                Struct_bi_const_circulator j = ++i;
                CGAL_assertion(  i ==  j);
                if ( i != start) {
                    CGAL_assertion( (*i).key == (*j).key);
                }
            } while (i != start);  // Inequality and equality checked.
        }
        CGAL_assertion( i == start);  // Equality checked.
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = start;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, start)) {   // superfluous
            do {
                CGAL_assertion( k == (*i).key);
                su += (*i).key;
                ++k;
                Struct_bi_const_circulator j = i++;
                CGAL_assertion(  i !=  j);
                if ( i != start) {
                    CGAL_assertion( (*i).key == (*j).key + 1);
                }
            } while (i != start);
        }
        CGAL_assertion( i == start);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_bidirectional_category(start);
        CGAL::Assert_is_at_least_bidirectional_category(start);
        // Loop backwards and pre-decrement.
        Struct_bi_const_circulator i = start;
        int su = 0;
        int k  = 5;
        do {
            Struct_bi_const_circulator j = --i;
            CGAL_assertion(  i ==  j);
            CGAL_assertion( (*i).key == (*j).key);
            CGAL_assertion( k == (*i).key);
            su += (*i).key;
            --k;
        } while (i != start);
        CGAL_assertion( i == start);
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = start;
        su = 0;
        k  = 5;
        // Loop with post-decrement.
        do {
            Struct_bi_const_circulator j = i--;
            CGAL_assertion(  i !=  j);
            if ( j != start) {
                CGAL_assertion( (*i).key == (*j).key - 1);
            }
            CGAL_assertion( k == (*i).key);
            su += (*i).key;
            --k;
        } while (i != start);
        CGAL_assertion( i == start);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_circulator( start);
        CGAL::Assert_circulator( start);
    
        // Check the local type parameters.
        Struct_bi_const_circulator::value_type      k1;
        k1.key = 1;
        Struct_bi_const_circulator::reference       k2 = k1;
        (void)k2;
        CGAL_assertion( k2.key == 1);
        Struct_bi_const_circulator::pointer         p1 = &k1;
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
        Struct_bi_const_circulator::size_type s = 5;
        CGAL_assertion( s == 5);
        Struct_bi_const_circulator::difference_type d = -5;
        CGAL_assertion( d == -5);
    
        // Check tests for empty data structures.
        Struct_bi_const_circulator z = Struct_bi_const_circulator();
        CGAL_assertion(   z == CGAL_CIRC_NULL);
        CGAL_assertion( ! (z != CGAL_CIRC_NULL));
        Struct_bi_const_circulator i = start;
        CGAL_assertion( ! (i == CGAL_CIRC_NULL));
        CGAL_assertion(   i != CGAL_CIRC_NULL);
        CGAL_assertion( i == start);
        CGAL_assertion( i == start);
        // Do I reach myself.
        ++i;
        Struct_bi_const_circulator j = i;
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
        Struct_bi_const_circulator i = start;
        ++i;
        Struct_bi_const_circulator j = i;
        int k = 0;
        do {
            CGAL_assertion( k < 5);
            ++k;
            --i;
        } while( i != j);
        CGAL_assertion( k == 5);
    }
  }
  delete_nodes(data_struct);
}
void test_class(){
  CNode* data_struct = generate_cnodes( 5);
  {
    Class_circulator start(data_struct);
    Assert_forward_category(start);
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(start);
        CGAL::Assert_circulator_or_iterator(start);
        CGAL::Assert_is_at_least_forward_category(start);
        CGAL::Assert_is_at_least_forward_category(start);
        typedef std::iterator_traits<Class_circulator> I_Traits;
        typedef I_Traits::value_type      I_value_type;
        typedef I_Traits::difference_type I_difference_type;
        CGAL_assertion(1==test_value_type( (I_value_type*)(0)));
        CGAL_assertion(1==test_difference_type( (I_difference_type*)(0)));
    
        // Default constructor.
        Class_circulator z = Class_circulator();
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Class_circulator i = start;
    
        // Check general support for circulators and iterators.
        CGAL_assertion( CGAL::is_empty_range( z, z));
        CGAL_assertion( ! CGAL::is_empty_range( i, start));
    
        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, start)) {   // superfluous
            do {
                CGAL_assertion( k == (*i).key);
                su += (*i).key;
                ++k;
                Class_circulator j = ++i;
                CGAL_assertion(  i ==  j);
                if ( i != start) {
                    CGAL_assertion( (*i).key == (*j).key);
                }
            } while (i != start);  // Inequality and equality checked.
        }
        CGAL_assertion( i == start);  // Equality checked.
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = start;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, start)) {   // superfluous
            do {
                CGAL_assertion( k == (*i).key);
                su += (*i).key;
                ++k;
                Class_circulator j = i++;
                CGAL_assertion(  i !=  j);
                if ( i != start) {
                    CGAL_assertion( (*i).key == (*j).key + 1);
                }
            } while (i != start);
        }
        CGAL_assertion( i == start);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        // Change three elements and check post-/pre-increment.
        Class_circulator i = start;
        (*i++).key = 4;
        CGAL_assertion( 4 == (*start).key);
        CGAL_assertion( 2 == (*i).key);
        (*i++).key = 3;
        CGAL_assertion( 3 == (*i).key);
        (*++i).key = 7;
        CGAL_assertion( 7 == (*i).key);
    
        // Check the setting and reset these elements
        // to their original values.
        i = start;
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
        i = start;
        int k = 1;
        do {
            CGAL_assertion( k == (*i).key);
            ++i;
            ++k;
        } while (i != start);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_circulator( start);
        CGAL::Assert_circulator( start);
    
        // Check the local type parameters.
        Class_circulator::value_type      k1;
        k1.key = 1;
        Class_circulator::reference       k2 = k1;
        (void)k2;
        CGAL_assertion( k2.key == 1);
        Class_circulator::pointer         p1 = &k1;
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
        Class_circulator::size_type s = 5;
        CGAL_assertion( s == 5);
        Class_circulator::difference_type d = -5;
        CGAL_assertion( d == -5);
    
        // Check tests for empty data structures.
        Class_circulator z = Class_circulator();
        CGAL_assertion(   z == CGAL_CIRC_NULL);
        CGAL_assertion( ! (z != CGAL_CIRC_NULL));
        Class_circulator i = start;
        CGAL_assertion( ! (i == CGAL_CIRC_NULL));
        CGAL_assertion(   i != CGAL_CIRC_NULL);
        CGAL_assertion( i == start);
        CGAL_assertion( i == start);
        // Do I reach myself.
        ++i;
        Class_circulator j = i;
        int k = 0;
        do {
            CGAL_assertion( k < 5);
            ++k;
            ++i;
        } while( i != j);
        CGAL_assertion( k == 5);
    }
  }{
    Class_const_circulator start(data_struct);
    Assert_forward_category(start);
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(start);
        CGAL::Assert_circulator_or_iterator(start);
        CGAL::Assert_is_at_least_forward_category(start);
        CGAL::Assert_is_at_least_forward_category(start);
        typedef std::iterator_traits<Class_const_circulator> I_Traits;
        typedef I_Traits::value_type      I_value_type;
        typedef I_Traits::difference_type I_difference_type;
        CGAL_assertion(1==test_value_type( (I_value_type*)(0)));
        CGAL_assertion(1==test_difference_type( (I_difference_type*)(0)));
    
        // Default constructor.
        Class_const_circulator z = Class_const_circulator();
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Class_const_circulator i = start;
    
        // Check general support for circulators and iterators.
        CGAL_assertion( CGAL::is_empty_range( z, z));
        CGAL_assertion( ! CGAL::is_empty_range( i, start));
    
        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, start)) {   // superfluous
            do {
                CGAL_assertion( k == (*i).key);
                su += (*i).key;
                ++k;
                Class_const_circulator j = ++i;
                CGAL_assertion(  i ==  j);
                if ( i != start) {
                    CGAL_assertion( (*i).key == (*j).key);
                }
            } while (i != start);  // Inequality and equality checked.
        }
        CGAL_assertion( i == start);  // Equality checked.
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = start;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, start)) {   // superfluous
            do {
                CGAL_assertion( k == (*i).key);
                su += (*i).key;
                ++k;
                Class_const_circulator j = i++;
                CGAL_assertion(  i !=  j);
                if ( i != start) {
                    CGAL_assertion( (*i).key == (*j).key + 1);
                }
            } while (i != start);
        }
        CGAL_assertion( i == start);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_circulator( start);
        CGAL::Assert_circulator( start);
    
        // Check the local type parameters.
        Class_const_circulator::value_type      k1;
        k1.key = 1;
        Class_const_circulator::reference       k2 = k1;
        (void)k2;
        CGAL_assertion( k2.key == 1);
        Class_const_circulator::pointer         p1 = &k1;
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
        Class_const_circulator::size_type s = 5;
        CGAL_assertion( s == 5);
        Class_const_circulator::difference_type d = -5;
        CGAL_assertion( d == -5);
    
        // Check tests for empty data structures.
        Class_const_circulator z = Class_const_circulator();
        CGAL_assertion(   z == CGAL_CIRC_NULL);
        CGAL_assertion( ! (z != CGAL_CIRC_NULL));
        Class_const_circulator i = start;
        CGAL_assertion( ! (i == CGAL_CIRC_NULL));
        CGAL_assertion(   i != CGAL_CIRC_NULL);
        CGAL_assertion( i == start);
        CGAL_assertion( i == start);
        // Do I reach myself.
        ++i;
        Class_const_circulator j = i;
        int k = 0;
        do {
            CGAL_assertion( k < 5);
            ++k;
            ++i;
        } while( i != j);
        CGAL_assertion( k == 5);
    }
  }{
    Class_bi_circulator start(data_struct);
    Assert_bidirectional_category(start);
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(start);
        CGAL::Assert_circulator_or_iterator(start);
        CGAL::Assert_is_at_least_forward_category(start);
        CGAL::Assert_is_at_least_forward_category(start);
        typedef std::iterator_traits<Class_bi_circulator> I_Traits;
        typedef I_Traits::value_type      I_value_type;
        typedef I_Traits::difference_type I_difference_type;
        CGAL_assertion(1==test_value_type( (I_value_type*)(0)));
        CGAL_assertion(1==test_difference_type( (I_difference_type*)(0)));
    
        // Default constructor.
        Class_bi_circulator z = Class_bi_circulator();
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Class_bi_circulator i = start;
    
        // Check general support for circulators and iterators.
        CGAL_assertion( CGAL::is_empty_range( z, z));
        CGAL_assertion( ! CGAL::is_empty_range( i, start));
    
        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, start)) {   // superfluous
            do {
                CGAL_assertion( k == (*i).key);
                su += (*i).key;
                ++k;
                Class_bi_circulator j = ++i;
                CGAL_assertion(  i ==  j);
                if ( i != start) {
                    CGAL_assertion( (*i).key == (*j).key);
                }
            } while (i != start);  // Inequality and equality checked.
        }
        CGAL_assertion( i == start);  // Equality checked.
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = start;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, start)) {   // superfluous
            do {
                CGAL_assertion( k == (*i).key);
                su += (*i).key;
                ++k;
                Class_bi_circulator j = i++;
                CGAL_assertion(  i !=  j);
                if ( i != start) {
                    CGAL_assertion( (*i).key == (*j).key + 1);
                }
            } while (i != start);
        }
        CGAL_assertion( i == start);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        // Change three elements and check post-/pre-increment.
        Class_bi_circulator i = start;
        (*i++).key = 4;
        CGAL_assertion( 4 == (*start).key);
        CGAL_assertion( 2 == (*i).key);
        (*i++).key = 3;
        CGAL_assertion( 3 == (*i).key);
        (*++i).key = 7;
        CGAL_assertion( 7 == (*i).key);
    
        // Check the setting and reset these elements
        // to their original values.
        i = start;
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
        i = start;
        int k = 1;
        do {
            CGAL_assertion( k == (*i).key);
            ++i;
            ++k;
        } while (i != start);
    }
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(start);
        CGAL::Assert_circulator_or_iterator(start);
        CGAL::Assert_is_at_least_forward_category(start);
        CGAL::Assert_is_at_least_forward_category(start);
        typedef std::iterator_traits<Class_bi_circulator> I_Traits;
        typedef I_Traits::value_type      I_value_type;
        typedef I_Traits::difference_type I_difference_type;
        CGAL_assertion(1==test_value_type( (I_value_type*)(0)));
        CGAL_assertion(1==test_difference_type( (I_difference_type*)(0)));
    
        // Default constructor.
        Class_bi_circulator z = Class_bi_circulator();
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Class_bi_circulator i = start;
    
        // Check general support for circulators and iterators.
        CGAL_assertion( CGAL::is_empty_range( z, z));
        CGAL_assertion( ! CGAL::is_empty_range( i, start));
    
        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, start)) {   // superfluous
            do {
                CGAL_assertion( k == (*i).key);
                su += (*i).key;
                ++k;
                Class_bi_circulator j = ++i;
                CGAL_assertion(  i ==  j);
                if ( i != start) {
                    CGAL_assertion( (*i).key == (*j).key);
                }
            } while (i != start);  // Inequality and equality checked.
        }
        CGAL_assertion( i == start);  // Equality checked.
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = start;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, start)) {   // superfluous
            do {
                CGAL_assertion( k == (*i).key);
                su += (*i).key;
                ++k;
                Class_bi_circulator j = i++;
                CGAL_assertion(  i !=  j);
                if ( i != start) {
                    CGAL_assertion( (*i).key == (*j).key + 1);
                }
            } while (i != start);
        }
        CGAL_assertion( i == start);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_bidirectional_category(start);
        CGAL::Assert_is_at_least_bidirectional_category(start);
        // Loop backwards and pre-decrement.
        Class_bi_circulator i = start;
        int su = 0;
        int k  = 5;
        do {
            Class_bi_circulator j = --i;
            CGAL_assertion(  i ==  j);
            CGAL_assertion( (*i).key == (*j).key);
            CGAL_assertion( k == (*i).key);
            su += (*i).key;
            --k;
        } while (i != start);
        CGAL_assertion( i == start);
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = start;
        su = 0;
        k  = 5;
        // Loop with post-decrement.
        do {
            Class_bi_circulator j = i--;
            CGAL_assertion(  i !=  j);
            if ( j != start) {
                CGAL_assertion( (*i).key == (*j).key - 1);
            }
            CGAL_assertion( k == (*i).key);
            su += (*i).key;
            --k;
        } while (i != start);
        CGAL_assertion( i == start);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_circulator( start);
        CGAL::Assert_circulator( start);
    
        // Check the local type parameters.
        Class_bi_circulator::value_type      k1;
        k1.key = 1;
        Class_bi_circulator::reference       k2 = k1;
        (void)k2;
        CGAL_assertion( k2.key == 1);
        Class_bi_circulator::pointer         p1 = &k1;
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
        Class_bi_circulator::size_type s = 5;
        CGAL_assertion( s == 5);
        Class_bi_circulator::difference_type d = -5;
        CGAL_assertion( d == -5);
    
        // Check tests for empty data structures.
        Class_bi_circulator z = Class_bi_circulator();
        CGAL_assertion(   z == CGAL_CIRC_NULL);
        CGAL_assertion( ! (z != CGAL_CIRC_NULL));
        Class_bi_circulator i = start;
        CGAL_assertion( ! (i == CGAL_CIRC_NULL));
        CGAL_assertion(   i != CGAL_CIRC_NULL);
        CGAL_assertion( i == start);
        CGAL_assertion( i == start);
        // Do I reach myself.
        ++i;
        Class_bi_circulator j = i;
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
        Class_bi_circulator i = start;
        ++i;
        Class_bi_circulator j = i;
        int k = 0;
        do {
            CGAL_assertion( k < 5);
            ++k;
            --i;
        } while( i != j);
        CGAL_assertion( k == 5);
    }
  }{
    Class_bi_const_circulator start(data_struct);
    Assert_bidirectional_category(start);
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(start);
        CGAL::Assert_circulator_or_iterator(start);
        CGAL::Assert_is_at_least_forward_category(start);
        CGAL::Assert_is_at_least_forward_category(start);
        typedef std::iterator_traits<Class_bi_const_circulator> I_Traits;
        typedef I_Traits::value_type      I_value_type;
        typedef I_Traits::difference_type I_difference_type;
        CGAL_assertion(1==test_value_type( (I_value_type*)(0)));
        CGAL_assertion(1==test_difference_type( (I_difference_type*)(0)));
    
        // Default constructor.
        Class_bi_const_circulator z = Class_bi_const_circulator();
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Class_bi_const_circulator i = start;
    
        // Check general support for circulators and iterators.
        CGAL_assertion( CGAL::is_empty_range( z, z));
        CGAL_assertion( ! CGAL::is_empty_range( i, start));
    
        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, start)) {   // superfluous
            do {
                CGAL_assertion( k == (*i).key);
                su += (*i).key;
                ++k;
                Class_bi_const_circulator j = ++i;
                CGAL_assertion(  i ==  j);
                if ( i != start) {
                    CGAL_assertion( (*i).key == (*j).key);
                }
            } while (i != start);  // Inequality and equality checked.
        }
        CGAL_assertion( i == start);  // Equality checked.
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = start;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, start)) {   // superfluous
            do {
                CGAL_assertion( k == (*i).key);
                su += (*i).key;
                ++k;
                Class_bi_const_circulator j = i++;
                CGAL_assertion(  i !=  j);
                if ( i != start) {
                    CGAL_assertion( (*i).key == (*j).key + 1);
                }
            } while (i != start);
        }
        CGAL_assertion( i == start);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_bidirectional_category(start);
        CGAL::Assert_is_at_least_bidirectional_category(start);
        // Loop backwards and pre-decrement.
        Class_bi_const_circulator i = start;
        int su = 0;
        int k  = 5;
        do {
            Class_bi_const_circulator j = --i;
            CGAL_assertion(  i ==  j);
            CGAL_assertion( (*i).key == (*j).key);
            CGAL_assertion( k == (*i).key);
            su += (*i).key;
            --k;
        } while (i != start);
        CGAL_assertion( i == start);
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = start;
        su = 0;
        k  = 5;
        // Loop with post-decrement.
        do {
            Class_bi_const_circulator j = i--;
            CGAL_assertion(  i !=  j);
            if ( j != start) {
                CGAL_assertion( (*i).key == (*j).key - 1);
            }
            CGAL_assertion( k == (*i).key);
            su += (*i).key;
            --k;
        } while (i != start);
        CGAL_assertion( i == start);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_circulator( start);
        CGAL::Assert_circulator( start);
    
        // Check the local type parameters.
        Class_bi_const_circulator::value_type      k1;
        k1.key = 1;
        Class_bi_const_circulator::reference       k2 = k1;
        (void)k2;
        CGAL_assertion( k2.key == 1);
        Class_bi_const_circulator::pointer         p1 = &k1;
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
        Class_bi_const_circulator::size_type s = 5;
        CGAL_assertion( s == 5);
        Class_bi_const_circulator::difference_type d = -5;
        CGAL_assertion( d == -5);
    
        // Check tests for empty data structures.
        Class_bi_const_circulator z = Class_bi_const_circulator();
        CGAL_assertion(   z == CGAL_CIRC_NULL);
        CGAL_assertion( ! (z != CGAL_CIRC_NULL));
        Class_bi_const_circulator i = start;
        CGAL_assertion( ! (i == CGAL_CIRC_NULL));
        CGAL_assertion(   i != CGAL_CIRC_NULL);
        CGAL_assertion( i == start);
        CGAL_assertion( i == start);
        // Do I reach myself.
        ++i;
        Class_bi_const_circulator j = i;
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
        Class_bi_const_circulator i = start;
        ++i;
        Class_bi_const_circulator j = i;
        int k = 0;
        do {
            CGAL_assertion( k < 5);
            ++k;
            --i;
        } while( i != j);
        CGAL_assertion( k == 5);
    }
  }
  delete_cnodes(data_struct);
}
void test_array() {
  {
    typedef Circulator_over_array<
        std::vector<int>,
        std::vector<int>::value_type,
        std::vector<int>::size_type,
        std::vector<int>::difference_type
    > Circulator;
    Circulator start( V, V.size());
    Assert_random_access_category(start);
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(start);
        CGAL::Assert_circulator_or_iterator(start);
        CGAL::Assert_is_at_least_forward_category(start);
        CGAL::Assert_is_at_least_forward_category(start);
        typedef std::iterator_traits<Circulator> I_Traits;
        typedef I_Traits::value_type      I_value_type;
        typedef I_Traits::difference_type I_difference_type;
        CGAL_assertion(1==test_value_type( (I_value_type*)(0)));
        CGAL_assertion(1==test_difference_type( (I_difference_type*)(0)));
    
        // Default constructor.
        Circulator z = Circulator();
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Circulator i = start;
    
        // Check general support for circulators and iterators.
        CGAL_assertion( CGAL::is_empty_range( z, z));
        CGAL_assertion( ! CGAL::is_empty_range( i, start));
    
        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, start)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                Circulator j = ++i;
                CGAL_assertion(  i ==  j);
                if ( i != start) {
                    CGAL_assertion( (*i) == (*j));
                }
            } while (i != start);  // Inequality and equality checked.
        }
        CGAL_assertion( i == start);  // Equality checked.
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = start;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, start)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                Circulator j = i++;
                CGAL_assertion(  i !=  j);
                if ( i != start) {
                    CGAL_assertion( (*i) == (*j) + 1);
                }
            } while (i != start);
        }
        CGAL_assertion( i == start);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        // Change three elements and check post-/pre-increment.
        Circulator i = start;
        (*i++) = 4;
        CGAL_assertion( 4 == (*start));
        CGAL_assertion( 2 == (*i));
        (*i++) = 3;
        CGAL_assertion( 3 == (*i));
        (*++i) = 7;
        CGAL_assertion( 7 == (*i));
    
        // Check the setting and reset these elements
        // to their original values.
        i = start;
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
        i = start;
        int k = 1;
        do {
            CGAL_assertion( k == (*i));
            ++i;
            ++k;
        } while (i != start);
    }
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(start);
        CGAL::Assert_circulator_or_iterator(start);
        CGAL::Assert_is_at_least_forward_category(start);
        CGAL::Assert_is_at_least_forward_category(start);
        typedef std::iterator_traits<Circulator> I_Traits;
        typedef I_Traits::value_type      I_value_type;
        typedef I_Traits::difference_type I_difference_type;
        CGAL_assertion(1==test_value_type( (I_value_type*)(0)));
        CGAL_assertion(1==test_difference_type( (I_difference_type*)(0)));
    
        // Default constructor.
        Circulator z = Circulator();
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Circulator i = start;
    
        // Check general support for circulators and iterators.
        CGAL_assertion( CGAL::is_empty_range( z, z));
        CGAL_assertion( ! CGAL::is_empty_range( i, start));
    
        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, start)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                Circulator j = ++i;
                CGAL_assertion(  i ==  j);
                if ( i != start) {
                    CGAL_assertion( (*i) == (*j));
                }
            } while (i != start);  // Inequality and equality checked.
        }
        CGAL_assertion( i == start);  // Equality checked.
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = start;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, start)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                Circulator j = i++;
                CGAL_assertion(  i !=  j);
                if ( i != start) {
                    CGAL_assertion( (*i) == (*j) + 1);
                }
            } while (i != start);
        }
        CGAL_assertion( i == start);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_bidirectional_category(start);
        CGAL::Assert_is_at_least_bidirectional_category(start);
        // Loop backwards and pre-decrement.
        Circulator i = start;
        int su = 0;
        int k  = 5;
        do {
            Circulator j = --i;
            CGAL_assertion(  i ==  j);
            CGAL_assertion( (*i) == (*j));
            CGAL_assertion( k == (*i));
            su += (*i);
            --k;
        } while (i != start);
        CGAL_assertion( i == start);
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = start;
        su = 0;
        k  = 5;
        // Loop with post-decrement.
        do {
            Circulator j = i--;
            CGAL_assertion(  i !=  j);
            if ( j != start) {
                CGAL_assertion( (*i) == (*j) - 1);
            }
            CGAL_assertion( k == (*i));
            su += (*i);
            --k;
        } while (i != start);
        CGAL_assertion( i == start);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(start);
        CGAL::Assert_circulator_or_iterator(start);
        CGAL::Assert_is_at_least_forward_category(start);
        CGAL::Assert_is_at_least_forward_category(start);
        typedef std::iterator_traits<Circulator> I_Traits;
        typedef I_Traits::value_type      I_value_type;
        typedef I_Traits::difference_type I_difference_type;
        CGAL_assertion(1==test_value_type( (I_value_type*)(0)));
        CGAL_assertion(1==test_difference_type( (I_difference_type*)(0)));
    
        // Default constructor.
        Circulator z = Circulator();
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Circulator i = start;
    
        // Check general support for circulators and iterators.
        CGAL_assertion( CGAL::is_empty_range( z, z));
        CGAL_assertion( ! CGAL::is_empty_range( i, start));
    
        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, start)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                Circulator j = ++i;
                CGAL_assertion(  i ==  j);
                if ( i != start) {
                    CGAL_assertion( (*i) == (*j));
                }
            } while (i != start);  // Inequality and equality checked.
        }
        CGAL_assertion( i == start);  // Equality checked.
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = start;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, start)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                Circulator j = i++;
                CGAL_assertion(  i !=  j);
                if ( i != start) {
                    CGAL_assertion( (*i) == (*j) + 1);
                }
            } while (i != start);
        }
        CGAL_assertion( i == start);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_bidirectional_category(start);
        CGAL::Assert_is_at_least_bidirectional_category(start);
        // Loop backwards and pre-decrement.
        Circulator i = start;
        int su = 0;
        int k  = 5;
        do {
            Circulator j = --i;
            CGAL_assertion(  i ==  j);
            CGAL_assertion( (*i) == (*j));
            CGAL_assertion( k == (*i));
            su += (*i);
            --k;
        } while (i != start);
        CGAL_assertion( i == start);
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = start;
        su = 0;
        k  = 5;
        // Loop with post-decrement.
        do {
            Circulator j = i--;
            CGAL_assertion(  i !=  j);
            if ( j != start) {
                CGAL_assertion( (*i) == (*j) - 1);
            }
            CGAL_assertion( k == (*i));
            su += (*i);
            --k;
        } while (i != start);
        CGAL_assertion( i == start);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_random_access_category(start);
        CGAL::Assert_is_at_least_random_access_category(start);
        // Random access.
        int k;
        for( k = 0; k < 5; k++) {
            CGAL_assertion( 1+k == start[k]);
        }
        int su = start[0]
               + start[1]
               + start[2]
               + start[3]
               + start[4];
        CGAL_assertion( su == 15);
    
        // Jump around.
        Circulator i = start;
        i += 3;
        CGAL_assertion( 4 == (*i));
        i -= 2;
        CGAL_assertion( 2 == (*i));
        i += 3;
        CGAL_assertion( 5 == (*i));
        i -= 4;
        CGAL_assertion( 1 == (*i));
        CGAL_assertion( i == start);
        Circulator j = i + 3;
        CGAL_assertion( 4 == (*j));
        Circulator jj = j - 2;
        CGAL_assertion( 2 == (*jj));
        jj = 4 + jj;
        CGAL_assertion( jj == start);
        Circulator ij = jj - 5;
        CGAL_assertion( ij == start);
    
        // Difference test.
        CGAL_assertion( jj - i == 5  ||  jj - i == 0);
        CGAL_assertion( i + (j-i) == j);
        CGAL_assertion( (j-i) + i == j);
    }
    { // Open own scope to hide local variables.
        Circulator i = start;
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
        i = start;
        int k = 1;
        do {
            CGAL_assertion( k == (*i));
            ++i;
            ++k;
        } while (i != start);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_circulator( start);
        CGAL::Assert_circulator( start);
    
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
        Circulator i = start;
        CGAL_assertion( ! (i == CGAL_CIRC_NULL));
        CGAL_assertion(   i != CGAL_CIRC_NULL);
        CGAL_assertion( i == start);
        CGAL_assertion( i == start);
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
        Circulator i = start;
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
        Circulator::difference_type d = start - start;
        CGAL_assertion( d == 0);
        d = start - start;
        CGAL_assertion( d == 0);
        Circulator i = start + 1;
        CGAL_assertion( start - i == 1 ||  start - i == -1);
        CGAL_assertion( i - start == 1 ||  i - start == -1);
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
  }{
    typedef Const_circulator_over_array<
        std::vector<int>,
        std::vector<int>::value_type,
        std::vector<int>::size_type,
        std::vector<int>::difference_type
    > Circulator;
    const std::vector<int>& W = V;
    Circulator start( W, W.size());
    Assert_random_access_category(start);
    { // Open own scope to hide local variables.
        // Check generally correct parameter properties.
        CGAL::Assert_circulator_or_iterator(start);
        CGAL::Assert_circulator_or_iterator(start);
        CGAL::Assert_is_at_least_forward_category(start);
        CGAL::Assert_is_at_least_forward_category(start);
        typedef std::iterator_traits<Circulator> I_Traits;
        typedef I_Traits::value_type      I_value_type;
        typedef I_Traits::difference_type I_difference_type;
        CGAL_assertion(1==test_value_type( (I_value_type*)(0)));
        CGAL_assertion(1==test_difference_type( (I_difference_type*)(0)));
    
        // Default constructor.
        Circulator z = Circulator();
        CGAL::Assert_circulator_or_iterator(z);
        // Copy constructor.
        Circulator i = start;
    
        // Check general support for circulators and iterators.
        CGAL_assertion( CGAL::is_empty_range( z, z));
        CGAL_assertion( ! CGAL::is_empty_range( i, start));
    
        int su = 0;
        int k  = 1;
        // Check general loop, pre-increment, dereference.
        if (! CGAL::is_empty_range( i, start)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                Circulator j = ++i;
                CGAL_assertion(  i ==  j);
                if ( i != start) {
                    CGAL_assertion( (*i) == (*j));
                }
            } while (i != start);  // Inequality and equality checked.
        }
        CGAL_assertion( i == start);  // Equality checked.
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = start;
        su = 0;
        k  = 1;
        // Loop with post increment.
        if (! CGAL::is_empty_range( i, start)) {   // superfluous
            do {
                CGAL_assertion( k == (*i));
                su += (*i);
                ++k;
                Circulator j = i++;
                CGAL_assertion(  i !=  j);
                if ( i != start) {
                    CGAL_assertion( (*i) == (*j) + 1);
                }
            } while (i != start);
        }
        CGAL_assertion( i == start);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_bidirectional_category(start);
        CGAL::Assert_is_at_least_bidirectional_category(start);
        // Loop backwards and pre-decrement.
        Circulator i = start;
        int su = 0;
        int k  = 5;
        do {
            Circulator j = --i;
            CGAL_assertion(  i ==  j);
            CGAL_assertion( (*i) == (*j));
            CGAL_assertion( k == (*i));
            su += (*i);
            --k;
        } while (i != start);
        CGAL_assertion( i == start);
        CGAL_assertion( su == 15);
    
        // Assignment.
        i = start;
        su = 0;
        k  = 5;
        // Loop with post-decrement.
        do {
            Circulator j = i--;
            CGAL_assertion(  i !=  j);
            if ( j != start) {
                CGAL_assertion( (*i) == (*j) - 1);
            }
            CGAL_assertion( k == (*i));
            su += (*i);
            --k;
        } while (i != start);
        CGAL_assertion( i == start);
        CGAL_assertion( su == 15);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_is_at_least_random_access_category(start);
        CGAL::Assert_is_at_least_random_access_category(start);
        // Random access.
        int k;
        for( k = 0; k < 5; k++) {
            CGAL_assertion( 1+k == start[k]);
        }
        int su = start[0]
               + start[1]
               + start[2]
               + start[3]
               + start[4];
        CGAL_assertion( su == 15);
    
        // Jump around.
        Circulator i = start;
        i += 3;
        CGAL_assertion( 4 == (*i));
        i -= 2;
        CGAL_assertion( 2 == (*i));
        i += 3;
        CGAL_assertion( 5 == (*i));
        i -= 4;
        CGAL_assertion( 1 == (*i));
        CGAL_assertion( i == start);
        Circulator j = i + 3;
        CGAL_assertion( 4 == (*j));
        Circulator jj = j - 2;
        CGAL_assertion( 2 == (*jj));
        jj = 4 + jj;
        CGAL_assertion( jj == start);
        Circulator ij = jj - 5;
        CGAL_assertion( ij == start);
    
        // Difference test.
        CGAL_assertion( jj - i == 5  ||  jj - i == 0);
        CGAL_assertion( i + (j-i) == j);
        CGAL_assertion( (j-i) + i == j);
    }
    { // Open own scope to hide local variables.
        CGAL::Assert_circulator( start);
        CGAL::Assert_circulator( start);
    
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
        Circulator i = start;
        CGAL_assertion( ! (i == CGAL_CIRC_NULL));
        CGAL_assertion(   i != CGAL_CIRC_NULL);
        CGAL_assertion( i == start);
        CGAL_assertion( i == start);
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
        Circulator i = start;
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
        Circulator::difference_type d = start - start;
        CGAL_assertion( d == 0);
        d = start - start;
        CGAL_assertion( d == 0);
        Circulator i = start + 1;
        CGAL_assertion( start - i == 1 ||  start - i == -1);
        CGAL_assertion( i - start == 1 ||  i - start == -1);
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
  }
}


int main(){
    init_global_data();
    test_struct();
    test_class();
    test_array();
    clean_global_data();
    return 0;
}
// EOF //
