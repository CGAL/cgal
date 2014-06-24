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
// file          : src/test_vector.C
// package       : $CGAL_Package: STL_Extension $
// chapter       : $CGAL_Chapter: STL Extensions for CGAL $
//
// revision      : $Id$
// revision_date : $Date$
//
// author(s)     : Andreas Fabri <Andreas.Fabri@sophia.inria.fr>
//                 Lutz Kettner <kettner@mpi-sb.mpg.de>
// maintainer    : Andreas Fabri <Andreas.Fabri@sophia.inria.fr>
// coordinator   : ETH Zurich
//
// A vector container class: internal replacement that works with declared
// but not yet defined types as they occur with cyclic type dependencies
// with templates.
// ============================================================================

#include <CGAL/vector.h>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <list>
#include <CGAL/use.h>

class X {
public:
    X(): x_(-1) {}
    X(int i) : x_(i) {}
    int x() const { return x_; }
private:
    int x_;
};

bool operator==( const X& x1, const X& x2) { return x1.x() == x2.x(); }
bool operator!=( const X& x1, const X& x2) { return x1.x() != x2.x(); }
bool operator< ( const X& x1, const X& x2) { return x1.x() <  x2.x(); }
bool operator<=( const X& x1, const X& x2) { return x1.x() <=  x2.x(); }
bool operator> ( const X& x1, const X& x2) { return x1.x() >  x2.x(); }
bool operator>=( const X& x1, const X& x2) { return x1.x() >= x2.x(); }

std::ostream& operator<< ( std::ostream& s, const X& x) {
    s << x.x();
    return s; 
}

int main() {
    typedef std::list<X>                      list;
    typedef CGAL::internal::vector<X>            vector;
    typedef vector::iterator                  iterator;
    typedef vector::const_iterator            const_iterator;
    typedef vector::reverse_iterator          reverse_iterator;
    typedef vector::const_reverse_iterator    const_reverse_iterator;
    typedef vector::value_type                value_type;
    typedef vector::pointer                   pointer;
    typedef vector::reference                 reference;
    typedef vector::allocator_type            allocator_type;
    typedef vector::difference_type           difference_type;
    typedef vector::size_type                 size_type;

    typedef iterator::value_type              I_value_type;
    typedef iterator::pointer                 I_pointer;
    typedef iterator::reference               I_reference;
    typedef iterator::difference_type         I_difference_type;
    typedef iterator::iterator_category       I_iterator_category;

    typedef const_iterator::value_type        C_value_type;
    typedef const_iterator::pointer           C_pointer;
    typedef const_iterator::reference         C_reference;
    typedef const_iterator::difference_type   C_difference_type;
    typedef const_iterator::iterator_category C_iterator_category;

    CGAL_USE_TYPE(const_reverse_iterator);
    CGAL_USE_TYPE(value_type);
    CGAL_USE_TYPE(pointer);
    CGAL_USE_TYPE(reference);
    CGAL_USE_TYPE(allocator_type);
    CGAL_USE_TYPE(difference_type);
    CGAL_USE_TYPE(size_type);
    CGAL_USE_TYPE(I_value_type);
    CGAL_USE_TYPE(I_pointer);
    CGAL_USE_TYPE(I_reference);
    CGAL_USE_TYPE(I_difference_type);
    CGAL_USE_TYPE(I_iterator_category);
    CGAL_USE_TYPE(C_value_type);
    CGAL_USE_TYPE(C_pointer);
    CGAL_USE_TYPE(C_reference);
    CGAL_USE_TYPE(C_difference_type);
    CGAL_USE_TYPE(C_iterator_category);


    vector V;
    std::cout << V.capacity() << ", " << V.size()<< std::endl;
    V.reserve(10);
    assert(V.empty());

    V.push_back(X(0));

    assert(! V.empty());

    V.push_back(X(1));
    V.push_back(X(2));

    X x0(0);
    X x2(2);
    assert( V.front() == x0);
    assert( V.back() == x2);
    iterator it = V.begin();
    assert( *it == X(0));
    it = V.end();
    it--;
    assert(*it == X(2));

    const_iterator cit = V.begin();
    assert( *cit == X(0));
    cit = V.end();
    cit--;
    assert(*cit == X(2));

    reverse_iterator rit = V.rbegin();
    assert( *rit == X(2));
    rit = V.rend();
    rit--;
    assert(*rit == X(0));

    std::copy( V.begin(), V.end(), std::ostream_iterator<X>( std::cout, ", "));

    vector V2(V);
    vector V3 = V;
    
  
    assert( V3[0] == X(0));
    assert( V3[1] == X(1));
    assert( V3[2] == X(2));

    assert( V == V2);
    assert( V == V3);

    V3[0] = X();
    assert( ! (V == V3));

    V.pop_back();
    V.pop_back();
    V.pop_back();
    assert( V.empty());
  
    V2.clear();
    assert( V2.empty());

    V3.insert( V3.begin(), X(5));
    assert( V3.size() == 4);
    V3.erase( V3.begin());
    assert( V3.size() == 3);
    V3.insert( V3.begin());
    assert( V3[0] == X());
    V3.erase( V3.begin());
    V3.erase( V3.begin());
    assert( V3[0] == X(1));
    V3.erase( V3.begin());
    assert( V3[0] == X(2));
    V3.erase( V3.begin());
    assert( V3.empty());

    vector V5(3);
    assert( V5.size() == 3);
    unsigned int i = 0;
    for( ; i < V5.size(); i++) {
        assert( V5[i] == X());
    }
    vector V6(3, X(7812));
    assert(V6.size() == 3);
    for( i = 0; i < V6.size(); i++) {
        assert( V6[i] == X(7812));
    }
    list L;
    L.push_back( X(0));
    L.push_back( X(1));
  
    vector V4(L.begin(), L.end());
    assert( V4.size() == 2);

    vector V7(7, X(7));
    vector V8(7, X(8));
    V7.swap(V8);
    for( i = 0; i< V7.size(); i++) {
        assert( V7[i] == X(8));
    }
    for( i = 0; i< V8.size(); i++) {
        assert( V8[i] == X(7));
    }
    assert( V7.size() == 7);
    it = V7.begin();
    it++;
    V7.insert( it, V8.begin(), V8.end());
    assert( V7.size() == 14);
    V7.insert( V7.begin(), vector::size_type(2), X(9));
    assert( V7.size() == 16);
    std::cout << V7[0] << " " << V7[1] << " " << V7[2] << " " << V7[3]
              << std::endl;
    assert( V7[0] == X(9));
    assert( V7[1] == X(9));
    assert( V7.size() == 16);
    std::cout << "done" << std::endl;

    // test remove_const
    const_iterator mycit=V3.begin();
    iterator myit=cit.remove_const();
    iterator myit2=myit.remove_const();
    CGAL_USE(mycit);
    CGAL_USE(myit);
    CGAL_USE(myit2);

    return 0;
}
