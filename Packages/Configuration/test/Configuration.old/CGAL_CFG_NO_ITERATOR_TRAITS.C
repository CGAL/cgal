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
// release       : $CGAL_Revision: CGAL-0.9-I-03 $
// release_date  : $CGAL_Date: 1997/12/01 $
//
// file          : config/testfiles/CGAL_CFG_NO_ITERATOR_TRAITS.C
// source        :
// revision      : 1.11
// revision_date : 29 Mar 1998
// author(s)     : various
//
// coordinator   : Utrecht University
//
// ============================================================================

// CGAL_CFG_NO_ITERATOR_TRAITS.C
// ---------------------------------------------------------------------
// A short test program to evaluate a C++ compiler.
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| Iterator traits are documented in the Dec. 1996 C++ Standard draft.
//| The following definition is set if iterator are not fully supported 
//| including their use in a template class, as a default template
//| argument and as a return type of global function.

#include <assert.h>
#include <iterator.h>
#include <vector.h>

// This class implements an iterator adaptor that forwards all
// member function calls to its template argument. It uses 
// iterator traits to derive correct types and iterator category.

template < class I, class category = iterator_traits<I>::iterator_category>
class Adaptor {
    I _i;
public:
    typedef iterator_traits<I>::value_type       value_type;
    typedef iterator_traits<I>::difference_type  difference_type;
    typedef iterator_traits<I>::reference        reference;
    typedef iterator_traits<I>::pointer          pointer;
    typedef category                             iterator_category;
    Adaptor( I i) : _i(i) {}
    Adaptor<I, category>&
    operator++() {
	++_i;
	return *this;
    }
    reference
    operator*() const {
	return *_i;
    }
};

// A global function to extract the iterator category.

template < class I> inline
iterator_traits<I>::iterator_category
query( I i) {
    return iterator_traits<I>::iterator_category();
}

// A function to match bidirectional iterators.
inline
int discr( bidirectional_iterator_tag tag) { return 42; }

int main() {
    vector<int> v;
    v.push_back(32);
    v.push_back(33);
    v.push_back(42);
    Adaptor< vector<int>::iterator> i( v.begin());
    ++i;
    assert( *i == 33);
    ++i;
    assert( *i == 42);
    assert( discr( query( i)) == 42);
    return 0;
}

// EOF //
