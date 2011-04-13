// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  :
//
// file          : config/testfiles/CGAL_CFG_NO_ITERATOR_TRAITS.C
// source        :
// revision      : 1.11
// revision_date : 29 Mar 1998
// author(s)     : various
//
// coordinator   : Utrecht University
//
// ======================================================================

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

#include <iterator>
#include <vector>

// This class implements an iterator adaptor that forwards all
// member function calls to its template argument. It uses 
// iterator traits to derive correct types and iterator category.

#ifdef  _MSC_VER
using namespace std; // MSC hates "using std::{blah};"....
#  if _MSC_VER < 1300
#    define typename     // preventing MSVC 6.0 "error C2899:
                     // typename cannot be used outside a template
#  endif  // MSVC 6.0
#else
using std::iterator_traits;
#endif // _MSC_VER

template < class I,
           class category = typename iterator_traits<I>::iterator_category >
class Adaptor {
    I _i;
public:
    typedef typename iterator_traits<I>::value_type       value_type;
    typedef typename iterator_traits<I>::difference_type  difference_type;
    typedef typename iterator_traits<I>::reference        reference;
    typedef typename iterator_traits<I>::pointer          pointer;
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
typename iterator_traits<I>::iterator_category
query( I) {
    return iterator_traits<I>::iterator_category();
}

// A function to match bidirectional iterators.
inline
int discr( std::bidirectional_iterator_tag ) { return 42; }

bool all_assertions_correct = true;

int main() {
    std::vector<int> v;
    v.push_back(32);
    v.push_back(33);
    v.push_back(42);
    Adaptor< std::vector<int>::iterator> i( v.begin());
    ++i;
    all_assertions_correct &= ( *i == 33);
    ++i;
    all_assertions_correct &= ( *i == 42);
    all_assertions_correct &= ( discr( query( i)) == 42);
    return !all_assertions_correct;
}

