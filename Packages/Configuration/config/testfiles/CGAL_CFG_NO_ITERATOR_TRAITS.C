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

//| The class std::iterator_traits is part of the std library.
//| It is used to access certain properties of iterators, such as
//| their value type or iterator category (forward, bidirectional, etc.).
//| The macro CGAL_CFG_NO_ITERATOR_TRAITS is set if std::iterator_traits
//| is not fully supported.

#include <iterator>
#include <vector>

#if defined(__sun) && defined(__SUNPRO_CC)
// For sunpro 5.3 we fake it, since it can do partial specialization
// but ships a non-compliant std library for backwards compatibility.
namespace std {
  template <class Iterator> struct iterator_traits
  {
    typedef _TYPENAME Iterator::value_type value_type;
    typedef _TYPENAME Iterator::difference_type difference_type;
    typedef _TYPENAME Iterator::pointer pointer;
    typedef _TYPENAME Iterator::reference reference;
    typedef _TYPENAME Iterator::iterator_category iterator_category;
  };
  template <class T> struct iterator_traits<T*>
  {
    typedef T value_type;
    typedef ptrdiff_t difference_type;
    typedef T* pointer;
    typedef T& reference;
    typedef random_access_iterator_tag iterator_category;
  };
  template <class T> struct iterator_traits<const T*>
  {
    typedef T value_type;
    typedef ptrdiff_t difference_type;
    typedef const T* pointer;
    typedef const T& reference;
    typedef random_access_iterator_tag iterator_category;
  };
}
#endif // defined(__sun) && defined(__SUNPRO_CC)

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
    return typename iterator_traits<I>::iterator_category();
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

