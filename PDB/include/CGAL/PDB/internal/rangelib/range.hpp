// Boost. Iterable Range Library (rangelib)
//
// Copyright 2003-2004 John Torjo (john@torjo.com) and Matthew Wilson (matthew@synesis.com.au)
//
// Permission to copy, use, sell and distribute this software is granted
// provided this copyright notice appears in all copies.
// Permission to modify the code and to distribute modified code is granted
// provided this copyright notice appears in all copies, and a notice
// that the code was modified is included with the copyright notice.
//
// This software is provided "as is" without express or implied warranty,
// and with no claim as to its suitability for any purpose.
 
// See http://www.boost.org for updates, documentation, and revision history.


#ifndef CGAL_PDB_BOOST_RTL_HPP_INCLUDED
#define CGAL_PDB_BOOST_RTL_HPP_INCLUDED

#include <CGAL/PDB/internal/rangelib/priv/traits.hpp>
#include <utility> // for std::pair


// rangelib = Range Library
namespace CGAL { namespace PDB { namespace internal { namespace rangelib { 

/*
    Purpose:
    To save you as many keystrokes as possible, when dealing with 
    iterators/containers/algorithms,
    so that you will write compact code that is easy and straightforward.
*/
    
// FIXME - create helpers for containers, like, to very simply access key-only or value-only


// range for containers
template< class container_type>
struct crange : public detail::container_traits< container_type> {

private:
	typedef detail::container_traits< container_type>	parent_class_type;
public:

	typedef typename parent_class_type::iterator	iterator;
	typedef typename parent_class_type::reference	reference;
	typedef typename parent_class_type::pointer		pointer;
    typedef crange< container_type>					self_type; // this class

    // the type does not matter; I just want to know if a given type is a range or not
    struct is_range_type {};

    // * Construction
    template< class cont>
    crange( cont & c) 
        : m_first( c.begin() ), m_last( c.end() ) {
    }
#ifndef CGAL_PDB_BOOST_RTL_WORKAROUND_VC6
    template< class cont>
    crange( const cont & c) 
        : m_first( c.begin() ), m_last( c.end() ) {
    }
#endif

    crange( iterator first, iterator last) 
        : m_first( first), m_last( last) {
    }

    crange( ::std::pair<iterator,iterator> val) 
        : m_first( val.first), m_last( val.second) {
    }


    // conversion
    template< class r>
	crange( const crange<r> & other) 
        : m_first( other.begin()), m_last( other.end()) {
    }


    // * Range access

    iterator begin() const { return m_first; }
    iterator end() const { return m_last; }
    void begin( iterator first) { m_first = first; }
    void end( iterator last) { m_last = last; }

    /*
    Use rng::advanced() instead!

    // useful for algorithms that require some value in the middle of the range,
    // like, nth_element or partial_sort
    template< class Dist> iterator advanced(Dist n) const {
        iterator i;
        ::std::advance( i, n);
        return i;
    }
    */

    // * Operators

    self_type & operator=( const self_type & other) {
        m_first = other.m_first;
        m_last = other.m_last;
        return *this;
    }

    reference operator*() const { return *m_first; }
    pointer operator->() const { return &*m_first; }

    // FIXME - in the future, maybe allow checks on ranges (are they valid?)

    // is the range valid (or we've reached the end) ?
    operator bool() const {
        return m_first != m_last;
    }

    self_type & operator++() { 
        ++m_first;
        return *this;
    }

    self_type operator++(int) {
        self_type tmp = *this;
        ++*this;
        return tmp;
    }

    self_type & operator--() { 
        --m_first;
        return *this;
    }

    self_type operator--(int) {
        self_type tmp = *this;
        --*this;
        return tmp;
    }

    self_type & operator+=( size_t n) {
        m_first += n;
        return *this;
    }
    self_type & operator-=( size_t n) {
        m_first += n;
        return *this;
    }

    reference operator[]( size_t n) {
        return m_first[n];
    }

    // * Operators not implemented

    // Note: it's dangerous to compare two ranges, and it makes no point
    // in doing so. 
    // * it's not very straightforward when two ranges should be equal.
    //   (should r1.begin() == r2.begin() suffice?)
    // * if we compare r1.begin() == r2.begin() && r1.end() == r2.end(),
    //   it might be inefficient. 
    //
    // You can, however, compare two begin() iterators.
    // You can do so, at will, very simple:
    // if ( r1.begin() == r2.begin() )...
    // (same for end)
    bool operator !=( const self_type&) const;
    bool operator ==( const self_type&) const;
    bool operator <( const self_type&) const;
    bool operator <=( const self_type&) const;
    bool operator >( const self_type&) const;
    bool operator >=( const self_type&) const;
    void operator -( const self_type&) const;
    
private:
    iterator m_first, m_last;
};

// * non-member Operators
template< class r> inline crange<r> operator+( crange<r> val, size_t n) {
    return val += n;
}
template< class r> inline crange<r> operator+( size_t n, crange<r> val) {
    return val += n;
}

template< class r> inline crange<r> operator-( crange<r> val, size_t n) {
    return val -= n;
}

template< class r> inline crange<r> operator-( size_t n, crange<r> val) {
    return val -= n;
}


namespace detail {
    template< class iterator_type>
    struct iterator_pair : public ::CGAL::PDB::internal::rangelib::detail::iterator_traits< iterator_type> {

private:
	typedef ::CGAL::PDB::internal::rangelib::detail::iterator_traits< iterator_type>	parent_class_type;
public:

	typedef typename parent_class_type::iterator	iterator;

        iterator begin() const { return m_first; }
        iterator end() const { return m_last; }
        iterator_pair( iterator first, iterator last) : m_first( first), m_last( last) {}
    private:
        iterator m_first, m_last;
    };
} // detail


// range for iterators
template< class iterator_type> 
struct irange : public crange< detail::iterator_pair< iterator_type> > {
    typedef detail::iterator_pair< iterator_type> pair_type;
    typedef crange< pair_type > base;
    typedef irange< iterator_type> self_type;
    // end has a default value, for facilitating the 'default iterator is passed-end iterator' idiom
    irange( iterator_type begin, iterator_type end = iterator_type()) 
        : base( begin,end ) {}

#if !defined(BOOST_MSVC) || \
	_MSC_VER >= 1300
	template <typename T, size_t N>
    irange( T (&ar)[N]) 
        : base( pair_type(&ar[0], &ar[N]) )
	{}
#endif /* _MSC_VER >= 1300 */

    irange( const ::std::pair<iterator_type,iterator_type> & both) 
        : base( pair_type(both.first,both.second)) {
    }

    self_type & operator=( const self_type & other) {
        base::operator=( other);
        return *this;
    }
};


// * converter functions - c_range & i_range
//   (so that you don't have to create an irange<> unless necesary)

template< class container> 
inline crange< container> c_range( container & c) {
    return crange< container>( c);
}

template< class container> 
inline crange< const container> c_range( const container & c) {
    return crange< const container>( c);
}

template< class iterator> 
inline irange< iterator> i_range( iterator first, iterator last) {
    return irange< iterator>( first, last);
}

// C-like range
template< class iterator> 
inline irange< iterator> i_range( iterator first, int n) {
    // only works for random iterators (where operator + is implemented)
    return irange< iterator>( first, first + n);
}




}}}}


#endif // CGAL_PDB_BOOST_RTL_HPP_INCLUDED
