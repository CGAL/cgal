
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


#ifndef CGAL_PDB_BOOST_RTL_RESORT_HPP_INCLUDED
#define CGAL_PDB_BOOST_RTL_RESORT_HPP_INCLUDED

#include <boost/config.hpp>

#include <CGAL/PDB/internal/rangelib/priv/defs.hpp>
#include <CGAL/PDB/internal/rangelib/priv/traits.hpp>
#include <CGAL/PDB/internal/rangelib/range.hpp>
#include <CGAL/PDB/internal/rangelib/algo.hpp>
#include <CGAL/PDB/internal/rangelib/indirect.hpp>

#include <boost/shared_ptr.hpp>
#include <vector>
#include <functional>

// rangelib = Range Library
namespace CGAL { namespace PDB { namespace internal { namespace rangelib {

namespace detail {
    template< class r> 
    struct resort_iterator_finder {
        typedef typename ::CGAL::PDB::internal::rangelib::range_finder<r>::range_type r_type;
        typedef typename r_type::iterator i_type;
        typedef typename r_type::value_type v_type;
        typedef typename r_type::pointer pointer;
        typedef typename r_type::reference reference;
        
        typedef typename ::std::vector<pointer>::iterator indirected_i;
        typedef ::boost::indirect_iterator<indirected_i > resort_type;
    };
}


// FIXME - I should find a way to disallow shrinking;
// something like, disallow this:
// crange<whatever> r = resorted(x);
//
// the above is a bug, because the iterators returned from resorted range point to 
// an INTERNAL vector, which is destructed at resorted_range's destruction
//
// TOTHINK - perhaps the best way is to simply disallow copy-construction???


// FIXME - the range should at least contain forward iterators
template< class r> 
struct resorted_range : public irange< typename ::CGAL::PDB::internal::rangelib::detail::resort_iterator_finder<r>::resort_type > {

    typedef irange< typename ::CGAL::PDB::internal::rangelib::detail::resort_iterator_finder<r>::resort_type> base;
    typedef typename ::CGAL::PDB::internal::rangelib::detail::resort_iterator_finder<r>::i_type old_iterator;
    typedef typename ::CGAL::PDB::internal::rangelib::detail::resort_iterator_finder<r>::resort_type new_iterator;

    typedef typename ::CGAL::PDB::internal::rangelib::detail::resort_iterator_finder<r>::v_type v_type;
    typedef typename ::CGAL::PDB::internal::rangelib::detail::resort_iterator_finder<r>::pointer pointer;
    typedef typename ::CGAL::PDB::internal::rangelib::detail::resort_iterator_finder<r>::reference reference;
    typedef std::vector< pointer> array;

    resorted_range( old_iterator first, old_iterator last)
        : base( new_iterator(), new_iterator() ),
          m_resorted( new array) {
        typedef ::std::less< v_type> default_sort;
        resort( first, last, default_sort() );
    }

    template< class pred>
    resorted_range( old_iterator first, old_iterator last, pred p)
        : base( new_iterator(), new_iterator() ),
          m_resorted( new array) {
        resort( first, last, p );
    }
    
    resorted_range( r & rng) 
        : base( new_iterator(), 
                new_iterator() ),
          m_resorted( new array) {
        typedef ::std::less< v_type> default_sort;

        resort( rng.begin(), rng.end(), default_sort());
    }

    template< class pred>
    resorted_range( r & rng, pred p ) 
        : base( new_iterator(), 
                new_iterator() ),
          m_resorted( new array) {
        resort( rng.begin(), rng.end(), p );
    }

private:
    struct resorted_inserter {
        resorted_inserter( array & v) : m_v( v) {}
        void operator()( reference val) {
            m_v.push_back( &val);
        }
    private:
        array & m_v;
    };

    template< class pred>
    struct dereference_ptr {
    dereference_ptr( pred p ) : m_p(p) {}
    bool operator() (pointer first, pointer second) {
        return m_p( *first, *second);
    }
    private:
        pred m_p;
    };

    template<class pred>
    void resort( old_iterator first, old_iterator last, pred p) {
        m_resorted->reserve( ::std::distance(first, last) );
        ::std::for_each( first, last, resorted_inserter( *m_resorted) );
        ::CGAL::PDB::internal::rangelib::rng::sort( *m_resorted, dereference_ptr<pred>(p) );

        // reset the iterators
        begin( m_resorted->begin() );
        end( m_resorted->end() ); 
    }
private:
    ::boost::shared_ptr<array> m_resorted;
};


// resorts using "operator<"
#ifndef CGAL_PDB_BOOST_RTL_WORKAROUND_VC6
template< class r> inline resorted_range<const r> 
resorted( const r & rng) {
    return resorted_range<const r>( rng.begin(), rng.end() );
}
#endif

template< class r> inline resorted_range<r> 
resorted( r & rng) {
    return resorted_range<r>( rng.begin(), rng.end() );
}


// using binary predicate
#ifndef CGAL_PDB_BOOST_RTL_WORKAROUND_VC6
template< class r, class pred> inline resorted_range<const r> 
resorted( const r & rng, pred p) {
    return resorted_range<const r>( rng.begin(), rng.end(), p );
}
#endif

template< class r, class pred> inline resorted_range<r> 
resorted( r & rng, pred p) {
    return resorted_range<r>( rng.begin(), rng.end(), p );
}




}}}}


#endif
