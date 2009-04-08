
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


#ifndef CGAL_PDB_BOOST_RTL_REVERSE_HPP_INCLUDED
#define CGAL_PDB_BOOST_RTL_REVERSE_HPP_INCLUDED

#include <CGAL/PDB/internal/rangelib/priv/defs.hpp>
#include <CGAL/PDB/internal/rangelib/priv/traits.hpp>
#include <CGAL/PDB/internal/rangelib/range.hpp>

#include <boost/iterator/reverse_iterator.hpp>


// rangelib = Range Library
namespace CGAL { namespace PDB { namespace internal { namespace rangelib {

namespace detail {

    template< class r> 
    struct reverse_iterator_finder {
        typedef typename ::CGAL::PDB::internal::rangelib::range_finder<r>::range_type r_type;
        typedef typename r_type::iterator i_type;

        typedef ::boost::reverse_iterator<i_type> reverse_type;
    };
}


template< class r > 
struct reversed_range : public irange< typename ::CGAL::PDB::internal::rangelib::detail::reverse_iterator_finder<r>::reverse_type > {
    typedef irange< typename ::CGAL::PDB::internal::rangelib::detail::reverse_iterator_finder<r>::reverse_type > base;
    typedef typename ::CGAL::PDB::internal::rangelib::detail::reverse_iterator_finder<r>::i_type old_iterator;
    typedef typename ::CGAL::PDB::internal::rangelib::detail::reverse_iterator_finder<r>::reverse_type new_iterator;

    reversed_range( old_iterator first, old_iterator last)
        : base( new_iterator(last), new_iterator(first) ) {
    }

    reversed_range( r & rng) 
        : base( new_iterator( rng.end() ), 
                new_iterator( rng.begin() ) ) {
    }
};


#ifndef CGAL_PDB_BOOST_RTL_WORKAROUND_VC6
template< class r> inline reversed_range<const r> 
reversed( const r & rng) {
    return reversed_range<const r>( rng.begin(), rng.end() );
}
#endif

template< class r> inline reversed_range<r> 
reversed( r & rng) {
    return reversed_range<r>( rng.begin(), rng.end() );
}



}}}}


#endif

