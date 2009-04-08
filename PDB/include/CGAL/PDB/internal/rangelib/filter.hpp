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


#ifndef CGAL_PDB_BOOST_RTL_FILTER_HPP_INCLUDED
#define CGAL_PDB_BOOST_RTL_FILTER_HPP_INCLUDED

#include <CGAL/PDB/internal/rangelib/priv/defs.hpp>
#include <CGAL/PDB/internal/rangelib/priv/traits.hpp>
#include <CGAL/PDB/internal/rangelib/range.hpp>

#include <boost/iterator/filter_iterator.hpp>


// rangelib = Range Library
namespace CGAL { namespace PDB { namespace internal { namespace rangelib {

namespace detail {
    template< class r, class pred> 
    struct filter_iterator_finder {
        typedef typename ::CGAL::PDB::internal::rangelib::range_finder<r>::range_type r_type;
        typedef typename r_type::iterator i_type;
        typedef ::boost::filter_iterator<pred, i_type> filter_type;
    };
}


// filters a range or a container
template< class r, class pred>
struct filtered_range : public irange< typename ::CGAL::PDB::internal::rangelib::detail::filter_iterator_finder<r, pred>::filter_type > {
    typedef irange< typename ::CGAL::PDB::internal::rangelib::detail::filter_iterator_finder<r, pred>::filter_type > base;
    typedef typename ::CGAL::PDB::internal::rangelib::detail::filter_iterator_finder<r, pred>::i_type old_iterator;
    typedef typename ::CGAL::PDB::internal::rangelib::detail::filter_iterator_finder<r, pred>::filter_type new_iterator;

    filtered_range( old_iterator first, old_iterator last, pred p = pred() )
        : base( new_iterator(p,first,last), new_iterator(p,last,last) ) {
    }

    filtered_range( r & rng, pred p = pred() ) 
        : base( new_iterator(p,rng.begin(),rng.end()), 
                new_iterator(p,rng.end(),rng.end() ) ) {
    }
};


#ifndef CGAL_PDB_BOOST_RTL_WORKAROUND_VC6
template< class r, class pred> inline filtered_range<const r,pred> 
filtered( const r & rng, pred p) {
    return filtered_range<const r,pred>( rng.begin(), rng.end(), p);
}
#endif

template< class r, class pred> inline filtered_range<r,pred> 
filtered( r & rng, pred p) {
    return filtered_range<r,pred>( rng.begin(), rng.end(), p);
}


}}}}


#endif
