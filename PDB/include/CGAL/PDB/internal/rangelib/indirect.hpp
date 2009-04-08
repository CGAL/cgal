
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


#ifndef CGAL_PDB_BOOST_RTL_INDIRECT_HPP_INCLUDED
#define CGAL_PDB_BOOST_RTL_INDIRECT_HPP_INCLUDED

#include <CGAL/PDB/internal/rangelib/priv/defs.hpp>
#include <CGAL/PDB/internal/rangelib/priv/traits.hpp>
#include <CGAL/PDB/internal/rangelib/range.hpp>

#include <boost/iterator/indirect_iterator.hpp>


// rangelib = Range Library
namespace CGAL { namespace PDB { namespace internal { namespace rangelib {

namespace detail {
    struct no_type {};
    template< class x> struct is_no_type_finder { BOOST_STATIC_CONSTANT(bool, is = false); };
    template<> struct is_no_type_finder<no_type> { BOOST_STATIC_CONSTANT(bool, is = true); };

    template< class r, class v_type> 
    struct indirect_iterator_finder {
        BOOST_STATIC_CONSTANT(bool, is_notype = ::CGAL::PDB::internal::rangelib::detail::is_no_type_finder<v_type>::is ); 
        
        typedef typename ::CGAL::PDB::internal::rangelib::range_finder<r>::range_type r_type;
        typedef typename r_type::iterator i_type;

        typedef typename ::boost::mpl::if_c< is_notype, 
            ::boost::indirect_iterator<i_type>,
            ::boost::indirect_iterator<i_type, v_type>
            >::type indirect_type;

    };
}


template< class r, class v_type = ::CGAL::PDB::internal::rangelib::detail::no_type > 
struct indirected_range : public irange< typename ::CGAL::PDB::internal::rangelib::detail::indirect_iterator_finder<r,v_type>::indirect_type > {
    typedef irange< typename ::CGAL::PDB::internal::rangelib::detail::indirect_iterator_finder<r,v_type>::indirect_type > base;
    typedef typename ::CGAL::PDB::internal::rangelib::detail::indirect_iterator_finder<r,v_type>::i_type old_iterator;
    typedef typename ::CGAL::PDB::internal::rangelib::detail::indirect_iterator_finder<r,v_type>::indirect_type new_iterator;


    indirected_range( old_iterator first, old_iterator last)
        : base( new_iterator(first), new_iterator(last) ) {
    }

    indirected_range( r & rng) 
        : base( new_iterator( rng.begin() ), 
                new_iterator( rng.end() ) ) {
    }
};


#ifndef CGAL_PDB_BOOST_RTL_WORKAROUND_VC6
template< class r> inline indirected_range<const r> 
indirected( const r & rng) {
    return indirected_range<const r>( rng.begin(), rng.end());
}
#endif

template< class r> inline indirected_range<r> 
indirected( r & rng) {
    return indirected_range<r>( rng.begin(), rng.end());
}

// here, you specify the value_type directly
#ifndef CGAL_PDB_BOOST_RTL_WORKAROUND_VC6
template< class v_type, class r> inline indirected_range<const r, v_type> 
indirected( const r & rng, v_type* = 0) {
    return indirected_range<const r, v_type>( rng.begin(), rng.end());
}
#endif

template< class v_type, class r> inline indirected_range<r, v_type> 
indirected( r & rng, v_type* = 0) {
    return indirected_range<r, v_type>( rng.begin(), rng.end());
}



}}}}


#endif
