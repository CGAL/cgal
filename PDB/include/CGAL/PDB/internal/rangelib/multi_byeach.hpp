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

#ifndef CGAL_PDB_BOOST_RTL_MULTI_BYEACH_HPP_INCLUDED
#define CGAL_PDB_BOOST_RTL_MULTI_BYEACH_HPP_INCLUDED

#include <CGAL/PDB/internal/rangelib/priv/defs.hpp>
#include <CGAL/PDB/internal/rangelib/priv/traits.hpp>
#include <boost/iterator.hpp>

#include <CGAL/PDB/internal/rangelib/priv/rng_fun.hpp>

// rangelib = Range Library
namespace CGAL { namespace PDB { namespace internal { namespace rangelib { 

namespace detail {

    template< class r, class multirange_function>
    struct multi_byeach_iterator : public ::boost::iterator< 
                                        ::std::input_iterator_tag, 
                                        typename first_argument_finder<multirange_function>::argument_type > {

        typedef multi_byeach_iterator<r,multirange_function> self_type;

        typedef typename range_finder<r>::range_type range;
        typedef typename range::iterator iterator;
        typedef typename range::value_type src_value_type;
        typedef typename first_argument_finder<multirange_function>::argument_type dest_value_type;

        multi_byeach_iterator( iterator first, iterator last, multirange_function f)
            : m_f( f), 
              m_first( first), m_last( last), 
              m_val(), 
              m_idx( 0),
              m_processed_last( first == last) {
            compute_value();
        }

        self_type& operator++() {
            compute_value();
            return *this;
        }
        self_type operator++(int) {
            self_type tmp( *this);
            ++*this;
            return tmp;
        }
        const dest_value_type& operator*() const { return m_val; }
        const dest_value_type* operator->() const { return &m_val; }

        bool more() const { return (m_first != m_last) || !m_processed_last; }
    private:

        void compute_value() {
            if ( m_first == m_last) {
                m_processed_last = true;
                return;
            }

            if ( m_f( m_val, *m_first, m_idx)) {
                ++m_idx;
                return;
            }
            else
                ++m_first;
        }
    private:
        multirange_function m_f;
        iterator m_first, m_last;
        
        dest_value_type m_val;
        int m_idx;

        // has the last element been processed?
        bool m_processed_last;
    };


    template<class r, class f>
    inline bool operator==( const multi_byeach_iterator<r,f> & first, const multi_byeach_iterator<r,f> & second) {
        return first.more() == second.more();
    }
    template<class r, class f>
    inline bool operator!=( const multi_byeach_iterator<r,f> & first, const multi_byeach_iterator<r,f> & second) {
        return first.more() != second.more();
    }
}


template< class r, class multirange_function>
struct multi_byeach_range : public irange< ::CGAL::PDB::internal::rangelib::detail::multi_byeach_iterator<r,multirange_function> > {
    typedef ::CGAL::PDB::internal::rangelib::detail::multi_byeach_iterator<r,multirange_function> new_iterator;
    typedef irange<new_iterator> base;

    multi_byeach_range( r & rng, multirange_function f)
        : base( new_iterator( rng.begin(), rng.end(), f),
                new_iterator( rng.end(), rng.end(), f) ) {}
};


#ifndef CGAL_PDB_BOOST_RTL_WORKAROUND_VC6
template<class r, class f>
multi_byeach_range<const r,f> multied_byeach( const r & rng, f func) {
    return multi_byeach_range<const r,f>(rng,func);
}
#endif

template<class r, class f>
multi_byeach_range<r,f> multied_byeach( r & rng, f func) {
    return multi_byeach_range<r,f>(rng,func);
}




}}}}

#endif
