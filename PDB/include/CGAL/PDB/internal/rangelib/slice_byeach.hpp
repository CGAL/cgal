

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


#ifndef CGAL_PDB_BOOST_RTL_SLICEBYEACH_HPP_INCLUDED
#define CGAL_PDB_BOOST_RTL_SLICEBYEACH_HPP_INCLUDED

#include "boost/config.hpp"
#include <CGAL/PDB/internal/rangelib/priv/defs.hpp>
#include <CGAL/PDB/internal/rangelib/priv/traits.hpp>
#include <CGAL/PDB/internal/rangelib/range.hpp>

#include <boost/iterator.hpp>
#include <functional>

#include <CGAL/PDB/internal/rangelib/priv/rng_fun.hpp>

// rangelib = Range Library
namespace CGAL { namespace PDB { namespace internal { namespace rangelib {


namespace detail {

    // slices, by querying each element once at a time;
    // useful for input iterator ranges (where you cannot access 
    // two elements at a time)
    //
    // VERY IMPORTANT:
    // You should understand that each sclice is made of at least ONE element from the 
    // original slice. In case you want intermediate "slices", you should use multi_range
    //
    //
    //
    // note: the byeach_function should be a binary functor, taking 2 arguments:
    //       arg1 - the result type (what we return for a slice of the range)
    //       arg2 - one value from the slice 
    //       It has to have the 'first_argument_type' defined (a typedef, most likely)
    //
    //  if f is a byeach_function, then
    //  f( result, val);
    // will recompute the result based on val 
    // For each value 'val' from a slice, we'll have a call to f(result,val);
    template< class r, class byeach_function, class slicer>
    struct slice_byeach_iterator : public ::boost::iterator< 
                                        ::std::input_iterator_tag, 
                                        typename first_argument_finder<byeach_function>::argument_type > {
        typedef slice_byeach_iterator<r,byeach_function, slicer> self_type;

        typedef typename range_finder<r>::range_type range;
        typedef typename range::iterator iterator;
        typedef typename range::value_type src_value_type;
        typedef typename first_argument_finder<byeach_function>::argument_type dest_value_type;

        slice_byeach_iterator( iterator first, iterator last, byeach_function f, slicer s, dest_value_type init)
            : m_f( f), m_slicer( s), 
              m_first( first), m_last( last), 
              m_init( init), m_val( init),
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

            // FIXME for iterator categories > input_iterator
            // don't use any temporary
            m_val = m_init;
            while ( true) {
                src_value_type tmp = *m_first;
                m_f( m_val, tmp);
                ++m_first;
                if ( m_first != m_last) {
                    if ( !m_slicer( tmp, *m_first))
                        // got to a new slice
                        break;
                }
                else
                    // got to the end
                    break;
            }
        }
    private:
        byeach_function m_f;
        slicer m_slicer;
        iterator m_first, m_last;
        dest_value_type m_init, m_val;
        // has the last sliced element been processed?
        bool m_processed_last;
    };


    template<class r, class f, class s>
    inline bool operator==( const slice_byeach_iterator<r,f,s> & first, const slice_byeach_iterator<r,f,s> & second) {
        return first.more() == second.more();
    }
    template<class r, class f, class s>
    inline bool operator!=( const slice_byeach_iterator<r,f,s> & first, const slice_byeach_iterator<r,f,s> & second) {
        return first.more() != second.more();
    }

}


template< class r, class byeach_function, class slicer>
struct sliced_byeach_range : public irange< ::CGAL::PDB::internal::rangelib::detail::slice_byeach_iterator<r,byeach_function,slicer> > {
    typedef ::CGAL::PDB::internal::rangelib::detail::slice_byeach_iterator<r,byeach_function,slicer> iterator_type;
    typedef irange< iterator_type> base;

    typedef typename iterator_type::dest_value_type v_type;

    sliced_byeach_range( r & rng, byeach_function f, slicer s, v_type init = v_type() )
        : base( iterator_type( rng.begin(), rng.end(), f, s, init),
                iterator_type( rng.end(), rng.end(), f, s, init)) {}

};


#ifndef CGAL_PDB_BOOST_RTL_WORKAROUND_VC6
template< class r, class f, class s, class initializer> inline sliced_byeach_range<const r,f,s> 
sliced_byeach( const r & rng, f func, s slicer, initializer init) {
    return sliced_byeach_range<const r,f,s>( rng, func, slicer, init);
}
template< class r, class f, class s> inline sliced_byeach_range<const r,f,s> 
sliced_byeach( const r & rng, f func, s slicer) {
    return sliced_byeach_range<const r,f,s>( rng, func, slicer);
}
#endif

template< class r, class f, class s, class initializer> inline sliced_byeach_range<r,f,s> 
sliced_byeach( r & rng, f func, s slicer, initializer init) {
    return sliced_byeach_range<r,f,s>( rng, func, slicer, init);
}
template< class r, class f, class s> inline sliced_byeach_range<r,f,s> 
sliced_byeach( r & rng, f func, s slicer) {
    return sliced_byeach_range<r,f,s>( rng, func, slicer);
}







}}}}


#endif
