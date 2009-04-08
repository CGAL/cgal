

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


#ifndef CGAL_PDB_BOOST_RTL_SLICEBRNG_HPP_INCLUDED
#define CGAL_PDB_BOOST_RTL_SLICEYRNG_HPP_INCLUDED

#include "boost/config.hpp"
#include <CGAL/PDB/internal/rangelib/priv/defs.hpp>
#include <CGAL/PDB/internal/rangelib/priv/traits.hpp>
#include <CGAL/PDB/internal/rangelib/range.hpp>

#include <boost/iterator/transform_iterator.hpp>
#include <functional>

#include <CGAL/PDB/internal/rangelib/priv/rng_fun.hpp>

// rangelib = Range Library
namespace CGAL { namespace PDB { namespace internal { namespace rangelib {


/*
    FIXME
    sliced() - allow for a transformation that will do this:
    from a given range, slice multiple consecutive records, into only ONE record
    (of a certain type).
    TOTHINK

  allow for two types of slices:
  One that allows to work on one element at a time, and 
  one that allows to work on a range at a time (takes a range, and returns one value)
  Example: you take a string, and return its words (also, quoted text is treated as one word).

  [TOTHINK]
  also, allow to somehow, for a range, to return multiple values.
  For example, if we're generating a daily report for one month.
  For some days, we have no activity (therefore, no records). We would like to
  create some "empty" records to specify that.
  Maybe this could just be another type of range. YUP:
  We'll have a range that takes two elements. Given these, it returns as many other
  elements as possible.
*/

namespace detail {

    // slices, by querying each element once at a time;
    // useful for input iterator ranges (where you cannot access 
    // two elements at a time)
    //
    // VERY IMPORTANT:
    // You should understand that each slice is made of at least ONE element from the 
    // original slice. In case you want intermediate "slices", you should use multirng
    //
    //
    //
    // note: the byrng_function should be a functor taking three arguments
    //       arg1 - the result type (what we return for a slice of the range)
    //       arg2 & arg3 - the iterators that represent a slice
    //       It has to have 'result_type' defined (what you're returning)
    //       result_type must be default-constructible
    template< class r, class byrng_function, class slicer>
    struct slice_byrange_iterator : public ::boost::iterator< 
                                        ::std::input_iterator_tag, 
                                        typename first_argument_finder<byrng_function>::argument_type > {
        typedef slice_byrange_iterator<r,byrng_function, slicer> self_type;

        typedef typename range_finder<r>::range_type range_type;
        typedef typename range_type::iterator iterator_type;
        typedef typename first_argument_finder<byrng_function>::argument_type result_type;

        slice_byrange_iterator( iterator_type first, iterator_type last, byrng_function f, slicer s)
            : m_f( f), m_slicer( s), 
              m_first( first), m_last( last), 
              m_val( result_type() ),
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
        const result_type& operator*() const { return m_val; }
        const result_type* operator->() const { return &m_val; }

        bool more() const { return (m_first != m_last) || !m_processed_last; }
    private:
        void compute_value() {
            if ( m_first == m_last) {
                m_processed_last = true;
                return;
            }

            iterator_type begin = m_first;
            while ( m_first != m_last) {
                iterator_type prev = m_first;
                ++m_first;
				if ( m_first == m_last)
					break;
                if ( !m_slicer( *prev, *m_first))
                    break;
            }
            m_f( m_val, begin, m_first);
        }
    private:
        byrng_function m_f;
        slicer m_slicer;
        iterator_type m_first, m_last;
        result_type m_val;
        // has the last sliced element been processed?
        bool m_processed_last;
    };


    template<class r, class f, class s>
    inline bool operator==( const slice_byrange_iterator<r,f,s> & first, const slice_byrange_iterator<r,f,s> & second) {
        return first.more() == second.more();
    }
    template<class r, class f, class s>
    inline bool operator!=( const slice_byrange_iterator<r,f,s> & first, const slice_byrange_iterator<r,f,s> & second) {
        return first.more() != second.more();
    }

}


template< class r, class byrng_function, class slicer>
struct sliced_byrange_range : public irange< ::CGAL::PDB::internal::rangelib::detail::slice_byrange_iterator<r,byrng_function,slicer> > {
    typedef ::CGAL::PDB::internal::rangelib::detail::slice_byrange_iterator<r,byrng_function,slicer> iterator_type;
    typedef irange< iterator_type> base;

    sliced_byrange_range( r & rng, byrng_function f, slicer s )
    : base( iterator_type( rng.begin(), rng.end(), f, s),
            iterator_type( rng.end(), rng.end(), f, s)) {}

};


#ifndef CGAL_PDB_BOOST_RTL_WORKAROUND_VC6
template< class r, class f, class s> inline sliced_byrange_range<const r,f,s> 
sliced_byrange( const r & rng, f func, s slicer) {
    return sliced_byrange_range<const r,f,s>( rng, func, slicer);
}
#endif

template< class r, class f, class s> inline sliced_byrange_range<r,f,s> 
sliced_byrange( r & rng, f func, s slicer) {
    return sliced_byrange_range<r,f,s>( rng, func, slicer);
}





}}}}


#endif
