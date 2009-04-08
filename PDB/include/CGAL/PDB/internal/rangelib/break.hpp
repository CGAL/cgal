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


#ifndef CGAL_PDB_BOOST_RTL_break_HPP_INCLUDED
#define CGAL_PDB_BOOST_RTL_break_HPP_INCLUDED

#include <CGAL/PDB/internal/rangelib/priv/defs.hpp>
#include <CGAL/PDB/internal/rangelib/priv/traits.hpp>
#include <CGAL/PDB/internal/rangelib/range.hpp>

#include <boost/iterator.hpp>

// FIXME port everything to work with Thorsten's range_traits lib.

// rangelib = Range Library
namespace CGAL { namespace PDB { namespace internal { namespace rangelib {

namespace detail {
    template<class iterator_type, class pred> 
    struct break_iterator :  public ::boost::iterator< 
                                        ::std::input_iterator_tag, 
                                        typename ::CGAL::PDB::internal::rangelib::detail::iterator_traits<iterator_type>::value_type> {
        typedef typename ::CGAL::PDB::internal::rangelib::detail::iterator_traits<iterator_type>::iterator i_type;
        typedef typename ::CGAL::PDB::internal::rangelib::detail::iterator_traits<iterator_type>::value_type v_type;
        typedef typename ::CGAL::PDB::internal::rangelib::detail::iterator_traits<iterator_type>::reference r_type;
        typedef typename ::CGAL::PDB::internal::rangelib::detail::iterator_traits<iterator_type>::pointer p_type;

        typedef break_iterator<iterator_type,pred> self_type;

        break_iterator(i_type iter, pred p, bool breaked) : m_iterator(iter), m_p(p), m_breaked(breaked) {
            if ( !m_breaked)
                check_breaked();
        }

        self_type& operator++() {
            ++m_iterator;
            check_breaked();
            return *this;
        }
        self_type operator++(int) {
            self_type tmp( *this);
            ++*this;
            return tmp;
        }
        r_type operator*() { return *m_iterator; }
        p_type operator->() { return &*m_iterator; }


        inline bool operator==( const self_type & other) const {
            if ( m_iterator == other.m_iterator) return true;
            return m_breaked == other.m_breaked; // see if they reached the end...
        }

    private:
        // see if we've breaked...
        void check_breaked() {
            if ( m_breaked) return;
            m_breaked = !m_p(*m_iterator);
        }

    private:
        pred m_p;
        i_type m_iterator;
        bool m_breaked;
    };

    template<class i, class p>
    inline bool operator!=( const break_iterator<i,p> & first, const break_iterator<i,p> & second) {
        return !(first == second);
    }




    template< class r, class pred> struct break_iterator_finder {
        typedef typename ::CGAL::PDB::internal::rangelib::range_finder<r>::range_type r_type;
        typedef typename r_type::iterator i_type;
        typedef ::CGAL::PDB::internal::rangelib::detail::break_iterator<i_type, pred> break_type;
    };
}


/* 
    Wraps a range or a container
    When a certain functor returns false, it "breaks".

    It is equivalent to:
    for ( some_range r(...); r; ++r)
        if ( some_functor(*r) )
            do_something( *r);
        else
            break;  // **********

*/
template< class r, class pred>
struct breaked_range : public irange< typename ::CGAL::PDB::internal::rangelib::detail::break_iterator_finder<r, pred>::break_type > {
    typedef irange< typename ::CGAL::PDB::internal::rangelib::detail::break_iterator_finder<r, pred>::break_type > base;
    typedef typename ::CGAL::PDB::internal::rangelib::detail::break_iterator_finder<r, pred>::i_type old_iterator;
    typedef typename ::CGAL::PDB::internal::rangelib::detail::break_iterator_finder<r, pred>::break_type new_iterator;

    breaked_range( old_iterator first, old_iterator last, pred p = pred() )
        : base( new_iterator(first,p,false), new_iterator(last,p,true) ) {
    }

    breaked_range( r & rng, pred p = pred() ) 
        : base( new_iterator( rng.begin(),p,false), 
                new_iterator( rng.end(),p,true) ) {
    }
};


#ifndef CGAL_PDB_BOOST_RTL_WORKAROUND_VC6
template< class r, class pred> inline breaked_range<const r,pred> 
breaked( const r & rng, pred p) {
    return breaked_range<const r,pred>( rng.begin(), rng.end(), p);
}
#endif

template< class r, class pred> inline breaked_range<r,pred> 
breaked( r & rng, pred p) {
    return breaked_range<r,pred>( rng.begin(), rng.end(), p);
}


// helper for breaked() - it breaks after n elements.
struct bk_after {
    bk_after(int n) : m_n(n), m_idx(0) {}
    template<class value_type> bool operator()(const value_type&) {
        return m_idx++ < m_n;
    }
private:
    int m_n, m_idx;
};
typedef bk_after break_after; // longer name, if you wish ;)


}}}}


#endif


/**
FIXME
have predicate for stoppping.
break_at()

- make it simple to specify something like break_after(n) - breaks after n elements

*/