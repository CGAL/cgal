
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


#ifndef CGAL_PDB_BOOST_RTL_GENERATE_HPP_INCLUDED
#define CGAL_PDB_BOOST_RTL_GENERATE_HPP_INCLUDED

#include <CGAL/PDB/internal/rangelib/priv/defs.hpp>
#include <CGAL/PDB/internal/rangelib/priv/traits.hpp>
#include <CGAL/PDB/internal/rangelib/range.hpp>

#include <boost/iterator.hpp>
#include <functional>

// rangelib = Range Library
namespace CGAL { namespace PDB { namespace internal { namespace rangelib {

namespace detail {
    struct end_iterator {};

    template< class res>
    struct generator_functor {
        typedef res (*func)();
        typedef res result_type;
        generator_functor( func f) : m_f(f) {}
        result_type operator() () { return m_f(); }
    private:
        func m_f;
    };

    // note: result_type must be Assignable and copy-constructible
    //
    // the stopper will return false when we should stop
    template< class generator, class stopper>
    struct generated_iterator : public ::boost::iterator< 
                                        ::std::input_iterator_tag, typename generator::result_type> {
        typedef generated_iterator<generator, stopper> self_type;
        typedef typename generator::result_type result_type;

        generated_iterator( generator g, stopper s)
            : m_g( g), m_s( s), m_more( true), 
              m_val( m_g() ) {
            check_more();
        }

        generated_iterator( generator g, stopper s, const end_iterator & )
            : m_g( g), m_s( s), m_more( false), 
              // initializing to this, since the result_type is not required to be
              // default constructible
              m_val( m_g() ) {
        }

        self_type& operator++() {
            m_val = m_g();
            check_more();
            return *this;
        }
        self_type operator++(int) {
            self_type tmp( *this);
            ++*this;
            return tmp;
        }
        const result_type& operator*() const { return m_val; }
        const result_type* operator->() const { return &m_val; }

        bool more() const { return m_more; }
    private:
        void check_more() {
            m_more = m_s( m_val); // see if reached the end
        }

    private:
        generator m_g;
        stopper m_s;
        bool m_more;
        // must be INITIALIZED Last (because it depends on m_g)
        result_type m_val;
    };

    template<class g, class s>
    inline bool operator==( const generated_iterator<g,s> & first, const generated_iterator<g,s> & second) {
        return first.more() == second.more();
    }
    template<class g, class s>
    inline bool operator!=( const generated_iterator<g,s> & first, const generated_iterator<g,s> & second) {
        return first.more() != second.more();
    }
}


// represents a range, generated from a functor (generated)
//
// see the ::std::generate/ ::std::generate_n algorithms
template< class generator, class stopper>
struct generated_range : public irange< ::CGAL::PDB::internal::rangelib::detail::generated_iterator<generator,stopper> > {
    typedef ::CGAL::PDB::internal::rangelib::detail::generated_iterator<generator,stopper> iterator_type;
    typedef irange< iterator_type> base;
    generated_range( generator g, stopper s)
        : base( iterator_type(g, s),
                iterator_type(g, s, ::CGAL::PDB::internal::rangelib::detail::end_iterator()) ) {
    }
};


template< class generator, class stopper>
generated_range<generator,stopper> generated( generator g, stopper s) {
    return generated_range<generator,stopper>( g,s);
}

template< class res, class stopper>
generated_range< ::CGAL::PDB::internal::rangelib::detail::generator_functor<res>, stopper> generated_f( res (*g)(), stopper s) {
    typedef ::CGAL::PDB::internal::rangelib::detail::generator_functor<res> g_type;
    return generated_range<g_type,stopper>( g_type(g), s);
}


///////////////////////////////////////////////////////////////////
// Stoppers

// returns the first n elements from the range
struct gen_n {
    gen_n( int n) : m_remaining(n) {}

    template<class value_type> bool operator() ( const value_type &) {
        return m_remaining-- > 0;
    }
private:
    int m_remaining;
};


// generates up to a given value (EXCLUDING that value)
template< class value_type>
struct gen_upto_t {
    gen_upto_t( value_type val) : m_val( val) {}
    bool operator()( const value_type & other) const { 
        return other < m_val;
    }
private:
    value_type m_val;
};

template< class value_type> inline gen_upto_t<value_type>
gen_upto( value_type val) {
    return gen_upto_t<value_type>( val);
}

// generates down to a given value (EXCLUDING that value)
template< class value_type>
struct gen_downto_t {
    gen_downto_t( value_type val) : m_val( val) {}
    bool operator()( const value_type & other) const { 
        return m_val < other;
    }
private:
    value_type m_val;
};

template< class value_type> inline gen_downto_t<value_type>
gen_downto( value_type val) {
    return gen_downto_t<value_type>( val);
}


}}}}

#endif
