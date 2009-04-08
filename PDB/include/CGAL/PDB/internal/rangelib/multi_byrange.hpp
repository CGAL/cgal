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

#ifndef CGAL_PDB_BOOST_RTL_MULTI_BYRNG_HPP_INCLUDED
#define CGAL_PDB_BOOST_RTL_MULTI_BYRNG_HPP_INCLUDED

#include <CGAL/PDB/internal/rangelib/priv/defs.hpp>
#include <CGAL/PDB/internal/rangelib/priv/traits.hpp>
#include <boost/iterator.hpp>

#include <CGAL/PDB/internal/rangelib/priv/rng_fun.hpp>
#include <vector>
#include <iterator>
#include <functional>

// rangelib = Range Library
namespace CGAL { namespace PDB { namespace internal { namespace rangelib { 

namespace detail {
	template< class type> struct copy_ptr
	{
	public:
		copy_ptr( const copy_ptr<type> & other) : m_pType( other.copy()) {
		}
		copy_ptr< type> & operator=( const copy_ptr< type> & other) {
			type * copy = m_pType;
			m_pType = 0;
			delete copy;
			m_pType = other.copy();
		}

		copy_ptr( type * pType = 0) : m_pType( pType) {}
		~copy_ptr() { delete m_pType; }

		type * get() const { return m_pType; }
		
	private:
		type * copy() const { return m_pType ? new type(*m_pType) : 0; }
	private:
		type * m_pType;
	};

	
	template< class f> struct multirange_function_finder {
		typedef typename f::result_type result_type;
		typedef typename result_type::value_type value_type;
	};

	// FIXME recalculate ALL iterator tags; for instance, in this case, we can have a forward iterator.
    template< class r, class multirange_function>
    struct multi_byrange_iterator : public ::boost::iterator< 
												::std::input_iterator_tag, 
												typename multirange_function_finder<multirange_function>::value_type > {

        typedef multi_byrange_iterator<r,multirange_function> self_type;

        typedef typename multirange_function_finder<multirange_function>::result_type result_range_type;
		typedef typename result_range_type::value_type result_value_type;

		typedef typename range_finder<r>:: range_type source_range_type;
        typedef typename source_range_type::iterator iterator_type;

        multi_byrange_iterator( iterator_type first, iterator_type last, multirange_function f)
            : m_f( f), 
			  m_first( first), m_last( last),
			  m_result(first != last ? new result_range_type( m_f(*m_first++)) : 0) {

			if ( m_result.get())
				if ( !*m_result.get() )
					// the first element contained an empty range
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
        const result_value_type& operator*() const { return **m_result.get(); }
        const result_value_type* operator->() const { return &**m_result.get(); }

        bool more() const { 
			return m_result.get() && *m_result.get(); 
		}
    private:

        void compute_value() {
			if ( m_result.get() == 0)
				return; // the original range was empty

			if ( *m_result.get())
				++*m_result.get();

			while ( !*m_result.get() && (m_first != m_last)) 
				*m_result.get() = m_f( *m_first++);
        }
    private:
		// IMPORTANT: order counts, see initialization in constructor
        multirange_function m_f;
        iterator_type m_first, m_last;

		copy_ptr<result_range_type> m_result;
    };


    template<class r, class f>
    inline bool operator==( const multi_byrange_iterator<r,f> & first, const multi_byrange_iterator<r,f> & second) {
        return first.more() == second.more();
    }
    template<class r, class f>
    inline bool operator!=( const multi_byrange_iterator<r,f> & first, const multi_byrange_iterator<r,f> & second) {
        return first.more() != second.more();
    }
}


template< class r, class multirange_function>
struct multi_byrange_range : public irange< ::CGAL::PDB::internal::rangelib::detail::multi_byrange_iterator<r,multirange_function> > {
    typedef ::CGAL::PDB::internal::rangelib::detail::multi_byrange_iterator<r,multirange_function> new_iterator;
    typedef irange<new_iterator> base;

    multi_byrange_range( r & rng, multirange_function f)
        : base( new_iterator( rng.begin(), rng.end(), f),
                new_iterator( rng.end(), rng.end(), f) ) {}
};


#ifndef CGAL_PDB_BOOST_RTL_WORKAROUND_VC6
template<class r, class f>
multi_byrange_range<const r,f> multied_byrange( const r & rng, f func) {
    return multi_byrange_range<const r,f>(rng,func);
}

template<class r, class res, class arg>
multi_byrange_range<const r, ::std::pointer_to_unary_function<arg,res> > multied_byrange_f( const r & rng, res (*func)(arg) ) {
	typedef ::std::pointer_to_unary_function<arg,res> ptr;
    return multi_byrange_range<const r, ptr>(rng, ptr(func) );
}

#endif

template<class r, class f>
multi_byrange_range<r,f> multied_byrange( r & rng, f func) {
    return multi_byrange_range<r,f>(rng,func);
}


// stupid VC6 again, can't find the right overload, for a function
template<class r, class res, class arg>
multi_byrange_range<r, ::std::pointer_to_unary_function<arg,res> > multied_byrange_f( r & rng, res (*func)(arg) ) {
	typedef ::std::pointer_to_unary_function<arg,res> ptr;
    return multi_byrange_range<r, ptr>(rng, ptr(func) );
}



}}}}

#endif
