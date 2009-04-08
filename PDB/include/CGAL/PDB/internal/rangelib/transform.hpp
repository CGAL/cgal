
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


#ifndef CGAL_PDB_BOOST_RTL_TRANSFORM_HPP_INCLUDED
#define CGAL_PDB_BOOST_RTL_TRANSFORM_HPP_INCLUDED

#include <CGAL/PDB/internal/rangelib/priv/defs.hpp>
#include <CGAL/PDB/internal/rangelib/priv/traits.hpp>
#include <CGAL/PDB/internal/rangelib/range.hpp>

#include <boost/iterator/transform_iterator.hpp>
#include <functional>


// rangelib = Range Library
namespace CGAL { namespace PDB { namespace internal { namespace rangelib {

namespace detail {

    template< class r, class function> 
    struct transform_iterator_finder {
        typedef typename ::CGAL::PDB::internal::rangelib::range_finder<r>::range_type r_type;
        typedef typename r_type::iterator i_type;

        typedef ::boost::transform_iterator<function, i_type> transform_type;
    };
}


template< class r, class function > 
struct transformed_range : public irange< typename ::CGAL::PDB::internal::rangelib::detail::transform_iterator_finder<r,function>::transform_type > {
    typedef irange< typename ::CGAL::PDB::internal::rangelib::detail::transform_iterator_finder<r,function>::transform_type> base;
    typedef typename ::CGAL::PDB::internal::rangelib::detail::transform_iterator_finder<r,function>::i_type old_iterator;
    typedef typename ::CGAL::PDB::internal::rangelib::detail::transform_iterator_finder<r,function>::transform_type new_iterator;

//    transformed_range( old_iterator first, old_iterator last, function f)
  //      : base( new_iterator(first,f), new_iterator(last,f) ) {
    //}

  typedef typename base::iterator const_iterator;

    transformed_range( r & rng, function f ) 
        : base( new_iterator( rng.begin(),f ), 
                new_iterator( rng.end(),f ) ) {
    }
};


#ifndef CGAL_PDB_BOOST_RTL_WORKAROUND_VC6
template< class r, class function> inline transformed_range<const r, function> 
transformed( const r & rng, function f) {
//    return transformed_range<const r, function>( rng.begin(), rng.end(), f);
    return transformed_range<const r, function>( rng, f);
}
template< class r, class arg, class res> inline transformed_range<const r, ::std::pointer_to_unary_function<arg,res> > 
transformed_f( const r & rng, res (*f)(arg) ) {
//    return transformed_range<const r, ::std::pointer_to_unary_function<arg,res> >( rng.begin(), rng.end(), 
  //      ::std::pointer_to_unary_function<arg,res>(f) );
    return transformed_range<const r, ::std::pointer_to_unary_function<arg,res> >( rng,
        ::std::pointer_to_unary_function<arg,res>(f) );
}
#endif


template< class r, class function> inline transformed_range<r, function> 
transformed( r & rng, function f) {
//    return transformed_range<r, function>( rng.begin(), rng.end(), f);
    return transformed_range<r, function>( rng, f);
}

// VC6 chokes if both these functions have the same name
template< class r, class arg, class res> inline transformed_range<r, ::std::pointer_to_unary_function<arg,res> > 
transformed_f( r & rng, res (*f)(arg) ) {
//    return transformed_range<r, ::std::pointer_to_unary_function<arg,res> >( rng.begin(), rng.end(), 
  //      ::std::pointer_to_unary_function<arg,res>(f) );
    return transformed_range<r, ::std::pointer_to_unary_function<arg,res> >( rng,
        ::std::pointer_to_unary_function<arg,res>(f) );
}




// FIXME explain usage
template< class r_and_f_finder> inline transformed_range<typename r_and_f_finder::range_type, typename r_and_f_finder::function_type> 
transformed( const r_and_f_finder & val) {
    return transformed( val.range(), val.function() );
}



// FIXME maybe select1st/select2nd from MPL can be used
template<class pair_type>
struct select_2nd {
    BOOST_STATIC_CONSTANT( bool, is_const = ::boost::is_const<pair_type>::value);
    typedef typename pair_type::second_type result_type_raw;
    typedef typename ::boost::mpl::if_c<is_const, const result_type_raw, result_type_raw>::type result_type;

    result_type & operator() ( pair_type & val) const { return val.second; }
};

template<class pair_type>
struct select_1st {
    BOOST_STATIC_CONSTANT( bool, is_const = ::boost::is_const<pair_type>::value);
    typedef typename pair_type::first_type result_type_raw;
    typedef typename ::boost::mpl::if_c<is_const, const result_type_raw, result_type_raw>::type result_type;

    result_type & operator() ( pair_type & val) const { return val.first; }
};


namespace detail {
    template< class r_type, class transformer_type>
    struct transform_helper {
        typedef r_type range_type;
        // this is so that, in case you're using the wrong (transformer) type,
        // to get a compile-time error ASAP
        typedef typename transformer_type::result_type result_type;
        typedef transformer_type function_type;
        range_type & range() const { return m_r; }
        function_type function() const { return m_f; }

        transform_helper( range_type & r, function_type f) : m_r(r), m_f(f) {}
    private:
        mutable range_type & m_r;
        function_type m_f;
    };

    template< class collection> struct transform_key_finder {
        BOOST_STATIC_CONSTANT( bool, is_const = ::boost::is_const<collection>::value);
        typedef typename collection::value_type value_type_raw;
        typedef typename ::boost::mpl::if_c<is_const,const value_type_raw, value_type_raw>::type value_type;

        typedef select_1st<value_type> key_finder;
        typedef transform_helper<collection,key_finder> type;
    };

    template< class collection> struct transform_value_finder {
        BOOST_STATIC_CONSTANT( bool, is_const = ::boost::is_const<collection>::value);
        typedef typename collection::value_type value_type_raw;
        typedef typename ::boost::mpl::if_c<is_const,const value_type_raw, value_type_raw>::type value_type;

        typedef select_2nd<value_type> value_finder;
        typedef transform_helper<collection,value_finder> type;
    };
}


// FIXME explain this in the docs ;)

/*
    given a collection (or a sequence of pairs),
    these allow you to easily grab either the first (the key) or the second,
    and return a transformed range.
    (a range of keys or a range of values)
*/
#ifndef CGAL_PDB_BOOST_RTL_WORKAROUND_VC6
template<class collection> inline typename ::CGAL::PDB::internal::rangelib::detail::transform_key_finder<const collection>::type
pair_1st( const collection & c) {
    typedef typename ::CGAL::PDB::internal::rangelib::detail::transform_key_finder<const collection> finder;
    typedef typename finder::type result_type;
    typedef typename finder::key_finder key_finder;
    return result_type(c, key_finder() );
}
template<class collection> inline typename ::CGAL::PDB::internal::rangelib::detail::transform_value_finder<const collection>::type
pair_2nd( const collection & c) {
    typedef typename ::CGAL::PDB::internal::rangelib::detail::transform_value_finder<const collection> finder;
    typedef typename finder::type result_type;
    typedef typename finder::value_finder value_finder;
    return result_type(c, value_finder() );
}
#endif

template<class collection> inline typename ::CGAL::PDB::internal::rangelib::detail::transform_key_finder<collection>::type
pair_1st( collection & c) {
    typedef typename ::CGAL::PDB::internal::rangelib::detail::transform_key_finder<collection> finder;
    typedef typename finder::type result_type;
    typedef typename finder::key_finder key_finder;
    return result_type(c, key_finder() );
}
template<class collection> inline typename ::CGAL::PDB::internal::rangelib::detail::transform_value_finder<collection>::type
pair_2nd( collection & c) {
    typedef typename ::CGAL::PDB::internal::rangelib::detail::transform_value_finder<collection> finder;
    typedef typename finder::type result_type;
    typedef typename finder::value_finder value_finder;
    return result_type(c, value_finder() );
}

// FIXME - make it possible to transform using a member

}}}}


#endif
