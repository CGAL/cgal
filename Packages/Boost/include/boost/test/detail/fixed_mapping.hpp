//  (C) Copyright Gennadiy Rozental 2001-2004.
//  Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE_1_0.txt or copy at 
//  http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org/libs/test for the library home page.
//
//  File        : $RCSfile$
//
//  Version     : $Revision$
//
//  Description : 
// ***************************************************************************

#ifndef BOOST_TEST_FIXED_MAPPING_HPP_071894GER
#define BOOST_TEST_FIXED_MAPPING_HPP_071894GER

// Boost
#include <boost/preprocessor/repetition/repeat.hpp>
#include <boost/preprocessor/arithmetic/add.hpp>
#include <boost/call_traits.hpp>
#include <boost/detail/binary_search.hpp>

// STL
#include <vector>
#include <functional>
#include <algorithm>
#include <utility>

namespace boost {

namespace unit_test {

// configurable maximum fixed sized mapping size supported by this header.
// You could redefine it before inclusion of this file.
#ifndef MAX_MAP_SIZE
#define MAX_MAP_SIZE 15
#endif

#define CONSTR_DECL_MID( z, i, dummy1 ) key_param_type key##i, value_param_type v##i,
#define CONSTR_BODY_MID( z, i, dummy1 ) add_pair( key##i, v##i );

#define CONSTR_DECL( z, n, dummy1 )                                 \
    fixed_mapping( BOOST_PP_REPEAT_ ## z( n, CONSTR_DECL_MID, "" )  \
                         value_param_type invalid_value )           \
    : m_invalid_value( invalid_value )                              \
    {                                                               \
        BOOST_PP_REPEAT_ ## z( n, CONSTR_BODY_MID, "" )             \
        init();                                                     \
    }                                                               \
/**/

#define CONTRUCTORS( n ) BOOST_PP_REPEAT( n, CONSTR_DECL, "" )

template<typename Key, typename Value, typename Compare = std::less<Key> >
class fixed_mapping
{
    typedef std::pair<Key,Value>                            elem_type;
    typedef std::vector<elem_type>                          map_type;
    typedef typename std::vector<elem_type>::const_iterator iterator;

    typedef typename call_traits<Key>::param_type           key_param_type;
    typedef typename call_traits<Value>::param_type         value_param_type;
    typedef typename call_traits<Value>::const_reference    value_ref_type;

    // bind( Compare(), bind(select1st<elem_type>(), _1),  bind(identity<Key>(), _2) )
    struct p1 : public std::binary_function<elem_type,Key,bool>
    {
        bool operator()( elem_type const& x, Key const& y ) const { return Compare()( x.first, y ); }
    };

    // bind( Compare(), bind(select1st<elem_type>(), _1), bind(select1st<elem_type>(), _2) )
    struct p2 : public std::binary_function<elem_type,elem_type,bool>
    {
        bool operator()( elem_type const& x, elem_type const& y ) const { return Compare()( x.first, y.first ); }
    };

public:
    // Constructors
    CONTRUCTORS( BOOST_PP_ADD( MAX_MAP_SIZE, 1 ) )

    // key -> value access
    value_ref_type  operator[]( key_param_type key ) const
    {
        iterator it = boost::detail::lower_bound( m_map.begin(), m_map.end(), key, p1() );

        return (it == m_map.end() || Compare()( key, it->first ) ) ? m_invalid_value : it->second;
    }

private:
    // Implementation
    void            init()                                                  { std::sort( m_map.begin(), m_map.end(), p2() ); }
    void            add_pair( key_param_type key, value_param_type value )  { m_map.push_back( elem_type( key, value ) ); }

    // Data members
    Value           m_invalid_value;
    map_type        m_map;
};

} // namespace unit_test

} // namespace boost

#undef MAX_MAP_SIZE
#undef CONSTR_DECL_MID
#undef CONSTR_BODY_MID
#undef CONSTR_DECL
#undef CONTRUCTORS

// ***************************************************************************
//  Revision History :
//  
//  $Log$
//  Revision 1.1  2004/11/20 10:52:19  spion
//  Initial revision
//
//  Revision 1.4  2004/07/26 05:35:53  david_abrahams
//  Use boost::detail::lower_bound, which is actually guaranteed to work
//  with heterogeneous comparisons.  C++98 binary searches are not, and
//  VC++ 8.0 beta rejects the code using std::lower_bound.
//
//  Revision 1.3  2004/07/19 12:21:08  rogeeff
//  guard rename
//
//  Revision 1.2  2004/05/21 06:19:35  rogeeff
//  licence update
//
//  Revision 1.1  2004/05/13 09:06:48  rogeeff
//  added fixed_mapping
//
// ***************************************************************************

#endif // BOOST_TEST_FIXED_MAPPING_HPP_071894GER

