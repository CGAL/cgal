//  (C) Copyright Gennadiy Rozental 2001-2003.
//  Use, modification, and distribution are subject to the 
//  Boost Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org/libs/test for the library home page.
//
//  File        : $RCSfile$
//
//  Version     : $Revision$
//
//  Description : provides extention for unit test framework that allows usage
//  boost::function as a test case base function.
// ***************************************************************************

#ifndef BOOST_UNIT_TEST_SUITE_EX_HPP
#define BOOST_UNIT_TEST_SUITE_EX_HPP

// Boost.Test
#include <boost/test/unit_test_suite.hpp>
#include <boost/test/detail/unit_test_config.hpp>

// BOOST
#include <boost/function/function0.hpp>
#include <boost/function/function1.hpp>

// STL
#include <string>  // for std::string

namespace boost {

namespace unit_test_framework {

// ************************************************************************** //
// **************            boost_function_test_case             ************** //
// ************************************************************************** //

class boost_function_test_case : public test_case {
public:
    typedef function0<void> function_type;

    // Constructor
    boost_function_test_case( function_type f_, std::string const& name_ )
    : test_case( name_, true, 1 ), m_function( f_ ) {}

protected:
    // test case implementation
    void                do_run()        { m_function(); }

private:
    // Data members
    function_type       m_function;
};

// ************************************************************************** //
// **************        parametrized_function_test_case       ************** //
// ************************************************************************** //

template <typename ParamIterator, typename ParameterType>
class parametrized_boost_function_test_case : public test_case {
public:
    typedef function1<void,ParameterType> function_type;

    // Constructor
    parametrized_boost_function_test_case( function_type f_, std::string const& name_,
                                        ParamIterator const& begin_, ParamIterator const& end_ )
    : test_case( name_, true, 0 ), m_first_parameter( begin_ ), m_last_parameter( end_ ), m_function( f_ )
    {
        p_stages_amount.set( detail::distance( begin_, end_ ) );
    }

    // test case implementation
    void                do_init()       { m_curr_parameter = m_first_parameter; }
    void                do_run()        { m_function( *m_curr_parameter ); ++m_curr_parameter; }

private:
    // Data members
    ParamIterator       m_first_parameter;
    ParamIterator       m_last_parameter;
    ParamIterator       m_curr_parameter;

    function_type       m_function;
};

// ************************************************************************** //
// **************               object generators              ************** //
// ************************************************************************** //

inline test_case*
create_test_case( function0<void> const& fct_, std::string name_ )
{
    return new boost_function_test_case( fct_, detail::normalize_test_case_name( name_ ) );
}

template<typename ParamIterator, typename ParameterType>
inline test_case*
create_test_case( function1<void,ParameterType> const& fct_, std::string name_, 
                  ParamIterator const& begin_, ParamIterator const& end_ )
{
    return new parametrized_boost_function_test_case<ParamIterator,ParameterType>(
                    fct_, detail::normalize_test_case_name( name_ ), begin_, end_ );
}

} // unit_test_framework

} // namespace boost

// ***************************************************************************
//  Revision History :
//  
//  $Log$
//  Revision 1.1  2004/05/23 10:51:38  spion
//  Initial revision
//
//  Revision 1.15  2003/12/01 00:41:56  rogeeff
//  prerelease cleaning
//

// ***************************************************************************

#endif // BOOST_UNIT_TEST_SUITE_EX_HPP
