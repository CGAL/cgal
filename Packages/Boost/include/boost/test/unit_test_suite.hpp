//  (C) Copyright Gennadiy Rozental 2001-2003.
//  (C) Copyright Ullrich Koethe 2001.
//  Use, modification, and distribution are subject to the 
//  Boost Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org/libs/test for the library home page.
//
//  File        : $RCSfile$
//
//  Version     : $Revision$
//
//  Description : defines all classes in test_case hierarchy and object generators
//  for them.
// ***************************************************************************

#ifndef BOOST_UNIT_TEST_SUITE_HPP
#define BOOST_UNIT_TEST_SUITE_HPP

// Boost.Test
#include <boost/test/detail/unit_test_monitor.hpp>
#include <boost/test/detail/unit_test_config.hpp>
#include <boost/test/detail/class_properties.hpp>

// BOOST
#include <boost/shared_ptr.hpp>

// STL
#include <string>  // for std::string

//____________________________________________________________________________//

#define BOOST_TEST_CASE( function ) \
boost::unit_test_framework::create_test_case((function), #function )
#define BOOST_CLASS_TEST_CASE( function, tc_instance ) \
boost::unit_test_framework::create_test_case((function), #function, tc_instance )
#define BOOST_PARAM_TEST_CASE( function, begin, end ) \
boost::unit_test_framework::create_test_case((function), #function, (begin), (end) )
#define BOOST_PARAM_CLASS_TEST_CASE( function, tc_instance, begin, end ) \
boost::unit_test_framework::create_test_case((function), #function, tc_instance, (begin), (end) )
#define BOOST_TEST_SUITE( testsuite_name ) \
( new boost::unit_test_framework::test_suite( testsuite_name ) )

namespace boost {

namespace unit_test_framework {

// ************************************************************************** //
// **************                   test_case                  ************** //
// ************************************************************************** //

class test_case {
public:
    typedef detail::unit_test_monitor::error_level error_level_type;

    // post creation configuration
    void                depends_on( test_case const* rhs );

    // Destructor
    virtual             ~test_case()    {}

    // total number of test cases
    virtual unit_test_counter size() const;

    // execute this method to run the test case
    void                run();

    // status
    bool                has_passed() const;

    // public properties
    BOOST_READONLY_PROPERTY( int, 2, (test_case,test_suite) )
                        p_timeout;                  // timeout for the excecution monitor
    BOOST_READONLY_PROPERTY( unit_test_counter, 1, (test_suite) )
                        p_expected_failures;        // number of assertions that are expected to fail in this test case
    BOOST_READONLY_PROPERTY( bool, 0, () )
                        p_type;                     // true = test case, false - test suite
    BOOST_READONLY_PROPERTY( std::string, 0, () )
                        p_name;                     // name for this test case


protected:
    // protected properties
    BOOST_READONLY_PROPERTY( bool, 2, (test_case,test_suite) )
                        p_compound_stage;           // used to properly manage progress report
    BOOST_READWRITE_PROPERTY( unit_test_counter )
                        p_stages_amount;            // number of stages this test consist of; stage could be another test case
                                                    // like with test_suite, another parameterized test for parameterized_test_case
                                                    // or 1 stage that reflect single test_case behaviour

    // access methods
    void                curr_stage_is_compound();

    // Constructor
    explicit            test_case( std::string const&   name_,
                                   bool                 type_,
                                   unit_test_counter    stages_amount_,
                                   bool                 monitor_run_    = true );

    // test case implementation hooks to be called with unit_test_monitor or alone
    virtual void        do_init()       {}
    virtual void        do_run()        {}
    virtual void        do_destroy()    {}

private:
    // Data members
    struct Impl;
    boost::shared_ptr<Impl> m_pimpl;
};

//____________________________________________________________________________//

extern detail::unit_test_monitor the_monitor;

template<typename Exception, typename ExceptionTranslator>
void
register_exception_translator( ExceptionTranslator const& tr, boost::type<Exception>* d = 0 )
{
    the_monitor.register_exception_translator( tr, d );
}

//____________________________________________________________________________//

// ************************************************************************** //
// **************              function_test_case              ************** //
// ************************************************************************** //

class function_test_case : public test_case {
public:
    typedef void  (*function_type)();

    // Constructor
    function_test_case( function_type f_, std::string const& name_ )
    : test_case( name_, true, 1 ), m_function( f_ ) {}

protected:
    // test case implementation
    void                do_run()        { m_function(); }

private:
    // Data members
    function_type       m_function;
};

// ************************************************************************** //
// **************                class_test_case               ************** //
// ************************************************************************** //

template<class UserTestCase>
class class_test_case : public test_case {
public:
    typedef void  (UserTestCase::*function_type)();

    // Constructor
    class_test_case( function_type f_, std::string const& name_, boost::shared_ptr<UserTestCase> const& user_test_case_ )
    : test_case( name_, true, 1 ), m_user_test_case( user_test_case_ ), m_function( f_ ) 
    {}

private:
    // test case implementation
    void                do_run()
    { 
        if( (!!m_user_test_case) && m_function ) 
            ((*m_user_test_case).*m_function)();
    }
    void                do_destroy()
    { 
        m_user_test_case.reset(); // t ofree the reference to the shared use test case instance
    }

    // Data members
    boost::shared_ptr<UserTestCase> m_user_test_case;
    function_type       m_function;
};

// ************************************************************************** //
// **************        parametrized_function_test_case       ************** //
// ************************************************************************** //

template <typename ParamIterator, typename ParameterType>
class parametrized_function_test_case : public test_case {
public:
    typedef void  (*function_type)( ParameterType );

    // Constructor
    parametrized_function_test_case( function_type f_, std::string const& name_,
                                     ParamIterator const& par_begin_, ParamIterator const& par_end_ )
    : test_case( name_, true, 0 ), m_first_parameter( par_begin_ ), m_last_parameter( par_end_ ), m_function( f_ )
    {
       // the typecasts are here to keep Borland C++ Builder 5 happy, for other compilers they have no effect:
       p_stages_amount.set( detail::distance( (ParamIterator)par_begin_, (ParamIterator)par_end_ ) );
    }

    // test case implementation
    void                do_init()       { m_curr_parameter = m_first_parameter; }
    void                do_run()        { m_function( *m_curr_parameter++ ); }

private:
    // Data members
    ParamIterator       m_first_parameter;
    ParamIterator       m_last_parameter;
    ParamIterator       m_curr_parameter;

    function_type       m_function;
};

// ************************************************************************** //
// **************         parametrized_class_test_case         ************** //
// ************************************************************************** //

template <class UserTestCase, class ParamIterator, typename ParameterType>
class parametrized_class_test_case : public test_case {
public:
    typedef void  (UserTestCase::*function_type)( ParameterType );

    // Constructor
    parametrized_class_test_case( function_type f_, std::string const& name_, boost::shared_ptr<UserTestCase>const & user_test_case_,
                                  ParamIterator const& par_begin_, ParamIterator const& par_end_ )
    : test_case( name_, true, 0 ), m_first_parameter( par_begin_ ), m_last_parameter( par_end_ ),
      m_user_test_case( user_test_case_ ), m_function( f_ )
    {
       // the typecasts are here to keep Borland C++ Builder 5 happy, for other compilers they have no effect:
       p_stages_amount.set( detail::distance( (ParamIterator)par_begin_, (ParamIterator)par_end_ ) );
    }

    // test case implementation
    void                do_init()       { m_curr_parameter = m_first_parameter; }
    void                do_run()        { ((*m_user_test_case).*m_function)( *m_curr_parameter++ ); }
    void                do_destroy()    { m_user_test_case.reset(); }

private:
    // Data members
    ParamIterator       m_first_parameter;
    ParamIterator       m_last_parameter;
    ParamIterator       m_curr_parameter;

    boost::shared_ptr<UserTestCase> m_user_test_case;
    function_type       m_function;
};

// ************************************************************************** //
// **************                  test_suite                  ************** //
// ************************************************************************** //

class test_suite : public test_case {
public:
    // Constructor
    explicit test_suite( std::string const& name_ = "Master" );

    // Destructor
    virtual             ~test_suite();

    // test case list management
    void                add( test_case* tc_, unit_test_counter expected_failures_ = 0, int timeout_ = 0 );

    // access methods
    unit_test_counter   size() const;

    // test case implementation
    void                do_init();
    void                do_run();

private:
    // Data members
    struct Impl;
    boost::shared_ptr<Impl> m_pimpl;
};

// ************************************************************************** //
// **************               object generators              ************** //
// ************************************************************************** //

namespace detail {

std::string const& normalize_test_case_name( std::string& name_ );

} // namespace detail

//____________________________________________________________________________//

inline test_case*
create_test_case( void (*fct_)(), std::string name_ )
{
    return new function_test_case( fct_, detail::normalize_test_case_name( name_ ) );
}

//____________________________________________________________________________//

template<class UserTestCase>
inline test_case*
create_test_case( void (UserTestCase::*fct_)(), std::string name_, boost::shared_ptr<UserTestCase> const& user_test_case_ )
{
    return new class_test_case<UserTestCase>( fct_, detail::normalize_test_case_name( name_ ), user_test_case_ );
}

//____________________________________________________________________________//

template<typename ParamIterator, typename ParamType>
inline test_case*
create_test_case( void (*fct_)( ParamType ), std::string name_, ParamIterator const& par_begin_, ParamIterator const& par_end_ )
{
    return new parametrized_function_test_case<ParamIterator,ParamType>(
        fct_, detail::normalize_test_case_name( name_ ), par_begin_, par_end_ );
}

//____________________________________________________________________________//

template<class UserTestCase, typename ParamIterator, typename ParamType>
inline test_case*
create_test_case( void (UserTestCase::*fct_)( ParamType ), std::string name_, boost::shared_ptr<UserTestCase> const& user_test_case_,
                  ParamIterator const& par_begin_, ParamIterator const& par_end_ )
{
    return new parametrized_class_test_case<UserTestCase,ParamIterator,ParamType>(
        fct_, detail::normalize_test_case_name( name_ ), user_test_case_, par_begin_, par_end_ );
}

//____________________________________________________________________________//

} // unit_test_framework

} // namespace boost

// ***************************************************************************
//  Revision History :
//  
//  $Log$
//  Revision 1.1  2004/05/23 10:51:38  spion
//  Initial revision
//
//  Revision 1.18  2003/12/01 00:41:56  rogeeff
//  prerelease cleaning
//

// ***************************************************************************

#endif // BOOST_UNIT_TEST_SUITE_HPP

