//  (C) Copyright Gennadiy Rozental 2003.
//  Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE_1_0.txt or copy at 
//  http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org/libs/test for the library home page.
//
//  File        : $RCSfile$
//
//  Version     : $Revision$
//
//  Description : implements support for test cases templates instantiated with 
//                sequence of test types
// ***************************************************************************

#ifndef BOOST_TEST_CASE_TEMPLATE_HPP_071894GER
#define BOOST_TEST_CASE_TEMPLATE_HPP_071894GER

// Boost.Test
#include <boost/test/unit_test_suite.hpp>

// Boost
#include "boost/mpl/size.hpp"
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/identity.hpp>

//____________________________________________________________________________//

#define BOOST_META_FUNC_TEST_CASE( the_function )   \
struct meta_ ## the_function {                      \
    template<typename T>                            \
    static void execute( T* = 0 )                   \
    {                                               \
        the_function<T>();                          \
    }                                               \
}                                                   \
/**/

#define BOOST_FUNC_TEMPLATE_TEST_CASE( the_function, typelist ) \
boost::unit_test::create_test_case_template( meta_ ## the_function(), typelist(), #the_function )
    
namespace boost {
namespace unit_test {
namespace ut_detail {

// ************************************************************************** //
// **************          test_case_template_instance         ************** //
// ************************************************************************** //
// Generate test case by supplied test case template and test type

template<typename TestCaseTemplate,typename TestType>
class test_case_template_instance : public test_case {
    typedef TestType*   test_type_ptr;
public:
    explicit            test_case_template_instance( const_string template_name_ )
    : test_case( template_name_, true, 1 )  {}
    
protected:
    // test case implementation
    void                do_run()            { TestCaseTemplate::execute( test_type_ptr() ); }

};

//____________________________________________________________________________//

// ************************************************************************** //
// **************           test_case_instance_runner          ************** //
// ************************************************************************** //
// Instantiate generated test case and run it.

template<typename TestCaseTemplate>
struct test_case_instance_runner {
    explicit            test_case_instance_runner( const_string template_name_ )
    : m_template_name( template_name_ ) {}

    template<typename TestType>
    void                operator()( ::boost::mpl::identity<TestType> )
    {
        test_case_template_instance<TestCaseTemplate,TestType> the_instance( m_template_name ); //!! could this throw?

        the_instance.run();
    }

    const_string    m_template_name;
};

} // namespace ut_detail

// ************************************************************************** //
// **************              test_case_template              ************** //
// ************************************************************************** //

template<typename TestCaseTemplate,typename TestTypesList>
class test_case_template : public test_case {
public:
    // Constructor
    explicit            test_case_template( const_string name_ )
    : test_case( name_, false, 1, false ), m_template_holder( p_name.get() ) {}

    // access methods
    unit_test_counter   size() const    { return ::boost::mpl::size<TestTypesList>::value; }

protected:
    
    // test case implementation
    void                do_run()
    {
        ::boost::mpl::for_each<TestTypesList,mpl::make_identity<boost::mpl::_> >( m_template_holder );
    }

    // Data members
    ut_detail::test_case_instance_runner<TestCaseTemplate> m_template_holder; // need instance to match for_each interface
};

//____________________________________________________________________________//

// ************************************************************************** //
// **************               object generators              ************** //
// ************************************************************************** //

template<typename TestCaseTemplate, typename TestTypesList>
inline test_case*
create_test_case_template( TestCaseTemplate, TestTypesList, std::string name_ )
{
    return new test_case_template<TestCaseTemplate,TestTypesList>( ut_detail::normalize_test_case_name( name_ ) );
}

//____________________________________________________________________________//

} // unit_test
} // namespace boost

// ***************************************************************************
//  Revision History :
//  
//  $Log$
//  Revision 1.1.1.2  2004/11/20 10:52:14  spion
//  Import of Boost v. 1.32.0
//
//  Revision 1.8  2004/07/19 12:14:34  rogeeff
//  guard rename
//
//  Revision 1.7  2004/06/07 07:33:49  rogeeff
//  detail namespace renamed
//
//  Revision 1.6  2004/05/21 06:19:35  rogeeff
//  licence update
//
//  Revision 1.5  2004/05/11 11:00:35  rogeeff
//  basic_cstring introduced and used everywhere
//  class properties reworked
//
//  Revision 1.4  2003/12/01 00:41:56  rogeeff
//  prerelease cleaning
//
// ***************************************************************************

#endif // BOOST_TEST_CASE_TEMPLATE_HPP_071894GER

