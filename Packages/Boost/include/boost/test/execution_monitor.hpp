//  (C) Copyright Gennadiy Rozental 2001-2003.
//  (C) Copyright Beman Dawes 2001.
//  Use, modification, and distribution are subject to the 
//  Boost Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org/libs/test for the library home page.
//
//  File        : $RCSfile$
//
//  Version     : $Revision$
//
//  Description : defines abstract monitor interfaces and implements execution excepiton
//  The original Boost Test Library included an implementation detail function
//  named catch_exceptions() which caught otherwise uncaught C++ exceptions.
//  It was derived from an existing test framework by Beman Dawes.  The
//  intent was to expand later to catch other detectable but platform dependent
//  error events like Unix signals or Windows structured C exceptions.
//
//  Requests from early adopters of the Boost Test Library included
//  configurable levels of error message detail, elimination of templates,
//  separation of error reporting, and making the catch_exceptions() facilities
//  available as a public interface.  Support for unit testing also stretched
//  the function based design.  Implementation within the header became less
//  attractive due to the need to include many huge system dependent headers,
//  although still preferable in certain cases.
//
//  All those issues have been addressed by introducing the class-based
//  design presented here.
// ***************************************************************************

#ifndef BOOST_EXECUTION_MONITOR_HPP
#define BOOST_EXECUTION_MONITOR_HPP

// BOOST TEST
#include <boost/test/detail/unit_test_config.hpp>

// BOOST
#include <boost/scoped_ptr.hpp>
#include <boost/type.hpp>

namespace boost {

class execution_monitor;

namespace detail {

// ************************************************************************** //
// **************       detail::translate_exception_base       ************** //
// ************************************************************************** //

class translate_exception_base {
public:
    // Constructor
    explicit    translate_exception_base( boost::scoped_ptr<translate_exception_base>& next )
    {
        next.swap( m_next );
    }

    // Destructor
    virtual     ~translate_exception_base() {}

    virtual int operator()( boost::execution_monitor& mon ) = 0;

protected:
    // Data members
    boost::scoped_ptr<translate_exception_base> m_next;
};

}

// ************************************************************************** //
// **************              execution_exception             ************** //
// ************************************************************************** //
    
//  design rationale: fear of being out (or nearly out) of memory.
    
class execution_exception {
    typedef unit_test_framework::c_string_literal c_string_literal;
public:
    enum error_code {
        //  These values are sometimes used as program return codes.
        //  The particular values have been choosen to avoid conflicts with
        //  commonly used program return codes: values < 100 are often user
        //  assigned, values > 255 are sometimes used to report system errors.
        //  Gaps in values allow for orderly expansion.
        
        no_error               = 0,   // for completeness only; never returned
        user_error             = 200, // user reported non-fatal error
        cpp_exception_error    = 205, // see note (1) below
        system_error           = 210, // see note (2) below
        timeout_error          = 215, // only detectable on certain platforms
        user_fatal_error       = 220, // user reported fatal error
        system_fatal_error     = 225  // see note (2) below
        
        //  Note 1: Only uncaught C++ exceptions are treated as errors.
        //  If the application catches a C++ exception, it will never reach
        //  the execution_monitor.
        
        //  Note 2: These errors include Unix signals and Windows structured
        //  exceptions.  They are often initiated by hardware traps.
        //
        //  The implementation decides what's a fatal_system_exception and what's
        //  just a system_exception.  Fatal errors are so likely to have corrupted
        //  machine state (like a stack overflow or addressing exception) that it
        //  is unreasonable to continue execution.
    };
    
    // Constructor
    execution_exception( error_code ec_, c_string_literal what_msg_ ) // max length 256 inc '\0'
    : m_error_code( ec_ ), m_what( what_msg_ ) {}

    // access methods
    error_code          code() const { return m_error_code; }
    c_string_literal    what() const { return m_what; }

private:
    // Data members
    error_code          m_error_code;
    c_string_literal    m_what;
}; // execution_exception

// ************************************************************************** //
// **************               execution_monitor              ************** //
// ************************************************************************** //

class execution_monitor {
public:
    // Destructor
    virtual ~execution_monitor()    {}
    
    int execute( bool catch_system_errors = true, int timeout_ = 0 );  // timeout is in seconds
    //  The catch_system_errors parameter specifies whether the monitor should 
    //  try to catch system errors/exeptions that would cause program to crash 
    //  in regular case
    //  The timeout argument specifies the seconds that elapse before
    //  a timer_error occurs.  May be ignored on some platforms.
    //
    //  Returns:  Value returned by function().
    //
    //  Effects:  Calls run_function() inside a try/catch block which also may
    //  include other unspecified platform dependent error detection code.
    //
    //  Throws: execution_exception on an uncaught C++ exception,
    //  a hardware or software signal, trap, or other exception.
    //
    //  Note: execute() doesn't consider it an error for function() to
    //  return a non-zero value.
    
    virtual int function() = 0;
    //  user supplied function called by run_function()
    
    int         run_function();
    // call function() and translate user exceptions with translators registered

    template<typename Exception, typename ExceptionTranslator>
    void        register_exception_translator( ExceptionTranslator const& tr, boost::type<Exception>* = 0 );

private:
    // Data members
    boost::scoped_ptr<detail::translate_exception_base> m_custom_translators;

}; // execution_monitor

namespace detail {

// ************************************************************************** //
// **************         detail::translate_exception          ************** //
// ************************************************************************** //

template<typename Exception, typename ExceptionTranslator>
class translate_exception : public translate_exception_base
{
    typedef boost::scoped_ptr<translate_exception_base> base_ptr;
public:
    explicit    translate_exception( ExceptionTranslator const& tr, base_ptr& next )
    : translate_exception_base( next ), m_translator( tr ) {}

    virtual int operator()( boost::execution_monitor& mon )
    {
        try {
            return m_next ? (*m_next)( mon ) : mon.function();
        }
        catch( Exception const& e )
        {
            m_translator( e );
            return true;
        }
    }

private:
    // Data members
    ExceptionTranslator m_translator;
};

} // namespace detail

template<typename Exception, typename ExceptionTranslator>
void
execution_monitor::register_exception_translator( ExceptionTranslator const& tr, boost::type<Exception>* )
{
    m_custom_translators.reset( 
        new detail::translate_exception<Exception,ExceptionTranslator>( tr,m_custom_translators ) );
}

}  // namespace boost

// ***************************************************************************
//  Revision History :
//  
//  $Log$
//  Revision 1.1  2004/05/23 10:51:33  spion
//  Initial revision
//
//  Revision 1.14  2003/12/01 00:41:56  rogeeff
//  prerelease cleaning
//

// ***************************************************************************

#endif

