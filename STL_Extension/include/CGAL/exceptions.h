// Copyright (c) 2006  INRIA Sophia-Antipolis (France) and
// Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Lutz Kettner, Sylvain Pion

#ifndef CGAL_EXCEPTIONS_H
#define CGAL_EXCEPTIONS_H

#include <CGAL/config.h>
#include <stdexcept>
#include <string>

// Address the warning C4003: not enough actual parameters for macro 'BOOST_PP_SEQ_DETAIL_IS_NOT_EMPTY'
// lexical_cast.hpp includes files from boost/preprocessor
// This concerns boost 1_67_0
#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable: 4003)
#endif
#include <boost/lexical_cast.hpp>
#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif


namespace CGAL {

// [Sylvain] This was originaly written in the Exacus library.
// I kept most doxygen comments.


/*! \defgroup CGAL_assert Assertions
    \brief <tt>\#include <CGAL/basic.h></tt> for pre- and postconditions,
    assertions, warnings, and error handler.

    \c CGAL/basic.h provides macros for pre- and postconditions, assertions,
    warnings, and errors that are active by default. There are variants
    for expensive checks that are inactive by default.

    A failed precondition, postcondition or assertion is a fatal error
    and the library aborts by default the program execution. For warnings
    the program continues to run. Both can be changed with the respective
    \c error_behavior or \c warning_behavior, or with the corresponding
    \c error_handler or \c warning_handler.

    However, if the checks are disabled with corresponding
    preprocessor symbols or \c NDEBUG, none of this error handling
    will happen and the program execution continues with potentially
    disastrous results. So, these checks are not for normal error
    handling, for example testing proper user input or file
    formats. Furthermore, it is recommended to keep the checks also in
    the production version of the program, keeping in mind that they
    should not cost too much runtime overhead, maybe 10%.

    \see \c CGAL_test for checks that cannot be disabled and
         \c CGAL_error_msg for an unconditional error message that cannot
         be disabled with macros.
*/

// Failure types and exceptions
// ============================

//! Exception base class for all failure types of assertions etc.
class Failure_exception : public std::logic_error {
    std::string      m_lib;
    std::string      m_expr; // can be empty
    std::string      m_file;
    int              m_line;
    std::string      m_msg;  // can be empty
public:
    //! initializes local members and the <tt>std::logic_error</tt> with
    //! a suitable message.
    Failure_exception( std::string   lib,
                       std::string   expr,
                       std::string   file,
                       int           line,
                       std::string   msg,
                       std::string   kind = "Unknown kind") :
    std::logic_error( lib + std::string( " ERROR: ") + kind + std::string( "!")
        + ((expr.empty()) ? (std::string("")) : (std::string("\nExpr: ")+expr))
        + std::string( "\nFile: ") + file
        + std::string( "\nLine: ") + boost::lexical_cast<std::string>(line)
        + ((msg.empty()) ? (std::string(""))
                         : (std::string("\nExplanation: ") + msg))),
    m_lib( lib),
    m_expr( expr),
    m_file( file),
    m_line( line),
    m_msg( msg)
    {}

    //! the name of the library that issues this message.
    std::string      library() const { return m_lib; }

    //! expression that failed in assertion, pre-, or postcondition.
    std::string      expression() const { return m_expr; }

    //! source code filename where the failure was detected.
    std::string      filename() const { return m_file; }

    //! line number in source code file where the failure was detected.
    int              line_number() const { return m_line; }

    //! an optional message explaining the kind of failure.
    std::string      message() const { return m_msg; }
};

//! Exception thrown for \c CGAL_error_msg.
class Error_exception : public Failure_exception {
public:
    Error_exception( std::string lib,
                     std::string msg,
                     std::string file,
                     int line)
        : Failure_exception( lib, "", file, line, msg, "failure" ) {}
};

//! Exception thrown for \c CGAL_precond.
class Precondition_exception : public Failure_exception {
public:
    Precondition_exception( std::string lib,
                       std::string expr,
                       std::string file,
                       int line,
                       std::string msg)
        : Failure_exception( lib, expr, file, line, msg,
                             "precondition violation") {}
};

//! Exception thrown for \c CGAL_postcond.
class Postcondition_exception : public Failure_exception {
public:
    Postcondition_exception( std::string lib,
                        std::string expr,
                        std::string file,
                        int line,
                        std::string msg)
        : Failure_exception( lib, expr, file, line, msg,
                             "postcondition violation") {}
};

//! Exception thrown for \c CGAL_assert.
class Assertion_exception : public Failure_exception {
public:
    Assertion_exception( std::string lib,
                      std::string expr,
                      std::string file,
                      int line,
                      std::string msg)
        : Failure_exception( lib, expr, file, line, msg,
                             "assertion violation") {}
};

//! Exception thrown for \c CGAL_test.
class Test_exception : public Failure_exception {
public:
    Test_exception( std::string lib,
                    std::string expr,
                    std::string file,
                    int line,
                    std::string msg)
        : Failure_exception( lib, expr, file, line, msg,
                             "test in test-suite violation") {}
};

//! Exception thrown for \c CGAL_warning.
class Warning_exception : public Failure_exception {
public:
    Warning_exception( std::string lib,
                       std::string expr,
                       std::string file,
                       int line,
                       std::string msg)
        : Failure_exception( lib, expr, file, line, msg,
                             "warning condition failed") {}
};

namespace internal {

// The following classes are useful to create output iterators (with the help
// of boost::function_output_iterator) that will throw as soon as something is being written.
class Throw_at_output_exception : public std::exception { };

struct Throw_at_output
{
  template<class T>
  void operator()(const T& /* t */) const { throw Throw_at_output_exception(); }
};

} // namespace internal
} //namespace CGAL

#endif // CGAL_EXCEPTIONS_H
