/*=============================================================================
    Boost.Wave: A Standard compliant C++ preprocessor library

    http://www.boost.org/

    Copyright (c) 2001-2005 Hartmut Kaiser. Distributed under the Boost
    Software License, Version 1.0. (See accompanying file
    LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
=============================================================================*/

#if !defined(CPP_EXCEPTIONS_HPP_5190E447_A781_4521_A275_5134FF9917D7_INCLUDED)
#define CPP_EXCEPTIONS_HPP_5190E447_A781_4521_A275_5134FF9917D7_INCLUDED

#include <exception>
#include <string>

#include <boost/assert.hpp>
#include <boost/config.hpp>

#include <boost/wave/wave_config.hpp>

///////////////////////////////////////////////////////////////////////////////
// helper macro for throwing exceptions
#if !defined(BOOST_WAVE_THROW)
#ifdef BOOST_NO_STRINGSTREAM
#include <strstream>
#define BOOST_WAVE_THROW(cls, code, msg, act_pos)                             \
    {                                                                         \
    using namespace boost::wave;                                              \
    std::strstream stream;                                                    \
        stream << cls::severity_text(cls::code) << ": "                       \
        << cls::error_text(cls::code);                                        \
    if ((msg)[0] != 0) stream << ": " << (msg);                               \
    stream << std::ends;                                                      \
    std::string throwmsg = stream.str(); stream.freeze(false);                \
    throw cls(throwmsg.c_str(), cls::code, (act_pos).get_line(),              \
        (act_pos).get_column(), (act_pos).get_file().c_str());                \
    }                                                                         \
    /**/
#else
#include <sstream>
#define BOOST_WAVE_THROW(cls, code, msg, act_pos)                             \
    {                                                                         \
    using namespace boost::wave;                                              \
    std::stringstream stream;                                                 \
        stream << cls::severity_text(cls::code) << ": "                       \
        << cls::error_text(cls::code);                                        \
    if ((msg)[0] != 0) stream << ": " << (msg);                               \
    stream << std::ends;                                                      \
    throw cls(stream.str().c_str(), cls::code, (act_pos).get_line(),          \
        (act_pos).get_column(), (act_pos).get_file().c_str());                \
    }                                                                         \
    /**/
#endif // BOOST_NO_STRINGSTREAM
#endif // BOOST_WAVE_THROW

///////////////////////////////////////////////////////////////////////////////
namespace boost {
namespace wave {

///////////////////////////////////////////////////////////////////////////////
// exception severity
namespace util {

    enum severity {
        severity_remark = 0,
        severity_warning,
        severity_error,
        severity_fatal,
        severity_commandline_error
    };
    
    inline char const *
    get_severity(severity level) 
    {
        static char const *severity_text[] = 
        {
            "remark",           // severity_remark
            "warning",          // severity_warning
            "error",            // severity_error
            "fatal error",      // severity_fatal
            "command line error"    // severity_commandline_error
        };
        BOOST_ASSERT(severity_remark <= level && 
            level <= severity_commandline_error);
        return severity_text[level];
    }
}

///////////////////////////////////////////////////////////////////////////////
//  cpp_exception, the base class for all specific C preprocessor exceptions 
class cpp_exception
:   public std::exception
{
public:
    cpp_exception(int line_, int column_, char const *filename_) throw() 
    :   line(line_), column(column_) 
    {
        unsigned int off = 0;
        while (off < sizeof(filename) && *filename_)
            filename[off++] = *filename_++;
        filename[off] = 0;
    }
    ~cpp_exception() throw() {}
    
    virtual char const *what() const throw() = 0;   // to be overloaded
    virtual char const *description() const throw() = 0;
    
    int line_no() const throw() { return line; }
    int column_no() const throw() { return column; }
    char const *file_name() const throw() { return filename; }
    
protected:
    char filename[512];
    int line;
    int column;
};

///////////////////////////////////////////////////////////////////////////////
// preprocessor error
class preprocess_exception :
    public cpp_exception
{
public:
    enum error_code {
        unexpected_error = 0,
        macro_redefinition,
        macro_insertion_error,
        bad_include_file,
        bad_include_statement,
        ill_formed_directive,
        error_directive,
        warning_directive,
        ill_formed_expression,
        missing_matching_if,
        missing_matching_endif,
        ill_formed_operator,
        bad_define_statement,
        too_few_macroarguments,
        too_many_macroarguments,
        empty_macroarguments,
        improperly_terminated_macro,
        bad_line_statement,
        bad_undefine_statement,
        bad_macro_definition,
        illegal_redefinition,
        duplicate_parameter_name,
        invalid_concat,
        last_line_not_terminated,
        ill_formed_pragma_option,
        include_nesting_too_deep,
        misplaced_operator,
        alreadydefined_name,
        undefined_macroname,
        invalid_macroname,
        unexpected_qualified_name,
        division_by_zero,
        integer_overflow,
        illegal_operator_redefinition,
        ill_formed_integer_literal,
        ill_formed_character_literal,
        unbalanced_if_endif,
        character_literal_out_of_range
    };

    preprocess_exception(char const *what_, error_code code, int line_, 
        int column_, char const *filename_) throw() 
    :   cpp_exception(line_, column_, filename_), level(severity_level(code))
    {
        unsigned int off = 0;
        while (off < sizeof(buffer) && *what_)
            buffer[off++] = *what_++;
        buffer[off] = 0;
    }
    ~preprocess_exception() throw() {}
    
    virtual char const *what() const throw()
    {
        return "boost::wave::preprocess_exception";
    }
    virtual char const *description() const throw()
    {
        return buffer;
    }
    util::severity get_severity()
    {
        return level;
    }

    static char const *error_text(int code)
    {
    // error texts in this array must appear in the same order as the items in
    // the error enum above
        static char const *preprocess_exception_errors[] = {
            "unexpected error (should not happen)",     // unexpected_error
            "illegal macro redefinition",               // macro_redefinition
            "macro definition failed (out of memory?)", // macro_insertion_error
            "could not find include file",              // bad_include_file
            "ill formed #include directive",            // bad_include_statement
            "ill formed preprocessor directive",        // ill_formed_directive
            "encountered #error directive or #pragma wave stop()", // error_directive
            "encountered #warning directive",           // warning_directive
            "ill formed preprocessor expression",       // ill_formed_expression
            "the #if for this directive is missing",    // missing_matching_if
            "detected at least one missing #endif directive",   // missing_matching_endif
            "ill formed preprocessing operator",        // ill_formed_operator
            "ill formed #define directive",             // bad_define_statement
            "too few macro arguments",                  // too_few_macroarguments
            "too many macro arguments",                 // too_many_macroarguments
            "empty macro arguments are not supported in pure C++ mode, "
            "use variadics mode to allow these",        // empty_macroarguments
            "improperly terminated macro invocation "
            "or replacement-list terminates in partial "
            "macro expansion (not supported yet)",      // improperly_terminated_macro
            "ill formed #line directive",               // bad_line_statement
            "#undef may not be used on this predefined name",   // bad_undefine_statement
            "invalid macro definition",                 // bad_macro_definition
            "this predefined name may not be redefined",    // illegal_redefinition
            "duplicate macro parameter name",           // duplicate_parameter_name
            "pasting the following two tokens does not "
            "give a valid preprocessing token",         // invalid_concat
            "last line of file ends without a newline", // last_line_not_terminated
            "unknown or illformed pragma option",       // ill_formed_pragma_option
            "include files nested too deep",            // include_nesting_too_deep
            "misplaced operator defined()",             // misplaced_operator
            "the name is already used in this scope as "
            "a macro or scope name",                    // alreadydefined_name
            "undefined macro or scope name may not be imported", // undefined_macroname
            "ill formed macro name",                    // invalid_macroname
            "qualified names are supported in C++0x mode only",  // unexpected_qualified_name
            "division by zero in preprocessor expression",       // division_by_zero
            "integer overflow in preprocessor expression",       // integer_overflow
            "this macro name cannot be used as a as it is an operator in C++",  // illegal_operator_redefinition
            "ill formed integer literal or integer constant too large",   // ill_formed_integer_literal
            "ill formed character literal",             // ill_formed_character_literal
            "unbalanced #if/#endif in include file",    // unbalanced_if_endif
            "character literal out of range"            // character_literal_out_of_range
        };
        BOOST_ASSERT(unexpected_error <= code && 
            code <= character_literal_out_of_range);
        return preprocess_exception_errors[code];
    }

    static util::severity severity_level(int code)
    {
        static util::severity preprocess_exception_severity[] = {
            util::severity_fatal,              // unexpected_error
            util::severity_warning,            // macro_redefinition
            util::severity_fatal,              // macro_insertion_error
            util::severity_error,              // bad_include_file
            util::severity_error,              // bad_include_statement
            util::severity_error,              // ill_formed_directive
            util::severity_fatal,              // error_directive
            util::severity_warning,            // warning_directive
            util::severity_error,              // ill_formed_expression
            util::severity_error,              // missing_matching_if
            util::severity_error,              // missing_matching_endif
            util::severity_error,              // ill_formed_operator
            util::severity_error,              // bad_define_statement
            util::severity_warning,            // too_few_macroarguments
            util::severity_warning,            // too_many_macroarguments
            util::severity_warning,            // empty_macroarguments
            util::severity_error,              // improperly_terminated_macro
            util::severity_warning,            // bad_line_statement
            util::severity_warning,            // bad_undefine_statement
            util::severity_commandline_error,  // bad_macro_definition
            util::severity_warning,            // illegal_redefinition
            util::severity_error,              // duplicate_parameter_name
            util::severity_error,              // invalid_concat
            util::severity_warning,            // last_line_not_terminated
            util::severity_warning,            // ill_formed_pragma_option
            util::severity_fatal,              // include_nesting_too_deep
            util::severity_error,              // misplaced_operator
            util::severity_error,              // alreadydefined_name
            util::severity_error,              // undefined_macroname
            util::severity_error,              // invalid_macroname
            util::severity_error,              // unexpected_qualified_name
            util::severity_fatal,              // division_by_zero
            util::severity_error,              // integer_overflow
            util::severity_error,              // illegal_operator_redefinition
            util::severity_error,              // ill_formed_integer_literal
            util::severity_error,              // ill_formed_character_literal
            util::severity_warning,            // unbalanced_if_endif
            util::severity_warning             // character_literal_out_of_range
        };
        BOOST_ASSERT(unexpected_error <= code && 
            code <= character_literal_out_of_range);
        return preprocess_exception_severity[code];
    }
    static char const *severity_text(int code)
    {
        return util::get_severity(severity_level(code));
    }

private:
    char buffer[512];
    util::severity level;
};

///////////////////////////////////////////////////////////////////////////////
}   // namespace wave
}   // namespace boost

#endif // !defined(CPP_EXCEPTIONS_HPP_5190E447_A781_4521_A275_5134FF9917D7_INCLUDED)
