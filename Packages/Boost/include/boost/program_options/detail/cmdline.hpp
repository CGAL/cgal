// Copyright Vladimir Prus 2002-2004.
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)


#ifndef BOOST_CMDLINE_VP_2003_05_19
#define BOOST_CMDLINE_VP_2003_05_19

#include <boost/program_options/config.hpp>
#include <boost/program_options/errors.hpp>
#include <boost/program_options/cmdline.hpp>

#include <boost/detail/workaround.hpp>

#include <boost/function.hpp>

#include <string>
#include <vector>

namespace boost { namespace program_options { namespace detail {

    /** Command line parser class. Main requirements were:
        - Powerful enough to support all common uses.
        - Simple and easy to learn/use.
        - Minimal code size and external dependencies.
        - Extensible for custom syntaxes.

        First all options are registered. After that, elements of command line
        are extracted using operator++. 

        For each element, user can find
        - if it's an option or an argument
        - name of the option
        - index of the option
        - option value(s), if any
        
        Sometimes the registered option name is not equal to the encountered
        one, for example, because name abbreviation is supported.  Therefore
        two option names can be obtained: 
        - the registered one 
        - the one found at the command line

        There are lot of style options, which can be used to tune the command
        line parsing. In addition, it's possible to install additional parser
        which will process custom option styles.

        @todo mininal match length for guessing?
    */
    class BOOST_PROGRAM_OPTIONS_DECL cmdline {
    public:

        typedef ::boost::program_options::command_line_style::style_t style_t;

        typedef function1<std::pair<std::string, std::string>, 
                          const std::string&> 
            additional_parser;
        
        /** Constructs a command line parser for (argc, argv) pair. Uses
            style options passed in 'style', which should be binary or'ed values
            of style_t enum. It can also be zero, in which case a "default"
            style will be used. If 'allow_unregistered' is true, then allows 
            unregistered options. They will be assigned index 1 and are
            assumed to have optional parameter.
        */
        cmdline(const std::vector<std::string>& args, int style,
                bool allow_unregistered = false);

        /** @overload */
        cmdline(int argc, const char*const * argv, int style, 
                bool allow_unregistered = false);

        /** Set additional parser. This will be called for each token
            of command line. If first string in pair is not empty,
            then the token is considered matched by this parser,
            and the first string will be considered an option name
            (which can be long or short), while the second will be
            option's parameter (if not empty). 
            Note that additional parser can match only one token.
        */
        void set_additional_parser(additional_parser p);

        /** Registers a new option.
            @param long_name Long name to use. When ending '*', symbols up to
            it give an allowed prefix -- all options starting with it will be
            allowed. The first character may not be '-'. Empty string means no
            long name.
            @param short_name Short name to use. Value of '\0' means no short
            name.
            @param properties Tell about possible parameters
               '|' -- no parameter
               '?' -- optional parameter
               ':' -- required parameter
               '*' -- 0 or more parameters
               '+' -- 1 or more parameters
            @param index A distinguishing value for the option.
            The index will be returned by 'option_index' member function. 
            Indices need not be unique -- e.g. client can set all indices to 0, 
            and use the string value to recognize options.
        */
        void add_option(const std::string& long_name, char short_name,
                        char properties = '|', int index = 0);

        /** @overload */
        void add_option(const char* long_name, char short_name,
                        char properties = '|', int index = 0);


        /** Returns false when nothing more can be extracted */
        operator bool() const;

        /** Advances to the next element on command line. 
            When called for the first time, positions on the first element. */
        cmdline& operator++();

        /// Tells if the current element is option.
        bool at_option() const;

        /// Tells if the current element is argument.
        bool at_argument() const;

        /** Returns the option name. If there's long option name associated with
            this option, it is returned, even if short name was used in command 
            line.Otherwise, the short name given to 'add_option' is returned 
            with a '-' prefix. 
            For purposes of simplicity, '-' is used even when dos-style short 
            option was found.
        */
        const std::string& option_name() const;
        /** Returns the index for the current option. */
        int option_index() const;
        /** Returns the option name as found on the command line. Any symbols
            that introduce the option, or delimit its parameter will be
            stripped. This function allows to work with allowed prefixes, in
            which case 'option_name' will return the prefix specification, and
            full option name should be queried explicitly.
            */
        const std::string& raw_option_name() const;
        /** Returns the first of option values. If there's more than one option
            throws multiple_values. If there are no options, returns empty 
            string. */
        const std::string& option_value() const;
        /** Returns all option values. */
        const std::vector<std::string>& option_values() const;
        /** Returns the argument. */
        const std::string& argument() const;

        /** Returns all arguments read by this command line parser. */
        const std::vector<std::string>& arguments() const;

        /** Returns the token that was current when 'operator++' was
            last invoked. */
        const std::string& last() const;

    private:

        // Properties of an option.
        enum properties_t {
            no_parameter = 1,
            /// 0 or 1 parameter
            allow_parameter,
            /// exactly 1 parameter
            require_parameter,
            /// 0 or more parameters
            allow_parameters,
            /// 1 or more parameters 
            require_parameters
        };

        enum element_kind_t {
            ek_option = 0,
            ek_argument
        };

        // General error status.
        enum error_type_t {
            no_error = 0,
            unknown_option,
            ambiguous_option,
            invalid_syntax
        };

        // Detailed error status.
        enum error_description_t {
            ed_success = 0,
            ed_long_not_allowed,
            ed_long_adjacent_not_allowed,
            ed_short_adjacent_not_allowed,
            ed_empty_adjacent_parameter,
            ed_missing_parameter,
            ed_extra_parameter,
            ed_unknown_option,
            ed_ambiguous_option
        };

        // The standard say that nested classes has no access to
        // private member of enclosing class. However, most compilers
        // allow that and it's likely be to allowed in future:
        // http://std.dkuug.dk/jtc1/sc22/wg21/docs/cwg_defects.html#45  
        // For Sun compiler, try using friend declaration.
#if BOOST_WORKAROUND(__SUNPRO_CC, BOOST_TESTED_AT(0x560))
        struct option;
        friend struct option;
#endif

        struct option {
            option(const std::string& long_name, char short_name,
                   properties_t properties, int index)
            : long_name(long_name), short_name(short_name), index(index),
            properties(properties)
            {}

            std::string long_name;
            char short_name;
            int index;
            properties_t properties;
        };


        std::vector<option> options;

        void init(const std::vector<std::string>& args, int style,
                  bool allow_unregistered);
        void check_style(int style) const;
       
        void next();

        const option* find_long_option(const char* name);
        const option* find_short_option(char name);

        enum option_kind { error_option, no_option, long_option, short_option,
                           dos_option };
        option_kind is_option(const char* s);
        // All handle_* member functions take string without any "--" or "-" or "/"
        // They return true and success and set m_num_tokens to the number of
        // tokens that were consumed.
        bool handle_long_option(const char* s);
        bool handle_short_option(const char* s, bool ignore_sticky = false);
        bool handle_dos_option(const char* s);
        // Attempts to apply additional parser to 's'.
        bool handle_additional_parser(const std::pair<std::string, std::string>& p);

        bool process_parameter(const option* opt, bool adjacent_parameter,
                              bool next_parameter);
        void advance(int count);

        /// Converts parameter property character into enum value.
        properties_t translate_property(char p);

        /** Clears the error state. If there were an error, throws appripriate
            exception. */
        void clear_error();

        // Copies of input.
        std::vector<std::string> args;
        style_t style;
        bool allow_unregistered;

        // Current state of parsing.
        unsigned index;
        const char* m_current;
        const char* m_next;
        const char* pending_short_option;
        bool m_no_more_options;
        error_description_t m_error_description;
        element_kind_t m_element_kind;
        int m_option_index;

        // Results of parsing the last option
        std::string m_last;
        const option* m_opt;
        std::string  m_option_name, m_raw_option_name, m_argument;
        std::vector<std::string> m_option_values; 
        int m_num_tokens;
        bool m_disguised_long;

        std::vector<std::string> m_arguments;

        additional_parser m_additional_parser;
    };
    
    void test_cmdline_detail();
    
}}}

#endif

