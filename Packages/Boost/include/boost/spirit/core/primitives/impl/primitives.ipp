/*=============================================================================
    Copyright (c) 1998-2003 Joel de Guzman
    Copyright (c) 2003 Martin Wille
    http://spirit.sourceforge.net/

    Use, modification and distribution is subject to the Boost Software
    License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
    http://www.boost.org/LICENSE_1_0.txt)
=============================================================================*/
#if !defined(BOOST_SPIRIT_PRIMITIVES_IPP)
#define BOOST_SPIRIT_PRIMITIVES_IPP

// This should eventually go to a config file.
#if defined(__GNUC__) && (__GNUC__ < 3) && !defined(_STLPORT_VERSION)
#  ifndef BOOST_SPIRIT_NO_CHAR_TRAITS
#    define BOOST_SPIRIT_NO_CHAR_TRAITS
#  endif
#endif

#include <cctype>
#if !defined(BOOST_NO_CWCTYPE)
#include <cwctype>
#endif

#ifndef BOOST_SPIRIT_NO_CHAR_TRAITS
#  include <string> // char_traits
#endif

#if defined(BOOST_MSVC)
#  pragma warning(disable:4800)
#endif

namespace boost { namespace spirit {

    template <typename DrivedT> struct char_parser;

    namespace impl 
    {
        template <typename IteratorT>
        inline IteratorT
        get_last(IteratorT first)
        {
            while (*first)
                first++;
            return first;
        }

        template<
            typename RT,
            typename IteratorT,
            typename ScannerT>
        inline RT
        string_parser_parse(
            IteratorT str_first,
            IteratorT str_last,
            ScannerT& scan)
        {
            typedef typename ScannerT::iterator_t iterator_t;
            iterator_t saved = scan.first;
            std::size_t slen = str_last - str_first;
    
            while (str_first != str_last)
            {
                if (scan.at_end() || (*str_first != *scan))
                    return scan.no_match();
                ++str_first;
                ++scan;
            }
    
            return scan.create_match(slen, nil_t(), saved, scan.first);
        }

    ///////////////////////////////////////////////////////////////////////////
    //
    // Conversion from char_type to int_type
    //
    ///////////////////////////////////////////////////////////////////////////

#ifndef BOOST_SPIRIT_NO_CHAR_TRAITS
#  define BOOST_SPIRIT_CHAR_TRAITS_NAMESPACE std
#else

        template <typename CharT>
        struct char_traits
        {
            typedef CharT int_type;
            typedef CharT char_type;
        };

        template<>
        struct char_traits<char>
        {
            typedef int int_type;
            typedef char char_type;
    
            static char_type 
            to_char_type(int_type c)
            { 
                return static_cast<char_type>(c); 
            }
    
            static int 
            to_int_type(char c)
            { 
                return static_cast<unsigned char>(c); 
            }
        };

        template<>
        struct char_traits<unsigned char>
        {
            typedef int int_type;
            typedef unsigned char char_type;
    
            static char_type 
            to_char_type(int_type c)
            { 
                return static_cast<char_type>(c); 
            }
    
            static int 
            to_int_type(unsigned char c)
            { 
                return c; 
            }
        };

#  define BOOST_SPIRIT_CHAR_TRAITS_NAMESPACE impl
#  ifndef BOOST_NO_CWCHAR

        template<>
        struct char_traits<wchar_t>
        {
            typedef wint_t int_type;
            typedef wchar_t char_type;
            
            static char_type 
            to_char_type(int_type c)
            { 
                return static_cast<char_type>(c); 
            }
    
            static wint_t 
            to_int_type(wchar_t c)
            { 
                return c; 
            }
        };

#  endif
#endif // BOOST_SPIRIT_NO_CHAR_TRAITS

        template <typename CharT>
        inline typename
        BOOST_SPIRIT_CHAR_TRAITS_NAMESPACE::char_traits<CharT>::int_type
        to_int_type(CharT c)
        {
            return BOOST_SPIRIT_CHAR_TRAITS_NAMESPACE
                ::char_traits<CharT>::to_int_type(c);
        }
    
        template <typename CharT>
        inline CharT
        to_char_type(typename 
            BOOST_SPIRIT_CHAR_TRAITS_NAMESPACE::char_traits<CharT>::int_type c)
        {
            return BOOST_SPIRIT_CHAR_TRAITS_NAMESPACE
                ::char_traits<CharT>::to_char_type(c);
        }

        ///////////////////////////////////////////////////////////////////////
        //
        //  Convenience functions
        //
        ///////////////////////////////////////////////////////////////////////

        inline bool 
        isalnum_(char c)
        { 
            using namespace std; 
            return isalnum(to_int_type(c)); 
        }
    
        inline bool 
        isalpha_(char c)
        { 
            using namespace std; 
            return isalpha(to_int_type(c)); 
        }
    
        inline bool 
        iscntrl_(char c)
        { 
            using namespace std; 
            return iscntrl(to_int_type(c)); 
        }
    
        inline bool 
        isdigit_(char c)
        { 
            using namespace std; 
            return isdigit(to_int_type(c)); 
        }
    
        inline bool 
        isgraph_(char c)
        { 
            using namespace std; 
            return isgraph(to_int_type(c)); 
        }
    
        inline bool 
        islower_(char c)
        { 
            using namespace std; 
            return islower(to_int_type(c)); 
        }
    
        inline bool 
        isprint_(char c)
        { 
            using namespace std; 
            return isprint(to_int_type(c)); 
        }
    
        inline bool 
        ispunct_(char c)
        { 
            using namespace std; 
            return ispunct(to_int_type(c)); 
        }
    
        inline bool 
        isspace_(char c)
        { 
            using namespace std; 
            return isspace(to_int_type(c)); 
        }
    
        inline bool 
        isupper_(char c)
        { 
            using namespace std; 
            return isupper(to_int_type(c)); }
    
        inline bool 
        isxdigit_(char c)
        { 
            using namespace std; 
            return isxdigit(to_int_type(c)); 
        }
    
        inline bool 
        isblank_(char c)
        { 
            return (c == ' ' || c == '\t'); 
        }
        
        inline char 
        tolower_(char c)
        { 
            using namespace std; 
            return to_char_type<char>(tolower(to_int_type(c))); 
        }
    
        inline char 
        toupper_(char c)
        { 
            using namespace std; 
            return to_char_type<char>(toupper(to_int_type(c))); 
        }

#if !defined(BOOST_NO_CWCTYPE)

        inline bool 
        isalnum_(wchar_t c)
        { 
            using namespace std; 
            return iswalnum(to_int_type(c)); 
        }
    
        inline bool 
        isalpha_(wchar_t c)
        { 
            using namespace std; 
            return iswalpha(to_int_type(c)); 
        }
    
        inline bool 
        iscntrl_(wchar_t c)
        { 
            using namespace std; 
            return iswcntrl(to_int_type(c)); 
        }
    
        inline bool 
        isdigit_(wchar_t c)
        { 
            using namespace std; 
            return iswdigit(to_int_type(c)); 
        }
    
        inline bool 
        isgraph_(wchar_t c)
        { 
            using namespace std; 
            return iswgraph(to_int_type(c)); 
        }
    
        inline bool 
        islower_(wchar_t c)
        { 
            using namespace std; 
            return iswlower(to_int_type(c)); 
        }
    
        inline bool 
        isprint_(wchar_t c)
        { 
            using namespace std; 
            return iswprint(to_int_type(c)); 
        }
    
        inline bool 
        ispunct_(wchar_t c)
        { 
            using namespace std; 
            return iswpunct(to_int_type(c)); 
        }
    
        inline bool 
        isspace_(wchar_t c)
        { 
            using namespace std; 
            return iswspace(to_int_type(c)); 
        }
    
        inline bool 
        isupper_(wchar_t c)
        { 
            using namespace std; 
            return iswupper(to_int_type(c)); 
        }
    
        inline bool 
        isxdigit_(wchar_t c)
        { 
            using namespace std; 
            return iswxdigit(to_int_type(c)); 
        }
    
        inline bool 
        isblank_(wchar_t c)
        { 
            return (c == L' ' || c == L'\t'); 
        } 
    
        inline wchar_t 
        tolower_(wchar_t c)
        { 
            using namespace std; 
            return to_char_type<wchar_t>(towlower(to_int_type(c))); 
        }
    
        inline wchar_t 
        toupper_(wchar_t c)
        { 
            using namespace std; 
            return to_char_type<wchar_t>(towupper(to_int_type(c))); 
        }

#endif // !defined(BOOST_NO_CWCTYPE)

}}} // namespace boost::spirit::impl

#endif
