// Copyright David Abrahams 2002. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef STRING_LITERAL_DWA2002629_HPP
# define STRING_LITERAL_DWA2002629_HPP

# include <cstddef>
# include <boost/type.hpp>
# include <boost/type_traits/array_traits.hpp>
# include <boost/type_traits/same_traits.hpp>
# include <boost/mpl/bool.hpp>

namespace boost { namespace python { namespace detail { 

# ifndef BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION
template <class T>
struct is_string_literal : mpl::false_
{
};

#  if !defined(__MWERKS__) || __MWERKS__ > 0x2407
template <std::size_t n>
struct is_string_literal<char const[n]> : mpl::true_
{
};

#   if (defined(__DECCXX_VER) && __DECCXX_VER <= 60590031) \
  || (defined(__sgi) && defined(_COMPILER_VERSION) && _COMPILER_VERSION <= 730)
// This compiler mistakenly gets the type of string literals as char*
// instead of char[NN].
template <>
struct is_string_literal<char* const> : mpl::true_
{
};
#   endif

#  else

// CWPro7 has trouble with the array type deduction above
template <class T, std::size_t n>
struct is_string_literal<T[n]>
    : is_same<T, char const>
{
};
#  endif 
# else
template <bool is_array = true>
struct string_literal_helper
{
    typedef char (&yes_string_literal)[1];
    typedef char (&no_string_literal)[2];

    template <class T>
    struct apply
    {
        typedef apply<T> self;
        static T x;
        static yes_string_literal check(char const*);
        static no_string_literal check(char*);
        static no_string_literal check(void const volatile*);
        
        BOOST_STATIC_CONSTANT(
            bool, value = sizeof(self::check(x)) == sizeof(yes_string_literal));
        typedef mpl::bool_<value> type;
    };
};

template <>
struct string_literal_helper<false>
{
    template <class T>
    struct apply : mpl::false_
    {
    };
};

template <class T>
struct is_string_literal
    : string_literal_helper<is_array<T>::value>::apply<T>
{
};
# endif

}}} // namespace boost::python::detail

#endif // STRING_LITERAL_DWA2002629_HPP
