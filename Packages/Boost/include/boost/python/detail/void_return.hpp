// Copyright David Abrahams 2002. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef VOID_RETURN_DWA200274_HPP
# define VOID_RETURN_DWA200274_HPP

# include <boost/config.hpp>

namespace boost { namespace python { namespace detail { 

struct void_return
{
    void_return() {}
 private: 
    void operator=(void_return const&);
};

template <class T>
struct returnable
{
    typedef T type;
};

# ifdef BOOST_NO_VOID_RETURNS
template <>
struct returnable<void>
{
    typedef void_return type;
};

#  ifndef BOOST_NO_CV_VOID_SPECIALIZATIONS
template <> struct returnable<const void> : returnable<void> {};
template <> struct returnable<volatile void> : returnable<void> {};
template <> struct returnable<const volatile void> : returnable<void> {};
#  endif

# endif // BOOST_NO_VOID_RETURNS

}}} // namespace boost::python::detail

#endif // VOID_RETURN_DWA200274_HPP
