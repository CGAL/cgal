// Copyright David Abrahams 2002. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef COPY_CONST_REFERENCE_DWA2002131_HPP
# define COPY_CONST_REFERENCE_DWA2002131_HPP

# include <boost/python/detail/prefix.hpp>
# include <boost/python/detail/indirect_traits.hpp>
# include <boost/mpl/if.hpp>
# include <boost/python/to_python_value.hpp>

namespace boost { namespace python { 

namespace detail
{
  template <class R>
  struct copy_const_reference_expects_a_const_reference_return_type
# if defined(__GNUC__) && __GNUC__ >= 3 || defined(__EDG__)
  {}
# endif
  ;
}

template <class T> struct to_python_value;

struct copy_const_reference
{
    template <class T>
    struct apply
    {
        typedef typename mpl::if_c<
            detail::is_reference_to_const<T>::value
            , to_python_value<T>
            , detail::copy_const_reference_expects_a_const_reference_return_type<T>
        >::type type;
    };
};


}} // namespace boost::python

#endif // COPY_CONST_REFERENCE_DWA2002131_HPP
