// Copyright David Abrahams 2004. Distributed under the Boost
// Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#ifndef UNWRAP_WRAPPER_DWA2004723_HPP
# define UNWRAP_WRAPPER_DWA2004723_HPP

# include <boost/python/detail/prefix.hpp>
# include <boost/python/detail/is_wrapper.hpp>
#  if defined(BOOST_PYTHON_NO_SFINAE)
#   include <boost/mpl/eval_if.hpp>
#   include <boost/mpl/identity.hpp>
#  else 
#   include <boost/python/detail/enable_if.hpp>
#  endif 

namespace boost { namespace python { namespace detail { 

#  if defined(BOOST_PYTHON_NO_SFINAE)
template <class T>
struct unwrap_wrapper_helper
{
    typedef typename T::_wrapper_wrapped_type_ type;
};
  
template <class T>
typename mpl::eval_if<is_wrapper<T>,unwrap_wrapper_helper<T>,mpl::identity<T> >::type*
unwrap_wrapper(T*)
{
    return 0;
}
#  else 
template <class T>
typename disable_if_ret<is_wrapper<T>,T*>::type
unwrap_wrapper(T*)
{
    return 0;
}

template <class T>
T* unwrap_wrapper(wrapper<T>*)
{
    return 0;
}
#  endif 

}}} // namespace boost::python::detail

#endif // UNWRAP_WRAPPER_DWA2004723_HPP
