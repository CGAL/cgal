// Copyright David Abrahams 2003. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef IS_AUTO_PTR_DWA2003224_HPP
# define IS_AUTO_PTR_DWA2003224_HPP

# ifndef BOOST_NO_AUTO_PTR
#  include <boost/python/detail/is_xxx.hpp>
#  include <memory>
# endif

namespace boost { namespace python { namespace detail { 

# if !defined(BOOST_NO_AUTO_PTR)

BOOST_PYTHON_IS_XXX_DEF(auto_ptr, std::auto_ptr, 1)

# else

template <class T>
struct is_auto_ptr : mpl::false_
{
};

# endif
    
}}} // namespace boost::python::detail

#endif // IS_AUTO_PTR_DWA2003224_HPP
