// Copyright David Abrahams 2003. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef IS_SHARED_PTR_DWA2003224_HPP
# define IS_SHARED_PTR_DWA2003224_HPP

# include <boost/python/detail/is_xxx.hpp>
# include <boost/shared_ptr.hpp>

namespace boost { namespace python { namespace detail { 

BOOST_PYTHON_IS_XXX_DEF(shared_ptr, shared_ptr, 1)
    
}}} // namespace boost::python::detail

#endif // IS_SHARED_PTR_DWA2003224_HPP
