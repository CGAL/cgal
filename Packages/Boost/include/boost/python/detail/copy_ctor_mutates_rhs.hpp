// Copyright David Abrahams 2003. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef COPY_CTOR_MUTATES_RHS_DWA2003219_HPP
# define COPY_CTOR_MUTATES_RHS_DWA2003219_HPP

#include <boost/python/detail/is_auto_ptr.hpp>
#include <boost/mpl/bool.hpp>

namespace boost { namespace python { namespace detail { 

template <class T>
struct copy_ctor_mutates_rhs
    : is_auto_ptr<T>
{
};

}}} // namespace boost::python::detail

#endif // COPY_CTOR_MUTATES_RHS_DWA2003219_HPP
