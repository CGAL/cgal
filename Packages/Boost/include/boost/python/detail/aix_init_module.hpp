// Copyright David Abrahams 2002. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef AIX_INIT_MODULE_DWA2002529_HPP
# define AIX_INIT_MODULE_DWA2002529_HPP
# ifdef _AIX
# include <boost/python/detail/prefix.hpp>
# include <cstdio>
# ifdef __KCC
#  include <iostream> // this works around a problem in KCC 4.0f
# endif 

namespace boost { namespace python { namespace detail { 

extern "C"
{
    typedef PyObject* (*so_load_function)(char*,char*,FILE*);
}

void aix_init_module(so_load_function, char const* name, void (*init_module)());

}}} // namespace boost::python::detail
# endif

#endif // AIX_INIT_MODULE_DWA2002529_HPP
