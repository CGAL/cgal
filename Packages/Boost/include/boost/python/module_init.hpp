// Copyright David Abrahams 2002. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef MODULE_INIT_DWA20020722_HPP
# define MODULE_INIT_DWA20020722_HPP

# include <boost/python/detail/prefix.hpp>

# ifndef BOOST_PYTHON_MODULE_INIT

namespace boost { namespace python { namespace detail {

BOOST_PYTHON_DECL void init_module(char const* name, void(*)());

}}}

#  if (defined(_WIN32) || defined(__CYGWIN__)) && !defined(BOOST_PYTHON_STATIC_MODULE)

#   define BOOST_PYTHON_MODULE_INIT(name)               \
void init_module_##name();                              \
extern "C" __declspec(dllexport) void init##name()      \
{                                                       \
    boost::python::detail::init_module(                 \
        #name,&init_module_##name);                     \
}                                                       \
void init_module_##name()

#  elif defined(_AIX) && !defined(BOOST_PYTHON_STATIC_MODULE)

#   include <boost/python/detail/aix_init_module.hpp>
#   define BOOST_PYTHON_MODULE_INIT(name)                               \
void init_module_##name();                                              \
extern "C"                                                              \
{                                                                       \
    extern PyObject* _PyImport_LoadDynamicModule(char*, char*, FILE *); \
    void init##name()                                                   \
    {                                                                   \
        boost::python::detail::aix_init_module(                         \
            _PyImport_LoadDynamicModule, #name, &init_module_##name);   \
    }                                                                   \
}                                                                       \
void init_module_##name()

# else

#   define BOOST_PYTHON_MODULE_INIT(name)                               \
void init_module_##name();                                              \
extern "C"  void init##name()                                           \
{                                                                       \
    boost::python::detail::init_module(#name, &init_module_##name);     \
}                                                                       \
void init_module_##name()

#  endif

# endif 

#endif // MODULE_INIT_DWA20020722_HPP
