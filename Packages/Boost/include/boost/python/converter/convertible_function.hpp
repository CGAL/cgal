// Copyright David Abrahams 2002. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef CONVERTIBLE_FUNCTION_DWA200278_HPP
# define CONVERTIBLE_FUNCTION_DWA200278_HPP

namespace boost { namespace python { namespace converter { 

typedef void* (*convertible_function)(PyObject*);
    
}}} // namespace boost::python::converter

#endif // CONVERTIBLE_FUNCTION_DWA200278_HPP
