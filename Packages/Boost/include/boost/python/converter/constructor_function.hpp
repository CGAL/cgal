// Copyright David Abrahams 2002. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef CONSTRUCTOR_FUNCTION_DWA200278_HPP
# define CONSTRUCTOR_FUNCTION_DWA200278_HPP

namespace boost { namespace python { namespace converter { 

// Declares the type of functions used to construct C++ objects for
// rvalue from_python conversions.
struct rvalue_from_python_stage1_data;
typedef void (*constructor_function)(PyObject* source, rvalue_from_python_stage1_data*);

}}} // namespace boost::python::converter

#endif // CONSTRUCTOR_FUNCTION_DWA200278_HPP
