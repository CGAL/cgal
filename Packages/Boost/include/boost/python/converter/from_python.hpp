// Copyright David Abrahams 2002. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef FIND_FROM_PYTHON_DWA2002223_HPP
# define FIND_FROM_PYTHON_DWA2002223_HPP

# include <boost/python/detail/prefix.hpp>
# include <boost/python/converter/rvalue_from_python_data.hpp>

namespace boost { namespace python { namespace converter { 

struct registration;


BOOST_PYTHON_DECL void* get_lvalue_from_python(
    PyObject* source, registration const&);

BOOST_PYTHON_DECL bool implicit_rvalue_convertible_from_python(
    PyObject* source, registration const&);

BOOST_PYTHON_DECL rvalue_from_python_stage1_data rvalue_from_python_stage1(
    PyObject* source, registration const&);

BOOST_PYTHON_DECL void* rvalue_from_python_stage2(
    PyObject* source, rvalue_from_python_stage1_data&, registration const&);

BOOST_PYTHON_DECL void* rvalue_result_from_python(
    PyObject*, rvalue_from_python_stage1_data&);

BOOST_PYTHON_DECL void* reference_result_from_python(PyObject*, registration const&);
BOOST_PYTHON_DECL void* pointer_result_from_python(PyObject*, registration const&);

BOOST_PYTHON_DECL void void_result_from_python(PyObject*);

BOOST_PYTHON_DECL void throw_no_pointer_from_python(PyObject*, registration const&);
BOOST_PYTHON_DECL void throw_no_reference_from_python(PyObject*, registration const&);

}}} // namespace boost::python::converter

#endif // FIND_FROM_PYTHON_DWA2002223_HPP
