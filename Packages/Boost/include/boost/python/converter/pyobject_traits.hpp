// Copyright David Abrahams 2002. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef PYOBJECT_TRAITS_DWA2002720_HPP
# define PYOBJECT_TRAITS_DWA2002720_HPP

# include <boost/python/detail/prefix.hpp>
# include <boost/python/converter/pyobject_type.hpp>

namespace boost { namespace python { namespace converter { 

template <class> struct pyobject_traits;

template <>
struct pyobject_traits<PyObject>
{
    // All objects are convertible to PyObject
    static bool check(PyObject*) { return true; }
    static PyObject* checked_downcast(PyObject* x) { return x; }
};

//
// Specializations
//

# define BOOST_PYTHON_BUILTIN_OBJECT_TRAITS(T)                  \
    template <> struct pyobject_traits<Py##T##Object>           \
        : pyobject_type<Py##T##Object, &Py##T##_Type> {}

// This is not an exhaustive list; should be expanded.
BOOST_PYTHON_BUILTIN_OBJECT_TRAITS(Type);
BOOST_PYTHON_BUILTIN_OBJECT_TRAITS(List);
BOOST_PYTHON_BUILTIN_OBJECT_TRAITS(Int);
BOOST_PYTHON_BUILTIN_OBJECT_TRAITS(Long);
BOOST_PYTHON_BUILTIN_OBJECT_TRAITS(Dict);
BOOST_PYTHON_BUILTIN_OBJECT_TRAITS(Tuple);

}}} // namespace boost::python::converter

#endif // PYOBJECT_TRAITS_DWA2002720_HPP
