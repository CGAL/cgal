// Copyright David Abrahams 2002. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef PYOBJECT_TYPE_DWA2002720_HPP
# define PYOBJECT_TYPE_DWA2002720_HPP

# include <boost/python/cast.hpp>

namespace boost { namespace python { namespace converter { 

BOOST_PYTHON_DECL PyObject* checked_downcast_impl(PyObject*, PyTypeObject*);

// Used as a base class for specializations which need to provide
// Python type checking capability.
template <class Object, PyTypeObject* pytype>
struct pyobject_type
{
    static bool check(PyObject* x)
    {
        return ::PyObject_IsInstance(x, (PyObject*)pytype);
    }

    static Object* checked_downcast(PyObject* x)
    {
        return python::downcast<Object>(
            (checked_downcast_impl)(x, pytype)
            );
    }
};

}}} // namespace boost::python::converter

#endif // PYOBJECT_TYPE_DWA2002720_HPP
