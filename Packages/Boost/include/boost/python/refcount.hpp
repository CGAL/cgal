// Copyright David Abrahams 2002. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef REFCOUNT_DWA2002615_HPP
# define REFCOUNT_DWA2002615_HPP

# include <boost/python/detail/prefix.hpp>
# include <boost/python/cast.hpp>

namespace boost { namespace python { 

template <class T>
inline T* incref(T* p)
{
    Py_INCREF(python::upcast<PyObject>(p));
    return p;
}

template <class T>
inline T* xincref(T* p)
{
    Py_XINCREF(python::upcast<PyObject>(p));
    return p;
}

template <class T>
inline void decref(T* p)
{
    Py_DECREF(python::upcast<PyObject>(p));
}

template <class T>
inline void xdecref(T* p)
{
    Py_XDECREF(python::upcast<PyObject>(p));
}

}} // namespace boost::python

#endif // REFCOUNT_DWA2002615_HPP
