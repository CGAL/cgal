// Copyright David Abrahams 2002. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef DECREF_GUARD_DWA20021220_HPP
# define DECREF_GUARD_DWA20021220_HPP

namespace boost { namespace python { namespace detail { 

struct decref_guard
{
    decref_guard(PyObject* o) : obj(o) {}
    ~decref_guard() { Py_XDECREF(obj); }
    void cancel() { obj = 0; }
 private:
    PyObject* obj;
};

}}} // namespace boost::python::detail

#endif // DECREF_GUARD_DWA20021220_HPP
