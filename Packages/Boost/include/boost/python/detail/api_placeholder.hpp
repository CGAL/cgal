// Copyright David Abrahams 2002. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef BOOST_PYTHON_API_PLACE_HOLDER_HPP
#define BOOST_PYTHON_API_PLACE_HOLDER_HPP

namespace boost { namespace python {

    inline long len(object const& obj)
    {
        long result = PyObject_Length(obj.ptr());
        if (PyErr_Occurred()) throw_error_already_set();
        return result;
    }
}} // namespace boost::python

#endif // BOOST_PYTHON_API_PLACE_HOLDER_HPP
