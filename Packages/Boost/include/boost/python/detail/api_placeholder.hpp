// Copyright David Abrahams 2002.
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
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
