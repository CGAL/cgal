// Copyright Gottfried Ganﬂauge 2003. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.

# ifndef BOOST_PYTHON_DETAIL_DEALLOC_HPP_
# define BOOST_PYTHON_DETAIL_DEALLOC_HPP_
namespace boost { namespace python { namespace detail {
    extern "C"
    {
        inline void dealloc(PyObject* self)
        {
          PyObject_Del(self);
        }
    }
}}} // namespace boost::python::detail
# endif    // BOOST_PYTHON_DETAIL_DEALLOC_HPP_
