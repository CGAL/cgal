//  (C) Copyright David Abrahams 2000. Permission to copy, use, modify, sell and
//  distribute this software is granted provided this copyright notice appears
//  in all copies. This software is provided "as is" without express or implied
//  warranty, and with no claim as to its suitability for any purpose.
//
//  The author gratefully acknowleges the support of Dragon Systems, Inc., in
//  producing this work.

#ifndef NONE_DWA_052000_H_
# define NONE_DWA_052000_H_

# include <boost/python/detail/prefix.hpp>

namespace boost { namespace python { namespace detail {

inline PyObject* none() { Py_INCREF(Py_None); return Py_None; }
    
}}} // namespace boost::python::detail

#endif // NONE_DWA_052000_H_
