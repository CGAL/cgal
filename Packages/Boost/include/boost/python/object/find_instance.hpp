// Copyright David Abrahams 2002. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef FIND_INSTANCE_DWA2002312_HPP
# define FIND_INSTANCE_DWA2002312_HPP

# include <boost/python/type_id.hpp>

namespace boost { namespace python { namespace objects { 

// Given a type_id, find the instance data which corresponds to it, or
// return 0 in case no such type is held.
BOOST_PYTHON_DECL void* find_instance_impl(PyObject*, type_info);

}}} // namespace boost::python::objects

#endif // FIND_INSTANCE_DWA2002312_HPP
