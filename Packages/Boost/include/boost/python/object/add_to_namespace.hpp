// Copyright David Abrahams 2002. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef ADD_TO_NAMESPACE_DWA200286_HPP
# define ADD_TO_NAMESPACE_DWA200286_HPP

# include <boost/python/object_fwd.hpp>

namespace boost { namespace python { namespace objects { 

//
// A setattr that's "smart" about function overloading (and docstrings).
//
BOOST_PYTHON_DECL void add_to_namespace(
    object const& name_space, char const* name, object const& attribute);

BOOST_PYTHON_DECL void add_to_namespace(
    object const& name_space, char const* name, object const& attribute, char const* doc);

}}} // namespace boost::python::objects

#endif // ADD_TO_NAMESPACE_DWA200286_HPP
