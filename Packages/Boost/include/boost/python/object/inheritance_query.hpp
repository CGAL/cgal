// Copyright David Abrahams 2003. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef INHERITANCE_QUERY_DWA2003520_HPP
# define INHERITANCE_QUERY_DWA2003520_HPP

# include <boost/python/type_id.hpp>

namespace boost { namespace python { namespace objects {

BOOST_PYTHON_DECL void* find_static_type(void* p, type_info src, type_info dst);
BOOST_PYTHON_DECL void* find_dynamic_type(void* p, type_info src, type_info dst);

}}} // namespace boost::python::object

#endif // INHERITANCE_QUERY_DWA2003520_HPP
