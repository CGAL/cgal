// Copyright David Abrahams 2002. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef CLASS_DETAIL_DWA200295_HPP
# define CLASS_DETAIL_DWA200295_HPP

# include <boost/python/handle.hpp>
# include <boost/python/type_id.hpp>

namespace boost { namespace python { namespace objects { 

BOOST_PYTHON_DECL type_handle registered_class_object(type_info id);
BOOST_PYTHON_DECL type_handle class_metatype();
BOOST_PYTHON_DECL type_handle class_type();

}}} // namespace boost::python::object

#endif // CLASS_DETAIL_DWA200295_HPP
