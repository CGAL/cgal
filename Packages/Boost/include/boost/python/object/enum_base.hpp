// Copyright David Abrahams 2002. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef ENUM_BASE_DWA200298_HPP
# define ENUM_BASE_DWA200298_HPP

# include <boost/python/object_core.hpp>
# include <boost/python/type_id.hpp>
# include <boost/python/converter/to_python_function_type.hpp>
# include <boost/python/converter/convertible_function.hpp>
# include <boost/python/converter/constructor_function.hpp>

namespace boost { namespace python { namespace objects { 

struct BOOST_PYTHON_DECL enum_base : python::api::object
{
 protected:
    enum_base(
        char const* name
        , converter::to_python_function_t
        , converter::convertible_function
        , converter::constructor_function
        , type_info);

    void add_value(char const* name, long value);
    void export_values();
    
    static PyObject* to_python(PyTypeObject* type, long x);
};

}}} // namespace boost::python::object

#endif // ENUM_BASE_DWA200298_HPP
