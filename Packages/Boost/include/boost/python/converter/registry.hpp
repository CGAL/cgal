//  Copyright David Abrahams 2001. Permission to copy, use,
//  modify, sell and distribute this software is granted provided this
//  copyright notice appears in all copies. This software is provided
//  "as is" without express or implied warranty, and with no claim as
//  to its suitability for any purpose.
#ifndef REGISTRY_DWA20011127_HPP
# define REGISTRY_DWA20011127_HPP
# include <boost/python/type_id.hpp>
# include <boost/python/converter/to_python_function_type.hpp>
# include <boost/python/converter/rvalue_from_python_data.hpp>
# include <boost/python/converter/constructor_function.hpp>
# include <boost/python/converter/convertible_function.hpp>

namespace boost { namespace python { namespace converter {

struct registration;

// This namespace acts as a sort of singleton
namespace registry
{
  // Get the registration corresponding to the type, creating it if neccessary
  BOOST_PYTHON_DECL registration const& lookup(type_info);

  // Return a pointer to the corresponding registration, if one exists
  BOOST_PYTHON_DECL registration const* query(type_info);
  
  BOOST_PYTHON_DECL void insert(to_python_function_t, type_info);

  // Insert an lvalue from_python converter
  BOOST_PYTHON_DECL void insert(void* (*convert)(PyObject*), type_info);

  // Insert an rvalue from_python converter
  BOOST_PYTHON_DECL void insert(
      convertible_function
      , constructor_function
      , type_info
      );
  
  // Insert an rvalue from_python converter at the tail of the
  // chain. Used for implicit conversions
  BOOST_PYTHON_DECL void push_back(
      convertible_function
      , constructor_function
      , type_info
      );
}

}}} // namespace boost::python::converter

#endif // REGISTRY_DWA20011127_HPP
