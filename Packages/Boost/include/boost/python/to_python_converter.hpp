// Copyright David Abrahams 2002. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef TO_PYTHON_CONVERTER_DWA200221_HPP
# define TO_PYTHON_CONVERTER_DWA200221_HPP

# include <boost/python/detail/prefix.hpp>

# include <boost/python/converter/registry.hpp>
# include <boost/python/converter/as_to_python_function.hpp>
# include <boost/python/type_id.hpp>

namespace boost { namespace python { 

template <class T, class Conversion>
struct to_python_converter
{
    to_python_converter();
};

//
// implementation
//

template <class T, class Conversion>
to_python_converter<T,Conversion>::to_python_converter()
{
    typedef converter::as_to_python_function<
        T, Conversion
        > normalized;
        
    converter::registry::insert(
        &normalized::convert
        , type_id<T>());
}

}} // namespace boost::python

#endif // TO_PYTHON_CONVERTER_DWA200221_HPP
