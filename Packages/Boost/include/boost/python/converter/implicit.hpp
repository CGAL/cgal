// Copyright David Abrahams 2002. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef IMPLICIT_DWA2002326_HPP
# define IMPLICIT_DWA2002326_HPP

# include <boost/python/converter/rvalue_from_python_data.hpp>
# include <boost/python/converter/registrations.hpp>
# include <boost/python/converter/registered.hpp>

# include <boost/python/extract.hpp>

namespace boost { namespace python { namespace converter { 

template <class Source, class Target>
struct implicit
{
    static void* convertible(PyObject* obj)
    {
        // Find a converter which can produce a Source instance from
        // obj. The user has told us that Source can be converted to
        // Target, and instantiating construct() below, ensures that
        // at compile-time.
        return implicit_rvalue_convertible_from_python(obj, registered<Source>::converters)
            ? obj : 0;
    }
      
    static void construct(PyObject* obj, rvalue_from_python_stage1_data* data)
    {
        void* storage = ((rvalue_from_python_storage<Target>*)data)->storage.bytes;
        
        new (storage) Target(extract<Source>(obj)());
        
        // record successful construction
        data->convertible = storage;
    }
};

}}} // namespace boost::python::converter

#endif // IMPLICIT_DWA2002326_HPP
