// Copyright David Abrahams 2002. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef SHARED_PTR_FROM_PYTHON_DWA20021130_HPP
# define SHARED_PTR_FROM_PYTHON_DWA20021130_HPP

# include <boost/python/handle.hpp>
# include <boost/python/converter/shared_ptr_deleter.hpp>

namespace boost { namespace python { namespace converter { 

template <class T>
struct shared_ptr_from_python
{
    shared_ptr_from_python()
    {
        converter::registry::insert(&convertible, &construct, type_id<shared_ptr<T> >());
    }

 private:
    static void* convertible(PyObject* p)
    {
        if (p == Py_None)
            return p;
        
        return converter::get_lvalue_from_python(p, registered<T>::converters);
    }
    
    static void construct(PyObject* source, rvalue_from_python_stage1_data* data)
    {
        void* const storage = ((converter::rvalue_from_python_storage<shared_ptr<T> >*)data)->storage.bytes;
        // Deal with the "None" case.
        if (data->convertible == source)
            new (storage) shared_ptr<T>();
        else
            new (storage) shared_ptr<T>(
                static_cast<T*>(data->convertible),
                shared_ptr_deleter(handle<>(borrowed(source)))
                );
        
        data->convertible = storage;
    }
};

}}} // namespace boost::python::converter

#endif // SHARED_PTR_FROM_PYTHON_DWA20021130_HPP
