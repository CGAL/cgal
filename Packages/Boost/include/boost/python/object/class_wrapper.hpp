// Copyright David Abrahams 2001. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef CLASS_WRAPPER_DWA20011221_HPP
# define CLASS_WRAPPER_DWA20011221_HPP

# include <boost/python/to_python_converter.hpp>
# include <boost/ref.hpp>

namespace boost { namespace python { namespace objects { 

//
// These two classes adapt the static execute function of a class
// MakeInstance execute() function returning a new PyObject*
// reference. The first one is used for class copy constructors, and
// the second one is used to handle smart pointers.
//

template <class Src, class MakeInstance>
struct class_cref_wrapper
    : to_python_converter<Src,class_cref_wrapper<Src,MakeInstance> >
{
    static PyObject* convert(Src const& x)
    {
        return MakeInstance::execute(boost::ref(x));
    }
};

template <class Src, class MakeInstance>
struct class_value_wrapper
    : to_python_converter<Src,class_value_wrapper<Src,MakeInstance> >
{
    static PyObject* convert(Src x)
    {
        return MakeInstance::execute(x);
    }
};

}}} // namespace boost::python::objects

#endif // CLASS_WRAPPER_DWA20011221_HPP
