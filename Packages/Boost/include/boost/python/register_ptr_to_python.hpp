// Copyright David Abrahams 2002. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.     
#ifndef REGISTER_PTR_TO_PYTHON_HPP
#define REGISTER_PTR_TO_PYTHON_HPP

#include <boost/python/pointee.hpp>
#include <boost/python/object.hpp>
#include <boost/python/object/class_wrapper.hpp>

namespace boost { namespace python {
    
template <class P>
void register_ptr_to_python(BOOST_EXPLICIT_TEMPLATE_TYPE(P))
{
    typedef typename boost::python::pointee<P>::type X;
    objects::class_value_wrapper<
        P
      , objects::make_ptr_instance<
            X
          , objects::pointer_holder<P,X>
        >
    >();
}           

}} // namespace boost::python

#endif // REGISTER_PTR_TO_PYTHON_HPP


