// Copyright David Abrahams 2002. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef VOID_PTR_DWA200239_HPP
# define VOID_PTR_DWA200239_HPP

# include <boost/type_traits/remove_cv.hpp>

namespace boost { namespace python { namespace detail { 

template <class U>
inline U& void_ptr_to_reference(void const volatile* p, U&(*)())
{
    return *(U*)p;
}

template <class T>
inline void write_void_ptr(void const volatile* storage, void* ptr, T*)
{
    *(T**)storage = (T*)ptr;
}

// writes U(ptr) into the storage
template <class U>
inline void write_void_ptr_reference(void const volatile* storage, void* ptr, U&(*)())
{
    // stripping CV qualification suppresses warnings on older EDGs
    typedef typename remove_cv<U>::type u_stripped; 
    write_void_ptr(storage, ptr, u_stripped(0));
}

}}} // namespace boost::python::detail

#endif // VOID_PTR_DWA200239_HPP
