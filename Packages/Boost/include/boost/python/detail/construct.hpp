// Copyright David Abrahams 2002. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef CONSTRUCT_REFERENCE_DWA2002716_HPP
# define CONSTRUCT_REFERENCE_DWA2002716_HPP

namespace boost { namespace python { namespace detail { 

template <class T, class Arg>
void construct_pointee(void* storage, Arg& x
# if !defined(BOOST_MSVC) || BOOST_MSVC > 1300
                       , T const volatile*
# else 
                       , T const*
# endif 
    )
{
    new (storage) T(x);
}

template <class T, class Arg>
void construct_referent_impl(void* storage, Arg& x, T&(*)())
{
    construct_pointee(storage, x, (T*)0);
}

template <class T, class Arg>
void construct_referent(void* storage, Arg const& x, T(*tag)() = 0)
{
    construct_referent_impl(storage, x, tag);
}

template <class T, class Arg>
void construct_referent(void* storage, Arg& x, T(*tag)() = 0)
{
    construct_referent_impl(storage, x, tag);
}

}}} // namespace boost::python::detail

#endif // CONSTRUCT_REFERENCE_DWA2002716_HPP
