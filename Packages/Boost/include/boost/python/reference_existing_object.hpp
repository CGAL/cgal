// Copyright David Abrahams 2002. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef REFERENCE_EXISTING_OBJECT_DWA200222_HPP
# define REFERENCE_EXISTING_OBJECT_DWA200222_HPP

# include <boost/python/detail/prefix.hpp>
# include <boost/python/detail/indirect_traits.hpp>
# include <boost/mpl/if.hpp>
# include <boost/python/to_python_indirect.hpp>
# include <boost/type_traits/composite_traits.hpp>

namespace boost { namespace python { 

namespace detail
{
  template <class R>
  struct reference_existing_object_requires_a_pointer_or_reference_return_type
# if defined(__GNUC__) && __GNUC__ >= 3 || defined(__EDG__)
  {}
# endif
  ;
}

template <class T> struct to_python_value;

struct reference_existing_object
{
    template <class T>
    struct apply
    {
        BOOST_STATIC_CONSTANT(
            bool, ok = is_pointer<T>::value || is_reference<T>::value);
        
        typedef typename mpl::if_c<
            ok
            , to_python_indirect<T, detail::make_reference_holder>
            , detail::reference_existing_object_requires_a_pointer_or_reference_return_type<T>
        >::type type;
    };
};

}} // namespace boost::python

#endif // REFERENCE_EXISTING_OBJECT_DWA200222_HPP
