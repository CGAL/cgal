// Copyright David Abrahams 2002. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef MANAGE_NEW_OBJECT_DWA200222_HPP
# define MANAGE_NEW_OBJECT_DWA200222_HPP

# include <boost/python/detail/prefix.hpp>
# include <boost/python/detail/indirect_traits.hpp>
# include <boost/mpl/if.hpp>
# include <boost/python/to_python_indirect.hpp>
# include <boost/type_traits/composite_traits.hpp>

namespace boost { namespace python { 

namespace detail
{
  template <class R>
  struct manage_new_object_requires_a_pointer_return_type
# if defined(__GNUC__) && __GNUC__ >= 3 || defined(__EDG__)
  {}
# endif
  ;
}

struct manage_new_object
{
    template <class T>
    struct apply
    {
        typedef typename mpl::if_c<
            boost::is_pointer<T>::value
            , to_python_indirect<T, detail::make_owning_holder>
            , detail::manage_new_object_requires_a_pointer_return_type<T>
        >::type type;
    };
};

}} // namespace boost::python

#endif // MANAGE_NEW_OBJECT_DWA200222_HPP
